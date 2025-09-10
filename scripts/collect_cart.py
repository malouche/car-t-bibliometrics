#!/usr/bin/env python3
"""
Collect CAR-T therapy in lymphoma (2009–2025) records from PubMed,
enrich with OpenAlex, and write a CSV in a fixed schema.

Env vars expected:
- NCBI_API_KEY      (required for polite PubMed rate limits)
- CONTACT_EMAIL     (optional; used in User-Agent headers)

Output:
- data/cart_lymphoma_2009_2025.csv
"""

import os
import sys
import time
import math
import json
import re
import csv
import html
import logging
from datetime import datetime, timezone
from urllib.parse import urlencode, quote

import requests
import pandas as pd
from xml.etree import ElementTree as ET

# ----------------------------- Config ---------------------------------

OUT_PATH = "data/cart_lymphoma_2009_2025.csv"
MAX_PUBMED = 5000               # hard cap to avoid runaway jobs
BATCH_SIZE = 200                # efetch batch size (<= 200 recommended)
ES_RETMAX = 10000               # esearch retmax (PubMed allows up to 100k)
SLEEP = 0.34                    # ~3 requests/sec when no API key; we'll be polite anyway

NCBI_API_KEY = os.getenv("NCBI_API_KEY", "")
CONTACT_EMAIL = os.getenv("CONTACT_EMAIL", "")

USER_AGENT = f"car-t-bibliometrics/1.0 (mailto:{CONTACT_EMAIL})" if CONTACT_EMAIL else "car-t-bibliometrics/1.0"

PUBMED_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OPENALEX_BASE = "https://api.openalex.org"

# PubMed query (lymphoma focus, 2009–2025)
PUBMED_QUERY = r'''(
  "chimeric antigen receptor t"[Title/Abstract] OR "CAR T"[Title/Abstract] OR "CAR-T"[Title/Abstract]
  OR "CAR T-cell*"[Title/Abstract] OR "CAR-T-cell*"[Title/Abstract]
  OR kymriah[Title/Abstract] OR tisagenlecleucel[Title/Abstract]
  OR yescarta[Title/Abstract] OR "axicabtagene ciloleucel"[Title/Abstract]
  OR tecartus[Title/Abstract] OR "brexucabtagene autoleucel"[Title/Abstract]
  OR breyanzi[Title/Abstract] OR "lisocabtagene maraleucel"[Title/Abstract]
)
AND
(
  lymphoma[Title/Abstract] OR "non-hodgkin lymphoma"[Title/Abstract] OR "hodgkin lymphoma"[Title/Abstract]
  OR "diffuse large b-cell lymphoma"[Title/Abstract] OR DLBCL[Title/Abstract]
  OR "follicular lymphoma"[Title/Abstract] OR "mantle cell lymphoma"[Title/Abstract]
  OR "marginal zone lymphoma"[Title/Abstract] OR MZL[Title/Abstract]
  OR "primary mediastinal b-cell lymphoma"[Title/Abstract] OR PMBCL[Title/Abstract]
  OR "small lymphocytic lymphoma"[Title/Abstract] OR SLL[Title/Abstract]
  OR "burkitt lymphoma"[Title/Abstract] OR "peripheral t-cell lymphoma"[Title/Abstract] OR PTCL[Title/Abstract]
  OR "anaplastic large cell lymphoma"[Title/Abstract] OR ALCL[Title/Abstract]
  OR "cutaneous t-cell lymphoma"[Title/Abstract] OR CTCL[Title/Abstract]
)
AND ("2009/01/01"[Date - Publication] : "2025/12/31"[Date - Publication])
AND (english[lang])
AND (Humans[Mesh] OR Humans[Title/Abstract])
AND (Review[ptyp] OR Journal Article[ptyp])'''

# Tagging regexes
PRODUCT_RE = re.compile(
    r"(kymriah|tisagenlecleucel|yescarta|axicabtagene\s+ciloleucel|tecartus|brexucabtagene\s+autoleucel|breyanzi|lisocabtagene\s+maraleucel)",
    re.IGNORECASE,
)
ANTIGEN_MAP = [
    ("CD19", re.compile(r"\bCD19\b", re.IGNORECASE)),
    ("CD20", re.compile(r"\bCD20\b", re.IGNORECASE)),
    ("CD22", re.compile(r"\bCD22\b", re.IGNORECASE)),
    ("CD30", re.compile(r"\bCD30\b", re.IGNORECASE)),
    ("CD7",  re.compile(r"\bCD7\b",  re.IGNORECASE)),
]
DISEASE_MAP = [
    ("DLBCL", re.compile(r"diffuse large b-?cell lymphoma|DLBCL", re.IGNORECASE)),
    ("FL",    re.compile(r"follicular lymphoma", re.IGNORECASE)),
    ("MCL",   re.compile(r"mantle cell lymphoma", re.IGNORECASE)),
    ("MZL",   re.compile(r"marginal zone lymphoma|MZL", re.IGNORECASE)),
    ("PMBCL", re.compile(r"primary mediastinal b-?cell lymphoma|PMBCL", re.IGNORECASE)),
    ("HL",    re.compile(r"hodgkin lymphoma", re.IGNORECASE)),
    ("PTCL",  re.compile(r"peripheral t-?cell lymphoma|PTCL", re.IGNORECASE)),
    ("ALCL",  re.compile(r"anaplastic large cell lymphoma|ALCL", re.IGNORECASE)),
    ("CTCL",  re.compile(r"cutaneous t-?cell lymphoma|CTCL", re.IGNORECASE)),
    ("BL",    re.compile(r"burkitt lymphoma", re.IGNORECASE)),
    ("SLL",   re.compile(r"small lymphocytic lymphoma|SLL", re.IGNORECASE)),
]

CSV_COLUMNS = [
    "doi","pmid","openalex_id","title","abstract","year","month","journal","issn","publisher",
    "document_type","language","authors","author_affiliations","author_countries",
    "corresponding_author","corresponding_author_country",
    "author_keywords","mesh_terms","funding_text","grant_numbers",
    "cited_by_count","reference_count","clinical_trial_flag",
    "product_tag","target_antigen","disease_tag","study_type",
    "source_record","search_date_utc","notes"
]

# ----------------------------- Utils ----------------------------------

def http_get(url, params=None, headers=None, timeout=30):
    hdrs = {"User-Agent": USER_AGENT}
    if headers:
        hdrs.update(headers)
    r = requests.get(url, params=params, headers=hdrs, timeout=timeout)
    r.raise_for_status()
    return r

def chunked(iterable, size):
    for i in range(0, len(iterable), size):
        yield iterable[i:i+size]

def strip_or_na(x):
    return x.strip() if isinstance(x, str) and x.strip() else "NA"

def join_or_na(xs, sep="; "):
    xs = [x for x in xs if x and str(x).strip()]
    return sep.join(xs) if xs else "NA"

# --------------------------- PubMed layer ------------------------------

def pubmed_esearch(query, retmax=ES_RETMAX):
    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": retmax,
        "usehistory": "y",
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    r = http_get(f"{PUBMED_BASE}/esearch.fcgi", params=params)
    data = r.json()
    count = int(data["esearchresult"]["count"])
    webenv = data["esearchresult"].get("webenv")
    query_key = data["esearchresult"].get("querykey")
    idlist = data["esearchresult"].get("idlist", [])
    return count, idlist, webenv, query_key

def pubmed_efetch_ids(id_list):
    """Fetch PubMed records for a list of PMIDs (<= 200). Returns XML root."""
    params = {
        "db": "pubmed",
        "retmode": "xml",
        "id": ",".join(id_list),
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    r = http_get(f"{PUBMED_BASE}/efetch.fcgi", params=params)
    return ET.fromstring(r.text)

def get_text(elem, path):
    node = elem.find(path)
    return node.text if node is not None and node.text is not None else ""

def parse_pubmed_article(article):
    """Extract fields from a PubMedArticle element."""
    med = article.find("./MedlineCitation")
    art = med.find("./Article") if med is not None else None
    if med is None or art is None:
        return None

    pmid = get_text(med, "./PMID")
    title = get_text(art, "./ArticleTitle")
    # Abstract may have multiple parts
    abstract_texts = []
    for ab in art.findall("./Abstract/AbstractText"):
        t = (ab.text or "")
        if "Label" in ab.attrib and ab.attrib["Label"]:
            t = f"{ab.attrib['Label']}: {t}"
        abstract_texts.append(t)
    abstract = " ".join(abstract_texts).strip()

    # Journal info
    journal = get_text(art, "./Journal/Title") or get_text(art, "./Journal/ISOAbbreviation")
    year = get_text(art, "./Journal/JournalIssue/PubDate/Year")
    month = get_text(art, "./Journal/JournalIssue/PubDate/Month")

    # DOI
    doi = ""
    for eid in art.findall("./ELocationID"):
        if eid.attrib.get("EIdType","").lower() == "doi":
            doi = (eid.text or "").strip().lower()
            break

    # Language
    language = get_text(art, "./Language") or "NA"

    # Publication types → document_type + clinical_trial_flag + study_type
    pubtypes = [ (pt.text or "").lower() for pt in art.findall("./PublicationTypeList/PublicationType") ]
    document_type = "Review" if any("review" in x for x in pubtypes) else "Article"

    clinical_trial_flag = "Y" if any("clinical trial" in x for x in pubtypes) else "N"
    study_type = "review" if "review" in " ".join(pubtypes) else ("trial" if "clinical trial" in " ".join(pubtypes) else "real-world")

    # Authors and affiliations
    authors, affils, countries = [], [], []
    for au in art.findall("./AuthorList/Author"):
        last = get_text(au, "./LastName")
        fore = get_text(au, "./ForeName") or get_text(au, "./Initials")
        name = " ".join([fore, last]).strip() if (fore or last) else ""
        if name:
            # Convert to "Last, First"
            if last and fore:
                authors.append(f"{last}, {fore}")
            else:
                authors.append(name)
        aff = [ (a.text or "") for a in au.findall("./AffiliationInfo/Affiliation") ]
        aff_text = "; ".join(aff) if aff else ""
        affils.append(aff_text if aff_text else "NA")
        # Crude country guess: last token after comma
        if aff_text and "," in aff_text:
            cand = aff_text.split(",")[-1].strip()
            countries.append(cand if cand else "NA")
        else:
            countries.append("NA")

    authors_s = join_or_na(authors)
    affils_s = join_or_na(affils)
    countries_s = join_or_na(countries)

    # Corresponding author (best-effort: look for 'corresponding' or email)
    corr_name, corr_country = "NA", "NA"
    for a_name, a_aff in zip(authors, affils):
        if re.search(r"correspond", a_aff, re.IGNORECASE) or re.search(r"@", a_aff):
            corr_name = a_name
            if "," in a_aff:
                corr_country = a_aff.split(",")[-1].strip() or "NA"
            break
    if corr_name == "NA" and authors:
        corr_name = authors[0]
        corr_country = countries[0] if countries else "NA"

    # MeSH terms
    mesh_terms = [ (mh.find("DescriptorName").text or "") for mh in med.findall("./MeshHeadingList/MeshHeading") if mh.find("DescriptorName") is not None ]
    mesh_s = join_or_na(mesh_terms)

    # Funding / grants
    grant_texts, grant_ids = [], []
    for g in art.findall("./GrantList/Grant"):
        gid = get_text(g, "./GrantID")
        ag  = get_text(g, "./Agency")
        if gid or ag:
            grant_texts.append("; ".join([x for x in [ag, gid] if x]))
        if gid:
            grant_ids.append(gid)
    funding_text = join_or_na(grant_texts)
    grant_numbers = join_or_na(grant_ids)

    record = {
        "pmid": pmid,
        "title": html.unescape(title).strip(),
        "abstract": html.unescape(abstract).strip() if abstract else "NA",
        "year": year or "NA",
        "month": month or "NA",
        "journal": journal or "NA",
        "document_type": document_type,
        "language": language or "NA",
        "authors": authors_s,
        "author_affiliations": affils_s,
        "author_countries": countries_s,
        "corresponding_author": corr_name,
        "corresponding_author_country": corr_country,
        "mesh_terms": mesh_s,
        "funding_text": funding_text,
        "grant_numbers": grant_numbers,
        "pubtypes_raw": "; ".join(pubtypes),
        "doi": doi,
    }
    return record

def fetch_pubmed_records(pmids):
    results = []
    for batch in chunked(pmids, BATCH_SIZE):
        root = pubmed_efetch_ids(batch)
        for art in root.findall("./PubmedArticle"):
            rec = parse_pubmed_article(art)
            if rec:
                results.append(rec)
        time.sleep(SLEEP)
    return results

# --------------------------- OpenAlex layer ----------------------------

def openalex_enrich(doi):
    """Return dict: openalex_id, cited_by_count, reference_count, issn, publisher"""
    if not doi or doi == "NA":
        return {}
    # normalize doi (no leading https://doi.org/)
    doi_norm = doi.lower().strip()
    doi_norm = doi_norm.replace("https://doi.org/", "").replace("http://doi.org/", "")
    url = f"{OPENALEX_BASE}/works/doi:{quote(doi_norm)}"
    try:
        r = http_get(url, timeout=30)
        data = r.json()
    except requests.HTTPError as e:
        # try alternate route
        try:
            r = http_get(f"{OPENALEX_BASE}/works", params={"filter": f"doi:{doi_norm}"}, timeout=30)
            js = r.json()
            data = js["results"][0] if js.get("results") else {}
        except Exception:
            data = {}

    if not data:
        return {}

    # Fields changed a bit over time; handle both
    openalex_id = data.get("id", "NA")
    cited_by_count = data.get("cited_by_count", None)
    referenced_works = data.get("referenced_works", None)
    reference_count = len(referenced_works) if isinstance(referenced_works, list) else data.get("reference_count", None)

    issn = "NA"
    publisher = "NA"
    # new keys
    if "primary_location" in data and isinstance(data["primary_location"], dict):
        src = data["primary_location"].get("source", {}) or {}
        if isinstance(src, dict):
            # issn may be list
            issn_val = src.get("issn")
            if isinstance(issn_val, list) and issn_val:
                issn = "; ".join(issn_val)
            elif isinstance(issn_val, str) and issn_val:
                issn = issn_val
    if "host_venue" in data and isinstance(data["host_venue"], dict):
        publisher = data["host_venue"].get("publisher") or publisher

    return {
        "openalex_id": openalex_id or "NA",
        "cited_by_count": cited_by_count if cited_by_count is not None else "NA",
        "reference_count": reference_count if reference_count is not None else "NA",
        "issn": issn if issn else "NA",
        "publisher": publisher if publisher else "NA",
    }

# --------------------------- Tagging logic -----------------------------

def infer_product(text):
    m = PRODUCT_RE.search(text)
    return m.group(1) if m else "NA"

def infer_antigen(text):
    for label, rgx in ANTIGEN_MAP:
        if rgx.search(text):
            return label
    return "Other/NA"

def infer_disease(text):
    for label, rgx in DISEASE_MAP:
        if rgx.search(text):
            return label
    return "Other"

def infer_study_type(pubtypes_raw, abstract):
    if "review" in pubtypes_raw.lower():
        return "review"
    if "clinical trial" in pubtypes_raw.lower():
        return "trial"
    # crude cues
    if re.search(r"retrospective|real[- ]world|registry|cohort", abstract or "", re.IGNORECASE):
        return "real-world"
    return "preclinical" if re.search(r"mouse|murine|preclinical|in vivo|in vitro", abstract or "", re.IGNORECASE) else "real-world"

# ------------------------------ Main ----------------------------------

def main():
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
    logging.info("Starting PubMed esearch…")

    count, idlist, webenv, query_key = pubmed_esearch(PUBMED_QUERY, retmax=min(ES_RETMAX, MAX_PUBMED))
    logging.info(f"PubMed count: {count}, fetched IDs: {len(idlist)}")

    if not idlist:
        logging.error("No PubMed IDs found. Exiting.")
        sys.exit(1)

    # Limit total processed
    pmids = idlist[:min(len(idlist), MAX_PUBMED)]
    logging.info(f"Fetching details for {len(pmids)} PMIDs…")

    records = fetch_pubmed_records(pmids)

    # Filter: must have DOI, and document_type Article/Review already enforced by query, but recheck
    filt = []
    for r in records:
        if r.get("doi") and r.get("document_type") in {"Article", "Review"}:
            filt.append(r)
    logging.info(f"Records with DOI after filtering: {len(filt)}")

    # Deduplicate by DOI
    seen = set()
    dedup = []
    for r in filt:
        doi = r["doi"].lower()
        if doi not in seen:
            seen.add(doi)
            dedup.append(r)
    logging.info(f"After DOI dedup: {len(dedup)}")

    # Enrich with OpenAlex
    out_rows = []
    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    for i, r in enumerate(dedup, 1):
        text_for_tag = " ".join([
            r.get("title",""),
            r.get("abstract",""),
            r.get("mesh_terms",""),
        ])

        # OpenAlex enrichment
        try:
            oa = openalex_enrich(r["doi"])
            time.sleep(SLEEP)
        except Exception as e:
            logging.warning(f"OpenAlex enrich failed for DOI {r['doi']}: {e}")
            oa = {}

        # Study type refinement
        study_type = infer_study_type(r.get("pubtypes_raw",""), r.get("abstract",""))

        row = {
            "doi": r.get("doi","NA"),
            "pmid": r.get("pmid","NA"),
            "openalex_id": oa.get("openalex_id","NA"),
            "title": r.get("title","NA"),
            "abstract": r.get("abstract","NA"),
            "year": r.get("year","NA"),
            "month": r.get("month","NA"),
            "journal": r.get("journal","NA"),
            "issn": oa.get("issn","NA"),
            "publisher": oa.get("publisher","NA"),
            "document_type": r.get("document_type","NA"),
            "language": r.get("language","NA"),
            "authors": r.get("authors","NA"),
            "author_affiliations": r.get("author_affiliations","NA"),
            "author_countries": r.get("author_countries","NA"),
            "corresponding_author": r.get("corresponding_author","NA"),
            "corresponding_author_country": r.get("corresponding_author_country","NA"),
            "author_keywords": "NA",  # PubMed often lacks author keywords
            "mesh_terms": r.get("mesh_terms","NA"),
            "funding_text": r.get("funding_text","NA"),
            "grant_numbers": r.get("grant_numbers","NA"),
            "cited_by_count": oa.get("cited_by_count","NA"),
            "reference_count": oa.get("reference_count","NA"),
            "clinical_trial_flag": "Y" if "clinical trial" in r.get("pubtypes_raw","").lower() else "N",
            "product_tag": infer_product(text_for_tag),
            "target_antigen": infer_antigen(text_for_tag),
            "disease_tag": infer_disease(text_for_tag),
            "study_type": study_type,
            "source_record": "Both" if oa else "PubMed",
            "search_date_utc": now_utc,
            "notes": "NA",
        }
        out_rows.append(row)

        if i % 100 == 0:
            logging.info(f"Processed {i} records…")

    # Write CSV with fixed column order
    df = pd.DataFrame(out_rows, columns=CSV_COLUMNS)
    df.to_csv(OUT_PATH, index=False, quoting=csv.QUOTE_MINIMAL)
    logging.info(f"Wrote {len(df)} rows to {OUT_PATH}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(f"Fatal error: {e}")
        sys.exit(1)
