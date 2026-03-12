"""
Configuration for PubMed cost-ratio literature search.
Diseases from the Vaccine Impact Modelling Consortium (VIMC).
"""

# ── NCBI E-utilities ────────────────────────────────────────────────────────
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ESEARCH_URL = f"{EUTILS_BASE}/esearch.fcgi"
EFETCH_URL = f"{EUTILS_BASE}/efetch.fcgi"
RATE_LIMIT_DELAY = 0.34  # seconds between requests (< 3 req/s)
EFETCH_BATCH_SIZE = 50
MAX_RESULTS_PER_DISEASE = 200

# ── Disease definitions ─────────────────────────────────────────────────────
# Each entry: key → dict with "name" (display), "synonyms" (PubMed search terms)
DISEASES = {
    "cholera": {
        "name": "Cholera",
        "synonyms": ["cholera", "Vibrio cholerae"],
    },
    "covid": {
        "name": "COVID-19",
        "synonyms": ["COVID-19", "SARS-CoV-2", "coronavirus disease 2019"],
    },
    "typhoid": {
        "name": "Typhoid",
        "synonyms": ["typhoid fever", "Salmonella typhi", "enteric fever"],
    },
    "dengue": {
        "name": "Dengue",
        "synonyms": ["dengue", "dengue fever", "dengue hemorrhagic fever"],
    },
    "hib": {
        "name": "Haemophilus influenzae type b",
        "synonyms": [
            "Haemophilus influenzae type b",
            "Hib disease",
            "Hib meningitis",
        ],
    },
    "rotavirus": {
        "name": "Rotavirus",
        "synonyms": ["rotavirus", "rotavirus gastroenteritis"],
    },
    "pneumococcal": {
        "name": "Pneumococcal disease (PCV)",
        "synonyms": [
            "pneumococcal disease",
            "Streptococcus pneumoniae",
            "pneumococcal pneumonia",
            "pneumococcal meningitis",
        ],
    },
    "hpv": {
        "name": "Human papillomavirus (HPV)",
        "synonyms": [
            "human papillomavirus",
            "HPV",
            "cervical cancer HPV",
        ],
    },
    "japanese_encephalitis": {
        "name": "Japanese encephalitis",
        "synonyms": ["Japanese encephalitis", "JE virus"],
    },
    "malaria": {
        "name": "Malaria",
        "synonyms": ["malaria", "Plasmodium falciparum", "Plasmodium vivax"],
    },
    "measles": {
        "name": "Measles",
        "synonyms": ["measles", "measles virus", "rubeola"],
    },
    "meningitis": {
        "name": "Meningococcal meningitis",
        "synonyms": [
            "meningococcal disease",
            "Neisseria meningitidis",
            "meningococcal meningitis",
        ],
    },
    "rubella": {
        "name": "Rubella",
        "synonyms": ["rubella", "congenital rubella syndrome"],
    },
    "yellow_fever": {
        "name": "Yellow fever",
        "synonyms": ["yellow fever", "yellow fever virus"],
    },
}

# ── Cost / economic search terms ────────────────────────────────────────────
COST_TERMS = [
    "cost of illness",
    "cost-effectiveness",
    "cost effectiveness",
    "economic burden",
    "cost-benefit",
    "cost benefit",
    "economic evaluation",
    "direct medical cost",
    "indirect cost",
    "productivity loss",
    "disability-adjusted life year",
    "DALY",
    "cost per case",
    "hospitalization cost",
    "treatment cost",
    "vaccination cost",
    "vaccine cost",
    "immunization cost",
]

# ── Output paths (relative to project root) ─────────────────────────────────
RAW_CACHE_PATH = "data/pubmed_raw_results.json"
EXTRACTION_PATH = "data/pubmed_cost_extractions.json"
FULLTEXT_DIR = "data/fulltext_pdfs"
FULLTEXT_MAP_PATH = "data/pubmed_fulltext_map.json"
QMD_OUTPUT_PATH = "cost_ratio_literature.qmd"
