"""
Stage 2: Extract cost-related data from cached PubMed articles.

Uses regex heuristics to identify monetary amounts, study types,
countries, and cost components. When full-text PDFs have been downloaded
and extracted, uses fulltext; otherwise falls back to abstract only.
Designed for human/AI review — not intended as a definitive extraction.
"""

import json
import re
from pathlib import Path

from .config import RAW_CACHE_PATH, EXTRACTION_PATH

# ── Regex patterns ──────────────────────────────────────────────────────────

# Match dollar amounts: $1,234.56, US$50, USD 1.2 million, etc.
_MONEY_RE = re.compile(
    r"(?:US\s*\$|USD\s*|Int(?:ernational)?\s*\$|\$)\s*"
    r"([\d,]+(?:\.\d+)?)\s*"
    r"(thousand|million|billion|trillion)?",
    re.IGNORECASE,
)

# Study type keywords
_STUDY_TYPES = {
    "CEA": [
        r"cost[\s-]*effective",
        r"incremental cost[\s-]*effectiveness ratio",
        r"ICER",
        r"cost per DALY",
        r"cost per QALY",
    ],
    "COI": [
        r"cost of illness",
        r"economic burden",
        r"disease burden",
        r"financial burden",
    ],
    "CBA": [
        r"cost[\s-]*benefit",
        r"benefit[\s-]*cost ratio",
        r"net benefit",
    ],
    "BIA": [
        r"budget impact",
    ],
}

# Cost component keywords (searched in context around money mentions)
_COST_COMPONENTS = {
    "direct_medical": [
        r"direct\s+medical",
        r"hospitali[sz]ation",
        r"inpatient",
        r"outpatient",
        r"treatment cost",
        r"medical cost",
        r"clinical cost",
        r"hospital cost",
    ],
    "productivity_illness": [
        r"productivity loss",
        r"indirect cost",
        r"work\s*days?\s+lost",
        r"absenteeism",
        r"presenteeism",
    ],
    "productivity_death": [
        r"premature death",
        r"premature mortality",
        r"years? of life lost",
        r"YLL",
        r"mortality cost",
        r"human capital",
    ],
    "vaccine_cost": [
        r"vaccine\s+cost",
        r"vaccine\s+price",
        r"cost\s+per\s+dose",
        r"price\s+per\s+dose",
        r"immunization\s+cost",
    ],
    "delivery_cost": [
        r"delivery\s+cost",
        r"administration\s+cost",
        r"supply\s+chain",
        r"program\s*(?:me)?\s+cost",
        r"logistics\s+cost",
    ],
    "cost_per_case": [
        r"cost\s+per\s+case",
        r"cost\s+per\s+episode",
        r"per[\s-]+case\s+cost",
        r"per[\s-]+patient\s+cost",
    ],
    "cost_per_daly": [
        r"cost\s+per\s+DALY",
        r"\$/\s*DALY",
        r"per\s+DALY\s+averted",
    ],
}

# Common LMICs and regions (non-exhaustive, for rough extraction)
_COUNTRIES_RE = re.compile(
    r"\b(?:Bangladesh|India|Pakistan|Nigeria|Ethiopia|Kenya|Tanzania|"
    r"Uganda|Mozambique|DRC|Congo|Haiti|Somalia|Yemen|Sudan|"
    r"South\s+Africa|Ghana|Senegal|Mali|Niger|Chad|"
    r"Indonesia|Philippines|Vietnam|Thailand|Myanmar|Cambodia|Laos|"
    r"Brazil|Colombia|Peru|Mexico|Guatemala|Honduras|"
    r"China|Egypt|Iraq|Afghanistan|Nepal|Zimbabwe|Zambia|Malawi|"
    r"sub[\s-]*Saharan\s+Africa|South[\s-]*East\s+Asia|"
    r"Southeast\s+Asia|Latin\s+America|low[\s-]*income|"
    r"middle[\s-]*income|LMIC|LIC|United\s+States|US(?=\s|$)|UK|Europe|"
    r"global|worldwide)\b",
    re.IGNORECASE,
)


def classify_study_type(text: str) -> list[str]:
    """Classify study type(s) from abstract text."""
    found = []
    for stype, patterns in _STUDY_TYPES.items():
        for pat in patterns:
            if re.search(pat, text, re.IGNORECASE):
                found.append(stype)
                break
    return found if found else ["other"]


def extract_countries(text: str) -> list[str]:
    """Extract country/region mentions from text."""
    matches = _COUNTRIES_RE.findall(text)
    # Deduplicate preserving order
    seen = set()
    unique = []
    for m in matches:
        m_lower = m.lower()
        if m_lower not in seen:
            seen.add(m_lower)
            unique.append(m)
    return unique


def extract_money_mentions(text: str) -> list[dict]:
    """
    Extract monetary amounts with surrounding context.

    Returns list of dicts with keys: amount_raw, amount_usd, context, components.
    """
    mentions = []
    for match in _MONEY_RE.finditer(text):
        raw_num = match.group(1).replace(",", "")
        multiplier_str = (match.group(2) or "").lower()
        multipliers = {
            "": 1,
            "thousand": 1_000,
            "million": 1_000_000,
            "billion": 1_000_000_000,
            "trillion": 1_000_000_000_000,
        }
        try:
            amount = float(raw_num) * multipliers.get(multiplier_str, 1)
        except ValueError:
            amount = None

        # Get surrounding context (±150 chars)
        start = max(0, match.start() - 150)
        end = min(len(text), match.end() + 150)
        context = text[start:end]

        # Classify cost component from context
        components = []
        for comp, patterns in _COST_COMPONENTS.items():
            for pat in patterns:
                if re.search(pat, context, re.IGNORECASE):
                    components.append(comp)
                    break

        mentions.append(
            {
                "amount_raw": match.group(0).strip(),
                "amount_usd": amount,
                "context": context.strip(),
                "components": components,
            }
        )
    return mentions


def extract_from_article(
    article: dict, disease_key: str, fulltext: str | None = None,
) -> dict:
    """
    Extract structured cost information from a single article.

    Parameters
    ----------
    article : dict
        Article dict from PubMed cache (pmid, title, abstract, ...).
    disease_key : str
        The disease this article was found under.
    fulltext : str or None
        Full-text content extracted from PDF, if available.

    Returns
    -------
    dict
        Structured extraction with all identified cost data.
    """
    abstract_text = f"{article.get('title', '')} {article.get('abstract', '')}"
    # Use fulltext if available; extract from both but prefer fulltext for money
    has_fulltext = bool(fulltext and fulltext.strip())
    search_text = fulltext if has_fulltext else abstract_text

    study_types = classify_study_type(search_text)
    countries = extract_countries(search_text)
    money = extract_money_mentions(search_text)

    return {
        "disease_key": disease_key,
        "pmid": article.get("pmid", ""),
        "title": article.get("title", ""),
        "year": article.get("year", ""),
        "journal": article.get("journal", ""),
        "study_types": study_types,
        "countries": countries,
        "money_mentions": money,
        "n_money_mentions": len(money),
        "has_cost_per_case": any(
            "cost_per_case" in m["components"] for m in money
        ),
        "has_vaccine_cost": any(
            "vaccine_cost" in m["components"] or "delivery_cost" in m["components"]
            for m in money
        ),
        "source": "fulltext" if has_fulltext else "abstract",
        "abstract_snippet": article.get("abstract", "")[:500],
    }


def extract_all(
    cache_path: str = RAW_CACHE_PATH,
    output_path: str = EXTRACTION_PATH,
    fulltext_path: str | None = None,
    disease_keys: list[str] | None = None,
) -> list[dict]:
    """
    Run extraction on all cached articles and save JSON.

    Parameters
    ----------
    cache_path : str
        Path to PubMed raw results JSON cache.
    output_path : str
        Where to save extractions JSON.
    fulltext_path : str or None
        Path to fulltext JSON map ({pmid: text}). If provided, fulltext
        is used for extraction when available.
    disease_keys : list[str] or None
        Subset of diseases to extract. None = all in cache.

    Returns
    -------
    list[dict]
        All extractions.
    """
    with open(cache_path, "r", encoding="utf-8") as f:
        cache = json.load(f)

    # Load fulltext map if available
    fulltext_map = {}
    if fulltext_path:
        try:
            with open(fulltext_path, "r", encoding="utf-8") as f:
                fulltext_map = json.load(f)
            print(f"  Loaded fulltext for {len(fulltext_map)} articles.")
        except FileNotFoundError:
            print(f"  No fulltext file at {fulltext_path}, using abstracts only.")

    keys = disease_keys or list(cache.keys())
    extractions = []

    for key in keys:
        if key not in cache:
            print(f"  WARNING: '{key}' not in cache, skipping.")
            continue
        disease_data = cache[key]
        articles = disease_data.get("articles", [])
        print(f"  Extracting {key}: {len(articles)} articles ...", end=" ")

        n_ft = 0
        for article in articles:
            pmid = article.get("pmid", "")
            fulltext = fulltext_map.get(pmid)
            if fulltext:
                n_ft += 1
            ext = extract_from_article(article, key, fulltext=fulltext)
            if ext["n_money_mentions"] > 0 or ext["study_types"] != ["other"]:
                extractions.append(ext)

        n_kept = sum(1 for e in extractions if e["disease_key"] == key)
        ft_note = f" ({n_ft} with fulltext)" if n_ft else ""
        print(f"{n_kept} with cost data{ft_note}.")

    # Save
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(extractions, f, indent=2, ensure_ascii=False)
    print(f"Extractions saved to {output_path} ({len(extractions)} articles)")

    return extractions
