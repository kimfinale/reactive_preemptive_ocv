---
name: pubmed-search
description: Search PubMed for cost-of-illness and vaccination cost literature to estimate the cost ratio R = C_I/C_V across vaccine-preventable diseases.
user_invocable: true
---

# PubMed Cost-Ratio Literature Search

Search PubMed for economic/cost studies on vaccine-preventable diseases and extract cost data to estimate plausible ranges of the cost ratio R = C_I / C_V.

## Diseases covered

14 VIMC diseases: cholera, COVID-19, typhoid, dengue, Hib, rotavirus, pneumococcal (PCV), HPV, Japanese encephalitis, malaria, measles, meningitis, rubella, yellow fever.

## Pipeline

The tool runs a multi-stage pipeline:

1. **Search** (`--stage search`) — Query PubMed via NCBI E-utilities. Cache raw results as JSON.
2. **Fulltext download** (`--stage fulltext`) — Download PDFs via Sci-Hub (vendored client with mirror fallback). Stores in `data/fulltext_pdfs/{disease}/`.
3. **Fulltext extract** (`--stage fulltext_extract`) — Extract text from PDFs using PyMuPDF. Saves `{pmid: text}` JSON map.
4. **Cost extraction** (`--stage extract`) — Regex-based extraction of USD amounts, study types (CEA/COI/CBA/BIA), countries, cost components. Uses fulltext when available, falls back to abstract.
5. **Output** (`--stage output`) — Generate Quarto document with per-disease summary tables.

## Usage

```bash
# Install dependencies
pip install -r scripts/pubmed_search/requirements.txt

# Abstract-only pipeline (fast, no Sci-Hub)
python -m scripts.pubmed_search.main

# Full pipeline including Sci-Hub fulltext
python -m scripts.pubmed_search.main --stage all_with_fulltext

# Individual stages
python -m scripts.pubmed_search.main --stage search
python -m scripts.pubmed_search.main --stage fulltext --max-pdfs 10 --delay 5
python -m scripts.pubmed_search.main --stage fulltext_extract
python -m scripts.pubmed_search.main --stage extract
python -m scripts.pubmed_search.main --stage output

# Specific diseases
python -m scripts.pubmed_search.main --diseases cholera dengue typhoid
```

## Outputs

- `data/pubmed_raw_results.json` — Cached PubMed search results
- `data/fulltext_pdfs/{disease}/` — Downloaded PDFs (when using fulltext stage)
- `data/pubmed_fulltext_map.json` — Extracted fulltext keyed by PMID
- `data/pubmed_cost_extractions.json` — Structured cost extractions per article
- `cost_ratio_literature.qmd` — Quarto narrative document with summary tables

## Two-stage caching

JSON caches between stages allow re-running downstream stages without re-querying PubMed or re-downloading PDFs. For example, to improve regex patterns: run `--stage extract` then `--stage output`.

## Limitations

- Sci-Hub access may be blocked by CAPTCHAs or unavailable mirrors
- Regex extraction is approximate — designed for human/AI review
- No currency conversion or inflation adjustment
- Non-English articles are excluded
