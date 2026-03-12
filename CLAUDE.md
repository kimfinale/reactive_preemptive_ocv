# Project: Reactive vs Pre-emptive Vaccination

Mathematical framework and policy analysis for optimal vaccine allocation between
pre-emptive deployment and reactive outbreak response.

## Stack

- **R** — analytics, simulation, figures (Quarto chunks + `R/` scripts)
- **Python** — PubMed literature search pipeline (`scripts/pubmed_search/`)
- **Quarto** — documents rendered to HTML (`quarto render <file>.qmd`)

## Key commands

> These are reference commands — nothing here runs automatically.
> Rendered HTML files already exist alongside each .qmd; read the .qmd source
> or the .html for context rather than re-rendering (rendering is slow).

```bash
# Render a Quarto document (only when source has changed)
quarto render analytic_study.qmd
quarto render policy_report.qmd
quarto render cost_ratio_literature.qmd

# Run R simulations standalone
Rscript policy_simulations.R

# PubMed literature search pipeline (use uv, not pip)
uv pip install -r scripts/pubmed_search/requirements.txt
python -m scripts.pubmed_search.main                          # abstract-only (all 14 diseases)
python -m scripts.pubmed_search.main --stage search           # PubMed search only
python -m scripts.pubmed_search.main --stage fulltext --max-pdfs 30 --delay 5  # Sci-Hub PDFs
python -m scripts.pubmed_search.main --stage fulltext_extract # extract text from PDFs
python -m scripts.pubmed_search.main --stage extract          # cost extraction
python -m scripts.pubmed_search.main --stage output           # build cost_ratio_literature.qmd
python -m scripts.pubmed_search.main --stage all_with_fulltext --max-pdfs 30
python -m scripts.pubmed_search.main --diseases cholera dengue  # subset
```

## Core analytical framework

The central question: given vaccination capacity covering fraction **f** of
n subunits, what fraction **α** to allocate pre-emptively vs. hold for reactive use?

Key parameters:
- **R = C_I / C_V** — cost ratio (infection cost per subunit / vaccination cost per subunit)
- **p** — outbreak probability per subunit
- **r** — reactive effectiveness (fraction of cases averted by reactive campaign)
- **α** — pre-emptive fraction (α=0 pure reactive, α=1 pure pre-emptive)
- **ρ** — Spearman rank correlation between true risk and targeting score (targeting accuracy)
- **f** — vaccination capacity (fraction of subunits coverable)
- **n** — number of subunits (policy scenario: n=50, f=0.5)

Critical threshold: **p_crit** — above which pre-emptive is preferred.
- Single subunit: `p_crit^(1) = r / (r + R(ν - r))` where ν = 1 (full pre-emptive coverage)
- Multi-subunit (abundant capacity): `p_crit^(n) = f / (1 + R(fν - r))`

## Key source documents

| File | Content |
|------|---------|
| `analytic_study.qmd` | Single/multi-population cost models, p_crit derivations |
| `analytic_hetero.qmd` | Heterogeneous risk (Beta distribution), imperfect targeting |
| `simulation.qmd` | Endogenous reactive effectiveness from epidemic trajectory |
| `policy_report.qmd` | Policy-facing report for n=50, f=0.5 scenario |
| `cost_ratio_literature.qmd` | PubMed-derived R estimates across 14 VIMC diseases |
| `R/analytic_utils.R` | All core R functions |
| `policy_simulations.R` | 5 simulation experiments for n=50 scenario |

## PubMed pipeline data files

| File | Content |
|------|---------|
| `data/pubmed_raw_results.json` | Cached PubMed search results (14 diseases, ~2600 articles) |
| `data/pubmed_cost_extractions.json` | Structured cost extractions (~2019 articles with cost data) |
| `data/pubmed_fulltext_map.json` | Extracted fulltext keyed by PMID (where PDFs were downloaded) |
| `data/fulltext_pdfs/{disease}/` | Downloaded PDFs via Sci-Hub |

## PubMed search pipeline — notes

- Sci-Hub client is vendored in `scripts/pubmed_search/scihub_client.py`
- Requires DOI (not PMID) — sci-hub.kr does not support PMID lookup
- DOIs are now extracted during efetch (added to `pubmed_api.py`)
- Working mirror as of March 2026: `https://sci-hub.kr` (other mirrors may be blocked)
- Two-stage JSON cache: search cache → extraction cache → QMD output

## Diseases (VIMC list)

cholera, covid, typhoid, dengue, hib, rotavirus, pneumococcal, hpv,
japanese_encephalitis, malaria, measles, meningitis, rubella, yellow_fever

## Manuscript framing (three options identified)

1. **Stockpile allocation problem** — optimal pre-emptive/reactive split for ICG OCV
   stockpile; target *Nature Medicine* or *Lancet Global Health*
2. **Reactive effectiveness gap** — when reactive vaccination fails despite high VE
   due to delays + epidemic speed; target *Lancet Infectious Diseases* or *Science*
3. **Value of risk intelligence** — ρ-threshold for surveillance investment;
   target *Nature Medicine* or *Lancet Global Health*

## R packages used

`tidyverse`, `ggplot2`, `patchwork`, `knitr`, `copula` (for Gaussian copula targeting)

## Simulation outputs

Located in project root as `.rds` files:
`policy_sim1_optimal_alpha.rds` through `policy_sim5_robustness.rds`
