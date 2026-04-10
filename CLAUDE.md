# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Mathematical framework and policy analysis for optimal OCV (oral cholera vaccine) allocation between pre-emptive deployment and reactive outbreak response, using ADM2-level cholera surveillance data from 2010–2012.

## Stack

- **R** — analytics, simulation, figures (Quarto chunks + `R/` scripts)
- **Python** — PubMed literature search pipeline (`scripts/pubmed_search/`)
- **Quarto** — documents rendered to HTML

## Key commands

> Rendered HTML files already exist alongside each `.qmd`. Read the `.qmd` source rather than re-rendering (rendering is slow, especially `vaccine_impact_3yr.qmd`).

```bash
# Render a Quarto document (only when source has changed)
# Windows path — Rscript is at "C:/Program Files/R/R-4.5.3/bin/Rscript.exe"
# Quarto is at "C:/Users/jonghoon.kim/AppData/Local/Programs/Quarto/bin/quarto.exe"
"/c/Users/jonghoon.kim/AppData/Local/Programs/Quarto/bin/quarto.exe" render vaccine_impact_3yr.qmd
"/c/Users/jonghoon.kim/AppData/Local/Programs/Quarto/bin/quarto.exe" render analytic_study.qmd

# Run standalone R scripts
"/c/Program Files/R/R-4.5.3/bin/Rscript.exe" --no-save --no-restore some_script.R

# PubMed pipeline (use uv, not pip)
uv pip install -r scripts/pubmed_search/requirements.txt
python -m scripts.pubmed_search.main --stage search
python -m scripts.pubmed_search.main --stage extract
python -m scripts.pubmed_search.main --stage output
python -m scripts.pubmed_search.main --diseases cholera dengue
```

> **Warning**: `Rscript.exe` segfaults when loading `outbreak_shapefiles.rds` and then calling most functions on the sf object inline. Use a separate script file (`explore_shp.R`) and call `saveRDS()` to extract data rather than chaining operations in `-e` expressions.

## Core analytical framework

**Central question**: Given vaccination capacity covering fraction **f** of n subunits, what fraction **α** to allocate pre-emptively vs. hold for reactive use?

**Key parameters**:
- **R = C_I / C_V** — cost ratio (illness cost / vaccine cost per dose)
- **p** — outbreak probability per subunit (in `vaccine_impact_3yr.qmd`: `p_3yr` = fraction of all shapefile ADM2 with ≥1 outbreak in WIN_YEARS)
- **p_eff** — expected reactive events per site over shapefile universe (= `sum(n_events) / N_SHP`)
- **r** — reactive VE (`R_REACT = 0.40`)
- **α** — pre-emptive fraction (α=0 pure reactive, α=1 pure pre-emptive)
- **ρ** — Spearman rank correlation between true risk and targeting score
- **f** — vaccination capacity (`F_VAX = 0.60`)
- **MIN_REACT_GAP** — days of reactive protection per campaign (default 180); outbreaks within this window are merged into one event
- **MIN_PREEMPTIVE_GAP** — days of pre-emptive protection (default `T_WIN × 365` ≈ 3 years)

**Three-regime α* formula** (homogeneous case, from `alpha_star_revised`):
- `R < R_lower = (1.5 - 0.5·p_eff) / (q̄·(1-r))` → α* = 0 (pure reactive)
- `R_lower ≤ R ≤ R_upper` → α* = α_c = (f − p_eff) / (f·(1 − p_eff)) (interior mixed)
- `R > R_upper = p_eff / (q̄·(p_eff - r))` → α* = 1 (pure pre-emptive)

## Key source documents

| File | Content |
|------|---------|
| `vaccine_impact_3yr.qmd` | **Primary analysis** — empirical cost model using 2010–2012 ADM2 cholera data |
| `analytic_study.qmd` | Single/multi-population cost models, p_crit derivations |
| `analytic_hetero.qmd` | Heterogeneous risk (Beta distribution), imperfect targeting |
| `simulation.qmd` | Endogenous reactive effectiveness from epidemic trajectory |
| `policy_report.qmd` | Policy-facing report for n=50, f=0.5 scenario |
| `cost_ratio_literature.qmd` | PubMed-derived R estimates across 14 VIMC diseases |
| `R/analytic_utils.R` | Closed-form single/multi-population cost and threshold functions |
| `policy_simulations.R` | 5 simulation experiments for n=50 scenario |

## Data architecture in `vaccine_impact_3yr.qmd`

The document is self-contained; all data flows from two source RDS files:

1. **`data/outbreak_files_Oct2024/time_series_outbreak_extraction.rds`** — time-series outbreak records; filtered to ADM2 (`str_count(location, "::") - 1 == 2`) → `ts2`
2. **`data/outbreak_files_Oct2024/outbreak_shapefiles.rds`** — sf object (4381 rows) providing the full geographic universe; `lctn_pr` column (character, 4302 numeric + 79 "composite_loc_*") maps 1:1 to `ts2$location_period_id`; used only for `N_SHP` (row count as denominator for p_3yr and p_eff)

Key derived objects built in chunk order:
- `sim_df` — one row per ts2 ADM2 site: `pop_sum`, `cases_win`, `any_outbreak`, `q_i`, `n_events`, `cases_pre`, `q_pre_i`
- `event_df` — one row per reactive event (post `merge_gap_events` gap filter): `location`, `location_period_id`, `outbreak_number`, `date_event`, `cases`, `pop_sum`, `site_idx`
- `N_SHP` — count of all shapefile rows (4381), used as denominator for `p_3yr` and `p_eff`

The `site_idx` column in `event_df` is a row index into `sim_df` — these must stay aligned (never reorder `sim_df` after `event_df` is built).

## Helper functions in `vaccine_impact_3yr.qmd`

**`merge_gap_events(dates, cases, gap_days)`** — greedy within-location merger: outbreaks within `gap_days` of the last vaccination are absorbed (cases added to trigger); returns `list(keep, merged_cases)`.

**`cost_incidence(alpha, f, q, pop, score, evt, R, r, nu, q_star, q_pre)`** — total per-capita cost given an α; uses `DOSE_PRE`, `DOSE_REACT`, `DOSE_WASTE` globals; `q_pre` defaults to `q` (cases averted by pre-emptive limited to protection window).

**`cases_averted_breakdown(...)`** — same allocation logic, returns list of case counts instead of cost.

Both functions allocate pre-emptive sites by descending `score`, then reactive events chronologically by `date_event` until budget exhausted.

## PubMed pipeline

Two-stage JSON cache: search cache (`pubmed_raw_results.json`) → extraction cache (`pubmed_cost_extractions.json`) → QMD output. Sci-Hub client requires DOI (not PMID); working mirror as of March 2026: `https://sci-hub.kr`.

## Diseases (VIMC list)

cholera, covid, typhoid, dengue, hib, rotavirus, pneumococcal, hpv, japanese_encephalitis, malaria, measles, meningitis, rubella, yellow_fever

## Manuscript framing

1. **Stockpile allocation** — optimal pre-emptive/reactive split for ICG OCV stockpile; target *Nature Medicine* or *Lancet Global Health*
2. **Reactive effectiveness gap** — when reactive fails despite high VE due to delays; target *Lancet Infectious Diseases* or *Science*
3. **Value of risk intelligence** — ρ-threshold for surveillance investment; target *Nature Medicine* or *Lancet Global Health*
