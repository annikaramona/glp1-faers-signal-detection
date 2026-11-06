# GLP-1 Receptor Agonist Safety Signal Detection (FAERS / openFDA)

A reproducible R pipeline that queries the FDA Adverse Event Reporting System (FAERS) through the openFDA API and runs a **disproportionality analysis**. The goal is to identify possible safety signals for GLP-1 receptor agonists (GLP-1 RAs) and selected gastrointestinal and hepatobiliary adverse events.

## Research question

Are GLP-1 RAs (semaglutide, tirzepatide, liraglutide) reported with pancreatitis, gallbladder disease, intestinal obstruction, or impaired gastric emptying more often than expected? The analysis compares them against:

1. All other drugs in FAERS
2. An active comparator, **sitagliptin** (a DPP-4 inhibitor), which is also used to treat type 2 diabetes

## Data source

- **FAERS**, accessed via the [openFDA drug adverse event API](https://open.fda.gov/apis/drug/event/)
- Drugs are matched on the standardized generic name field (`patient.drug.openfda.generic_name`)
- Events are matched on MedDRA Preferred Terms (`patient.reaction.reactionmeddrapt`)
- The pipeline pulls aggregate **report counts** only and never downloads individual case records

## Methods

### Pipeline structure

| Step | What it does |
|---|---|
| **Parameters** | Drugs, events, comparator, and years are defined in one section, so the analysis can be repointed to a different drug class by editing a few lines |
| **Extract** | A reusable function, `count_reports()`, queries the API and returns a report count. It is looped over every drug-event pair (`purrr`) to build 2×2 tables, with retry logic and rate limiting |
| **Transform** | Computes disproportionality metrics against all other drugs and against the active comparator, plus quarterly reporting trends |
| **Output** | Writes CSV tables and plots to `output/` |

### Disproportionality metrics

For each drug-event pair, a 2×2 table of **reports**:

|  | Event | Other events |
|---|---|---|
| **Drug of interest** | a | b |
| **All other drugs** | c | d |

- **Reporting odds ratio (ROR)** = (a/b) / (c/d), with a 95% CI on the log scale
- **Proportional reporting ratio (PRR)** = [a/(a+b)] / [c/(c+d)], with a 95% CI
- **Signal criterion:** at least 3 reports and a lower 95% CI bound of the ROR above 1

### Active comparator

RORs are also computed by restricting the comparison to GLP-1 RA reports vs. sitagliptin reports. Both drugs are used for type 2 diabetes, so this partly addresses **confounding by indication**: events that track with diabetes itself, rather than with the drug, are less likely to show up as false signals.

### Reporting trends

Report counts per quarter (2018–2025) are pulled for each drug to describe how reporting volume changed as GLP-1 RA use expanded.

## Outputs

| File | Description |
|---|---|
| `output/raw_counts.csv` | 2×2 counts for every drug-event pair |
| `output/disproportionality_all_drugs.csv` | ROR, PRR, 95% CIs, and signal flag vs. all other drugs |
| `output/disproportionality_vs_comparator.csv` | ROR and 95% CI vs. sitagliptin |
| `output/quarterly_trends.csv` | Reports received per quarter, by drug |
| `output/forest_all_drugs.png` | Forest plot of RORs vs. all other drugs |
| `output/forest_vs_comparator.png` | Forest plot of RORs vs. sitagliptin |
| `output/quarterly_trends.png` | Quarterly reporting trend plot |

## Results

<!-- Fill in after running the pipeline. Suggested structure:
- Which drug-event pairs met the signal criterion vs. all other drugs
- Whether those signals held up vs. the active comparator (sitagliptin)
- What the quarterly trend shows about reporting volume over time
- Embed plots, e.g. ![Forest plot](output/forest_all_drugs.png)
-->

*Results to be added.*

## Limitations

- **No denominator.** FAERS contains reports, not exposed patients, so results cannot estimate incidence or absolute or relative risk.
- **Reporting bias.** Reporting is voluntary and incomplete. Media attention, regulatory communications, and litigation can produce **stimulated reporting**, which is especially relevant for GLP-1 RAs given their recent high profile.
- **Duplicate reports.** The same case can be submitted more than once (by patients, clinicians, and manufacturers).
- **Drug role.** Counts include reports where the drug is listed, not only those where it is the primary suspect.
- **Confounding.** The active comparator reduces, but does not remove, confounding by indication. Differences in patient populations (e.g., weight-loss vs. diabetes use) remain.
- **Hypothesis-generating only.** Disproportionality identifies signals for further study and does not establish causation.

## Next steps

- Restrict to reports where the drug is the primary suspect
- Stratify by indication, age, and sex
- Test signals in real-world data (claims or EHR) using an **active-comparator, new-user design** with **propensity score** methods
- Add an LLM step that drafts a plain-language signal summary from the results table, with human review

## How to run

**Requirements:** R ≥ 4.4 and the packages `httr2`, `dplyr`, `purrr`, `tidyr`, `ggplot2`, and `readr`. The script installs any that are missing.

```bash
Rscript glp1_signal_pipeline.R
```

Or open `glp1_signal_pipeline.R` in RStudio, set the working directory to the script's location, and click **Source**. A full run takes a few minutes because of the API calls.

**Optional:** a free [openFDA API key](https://open.fda.gov/apis/authentication/) raises the daily request limit. Add `OPENFDA_API_KEY=yourkey` to `~/.Renviron` and restart R.

**Customizing:** edit the `drugs`, `comparator`, `events`, and `trend_years` objects in section 1. Event names must match MedDRA Preferred Terms as they appear in FAERS; use `top_terms("drugname")` to check spellings.

## Tools

R · openFDA REST API · httr2 · tidyverse (dplyr, purrr, tidyr, ggplot2, readr)
