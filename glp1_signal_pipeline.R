# =====================================================================
# GLP-1 RA Safety Signal Detection Pipeline (openFDA / FAERS)
# ---------------------------------------------------------------------
# Question: Do GLP-1 receptor agonists show disproportionate reporting of
# selected GI/hepatobiliary adverse events in FDA spontaneous reports,
# vs. all other drugs and vs. an active comparator (sitagliptin, DPP-4)?
#
# Run it:  RStudio -> open this file -> click "Source"
#          or terminal -> Rscript glp1_signal_pipeline.R
# Outputs are written to an "output/" folder next to this script.
#
# Key limitations (know these for interviews):
#  - Spontaneous reports: no denominator (can't estimate incidence/risk)
#  - Reporting bias & stimulated reporting (media attention, litigation)
#  - Duplicate reports; drug may be listed but not the suspect drug
#  - Disproportionality = hypothesis-generating signal, not causation
# =====================================================================


# ---- 0. Setup --------------------------------------------------------
pkgs <- c("httr2", "dplyr", "purrr", "tidyr", "ggplot2", "readr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

dir.create("output", showWarnings = FALSE)


# ---- 1. Study parameters (edit these) --------------------------------
drugs      <- c("semaglutide", "tirzepatide", "liraglutide", "sitagliptin")
comparator <- "sitagliptin"   # active comparator (DPP-4 inhibitor)

# MedDRA Preferred Terms, uppercase, exactly as they appear in FAERS.
# Verify spellings with top_terms() below before trusting results.
events <- c("PANCREATITIS",
            "CHOLELITHIASIS",
            "CHOLECYSTITIS",
            "INTESTINAL OBSTRUCTION",
            "IMPAIRED GASTRIC EMPTYING")

trend_years <- 2018:2025

# Optional: a free openFDA API key raises your daily request limit.
# Put OPENFDA_API_KEY=yourkey in ~/.Renviron, then restart R.
api_key <- Sys.getenv("OPENFDA_API_KEY")

base_url <- "https://api.fda.gov/drug/event.json"


# ---- 2. API helper functions -----------------------------------------

# Build search strings
drug_q  <- function(drug)  sprintf('patient.drug.openfda.generic_name:"%s"', drug)
event_q <- function(event) sprintf('patient.reaction.reactionmeddrapt.exact:"%s"', event)
date_q  <- function(start, end) sprintf("receivedate:[%s TO %s]",
                                        format(start, "%Y%m%d"), format(end, "%Y%m%d"))

# Return the number of reports matching a search (NULL = all reports).
# openFDA returns 404 when nothing matches, so we treat that as 0.
count_reports <- function(search = NULL) {
  req <- request(base_url) |>
    req_url_query(limit = 1) |>
    req_error(is_error = function(resp) FALSE) |>
    req_retry(max_tries = 3)
  if (!is.null(search)) req <- req_url_query(req, search = search)
  if (nzchar(api_key))  req <- req_url_query(req, api_key = api_key)

  resp <- req_perform(req)
  Sys.sleep(0.3)  # be polite to the API

  status <- resp_status(resp)
  if (status == 404) return(0)
  if (status != 200) stop("openFDA returned ", status, " for search: ", search)
  resp_body_json(resp)$meta$results$total
}

# Look up the most common reaction terms for a drug (use to check spellings)
# e.g. grep("GASTRIC", top_terms("semaglutide")$term, value = TRUE)
top_terms <- function(drug, n = 1000) {
  resp <- request(base_url) |>
    req_url_query(search = drug_q(drug),
                  count  = "patient.reaction.reactionmeddrapt.exact",
                  limit  = n) |>
    req_perform()
  resp_body_json(resp, simplifyVector = TRUE)$results
}


# ---- 3. EXTRACT: pull the 2x2 table counts ---------------------------
#                 event   no event
#   drug            a        b
#   other drugs     c        d
message("Pulling total report count...")
N <- count_reports()

message("Pulling drug totals...")
drug_totals <- tibble(drug = drugs) |>
  mutate(drug_total = map_dbl(drug, \(d) count_reports(drug_q(d))))

message("Pulling event totals...")
event_totals <- tibble(event = events) |>
  mutate(event_total = map_dbl(event, \(e) count_reports(event_q(e))))

message("Pulling drug-event pair counts...")
counts <- expand_grid(drug = drugs, event = events) |>
  mutate(a = map2_dbl(drug, event,
                      \(d, e) count_reports(paste(drug_q(d), "AND", event_q(e))))) |>
  left_join(drug_totals,  by = "drug") |>
  left_join(event_totals, by = "event") |>
  mutate(N = N,
         b = drug_total  - a,
         c = event_total - a,
         d = N - a - b - c)

write_csv(counts, "output/raw_counts.csv")


# ---- 4. TRANSFORM: disproportionality vs. all other drugs ------------
z <- qnorm(0.975)

results <- counts |>
  mutate(
    ror       = ifelse(a > 0, (a / b) / (c / d), NA),
    ror_se    = sqrt(1/a + 1/b + 1/c + 1/d),
    ror_lower = exp(log(ror) - z * ror_se),
    ror_upper = exp(log(ror) + z * ror_se),
    prr       = ifelse(a > 0, (a / (a + b)) / (c / (c + d)), NA),
    prr_se    = sqrt(1/a - 1/(a + b) + 1/c - 1/(c + d)),
    prr_lower = exp(log(prr) - z * prr_se),
    prr_upper = exp(log(prr) + z * prr_se),
    # Simple signal rule: at least 3 reports and lower 95% CI of ROR > 1
    signal    = a >= 3 & !is.na(ror_lower) & ror_lower > 1
  )

write_csv(results, "output/disproportionality_all_drugs.csv")


# ---- 5. TRANSFORM: ROR vs. active comparator -------------------------
# Restricts the comparison to GLP-1 RA vs. sitagliptin reports, which
# partly addresses confounding by indication (both used for diabetes).
comp <- counts |>
  filter(drug == comparator) |>
  select(event, a_comp = a, total_comp = drug_total)

active <- counts |>
  filter(drug != comparator) |>
  left_join(comp, by = "event") |>
  mutate(
    no_event      = drug_total - a,
    no_event_comp = total_comp - a_comp,
    ror_ac        = ifelse(a > 0 & a_comp > 0,
                           (a / no_event) / (a_comp / no_event_comp), NA),
    se_ac         = sqrt(1/a + 1/no_event + 1/a_comp + 1/no_event_comp),
    ror_ac_lower  = exp(log(ror_ac) - z * se_ac),
    ror_ac_upper  = exp(log(ror_ac) + z * se_ac)
  )

write_csv(active, "output/disproportionality_vs_comparator.csv")



# ---- 6. EXTRACT + TRANSFORM: quarterly reporting trends --------------
message("Pulling quarterly trends (takes a minute)...")
quarters <- expand_grid(year = trend_years, q = 1:4) |>
  mutate(
    start = as.Date(sprintf("%d-%02d-01", year, (q - 1) * 3 + 1)),
    end   = as.Date(ifelse(q == 4,
                           sprintf("%d-01-01", year + 1),
                           sprintf("%d-%02d-01", year, q * 3 + 1))) - 1
  )

trends <- expand_grid(drug = drugs, quarters) |>
  mutate(reports = pmap_dbl(list(drug, start, end),
                            \(d, s, e) count_reports(paste(drug_q(d), "AND", date_q(s, e)))))

write_csv(trends, "output/quarterly_trends.csv")


# ---- 7. OUTPUT: plots ------------------------------------------------
p_forest <- results |>
  filter(!is.na(ror)) |>
  ggplot(aes(x = event, y = ror, ymin = ror_lower, ymax = ror_upper, color = drug)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = 0.6)) +
  scale_y_log10() +
  coord_flip() +
  labs(title = "Reporting odds ratios vs. all other drugs (FAERS)",
       x = NULL, y = "ROR (95% CI, log scale)", color = "Drug") +
  theme_minimal()
ggsave("output/forest_all_drugs.png", p_forest, width = 8, height = 5, bg = "white")

p_active <- active |>
  filter(!is.na(ror_ac)) |>
  ggplot(aes(x = event, y = ror_ac, ymin = ror_ac_lower, ymax = ror_ac_upper, color = drug)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = 0.6)) +
  scale_y_log10() +
  coord_flip() +
  labs(title = paste("Reporting odds ratios vs. active comparator:", comparator),
       x = NULL, y = "ROR (95% CI, log scale)", color = "Drug") +
  theme_minimal()
ggsave("output/forest_vs_comparator.png", p_active, width = 8, height = 5, bg = "white")

p_trend <- ggplot(trends, aes(x = start, y = reports, color = drug)) +
  geom_line(linewidth = 0.9) +
  labs(title = "FAERS reports received per quarter",
       x = NULL, y = "Reports", color = "Drug") +
  theme_minimal()
ggsave("output/quarterly_trends.png", p_trend, width = 8, height = 5, bg = "white")


# ---- 8. Summary ------------------------------------------------------
message("\nDone. Signals (>=3 reports, ROR lower CI > 1, vs. all drugs):")
results |>
  filter(signal) |>
  select(drug, event, a, ror, ror_lower, ror_upper) |>
  mutate(across(c(ror, ror_lower, ror_upper), \(x) round(x, 2))) |>
  print(n = Inf)
message("Files saved to: ", normalizePath("output"))

