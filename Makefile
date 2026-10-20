# Reproduce every table and figure:   make all      (~40 min cold, ~10 min with cache)
# Fast smoke run (small Monte-Carlo):  make quick
# Package unit tests / checks:         make test | make check | make lint
# Recompute cached results:            make refresh

R      ?= Rscript
STEPS  := 01_prepare_data 02_marginal_gev 03_dependence_eda \
          04_simulation_study 05_archimax_pairwise 06_archimax_6d \
          07_key_findings 08_sensitivity_min_days 09_bootstrap_ci

.PHONY: all quick test check lint docs sensitivity refresh clean start metadata $(STEPS)

all: start $(STEPS) sensitivity metadata

start:
	@rm -f results/run_steps.csv

$(STEPS):
	$(R) analysis/$@.R

sensitivity:            ## Gumbel-standardised pseudo-observations (robustness)
	POBS=gumbel $(R) analysis/05_archimax_pairwise.R
	POBS=gumbel $(R) analysis/06_archimax_6d.R

metadata:
	$(R) analysis/99_run_metadata.R

quick:
	QUICK=1 $(MAKE) all

refresh:
	FORCE=1 $(MAKE) all

test:
	$(R) -e 'devtools::test(stop_on_failure = TRUE)'

check:
	$(R) -e 'devtools::check(document = FALSE, manual = FALSE, error_on = "warning")'

lint:
	$(R) -e 'l <- c(lintr::lint_package(), lintr::lint_dir("analysis")); print(l); quit(status = length(l) > 0)'

docs:
	$(R) -e 'roxygen2::roxygenise()'

clean:
	rm -rf data/processed/* results/figures/* results/tables/* results/cache/*
