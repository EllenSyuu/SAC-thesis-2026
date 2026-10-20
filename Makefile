# Reproduce every table and figure:   make all        (~20 min on 2 cores)
# Fast smoke run (small Monte-Carlo):  make quick
# Package unit tests:                  make test

R      ?= Rscript
STEPS  := 01_prepare_data 02_marginal_gev 03_dependence_eda \
          04_simulation_study 05_archimax_pairwise 06_archimax_6d 07_key_findings

.PHONY: all quick test lint docs sensitivity clean $(STEPS)

all: $(STEPS) sensitivity

$(STEPS):
	$(R) analysis/$@.R

sensitivity:            ## Gumbel-standardised pseudo-observations (robustness)
	POBS=gumbel $(R) analysis/05_archimax_pairwise.R
	POBS=gumbel $(R) analysis/06_archimax_6d.R

quick:
	QUICK=1 $(MAKE) all

test:
	$(R) -e 'devtools::test()'

lint:
	$(R) -e 'lintr::lint_dir("R"); lintr::lint_dir("analysis")'

docs:
	$(R) -e 'roxygen2::roxygenise()'

clean:
	rm -rf data/processed/* results/figures/* results/tables/*
