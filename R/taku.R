# Load the Taku analysis functions.
#
# This is the only R file that readers need to source directly.  The files in
# R/internal contain implementation details and are kept separate only so that
# model code, data generation, and estimators can be tested independently.
taku_source <- function(path) {
  source(file.path("R", "internal", path), local = FALSE)
}

for (file in c(
  "data.R",
  "scenarios.R",
  "models.R",
  "simulation_data.R",
  "frequentist.R",
  "bayesian.R",
  "workflow.R"
)) {
  taku_source(file)
}

rm(taku_source, file)
