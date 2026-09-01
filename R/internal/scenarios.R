# Define the complete set of Bayesian scenarios used in the simulation study.
#
# Each row combines one observation model with one prior for the latent weekly
# stock composition. Keeping this design in one table prevents scenario names,
# model choices, and output-folder mappings from being repeated across scripts.
#
# Columns:
#   scenario_id: Human-readable identifier accepted by the command-line runner.
#   observation_model: "rdm" or "mmd"; interpreted by analysis_model_code().
#   prior: "ar1" or "dirichlet"; interpreted by prior_model_code().
# Returns:
#   A four-row data frame, with one row per supported scenario.
taku_scenarios <- function() {
  data.frame(
    scenario_id = c(
      "rdm_ar",
      "rdm_dirichlet",
      "mmd_ar",
      "mmd_dirichlet"
    ),
    observation_model = c("rdm", "rdm", "mmd", "mmd"),
    prior = c("ar1", "dirichlet", "ar1", "dirichlet"),
    display_name = c(
      "RDM + AR(1)",
      "RDM + Dirichlet",
      "MMD + AR(1)",
      "MMD + Dirichlet"
    ),
    stringsAsFactors = FALSE
  )
}

# Define the fixed category order used by the Taku inputs.
#
# The numeric order is not alphabetical; it follows the labels in the regional
# composition data. Keeping the mapping in code avoids scattering unexplained
# column numbers across the data generator and abundance calculations.
#
# Returns:
#   A four-row data frame containing the one-based category index, reporting
#   unit name, and lake/river type.
taku_reporting_units <- function() {
  data.frame(
    category = 1:4,
    reporting_unit = c(
      "Taku River",
      "Other Lakes",
      "Tatsamenie Lake",
      "U.S. River"
    ),
    type = c("river", "lake", "lake", "river"),
    stringsAsFactors = FALSE
  )
}

# Return the category indices that make up the lake-type abundance denominator.
taku_lake_categories <- function() {
  taku_reporting_units()$category[taku_reporting_units()$type == "lake"]
}

# Look up and validate a single scenario identifier.
#
# Arguments:
#   scenario_id: One value from taku_scenarios()$scenario_id.
#
# Returns:
#   A one-row data frame. Keeping the result as a data frame makes its named
#   fields convenient to pass to the model-building and pipeline functions.
get_taku_scenario <- function(scenario_id) {
  scenarios <- taku_scenarios()
  match_index <- match(scenario_id, scenarios$scenario_id)

  if (is.na(match_index)) {
    stop(
      "unknown scenario_id: ", scenario_id,
      "; choose one of: ", paste(scenarios$scenario_id, collapse = ", ")
    )
  }

  scenarios[match_index, , drop = FALSE]
}
