# Analysis of flow scenarios for three target species and carp at target
#    waterways to demonstrate the influence of species interactions
#    on population outcomes
#
# Target species: 
#     Murray Cod (Maccullochella peelii)
#     River Blackfish (Gadopsis marmoratus)
#     Murray-Darling Rainbowfish (Melanotaenia fluviatilis)
#
# Target waterways:
#     TBD
#
# Flow scenarios consider observed daily flows with and without 
#    e-water allocations and near-term (to 2025) forecasts of flows under
#    different plausible climates and flow management strategies

# Author: Jian Yen (jian.yen [at] deeca.vic.gov.au)
# Last updated: 8 May 2024

# load some packages
library(qs)
library(dplyr)
library(tidyr)
library(lubridate)
library(aae.db)
library(aae.hydro)
library(aae.pop.templates)
library(sf)
library(ggplot2)
library(ragg)
library(rstanarm)
library(bayesplot)

# and some load helpers
source("R/utils.R")
source("R/fish.R")
source("R/flow.R")
source("R/validation.R")

# settings
set.seed(2024-05-01)
nsim <- 1000
nburnin <- 10
simulate_again <- TRUE

# load data
cpue <- fetch_fish(recompile = FALSE)
cpue_exotic <- cpue |> 
  filter(grepl("perca|gambusia", scientific_name, ignore.case = TRUE))
cpue_exotic <- cpue_exotic |>
  mutate(waterbody = paste(
    tolower(gsub(" ", "_", waterbody)), reach_no, sep = "_r")
  ) |>
  group_by(scientific_name, waterbody, survey_year) |>
  summarise(presence = ifelse(sum(cpue) > 0, 1, 0))
cpue <- cpue |> 
  filter(!grepl("perca|gambusia", scientific_name, ignore.case = TRUE))

# download CPUE of recruits and adults
cpue_recruits <- fetch_fish(recruit = TRUE, recompile = FALSE)
cpue_adults <- fetch_fish(recruit = FALSE, adult = TRUE, recompile = FALSE)

# include Broken Creek Reach 5 in analyses (reach 5 runs to Rices Weir)
#   as well as MacKenzie River Reach 2
cpue <- cpue |>
  mutate(
    reach_no = ifelse(waterbody == "Broken Creek" & reach_no == 5, 4, reach_no),
    reach_no = ifelse(waterbody == "MacKenzie River" & reach_no == 2, 3, reach_no)
  )
cpue_recruits <- cpue_recruits |>
  mutate(
    reach_no = ifelse(waterbody == "Broken Creek" & reach_no == 5, 4, reach_no),
    reach_no = ifelse(waterbody == "MacKenzie River" & reach_no == 2, 3, reach_no)
  )
cpue_adults <- cpue_adults |>
  mutate(
    reach_no = ifelse(waterbody == "Broken Creek" & reach_no == 5, 4, reach_no),
    reach_no = ifelse(waterbody == "MacKenzie River" & reach_no == 2, 3, reach_no)
  )

# load stocking info
stocking <- read.csv("data/stocking.csv")

# fetch hydrological data
flow <- fetch_flow(start = 2009, end = 2024, recompile = FALSE)

# specify futures (individual years and events; can chop and change these to 
#    create specific scenarios)
flow_futures <- specify_flow_future(flow)

# and add some extra information on hypoxia risk and K by species (expands
#   over species)
flow_futures <- flow_futures |>
  left_join(
    fetch_hypoxia_risk(),
    by = "waterbody"
  ) |>
  left_join(
    fetch_carrying_capacity(), 
    by = "waterbody",
    relationship = "many-to-many"
  )

# calculate flow metrics
metrics <- calculate_metrics(flow_futures, recompile = FALSE)

# list all possible future scenarios 
#   (81 combinations per species and waterbody)
scenario_options <- flow_futures |>
  distinct(species, waterbody, future, scenario) |>
  mutate(future_next = future, scenario_next = scenario) |>
  complete(nesting(species, waterbody), future, future_next, scenario, scenario_next)

# for each future, grab the correct 2024 metrics and then grab the
#   2025 metrics based on the *_next settings, but replace antecedent
#   for that scenario with the appropriate 2024 value
metrics_2024 <- metrics |>
  filter(water_year == max(water_year))
metrics_2025 <- metrics_2024 |> mutate(water_year = water_year + 1)
metrics_future <- scenario_options |>
  left_join(
    metrics_2024,
    by = c("species", "waterbody", "future", "scenario")
  ) |>
  left_join(
    metrics_2025,
    by = c("species", "waterbody", "future_next" = "future", "scenario_next" = "scenario"),
    suffix = c("_2024", "_2025")
  ) |>
  mutate(
    proportional_antecedent_flow_2025 = proportional_annual_flow_2024,
    proportional_max_antecedent_2025 = proportional_max_annual_2024
  )
metrics_2024 <- metrics_future |>
  select(species, waterbody, future, future_next, scenario, scenario_next, contains("2024")) |>
  rename(water_year = water_year_2024) |>
  rename_with(\(x) gsub("_2024", "", x), contains("2024"))
metrics_2025 <- metrics_future |>
  select(species, waterbody, future, future_next, scenario, scenario_next, contains("2025")) |>
  rename_with(\(x) gsub("_2025", "", x), contains("2025"))
metrics_future <- bind_rows(metrics_2024, metrics_2025)

# load data wihtout e-flows (check reaches and filter to those?)
# Ovens, Upper Loddon, and Broken R Reach 3 have no e-water (very small amounts)
counterfactual_lu <- c(
  "lowerbrokenc_r4" = "broken_creek_r4",
  "campaspe_r4" = "campaspe_river_r4",
  "glenelg_r1b" = "glenelg_river_r1",
  "glenelg_r2" = "glenelg_river_r2",
  "glenelg_r3" = "glenelg_river_r3",
  "goulburn_r4" = "goulburn_river_r4",
  "loddon_r4" = "loddon_river_r4",
  "macalister_r1" = "macalister_river_r1",
  "mackenzie_r3" = "mackenzie_river_r3",
  "moorabool_r3b" = "moorabool_river_r3",
  "thomson_r3" = "thomson_river_r3"
)

# load counterfactual data and filter to target systems
counterfactual <- read.csv("data/vic-ewater-collated-long.csv")

# reformat counterfactuals to match flow_futures and remove negative
#   counterfactual flow values
counterfactual <- counterfactual |>
  mutate(
    date_formatted = parse_date_time(date_formatted, orders = c("ymd_HMS", "ymd")),
    waterbody = counterfactual_lu[system],
    scenario = ifelse(scenario == "No e-water", "counterfactual", scenario),
    water_temperature_c = NA,
    future = NA,
    stream_discharge_mld = discharge_mld,
    stream_discharge_mld = ifelse(stream_discharge_mld < 0, 0, stream_discharge_mld)
  ) |>
  filter(scenario == "counterfactual", waterbody %in% counterfactual_lu) |>
  as_tibble() |>
  left_join(
    fetch_hypoxia_risk(),
    by = "waterbody"
  ) |>
  left_join(
    fetch_carrying_capacity(), 
    by = "waterbody",
    relationship = "many-to-many"
  ) |>
  select(all_of(colnames(flow_futures)))

# re-calculate flow metrics for counterfactual scenarios
metrics_counterfactual <- calculate_metrics(
  counterfactual,
  reference = flow_futures,
  recompile = FALSE,
  suffix = "counterfactual"
)

# then left_join this to the pre-2024 metrics set and fill NAs 
#    with observed values (i.e., no counterfactual implies no e-water)
metrics_observed <- metrics |>
  filter(future == "ave", scenario == "none", water_year < 2024) |>
  select(-future, -scenario) |>
  left_join(
    metrics_counterfactual |> select(-future, -scenario),
    by = c("species", "waterbody", "water_year"),
    suffix = c("", "_counterfactual")
  ) |>
  mutate(
    across(
      contains("_counterfactual"),
      .fns = \(x) ifelse(is.na(x), get(gsub("_counterfactual", "", cur_column())), x)
    ),
    hypoxia_risk_temp_counterfactual = hypoxia_risk_temp,
    hypoxia_risk_counterfactual = ifelse(
      hypoxia_risk_discharge_counterfactual & hypoxia_risk_temp_counterfactual,
      1,
      0
    )
  )
metrics_counterfactual <- metrics_observed |>
  select(species, waterbody, water_year, contains("_counterfactual")) |>
  rename_with(\(x) gsub("_counterfactual", "", x), contains("_counterfactual"))
metrics_observed <- metrics_observed |>
  select(!contains("_counterfactual"))

# add extra metrics for rainbowfish (redfin, gambusia, instream_cover)
#   and for blackfish (overhanging veg cover, instream_cover)
veg_overhang <- read.csv("data/eflow_veg_overhang_231113.csv")
veg_overhang <- veg_overhang |>
  as_tibble() |>
  mutate(
    waterbody = gsub(" ", "_", tolower(waterbody)),
    waterbody = paste(waterbody, vewh_reach, sep = "_r")
  ) |>
  filter(waterbody %in% unique(metrics_observed$waterbody)) |>
  mutate(proportion_overhang = overhang_area_m2 / buff5m_area_m2)
instream_cover <- read.csv("data/vefmap_habitat_site_metrics_231208.csv")
instream_cover <- instream_cover |>
  left_join(
    fetch_site_info(instream_cover) |> 
      select(id_site, reach_no) |> 
      collect(),
    by = "id_site"
  )
instream_cover <- instream_cover |>
  group_by(waterbody, reach_no) |>
  summarise(
    iwh = median(cur_pred_m3_m2, na.rm = TRUE),
    proportion_overhang = median(mean_vo_percentage, na.rm = TRUE) / 100
  ) |>
  ungroup() |>
  mutate(
    waterbody = gsub(" ", "_", tolower(waterbody)),
    waterbody = paste(waterbody, reach_no, sep = "_r")
  )
veg_overhang <- veg_overhang |>
  select(waterbody, proportion_overhang) |>
  full_join(
    instream_cover |> select(-reach_no),
    by = "waterbody",
    suffix = c("", "_hab_metrics")
  ) |>
  mutate(
    proportion_overhang_hab_metrics = ifelse(
      is.na(proportion_overhang_hab_metrics),
      1.2 * proportion_overhang, 
      proportion_overhang_hab_metrics
    ),
    iwh = ifelse(waterbody == "loddon_river_r2", 0.0199, iwh),
    iwh = ifelse(waterbody == "macalister_river_r1", 0.00716, iwh),
    iwh = ifelse(waterbody == "mackenzie_river_r3", 0.00886, iwh),
    iwh = ifelse(waterbody == "ovens_river_r5", max(iwh, na.rm = TRUE), iwh),
    iwh = ifelse(waterbody == "moorabool_river_r0", 0.0145, iwh),
  ) |>
  bind_rows(
    tibble(
      waterbody = c("broken_creek_r4", "loddon_river_r4"),
      proportion_overhang = rep(NA, 2),
      iwh = c(0.0128, 0.00471),
      proportion_overhang_hab_metrics = c(0.797, 0.641)
    )
  ) |>
  mutate(
    iwh = rescale_fn(iwh, base = 0.75, na.rm = TRUE),
    proportion_overhang_hab_metrics = rescale_fn(
      proportion_overhang_hab_metrics,
      base = 0.85,
      na.rm = TRUE
    )
  )

# add overhang/iwh and exotic species to metrics data sets
metrics_observed <- metrics_observed |>
  left_join(
    cpue_exotic |>
      mutate(
        species = "redfin",
        species = ifelse(
          grepl("gambusia", scientific_name, ignore.case = TRUE),
          "gambusia", 
          species
        )
      ) |>
      pivot_wider(
        id_cols = c(waterbody, survey_year),
        names_from = species,
        values_from = presence,
        names_prefix = "presence_"
      ),
    by = c("waterbody", "water_year" = "survey_year")
  ) |>
  mutate(
    presence_redfin = ifelse(is.na(presence_redfin), 0, presence_redfin),
    presence_gambusia = ifelse(is.na(presence_gambusia), 0, presence_gambusia)
  ) |>
  left_join(
    veg_overhang |>
      select(waterbody, iwh, proportion_overhang_hab_metrics) |>
      rename(
        instream_cover = iwh, 
        veg_overhang = proportion_overhang_hab_metrics
      ),
    by = "waterbody"
  )
metrics_counterfactual <- metrics_counterfactual |>
  left_join(
    metrics_observed |> 
      select(species, waterbody, water_year, contains("presence_"), instream_cover, veg_overhang),
    by = c("species", "waterbody", "water_year")
  )

# repeat for future metrics, but just use most recent observations to project
#   future occurrences
max_year <- cpue_exotic |>
  group_by(waterbody) |>
  summarise(max_year = max(survey_year))
metrics_future <- metrics_future |>
  left_join(
    cpue_exotic |>
      mutate(
        species = "redfin",
        species = ifelse(
          grepl("gambusia", scientific_name, ignore.case = TRUE),
          "gambusia", 
          species
        )
      ) |>
      pivot_wider(
        id_cols = c(waterbody, survey_year),
        names_from = species,
        values_from = presence,
        names_prefix = "presence_"
      ) |>
      left_join(max_year, by = "waterbody") |>
      filter(survey_year == max_year) |>
      select(contains("presence_"), waterbody),
    by = c("waterbody")
  ) |>
  mutate(
    presence_redfin = ifelse(is.na(presence_redfin), 0, presence_redfin),
    presence_gambusia = ifelse(is.na(presence_gambusia), 0, presence_gambusia)
  ) |>
  left_join(
    veg_overhang |>
      select(waterbody, iwh, proportion_overhang_hab_metrics) |>
      rename(
        instream_cover = iwh, 
        veg_overhang = proportion_overhang_hab_metrics
      ),
    by = "waterbody"
  )



# extract carrying capacity for species and reach
metrics_tmp <- metrics_observed |> filter(waterbody == "goulburn_river_r4")
k_mc <- flow_futures |> 
  filter(species == "maccullochella_peelii", waterbody == "goulburn_river_r4") |>
  pull(carrying_capacity) |>
  unique()
k_rb <- flow_futures |> 
  filter(species == "melanotaenia_fluviatilis", waterbody == "goulburn_river_r4") |>
  pull(carrying_capacity) |>
  unique()
k_cc <- flow_futures |> 
  filter(species == "cyprinus_carpio", waterbody == "goulburn_river_r4") |>
  pull(carrying_capacity) |>
  unique()

# specify initial conditions
initial <- specify_initial_conditions(
  species = "maccullochella_peelii",
  waterbody = "goulburn_river_r4",
  cpue = cpue,
  start = min(metrics_tmp$water_year),
  nsim = nsim,
  k = k_mc
)

# grab stocking info if required (default to zero, otherwise)
n_stocked <- rep(0, nrow(metrics_tmp))
system_lu <- c(
  "broken_creek_r4" = "Broken Creek",
  "broken_river_r3" = "Broken River",
  "campaspe_river_r4" = "Campaspe River",
  "goulburn_river_r4" = "Goulburn River",
  "loddon_river_r4" = "Loddon River",
  "ovens_river_r5" = "Ovens River"
)
stocking_rates <- stocking |>
  filter(
    Species == "Murray Cod",
    System == system_lu["goulburn_river_r4"]
  ) |>
  select(Year, Number) |>
  mutate(Year = Year + 1) |>
  rename(water_year = Year, number_stocked = Number)
n_stocked <- metrics_tmp |>
  filter(species == "maccullochella_peelii") |>
  select(water_year) |>
  left_join(stocking_rates, by = c("water_year")) |>
  mutate(number_stocked = ifelse(is.na(number_stocked), 0, number_stocked)) |>
  pull(number_stocked)

# initialise population model for this species and waterbody
mc <- specify_pop_model(
  species = "maccullochella_peelii",
  waterbody = "goulburn_river_r4",
  ntime = nrow(metrics_tmp), 
  nstocked = n_stocked,
  k = k_mc
)
rb <- specify_pop_model(
  species = "melanotaenia_fluviatilis",
  waterbody = "goulburn_river_r4",
  ntime = nrow(metrics_tmp), 
  k = k_rb
)
cc <- specify_pop_model(
  species = "cyprinus_carpio",
  waterbody = "goulburn_river_r4",
  ntime = nrow(metrics_tmp), 
  k = k_cc
)

# carp negatively affect smaller MC
mask_lt4 <- transition(mc$dynamics$matrix, dim = 1:4)
fun_lt4 <- function(x, n) {
  # n is the population vector of carp
  x * exp(-sum(n[3:28]) / 100000)
}

# carp benefit from smaller MC
mask_all <- transition(cc$dynamics$matrix)
fun_all <- function(x, n) {
  # n is the population vector of MC
  x / (1 + x * sum(n[1:4]) / 100000)
}

# carp reduce MC recruitment
mask_rec <- reproduction(mc$dynamics$matrix)
fun_rec <- function(x, n) {
  # n is the population vector of carp
  x * exp(-sum(n[3:28]) / 100000)
}

# carp reduce RB abundance
mask_rb_cc <- transition(rb$dynamics$matrix)
fun_rb_cc <- function(x, n) {
  # n is the population vector of carp
  x * exp(-sum(n[3:28]) / 100000)
}

# MC reduce RB abundance
mask_rb_mc <- transition(rb$dynamics$matrix)
fun_rb_mc <- function(x, n) {
  # n is the population vector of cod
  x * exp(-sum(n[3:50]) / 20000)
}

# big MC reduce survival of small carp
mask_lt5 <- transition(cc$dynamics$matrix, dim = 1:5)
fun_lt5 <- function(x, n) {
  # n is the population vector of MC
  x * exp(-sum(n[3:50]) / 20000)
}

# combine masks and functions into pairwise_interaction objects
mc_lt4 <- pairwise_interaction(mc$dynamics, cc$dynamics, mask_lt4, fun_lt4)
cc_all <- pairwise_interaction(cc$dynamics, mc$dynamics, mask_all, fun_all)
mc_rec <- pairwise_interaction(mc$dynamics, cc$dynamics, mask_rec, fun_rec)
rb_cc <- pairwise_interaction(rb$dynamics, cc$dynamics, mask_rb_cc, fun_rb_cc)
rb_mc <- pairwise_interaction(rb$dynamics, mc$dynamics, mask_rb_mc, fun_rb_mc)
cc_lt5 <- pairwise_interaction(cc$dynamics, mc$dynamics, mask_lt5, fun_lt5)

# compile a multispecies dynamics object
mc_cc_rb_dyn <- multispecies(mc_lt4, cc_all, mc_rec, rb_cc, rb_mc, cc_lt5)

# add covariates
metrics_tmp <- metrics_tmp |>
  mutate(
    flow_variability = spawning_flow_variability,
    floodplain_access = ifelse(proportional_spring_flow > 5, 1, 0)
  )

# rename a few metrics for some species
metrics_tmp <- metrics_tmp |>
  mutate(antecedent_flow = proportional_antecedent_flow)
metrics_tmp <- metrics_tmp |>
  mutate(blackwater_risk = hypoxia_risk)
metrics_tmp <- metrics_tmp |>
  mutate(redfin = presence_redfin, gambusia = presence_gambusia)

# TODO: wrap arg formatting in a function
# TODO: add stocking step to the MC model correctly
# TODO: work out how to set initial conditions appropriately
#  (is a list of inits, one for each species -- check order)
# TODO: setup args so they match hex codes of pop dynamics object
# TODO: split metrics_tmp into each species term
# TODO: set up a loop over all waterbodies
# TODO: set up blackfish variant

metrics_tmp_mc <- metrics_tmp |>
  filter(species == "maccullochella_peelii") 
mc_args <- list(
  density_dependence_n = mc$arguments$density_dependence_n,
  covariates = c(
    format_covariates(
      metrics_tmp |>
        filter(species == "maccullochella_peelii") |>
        select(all_of(get_metric_names("maccullochella_peelii")))
    ),
    list(threshold = 0.05),
    list(coefs = get_coefs("maccullochella_peelii", "goulburn_river_r4"))
  ),
  density_dependence = list(
    kdyn = lapply(
      seq_len(nrow(metrics_tmp_mc)),
      \(.x) metrics_tmp_mc$carrying_capacity[.x]
    )
  )
)
cc_args <- list(
  density_dependence_n = cc$arguments$density_dependence_n,
  covariates = c(
    format_covariates(
      metrics_tmp |>
        filter(species == "cyprinus_carpio") |>
        select(all_of(get_metric_names("cyprinus_carpio")))
    ),
    list(coefs = get_coefs("cyprinus_carpio", "goulburn_river_r4"))
  )
)
rb_args <- list(
  covariates = c(
    format_covariates(
      metrics_tmp |>
        filter(species == "melanotaenia_fluviatilis") |>
        select(all_of(get_metric_names("melanotaenia_fluviatilis")))
    ),
    list(coefs = get_coefs("melanotaenia_fluviatilis", "goulburn_river_r4")[1:3]),
    list(warmwater_coefficient = get_coefs("melanotaenia_fluviatilis", "goulburn_river_r4")[4]),
    list(coldwinter_coefficient = get_coefs("melanotaenia_fluviatilis", "goulburn_river_r4")[5])
  )
)

# simulate population dynamics
sims <- simulate(
  mc_cc_rb_dyn,
  nsim = nsim,
  args = list(
    mc_args,
    cc_args,
    rb_args
  ),
  options = list(
    ntime = nrow(metrics_tmp),
    update = update_binomial_leslie,
    tidy_abundances = floor
  )
)


# TODO: use below approach to initalise/burn-in the model inits
# initial <- simulate_scenario(
#   species = species_list[i],
#   x = pop, 
#   nsim = nsim, 
#   init = initial,
#   metrics = metrics_observed_wb[1, ],
#   coefs = get_coefs(species_list[i], waterbodies[j]),
#   nburnin = nburnin - 1
# )
# sims_observed <- simulate_scenario(
#   species = species_list[i],
#   x = pop, 
#   nsim = nsim, 
#   init = initial[, , dim(initial)[3]],
#   metrics = metrics_observed_wb,
#   coefs = get_coefs(species_list[i], waterbodies[j]),
#   nburnin = 0
# )
