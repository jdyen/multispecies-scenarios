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

# loop over all waterbodies and fit a model for each
wb_list <- metrics_observed |> pull(waterbody) |> unique()
northern_rivers <- c(
  "broken_creek_r4", "broken_river_r3", "campaspe_river_r4", 
  "goulburn_river_r4", "loddon_river_r4", "ovens_river_r5"
)
for (i in seq_along(wb_list)) {
  
  # pull out metrics for the target waterbody
  metrics_wb <- metrics_observed |> filter(waterbody == !!wb_list[i])
  metrics_wb_cf<- metrics_counterfactual |> filter(waterbody == !!wb_list[i])
  
  # create base population dynamics objects
  pops <- prepare_pop(
    metrics_wb, 
    waterbody = wb_list[i],
    stocking = stocking
  )
  
  # pull out system-appropriate populations/species 
  if (wb_list[i] %in% northern_rivers) {
    pop_list <- list(pops$mc, pops$cc, pops$rb)
  } else {
    pop_list <- list(pops$bf, pops$cc)
  }
  
  # compile a multispecies dynamics object
  interactions <- specify_interactions(pop_list)
  mspop <- do.call(multispecies, interactions)
  
  # check the order of species in the ms object
  species_order <- match_pops(mspop, pop_list)
  
  # prepare args and inits
  args_burnin <- prepare_args(
    metrics = metrics_wb,
    waterbody = wb_list[i],
    pops = pop_list,
    nburnin = 10
  )
  args_observed <- prepare_args(
    metrics = metrics_wb,
    waterbody = wb_list[i],
    pops = pop_list
  )
  args_counterfactual <- prepare_args(
    metrics = metrics_wb_cf,
    waterbody = wb_list[i],
    pops = pop_list
  )
  inits <- prepare_inits(
    waterbody = wb_list[i],
    cpue = cpue,
    metrics = metrics_wb,
    pops = pop_list,
    nsim = nsim
  )
  
  # simulate population dynamics for target system
  sims <- simulate_scenario(
    x = mspop,
    nsim = nsim,
    init = inits, 
    args = args_observed,
    args2 = args_counterfactual,
    args_burnin = args_burnin
  )
  
  # save outputs
  qsave(sims, file = "outputs/simulations/sims-", wb_list[i], ".qs")
  
}

