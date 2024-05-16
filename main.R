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
#     As for single-species scenarios
#
# Flow scenarios consider observed daily flows with and without 
#    e-water allocations and near-term (to 2025) forecasts of flows under
#    different plausible climates and flow management strategies

# Author: Jian Yen (jian.yen [at] deeca.vic.gov.au)
# Last updated: 17 May 2024

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
library(patchwork)

# and some load helpers
source("R/utils.R")
source("R/interactions.R")
source("R/fish.R")
source("R/flow.R")
source("R/validation.R")

# settings
set.seed(2024-05-01)
nsim <- 100
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

# loop over all waterbodies and fit a model for each if required
wb_list <- metrics_observed |> pull(waterbody) |> unique()
northern_rivers <- c(
  "broken_creek_r4", "broken_river_r3", "campaspe_river_r4", 
  "goulburn_river_r4", "loddon_river_r4", "ovens_river_r5"
)
sims_exist <- any(grepl("^sims-", dir("outputs/simulations/")))
if (simulate_again | !sims_exist) {
  for (i in seq_along(wb_list)) {
    
    # pull out metrics for the target waterbody
    metrics_wb <- metrics_observed |> filter(waterbody == !!wb_list[i])
    metrics_wb_cf<- metrics_counterfactual |> filter(waterbody == !!wb_list[i])
    
    # check for NA values in the metrics, fill with waterbody median if needed
    metrics_wb <- metrics_wb |>
      mutate(across(contains("nday_"), fill_na))
    metrics_wb_cf <- metrics_wb_cf |>
      mutate(across(contains("nday_"), fill_na))
    
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
    
    # compile a multispecies dynamics object for three levels of interaction
    #   strength (note: scale affects carrying capacities, so halving scale 
    #   doubles the strength of interactions [approximately])
    interactions <- specify_interactions(pop_list, scale = 2)
    interactions_double <- specify_interactions(pop_list, scale = 1)
    interactions_half <- specify_interactions(pop_list, scale = 4)
    mspop <- do.call(multispecies, interactions)
    mspop_double <- do.call(multispecies, interactions_double)
    mspop_half <- do.call(multispecies, interactions_half)
    
    # check the order of species in the ms object
    species_order <- match_pops(mspop, pop_list)
    
    # prepare args and inits
    args_burnin <- prepare_args(
      metrics = metrics_wb,
      waterbody = wb_list[i],
      pops = pop_list,
      nburnin = nburnin
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
    sims_double <- simulate_scenario(
      x = mspop_double,
      nsim = nsim,
      init = inits, 
      args = args_observed,
      args2 = args_counterfactual,
      args_burnin = args_burnin
    )
    sims_half <- simulate_scenario(
      x = mspop_half,
      nsim = nsim,
      init = inits, 
      args = args_observed,
      args2 = args_counterfactual,
      args_burnin = args_burnin
    )
    
    # save outputs
    qsave(sims, file = paste0("outputs/simulations/sims-", wb_list[i], ".qs"))
    qsave(sims_double, file = paste0("outputs/simulations/sensitivity-double-sims-", wb_list[i], ".qs"))
    qsave(sims_half, file = paste0("outputs/simulations/sensitivity-half-sims-", wb_list[i], ".qs"))
    
    # extract initial conditions for forecasts from sims_observed
    #  but use the no-interactions scenario because interactions
    #  has n = 0 for many cases
    init_future <- lapply(sims[[2]], \(x) x[, , dim(x)[3]])
    
    # simulate futures
    future_sub <- metrics_future |>
      filter(waterbody == wb_list[i]) |>
      distinct(future, future_next, scenario, scenario_next)
    
    # run through each future and generate forecasts
    for (j in seq_len(nrow(future_sub))) {
      
      # pull out metrics for a given scenario and calculate arguments
      metrics_future_sub <- metrics_future |>
        filter(
          future == future_sub$future[j],
          future_next == future_sub$future_next[j],
          scenario == future_sub$scenario[j],
          scenario_next == future_sub$scenario_next[j]
        )
      args_future <- prepare_args(
        metrics = metrics_future_sub,
        waterbody = wb_list[i],
        pops = pop_list
      )
      
      # simulate under the specific scenario
      sims_future <- simulate_scenario(
        x = mspop,
        nsim = nsim,
        init = init_future, 
        args = args_future,
        args2 = NULL,
        args_burnin = NULL
      )
      
      # and save output
      future_name <- paste(
        future_sub$future[j],
        future_sub$future_next[j],
        future_sub$scenario[j],
        future_sub$scenario_next[j],
        sep = "_"
      )
      qsave(
        sims_future, 
        file = paste0(
          "outputs/simulations/future-", wb_list[i], "-", future_name, ".qs"
        )
      )
      
    }
    
  } 
  
}

# load main sims files
sim_files <- dir("outputs/simulations/")
sim_files <- sim_files[grepl("^sims-", sim_files)]
sims_main <- lapply(sim_files, \(x) qread(paste0("outputs/simulations/", x)))

# model CPUE using an AR1 model to estimate values (simple AR1 model with
#   random terms to soak up variation)
iter <- 2000
warmup <- 1000
chains <- 3
cores <- 3
use_cached <- TRUE
cpue_mc <- estimate_cpue(
  x = cpue, 
  use_cached = use_cached,
  species = "Maccullochella peelii",
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores
)                       
cpue_recruit_mc <- estimate_cpue(
  x = cpue_recruits, 
  recruit = TRUE,
  use_cached = use_cached,
  species = "Maccullochella peelii",
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores
)                       
cpue_adult_mc <- estimate_cpue(
  x = cpue_adults, 
  adult = TRUE,
  use_cached = use_cached,
  species = "Maccullochella peelii",
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores
)                       
cpue_bf <- estimate_cpue(
  x = cpue, 
  use_cached = use_cached,
  species = "Gadopsis marmoratus",
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores
)                       
cpue_recruit_bf <- estimate_cpue(
  x = cpue_recruits, 
  recruit = TRUE,
  use_cached = use_cached,
  species = "Gadopsis marmoratus",
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores
)       
cpue_cc <- estimate_cpue(
  x = cpue, 
  use_cached = use_cached,
  species = "Cyprinus carpio",
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores
)                       
cpue_recruit_cc <- estimate_cpue(
  x = cpue_recruits, 
  recruit = TRUE,
  use_cached = use_cached,
  species = "Cyprinus carpio",
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores
)                       
cpue_rb <- estimate_cpue(
  x = cpue, 
  use_cached = use_cached,
  species = "Melanotaenia fluviatilis",
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores
)                       

# posterior checks (saved to figures)
pp_mc <- pp_check(cpue_mc) + scale_x_log10() + 
  xlab("Catch (total)") + ylab("Density") + 
  theme(legend.position = "none")
pp_bf <- pp_check(cpue_bf) + scale_x_log10() + 
  xlab("Catch (total)") + ylab("Density") +
  theme(legend.position = "none")
pp_mc_adult <- pp_check(cpue_adult_mc) + scale_x_log10() +
  xlab("Catch (adults)") + ylab("Density") +
  theme(legend.position = "none")
pp_mc_recruit <- pp_check(cpue_recruit_mc) + scale_x_log10() +
  xlab("Catch (young of year)") + ylab("Density") +
  theme(legend.position = "none")
pp_bf_recruit <- pp_check(cpue_recruit_bf) + scale_x_log10() +
  xlab("Catch (young of year)") + ylab("Density") +
  theme(legend.position = "none")
pp_cc <- pp_check(cpue_cc) + scale_x_log10() + 
  xlab("Catch (total)") + ylab("Density") + 
  theme(legend.position = "none")
pp_cc_recruit <- pp_check(cpue_recruit_cc) + scale_x_log10() +
  xlab("Catch (young of year)") + ylab("Density") + 
  theme(legend.position = "none")
pp_rb <- pp_check(cpue_rb) + scale_x_log10() + 
  xlab("Catch (total)") + ylab("Density") +
  theme(legend.position = "none")
pp_all <- (pp_mc | pp_mc_recruit) /
  (pp_mc_adult | pp_bf) /
  (pp_bf_recruit | pp_rb) +
  (pp_cc | pp_cc_recruit) +
  plot_annotation(tag_levels = "a")
ggsave(
  filename = "outputs/figures/pp-checks.png",
  plot = pp_all,
  device = ragg::agg_png,
  width = 6,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

# validation on main sims files (sims_main[[i]][[1]])
names(sims_main) <- gsub("sims-|\\.qs", "", sim_files)
mc_sim_metrics <- calculate_val_metrics(
  x = sims_main[northern_rivers],
  cpue_mod = cpue_mc,
  subset = 1:50, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 1
)
mc_sim_metrics_adult <- calculate_val_metrics(
  x = sims_main[northern_rivers],
  cpue_mod = cpue_adult_mc,
  subset = 5:50, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 1
)
mc_sim_metrics_recruit <- calculate_val_metrics(
  x = sims_main[northern_rivers],
  cpue_mod = cpue_recruit_mc,
  recruit = TRUE,
  subset = 1, 
  sim_years = (min(metrics_observed$water_year) - 1L):max(metrics_observed$water_year),
  sp_idx = 1
)
bf_sim_metrics <- calculate_val_metrics(
  x = sims_main[!names(sims_main) %in% northern_rivers],
  cpue = cpue_bf,
  subset = 1:11, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 1
)
bf_sim_metrics_recruit <- calculate_val_metrics(
  x = sims_main[!names(sims_main) %in% northern_rivers],
  cpue_mod = cpue_recruit_bf,
  recruit = TRUE,
  subset = 1, 
  sim_years = (min(metrics_observed$water_year) - 1L):max(metrics_observed$water_year),
  sp_idx = 1
)
cc_sim_metrics <- calculate_val_metrics(
  x = sims_main,
  cpue = cpue_cc,
  subset = 1:28,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 2
)
cc_sim_metrics_recruit <- calculate_val_metrics(
  x = sims_main,
  cpue_mod = cpue_recruit_cc,
  recruit = TRUE,
  subset = 1, 
  sim_years = (min(metrics_observed$water_year) - 1L):max(metrics_observed$water_year),
  sp_idx = 2
)
rb_sim_metrics <- calculate_val_metrics(
  x = sims_main[northern_rivers],
  cpue = cpue_rb,
  subset = 1:7, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 3
)
sim_metrics <- bind_rows(
  mc_sim_metrics |> mutate(species = "Murray Cod"),
  mc_sim_metrics_recruit |> mutate(species = "Murray Cod (young of year)"),
  mc_sim_metrics_adult |> mutate(species = "Murray Cod (adults)")
)
metrics_plot_mc <- plot_metric(sim_metrics)
metrics_plot_rb <- plot_metric(rb_sim_metrics |> mutate(species = "Murray-Darling Rainbowfish"))
metrics_plot_bf <- plot_metric(
  bind_rows(
    bf_sim_metrics |> mutate(species = "River Blackfish"),
    bf_sim_metrics_recruit |> mutate(species = "River Blackfish (young of year)")
  )
)
metrics_plot_cc <- plot_metric(
  bind_rows(
    cc_sim_metrics |> mutate(species = "Common Carp"),
    cc_sim_metrics_recruit |> mutate(species = "Common Carp (young of year)")
  )
)
ggsave(
  filename = "outputs/figures/metrics-mc.png",
  plot = metrics_plot_mc,
  device = ragg::agg_png,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/metrics-bf.png",
  plot = metrics_plot_bf,
  device = ragg::agg_png,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/metrics-rb.png",
  plot = metrics_plot_rb,
  device = ragg::agg_png,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/metrics-cc.png",
  plot = metrics_plot_cc,
  device = ragg::agg_png,
  width = 10,
  height = 6,
  units = "in",
  dpi = 600
)

# plot hindcasts for observed with and without interactions
mc_hindcasts <- plot_hindcasts(
  x = sims_main[northern_rivers],
  cpue = cpue_mc,
  subset = 1:50, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 1
) 
mc_adult_hindcasts <- plot_hindcasts(
  x = sims_main[northern_rivers],
  cpue = cpue_adult_mc,
  subset = 5:50, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 1
) 
mc_recruit_hindcasts <- plot_hindcasts(
  x = sims_main[northern_rivers],
  cpue = cpue_recruit_mc,
  recruit = TRUE,
  subset = 1, 
  sim_years = (min(metrics_observed$water_year) - 1L):max(metrics_observed$water_year),
  sp_idx = 1
) 
bf_hindcasts <- plot_hindcasts(
  x = sims_main[!names(sims_main) %in% northern_rivers],
  cpue = cpue_bf,
  subset = 1:11, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 1
) 
bf_recruit_hindcasts <- plot_hindcasts(
  x = sims_main[!names(sims_main) %in% northern_rivers],
  cpue = cpue_recruit_bf,
  recruit = TRUE,
  subset = 1, 
  sim_years = (min(metrics_observed$water_year) - 1L):max(metrics_observed$water_year),
  sp_idx = 1
) 
rb_hindcasts <- plot_hindcasts(
  x = sims_main[northern_rivers],
  cpue = cpue_rb,
  subset = 1:7, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 3
) 
cc_hindcasts <- plot_hindcasts(
  x = sims_main,
  cpue = cpue_cc,
  subset = 1:28, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 2
) 
cc_recruit_hindcasts <- plot_hindcasts(
  x = sims_main,
  cpue = cpue_recruit_cc,
  recruit = TRUE,
  subset = 1, 
  sim_years = (min(metrics_observed$water_year) - 1L):max(metrics_observed$water_year),
  sp_idx = 2
) 

# save the plots
ggsave(
  filename = "outputs/figures/hindcasts-mc.png",
  plot = mc_hindcasts,
  device = ragg::agg_png,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/hindcasts-mc-adult.png",
  plot = mc_adult_hindcasts,
  device = ragg::agg_png,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/hindcasts-mc-recruit.png",
  plot = mc_recruit_hindcasts,
  device = ragg::agg_png,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/hindcasts-bf.png",
  plot = bf_hindcasts,
  device = ragg::agg_png,
  width = 8,
  height = 7,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/hindcasts-bf-recruit.png",
  plot = bf_recruit_hindcasts,
  device = ragg::agg_png,
  width = 8,
  height = 7,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/hindcasts-cc.png",
  plot = cc_hindcasts,
  device = ragg::agg_png,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/hindcasts-cc-recruit.png",
  plot = cc_recruit_hindcasts,
  device = ragg::agg_png,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/hindcasts-rb.png",
  plot = rb_hindcasts,
  device = ragg::agg_png,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

# compare the four scenarios (obs, obs-noint, cf, cf-noint)
trajectories_mc_adult <- plot_trajectories(
  x = sims_main[northern_rivers], 
  subset = 5:50,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 1
)
trajectories_mc <- plot_trajectories(
  x = sims_main[northern_rivers], 
  subset = 1:50,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 1
)
trajectories_mc_recruit <- plot_trajectories(
  x = sims_main[northern_rivers], 
  subset = 1,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 1
)
trajectories_bf <- plot_trajectories(
  x = sims_main[!names(sims_main) %in% northern_rivers], 
  subset = 1:11,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 1
)
trajectories_bf_recruit <- plot_trajectories(
  x = sims_main[!names(sims_main) %in% northern_rivers], 
  subset = 1,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 1
)
trajectories_cc <- plot_trajectories(
  x = sims_main, 
  subset = 1:28,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 2
)
trajectories_cc_recruit <- plot_trajectories(
  x = sims_main, 
  subset = 1,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 2
)
trajectories_rb <- plot_trajectories(
  x = sims_main[northern_rivers[1:5]], 
  subset = 1:7,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 3
)

# save the plots
ggsave(
  filename = "outputs/figures/trajectories-mc.png",
  plot = trajectories_mc,
  device = ragg::agg_png,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/trajectories-mc-adult.png",
  plot = trajectories_mc_adult,
  device = ragg::agg_png,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/trajectories-mc-recruit.png",
  plot = trajectories_mc_recruit,
  device = ragg::agg_png,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/trajectories-bf.png",
  plot = trajectories_bf,
  device = ragg::agg_png,
  width = 8,
  height = 7,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/trajectories-bf-recruit.png",
  plot = trajectories_bf_recruit,
  device = ragg::agg_png,
  width = 8,
  height = 7,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/trajectories-cc.png",
  plot = trajectories_cc,
  device = ragg::agg_png,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/trajectories-cc-recruit.png",
  plot = trajectories_cc_recruit,
  device = ragg::agg_png,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/trajectories-rb.png",
  plot = trajectories_rb,
  device = ragg::agg_png,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)

# check sensitivity to interaction strength
sims_sens_files <- dir("outputs/simulations/")
sims_sens_files <- sims_sens_files[!grepl("future-", sims_sens_files)]
sims_sens <- lapply(sims_sens_files, \(x) qread(paste0("outputs/simulations/", x)))
names(sims_sens) <- gsub(
  "sims-", "", 
  gsub(
    "-sims-", "_",
    gsub(
      "sensitivity-|\\.qs", "", sims_sens_files
    )
  )
)
sens_double <- lapply(sims_sens[grepl("double_", names(sims_sens))], \(x) x[c(1, 3)])
sens_half <- lapply(sims_sens[grepl("half_", names(sims_sens))], \(x) x[c(1, 3)])
sens_main <- lapply(
  sims_sens[!(grepl("half_", names(sims_sens)) | grepl("double_", names(sims_sens)))],
  \(x) x[c(1, 3)]
)
sens_combined <- mapply(
  \(x, y, z) c(x, y, z),
  x = sens_main,
  y = sens_double,
  z = sens_half, 
  SIMPLIFY = FALSE
)

# plot trajectories under different interaction strengths
sensitivity_mc <- plot_trajectories(
  x = sens_combined[northern_rivers], 
  subset = 1:50,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 1,
  scn = 1:6,
  .names = paste0(
    rep(c("Baseline", "Double", "Half"), each = 2), 
    rep(c(" (with e-water)", " (without e-water)"), times = 3)
  )
)
sensitivity_mc_adult <- plot_trajectories(
  x = sens_combined[northern_rivers], 
  subset = 5:50,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 1,
  scn = 1:6,
  .names = paste0(
    rep(c("Baseline", "Double", "Half"), each = 2), 
    rep(c(" (with e-water)", " (without e-water)"), times = 3)
  )
)
sensitivity_mc_recruit <- plot_trajectories(
  x = sens_combined[northern_rivers], 
  subset = 1,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 1,
  scn = 1:6,
  .names = paste0(
    rep(c("Baseline", "Double", "Half"), each = 2), 
    rep(c(" (with e-water)", " (without e-water)"), times = 3)
  )
)
sensitivity_bf <- plot_trajectories(
  x = sens_combined[!names(sens_combined) %in% northern_rivers], 
  subset = 1:11,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 1,
  scn = 1:6,
  .names = paste0(
    rep(c("Baseline", "Double", "Half"), each = 2), 
    rep(c(" (with e-water)", " (without e-water)"), times = 3)
  )
)
sensitivity_bf_recruit <- plot_trajectories(
  x = sens_combined[!names(sens_combined) %in% northern_rivers], 
  subset = 1,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 1,
  scn = 1:6,
  .names = paste0(
    rep(c("Baseline", "Double", "Half"), each = 2), 
    rep(c(" (with e-water)", " (without e-water)"), times = 3)
  )
)
sensitivity_cc <- plot_trajectories(
  x = sens_combined, 
  subset = 1:28,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 2,
  scn = 1:6,
  .names = paste0(
    rep(c("Baseline", "Double", "Half"), each = 2), 
    rep(c(" (with e-water)", " (without e-water)"), times = 3)
  )
)
sensitivity_cc_recruit <- plot_trajectories(
  x = sens_combined, 
  subset = 1,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 2,
  scn = 1:6,
  .names = paste0(
    rep(c("Baseline", "Double", "Half"), each = 2), 
    rep(c(" (with e-water)", " (without e-water)"), times = 3)
  )
)
sensitivity_rb <- plot_trajectories(
  x = sens_combined[northern_rivers[1:5]], 
  subset = 1:7,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  sp_idx = 3,
  scn = 1:6,
  .names = paste0(
    rep(c("Baseline", "Double", "Half"), each = 2), 
    rep(c(" (with e-water)", " (without e-water)"), times = 3)
  )
)

# save the plots
ggsave(
  filename = "outputs/figures/sensitivity-mc.png",
  plot = sensitivity_mc,
  device = ragg::agg_png,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/sensitivity-mc-adult.png",
  plot = sensitivity_mc_adult,
  device = ragg::agg_png,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/sensitivity-mc-recruit.png",
  plot = sensitivity_mc_recruit,
  device = ragg::agg_png,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/sensitivity-bf.png",
  plot = sensitivity_bf,
  device = ragg::agg_png,
  width = 8,
  height = 7,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/sensitivity-bf-recruit.png",
  plot = sensitivity_bf_recruit,
  device = ragg::agg_png,
  width = 8,
  height = 7,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/sensitivity-cc.png",
  plot = sensitivity_cc,
  device = ragg::agg_png,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/sensitivity-cc-recruit.png",
  plot = sensitivity_cc_recruit,
  device = ragg::agg_png,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/sensitivity-rb.png",
  plot = sensitivity_rb,
  device = ragg::agg_png,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

# load futures
sims_future_files <- dir("outputs/simulations/")
sims_future_files <- sims_future_files[grepl("future-", sims_future_files)]
sims_future <- lapply(sims_future_files, \(x) qread(paste0("outputs/simulations/", x)))
names(sims_future) <- gsub("future-|\\.qs", "", sims_future_files)

# pull out names of each waterbody for matching systems/species
future_waterbody <- sapply(
  names(sims_future),
  \(x) strsplit(x, "-")[[1]][1]
)

# plot the forecasts for all systems for 2024
mc_forecasts <- summarise_forecasts(  
  x = sims_future[future_waterbody %in% northern_rivers], 
  subset = 1:50,
  sp_idx = 1,
  probs = c(0.1, 0.9),
  scn = 1:2,
  .names = c("With interactions", "Without interactions")
)
mc_forecasts_recruits <- summarise_forecasts(  
  x = sims_future[future_waterbody %in% northern_rivers], 
  subset = 1,
  sp_idx = 1,
  probs = c(0.1, 0.9),
  scn = 1:2,
  .names = c("With interactions", "Without interactions")
)
mc_forecasts_adults <- summarise_forecasts(  
  x = sims_future[future_waterbody %in% northern_rivers], 
  subset = 5:50,
  sp_idx = 1,
  probs = c(0.1, 0.9),
  scn = 1:2,
  .names = c("With interactions", "Without interactions")
)
bf_forecasts <- summarise_forecasts(  
  x = sims_future[!future_waterbody %in% northern_rivers], 
  subset = 1:11,
  sp_idx = 1,
  probs = c(0.1, 0.9),
  scn = 1:2,
  .names = c("With interactions", "Without interactions")
)
bf_forecasts_recruits <- summarise_forecasts(  
  x = sims_future[!future_waterbody %in% northern_rivers], 
  subset = 1,
  sp_idx = 1,
  probs = c(0.1, 0.9),
  scn = 1:2,
  .names = c("With interactions", "Without interactions")
)
cc_forecasts <- summarise_forecasts(  
  x = sims_future, 
  subset = 1:28,
  sp_idx = 2,
  probs = c(0.1, 0.9),
  scn = 1:2,
  .names = c("With interactions", "Without interactions")
)
cc_forecasts_recruits <- summarise_forecasts(  
  x = sims_future, 
  subset = 1,
  sp_idx = 2,
  probs = c(0.1, 0.9),
  scn = 1:2,
  .names = c("With interactions", "Without interactions")
)
rb_forecasts <- summarise_forecasts(  
  x = sims_future[future_waterbody %in% northern_rivers], 
  subset = 1:7,
  sp_idx = 3,
  probs = c(0.1, 0.9),
  scn = 1:2,
  .names = c("With interactions", "Without interactions")
)
mc_forecasts_plot <- plot_forecasts(
  x = mc_forecasts, 
  target = 2024
)
mc_forecasts_recruit_plot <- plot_forecasts(
  x = mc_forecasts, 
  target = 2024
)
mc_forecasts_adult_plot <- plot_forecasts(
  x = mc_forecasts_adults, 
  target = 2024
)
bf_forecasts_plot <- plot_forecasts(
  x = bf_forecasts, 
  target = 2024
)
bf_forecasts_recruit_plot <- plot_forecasts(
  x = bf_forecasts_recruits, 
  target = 2024
)
cc_forecasts_plot <- plot_forecasts(
  x = cc_forecasts, 
  target = 2024
)
cc_forecasts_recruit_plot <- plot_forecasts(
  x = cc_forecasts_recruits, 
  target = 2024
)
rb_forecasts_plot <- plot_forecasts(
  x = rb_forecasts, 
  target = 2024
)

# save plots
ggsave(
  filename = "outputs/figures/forecasts-mc.png",
  plot = mc_forecasts_plot,
  device = ragg::agg_png,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/forecasts-mc-recruit.png",
  plot = mc_forecasts_recruit_plot,
  device = ragg::agg_png,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/forecasts-mc-adult.png",
  plot = mc_forecasts_adult_plot,
  device = ragg::agg_png,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/forecasts-bf.png",
  plot = bf_forecasts_plot,
  device = ragg::agg_png,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/forecasts-bf-recruit.png",
  plot = bf_forecasts_recruit_plot,
  device = ragg::agg_png,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/forecasts-cc.png",
  plot = cc_forecasts_plot,
  device = ragg::agg_png,
  width = 10,
  height = 7,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/forecasts-cc-recruit.png",
  plot = cc_forecasts_recruit_plot,
  device = ragg::agg_png,
  width = 10,
  height = 7,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/forecasts-rb.png",
  plot = rb_forecasts_plot,
  device = ragg::agg_png,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
