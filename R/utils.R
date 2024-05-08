# function to tidy up lists by adding names as a column and then binding all rows
add_col_bind <- function(x, name = "scenario") {
  id <- names(x)
  out <- mapply(\(.x, .y) .x |> mutate(category = .y), .x = x, .y = id, SIMPLIFY = FALSE)
  bind_rows(out) |> rename(!!name := category)
}

# function to rescale values between [base, top]
rescale_fn <- function(x, base, top = 1, ...) {
  
  # shift to min = 0
  x <- x - min(x, ...)
  
  # scale to [0, 1]
  x <- x / max(x, ...)
  
  # shrink to [0, top - base]
  x <- (top - base) * x
  
  # and add base
  x <- base + x
  
  # return
  x
  
}

# function to set initial conditions based on assumed constant proportional
#   survival over all age classes, but modified based on survey data
#   to account for shifts to older/younger age classes
set_initial <- function(
    species,
    cpue,
    cpue_max,
    effort_h,
    n = 10000,
    nsim = 1,
    rescale = NULL
) {
  
  # define a population matrix and work out steady state age frequency
  if (species == "gadopsis_marmoratus") {
    adults <- 2:11
    mat <- river_blackfish(k = n)$dynamics$matrix
    initial_age_frequency <- Re(eigen(mat)$vectors[, 1])
  }
  if (species == "maccullochella_peelii") {
    adults <- 5:30
    mat <- murray_cod(k = n)$dynamics$matrix
    initial_age_frequency <- Re(eigen(mat)$vectors[, 1])
  }
  if (species == "melanotaenia_fluviatilis") {
    adults <- 2:5
    mat <- murray_rainbowfish(k = n)$dynamics$matrix
    initial_age_frequency <- Re(eigen(mat)$vectors[, 1])
  }
  
  # work out proportion of fish relative to max catch for a reach
  n_scaled <- 0
  if (cpue_max > 0)
    n_scaled <- n * (cpue / cpue_max)
  
  # rescale if required
  if (!is.null(rescale))
    n_scaled <- 0.6 * cpue * rescale  # current assumption is 60% of catch is adults
  
  # set a min value if cpue_adult = 0
  if (n_scaled == 0) {
    if (cpue_max > 0) {
      cpue_min <- min(c(cpue_max, 1)) / effort_h
      n_scaled <- n * (cpue_min / cpue_max)
    } else {
      n_scaled <- 1 / effort_h
    }
  }
  
  # standardise this for adult abundances only
  initial_age_frequency <- initial_age_frequency / sum(initial_age_frequency[adults])
  
  # return
  matrix(
    rpois(
      nsim * length(initial_age_frequency),
      lambda = n_scaled * initial_age_frequency
    ),
    nrow = nsim,
    byrow = TRUE
  )
  
}

# function to initialise populations
specify_initial_conditions <- function(species, waterbody, cpue, start, nsim, k, ...) {
  
  # stop if species isn't one of the three targets
  stopifnot(
    species %in% c("gadopsis_marmoratus", "maccullochella_peelii", "melanotaenia_fluviatilis")
  )
  
  # rename waterbody to avoid conflicts below
  wb <- waterbody
  
  # pull out CPUE for the target species
  cpue_sub <- cpue |> 
    mutate(
      sciname = tolower(gsub(" ", "_", scientific_name)),
      waterbody = paste(
        tolower(gsub(" ", "_", waterbody)),
        reach_no,
        sep = "_r"
      )
    ) |>
    filter(
      sciname == species,
      waterbody == wb
    ) |>
    mutate(cpue_max = max(catch / effort_h))
  
  # work out the minimum starting year that's actually in the data
  if (start %in% cpue_sub$survey_year) {
    cpue_sub <- cpue_sub |>
      filter(survey_year == start)
  } else {
    start <- abs(start - unique(cpue_sub$survey_year))
    start <- unique(cpue_sub$survey_year)[which.min(start)]
    cpue_sub <- cpue_sub |>
      filter(survey_year == start)
  }
  
  # calculate total CPUE per reach
  cpue_sub <- cpue_sub |>
    group_by(survey_year) |>
    summarise(
      cpue_max = unique(cpue_max),
      catch = sum(catch),
      effort_h = sum(effort_h)
    ) |>
    mutate(
      cpue = catch / effort_h
    )
  
  # pull out reach length for target waterbody
  reach_length <- .vewh_reach_lengths |> 
    filter(waterbody == wb) |>
    pull(reach_length) |>
    as.numeric()
  
  # k is 20% of capacity for river blackfish (why?)
  if (species == "gadopsis_marmoratus") 
    k <- 0.2 * k
  
  # work out the rescaling rate for catch to abundance based on a capture
  #    rate of 40%
  adult_rescale <- (1 / 0.4) * (reach_length / 100)
  init <- set_initial(
    species = species,
    cpue = cpue_sub$cpue,
    cpue_max = cpue_sub$cpue_max,
    effort_h = cpue_sub$effort_h,
    n = k,
    nsim = nsim,
    rescale = adult_rescale
  )
  
  # return
  init
  
}

# function to set up all three population models
specify_pop_model <- function(species, waterbody, ntime, nstocked, k, ...) {
  
  # stop if species isn't one of the three targets
  stopifnot(
    species %in% 
      c(
        "cyprinus_carpio",
        "gadopsis_marmoratus", 
        "maccullochella_peelii", 
        "melanotaenia_fluviatilis"
      )
  )
  
  # river blackfish model
  if (species == "gadopsis_marmoratus") {
    mod <- river_blackfish(
      k = k,
      ntime = ntime
    )
  }
  
  # murray cod model
  if (species == "maccullochella_peelii") {
    
    # use a lookup to define system from waterbody
    system <- switch(
      waterbody,
      "broken_creek_r4" = "broken_creek",
      "broken_river_r3" = "broken_river",
      "campaspe_river_r4" = "campaspe_river",
      "goulburn_river_r4" = "goulburn_river",
      "loddon_river_r4" = "campaspe_river",
      "ovens_river_r5" = "ovens_river",
      "murray_river"
    )
    
    mod <- murray_cod(
      k = k,
      system = system,
      n = list(
        # number stocked, accounting for fingerling mortality and 50:50 sex ratio
        0.5 * 0.38 * 0.31 * nstocked[seq_len(ntime)],
        rep(0, ntime),
        rep(0, ntime)
      ),
      ntime = ntime, 
      start = rep(1, 3), 
      end = rep(ntime, 3), 
      add = c(TRUE, TRUE, TRUE),
      p_capture = 0.1,   # 10% capture probability for any fish in slot
      slot = c(550, 750) # slot in mm
    )
  }
  
  # rainbowfish model
  if (species == "melanotaenia_fluviatilis") {
    mod <- murray_rainbowfish(
      k = k,
      ntime = ntime
    )
  }
  
  # carp model
  if (species == "cyprinus_carpio") {
    mod <- common_carp(
      k = k,
      system = "main_channel",
      ntime = ntime
    )
    
  }
  
  # return
  mod
  
}

# function to return coefficients for each species
get_coefs <- function(species, waterbody) {
  coefs <- list(
    "gadopsis_marmoratus" = list(
      "glenelg_river_r1" = c(-0.1, 50, 50, 100, 75, 0.1, 0.05),
      "glenelg_river_r2" = c(-0.1, 50, 50, 100, 50, 0.1, 0.05),
      "glenelg_river_r3" = c(-1, 50, 50, 100, 50, 0.1, 0.05),
      "loddon_river_r2" = c(-1, 50, 25, 70, 80, 0.1, 0.05),
      "macalister_river_r1" = c(-0.1, 50, 15, 10, 20, 0.1, 0.05),
      "mackenzie_river_r3" = c(-1, 50, 25, 70, 80, 0.1, 0.05),
      "moorabool_river_r3" = c(-0.1, 50, 40, 40, 30, 0.1, 0.05),
      "thomson_river_r3" = c(-1, 50, 20, 10, 15, 0.1, 0.05)
    ),
    "maccullochella_peelii" = list(
      "broken_creek_r4" = c(-10, 30, 6, -30, 100, 100),
      "broken_river_r3" = c(-15, 35, 45, -30, 80, 25),
      "campaspe_river_r4" = c(-60, 10, 20, -15, 45, 10),
      "goulburn_river_r4" = c(-5, 20, 6, -10, 30, 10),
      "loddon_river_r4" = c(-30, 20, 10, -10, 30, 10),
      "ovens_river_r5" = c(-10, 60, 20, -30, 60, 25)
    ),
    "melanotaenia_fluviatilis" = list(
      "broken_creek_r4" = c(-0.05, 40, 20, 0.2, 0.05),
      "broken_river_r3" = c(-0.05, 45, 30, 0.2, 0.05),
      "campaspe_river_r4" = c(-0.05, 10, 5, 0.2, 0.05),
      "goulburn_river_r4" = c(15, 5, 40, 0.2, 0.05),
      "loddon_river_r4" = c(-0.05, 10, 30, 0.2, 0.05),
      "ovens_river_r5" = c(-0.02, 40, 20, 0.2, 0.05)
    ),
    "cyprinus_carpio" = list(
      "broken_creek_r4" = c(-0.05),
      "broken_river_r3" = c(-0.05),
      "campaspe_river_r4" = c(-0.05),
      "goulburn_river_r4" = c(-0.05),
      "loddon_river_r4" = c(-0.05),
      "ovens_river_r5" = c(-0.05),
      "glenelg_river_r1" = c(-0.05),
      "glenelg_river_r2" = c(-0.05),
      "glenelg_river_r3" = c(-0.05),
      "loddon_river_r2" = c(-0.05),
      "macalister_river_r1" = c(-0.05),
      "mackenzie_river_r3" = c(-0.05),
      "moorabool_river_r3" = c(-0.05),
      "thomson_river_r3" = c(-0.05)
    )
  )
  coefs[[species]][[waterbody]]
}

# function to pull out names of all metrics for relevant species
get_metric_names <- function(species) {
  metric_list <- list(
    "gadopsis_marmoratus" = c(
      "spawning_flow_variability",
      "proportional_spring_flow",
      "proportional_summer_flow",
      "proportional_winter_flow",
      "antecedent_flow",
      "nday_lt5",
      "nday_gt16",
      "nday_lt18",
      "instream_cover",
      "veg_overhang"  # not included in template yet
    ),
    "maccullochella_peelii" = c(
      "spawning_flow_variability",
      "proportional_spring_flow",
      "proportional_max_antecedent",
      "proportional_summer_flow",
      "proportional_winter_flow",
      "spawning_temperature",
      "blackwater_risk",
      "minimum_daily_flow",
      "carrying_capacity"
    ),
    "melanotaenia_fluviatilis" = c(
      "nday_lt10",
      "nday_gt20",
      "redfin",
      "gambusia",
      "instream_cover",
      "spawning_flow_variability",
      "proportional_spring_flow",
      "proportional_summer_flow"
    ),
    "cyprinus_carpio" = c(
      "floodplain_access",
      "flow_variability"
    )
  )
  metric_list[[species]]
}

# function to simulate a pop model for a single species based on input
#    metrics and coefficients
simulate_scenario <- function(x, nsim, init, metrics, coefs, nburnin = 0, ...) {
  
  
  # BF args
  args = list(
    covariates = c(
      format_covariates(metrics[rep(1, nburnin), ]),
      list(
        coefs = coefs[1:5],
        temperature_coefficient = coefs[6],
        coldwater_coefficient = coefs[7]
      )
    )
  )
  
  # MC args
  args = list(
    covariates = c(
      format_covariates(metrics[rep(1, nburnin), ]),
      list(threshold = 0.05),
      list(coefs = coefs)
    ),
    density_dependence = list(
      kdyn = lapply(seq_len(ntime), function(.x) metrics[[i]]$kdyn[.x])
    ),
    density_dependence = list(
      kdyn = lapply(
        rep(1, nburnin), 
        function(.x) metrics$carrying_capacity[.x]
      )
    )
  )
  
  # RB args
  args = list(
    covariates = c(
      format_covariates(metrics[rep(1, nburnin), ]),
      list(
        coefs = coefs[1:3],
        warmwater_coefficient = coefs[4],
        coldwinter_coefficient = coefs[5]
      )
    )
  )
  
  # CC args
  
  # include a short loop to stabilise the initial conditions if 
  #   required
  if (nburnin > 0) {
    
    sims <- simulate(
      x,
      nsim = nsim,
      init = init,
      args = arg_list,
      options = list(
        update = update_binomial_leslie,
        tidy_abundances = floor
      )
    )
    
    # pull out final iteration as the starting point for main sims
    init <- sims[, , nburnin] 
    
  }
  
  # main simulation
  sims <- simulate(
    x,
    nsim = nsim,
    init = init,
    args = list(
      covariates = c(
        format_covariates(metrics),
        list(
          coefs = coefs[1:5],
          temperature_coefficient = coefs[6],
          coldwater_coefficient = coefs[7]
        )
      )
    ),
    options = list(
      update = update_binomial_leslie,
      tidy_abundances = floor
    )
  )
  # return
  sims
  
}

# function to load simulated scenarios and parse scenarios from filenames
load_simulated <- function(type, species) {
  
  # list all saved simulations
  sim_list <- dir("outputs/simulated/")
  
  # filter to target type and species
  sim_list_sub <- sim_list[
    grepl(paste(type, species, sep = "-"), sim_list, ignore.case = TRUE)
  ]
  out <- lapply(sim_list_sub, \(x) qread(paste0("outputs/simulated/", x)))
  
  # parse names
  waterbody <- sapply(strsplit(sim_list_sub, split = "-"), \(x) x[3])
  if (type == "future") {
    
    # work out the scenario as well as the system
    scenario <- sapply(strsplit(sim_list_sub, split = "-"), \(x) x[4])
    scenario <- gsub(".qs", "", scenario)
    
    # split it up by climate and flows, and add waterbody info
    naming_fn <- function(x, wb) {
      tibble(
        waterbody = wb,
        future = x[1],
        future_next = x[2],
        scenario = x[3],
        scenario_next = x[4]
      )
    }
    scenario <- strsplit(scenario, "_") |>
      mapply(FUN = naming_fn, wb = waterbody, SIMPLIFY = FALSE) |>
      bind_rows()
    
  } else {
    
    # just strip the file suffix from hte system name
    waterbody <- gsub(".qs", "", waterbody)
    
    # reformat
    scenario <- tibble(waterbody = waterbody)
    
  }
  
  # and return
  list(
    scenario = scenario,
    sims = out
  )
  
}

fetch_reach_lengths <- function(recompile = FALSE) {
  
  # check if vewh reach info already exists
  vewh_info_exists <- any(grepl("eflow-reaches", dir("data/")))
  if (!vewh_info_exists | recompile) {
    
    # download vewh reach info if not available
    vewh_reaches <- fetch_table("eflow_reaches_20171214", "projects") |>
      collect()
    st_geometry(vewh_reaches) <- st_as_sfc(vewh_reaches$geom, crs = 4283)
    
    # save to file
    qs::qsave(vewh_reaches, file = "data/eflow-reaches.qs")
    
  } else {
    
    # load from file
    vewh_reaches <- qs::qread("data/eflow-reaches.qs")
    
  }
  
  # calculate and return reach lengths  
  vewh_reaches |>
    select(eflowriver, shape_leng, vewh_reach) |>
    mutate(segment_length = st_length(vewh_reaches)) |>
    filter(
      grepl(
        "broken|glenelg|goulburn|mackenzie|campaspe|loddon|ovens|macalister|moorabool|thomson",
        eflowriver,
        ignore.case = TRUE
      )
    ) |>
    group_by(eflowriver, vewh_reach) |>
    summarise(surveyed_length_m = sum(segment_length)) |>
    ungroup() |>
    mutate(surveyed_length_km = surveyed_length_m / 1000)
  
}

plot_target_reaches <- function(recompile = FALSE) {
  
  # check if vewh reach info already exists
  vewh_info_exists <- any(grepl("eflow-reaches", dir("data/")))
  if (!vewh_info_exists | recompile) {
    
    # download vewh reach info if not available
    vewh_reaches <- fetch_table("eflow_reaches_20171214", "projects") |>
      collect()
    st_geometry(vewh_reaches) <- st_as_sfc(vewh_reaches$geom, crs = 4283)
    
    # save to file
    qs::qsave(vewh_reaches, file = "data/eflow-reaches.qs")
    
  } else {
    
    # load from file
    vewh_reaches <- qs::qread("data/eflow-reaches.qs")
    
  }
  
  # repeat for other spatial data sets
  spatial_exists <- any(grepl("watercourse-layer", dir("data/"))) &
    any(grepl("watercourse-sub-layer", dir("data/"))) &
    any(grepl("basin-layer", dir("data/"))) &
    any(grepl("vic-outline-layer", dir("data/")))
  if (!spatial_exists | recompile) {
    
    # download extra spatial layers to plot the target reaches
    watercourse_lm <- fetch_table("hy_watercourse_lm", schema = "spatial") |>
      collect()
    spatial_basin <- fetch_table("basin100", schema = "spatial") |>
      collect()
    vic_outline <- fetch_table("victoria_outline", schema = "spatial") |>
      collect()
    
    # set geometries properly
    st_geometry(watercourse_lm) <- st_as_sfc(watercourse_lm$geom, crs = 4283)
    st_geometry(spatial_basin) <- st_as_sfc(spatial_basin$geom, crs = 4283)
    st_geometry(vic_outline) <- st_as_sfc(vic_outline$geom, crs = 4283)
    
    # filter geofab to Victoria with a 5 km buffer
    vic_poly <- st_polygonize(vic_outline)
    vic_buffered <- vic_poly |>
      st_transform(7899) |>
      st_buffer(dist = 500) |>
      st_transform(4283)
    watercourse_vic <- watercourse_lm |>
      st_filter(vic_buffered, .predicate = st_intersects)
    
    # transform to GDA2020
    spatial_basin <- spatial_basin |> st_transform(7844)
    vic_outline <- vic_outline |> st_transform(7844)
    watercourse_vic <- watercourse_vic |> st_transform(7844)
    
    # create a watercourse layer that includes all major rivers in the target
    #    area
    watercourse_surveyed <- watercourse_vic |>
      filter(
        grepl(
          "broken|glenelg|goulburn|mackenzie|campaspe|loddon|ovens|macalister|moorabool|thomson|wimmera|burnt|kiewa|barwon|la trobe|wannon|pyramid|box|tullaroop",
          name,
          ignore.case = TRUE
        )
      )
    
    # save to file
    qs::qsave(watercourse_vic, file = "data/watercourse-layer.qs")
    qs::qsave(watercourse_surveyed, file = "data/watercourse-sub-layer.qs")
    qs::qsave(spatial_basin, file = "data/basin-layer.qs")
    qs::qsave(vic_outline, file = "data/vic-outline-layer.qs")
    
  } else {
    
    # load from file if it already exists
    watercourse_vic <- qs::qread("data/watercourse-layer.qs")
    watercourse_surveyed <- qs::qread("data/watercourse-sub-layer.qs")
    spatial_basin <- qs::qread("data/basin-layer.qs")
    vic_outline <- qs::qread("data/vic-outline-layer.qs")
    
  }
  
  # filter to target reaches
  .vewh_reach_list <- c(
    "Glenelg River: Reach 1a",
    "Glenelg River: Reach 1b",
    "Glenelg River: Reach 2",
    "Glenelg River: Reach 3",
    "Loddon River: Reach 2",
    "Loddon River: Reach 4a",
    # "Loddon River: Reach 4b", # Twelve Mile Ck
    "Loddon River: Reach 4c",
    "Loddon River: Reach 4d",
    "Macalister River: Reach 1",
    "MacKenzie River: Reach 3",
    # "Moorabool River: Reach 3a",  # West Branch from Lal Lal to confluence with East Branch
    "Moorabool River: Reach 3b",
    "Thomson River: Reach 3",
    "Lower Broken Creek: Reach 3",
    "Broken River: Reach 3",
    "Campaspe River: Reach 4",
    "Goulburn River: Reach 4",
    "Ovens River: Reach 5"
  )
  vewh_reaches <- vewh_reaches |> st_transform(7844) |>
    filter(
      grepl(
        "broken|glenelg|goulburn|mackenzie|campaspe|loddon|ovens|macalister|moorabool|thomson",
        eflowriver,
        ignore.case = TRUE
      )
    )
  
  # clean up reach names
  vewh_reaches <- vewh_reaches |>
    mutate(eflowriver_reach = paste(eflowriver, vewh_reach, sep = ": Reach ")) |>
    filter(eflowriver_reach %in% .vewh_reach_list) |>
    mutate(
      eflowriver_reach = gsub("a$|b$|c$|d$", "", eflowriver_reach),
      eflowriver_reach = gsub(
        "Lower Broken Creek: Reach 3", 
        "Broken Creek: Reach 4",
        eflowriver_reach
      )
    )
  vewh_reaches_mid <- st_centroid(
    vewh_reaches |> 
      filter(!duplicated(vewh_reaches$eflowriver_reach))
  ) 
  
  # order of reach labels (if adjustment required to locations below)  
  # [1] "Glenelg River: Reach 3"    "Glenelg River: Reach 1"    "Loddon River: Reach 4"     "Moorabool River: Reach 3" 
  # [5] "Ovens River: Reach 5"      "Macalister River: Reach 1" "Broken Creek: Reach 4"     "Broken River: Reach 3"    
  # [9] "Glenelg River: Reach 2"    "Campaspe River: Reach 4"   "Goulburn River: Reach 4"   "Loddon River: Reach 2"    
  # [13] "Thomson River: Reach 3"    "MacKenzie River: Reach 3" 
  
  # and plot it all
  vic_outline |>
    ggplot() + 
    geom_sf() + 
    geom_sf(data = spatial_basin) + 
    geom_sf(data = watercourse_surveyed, linewidth = 0.8, col = "gray30") +
    geom_sf(data = vewh_reaches, aes(col = eflowriver_reach), linewidth = 1.15) +
    geom_sf_label(
      data = vewh_reaches_mid, 
      aes(label = eflowriver_reach, col = eflowriver_reach),
      nudge_x = c(1.1, 1.4, -1.2, 0, 1.15, 0.5, 0.5, 0.7, 1.2, -0.5, 0.3, -0.2, 0, 0)[vewh_reaches_mid |> pull(eflowriver_reach) |> order()],
      nudge_y = c(-0.05, 0, 0.2, -0.25, 0, 0.2, 0.25, -0.2, -0.1, -0.2, -0.4, 0.2, -0.2, 0.3)[vewh_reaches_mid |> pull(eflowriver_reach) |> order()],
      size = 3.5,
      fill = scales::alpha("white", 0.95),
      label.size = 0
    ) +
    annotation_scale(location = "bl") +
    xlab("Longitude") +
    ylab("Latitude") +
    annotation_north_arrow(
      location = "bl", 
      which_north = "true", 
      pad_x = unit(0.05, "in"), 
      pad_y = unit(0.25, "in"),
      style = north_arrow_fancy_orienteering
    ) +
    theme(legend.position = "none", text = element_text(size = 8, face = "bold"))
  
}
