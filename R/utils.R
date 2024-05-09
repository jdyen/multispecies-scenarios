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
  if (species == "maccullochella_peelii") {
    adults <- 5:30
    mat <- murray_cod(k = n)$dynamics$matrix
    initial_age_frequency <- Re(eigen(mat)$vectors[, 1])
  }
  if (species == "cyprinus_carpio") {
    adults <- 3:28
    mat <- common_carp(k = n)$dynamics$matrix
    initial_age_frequency <- Re(eigen(mat)$vectors[, 1])
  }
  if (species == "melanotaenia_fluviatilis") {
    adults <- 2:5
    mat <- murray_rainbowfish(k = n)$dynamics$matrix
    initial_age_frequency <- Re(eigen(mat)$vectors[, 1])
  }
  if (species == "gadopsis_marmoratus") {
    adults <- 2:11
    mat <- river_blackfish(k = n)$dynamics$matrix
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
prepare_inits <- function(waterbody, cpue, metrics, pops, nsim, ...) {
  
  # extract carrying capacity for species and reach
  metrics_sub <- metrics |> filter(waterbody == !!waterbody)
  k_mc <- metrics_sub |> 
    filter(species == "maccullochella_peelii", waterbody == !!waterbody) |>
    pull(carrying_capacity) |>
    unique()
  k_cc <- metrics_sub |> 
    filter(species == "cyprinus_carpio", waterbody == !!waterbody) |>
    pull(carrying_capacity) |>
    unique()
  k_rb <- metrics_sub |> 
    filter(species == "melanotaenia_fluviatilis", waterbody == !!waterbody) |>
    pull(carrying_capacity) |>
    unique()
  k_bf <- metrics_sub |> 
    filter(species == "gadopsis_marmoratus", waterbody == !!waterbody) |>
    pull(carrying_capacity) |>
    unique()
  
  # pull out start year and carrying capacities
  start <- min(metrics$water_year)
  
  # pull out CPUE for the target species and add max values
  cpue_sub <- cpue |> 
    mutate(
      species = tolower(gsub(" ", "_", scientific_name)),
      waterbody = paste(
        tolower(gsub(" ", "_", waterbody)),
        reach_no,
        sep = "_r"
      )
    ) |>
    filter(waterbody == !!waterbody)
  cpue_max <- cpue_sub |>
    group_by(species) |>
    summarise(cpue_max = max(catch / effort_h)) |>
    ungroup()
  cpue_sub <- cpue_sub |>
    left_join(cpue_max, by = "species")
    
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
  cpue_mc <- cpue_sub |>
    filter(species == "maccullochella_peelii") |>
    group_by(survey_year) |>
    summarise(
      cpue_max = unique(cpue_max),
      catch = sum(catch),
      effort_h = sum(effort_h)
    ) |>
    mutate(
      cpue = catch / effort_h
    )
  cpue_cc <- cpue_sub |>
    filter(species == "cyprinus_carpio") |>
    group_by(survey_year) |>
    summarise(
      cpue_max = unique(cpue_max),
      catch = sum(catch),
      effort_h = sum(effort_h)
    ) |>
    mutate(
      cpue = catch / effort_h
    )
  cpue_rb <- cpue_sub |>
    filter(species == "melanotaenia_fluviatilis") |>
    group_by(survey_year) |>
    summarise(
      cpue_max = unique(cpue_max),
      catch = sum(catch),
      effort_h = sum(effort_h)
    ) |>
    mutate(
      cpue = catch / effort_h
    )
  cpue_bf <- cpue_sub |>
    filter(species == "gadopsis_marmoratus") |>
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
    filter(waterbody == !!waterbody) |>
    pull(reach_length) |>
    as.numeric()
  
  # work out the rescaling rate for catch to abundance based on a capture
  #    rate of 40%
  adult_rescale <- (1 / 0.4) * (reach_length / 100)
  init_mc <- init_cc <- init_rb <- init_bf <- NA
  if (length(k_mc) > 0) {
    init_mc <- set_initial(
      species = "maccullochella_peelii",
      cpue = cpue_mc$cpue,
      cpue_max = cpue_mc$cpue_max,
      effort_h = cpue_mc$effort_h,
      n = k_mc,
      nsim = nsim,
      rescale = adult_rescale
    )
  }
  if (length(k_cc) > 0) {
    init_cc <- set_initial(
      species = "cyprinus_carpio",
      cpue = cpue_cc$cpue,
      cpue_max = cpue_cc$cpue_max,
      effort_h = cpue_cc$effort_h,
      n = k_cc,
      nsim = nsim,
      rescale = adult_rescale
    )
  }
  if (length(k_rb) > 0) {
    init_rb <- set_initial(
      species = "melanotaenia_fluviatilis",
      cpue = cpue_rb$cpue,
      cpue_max = cpue_rb$cpue_max,
      effort_h = cpue_rb$effort_h,
      n = k_rb,
      nsim = nsim,
      rescale = adult_rescale
    )
  }
  if (length(k_bf) > 0) {
    init_bf <- set_initial(
      species = "gadopsis_marmoratus",
      cpue = cpue_bf$cpue,
      cpue_max = cpue_bf$cpue_max,
      effort_h = cpue_bf$effort_h,
      n = 0.2 * k_bf,
      nsim = nsim,
      rescale = adult_rescale
    )
  }
  
  # re-order based on inpts
  nstage <- sapply(pops, \(x) nrow(x$dynamics$matrix))
  stage_idx <- match(nstage, c(50, 28, 7, 11))
  stage_idx <- stage_idx[!is.na(stage_idx)]
  
  # return just those that match pops
  inits <- list(init_mc, init_cc, init_rb, init_bf)[stage_idx]
  inits[sapply(inits, length) > 0]
  
}

# function to prepare population model objects for each species
prepare_pop <- function(metrics, waterbody, stocking) {
  
  # extract carrying capacity for species and reach
  k_mc <- metrics |> 
    filter(species == "maccullochella_peelii", waterbody == !!waterbody) |>
    pull(carrying_capacity) |>
    unique()
  k_cc <- metrics |> 
    filter(species == "cyprinus_carpio", waterbody == !!waterbody) |>
    pull(carrying_capacity) |>
    unique()
  k_rb <- metrics |> 
    filter(species == "melanotaenia_fluviatilis", waterbody == !!waterbody) |>
    pull(carrying_capacity) |>
    unique()
  k_bf <- metrics |> 
    filter(species == "gadopsis_marmoratus", waterbody == !!waterbody) |>
    pull(carrying_capacity) |>
    unique()
  
  # grab stocking info if required (default to zero, otherwise)
  n_stocked <- rep(0, nrow(metrics))
  system_lu <- c(
    "broken_creek_r4" = "Broken Creek",
    "broken_river_r3" = "Broken River",
    "campaspe_river_r4" = "Campaspe River",
    "goulburn_river_r4" = "Goulburn River",
    "loddon_river_r4" = "Loddon River",
    "ovens_river_r5" = "Ovens River"
  )
  n_stocked <- metrics |>
    filter(species == "maccullochella_peelii") |>
    select(water_year) |>
    left_join(
      stocking |>
        filter(
          Species == "Murray Cod",
          System == system_lu[waterbody]
        ) |>
        select(Year, Number) |>
        mutate(Year = Year + 1) |>
        rename(water_year = Year, number_stocked = Number),
      by = c("water_year")
    ) |>
    mutate(number_stocked = ifelse(is.na(number_stocked), 0, number_stocked)) |>
    pull(number_stocked)
  
  # prepare populations and return
  mc <- specify_pop_model(
    species = "maccullochella_peelii",
    waterbody = waterbody,
    ntime = sum(metrics$species == "maccullochella_peelii"), 
    nstocked = n_stocked,
    k = k_mc
  )
  cc <- specify_pop_model(
    species = "cyprinus_carpio",
    waterbody = waterbody,
    ntime = sum(metrics$species == "cyprinus_carpio"), 
    k = k_cc
  )
  rb <- specify_pop_model(
    species = "melanotaenia_fluviatilis",
    waterbody = waterbody,
    ntime = sum(metrics$species == "melanotaenia_fluviatilis"), 
    k = k_rb
  )
  bf <- specify_pop_model(
    species = "gadopsis_marmoratus",
    waterbody = waterbody,
    ntime = sum(metrics$species == "gadopsis_marmoratus"), 
    k = k_bf
  )
  
  # return
  list(
    mc = mc,
    cc = cc,
    rb = rb,
    bf = bf
  )
  
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

# function to specify species interactions
specify_interactions <- function(pops, ...) {
  
  # check which type of river is being modelling
  sp_list <- sapply(pops, \(x) x$species)
  system <- "unspecified"
  if (all(c("Murray cod", "Common carp", "Murray rainbowfish") %in% sp_list))
    system <- "northern"
  if (all(c("Common carp", "River blackfish") %in% sp_list))
    system <- "coastal"
  if (system == "unspecified")
    stop("species lists must include complete species sets", call. = FALSE)
  
  # set up MC/CC/RB interactions for the northern system
  if (system == "northern") {
    
    # pull out target pop dynamics objects
    mc <- pops[[which(sp_list == "Murray cod")]]
    cc <- pops[[which(sp_list == "Common carp")]]
    rb <- pops[[which(sp_list == "Murray rainbowfish")]]
    
    # carp negatively affect smaller MC
    mask_lt4 <- transition(mc$dynamics$matrix, dim = 1:4)
    fun_lt4 <- function(x, n, interacting = TRUE) {
      # n is the population vector of carp
      out <- x
      if (interacting)
        out <- x * exp(-sum(n[3:28]) / 100000)
      out
    }
    
    # carp benefit from smaller MC
    mask_all <- transition(cc$dynamics$matrix)
    fun_all <- function(x, n, interacting = TRUE) {
      # n is the population vector of MC
      out <- x
      if (interacting)
        out <- x / (1 + x * sum(n[1:4]) / 100000)
      out
    }
    
    # carp reduce MC recruitment
    mask_rec <- reproduction(mc$dynamics$matrix)
    fun_rec <- function(x, n, interacting = TRUE) {
      # n is the population vector of carp
      out <- x
      if (interacting)
        out <- x * exp(-sum(n[3:28]) / 100000)
      out
    }
    
    # carp reduce RB abundance
    mask_rb_cc <- transition(rb$dynamics$matrix)
    fun_rb_cc <- function(x, n, interacting = TRUE) {
      # n is the population vector of carp
      out <- x
      if (interacting)
        out <- x * exp(-sum(n[3:28]) / 100000)
      out
    }
    
    # MC reduce RB abundance
    mask_rb_mc <- transition(rb$dynamics$matrix)
    fun_rb_mc <- function(x, n, interacting = TRUE) {
      # n is the population vector of cod
      out <- x
      if (interacting)
        out <- x * exp(-sum(n[3:50]) / 20000)
      out
    }
    
    # big MC reduce survival of small carp
    mask_lt5 <- transition(cc$dynamics$matrix, dim = 1:5)
    fun_lt5 <- function(x, n, interacting = TRUE) {
      # n is the population vector of MC
      out <- x
      if (interacting)
        out <- x * exp(-sum(n[3:50]) / 20000)
      out
    }
    
    # combine masks and functions into pairwise_interaction objects
    mc_lt4 <- pairwise_interaction(mc$dynamics, cc$dynamics, mask_lt4, fun_lt4)
    cc_all <- pairwise_interaction(cc$dynamics, mc$dynamics, mask_all, fun_all)
    mc_rec <- pairwise_interaction(mc$dynamics, cc$dynamics, mask_rec, fun_rec)
    rb_cc <- pairwise_interaction(rb$dynamics, cc$dynamics, mask_rb_cc, fun_rb_cc)
    rb_mc <- pairwise_interaction(rb$dynamics, mc$dynamics, mask_rb_mc, fun_rb_mc)
    cc_lt5 <- pairwise_interaction(cc$dynamics, mc$dynamics, mask_lt5, fun_lt5)
    
    # collate interactions
    interactions <- list(mc_lt4, cc_all, mc_rec, rb_cc, rb_mc, cc_lt5)
    
  } else {
    
    # pull out target dynamics objects
    cc <- pops[[which(sp_list == "Common carp")]]
    bf <- pops[[which(sp_list == "River blackfish")]]
    
    # carp negatively affect smaller BF
    mask_lt4 <- transition(bf$dynamics$matrix, dim = 1:4)
    fun_lt4 <- function(x, n, interacting = TRUE) {
      # n is the population vector of carp
      out <- x
      if (interacting)
        out <- x * exp(-sum(n[3:28]) / 100000)
      out
    }
    
    # carp benefit from smaller BF
    mask_all <- transition(cc$dynamics$matrix)
    fun_all <- function(x, n, interacting = TRUE) {
      # n is the population vector of MC
      out <- x
      if (interacting)
        out <- x / (1 + x * sum(n[1:4]) / 100000)
      out
    }
    
    # carp reduce BF recruitment
    mask_rec <- reproduction(bf$dynamics$matrix)
    fun_rec <- function(x, n, interacting = TRUE) {
      # n is the population vector of carp
      out <- x
      if (interacting)
        out <- x * exp(-sum(n[3:28]) / 100000)
      out
    }
    
    # big BF reduce survival of small carp
    mask_lt5 <- transition(cc$dynamics$matrix, dim = 1:5)
    fun_lt5 <- function(x, n, interacting = TRUE) {
      # n is the population vector of MC
      out <- x
      if (interacting)
        out <- x * exp(-sum(n[5:11]) / 2000)
      out
    }
    
    # combine masks and functions into pairwise_interaction objects
    bf_lt4 <- pairwise_interaction(bf$dynamics, cc$dynamics, mask_lt4, fun_lt4)
    cc_all <- pairwise_interaction(cc$dynamics, bf$dynamics, mask_all, fun_all)
    bf_rec <- pairwise_interaction(bf$dynamics, cc$dynamics, mask_rec, fun_rec)
    cc_lt5 <- pairwise_interaction(cc$dynamics, bf$dynamics, mask_lt5, fun_lt5)
    
    # collate interactions
    interactions <- list(bf_lt4, cc_all, bf_rec, cc_lt5)
    
  }
  
  # return
  interactions
  
}

# function to line up ms and single species objects
match_pops <- function(ms, pops) {
  
  # pull out unique ID for each pop dyanmics object
  ms_hex <- sapply(ms$dynamics, \(x) x$hex)
  pop_hex <- sapply(pops, \(x) x$dynamics$hex)
  
  # check that all are in there
  if (!all(ms_hex %in% pop_hex)) {
    stop(
      "multispecies object does not contain all dynamics objects", 
      call. = FALSE
    )
  }
  
  # work out species and return order in the ms object
  sp_names <- sapply(pops, \(x) x$species)
  sp_names[match(pop_hex, ms_hex)]
  
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

# function to prepare arguments and metrics for multispecies simulations
prepare_args <- function(metrics, waterbody, pops, nburnin = NULL) {
  
  # extract carrying capacity for species and reach
  metrics_sub <- metrics |> filter(waterbody == !!waterbody)
  
  # add or rename covariates
  metrics_sub <- metrics_sub |>
    mutate(
      flow_variability = spawning_flow_variability,
      floodplain_access = ifelse(proportional_spring_flow > 5, 1, 0),
      antecedent_flow = proportional_antecedent_flow,
      blackwater_risk = hypoxia_risk,
      redfin = presence_redfin, 
      gambusia = presence_gambusia
    )
  
  # filter metrics to each species
  metrics_sub_mc <- metrics_sub |>
    filter(species == "maccullochella_peelii") 
  metrics_sub_cc <- metrics_sub |>
    filter(species == "cyprinus_carpio") 
  metrics_sub_rb <- metrics_sub |>
    filter(species == "melanotaenia_fluviatilis") 
  metrics_sub_bf <- metrics_sub |>
    filter(species == "gadopsis_marmoratus") 
  
  # match pops with correct args
  nstage <- sapply(pops, \(x) nrow(x$dynamics$matrix))
  pop_idx <- match(c(50, 28, 7, 11), nstage)
  stage_idx <- match(nstage, c(50, 28, 7, 11))
  stage_idx <- stage_idx[!is.na(stage_idx)]
  
  # and collate appropriate args object for burnin or full scenarios
  if (is.null(nburnin)) {
    
    mc_args <- list(
      density_dependence_n = pops[[pop_idx[1]]]$arguments$density_dependence_n,
      covariates = c(
        format_covariates(
          metrics_sub_mc |>
            select(all_of(get_metric_names("maccullochella_peelii")))
        ),
        list(threshold = 0.05),
        list(coefs = get_coefs("maccullochella_peelii", waterbody))
      ),
      density_dependence = list(
        kdyn = lapply(
          seq_len(nrow(metrics_sub_mc)),
          \(.x) metrics_sub_mc$carrying_capacity[.x]
        )
      )
    )
    cc_args <- list(
      density_dependence_n = pops[[pop_idx[2]]]$arguments$density_dependence_n,
      covariates = c(
        format_covariates(
          metrics_sub_cc |>
            select(all_of(get_metric_names("cyprinus_carpio")))
        ),
        list(coefs = get_coefs("cyprinus_carpio", waterbody))
      )
    )
    rb_args <- list(
      covariates = c(
        format_covariates(
          metrics_sub_rb |>
            select(all_of(get_metric_names("melanotaenia_fluviatilis")))
        ),
        list(coefs = get_coefs("melanotaenia_fluviatilis", waterbody)[1:3]),
        list(warmwater_coefficient = get_coefs("melanotaenia_fluviatilis", waterbody)[4]),
        list(coldwinter_coefficient = get_coefs("melanotaenia_fluviatilis", waterbody)[5])
      )
    )
    bf_args <- list(
      covariates = c(
        format_covariates(
          metrics_sub_bf |>
            select(all_of(get_metric_names("gadopsis_marmoratus")))
        ),
        list(coefs = get_coefs("gadopsis_marmoratus", waterbody)[1:5]),
        list(temperature_coefficient = get_coefs("gadopsis_marmoratus", waterbody)[6]),
        list(coldwater_coefficient = get_coefs("gadopsis_marmoratus", waterbody)[7])
      )
    )
    
  } else {
    
    mtc_mc <- metrics_sub_mc |>
      select(all_of(get_metric_names("maccullochella_peelii")))
    mtc_cc <- metrics_sub_cc |>
      select(all_of(get_metric_names("cyprinus_carpio")))
    mtc_rb <- metrics_sub_rb |>
      select(all_of(get_metric_names("melanotaenia_fluviatilis")))
    mtc_bf <- metrics_sub_bf |>
      select(all_of(get_metric_names("gadopsis_marmoratus")))
    mc_args <- list(
      density_dependence_n = pops[[pop_idx[1]]]$arguments$density_dependence_n,
      covariates = c(
        format_covariates(mtc_mc[rep(1, nburnin), ]),
        list(threshold = 0.05),
        list(coefs = get_coefs("maccullochella_peelii", waterbody))
      ),
      density_dependence = list(
        kdyn = lapply(
          seq_len(nburnin),
          \(.x) metrics_sub_mc$carrying_capacity[1]
        )
      )
    )
    cc_args <- list(
      density_dependence_n = pops[[pop_idx[2]]]$arguments$density_dependence_n,
      covariates = c(
        format_covariates(mtc_cc[rep(1, nburnin), ]),
        list(coefs = get_coefs("cyprinus_carpio", waterbody))
      )
    )
    rb_args <- list(
      covariates = c(
        format_covariates(mtc_rb[rep(1, nburnin), ]),
        list(coefs = get_coefs("melanotaenia_fluviatilis", waterbody)[1:3]),
        list(warmwater_coefficient = get_coefs("melanotaenia_fluviatilis", waterbody)[4]),
        list(coldwinter_coefficient = get_coefs("melanotaenia_fluviatilis", waterbody)[5])
      )
    )
    bf_args <- list(
      covariates = c(
        format_covariates(mtc_bf[rep(1, nburnin), ]),
        list(coefs = get_coefs("gadopsis_marmoratus", waterbody)[1:5]),
        list(temperature_coefficient = get_coefs("gadopsis_marmoratus", waterbody)[6]),
        list(coldwater_coefficient = get_coefs("gadopsis_marmoratus", waterbody)[7])
      )
    )
    
  }
  
  # match to the models in pops and return in same order
  args <- list(mc_args, cc_args, rb_args, bf_args)[stage_idx]
  args[sapply(args, length) > 0]
  
}

# function to simulate a pop model for multiple species based on input
#    metrics and coefficients
simulate_scenario <- function(
    x, nsim, init, args, args2 = NULL, args_burnin = NULL, ...
) {
  
  # how many time steps?  
  ntime <- length(args[[1]]$covariates[[1]])
  
  # include a short loop to stabilise the initial conditions if required
  if (!is.null(args_burnin)) {
    
    # simulate nburnin times and store the final interation to initialise
    #   the full model run
    nburnin <- length(args_burnin[[1]]$covariates[[1]])
    sims <- simulate(
      x,
      nsim = nsim,
      args = args_burnin,
      init = inits,
      options = list(
        ntime = nburnin,
        update = update_binomial_leslie,
        tidy_abundances = floor
      )
    )
    
    # pull out final iteration as the starting point for main sims
    init <- lapply(sims, \(x) x[, , nburnin])
    
  }
  
  # main simulation
  sims <- simulate(
    x,
    nsim = nsim,
    init = init,
    args = args,
    options = list(
      update = update_binomial_leslie,
      tidy_abundances = floor
    )
  )
  
  # repeat but without interactions
  args_noint <- args
  for (i in seq_along(args))
    args_noint[[i]]$interactions <- list(interacting = FALSE)
  sims_noint <- simulate(
    x,
    nsim = nsim,
    init = init,
    args = args_noint,
    options = list(
      update = update_binomial_leslie,
      tidy_abundances = floor
    )
  )
  
  # collate two simulation outputs
  sims <- list(sims, sims_noint)
  
  # secondary simulation if required
  if (!is.null(args2)) {
    sims2 <- simulate(
      x,
      nsim = nsim,
      init = init,
      args = args2,
      options = list(
        update = update_binomial_leslie,
        tidy_abundances = floor
      )
    )
    
    # repeat but without interactions
    args2_noint <- args2
    for (i in seq_along(args2))
      args2_noint[[i]]$interactions <- list(interacting = FALSE)
    sims2_noint <- simulate(
      x,
      nsim = nsim,
      init = init,
      args = args2_noint,
      options = list(
        update = update_binomial_leslie,
        tidy_abundances = floor
      )
    )
    
    # collate all
    sims <- c(sims, list(sims2, sims2_noint))
    
  }
  
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
