# helper functions for validation summaries
# function to estimate CPUE from catch data
estimate_cpue <- function(
    x,
    species, 
    iter = 2000,
    warmup = 1000, 
    chains = 2, 
    cores = 1, 
    recruit = FALSE,
    adult = FALSE,
    use_cached = TRUE
) {
  
  # check species are in correct list
  if (!species %in% unique(x$scientific_name))
    stop("species must be represented in data set", call. = FALSE)
  
  # use saved version if available, otherwise refit the model
  species_clean <- gsub(" ", "_", tolower(species))
  if (recruit | adult) {
    if (recruit) {
      species_clean <- paste0("recruit-", species_clean)
    } else {
      species_clean <- paste0("adult-", species_clean)
    }
  }
  if (use_cached & any(grepl(species_clean, dir("outputs/fitted/")))) {
    
    # if exists, load fitted version
    mod <- qread(paste0("outputs/fitted/cpue-mod-", species_clean, ".qs"))
    
  } else {
    
    # subset data to target species and collate over all gear types and surveys
    #   within a year
    x <- x |>
      filter(scientific_name == species) |>
      group_by(id_site, waterbody, reach_no, survey_year) |>
      summarise(
        catch = sum(catch),
        effort_h = sum(effort_h)
      ) |>
      ungroup()
    
    # add a filter to drop out systems where a species isn't 
    #   detected
    include <- x |>
      group_by(waterbody, reach_no) |>
      summarise(include = sum(catch) > 0) |>
      ungroup()
    x <- x |>
      left_join(include, by = c("waterbody", "reach_no")) |>
      filter(include)
    
    # fit model
    if (recruit) {
      
      mod <- stan_glmer(
        catch ~ (1 | waterbody / reach_no) +
          (1 | id_site) +
          (1 | survey_year) +
          (1 | waterbody:survey_year) +
          offset(effort_h),
        family = poisson,
        data = x,
        iter = iter,
        warmup = warmup,
        chains = chains,
        cores = cores
      )
      
      qsave(mod, file = paste0("outputs/fitted/cpue-mod-", species_clean, ".qs"))
      
    } else {
      
      mod <- stan_glmer(
        catch ~ log_cpue_ym1 +
          (1 | waterbody / reach_no) +
          (1 | id_site) +
          (1 | survey_year) +
          (1 | waterbody:survey_year) +
          offset(effort_h),
        family = poisson,
        data = x |> 
          left_join(
            x |>
              mutate(
                survey_year = survey_year + 1,
                log_cpue_ym1 = log(catch + 1) - log(effort_h)
              ) |>
              select(id_site, survey_year, log_cpue_ym1),
            by = c("id_site", "survey_year")
          ) |>
          filter(!is.na(log_cpue_ym1)),
        iter = iter,
        warmup = warmup,
        chains = chains,
        cores = cores
      )
      
      qsave(mod, file = paste0("outputs/fitted/cpue-mod-", species_clean, ".qs"))
      
    }
    
  }
  
  # return
  mod
  
}  

# function to calculate mid, lower, and upper bounds from simulated
#   trajectories
summarise_sim <- function(x, y, subset, probs, growth_rate = TRUE, zscale = TRUE) {
  
  # pull out abundances for the subset of ages/stages
  abund <- apply(x[, subset, , drop = FALSE], c(1, 3), sum)
  
  # return growth rates if required (raw abundances returned otherwise)
  if (growth_rate) {
    
    # skip if all values are zero
    if (!all(abund == 0)) {
      
      # calculate the pop growth rate first, then summarise this if
      #   z-scaling
      # fill zeros with half of the min observed value
      abund[abund == 0] <- min(abund[abund > 0], na.rm = TRUE) / 2.0
      
      # calculate pop growth rate
      abund <- abund / abund[, c(1L, seq_len(ncol(abund) - 1L))]
      
      # calculate zscores
      if (zscale) {
        abund <- abund - mean(abund)
        abund <- abund / sd(abund)
      }
      
    }
    
    # remove first column (duplicated in previous to avoid errors when
    #    length(abund) == 1, suspect this will still error but not
    #    relevant to this analysis)
    abund <- abund[, -1]
    
  } else {
    
    # just calculate z-scores
    if (zscale) {
      abund <- abund - mean(abund)
      abund <- abund / sd(abund)
    }
    
  }
  
  # collate raw predicted values, dropping first column
  out <- tibble(
    waterbody = y,
    mid = apply(abund, 2, median, na.rm = TRUE),
    lower = apply(abund, 2, quantile, probs = probs[1], na.rm = TRUE),
    upper = apply(abund, 2, quantile, probs = probs[2], na.rm = TRUE)
  )
  
  # return
  out
  
}

# function to calculate pop growth rates from obsered data and compare
#   to simulated values
add_cpue <- function(
    sim,
    cpue_mod, 
    sim_years = 2010:2023,
    probs = c(0.1, 0.9)
) {
  
  #    generate new samples from the fitted posterior for each year/waterbody,
  #    setting previous cpue to 0 to estimate growth rate directly
  #    (no need to divide by catch_ym1)
  newdata <- cpue_mod$data |> 
    distinct(waterbody, reach_no, survey_year) |>
    filter(!is.na(reach_no)) |>
    mutate(
      log_cpue_ym1 = 0,
      effort_h = 1,
      id_site = "abc"
    )
  cpue_pred <- posterior_epred(
    cpue_mod, 
    newdata = newdata,
    re.form = ~ (1 | waterbody / reach_no) +
      (1 | survey_year)
  )
  cpue_ar1 <- tibble(
    newdata,
    cpue = apply(cpue_pred, 2, median),
    lower = apply(cpue_pred, 2, quantile, probs = probs[1]),
    upper = apply(cpue_pred, 2, quantile, probs = probs[2])
  )
  
  # add reach info and rename cpue field
  cpue_ar1 <- cpue_ar1 |>
    mutate(
      waterbody = paste0(
        tolower(gsub(" ", "_", waterbody)),
        "_r",
        reach_no
      )
    ) |>
    select(-reach_no) |>
    rename(growth_rate = cpue)
  
  # add in survey year info to simulated values
  nwaterbody <- sim |> pull(waterbody) |> unique() |> length()
  sim <- sim |> mutate(survey_year = rep(sim_years, nwaterbody))
  
  # z-scale it all
  cpue_std <- cpue_ar1 |>
    group_by(waterbody) |>
    summarise(
      center = mean(growth_rate, na.rm = TRUE),
      width = sd(growth_rate, na.rm = TRUE)
    )
  cpue_ar1 <- cpue_ar1 |>
    left_join(cpue_std, by = "waterbody") |>
    mutate(
      growth_rate_z = (growth_rate - center) / width,
      lower_z = (lower - center) / width,
      upper_z = (upper - center) / width
    ) |>
    select(waterbody, survey_year, growth_rate_z, lower_z, upper_z)
  
  # return this value joined to simulated pop growth rates
  sim |> 
    left_join(cpue_ar1, by = c("waterbody", "survey_year")) |>
    pivot_longer(
      cols = c(mid, growth_rate_z, lower, lower_z, upper, upper_z),
      values_to = "value",
      names_to = "type"
    ) |>
    mutate(
      category = ifelse(grepl("_z", type), "Observed", "Simulated"),
      type = gsub("_z", "", type),
      type = gsub("growth_rate", "mid", type)
    ) |>
    pivot_wider(
      id_cols = c(waterbody, survey_year, category),
      names_from = type,
      values_from = value
    )
  
}

# calculate summary metrics
calculate_val_metrics <- function(
    x, cpue_mod, 
    subset, 
    sim_years,
    sp_idx,
    probs = c(0.1, 0.9), 
    recruit = FALSE
) {
  
  # pull out relevant components of the sims output
  x_noint <- lapply(x, \(x) x[[2]][[sp_idx]])
  x <- lapply(x, \(x) x[[1]][[sp_idx]])
  
  # use functions above to summarise the simulated population trajectories
  x <- mapply(
    summarise_sim, 
    x = x,
    y = names(x),
    MoreArgs = list(
      subset = subset, probs = probs, growth_rate = !recruit
    ),
    SIMPLIFY = FALSE
  )
  x <- bind_rows(x)
  x_noint <- mapply(
    summarise_sim, 
    x = x_noint,
    y = names(x_noint),
    MoreArgs = list(
      subset = subset, probs = probs, growth_rate = !recruit
    ),
    SIMPLIFY = FALSE
  )
  x_noint <- bind_rows(x_noint)
  
  # add estimated CPUE
  x <- add_cpue(
    sim = x,
    cpue_mod = cpue_mod,
    sim_years = sim_years,
    probs = probs
  )
  x_noint <- add_cpue(
    sim = x_noint,
    cpue_mod = cpue_mod,
    sim_years = sim_years,
    probs = probs
  )
  
  # split out the modelled and observed values and calculate residuals
  x <- x |>
    select(waterbody, survey_year, category, mid) |>
    pivot_wider(
      id_cols = c(waterbody, survey_year),
      names_from = category,
      values_from = mid
    ) |>
    mutate(eps = Simulated - Observed)
  x_noint <- x_noint |>
    select(waterbody, survey_year, category, mid) |>
    pivot_wider(
      id_cols = c(waterbody, survey_year),
      names_from = category,
      values_from = mid
    ) |>
    mutate(eps = Simulated - Observed)
  
  # combine both sets of outputs
  x <- bind_rows(
    x |> mutate(type = "interactions"),
    x_noint |> mutate(type = "no interactions")
  )
  
  # calculate all the metrics
  x |>
    group_by(waterbody, type) |>
    summarise(
      r = ifelse(
        !all(is.na(Observed)), 
        cor(Simulated, Observed, use = "complete"),
        NA
      ),
      md = mean(eps, na.rm = TRUE),
      rmse = sqrt(mean(eps ^ 2, na.rm = TRUE)),
      sign = sum(sign(Observed) == sign(Simulated), na.rm = TRUE) / length(Observed)
    )
  
}

# helpers to tidy names in a plot
tidy_names <- function(x) {
  x <- gsub("_", " ", x)
  x <- strsplit(x, split = " ")
  y <- sapply(x, \(.x) gsub("r", "Reach ", .x[3]))
  x <- lapply(x, \(.x) paste0(toupper(substr(.x[1:2], 1, 1)), substr(.x[1:2], 2, nchar(.x))))
  paste(sapply(x, paste, collapse = " "), y, sep = ": ")
}
metric_names <- c(
  "r" = "r",
  "md" = "MD",
  "rmse" = "RMSE",
  "sign" = "Sign"
)

# function to plot validation metrics for one or more rivers at a time
plot_metric <- function(x) {
  
  # prepare data
  x <- x |>
    pivot_longer(
      cols = c(r, md, rmse, sign),
      names_to = "name",
      values_to = "value"
    ) |>
    mutate(
      waterbody = tidy_names(waterbody),
      metric = metric_names[name],
      metric = factor(metric, levels = c("r", "Sign", "RMSE", "MD")),
      species = factor(
        species,
        levels = c(
          "Murray Cod", "Murray Cod (young of year)", "Murray Cod (adults)",
          "River Blackfish", "River Blackfish (young of year)",
          "Murray-Darling Rainbowfish",
          "Common Carp", "Common Carp (young of year)"
        )
      ),
      type_num = ifelse(type == "interactions", 1, 0.25),
      type = factor(
        type,
        levels = c("interactions", "no interactions"),
        labels = c("With interactions", "Without interactions")
      )
    )
  
  # work out level and labels for text
  x <- x |>
    left_join(
      x |> 
        group_by(metric) |> 
        summarise(level = median(value, na.rm = TRUE)),
      by = "metric"
    ) |>
    mutate(label = ifelse(is.na(value), "*", ""))
  
  # set a width based on species
  width_set <- 0.45
  if (any(grepl("rainbowfish", x$species, ignore.case = TRUE)))
    width_set <- -0.675
  
  # create a dummy data set that adds lines at 1/0 for the different facets
  dummy <- tibble(
    metric = factor(c("r", "MD", "RMSE", "Sign"), levels = c("r", "Sign", "RMSE", "MD")),
    height = c(1, 0, 0, 1)
  )
  
  # plot 
  p <- x |>
    ggplot(aes(
      y = value,
      x = waterbody, 
      fill =  species,
      alpha = type
    )) + 
    geom_bar(
      position = position_dodge(width = 0.9, preserve = "single"),
      stat = "identity"
    ) +
    geom_hline(data = dummy, aes(yintercept = height), col = "gray30", linewidth = 1.25) +
    ylab("Value") +
    xlab("Waterbody") +
    facet_wrap( ~ metric, scales = "free") +
    scale_fill_brewer(palette = "Set2", name = "") +
    scale_alpha_manual(values = c(1, 0.4), name = "") +
    ggthemes::theme_hc() +
    theme(
      legend.text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 60, hjust = 1),
      panel.border = element_rect(fill = NA, colour = "gray30", linetype = 1),
      strip.background = element_rect(fill = "white")
    )
  
  # remove legend if just one species
  if (length(unique(x$species)) == 1)
    p <- p + theme(legend.position = "none")
  
  # return
  p
  
}

# river names lookup
.river_lookup <- c(
  "broken_creek_r4" = "Broken Creek (Reach 4)",
  "broken_river_r3" = "Broken River (Reach 3)",
  "campaspe_river_r4" = "Campaspe River (Reach 4)",
  "goulburn_river_r4" = "Goulburn River (Reach 4)",
  "loddon_river_r4" = "Loddon River (Reach 4)",
  "ovens_river_r5" = "Ovens River (Reach 5)",
  "glenelg_river_r1" = "Glenelg River (Reach 1)",
  "glenelg_river_r2" = "Glenelg River (Reach 2)",
  "glenelg_river_r3" = "Glenelg River (Reach 3)",
  "loddon_river_r2" = "Loddon River (Reach 2)",
  "macalister_river_r1" = "Macalister River (Reach 1)",
  "mackenzie_river_r3" = "MacKenzie River (Reach 3)",
  "moorabool_river_r3" = "Moorabool River (Reach 3)",
  "thomson_river_r3" = "Thomson River (Reach 3)"
)

# function to create abundance hindcast plots from simulated and observed data
plot_hindcasts <- function(
    x, 
    cpue, 
    subset, 
    sim_years, 
    sp_idx,
    probs = c(0.1, 0.9), 
    recruit = FALSE
) {
  
  # pull out relevant components of the sims output
  x_noint <- lapply(x, \(x) x[[2]][[sp_idx]])
  x <- lapply(x, \(x) x[[1]][[sp_idx]])
  
  # use functions above to summarise the simulated population trajectories
  x <- mapply(
    summarise_sim, 
    x = x,
    y = names(x),
    MoreArgs = list(
      subset = subset, probs = probs, growth_rate = !recruit
    ),
    SIMPLIFY = FALSE
  )
  x <- bind_rows(x)
  
  # add estimated CPUE
  x <- add_cpue(
    sim = x,
    cpue_mod = cpue,
    sim_years = sim_years,
    probs = probs
  )
  
  # repeat for version without interactions
  x_noint <- mapply(
    summarise_sim, 
    x = x_noint,
    y = names(x_noint),
    MoreArgs = list(
      subset = subset, probs = probs, growth_rate = !recruit
    ),
    SIMPLIFY = FALSE
  )
  x_noint <- bind_rows(x_noint)
  
  # add estimated CPUE
  x_noint <- add_cpue(
    sim = x_noint,
    cpue_mod = cpue,
    sim_years = sim_years,
    probs = probs
  )
  
  # bind both together
  x <- bind_rows(
    x,
    x_noint |>
      filter(category == "Simulated") |>
      mutate(category = "Simulated (no interactions)")
  )
  
  # remove groups without good data
  include <- x |>
    filter(category == "Observed") |>
    group_by(waterbody) |>
    summarise(include = any(!is.na(mid))) |>
    ungroup()
  x <- x |>
    left_join(include, by = "waterbody") |>
    filter(include)
  
  # set up base plot
  p <- x |>
    mutate(waterbody = .river_lookup[waterbody]) |>
    ggplot(aes(x = survey_year, y = mid, col = category, group = category)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_line(position = position_dodge(width = 0.2)) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper), 
      width = 0.2, 
      position = position_dodge(width = 0.2)
    ) +
    scale_color_brewer(
      name = "",
      palette = "Set2"
    ) +
    xlab("Water year") +
    ggthemes::theme_hc() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 8),
      panel.border = element_rect(fill = NA, colour = "gray30", linetype = 1),
      strip.background = element_rect(fill = "white")
    ) + 
    facet_wrap( ~ waterbody, scales = "free")
  
  if (recruit) {
    p <- p + ylab("Scaled recruitment")
  } else {
    p <- p + ylab("Scaled population growth rate")
  }
  
  # and return
  p
  
}

# function to extract trajectories from a simulation object
extract_trajectories <- function(
    scn, x, subset, probs, recruit, sp_idx, zscale, sim_years
) {
  
  # pull out relevant components of the sims output
  x <- lapply(x, \(x) x[[scn]][[sp_idx]])
  
  # use functions above to grab bounds from replicate simulated trajectories
  x <- mapply(
    summarise_sim, 
    x = x,
    y = names(x),
    MoreArgs = list(
      subset = subset, probs = probs, growth_rate = !recruit, zscale = zscale
    ),
    SIMPLIFY = FALSE
  )
  
  # flatten
  x <- bind_rows(x)
  
  # add survey year info
  nwaterbody <- x |> pull(waterbody) |> unique() |> length()
  x <- x |> mutate(survey_year = rep(c(sim_years[1] - 1, sim_years), nwaterbody))
  
  # and return
  x
  
}

# function to plot abundance trajectories from two scenarios
plot_trajectories <- function(
    x, 
    subset,
    sim_years,
    sp_idx,
    probs = c(0.1, 0.9),
    scn = 1:4,
    .names = c("With e-water", "With e-water (no interactions)",
               "Without e-water", "Without e-water (no interactions)")
) {
  
  # extract and summarise sims for all scenarios
  x <- lapply(
    scn,
    extract_trajectories,
    x = x,
    subset = subset, 
    probs = probs, 
    recruit = TRUE,
    zscale = FALSE,
    sp_idx = sp_idx,
    sim_years = sim_years
  )
  
  # and add scenario names  
  x <- mapply(
    \(x, y) x |> mutate(Scenario = y),
    x = x,
    y = .names,
    SIMPLIFY = FALSE
  )
  
  # combine both scenarios into a single object
  x <- do.call(bind_rows, x)
  
  # set up base plot
  p <- x |>
    mutate(waterbody = .river_lookup[waterbody]) |>
    ggplot(aes(x = survey_year, y = mid, col = Scenario, group = Scenario)) +
    geom_point(position = position_dodge(width = 0)) +
    geom_line(position = position_dodge(width = 0)) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper), 
      width = 0.2, 
      position = position_dodge(width = 0)
    ) +
    scale_color_brewer(
      name = "",
      palette = "Set2"
    ) +
    xlab("Water year") +
    ylab("Abundance") +
    ggthemes::theme_hc() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 8),
      panel.border = element_rect(fill = NA, colour = "gray30", linetype = 1),
      strip.background = element_rect(fill = "white")
    ) +
    facet_wrap( ~ waterbody, scales = "free_y")
  
  # and return
  p
  
}

# function to summarise forecasts for all systems for a single species
summarise_forecasts <- function(  
    x, 
    subset,
    sp_idx,
    probs = c(0.1, 0.9),
    scn = 1:2,
    .names = c("With interactions", "Without interactions")
) {
  
  # use functions above to summarise the simulated population trajectories
  x <- lapply(
    scn,
    extract_trajectories,
    x = x,
    subset = subset,
    probs = probs,
    recruit = TRUE,
    zscale = FALSE,
    sp_idx = sp_idx,
    sim_years = 2024:2025
  )
  
  # parse future scenario names into a useful format 
  x <- lapply(
    x,
    \(.x) .x |>
      mutate(
        scenario = sapply(waterbody, \(x) strsplit(x, "-")[[1]][2]),
        waterbody = sapply(waterbody, \(x) strsplit(x, "-")[[1]][1]),
        future = sapply(scenario, \(x) strsplit(x, "_")[[1]][1]),
        future_next = sapply(scenario, \(x) strsplit(x, "_")[[1]][2]),
        scenario_next = sapply(scenario, \(x) strsplit(x, "_")[[1]][4]),
        scenario = sapply(scenario, \(x) strsplit(x, "_")[[1]][3])
      )
  )    
  
  # and add scenario names
  x <- mapply(
    \(x, y) x |> mutate(Scenario = y),
    x = x,
    y = .names,
    SIMPLIFY = FALSE
  )
  
  # combine both scenarios into a single object
  do.call(bind_rows, x)
  
}

# function to plot near-term forecasts from start to final observed year
plot_forecasts <- function(x, system = NULL, target = NULL, climate = NULL) {
  
  # clean up variable values
  x <- x |>
    mutate(
      future = factor(
        future,
        levels = c("dry", "ave", "wet"), 
        labels = c("Dry (2023/2024)", "Ave. (2023/2024)", "Wet (2023/2024)")
      ),
      future_next = factor(
        future_next,
        levels = c("dry", "ave", "wet"),
        labels = c("Dry (2024/2025)", "Ave. (2024/2025)", "Wet (2024/2025)")
      ),
      scenario = factor(
        scenario,
        levels = c("none", "baseflow", "fresh"),
        labels = c("None", "Baseflows", "Freshes")
      ),
      scenario_next = factor(
        scenario_next,
        levels = c("none", "baseflow", "fresh"),
        labels = c("None", "Baseflows", "Freshes")
      )
    )
  
  # set target to latest year if not specified
  if (is.null(target))
    target <- max(x$survey_year)
  
  # and set a flag to work out plots below
  one_step_ahead <- target != max(x$survey_year)
  
  # filter to target water year
  x <- x |> filter(survey_year == target)
  
  # filter to a single system if required
  if (!is.null(system))
    x <- x |> filter(waterbody == system)
  
  # rename waterbodies for plotting
  x <- x |> mutate(waterbody = .river_lookup[waterbody])

  # if showing all climates, need a plot that expands out
  if (is.null(climate)) {
    
    # two options: simpler plot if only showing one step ahead
    if (one_step_ahead) {
      
      p <- x |>
        filter(
          future_next == "Ave. (2024/2025)",
          scenario_next == "None"
        ) |>
        mutate(future = factor(future, labels = c("Dry", "Ave.", "Wet"))) |>
        ggplot(aes(y = mid, x = future, fill = scenario, alpha = Scenario)) +
        geom_bar(position = position_dodge(0.9), stat = "identity") +
        geom_errorbar(
          aes(ymin = lower, ymax = upper),
          position = position_dodge(0.9),
          col = "black",
          width = 0.2
        ) +
        xlab("Climate state") +
        ylab("Abundance") +
        scale_alpha_manual(values = c(1, 0.4), name = "") +
        scale_fill_brewer(name = "Flow priority", palette = "Set2") +
        ggthemes::theme_hc() +
        theme(
          legend.position = "bottom",
          axis.text = element_text(size = 8),
          panel.border = element_rect(fill = NA, colour = "gray30", linetype = 1),
          strip.background = element_rect(fill = "white")
        )
      
      # and add a facet wrap if system is not specified
      if (is.null(system))
        p <- p + facet_wrap( ~ waterbody, scales = "free")
      
    } else {
      
      # must specify system or this plot gets unwieldy
      if (is.null(system))
        stop(
          "system must be specified when plotting multiple ",
          "time steps ahead",
          call. = FALSE
        )
      
      #  plot it
      p <- x |>
        mutate(survey_year = factor(survey_year)) |>
        ggplot(
          aes(y = mid, x = scenario_next, fill = scenario, alpha = Scenario)
        ) +
        geom_bar(position = position_dodge(0.9), stat = "identity") +
        geom_errorbar(
          aes(ymin = lower, ymax = upper),
          position = position_dodge(0.9),
          col = "black",
          width = 0.2
        ) +
        xlab("Flow priority (2024/2025)") +
        ylab("Abundance") +
        scale_alpha_manual(values = c(1, 0.4), name = "") +
        scale_fill_brewer(name = "Flow priority (2023/2024)", palette = "Set2") +
        facet_grid(future_next ~ future) +
        ggthemes::theme_hc() +
        theme(
          legend.position = "bottom",
          axis.text = element_text(size = 8),
          axis.text.x = element_text(angle = 60, hjust = 1),
          panel.border = element_rect(fill = NA, colour = "gray30", linetype = 1),
          strip.background = element_rect(fill = "white")
        )
      
    }
    
  } else {
    
    # show just a single climate for 2023/24 then all for 2024/25
    #  plot it
    p <- x |>
      filter(future == climate) |>
      mutate(survey_year = factor(survey_year)) |>
      ggplot(
        aes(y = mid, x = scenario_next, fill = scenario, alpha = Scenario)
      ) +
      geom_bar(position = position_dodge(0.9), stat = "identity") +
      geom_errorbar(
        aes(ymin = lower, ymax = upper),
        position = position_dodge(0.9),
        col = "black",
        width = 0.2
      ) +
      xlab("Flow priority (2024/2025)") +
      ylab("Abundance") +
      scale_alpha_manual(values = c(1, 0.4), name = "") +
      scale_fill_brewer(name = "Flow priority (2023/2024)", palette = "Set2") +
      facet_grid( ~ future_next) +
      ggthemes::theme_hc() +
      theme(
        legend.position = "bottom",
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.border = element_rect(fill = NA, colour = "gray30", linetype = 1),
        strip.background = element_rect(fill = "white")
      )
    
  }
  
  # and return
  p
  
}
