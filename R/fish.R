# specify species and locations
waterbody_list <- list(
  gadopsis_marmoratus = c(
    "Glenelg River",
    "Loddon River",
    "Mackenzie River",
    "Macalister River",
    "Moorabool River",
    "Thomson River"
  ),
  maccullochella_peelii = c(
    "Broken Creek",
    "Broken River",
    "Campaspe River",
    "Goulburn River", 
    "Loddon River",
    "Ovens River"
  ),
  melanotaenia_fluviatilis = c(
    "Broken Creek",
    "Broken River",
    "Campaspe River",
    "Goulburn River",
    "Loddon River",
    "Ovens River"
  ),
  cyprinus_carpio = c(
    "Broken Creek",
    "Broken River",
    "Campaspe River",
    "Goulburn River",
    "Loddon River",
    "Ovens River",
    "Glenelg River",
    "Mackenzie River",
    "Macalister River",
    "Moorabool River",
    "Thomson River"
  )
)
reach_list <- c(
  "Broken Creek_r4",
  "Broken River_r3",  # use 404224, very little e-water delivered here
  "Campaspe River_r4",
  "Glenelg River_r1", # use R1b gauge flows
  "Glenelg River_r2",
  "Glenelg River_r3",
  "Goulburn River_r4",
  "Loddon River_r2",  # no compliance assessed here, so no e-water recorded
  "Loddon River_r4",
  "Macalister River_r1",
  "Mackenzie River_r3",  # @ McKenzie Ck, too far down but no gauge in R1
  "Moorabool River_r3",
  "Ovens River_r5",  # very little e-water here
  "Thomson River_r3"
)

# grab reach info
.vewh_reach_list <- c(
  "Lower Broken Creek: Reach 3",
  "Broken River: Reach 3",
  "Campaspe River: Reach 4",
  "Glenelg River: Reach 1b",
  "Glenelg River: Reach 2",
  "Glenelg River: Reach 3",
  "Goulburn River: Reach 4",
  "Loddon River: Reach 2",
  "Loddon River: Reach 4c",
  "Macalister River: Reach 1",
  "MacKenzie River: Reach 3",
  "Moorabool River: Reach 3b",
  "Ovens River: Reach 5",
  "Thomson River: Reach 3"
)
.vewh_reach_lengths <- fetch_table("eflow_reaches_20171214", "projects", collect = TRUE)
st_geometry(.vewh_reach_lengths) <- st_as_sfc(.vewh_reach_lengths$geom, crs = 4283)
.vewh_reach_lengths <- .vewh_reach_lengths |>
  mutate(
    reach_length = st_length(geometry),
    waterbody_reach = paste(eflowriver, vewh_reach, sep = ": Reach ")
  ) |>
  filter(waterbody_reach %in% .vewh_reach_list) |>
  select(waterbody_reach, reach_length) |>
  arrange(waterbody_reach) |>
  mutate(
    waterbody = gsub(": Reach ", "_r", waterbody_reach),
    waterbody = gsub(" ", "_", waterbody),
    waterbody = tolower(waterbody),
    waterbody = gsub("lower_broken_creek_r3", "broken_creek_r4", waterbody),
    waterbody = gsub("_r4c", "_r4", waterbody),
    waterbody = gsub("_r3b", "_r3", waterbody),
    waterbody = gsub("_r1b", "_r1", waterbody)
  )

# function to return carrying capacity as a neat tibble
fetch_carrying_capacity <- function() {

  # set carrying capacity
  .carrying_capacity <- list(
    gadopsis_marmoratus = c(
      "glenelg_river_r1" = 50000,
      "glenelg_river_r2" = 100000,
      "glenelg_river_r3" = 200000,
      "loddon_river_r2" = 100000,
      "macalister_river_r1" = 100000,
      "mackenzie_river_r3" = 50000,
      "moorabool_river_r3" = 100000,
      "thomson_river_r3" = 200000
    ),
    maccullochella_peelii = c(
      "broken_creek_r4" = 50000,
      "broken_river_r3" = 200000,
      "campaspe_river_r4" = 50000,
      "goulburn_river_r4" = 100000,
      "loddon_river_r4" = 200000,
      "ovens_river_r5" = 100000
    ),
    melanotaenia_fluviatilis = c(
      "broken_creek_r4" = 1000,
      "broken_river_r3" = 5000,
      "campaspe_river_r4" = 5000,
      "goulburn_river_r4" = 10000,
      "loddon_river_r4" = 3000,
      "ovens_river_r5" = 10000
    ),
    cyprinus_carpio = c(
      "broken_creek_r4" = 500000,
      "broken_river_r3" = 400000,
      "campaspe_river_r4" = 250000,
      "goulburn_river_r4" = 300000,
      "loddon_river_r4" = 500000,
      "ovens_river_r5" = 100000,
      "glenelg_river_r1" = 25000,
      "glenelg_river_r2" = 50000,
      "glenelg_river_r3" = 100000,
      "loddon_river_r2" = 250000,
      "macalister_river_r1" = 500000,
      "mackenzie_river_r3" = 100000,
      "moorabool_river_r3" = 50000,
      "thomson_river_r3" = 500000
    )
  )
  
  # collate and return
  tibble(
    species = rep(names(.carrying_capacity), times = sapply(.carrying_capacity, length)),
    waterbody = unlist(lapply(.carrying_capacity, names)),
    carrying_capacity = unlist(.carrying_capacity)
  )
  
}

# internal function to download fish data from AAEDB
fetch_fish <- function(recruit = FALSE, adult = FALSE, recompile = FALSE) {
  
  # list target species and waterbodies
  .waterbody_list <- c(
    "Broken Creek",
    "Broken River", 
    "Campaspe River",
    "Glenelg River", 
    "Goulburn River",
    "Loddon River", 
    "Macalister River",
    "Moorabool River", 
    "Ovens River",
    "Mackenzie River",
    "MacKenzie River",
    "Thomson River"
  )
  .species_list <- c(
    "Gadopsis marmoratus",
    "Maccullochella peelii",
    "Melanotaenia fluviatilis",
    "Perca fluviatilis",
    "Gambusia holbrooki",
    "Cyprinus carpio"
  )
  
  # check if data exist
  if (recruit | adult) {
    fish_exists <- any(grepl("fish-compiled-adults.qs", dir("data/")))
    if (recruit)
      fish_exists <- any(grepl("fish-compiled-recruits.qs", dir("data/")))
  } else {
    fish_exists <- any(grepl("fish-compiled.qs", dir("data/")))
  }
  
  # if data exist and !recompile, load saved version. Re-extract otherwise
  if (fish_exists & !recompile) {
    
    # load data
    filename <- "fish-compiled"
    if (recruit | adult) {
      if (recruit) {
        filename <- paste0(filename, "-recruits")
      } else {
        filename <- paste0(filename, "-adults")
      }
    }
    cpue <- qread(paste0("data/", filename, ".qs"))
    
  } else {
    
    # grab cpue data from AAEDB, filtered to targets
    if (recruit | adult) {
      
      if (recruit) {
        
      # calculate CPUE for individuals below the estimated recruit length
      #   threshold for Murray cod and blackfish (rainbows are all "recruits")
      cpue_mc <- fetch_cpue(
        c(1, 2, 4, 6, 9, 10:13, 15:50),
        criterion = list(var = "length_cm", lower = 0, upper = 13.5)
      ) |>
        filter(
          scientific_name == "Maccullochella peelii",
          waterbody %in% !!.waterbody_list
        ) |>
        collect()
      cpue_bf <- fetch_cpue(
        c(1, 2, 4, 6, 9, 10:13, 15:50),
        criterion = list(var = "length_cm", lower = 0, upper = 8)
      ) |>
        filter(
          scientific_name == "Gadopsis marmoratus",
          waterbody %in% !!.waterbody_list
        ) |>
        collect()
      cpue_carp <- fetch_cpue(
        c(1, 2, 4, 6, 9, 10:13, 15:50),
        criterion = list(var = "length_cm", lower = 0, upper = 10.0)
      ) |>
        filter(
          scientific_name == "Cyprinus carpio",
          waterbody %in% !!.waterbody_list
        ) |>
        collect()
      
      # combine all three and add info as per full CPUE survey info
      cpue <- bind_rows(cpue_mc, cpue_bf, cpue_carp)
      
      } else {
        
        # grab values for MC, thresholded to adults 
        #    (5+, approx. 400 mm and larger)
        cpue <- fetch_cpue(
          c(1, 2, 4, 6, 9, 10:13, 15:50),
          criterion = list(var = "length_cm", lower = 50.0, upper = Inf)
        ) |>
          filter(
            scientific_name == "Maccullochella peelii",
            waterbody %in% !!.waterbody_list
          ) |>
          collect()
      }
      
    } else {
      
      cpue <- fetch_cpue(c(1, 2, 4, 6, 9, 10:13, 15:50)) |>
        filter(
          scientific_name %in% !!.species_list,
          waterbody %in% !!.waterbody_list
        ) |>
        collect()
      
    }
    
    # add some site info
    site_info <- cpue |> fetch_site_info() |> collect()
    st_geometry(site_info) <- st_as_sfc(site_info$geom_pnt, crs = 4283)
    
    # ignored for now, can use VEWH reach table to add reach info    
    vewh_reaches <- fetch_table("eflow_reaches_20171214", "projects") |>
      collect()
    st_geometry(vewh_reaches) <- st_as_sfc(vewh_reaches$geom, crs = 4283)
    site_info <- site_info |>
      st_join(vewh_reaches |> select(vewh_reach), join = st_within) |>
      mutate(reach_no = ifelse(is.na(reach_no), vewh_reach, reach_no)) |>
      select(-vewh_reach)
    
    # grab a few Ovens sites and duplicate for id_project = 9
    ovens_sub <- site_info |>
      filter(id_site %in% c(3160, 3162, 3163, 3183, 3185)) |>
      mutate(id_project = 9)
    site_info <- bind_rows(site_info, ovens_sub)
    
    # add reach info for unknown reaches and then append to CPUE data
    site_info <- site_info |>
      mutate(
        id_site = as.numeric(id_site),
        reach_no = ifelse(id_site %in% c(4468, 4066, 4068, 4069, 4073), 1, reach_no),
        reach_no = ifelse(id_site %in% c(4067, 4070:4072), 2, reach_no),
        reach_no = ifelse(id_site %in% c(3193), 2, reach_no),
        reach_no = ifelse(id_site %in% c(3194, 4266), 1, reach_no),
        reach_no = ifelse(id_site %in% c(4060, 4061), 1, reach_no),
        reach_no = ifelse(id_site %in% c(4382, 4383), 3, reach_no),
        reach_no = ifelse(id_site %in% c(4109, 4110, 4115), 4, reach_no),
        reach_no = ifelse(id_site %in% c(4116), 5, reach_no),
        reach_no = ifelse(id_site %in% c(3133), 0, reach_no),
        reach_no = ifelse(id_site %in% c(3134), 5, reach_no),
        reach_no = ifelse(id_site %in% c(1643:1646, 3164, 4225, 4229, 4231), 0, reach_no),
        reach_no = ifelse(id_site %in% c(3160:3163, 3182, 4172:4185, 4194:4197, 4199:4204, 4208:4212, 4217:4224, 4228, 4232:4241), 4, reach_no),
        reach_no = ifelse(id_site %in% c(3322, 3324, 3325), 2, reach_no),
        reach_no = ifelse(id_site %in% c(3167), 3, reach_no),
        reach_no = ifelse(id_site %in% c(3112, 3177:3180, 4451), 5, reach_no),
        reach_no = ifelse(id_site %in% c(3113, 3181), 4, reach_no),
        reach_no = ifelse(id_site %in% c(3171:3176), 6, reach_no),
        reach_no = ifelse(waterbody == "Broken Creek" & is.na(reach_no), 5, reach_no),
        reach_no = ifelse(waterbody == "Ovens River" & is.na(reach_no), 5, reach_no)
      )
    
    # add reach and lat/long info, filter to remove Buffalo sites and add reach info
    #   for Ovens sites in id_project 9
    cpue <- cpue |>
      left_join(
        site_info |>
          select(id_project, id_site, reach_no, latitude, longitude),
        by = c("id_project", "id_site")
      ) |>
      filter(!(id_site %in% c(4226, 4227, 4230))) |>
      mutate(waterbody = ifelse(id_site %in% c(3134), "Thomson River", waterbody))
    
    # filter out some upper reaches of the Moorabool and Macalister (0, 2, 0)
    cpue <- cpue |>
      filter(
        !(waterbody == "Macalister River" & reach_no == 0),
        !(waterbody == "Moorabool River" & reach_no %in% c(0, 2))
      )
    
    # filter out sites without geom information
    cpue <- cpue |>
      filter(
        id_site %in% !!(site_info |> filter(!is.na(geom_pnt)) |> pull(id_site) |> unique())
      )
    
    # save this
    filename <- "fish-compiled"
    if (recruit | adult) {
      filename <- paste0(filename, "-adults")
      if (recruit)
        filename <- paste0(filename, "-recruits")
    }
    qsave(cpue, file = paste0("data/", filename, ".qs"))
    
  }
  
  # return
  cpue
  
}
