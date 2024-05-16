#' @description
#'   function to specify an interaction between two species, either as a
#'   positive or negative effect of different strengths and for different
#'   stages of the source species
#' 
#' @param k carrying capacity of source species
#' @param theta strength of interaction. For negative terms, 0.7, 0.22, and
#'    0.05 give roughly 50%, 20%, and 5% reductions, respectively. For
#'    positive terms, 5, 0.85, and 0.2 give roughly 50%, 20%, and 5% increases
#'    at maximum values of the source
#' @param stage stages of the source species that influence the target
#' @param negative logical to specify whether the interaction is negative or
#'    positive in nature
interaction_fn <- function(
    k, theta, stage, negative = TRUE, cut = TRUE, ...
) {
  
  # evaluate pars so they don't get lost
  force(k)
  force(theta)
  force(stage)
  force(cut)
  
  # define a function with these parameters
  if (negative) {
    fn <- function(x, n, interacting = TRUE, ...) {
      
      out <- x
      
      if (interacting)
        out <- x * exp(-(theta * sum(n[stage])) / k)

      # truncate values to [0, 1] if cut == TRUE
      if (cut) {
        eps <- 1e-5
        out[out < 0] <- eps
        out[out > 1] <- 1 - eps
      }
      
      out
      
    }
  } else {
    fn <- function(x, n, interacting = TRUE, ...) {
      
      out <- x
      
      if (interacting)
        out <- x * (0.5 + (1 / (1 + exp(- (theta * sum(n[stage])) / k))))
      
      # truncate values to [0, 1] if cut == TRUE
      if (cut) {
        eps <- 1e-5
        out[out < 0] <- eps
        out[out > 1] <- 1 - eps
      }
      
      out
    }
  }
  
  # return
  fn
  
}

# function to specify species interactions
#' @param pops list of population dynamics objects for which interactions are
#'    required
#' @param scale scaling factor applied to maximum densities. Defaults to 1
#'    (no change), with decreases in this value increasing the strength of
#'    pairwise interactions (and vice versa)
specify_interactions <- function(pops, scale = 1, ...) {
  
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
    
    # specify age classes for each size class
    small_fish <- 1
    medium_fish <- 1:4
    large_mc <- 5:50
    large_carp <- 5:28
    rb_stages <- 1:7
    kmc <- scale * 20000
    kcarp <- scale * 10000
    krb <- scale * 1000
    
    # specify interactions one-by-one (reciprocating all unless explicitly not
    #    required)    
    
    #   - Strong negative effect of large carp on small MC
    mask_largecarp_eat_smallcod <- reproduction(mc$dynamics$matrix)
    fun_largecarp_eat_smallcod <- interaction_fn(
      k = kcarp, theta = 0.7, stage = large_carp, negative = TRUE, cut = FALSE
    ) 
    largecarp_eat_smallcod <- pairwise_interaction(
      target = mc$dynamics,
      source = cc$dynamics, 
      mask_largecarp_eat_smallcod, 
      fun_largecarp_eat_smallcod
    )
    
    #   - Strong positive effect of small MC on large carp
    mask_smallcod_feed_largecarp <- transition(
      cc$dynamics$matrix, dim = large_carp
    )
    fun_smallcod_feed_largecarp <- interaction_fn(
      k = kmc, theta = 5, stage = small_fish, negative = FALSE
    ) 
    smallcod_feed_largecarp <- pairwise_interaction(
      target = cc$dynamics,
      source = mc$dynamics, 
      mask_smallcod_feed_largecarp, 
      fun_smallcod_feed_largecarp
    )
    
    #   - Moderate negative effect of large MC on small carp
    mask_largecod_eat_smallcarp <- reproduction(cc$dynamics$matrix)
    fun_largecod_eat_smallcarp <- interaction_fn(
      k = kmc, theta = 0.22, stage = large_mc, negative = TRUE, cut = FALSE
    ) 
    largecod_eat_smallcarp <- pairwise_interaction(
      target = cc$dynamics,
      source = mc$dynamics, 
      mask_largecod_eat_smallcarp, 
      fun_largecod_eat_smallcarp
    )
    
    #   - Moderate positive effect of small carp on large MC
    mask_smallcarp_feed_largecod <- transition(
      mc$dynamics$matrix, dim = large_mc
    )
    fun_smallcarp_feed_largecod <- interaction_fn(
      k = kcarp, theta = 0.85, stage = small_fish, negative = FALSE
    ) 
    smallcarp_feed_largecod <- pairwise_interaction(
      target = mc$dynamics,
      source = cc$dynamics, 
      mask_smallcarp_feed_largecod, 
      fun_smallcarp_feed_largecod
    )
    
    #   - Weak negative effect of medium carp on small MC
    mask_medcarp_eat_smallcod <- reproduction(mc$dynamics$matrix)
    fun_medcarp_eat_smallcod <- interaction_fn(
      k = kcarp, theta = 0.05, stage = medium_fish, negative = TRUE, cut = FALSE
    ) 
    medcarp_eat_smallcod <- pairwise_interaction(
      target = mc$dynamics,
      source = cc$dynamics, 
      mask_medcarp_eat_smallcod, 
      fun_medcarp_eat_smallcod
    )
    
    #   - Weak positive effect of small MC on medium carp
    mask_smallcod_feed_medcarp <- transition(
      cc$dynamics$matrix, dim = medium_fish
    )
    fun_smallcod_feed_medcarp <- interaction_fn(
      k = kmc, theta = 0.2, stage = small_fish, negative = FALSE
    ) 
    smallcod_feed_medcarp <- pairwise_interaction(
      target = cc$dynamics,
      source = mc$dynamics, 
      mask_smallcod_feed_medcarp, 
      fun_smallcod_feed_medcarp
    )
    
    #   - Weak negative effect of medium MC on small carp
    mask_medcod_eat_smallcarp <- reproduction(cc$dynamics$matrix)
    fun_medcod_eat_smallcarp <- interaction_fn(
      k = kmc, theta = 0.05, stage = medium_fish, negative = TRUE, cut = FALSE
    ) 
    medcod_eat_smallcarp <- pairwise_interaction(
      target = mc$dynamics,
      source = cc$dynamics, 
      mask_medcod_eat_smallcarp, 
      fun_medcod_eat_smallcarp
    )
    
    #   - Weak positive effect of small carp on medium MC
    mask_smallcarp_feed_medcod <- transition(
      mc$dynamics$matrix, dim = medium_fish
    )
    fun_smallcarp_feed_medcod <- interaction_fn(
      k = kcarp, theta = 0.2, stage = small_fish, negative = FALSE
    ) 
    smallcarp_feed_medcod <- pairwise_interaction(
      target = mc$dynamics,
      source = cc$dynamics, 
      mask_smallcarp_feed_medcod, 
      fun_smallcarp_feed_medcod
    )
    
    #   - Moderate negative effect of large MC on rainbowfish
    mask_largecod_eat_rb <- transition(
      rb$dynamics$matrix, dim = rb_stages
    )
    fun_largecod_eat_rb <- interaction_fn(
      k = kmc, theta = 0.22, stage = large_mc, negative = TRUE
    ) 
    largecod_eat_rb <- pairwise_interaction(
      target = rb$dynamics,
      source = mc$dynamics, 
      mask_largecod_eat_rb, 
      fun_largecod_eat_rb
    )
    
    #   - Weak negative effect of medium MC on rainbowfish
    mask_medcod_eat_rb <- transition(
      rb$dynamics$matrix, dim = rb_stages
    )
    fun_medcod_eat_rb <- interaction_fn(
      k = kmc, theta = 0.05, stage = medium_fish, negative = TRUE
    ) 
    medcod_eat_rb <- pairwise_interaction(
      target = rb$dynamics,
      source = mc$dynamics, 
      mask_medcod_eat_rb, 
      fun_medcod_eat_rb
    )
    
    #   - Weak negative effect of large carp on rainbowfish
    mask_largecarp_eat_rb <- transition(
      rb$dynamics$matrix, dim = rb_stages
    )
    fun_largecarp_eat_rb <- interaction_fn(
      k = kcarp, theta = 0.05, stage = large_carp, negative = TRUE
    ) 
    largecarp_eat_rb <- pairwise_interaction(
      target = rb$dynamics,
      source = cc$dynamics, 
      mask_largecarp_eat_rb, 
      fun_largecarp_eat_rb
    )
    
    #   - Weak negative effect of medium carp on rainbowfish
    mask_medcarp_eat_rb <- transition(
      rb$dynamics$matrix, dim = rb_stages
    )
    fun_medcarp_eat_rb <- interaction_fn(
      k = kcarp, theta = 0.05, stage = medium_fish, negative = TRUE
    ) 
    medcarp_eat_rb <- pairwise_interaction(
      target = rb$dynamics,
      source = cc$dynamics, 
      mask_medcarp_eat_rb, 
      fun_medcarp_eat_rb
    )
    
    #   - Weak positive effect of rainbowfish on large MC
    mask_rb_feed_largecod <- transition(
      mc$dynamics$matrix, dim = large_mc
    )
    fun_rb_feed_largecod <- interaction_fn(
      k = krb, theta = 0.2, stage = rb_stages, negative = FALSE
    ) 
    rb_feed_largecod <- pairwise_interaction(
      target = mc$dynamics,
      source = rb$dynamics, 
      mask_rb_feed_largecod, 
      fun_rb_feed_largecod
    )
    
    #   - Weak positive effect of rainbowfish on medium MC
    mask_rb_feed_medcod <- transition(
      mc$dynamics$matrix, dim = medium_fish
    )
    fun_rb_feed_medcod <- interaction_fn(
      k = krb, theta = 0.2, stage = rb_stages, negative = FALSE
    ) 
    rb_feed_medcod <- pairwise_interaction(
      target = mc$dynamics,
      source = rb$dynamics, 
      mask_rb_feed_medcod, 
      fun_rb_feed_medcod
    )
    
    #   - Weak positive effect of rainbowfish on large carp
    mask_rb_feed_largecarp <- transition(
      cc$dynamics$matrix, dim = large_carp
    )
    fun_rb_feed_largecarp <- interaction_fn(
      k = krb, theta = 0.2, stage = rb_stages, negative = FALSE
    ) 
    rb_feed_largecarp <- pairwise_interaction(
      target = cc$dynamics,
      source = rb$dynamics, 
      mask_rb_feed_largecarp, 
      fun_rb_feed_largecarp
    )
    
    #   - Weak positive effect of rainbowfish on medium carp
    mask_rb_feed_medcarp <- transition(
      cc$dynamics$matrix, dim = medium_fish
    )
    fun_rb_feed_medcarp <- interaction_fn(
      k = krb, theta = 0.2, stage = rb_stages, negative = FALSE
    ) 
    rb_feed_medcarp <- pairwise_interaction(
      target = cc$dynamics,
      source = rb$dynamics, 
      mask_rb_feed_medcarp, 
      fun_rb_feed_medcarp
    )
    
    # collate interactions (16 in total)
    interactions <- list(
      largecarp_eat_smallcod,
      smallcod_feed_largecarp,
      largecod_eat_smallcarp,
      smallcarp_feed_largecod,
      medcarp_eat_smallcod,
      smallcod_feed_medcarp,
      medcod_eat_smallcarp,
      smallcarp_feed_medcod,
      largecod_eat_rb,
      rb_feed_largecod,
      largecarp_eat_rb,
      rb_feed_largecarp,
      medcod_eat_rb,
      rb_feed_medcod,
      medcarp_eat_rb,
      rb_feed_medcarp
    )
    
  } else {
    
    # pull out target dynamics objects
    cc <- pops[[which(sp_list == "Common carp")]]
    bf <- pops[[which(sp_list == "River blackfish")]]
    
    # specify age classes for each size class
    small_fish <- 1
    medium_fish <- 2
    large_bf <- 3:11 
    large_carp <- 3:28
    kbf <- scale * 10000
    kcarp <- scale * 100000
      
    # specify interactions one-by-one (reciprocating all unless explicitly not
    #    required)    
    
    #   - Strong negative effect of large carp on small blackfish
    mask_largecarp_eat_smallbf <- reproduction(bf$dynamics$matrix)
    fun_largecarp_eat_smallbf <- interaction_fn(
      k = kcarp, theta = 0.7, stage = large_carp, negative = TRUE, cut = FALSE
    ) 
    largecarp_eat_smallbf <- pairwise_interaction(
      target = bf$dynamics,
      source = cc$dynamics, 
      mask_largecarp_eat_smallbf, 
      fun_largecarp_eat_smallbf
    )
    
    #   - Strong positive effect of small blackfish on large carp
    mask_smallbf_feed_largecarp <- transition(
      cc$dynamics$matrix, dim = large_carp
    )
    fun_smallbf_feed_largecarp <- interaction_fn(
      k = kbf, theta = 5, stage = small_fish, negative = FALSE
    ) 
    smallbf_feed_largecarp <- pairwise_interaction(
      target = cc$dynamics,
      source = bf$dynamics, 
      mask_smallbf_feed_largecarp, 
      fun_smallbf_feed_largecarp
    )
    
    #   - Strong negative effect of large blackfish on small carp
    mask_largebf_eat_smallcarp <- reproduction(cc$dynamics$matrix)
    fun_largebf_eat_smallcarp <- interaction_fn(
      k = kbf, theta = 0.7, stage = large_bf, negative = TRUE, cut = FALSE
    ) 
    largebf_eat_smallcarp <- pairwise_interaction(
      target = cc$dynamics,
      source = bf$dynamics, 
      mask_largebf_eat_smallcarp, 
      fun_largebf_eat_smallcarp
    )
    
    #   - Strong positive effect of small carp on large blackfish
    mask_smallcarp_feed_largebf <- transition(
      bf$dynamics$matrix, dim = large_bf
    )
    fun_smallcarp_feed_largebf <- interaction_fn(
      k = kcarp, theta = 5, stage = small_fish, negative = FALSE
    ) 
    smallcarp_feed_largebf <- pairwise_interaction(
      target = bf$dynamics,
      source = cc$dynamics, 
      mask_smallcarp_feed_largebf, 
      fun_smallcarp_feed_largebf
    )
    
    #   - Moderate negative effect of medium carp on small blackfish
    mask_medcarp_eat_smallbf <- reproduction(bf$dynamics$matrix)
    fun_medcarp_eat_smallbf <- interaction_fn(
      k = kcarp, theta = 0.22, stage = medium_fish, negative = TRUE, cut = FALSE
    ) 
    medcarp_eat_smallbf <- pairwise_interaction(
      target = bf$dynamics,
      source = cc$dynamics, 
      mask_medcarp_eat_smallbf, 
      fun_medcarp_eat_smallbf
    )
    
    #   - Moderate positive effect of small blackfish on medium carp
    mask_smallbf_feed_medcarp <- transition(
      cc$dynamics$matrix, dim = medium_fish
    )
    fun_smallbf_feed_medcarp <- interaction_fn(
      k = kbf, theta = 0.85, stage = small_fish, negative = FALSE
    ) 
    smallbf_feed_medcarp <- pairwise_interaction(
      target = cc$dynamics,
      source = bf$dynamics, 
      mask_smallbf_feed_medcarp, 
      fun_smallbf_feed_medcarp
    )
    
    #   - Weak negative effect of medium blackfish on small carp
    mask_medbf_eat_smallcarp <- reproduction(cc$dynamics$matrix)
    fun_medbf_eat_smallcarp <- interaction_fn(
      k = kbf, theta = 0.05, stage = medium_fish, negative = TRUE, cut = FALSE
    ) 
    medbf_eat_smallcarp <- pairwise_interaction(
      target = cc$dynamics,
      source = bf$dynamics, 
      mask_medbf_eat_smallcarp, 
      fun_medbf_eat_smallcarp
    )
    
    #   - Weak positive effect of small carp on medium blackfish
    mask_smallcarp_feed_medbf <- transition(
      bf$dynamics$matrix, dim = medium_fish
    )
    fun_smallcarp_feed_medbf <- interaction_fn(
      k = kcarp, theta = 0.2, stage = small_fish, negative = FALSE
    ) 
    smallcarp_feed_medbf <- pairwise_interaction(
      target = bf$dynamics,
      source = cc$dynamics, 
      mask_smallcarp_feed_medbf, 
      fun_smallcarp_feed_medbf
    )
    
    # collate interactions (8 in total)
    interactions <- list(
      largecarp_eat_smallbf,
      smallbf_feed_largecarp,
      largebf_eat_smallcarp,
      smallcarp_feed_largebf,
      medcarp_eat_smallbf,
      smallbf_feed_medcarp,
      medbf_eat_smallcarp,
      smallcarp_feed_medbf
    )
    
  }
  
  # return
  interactions
  
}
