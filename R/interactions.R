# TODO: update to have multiple strengths of interactions to 
#   see where these start to have an impact
# TODO: add mean-field interactions into here somehow
#   Think through consequences of reciprocal interactions,
#    might appear as flat but actually be equal in both directions
#    (therefore don't matter?)

# Interactions to include:
#
# Definitions:
#   Strong = up to 50% change, Moderate = 20%, Weak = 5%
#   Test double these values too
#
# Northern (small = 0-135 mm, large = > 500 mm):
#  (small MC = YOY, medium MC = 1-4, large MC = 5+)
#  (small RB = all)
#  (small CC = YOY, medium CC = 1-4, large CC = 5+)
#   - Strong negative effect of large carp on small MC
#   - Moderate negative effect of large MC on small carp
#   - Weak negative effect of medium carp on small MC
#   - Weak negative effect of medium MC on small carp
#   - Moderate negative effect of large MC on rainbowfish
#   - Moderate negative effect of medium MC on rainbowfish
#   - Moderate negative effect of large carp on rainbowfish
#   - Moderate negative effect of medium carp on rainbowfish
#
# Coastal (small = 0-80 mm, large = > 160 mm):
#   (small BF = YOY, medium BF = 1, large BF = 2+)
#   (small CC = YOY, medium CC = 1, large CC = 2+)
#   - Strong negative effect of large carp on small blackfish
#   - Moderate negative effect of large blackfish on small carp
#   - Moderate negative effect of medium carp on small blackfish
#   - Weak negative effect of medium blackfish on small carp
#   - Neutral effect of large fish on each other

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
