#' How the population changes in a landscape.
#'
#' Pre-defined or custom functions to define population change during a simulation.
#' Please see the tutorial vignette titled "Creating custom *steps* functions"
#' for information on how to write custom functions for use in simulations.
#'
#' @name population_change_functions
#'
#' @seealso
#' \itemize{
#'   \item{\link[steps]{growth} is a default function for changing populations based on
#'   transition matrices and functions}
#'   }
NULL

#' Population growth
#' 
#' This function applies negative or positive growth to the population using matrix
#' multiplication. Stochasticity can be added to cell-based transition matrices or globally.
#' Users can also specify a built-in or custom function to modify the transition matrices
#' throughout a simulation. Please see the tutorial vignette titled "Creating custom
#' *steps* functions" for information on how to write custom functions for use in simulations.
#'
#' @param transition_matrix A symmetrical age-based (Leslie) or stage-based (Lefkovitch)
#'   population structure matrix.
#' @param global_stochasticity,local_stochasticity Either scalar values or
#' matrices (with the same dimension as \code{transition_matrix}) specifying
#' the variability in the transition matrix either for populations in all grid
#' cells (\code{global_stochasticity}) or for each grid cell population
#' separately (\code{local_stochasticity}). Values supplied here are the
#' standard deviation of a truncated normal distribution where the mean is the
#' value supplied in the transition matrix.
#' @param transition_function A function to specify or modify life-stage transitions
#'   at each timestep. See \link[steps]{transition_function}.
#' @param transition_order Order of transitions performed in growth function. This behaviour
#'   is only applied when demographic stochasticity is set to "full" (default) and transitions
#'   are applied sequentially. By default "fecundity" is performed first (calculating the
#'   number of new individuals to be added to the populations), then "survival" is applied.
#'   The final population is the sum of these. Users should be cautious of specifying
#'   "survival" to be performed first as typically survival of reproductive stages will already
#'   be accounted for in the fecundity values of the transition matrix.
#' @param two_sex Does the transition matrix include life stages for two sexes (i.e. male and
#'   female)? Default is FALSE which assumes a single sex matrix (e.g. females only).
#' 
#' @export
#' 
#' @examples
#' 
#' # Example of a growth function that changes the populations based on a transition matrix that
#' # is subject to global stochasticity. 
#' 
#' \dontrun{
#' stoch_growth <- growth(transition_matrix = egk_mat, global_stochasticity = egk_mat_stoch)
#' 
#' ls <- landscape(population = egk_pop, suitability = NULL, carrying_capacity = NULL)
#' 
#' pd <- population_dynamics(change = stoch_growth)
#' 
#' simulation(landscape = ls, population_dynamics = pd, habitat_dynamics = NULL, timesteps = 20)
#' }

growth <- function (transition_matrix,
                    global_stochasticity = 0,
                    local_stochasticity = 0,
                    transition_function = NULL,
                    transition_order = c("fecundity", "survival"),
                    two_sex = FALSE) { # updated 23.10.20
  
  transition_order <- match.arg(transition_order)
  is_function <- inherits(transition_function, "function")
  
  # if it's a list of functions, chain them together
  if (!is_function && inherits(transition_function, "list")) {
    
    # check the elements of this list are all functions, with the right arguments
    all_funs <- all(unlist(lapply(transition_function,
                                  function(x) inherits(x, "function"))))
    expected_params <- all(unlist(lapply(transition_function,
                                         function(x) names(formals(x)) %in% c("transition_array",
                                                                              "landscape",
                                                                              "timestep"))))
    
    if (!all_funs || !expected_params) {
      stop("A transition function list must contain only function objects that\n",
           "each include the arguments: transition_array, landscape, and timestep.\n",
           "Please check your inputs and re-run the simulation.")
    }
    
    transition_function_list <- transition_function
    
    
    transition_function <- function(transition_array, landscape, timestep) {
      
      for (fun in transition_function_list) {
        transition_array <- fun(transition_array, landscape, timestep)
      }
      
      transition_array
      
    }
    
    is_function <- TRUE
    
  }
  
  idx <- which(transition_matrix != 0)
  is_recruitment <- upper.tri(transition_matrix)[idx]
  if(two_sex == TRUE) { # added 23.10.20
    mat <- is.na(transition_matrix)
    mat[c(1, (nrow(transition_matrix) / 2) + 1), ] <- TRUE
    is_recruitment <- mat[idx]
  }
  upper <- ifelse(is_recruitment, Inf, 1)
  vals <- transition_matrix[idx]
  dim <- nrow(transition_matrix)
  
  if (is.matrix(global_stochasticity)) {
    stopifnot(identical(c(dim, dim), dim(global_stochasticity)))
    stopifnot(identical(which(global_stochasticity != 0), idx))
    global_stochasticity <- global_stochasticity[idx]
  }
  
  if (is.matrix(local_stochasticity)) {
    stopifnot(identical(c(dim, dim), dim(local_stochasticity)))
    stopifnot(identical(which(local_stochasticity != 0), idx))
    local_stochasticity <- local_stochasticity[idx]
  }
  
  pop_dynamics <- function (landscape, timestep) {
    
    population_raster <- landscape$population
    HO_raster <- landscape$Heterozygosity
    allele_stack <- landscape$allele_frequencies
    
    # get population as a matrix
    cell_idx <- which(!is.na(raster::values(population_raster[[1]])))
    population <- raster::extract(population_raster, cell_idx)
    Hcombine <- Hnew <- raster::extract(HO_raster, cell_idx)
    allele_matrix <- raster::extract(allele_stack, cell_idx)
    n_cells <- length(cell_idx)
    
    # calculate global and local noise
    global_noise <- stats::rnorm(length(idx), 0, global_stochasticity)
    local_noise <- stats::rnorm(length(idx) * n_cells, 0, local_stochasticity)
    total_noise <- global_noise + local_noise
    
    # pad the index to get corresponding elements in each slice  
    addition <- length(transition_matrix) * (seq_len(n_cells) - 1)
    idx_full <- as.numeric(outer(idx, addition, FUN = "+"))
    
    if (is_function) {
      
      # create transition array and fill with initial matrix values
      transition_array <- array(0, c(dim, dim, n_cells))
      transition_array[] <- transition_matrix[]
      
      # update the transition array
      transition_array <- transition_function(transition_array, landscape, timestep)
      
      values <- transition_array[idx_full] + total_noise
    } else {
      transition_array <- array(0, c(dim, dim, n_cells))
      values <- vals + total_noise
    }
    
    values <- pmax_zero(values)
    values <- pmin(values, rep(upper, n_cells))
    transition_array[idx_full] <- values
    
    if (steps_stash$demo_stochasticity == "full") {
      
      
      #H_new <- Q_new  <- allele_matrix
      
      #Loop through each loci
      #for (kk in 1:n_loci) {
        
        #outsurv <- genotypes <- split_genotype(Population=population, allele=allele_matrix[,kk])
        
        #Change the survival rates is modelling an adaptive loci
        #if (kk == n_loci & timestep > 150) {
          
          #rws <- dim(transition_array)[1]
          #resistance <- which(rowSums(genotypes[[3]])>0 & disease==1)
          #if (length(resistance) > 0){
            #transition_array[2:rws,,resistance] <- transition_array[2:rws,,resistance] * 1.00
          #}
        #}
        
        #Split by females and males and pair for breeding
        #females <- genotypes[,1:6,]
        #males <- genotypes[,7:12,]
        
    
        #for (v in 1:dim(genotypes)[1]) {
          #pp_females <- sum(females[v,,1]) 
          #pq_females <- sum(females[v,,2])
          #Ind_females <- data.frame("Females" = c(1:sum(females[v,,1:3])), "genotype" = NA, "Pair" = 0)
          #Ind_females[1:pp_females,2] <- "pp"
          #Ind_females[(pp_females+1):(pp_females + pq_females),2] <- "pq"
          #Ind_females[is.na(Ind_females)] <- "qq"
          
          #pp_males <- sum(males[v,,1]) 
          #pq_males <- sum(males[v,,2])
          #Ind_males <- data.frame("Males" = c(1:sum(males[v,,1:3])), "genotype" = NA)
          #Ind_males[1:pp_males,2] <- "pp"
          #Ind_males[(pp_males+1):(pp_males + pq_males),2] <- "pq"
          #Ind_males[is.na(Ind_males)] <- "qq"
          #Ind_females[,3] <- sample(Ind_males[,2],length(Ind_females[,2]),replace=TRUE)
        #}
        
        
        
        #loop through genotypes
        #for (ss in 1:dim(genotypes)[3]) { 
          n_loci <- nlayers(allele_stack)
          population2 <- population
          
          total_pop <- rowSums(population2)
          if (sum(total_pop) == 0) {next}
          has_pop <- total_pop > 0 
          
          # if two-sex model, assumes that matrix is setup in a way that the fecundity values are
          # in the first, and half the row count + 1, rows - added 23.10.20
          fec_rows <- 1
          if (two_sex) fec_rows <- c(1, (nrow(transition_matrix) / 2) + 1)
          
          if (transition_order == "fecundity" && sum(has_pop) >= 1) {
            # first step - perform fecundity to add individuals to the populations
            new_population <- add_offspring_edits(population2[has_pop, ], transition_array[ , , has_pop], fec_rows) # updated 23.10.20
            
            
            #For each loci, randomly sample allele frequency of offspring given #offspring, number parents and parent allele frequency
            offspring_both_sexes <- rowSums(new_population) #total male and female offspring
            allele_matrix_offspring <- allele_matrix
            for (kk in 1:n_loci) {
              offspring_p <- rbinom(offspring_both_sexes, total_pop[has_pop] * 2, allele_matrix[,kk][has_pop])
              offspring_p <- offspring_p/(total_pop[has_pop] * 2)
              allele_matrix_offspring[,kk][has_pop] <- offspring_p
            }
            
            
            # second step - perform survival on new population
            surv_population <- surviving_population(population2[has_pop, ], transition_array[ , , has_pop], fec_rows) # updated 23.10.20
            
            allele_matrix_new <- allele_matrix
            
            for (kk in 1:n_loci) {
              
              genotypes_offspring <- split_genotype(Population=new_population, allele=allele_matrix_offspring[,kk])
              genotypes_adults <- split_genotype(Population=surv_population, allele=allele_matrix[,kk])
              genotypes_adults[ , fec_rows,] <- genotypes_adults[ , fec_rows,] + genotypes_offspring
              allele_matrix_new[has_pop,kk] <- (2 * rowSums(genotypes_adults[,,1]) + rowSums(genotypes_adults[,,2]))/(2*rowSums(genotypes_adults[,,1:3])) 
              
            }
            
            
            # add new and surviving populations
            surv_population[ , fec_rows] <- surv_population[ , fec_rows] + new_population
            population2[has_pop, ] <- surv_population
            
          }
          
          if (transition_order == "survival" && sum(has_pop) >= 1) {
            # first step - perform survival on population
            surv_population <- surviving_population(population2[has_pop, ], transition_array[ , , has_pop], fec_rows) # updated 23.10.20
            # second step - perform fecundity to add individuals to the populations
            new_population <- add_offspring_edits(surv_population, transition_array[ , , has_pop], fec_rows) # updated 23.10.20
            # add new and surviving populations
            surv_population[ , fec_rows] <- surv_population[ , fec_rows] + new_population
            population2[has_pop, ] <- surv_population
          }
          
          #outsurv[,,ss] <- population2 
          
        #} #End of loop through different genotypes for allele i 
        
        #Calculate new p and q based on survival and fecundity
        #if (kk==1) {pop_new <- outsurv[,,1] + outsurv[,,2] + outsurv[,,3]}
        #finalpop <- outsurv[,,1] + outsurv[,,2] + outsurv[,,3]
        
        #Adaptive allele
        #pnew <- (2 * rowSums(outsurv[,,1]) + rowSums(outsurv[,,2]))/(2*rowSums(finalpop)) #Check equation
        #qnew <- 1 - pnew    
        #Hnew <- 1 - (pnew)^2
        
        #Q_new[,kk] <- qnew
        #H_new[,kk] <- Hnew
        
      #} #End loci loop
      
      #Hcombine <- 1/n_loci*(rowSums(H_new))
      
    } else {
      
      population <- sapply(seq_len(n_cells),
                           function(x) transition_array[ , , x] %*% matrix(population[x, ]))
      
      population <- t(population)
      
    }
    
    # put back in the raster
    population_raster[cell_idx] <- population2
    #HO_raster[cell_idx] <- Hcombine
    allele_stack[cell_idx] <- allele_matrix_new
    H0_raster <- calc_heterozygosity(frequencies = allele_stack)
    
    landscape$Heterozygosity <- HO_raster
    landscape$population <- population_raster
    landscape$allele_frequencies <- allele_stack
    
    landscape
    
    # import components from landscape object ----#NOTE CHANGE SO THAT RASTER STACK IS LOADED, NOT INDIVIDUAL RASTER LAYERS FOR EACH LOCI
    #population_raster <- landscape$population
    #disease_raster <- landscape$DFTD1 --- #NOT SURE WHY I NEED DISEASE RASTER LOADED IN
    #HO_raster <- landscape$Heterozygosity
    
    #alleleA_raster <- landscape$allele_a    
    #alleleN1_raster <- landscape$allele_n
    #alleleN2_raster <- landscape$allele_n2 
    #alleleN3_raster <- landscape$allele_n3 
    #alleleN4_raster <- landscape$allele_n4 
    #alleleN5_raster <- landscape$allele_n5
    #alleleN6_raster <- landscape$allele_n6 
    #alleleN7_raster <- landscape$allele_n7 
    #alleleN8_raster <- landscape$allele_n8 
    #alleleN9_raster <- landscape$allele_n9 
    #alleleN10_raster <- landscape$allele_n10 
    
    # get population as a matrix
    #cell_idx <- which(!is.na(raster::values(population_raster[[1]])))
    #population <- raster::extract(population_raster, cell_idx)
    #disease <- raster::extract(disease_raster, cell_idx)
    #Hcombine <- Hnew <- raster::extract(HO_raster, cell_idx)
    
    #alleleA <- raster::extract(alleleA_raster, cell_idx)    
    #alleleN1 <- raster::extract(alleleN1_raster, cell_idx)
    #alleleN2 <- raster::extract(alleleN2_raster, cell_idx)
    #alleleN3 <- raster::extract(alleleN3_raster, cell_idx)
    #alleleN4 <- raster::extract(alleleN4_raster, cell_idx)
    #alleleN5 <- raster::extract(alleleN5_raster, cell_idx)
    #alleleN6 <- raster::extract(alleleN6_raster, cell_idx)
    #alleleN7 <- raster::extract(alleleN7_raster, cell_idx)
    #alleleN8 <- raster::extract(alleleN8_raster, cell_idx)
    #alleleN9 <- raster::extract(alleleN9_raster, cell_idx)
    #alleleN10 <- raster::extract(alleleN10_raster, cell_idx)
    
    #n_cells <- length(cell_idx)
   
    #In genetic layers, turn all cells that are NA to 0
    #alleleA[is.na(alleleA)] <- 0  
    #alleleN1[is.na(alleleN1)] <- 0  
    #alleleN2[is.na(alleleN2)] <- 0  
    #alleleN3[is.na(alleleN3)] <- 0  
    #alleleN4[is.na(alleleN4)] <- 0  
    #alleleN5[is.na(alleleN5)] <- 0  
    #alleleN6[is.na(alleleN6)] <- 0  
    #alleleN7[is.na(alleleN7)] <- 0  
    #alleleN8[is.na(alleleN8)] <- 0  
    #alleleN9[is.na(alleleN9)] <- 0  
    #alleleN10[is.na(alleleN10)] <- 0  
    #Hcombine[is.na(Hcombine)] <- 0 
    
    # calculate global and local noise
    #global_noise <- stats::rnorm(length(idx), 0, global_stochasticity)
    #local_noise <- stats::rnorm(length(idx) * n_cells, 0, local_stochasticity)
    #total_noise <- global_noise + local_noise
    
    # pad the index to get corresponding elements in each slice  
    #addition <- length(transition_matrix) * (seq_len(n_cells) - 1)
    #idx_full <- as.numeric(outer(idx, addition, FUN = "+"))
    
   # if (is_function) {
      
   #   # create transition array and fill with initial matrix values
   #   transition_array <- array(0, c(dim, dim, n_cells))
   #   transition_array[] <- transition_matrix[]
      
      # update the transition array
   #   transition_array <- transition_function(transition_array, landscape, timestep)
      
   #   values <- transition_array[idx_full] + total_noise
   # } else {
   #   transition_array <- array(0, c(dim, dim, n_cells))
   #   values <- vals + total_noise
   # }
    
   # values <- pmax_zero(values)
   # values <- pmin(values, rep(upper, n_cells))
   # transition_array[idx_full] <- values
   
   # if (steps_stash$demo_stochasticity == "full") {
      
    ############################################################
    #Split populations into the three genotypes
    ############################################################
      
   # len <- 11 #ONCE I STACK THE LOCI, READ IN NUMBER OF LAYERS HERE SO ITS NOT FIXED AT 11
    
   # H_new <- Q_new  <- vector(mode = "list", length = len)
    
   # alleles <- list(alleleN1, alleleN2, alleleN3, alleleN4, alleleN5, alleleN6, alleleN7, alleleN8, alleleN9, alleleN10, alleleA)
      
   #   for (kk in 1:len) {
      
  #      outsurv <- genotypes <- split_genotype(Population=population, allele=alleles[[kk]])
        
    #    if (kk == len & timestep>15) {
          
   #       rws <- dim(transition_array)[1]
    #      resistance <- which(rowSums(genotypes[[3]])>0 & disease==1)
    #      if (length(resistance) > 0){
    #        transition_array[2:rws,,resistance] <- transition_array[2:rws,,resistance] * 1.00
     #     }
    #    }
        
    #  for (ss in 1:length(genotypes)) { #I've replaced every population with population2 everytime in this loop#
      
    #    population2 <- genotypes[[ss]]
      
    #    total_pop <- rowSums(population2)
    #    if (sum(total_pop) == 0) {next}
    #    has_pop <- total_pop > 0 
      
        # if two-sex model, assumes that matrix is setup in a way that the fecundity values are
        # in the first, and half the row count + 1, rows - added 23.10.20
    #    fec_rows <- 1
    #    if (two_sex) fec_rows <- c(1, (nrow(transition_matrix) / 2) + 1)
      
     #   if (transition_order == "fecundity" && sum(has_pop) >= 1) {
     #     # first step - perform fecundity to add individuals to the populations
     #     new_population <- add_offspring(population2[has_pop, ], transition_array[ , , has_pop], fec_rows) # updated 23.10.20
          # second step - perform survival on new population
     #     surv_population <- surviving_population(population2[has_pop, ], transition_array[ , , has_pop], fec_rows) # updated 23.10.20
          # add new and surviving populations
     #     surv_population[ , fec_rows] <- surv_population[ , fec_rows] + new_population
      #    population2[has_pop, ] <- surv_population
        
     #   }
      
     #   if (transition_order == "survival" && sum(has_pop) >= 1) {
     #     # first step - perform survival on population
     #     surv_population <- surviving_population(population2[has_pop, ], transition_array[ , , has_pop], fec_rows) # updated 23.10.20
          # second step - perform fecundity to add individuals to the populations
     #     new_population <- add_offspring(surv_population, transition_array[ , , has_pop], fec_rows) # updated 23.10.20
          # add new and surviving populations
      #    surv_population[ , fec_rows] <- surv_population[ , fec_rows] + new_population
      #    population2[has_pop, ] <- surv_population
      #  }
      
     #   outsurv[[ss]] <- population2 ####################################
      
    #  } #End of loop through different genotypes for allele i #############################
      
      #Calculate new p and q based on survival and fecundity
    #  finalpop <- outsurv[[1]] + outsurv[[2]] + outsurv[[3]]
      #Adaptive allele
    #  pnew <- (2 * rowSums(outsurv[[1]]) + rowSums(outsurv[[2]]))/(2*rowSums(finalpop))
    #  qnew <- 1 - pnew    # Frequency of C allele
    #  pnew[is.nan(pnew)] <- 0
    #  qnew[is.nan(qnew)] <- 0
   #   Hnew <- 1 - (pnew)^2
      
    #  Q_new[[kk]] <- qnew
    #  H_new[[kk]] <- Hnew
      
    #  } #End allele loop
      
    #  Hcombine <- 1/len*(H_new[[1]] + H_new[[2]] + H_new[[3]] + H_new[[4]] + H_new[[5]] + H_new[[6]] + H_new[[7]] + H_new[[8]] + H_new[[9]] + H_new[[10]])
      
      
   # } else {
      
    #  population <- sapply(seq_len(n_cells),
     #                      function(x) transition_array[ , , x] %*% matrix(population[x, ]))
      
      #population <- t(population)
      
    #}
    
    # put back in the raster
    #population <- sum(outsurv) ################################
  
   # population_raster[cell_idx] <- finalpop
    #alleleN1_raster[cell_idx] <- Q_new[[1]]
   # alleleN2_raster[cell_idx] <- Q_new[[2]]
   # alleleN3_raster[cell_idx] <- Q_new[[3]]
   # alleleN4_raster[cell_idx] <- Q_new[[4]]
    #alleleN5_raster[cell_idx] <- Q_new[[5]]
    #alleleN6_raster[cell_idx] <- Q_new[[6]]
    #alleleN7_raster[cell_idx] <- Q_new[[7]]
    #alleleN8_raster[cell_idx] <- Q_new[[8]]
    #alleleN9_raster[cell_idx] <- Q_new[[9]]
    #alleleN10_raster[cell_idx] <- Q_new[[10]]
    #alleleA_raster[cell_idx] <- Q_new[[len]]
    #HO_raster[cell_idx] <- Hcombine
    
    #In allele layers, turn all locally extinct cells from 0 to NA
    #extinct <- which(sum(population_raster)==0) #################################
    #alleleA_raster[sum(population_raster)==0] <- NA ########################################
    #alleleN1_raster[sum(population_raster)==0] <- NA ###################################
    #alleleN2_raster[sum(population_raster)==0] <- NA ######################################
    #alleleN3_raster[sum(population_raster)==0] <- NA ######################################
    #alleleN4_raster[sum(population_raster)==0] <- NA ######################################
    #alleleN5_raster[sum(population_raster)==0] <- NA ######################################
    #alleleN6_raster[sum(population_raster)==0] <- NA ######################################
    #alleleN7_raster[sum(population_raster)==0] <- NA ######################################
    #alleleN8_raster[sum(population_raster)==0] <- NA ######################################
    #alleleN9_raster[sum(population_raster)==0] <- NA ######################################
    #alleleN10_raster[sum(population_raster)==0] <- NA ######################################
    #HO_raster[sum(population_raster)==0] <- NA ######################################
    #HO_raster[sink.xy] <- 0.5
    
    #population_raster[cell_idx] <- population
    #print(class(population(raster)))
    #landscape$population <- population_raster
    #landscape$allele_a <- alleleA_raster
    #landscape$allele_n <- alleleN1_raster
    #landscape$allele_n2 <- alleleN2_raster
    #landscape$allele_n3 <- alleleN3_raster
    #landscape$allele_n4 <- alleleN4_raster
    #landscape$allele_n5 <- alleleN5_raster
    #landscape$allele_n6 <- alleleN6_raster
    #landscape$allele_n7 <- alleleN7_raster
    #landscape$allele_n8 <- alleleN8_raster
    #landscape$allele_n9 <- alleleN9_raster
    #landscape$allele_n10 <- alleleN10_raster
    #landscape$Heterozygosity <- HO_raster
    
    #landscape
    
    
  }
  
  result <- as.population_growth(pop_dynamics)
  #landscape
  result
}


##########################
### internal functions ###
##########################

as.population_growth <- function (simple_growth) {
  as_class(simple_growth, "population_growth", "function")
}

add_offspring <- function (population, transition_array, fec_rows) { # updated 23.10.20
  
  pops <- population
  
  # updated 23.10.20
  new_offspring <- sapply(fec_rows,
                          function(x) {
                            # get fecundities for all eligible stages
                            if (class(transition_array) == "matrix") {
                              fecundities <- t(transition_array[x, ])
                            } else {
                              fecundities <- t(transition_array[x, , ])
                            }
                            
                            # get expected number, then do a poisson draw about this
                            expected_offspring <- fecundities * pops
                            new_offspring_stochastic <- expected_offspring
                            new_offspring_stochastic[] <- stats::rpois(length(expected_offspring), expected_offspring[])
                            
                            # sum stage 1s created by all other stages
                            rowSums(new_offspring_stochastic)
                          })
  
  return(new_offspring) # updated 23.10.20
}

add_offspring_edits <- function (population, transition_array, fec_rows) { # updated 23.10.20
  
  pops <- population
  
  # updated 23.10.20
  new_offspring <- sapply(fec_rows,
                          function(x) {
                            # get fecundities for all eligible stages
                            if (class(transition_array) == "matrix") {
                              fecundities <- t(transition_array[x, ])
                            } else {
                              fecundities <- t(transition_array[x, , ])
                            }
                            
                            # get expected number, then do a poisson draw about this
                            expected_offspring <- fecundities * pops
                            new_offspring_stochastic <- expected_offspring
                            new_offspring_stochastic[] <- stats::rpois(length(expected_offspring), expected_offspring[])
                            
                            # sum stage 1s created by all other stages
                            rowSums(new_offspring_stochastic)
                          })
  
  return(new_offspring) # updated 23.10.20
}



surviving_population <- function (population, transition_array, fec_rows) { # updated 23.10.20
  survival_array <- transition_array
  
  if (class(transition_array) == "matrix") {
    survival_array[fec_rows, ] <- 0
  } else {
    survival_array[fec_rows, , ] <- 0
  }
  
  if(inherits(population, "numeric")) {
    n_stage <- length(population)
  } else {
    n_stage <- ncol(population)
  }
  
  # loop through stages, getting the stages to which they move (if they survive)
  if(inherits(population, "numeric")) {
    survival_stochastic <- matrix(0, 1, n_stage)
  } else {
    survival_stochastic <- matrix(0, nrow(population), n_stage)
  }
  
  for (stage in seq_len(n_stage)) {
    
    # get the populations that have any individuals of this stage
    if(inherits(population, "numeric")) {
      pops <- population[stage]
    } else {
      pops <- population[, stage]
    }
    
    # probability of transitioning to each other stage
    if (class(survival_array) == "matrix") {
      probs <- t(survival_array[ , stage])
    } else {
      probs <- t(survival_array[ , stage, ])
    }
    
    # add on probability of dying
    
    surv_prob <- rowSums(probs)
    
    # check for sensible values
    if(any(surv_prob > 1)) {
      stop("Survival values greater than one have been detected. Please check to ensure\n",
           "that your transition matrix values or stochasticity (standard deviation) values\n",
           "are sensible and re-run the simulation.")
    }
    
    probs <- cbind(probs, 1 - surv_prob)
    
    # loop through cells with population (rmultinom is not vectorised on probabilities)
    new_stages <- matrix(NA, length(pops), n_stage)
    idx <- seq_len(n_stage)
    for (i in seq_len(length(pops))) {
      new_stages[i, ] <- stats::rmultinom(1, pops[i], probs[i, ])[idx, ]
    }
    
    # update the population
    survival_stochastic <- survival_stochastic + new_stages
    
  }
  
  return(survival_stochastic)
}
