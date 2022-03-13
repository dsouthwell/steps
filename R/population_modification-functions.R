#' How the population is modified in a landscape.
#'
#' Pre-defined functions to define population modification (e.g. translocation) during a simulation.
#'
#' @name population_modification_functions
#' 
#' @seealso
#' \itemize{
#'   \item{\link[steps]{translocation} for specifying explicit spatial and temporal movements
#'   of populations}
#'   \item{\link[steps]{mortality} for specifying explicit spatial and temporal changes to populations}
#'   }
NULL

#' Translocate populations
#'
#' This function is used to move or introduce populations throughout a simulation. A user can
#' specify which life-stages will be affected (default is all) and in which timesteps the
#' translocations will take place. A warning will be generated if populations are not available
#' where specified to translocate from.
#'
#' @param origins_layer the name of a spatial layer in the landscape object with the locations
#'   and number of individuals to translocate from. Note, this layer will have only zero
#'   values if individuals are being introduced from outside the study area
#' @param destinations_layer the name of a spatial layer in the landscape object with the locations
#'   and number of individuals to translocate to. Note, this layer will have only zero
#'   values if individuals are being controlled (e.g. culling)
#' @param stages which life-stages are modified - default is all
#' @param effect_timesteps which timesteps in a single simulation do the translocations
#'   take place

#'
#' @export
#' 
#' @examples
#' 
#' # Modify populations in all life-stages using explicit layers of origin and destination populations
#' # in timesteps 5, 10, and 15.
#' 
#' \dontrun{
#' trans_pop <- translocation(origins_layer = "origins",
#'                            destinations_layer = "destinations",
#'                            stages = NULL,
#'                            effect_timesteps = c(5, 10, 15))
#' 
#' ls <- landscape(population = egk_pop,
#'                 suitability = NULL,
#'                 carrying_capacity = NULL,
#'                 "origins" = egk_origins,
#'                 "destinations" = egk_destinations)
#' 
#' pd <- population_dynamics(change = growth(egk_mat), modification = trans_pop)
#' 
#' simulation(landscape = ls, population_dynamics = pd, habitat_dynamics = NULL, timesteps = 20)
#' }

translocation <- function (origins_layer1, origins_layer2, origins_layer3, destinations_layer, random_release_sites, stages = NULL, effect_timesteps = 1) {
  
  pop_dynamics <- function (landscape, timestep) {
    
    if (timestep %in% effect_timesteps) {
      
      population_raster <- landscape$population
      nstages <- raster::nlayers(population_raster)
      
      # get population as a matrix
      idx <- which(!is.na(raster::getValues(population_raster[[1]])))
      population_matrix <- raster::extract(population_raster, idx)
      
      #Extract all origins layers
      origins1 <- raster::extract(landscape[[origins_layer1]][[timestep]], idx) #Maria Island
      origins2 <- raster::extract(landscape[[origins_layer2]][[timestep]], idx) #Forestier
      origins3 <- raster::extract(landscape[[origins_layer3]][[timestep]], idx) #Road collisions
      
      #Do binomial trial to determine road mortality
      origins3 <- ifelse(runif(length(origins3)) < origins3, 1, 0)
      
      #Now sum the different origins together
      origins <- origins1 + origins2 + origins3
      
      destinations <- raster::extract(landscape[[destinations_layer]][[timestep]], idx)
      
        if (random_release_sites == TRUE) {
        
        #Identify which cells dont have DFT1 or DFT2
        DFTD1_raster <- landscape$DFTD1
        DFTD2_raster <- landscape$DFTD2
        DFTF_vac <- which(raster::getValues(DFTD1_raster)==0 & raster::getValues(DFTD2_raster)==0)
        #Count how many release sites there were
        new_destinations <- sample(DFTF_vac, length(which(destinations > 0)), replace=FALSE)
        #Randomly select the same number of new release sites that aren't diseased
        n_releases <- max(destinations)
        destinations[] <- 0
        destinations[new_destinations] <- n_releases
        #Set the number of releases 
        
      }
      
      if (is.null(stages)) stages <- seq_len(nstages)
      
      for (stage in stages) {
        
        warn_once(any(population_matrix[ , stage] - origins < 0),
                  paste("The proposed number of translocated individuals do not exist for\nlife-stage",
                        stage,
                        "- only the maximum number of available\nindividuals in each cell will be translocated. Please check the\nspecified origins and destination layers."),
                  warning_name = paste0("translocated_individuals_", stage))
        
        population_matrix[ , stage] <- population_matrix[ , stage] + destinations - pmin(origins, population_matrix[ , stage])
        
      }
      
      #########################################
      #Reset heterozugosity for translocated populations
      
      allele_raster <- landscape$allele_frequencies
      H_raster <- landscape$Heterozygosity
      allele_matrix <- raster::extract(allele_raster, idx)
      H_matrix <- raster::extract(H_raster, idx)
      
      #Calculate weighted heterozygosity of source population (origins 1 and origins 2)
      
      #allele_source1 <- allele_source2 <- allele_source <- matrix(0, ncol=11, nrow=1)
      
      #if (sum(origins1) > 0) {
        #allele_source1 <- allele_matrix[which(origins1 > 0),]
        #allele_mean1 <- colMeans(allele_source1)
      #} else {allele_source1 <-  matrix(0, ncol=11, nrow=1)}
      
      #if (sum(origins2) > 0) {
        #allele_source2 <- allele_matrix[which(origins2 > 0),]
        #allele_mean2 <- colMeans(allele_source2)
      #} else {allele_source2 <-  matrix(0, ncol=11, nrow=1)}
      
     # wgts_source <- c(sum(origins1)/(sum(origins1) + sum(origins2)), sum(origins2)/(sum(origins1) + sum(origins2)))
      #for (ss in 1:length(allele_source1)) {
        #allele_source[ss] <- weighted.mean(c(allele_source1[ss], allele_source2[ss]), wgts_source)
      #}
      #allele_source <- colMeans(allele_source)
      
      #if (sum(destinations) > 0) {
        
        #sink_locations <- which(destinations > 0)
        #old_alleles <- allele_matrix[sink_locations,]
        #old_popsize <- rowSums(population_matrix)[sink_locations]
      
        #new_popsize <- old_popsize + destinations[sink_locations]
        #wgts_old <- old_popsize/new_popsize
       # wgts_trans <- 1 - wgts_old
      
       # new_alleles <- trans_alleles <- old_alleles
       # trans_alleles <- matrix(allele_source, nrow=nrow(old_alleles), ncol=ncol(new_alleles), byrow=TRUE) #rnorm(length(trans_alleles), 0.5, 0.05)
       # for (kk in 1:ncol(new_alleles)) {
       #   for (jj in 1:nrow(new_alleles)) {
       #     new_alleles[jj,kk] <- weighted.mean(c(old_alleles[jj,kk], trans_alleles[jj,kk]), c(wgts_old[jj], wgts_trans[jj]))
       #   }
       #}
        
        #H_new <- calc_heterozygosity(frequencies = new_alleles)
        #allele_matrix[sink_locations, ] <- new_alleles
        
        
        #allele_raster[idx] <- allele_matrix
        
        #H_raster <- calc_heterozygosity(frequencies = allele_raster)
        #landscape$allele_frequencies <- allele_raster
        
        #landscape$Heterozygosity <- H_raster
        
      #}
      
      #########################################
      #Reset inbreeding coefficient for translocated populations
      
      inbreeding_raster <- landscape$inbreeding
      diversity_raster <- landscape$Heterozygosity_relative
      inbreeding_matrix <- raster::extract(inbreeding_raster, idx)
      diversity_matrix <- raster::extract(diversity_raster, idx)
      
      F_source1 <- F_source2 <- F_source <- 0
      
      #Maria Island
      if (sum(origins1) > 0) {
        F_source1 <- mean(inbreeding_matrix[which(origins1 > 0)])
        F_source[] <- 0
      } 
      #Forestier
      if (sum(origins2) > 0) { F_source2 <- mean(inbreeding_matrix[which(origins2 > 0)])} 
      
      if (sum(origins1) > 0 | sum(origins2) > 0){
        wgts_source <- c(sum(origins1)/(sum(origins1) + sum(origins2)), sum(origins2)/(sum(origins1) + sum(origins2)))
        F_source <- weighted.mean(c(F_source1, F_source2), wgts_source)
      }
      
      if (sum(destinations) > 0) {
      sink_locations <- which(destinations > 0)
      old_F <- inbreeding_matrix[sink_locations]
      old_popsize <- rowSums(population_matrix)[sink_locations]
      
      new_popsize <- old_popsize + destinations[sink_locations]
      wgts_old <- old_popsize/new_popsize
      wgts_trans <- 1 - wgts_old
      
      new_F <- old_F
      for (ww in 1:length(old_F)){
        new_F[ww] <- weighted.mean(c(old_F[ww], F_source), c(wgts_old[ww], wgts_trans[ww]))
      }
      
      inbreeding_matrix[sink_locations] <- new_F
      diversity_matrix[sink_locations] <- 1 - new_F
      
      inbreeding_raster[idx] <- inbreeding_matrix
      diversity_raster[idx] <- diversity_matrix
      
      landscape$inbreeding <- inbreeding_raster
      landscape$Heterozygosity_relative <- diversity_raster
      }
      
  
      
      # put back in the raster
      population_raster[idx] <- population_matrix
      
      landscape$population <- population_raster
      
      
      
    }
    
    landscape
    
  }
  
  as.population_modification(pop_dynamics)
  
}


#' Directly affect populations
#' 
#' This function modifies a population by a mortality spatial layer included in a
#' steps landscape object. The mortality layer consists of values from 0???1 and
#' modifies the population by multiplying the population of a cell by the value of
#' the corresponding cell in a mortality layer. For example, a cell with ten
#' individuals before the mortality function is applied, and corresponding mortality
#' layer cell with a value of 0.2, would have two individuals remaining after
#' modification. Note, rounding also occurs after modification using a ceiling method
#' (i.e the largest whole integer is retained).

#' @param mortality_layer the name of spatial layer(s) in the landscape object with
#'   mortality proportions used to alter the populations for each timestep. If
#'   a stack of rasters is used then the number of layers must match the intended
#'   number of timesteps in the simulation.
#' @param stages which life-stages are modified - default is all
#'
#' @export
#' 
#' @examples
#' # Modify populations in all life-stages with fire intensity.
#' 
#' \dontrun{
#' fire_mortal <- mortality(mortality_layer = "fire", stages = NULL)
#' 
#' ls <- landscape(population = egk_pop,
#'                 suitability = egk_hab,
#'                 carrying_capacity = egk_k,
#'                 "fire" = egk_fire)
#' 
#' pd <- population_dynamics(change = growth(egk_mat), modification = fire_mortal)
#' 
#' simulation(landscape = ls, population_dynamics = pd, habitat_dynamics = NULL, timesteps = 20)
#' }

mortality <- function (mortality_layer, stages = NULL) {
  
  pop_dynamics <- function (landscape, timestep) {
    
    population_raster <- landscape$population
    nstages <- raster::nlayers(population_raster)
    
    # get population as a matrix
    idx <- which(!is.na(raster::getValues(population_raster[[1]])))
    population_matrix <- raster::extract(population_raster, idx)
    
    if (raster::nlayers(landscape[[mortality_layer]]) > 1) {
      mortality_prop <- raster::extract(landscape[[mortality_layer]][[timestep]], idx)
    } else {
      mortality_prop <- raster::extract(landscape[[mortality_layer]], idx)
    }
    
    if (is.null(stages)) stages <- seq_len(nstages)
    
    for (stage in stages) {
      
      population_matrix[ , stage] <- ceiling(population_matrix[ , stage] * mortality_prop)
      
    }
    
    # put back in the raster
    population_raster[idx] <- population_matrix
    
    landscape$population <- population_raster
    
    
    landscape
    
  }
  
  as.population_modification(pop_dynamics)
  
}

##########################
### internal functions ###
##########################

as.population_modification <- function (modification) {
  as_class(modification, "population_modification", "function")
}