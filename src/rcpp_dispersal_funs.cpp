#include "RcppArmadillo.h"
using namespace Rcpp;

//' dispersal function for dynamic metapopulation models
//' @param current_distribution raster of current population distribution.
//' @param habitat_suitability raster of habitat suitability that has been converted to carrying capacity
//' @param barrier_map raster of barriers to the dispersal, 1=barrier, 0=no barrier.
//' @export
// [[Rcpp::export]]

NumericMatrix a_dispersal_function(NumericMatrix current_population_state, //raster
						NumericMatrix potiential_carrying_capacity,
						NumericMatrix habitat_suitability_map, //raster
						NumericMatrix barriers_map, //raster
						int dispersal_steps,
						int dispersal_distance,
						int barrier_type,
						bool use_barriers = true,
						double dispersal_kernal,
						double dispersal_probability){

	arma::mat cps = as<arma::mat>(current_population_state);
	arma::mat pcc = as<arma::mat>(potiential_carrying_capacity);	
	arma::mat hsm = as<arma::mat>(habitat_suitability_map);
	arma::mat barriers = as<arma::mat>(barriers_map);
	int ncol = cdr.n_cols;
    int nrow = cdr.n_rows;
    arma::mat cca(nrow,ncol); // carrying capacity avaliable
    arma::mat fps(nrow,ncol);
    cca.fill(NA_REAL);
    int nr_step_colonized, nr_step_decolonized, loopID;

	// check how much carrying capacity is free per-cell - this will enable dispersal to these cells if needed.
     for(int i = 0; i < nrows; i++){
         for(int j = 0; j < ncols; j++){
			 if(arma::is_finite(cdr(i,j))){ 
	            cca(i,j) = pcc(i,j) - cps(i,j); // carrying capacity avaliable = potiential carrying capacity - current population state.
				//fps(i,j) = cps(i,j);
	     }
	   }
     }
     
      /* Filter the habitat suitability matrix in three ways:
      **  -> replace any value < 0 by 0 (this removes NoData).
      **  -> set habitat suitability to 0 where barrier = 1.
      **  -> set habitat suitability values to NoData where barrier = NoData. 
      ** and also  */
      cca_cleaned = clean_matrix(cca, barriers, true, true, true);
      
      
      /* "Unlimited" and "no dispersal" scenario pixel count. Here we compute the number of pixels
      ** that would be colonized if dispersal was unlimited or null. This is simply the sum of all
      ** potentially suitable habitats */
      n_dispersal_cells = total_dispersal_cells(cca_cleaned);
      update_no_dispersal_matrix(habitat_suitability, no_dispersal, &n_dispersal_cells);

      /* Reset number of decolonized cells within current dispersal step pixel counter */
	  nr_step_decolonized = 0;
	    
      /* Update for temporarily resilient pixels. */
      for(i = 0; i < nrows; i++){
	    for(j = 0; j < ncols; j++){
	      
	      /* Update non-suitable pixels. If a pixel turned unsuitable, we update its status to "Temporarily Resilient". */
	      if(hsr(i,j) == 0) && (cca_cleaned(i,j) > 0)){
	        
	        /* If the user selected TemporaryResilience==T, then the pixel is set to "Temporary Resilient" status. */
	        if(tempResilience == true){
	          current_state[i][j] = 29900;
	        }
	        else{
	          /* If not temporary resilience was specified, then the pixel is set to "decolonized" status. */
	          current_state[i][j] = -1 - loopID;
	        }
	        /* The number of decolonized cells within current step is increased by one */
	        nr_step_decolonized++;
	      }
	    }
      }

      
      /* *************************************** */
      /* Dispersal event step loop starts here.  */
      /* *************************************** */
      loopID = 0;
      for(dispersal_step = 1; dispersal_step <= dispersal_steps; dispersal_step++){
	    
	    /* Set the value of "loopID" for the current iteration of the dispersal loop. */
	    loopID++;
          
	    /* Reset pixel counters that count pixels within the current loop. */
	    nr_step_colonized = 0;
	    if(dispersal_step > 1) nr_step_decolonized = 0;

	    /* Source cell search: Can the sink pixel be colonized? There are four
	    ** conditions to be met for a sink pixel to become colonized:
	    **   1. Sink pixel is currently suitable. 
	    **	 2. Sink pixel has avaliable carrying capacity.
	    **   3. Sink pixel is within dispersal distance of an already colonised.
	    **   4. There is no obstacle (barrier) between the pixel to be colonised
	    **      (sink pixel) and the pixel that is already colonised (source
	    **      pixel).
	    **
	    ** Loop through the cellular automaton. */
	    for(i = 0; i < nrows; i++){
	      for(j = 0; j < ncols; j++){
	        
		    //The are boolian calls which indicate if habitat is suitable and if cell can disperse
	        habitat_is_suitable = false;
	        cell_in_dispersal_distance = false;
	        
	        /* 1. Test whether the pixel is a suitable sink (i.e., its habitat
	        **    is suitable, it's unoccupied and is not on a barrier or filter
	        **    pixel). */
	        if((hsm(i,j) > 0) && (cca_cleaned(i,j) > 0)) habitat_is_suitable = true;

	        /* 2. Test whether there is a source cell within the dispersal
	        **    distance. To be more time efficient, this code runs only if
	        **    the answer to the first question is positive. Additionally,
	        **    if there is a source cell within dispersion distance and if
	        **    a barrier was asked for, then we also check that there is no
	        **    barrier between this source cell and the sink cell (this
	        **    verification is carried out in the "SearchSourceCell"
	        **    function). */
	        if(habitat_is_suitable){
		      /* Now we search if there is a suitable source cell to colonize the sink cell. */
	          cell_in_dispersal_distance = can_source_cell_disperse(i, j, current_state, hsm, barriers));
	        }
	        
	        /* Update sink cell status. */
	        if(habitat_is_suitable &&  Rcpp::all(!Rcpp::is_na(cell_in_dispersal_distance))){
				
			  //sink cell = existing population + rbinom(1, population @ source cell, p=prob_of_dispersal)
			  //source cell = source cell population - rbinom realisation. 
	          
		      /* Only if the 2 conditions are fullfilled the cell's status is set to colonised. */
		      int source_x = cell_is_dispersal_distance[0];
		      int source_y = cell_is_dispersal_distance[1];
		      source_pop = cps(source_x,source_y);
		      source_pop_dispersed = rbinom(1,source_pop,dispersal_proportion);
	          fps(i,j) = cps(i,j) + source_pop_dispersed;
	          fps(sink_x,sink_y) = cps(source_x,source_y) - source_pop_dispersed;
	          nr_step_colonized++;
	        }
	      }
	    }
        

	    /* Update pixel counters. */
	    nrColonized = nrColonized + nrStepColonized - nrStepDecolonized;
	    nrAbsent = nrAbsent - nrStepColonized + nrStepDecolonized;
	    nrTotColonized += nrStepColonized;
	    nrTotDecolonized += nrStepDecolonized;
	    nrTotLDDSuccess += nrStepLDDSuccess;
	    	 
	    /* If the user has requested full output, also write the current state matrix to file. */
	    }
   return(wrap(fds));    /* end of dispersal */
}

/*
** can_source_cell_disperse: This function will search, for a given input location, for a
**            suitable "source" pixel that can lead to the colonization of
**            the input location (sink pixel).
**
** Parameters:
**   - i:        The row number of the sink pixel.
**   - j:        The column number of the sink pixel.
**   - current_population_state: The matrix that contains the current state of the cellular
**               automaton.
**   - initial_population_state:   A matrix giving the "age" of each colonized pixel.
** 	 - habitat_suitability 	
**   - loopID:   Indicates in which "loop" the simulation currently is. The
**               'loopID' enables to retrieve the 'environmental change' loop
**               and the 'dispersal' loop that the simulation is currently in.
**   - habSuit:  The habitat suitability of the current pixel.
**   - barriers: A pointer to the barriers matrix.
**
** Returns:
**   If a suitable source cell was found: true.
**   Otherwise:                           false.
*/

IntegerVector can_source_cell_disperse(int i, 
							  int j,
							  NumericMatrix current_population_state,
							  NumericMatrix initial_population,
							  NumericMatrix habitat_suitability_map,
							  NumericMatrix barriers_map
							  bool use_barrier=true,
							  int dispersal_loop_ID //which dispersal loop are we in.
							  int dispersal_distance,
							  NumericVector dispersal_kernal,
							  double dispersal_proportion){
								  
  	arma::mat cca = as<arma::mat>(carrying_capacity_avaliable); //A matrix of avaliable carrying capacity.
  	arma::mat cps = as<arma::mat>(current_population_state); // the current population (during population step)
  	arma::mat hsm = as<arma::mat>(habitat_suitability_map); // the habitat suitability
	arma::mat barriers = as<arma::mat>(barrier_map);
	int ncol = cps.n_cols;
    int nrow = cps.n_rows;								  
	int    k, l, real_distance;
	double prob_colonisation, rnd;
	IntegerVector source_found(2);
	source_found.fill(NA_REAL);

  /*
  ** For now let's set these paramters to fixed values. Later we can implement
  ** them as variables.
  */
  //source_found = false;
        
  /*
  ** Search for a potential source cell. i and j are the coordinates of the
  ** sink cell. k and l are the coordinates of the potential source cell.
  */
  for (k = i - dispersal_distance; k <= i + dispersal_distance; k++){
    for (l = j - dispersal_distance; l <= j + dispersal_distance; l++){
      /*
      ** 1. Test of basic conditions to see if a pixel could be a potential
      **    source cell:
      **    - The pixel must be within the limits of the matrix's extent.
      **    - The pixel must be colonized, but not during the current loop.
      ** 	- The pixel must have avaliable carrying capacity to allow recruitment.
      */
      if ((k >= 0) && (k < nrRows) && (l >= 0) && (l < nrCols)){
		  if ((cps(k,l) > 0) && (cps(k,l) != loopID)){
			  if (cca(k,l) > 0){
	    /*
	    ** 2. Compute the distance between sink and (potential) source pixel
	    **    and check if it is <= maximum dispersal distance. The distance
	    **    is computed in pixel units.
	    **    real_distance = round(Sqr((K - I) ^ 2 + (L - J) ^ 2) + 0.5)
	    */
	    real_distance = round(sqrt((k-i)*(k-i) + (l-j)*(l-j)));
	    if ((real_distance > 0) && (real_distance <= dispersal_distance)){
	      /*
	      ** 3. Compute the probability of colonization of the sink pixel.
	      **    This probability depends on several factors:
	      **    - Disance between source and sink cells.
	      */
		  prob_colonisation = disp_kernel[real_distance-1] * (hsm(k,l)/1000.0);
	      rnd = runif(1);
	      if (rnd < prob_colonisation || prob_colonisation == 1.0){
				/*
				** When we reach this stage, the last thing we need to check for
				** is whether there is a "barrier" obstacle between the source
				** and sink pixel. We check this last as it requires significant
				** computing time.
				*/
				if (use_barrier){
				  if (!barrier_to_dispersal(i, j, k, l, barriers)){
					source_found[0] = k;
					source_found[1] = l;
					return(source_found);
				  }
				} else {
					source_found[0] = k;
					source_found[1] = l;
				  return(source_found);
				}
			  }
			}
		  }
		}
      }
    }
  }
}


// This function cleans up the habitat martix before dispersal.

NumericMatrix clean_matrix(NumericMatrix in_matrix,
						   NumericMatrix barrier_matrix,
						   bool filter_no_data = true,
						   bool filter_barriers = true,
						   bool insert_no_data = true){
	int i, j;
	arma::mat inmat = as<arma::mat>(in_matrix);
	arma::mat barriers = as<arma::mat>(barrier_matrix);
	int ncols = inmat.n_cols;
	int nrows = inmat.n_rows;
	  // set any value < 0 to 0, removes nan data and data where carrying  capacity is */
	  if(filter_no_data){
		for (i = 0; i < nrows; i++){
			for (j = 0; j < ncols; j++){
				if (inmat(i,j) < 0) inmat(i,j) = 0;
			}
		 }
	  }

	  /* Turn barrier cells into zero. */
	  if(filter_barriers){
		for (i = 0; i < nrows; i++){
			for (j = 0; j < ncols; j++){
				if (barriers(i,j) == 1) inmat(i,j) = 0;
			}
		 }
	  }
	  
	  /* turn NA in barrier_matrix into NA . */
	  if(insert_no_data){
		for (i = 0; i < nrRows; i++){
			for (j = 0; j < nrCols; j++){
				if (barriers(i,j) == -9999) inmat(i,j) = -9999;
			}
		 }
	  }
  
}


//' is there a barrier to dispersal?

bool barrier_to_dispersal(int snkX, int snkY, int srcX, int srcY, int & barrier_map){
  int  dstX, dstY, i, pxlX, pxlY, distMax, barCounter;
  bool barrier_found;

  barrier_found = false;
  
  /*
  ** Calculate the distance in both dimensions between the source and sink
  ** pixels and take the largest of the two.
  */
  dstX = srcX - snkX;
  dstY = srcY - snkY;
  if (abs (dstX) >= abs (dstY))
  {
    distMax = abs (dstX);
  }
  else
  {
    distMax = abs (dstY);
  }

  /*
  ** Check the possible paths from source to sink and see if there is a path
  ** without barriers.
  */
  if (barrier_type == WEAK_BARRIER)
  {
    /*
    ** Weak barrier: If there is at least one free path we're good.
    **
    ** BARRIER MIDDLE
    */
    barrier_found = false;
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX + (1.0 * i / distMax * dstX));
      pxlY = (int)round (snkY + (1.0 * i / distMax * dstY));
      if (barriers[pxlX][pxlY] == 1)
      {
	barrier_found = true;
	break;
      }
    }
    if (!barrier_found)
    {
      goto End_of_Routine;
    }
    /*
    ** BARRIER TOP_LEFT
    */
    barrier_found = false;
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX - 0.49 + (1.0 * i / distMax * dstX));
      pxlY = (int)round (snkY - 0.49 + (1.0 * i / distMax * dstY));
      if (barriers[pxlX][pxlY] == 1)
      {
	barrier_found = true;
	break;
      }
    }
    if (!barrier_found)
    {
      goto End_of_Routine;
    }
    /*
    ** BARRIER TOP_RIGHT
    */
    barrier_found = false;
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX + 0.49 + (1.0 * i / distMax * dstX));
      pxlY = (int)round (snkY - 0.49 + (1.0 * i / distMax * dstY));
      if (barriers[pxlX][pxlY] == 1)
      {
	barrier_found = true;
	break;
      }
    }
    if (!barrier_found)
    {
      goto End_of_Routine;
    }
    /*
    ** Barrier DOWN_LEFT
    */
    barrier_found = false;
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX - 0.49 + (1.0 * i / distMax * dstX));
      pxlY = (int)round (snkY + 0.49 + (1.0 * i / distMax * dstY));
      if (barriers[pxlX][pxlY] == 1)
      {
	barrier_found = true;
	break;
      }
    }
    if (!barrier_found)
    {
      goto End_of_Routine;
    }
    /*
    ** Barrier DOWN_RIGHT
    */
    barrier_found = false;
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX + 0.49 + (1.0 * i / distMax * dstX));
      pxlY = (int)round (snkY + 0.49 + (1.0 * i / distMax * dstY));
      if (barriers[pxlX][pxlY] == 1)
      {
	barrier_found = true;
	break;
      }
    }
    if (!barrier_found)
    {
      goto End_of_Routine;
    }
  }
  else if (barrierType == STRONG_BARRIER)
  {
    /*
    ** Strong barrier: If more than one way is blocked by a barrier then
    **                 colonization fails.
    */
    barCounter = 0;
    barrier_found = false;
    /*
    ** BARRIER MIDDLE
    */
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX + (1.0 * i / distMax * dstX));
      pxlY = (int)round (snkY + (1.0 * i / distMax * dstY));
      if (barriers[pxlX][pxlY] == 1)
      {
	barCounter++;
	break;
      }
    }
    /*
    ** BARRIER TOP_LEFT
    */
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX - 0.49 + (((i-1.0) / distMax * dstX) +
					((1.0 / distMax * dstX) / 2.0)));
      pxlY = (int)round (snkY - 0.49 + (((i-1.0) / distMax * dstY) +
					((1.0 / distMax * dstY) / 2.0)));
      if (barriers[pxlX][pxlY] == 1)
      {
	barCounter++;
	break;
      }
    }
    if (barCounter > 1)
    {
      barrier_found = true;
      goto End_of_Routine;
    }
    /*    
    ** BARRIER TOP_RIGHT
    */
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX + 0.49 + (((i-1.0) / distMax * dstX) +
					((1.0 / distMax * dstX) / 2.0)));
      pxlY = (int)round (snkY - 0.49 + (((i-1.0) / distMax * dstY) +
					((1.0 / distMax * dstY) / 2.0)));
      if (barriers[pxlX][pxlY] == 1)
      {
	barCounter++;
	break;
      }
    }
    if (barCounter > 1)
    {
      barrier_found = true;
      goto End_of_Routine;
    }
    /*
    ** BARRIER DOWN_LEFT
    */
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX - 0.49 + (((i-1.0) / distMax * dstX) +
					((1.0 / distMax * dstX) / 2.0)));
      pxlY = (int)round (snkY + 0.49 + (((i-1.0) / distMax * dstY) +
					((1.0 / distMax * dstY) / 2.0)));
      if (barriers[pxlX][pxlY] == 1)
      {
	barCounter++;
	break;
      }
    }
    if (barCounter > 1)
    {
      barrier_found = true;
      goto End_of_Routine;
    }
    /*
    ** BARRIER DOWN_RIGHT
    */
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX + 0.49 + (((i-1.0) / distMax * dstX) +
					((1.0 / distMax * dstX) / 2.0)));
      pxlY = (int)round (snkY + 0.49 + (((i-1.0) / distMax * dstY) +
					((1.0 / distMax * dstY) / 2.0)));
      if (barriers[pxlX][pxlY] == 1)
      {
	barCounter++;
	break;
      }
    }
    if (barCounter > 1)
    {
      barrier_found = true;
      goto End_of_Routine;
    }
  }        

 End_of_Routine:
  /*
  ** Return the result.
  */
  return (barrier_found);
}


int total_dispersal_cells(NumericMatrix habitat_suitability_map){
  
	int i, j, count;
	arma::mat hs = as<arma::mat>(habitat_suitability_map);
	int ncols = hs.n_cols;
	int nrows = hs.n_rows;

   // Count the number of dispersable pixels.
	  count = 0;
	  for (i = 0; i < nrows; i++){
		for (j = 0; j < ncols; j++){
		  if (hs(i,j) > 0) count++;
		}
	  }
  return (count);
}




/*
** updateNoDispMat: This function updates the "NoDispersal_Matrix" with the
**                  habitat suitability values that are contained in the
**                  current "HS_Matrix".
**
** Parameters:
**   - hsMat:       A pointer to the habitat suitability matrix.
**   - niDispMat:   A pointer to the no-dispersal matrix.
**   - noDispCount: A pointer to the no-dispersal count variable (its value
**                  will be updated!).
*/

void update_no_dispersal_matrix(int **hsMat, int **noDispMat, int *noDispCount)
{                   
  int i, j;

  if (*noDispCount > 0){
    for(i = 0; i < nrRows; i++){
      for(j = 0; j < nrCols; j++){
	    if((noDispMat[i][j] == 1) && (hsMat[i][j] == 0)){
	      noDispMat[i][j] = 0;
	      (*noDispCount)--;
	    }
      }
    }
  }
}

