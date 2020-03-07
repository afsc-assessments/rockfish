///////////////////////////////////////////////////////////////////
// BSAI POP Model  
//	Has asymptotic fishery selectivity curve
//   
//
// Naming Conventions:
//
//  GENERAL:
//    styr, endyr begining year and ending year of model (catch data available)
//    nages       number of age groups considered
//    nyrs_       number of observations available to specific data set
//
//  DATA SPECIFIC:
//    catch_bio   Observed fishery catch biomass
//    nyrs_srv2   Number of years of slope survey index value (annual)
//    obs_srv2    Observed trawl slope survey index value (annual)
//    nyrs_srv3   Number of years of AI trawl survey index value (annual)
//    obs_srv3   Observed trawl AI trawl survey index value (annual)
//
//    oac_fish    ObsFselerved age comp from the fishery
//    oac_srv2 	  Observed age comp from the trawl survey
//
//    pred_catch  Predicted catch from the fishery
//    pred_srv2   Predicted index of the trawl survey
//
//    eac_fish    Expected age comp from fishery data
//    eac_srv2	  Expected age comp from trawl survey.
//
//    sel_fish    selectivity for fishery                
//    sel_srv2    selectivity for the trawl survey
//
//////////////////////////////////////////////////////////////////////////////
DATA_SECTION //------------------------------------------------------------------------------------------------
  init_int styr				// start year for model
  init_int styr_fish			// start year for the fishery
  init_int endyr			// end year for model and fishery
  init_int yrs_r      // number of years to peel for retrospective run
  init_int nages      // number of ages modeled
  init_int nages_dat // number of ages in data, can be less than the modeled ages
  init_ivector ages(1,nages)    // vector of ages (for model)
  init_ivector ages_dat(1,nages_dat)  // vector of ages (for data) 
  init_int nselages			// number of ages with estimated selectivities
  init_int rec_age			// recruitment age
  init_int nlen				// number of length groups
  init_ivector lengths(1,nlen) 		// vector of lengths
  init_int nyrs_fish			// number of years in the fishery
  init_ivector yrs_fish(1,nyrs_fish)    // vector of the index years in the fishery
  init_vector catch_bio(styr_fish,endyr)     // observed catch biomass
  //init_int nyrs_cpue			// number of years in  the CPUE index
  //init_ivector yrs_cpue(1,nyrs_cpue)	// vector of the index years in  the CPUE index
  //init_vector obs_cpue(1,nyrs_cpue)     // vector of observed CPUE values
  //init_int nyrs_srv2			// number of years in survey 2 (slope trawl survey)
  //init_ivector yrs_srv2(1,nyrs_srv2)	// vector of the index years in survey 2 (slope trawl survey)
  //init_vector obs_srv2(1,nyrs_srv2)     // vector of observed survey 2 biomasses
  //init_vector obs_srv2_sd(1,nyrs_srv2)   // vector of standard deviations of the survey 2 mean
  init_int nyrs_srv3			// number of years in survey 3 (wedge trawl survey)
  init_ivector yrs_srv3(1,nyrs_srv3)	// vector of the index years in survey 3 (wedge trawl survey)
  init_vector obs_srv3(1,nyrs_srv3)     // vector of observed survey 3 biomasses
  init_vector obs_srv3_sd(1,nyrs_srv3)   // vector of standard deviations of the survey 3 mean
  init_vector obs_srv3_lower(1,nyrs_srv3) // mean minus 2SD for survey 3
  init_vector obs_srv3_upper(1,nyrs_srv3) // mean plus 2SD for survey 3
  init_int nyrs_srv_ebs      // number of years in survey 2 (EBS slope survey)
  init_ivector yrs_srv_ebs(1,nyrs_srv_ebs)  // vector of the index years in survey 2 (EBS slope survey)
  init_vector obs_srv_ebs(1,nyrs_srv_ebs)     // vector of observed survey 2 biomasses
  init_vector obs_srv_ebs_sd(1,nyrs_srv_ebs)   // vector of standard deviations of the survey 2 mean
  init_vector obs_srv_ebs_lower(1,nyrs_srv_ebs) // mean minus 2SD for survey 2
  init_vector obs_srv_ebs_upper(1,nyrs_srv_ebs) // mean plus 2SD for survey 3
  init_matrix unbiasedages(1,nages,1,nages) // transition age error matrix for unbiased ages
  //init_matrix biasedages(1,nages,1,nages)   // transition age error matrix for biased ages  ****
  init_matrix translen(1,nlen,1,nages)     // transition matrix from ages to lengths
  //init_int nyrs_fish_biased_ac			// number of years with biased fishery age comps  ****
  //init_ivector yrs_fish_biased_ac(1,nyrs_fish_biased_ac) 	// vector of index years with biased fishery age comps  ****
  //init_matrix oac_fish_biased(1,nyrs_fish_biased_ac,1,nages) 	// observed biased fishery age comps  ****
  init_int nyrs_fish_unbiased_ac			// number of years with  unbiased fishery age comps
  init_ivector yrs_fish_unbiased_ac(1,nyrs_fish_unbiased_ac) 	// vector of index years with unbiased fishery age comps
  init_matrix oac_fish_unbiased(1,nyrs_fish_unbiased_ac,1,nages_dat) 	// observed unbiased fishery age comps
  init_int nyrs_fish_lc					// number of years with fishery length comps
  init_ivector yrs_fish_lc(1,nyrs_fish_lc) 		// vector of index years with fishery length comps
  init_matrix olc_fish(1,nyrs_fish_lc,1,nlen) 		// observed fishery length comps
  init_int nyrs_surv3_unbiased_ac			// number of years with  unbiased AI surv age comps
  init_ivector yrs_surv3_unbiased_ac(1,nyrs_surv3_unbiased_ac) 	// vector of index years with unbiased AI surv age comps
  init_matrix oac_surv3_unbiased(1,nyrs_surv3_unbiased_ac,1,nages_dat) 	// observed unbiased AI surv age comps
  init_int nyrs_surv3_lc					// number of years with AI survey length comps
  init_ivector yrs_surv3_lc(1,nyrs_surv3_lc) 		// vector of index years with AI survey length comps
  init_matrix olc_surv3(1,nyrs_surv3_lc,1,nlen) 		// observed AI survey length comps
  init_int nyrs_surv_ebs_unbiased_ac     // number of years with  unbiased EBS slope survey surv age comps
  init_ivector yrs_surv_ebs_unbiased_ac(1,nyrs_surv_ebs_unbiased_ac)  // vector of index years with unbiased EBS slope surv age comps
  init_matrix oac_surv_ebs_unbiased(1,nyrs_surv_ebs_unbiased_ac,1,nages_dat)  // observed unbiased EBS slope surv age comps
  init_int nyrs_surv_ebs_lc          // number of years with EBS survey length comps
  init_ivector yrs_surv_ebs_lc(1,nyrs_surv_ebs_lc)    // vector of index years with EBS survey length comps
  init_matrix olc_surv_ebs(1,nyrs_surv_ebs_lc,1,nlen)     // observed EBS survey length comps
  init_vector wt_pop(1,nages)   // population weight at age
  init_vector wt_fsh(1,nages)   // fishery weights at age
  vector wt_pop_bin(1,nages_dat)    // rescale weights to match the data plus group
  vector wt_fsh_bin(1,nages_dat)    // rescale weights to match the data plus group
  //init_vector maturity(1,nages) // Maturity at age
  //init_vector fecund(1,nages)	// predicted number of eggs by age per spawner; thousands
  //init_vector lvsrv(1,nages)	// predicted larval survival 
  init_number spawn_mo          // spawning month
  init_int fyear_ac_option	// first year age comp option
  init_int historic_catch	// historic catch level
  init_int sr_type		// option for S-R curve (avg or B-H)
  init_int fixedrec_yrs		// the number of years from the endyr in which we fix the recruitments
  init_number sigr 		// the sd of the recruitments 
  init_number priormeansrv3q	// prior mean of trawl survey q
  init_number priorcvsrv3q 	// prior cv of trawl survey q  
  init_vector prop_bio_ai(styr,endyr)  // the proportion of survey biomass in the AI survey area  
  init_number priormean_M	// prior mean of M
  init_number priorcv_M 	// prior cv of M
  init_int    sel_option  // option for selectivity functional form (1=logistic, 2=double logistic, 3= bicubic spline)
  init_int    nbins       // the number of bins
  init_int    sel_styr    // the start year for fitting selecitvity
  //init_int    sel_endyr    // the end year for fitting selectivity
  init_int    sel_fixedyrs    //the number of years from the endyr_r in which we fix the recruitments
  init_ivector binstart(1,nbins)  // the start year for each bin
  init_number sigma_aslp    // sigma for selectivity acsending slope deviations
  init_number sigma_a50     // sigma for selectivity a50 deviations
  init_number sigma_dslp    // sigma for selectivity decsending slope deviations
  init_number sigma_d50     // sigma for selectivity d50 deviations
  init_int    n_yr_nodes      // number of year nodes for bicubic spline selectivity
  init_int    n_age_nodes     // number of age nodes for bicubic spline selectivity
  int         isel_npar       // generalized number of nodes for sel_par
  int         jsel_npar       // generalized number of nodes for sep_par
  //  maturity data
  init_number nages_mat           // number of ages for maturity ogive
  init_vector matages(1,nages_mat) // maturity ages for ogive
  init_number nages_T             // number of ages, TenBrink
  init_vector ages_T(1,nages_T)   // ages, TenBrink
  init_vector T_n(1,nages_T)      // total number of fish observed by age by TenBrink (T_n)
  init_vector T_y(1,nages_T)      // number of mature fish observed by age by TenBrink (T_y)
  init_number nages_S             // number of ages, Shaw
  init_vector ages_S(1,nages_S)   // ages, Shaw
  init_vector S_n(1,nages_S)      // total number of fish observed by age by Shaw (S_n)
  init_vector S_y(1,nages_S)      // number of mature fish observed by age by Shaw (S_y)
  init_vector mat_lambda(1,nages_mat)  // weights for fitting the maturity data 
  !! cout << mat_lambda << endl; 


  // read in 2014 biomass for comparison plot
  init_vector biomass2016(1957,2016)
  !! cout << biomass2016 << endl; 

  //int    nbins                  // number of bins for fishery selectivity parameters
  //int    tmpnbins               // round the number of bins
  //!! cout << " olc_fish is " << olc_fish << endl;
   

  
  number spmo_frac		// spawning month as fraction of the year
  int styr_rec			// start year for recruits
  int styr_rec_dev  // start year for which we have recruitment deviations (some may be fixed)
  int lastyr_rec		// last year for which we estimate the recruitment
  int i                         // index for years 
  int j				// index for ages
  int k
  int l
  int bincount       // for looping through the number of bins                         
  int styr_fut			// Start year for projections
  int endyr_fut                 // End year for projections
  int num_proj_Fs 		// Number of Fs to evaluate in the future
  //vector cv_srv2(1,nyrs_srv2)	// vector of the cv's for the slope survey
  vector cv_srv3(1,nyrs_srv3)	// vector of the cv's for the wedge survey
  vector cv_srv_ebs(1,nyrs_srv_ebs) // vector of the cv's for the EBS  survey
  

  int phase_fydev            // phase for first year deviations
  int sr_phase                 // phase for estiamtion of rzero and steepness
 
  matrix tmp(1,nages,1,nages)   // tmp matrix for rescaling the age error matrix
  matrix tmp2(1,nages,1,nlen)   // tmp matrix for rescaling the transition matrix 

  


// things for getting the survey selectivity a10 (as modified by M)
  int firstage;              // the first age that exceeds 10% survey selection
  int excludeage;            // exclude ages at or below this

  int endyr_r   // end year for retrospective run
  int sel_endyr // end year for estimating selectivity

  int comp_yr_count  // for counting the number of years to be used for compositon data in retrospective run
  int surv_yr_count  // for counting the number of years to be used for survey biomass indices
  int nyrs_fish_unbiased_ac_r  // number of years to use for retrospecitve run  
  int nyrs_surv3_unbiased_ac_r  // number of years to use for retrospecitve run
  int nyrs_surv_ebs_unbiased_ac_r  // number of years to use for retrospecitve run
  int nyrs_fish_lc_r  // number of years to use for retrospecitve run
  int nyrs_surv3_lc_r  // number of years to use for retrospecitve run
  int nyrs_surv_ebs_lc_r  // number of years to use for retrospecitve run
  int nyrs_srv3_r  // number of years to use for retrospecitve run
  int nyrs_srv_ebs_r  // number of years to use for retrospecitve run

 LOCAL_CALCS
  // calculate endyr_r and sel_year for the surveys
  endyr_r = endyr - yrs_r;
  sel_endyr = endyr_r - sel_fixedyrs;

  // cout the number of years of age comps to be used in the retrospective run
  comp_yr_count = 0;
  for (i=1;i<=nyrs_fish_unbiased_ac;i++)
      if (yrs_fish_unbiased_ac(i) <= endyr_r)  comp_yr_count++;  
  nyrs_fish_unbiased_ac_r = comp_yr_count;

  comp_yr_count = 0;
  for (i=1;i<=nyrs_surv3_unbiased_ac;i++)
      if (yrs_surv3_unbiased_ac(i) <= endyr_r)  comp_yr_count++;  
  nyrs_surv3_unbiased_ac_r = comp_yr_count;

  comp_yr_count = 0;
  for (i=1;i<=nyrs_surv_ebs_unbiased_ac;i++)
      if (yrs_surv_ebs_unbiased_ac(i) <= endyr_r)  comp_yr_count++;  
  nyrs_surv_ebs_unbiased_ac_r = comp_yr_count;
  
  comp_yr_count = 0;
  for (i=1;i<=nyrs_fish_lc;i++)
      if (yrs_fish_lc(i) <= endyr_r)  comp_yr_count++;  
  nyrs_fish_lc_r = comp_yr_count;

  comp_yr_count = 0;
  for (i=1;i<=nyrs_surv3_lc;i++)
      if (yrs_surv3_lc(i) <= endyr_r)  comp_yr_count++;  
  nyrs_surv3_lc_r = comp_yr_count;

  comp_yr_count = 0;
  for (i=1;i<=nyrs_surv_ebs_lc;i++)
      if (yrs_surv_ebs_lc(i) <= endyr_r)  comp_yr_count++;  
  nyrs_surv_ebs_lc_r = comp_yr_count;

  // cout the number of years of survey indices to be used in the retrospective run
  surv_yr_count = 0;
  for (i=1;i<=nyrs_srv3;i++)
      if (yrs_srv3(i) <= endyr_r)  surv_yr_count++;  
  nyrs_srv3_r = surv_yr_count;

  surv_yr_count = 0;
  for (i=1;i<=nyrs_srv_ebs;i++)
      if (yrs_srv_ebs(i) <= endyr_r)  surv_yr_count++;  
  nyrs_srv_ebs_r = surv_yr_count;

 END_CALCS 
 
  matrix rescaled_sel_fish(styr,endyr_r,1,nages_dat)  // rescaled fishery selectivity matrix
  vector rescaled_F(styr,endyr_r)                 // rescaled F values
  vector recent_fish_sel(1,nages_dat)               // the recent fish selectivity (for spr calcs)

  vector yrs_fish_unbiased_ac_r(1,nyrs_fish_unbiased_ac_r)  // years to use for retrospecitve run  
  vector yrs_surv3_unbiased_ac_r(1,nyrs_surv3_unbiased_ac_r)  // years to use for retrospecitve run
  vector yrs_surv_ebs_unbiased_ac_r(1,nyrs_surv_ebs_unbiased_ac_r)  // years to use for retrospecitve run
  vector yrs_fish_lc_r(1,nyrs_fish_lc_r)  // years to use for retrospecitve run
  vector yrs_surv3_lc_r(1,nyrs_surv3_lc_r)  // years to use for retrospecitve run
  vector yrs_surv_ebs_lc_r(1,nyrs_surv_ebs_lc_r)  // years to use for retrospecitve run
  matrix oac_fish_unbiased_r(1,nyrs_fish_unbiased_ac_r,1,nages_dat)  // observed unbiased fishery age comps, retrospective run
  matrix olc_fish_r(1,nyrs_fish_lc_r,1,nlen)     // observed fishery length comps, retrospective run
  matrix oac_surv3_unbiased_r(1,nyrs_surv3_unbiased_ac_r,1,nages_dat)   // observed unbiased AI surv age comps, retrospective run
  matrix olc_surv3_r(1,nyrs_surv3_lc_r,1,nlen)     // observed AI survey length comps, retrospective run
  matrix oac_surv_ebs_unbiased_r(1,nyrs_surv_ebs_unbiased_ac_r,1,nages_dat)  // observed unbiased EBS slope surv age comps, retrospective run
  matrix olc_surv_ebs_r(1,nyrs_surv_ebs_lc_r,1,nlen)     // observed EBS survey length comps, retrospective run 



 LOCAL_CALCS
  // get the comp data and years for the retrospective run
  yrs_fish_unbiased_ac_r = yrs_fish_unbiased_ac(1,nyrs_fish_unbiased_ac_r);
  yrs_surv3_unbiased_ac_r = yrs_surv3_unbiased_ac(1,nyrs_surv3_unbiased_ac_r);
  yrs_surv_ebs_unbiased_ac_r = yrs_surv_ebs_unbiased_ac(1,nyrs_surv_ebs_unbiased_ac_r);
  yrs_fish_lc_r = yrs_fish_lc(1,nyrs_fish_lc_r);
  yrs_surv3_lc_r = yrs_surv3_lc(1,nyrs_surv3_lc_r);
  yrs_surv_ebs_lc_r = yrs_surv_ebs_lc(1,nyrs_surv_ebs_lc_r); 
  for (i=1;i<=nyrs_fish_unbiased_ac_r;i++)
     oac_fish_unbiased_r(i) = oac_fish_unbiased(i); 
  for (i=1;i<=nyrs_surv3_unbiased_ac_r;i++)
     oac_surv3_unbiased_r(i) = oac_surv3_unbiased(i); 
  for (i=1;i<=nyrs_surv_ebs_unbiased_ac_r;i++)
      oac_surv_ebs_unbiased_r(i) =  oac_surv_ebs_unbiased(i); 
  for (i=1;i<=nyrs_fish_lc_r;i++)
      olc_fish_r(i) =  olc_fish(i);  
  for (i=1;i<=nyrs_surv3_lc_r;i++)
      olc_surv3_r(i) =  olc_surv3(i);  
  for (i=1;i<=nyrs_surv_ebs_lc_r;i++)
      olc_surv_ebs_r(i) = olc_surv_ebs(i); 

  // bin the maturity and weight to match data plus group
  wt_pop_bin(1,nages_dat-1) = wt_pop(1,nages_dat-1);
  wt_pop_bin(nages_dat) = mean(wt_pop(nages_dat,nages)); 
  wt_fsh_bin(1,nages_dat-1) = wt_fsh(1,nages_dat-1);
  wt_fsh_bin(nages_dat) = mean(wt_fsh(nages_dat,nages)); 
  spmo_frac = (spawn_mo-1)/12.;
  num_proj_Fs = 5;
  styr_fut=endyr_r+1;
  endyr_fut = styr_fut+10;
  // define styr_rec +++++++++++++++++++++++++++++++++++++++

  if (fyear_ac_option == 1) // first year rec are combined with other recruitments  
  {
   styr_rec = styr-nages+1;
   styr_rec_dev = styr-nages_dat+1;    // some cohorts share a recruitment deviation
   phase_fydev = -1;
  }
  else if (fyear_ac_option == 2) // first year recruitments are in equilibrium with historic catch
  { 
   styr_rec = styr;
   styr_rec_dev = styr;
   phase_fydev = -1;
  }
  else if (fyear_ac_option == 3) // first year recruitments are stochastic, but separate from other recs
  { 
   styr_rec = styr;
   styr_rec_dev = styr;
   phase_fydev = 3;
  }  
  lastyr_rec = endyr_r - fixedrec_yrs;   // define the last year for which we estimate recruitment
  
  if (sr_type==1) sr_phase=-1;
    else sr_phase = 2; 

   tmp = trans(unbiasedages);
   for (i=1;i<=nages;i++) tmp(i) = tmp(i)/sum(tmp(i));
   unbiasedages = trans(tmp);
  
   tmp2 = trans(translen);
   for (i=1;i<=nages;i++) tmp2(i) = tmp2(i)/sum(tmp2(i));
   translen = trans(tmp2);

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // calculate cv's for the survey
  //cv_srv2=elem_div(obs_srv2_sd,obs_srv2);
  cv_srv3=elem_div(obs_srv3_sd,obs_srv3);
  cv_srv_ebs=elem_div(obs_srv_ebs_sd,obs_srv_ebs);
  
  // compute the number of bins
  //tmpnbins = (endyr - styr +1)/fselbinsize;
  //if ((endyr - styr +1)/fselbinsize - tmpnbins < 1e-3 ) nbins=tmpnbins;
  //else nbins=tmpnbins+1;
  //cout << "look here "<<endl;
  //cout << endyr - styr +1 <<" "<<fselbinsize<<" "<< tmpnbins <<" " <<nbins<<endl;  

  // start to read from the control file
  ad_comm::change_datafile_name("pop18.ctl");
 END_CALCS
  init_int phase_selcoff
  init_int phase_logist_sel
  init_int phase_f_sel_param
  init_int phase_s_sel_param_ai
  init_int phase_s_sel_param_ebs
  init_int phase_q_srv_ebs  
  init_int phase_srec
  init_int phase_f40
  init_int phase_proj
  init_int phase_historic_F
  init_int har_flag
  //init_vector fish_biased_ac_samp(1,nyrs_fish_biased_ac) //***
  init_vector raw_fish_unbiased_ac_samp(1,nyrs_fish_unbiased_ac)
  !! cout << " nyrs_fish_unbiased_ac " << nyrs_fish_unbiased_ac << endl;
  !! cout << " raw_fish_unbiased_ac_samp " << raw_fish_unbiased_ac_samp << endl;
  init_vector raw_fish_lc_samp(1,nyrs_fish_lc)
  !! cout << " raw_fish_lc_samp " << raw_fish_lc_samp << endl;
  init_vector raw_surv3_unbiased_ac_samp(1,nyrs_surv3_unbiased_ac)
  init_vector raw_surv3_lc_samp(1,nyrs_surv3_lc)
  init_vector raw_surv_ebs_unbiased_ac_samp(1,nyrs_surv_ebs_unbiased_ac)
  init_vector raw_surv_ebs_lc_samp(1,nyrs_surv_ebs_lc)

  init_vector lambda(1,15)
  !! cout << lambda << endl;

  // vectors for sample sizes for retrospective runs
  vector fish_unbiased_ac_samp_r(1,nyrs_fish_unbiased_ac_r)
  vector fish_lc_samp_r(1,nyrs_fish_lc_r)
  vector surv3_unbiased_ac_samp_r(1,nyrs_surv3_unbiased_ac_r)   
  vector surv3_lc_samp_r(1,nyrs_surv3_lc_r)
  vector surv_ebs_unbiased_ac_samp_r(1,nyrs_surv_ebs_unbiased_ac_r)
  vector surv_ebs_lc_samp_r(1,nyrs_surv_ebs_lc_r)




  LOCAL_CALCS 
  // start to read from the composition weight file
  // order is fac, flc, sac, slc
  ad_comm::change_datafile_name("compweights.ctl");
 END_CALCS
  init_vector compweights(1,6)
 
 LOCAL_CALCS 
    /*
    fish_unbiased_ac_samp = compweights(1)*raw_fish_unbiased_ac_samp;
    fish_lc_samp = compweights(2)*raw_fish_lc_samp;
    surv3_unbiased_ac_samp = compweights(3)*raw_surv3_unbiased_ac_samp;
    surv3_lc_samp = compweights(4)*raw_surv3_lc_samp;
    surv_ebs_unbiased_ac_samp = compweights(5)*raw_surv_ebs_unbiased_ac_samp; 
    surv_ebs_lc_samp = compweights(6)*raw_surv_ebs_lc_samp;
    */    
    fish_unbiased_ac_samp_r = compweights(1)*raw_fish_unbiased_ac_samp(1,nyrs_fish_unbiased_ac_r);
    fish_lc_samp_r = compweights(2)*raw_fish_lc_samp(1,nyrs_fish_lc_r);
    surv3_unbiased_ac_samp_r = compweights(3)*raw_surv3_unbiased_ac_samp(1,nyrs_surv3_unbiased_ac_r);
    surv3_lc_samp_r = compweights(4)*raw_surv3_lc_samp(1,nyrs_surv3_lc_r);
    surv_ebs_unbiased_ac_samp_r = compweights(5)*raw_surv_ebs_unbiased_ac_samp(1,nyrs_surv_ebs_unbiased_ac_r);
    surv_ebs_lc_samp_r = compweights(6)*raw_surv_ebs_lc_samp(1,nyrs_surv_ebs_lc_r);

 END_CALCS 
 

  int phase_f_sel_ascend;
  int phase_f_sel_descend;
  int phase_f_sel_par;
  int phase_a50_devs;
  int phase_aslope_devs;
  int phase_d50_devs;
  int phase_dslope_devs;
  //int phase_cpue;
  ivector binindex(sel_styr,sel_endyr);  
  //ivector binstart(1,nbins);           // the start year for each time bin
  vector scal_yr_nodes(1,n_yr_nodes);  // the yr nodes scaled from 0 to 1
  vector scal_age_nodes(1,n_age_nodes);  // the age nodes scaled from 0 to 1
 
 LOCAL_CALCS  // set phases for fishery and survey selectivity options, and time-varying deviations
  //  phase for the ascending parameters for logistic
  //  phase for the descending parameters for the logistic
  //  phase for the age and year nodes for the bicubic spline
  //  phase for ebs catchabilty and selectivity (turn off if they are not being fit)
  //  set phase to cpue to zero if not being fit
  
  // compute the years each bin starts
    //binstart(1) = styr;
    //for (i=2;i<=nbins;i++) 
    //    binstart(i) = binstart(i-1) + fselbinsize;

  // compute binindex for time-varying selectivity
    if (nbins > 1)
    {
      bincount = 1;
      for (i=sel_styr; i<=sel_endyr; i++)
          {
            if (i == binstart(bincount)) 
            {
             binindex(i) = bincount; 
             if( bincount < nbins ) bincount++;   
            }
            if (i < sel_endyr) binindex(i+1) = binindex(i);
          }
           cout << binindex << endl;
    



    // for (i=styr;i<=endyr;i++)
    // {
    //  k = (i - styr +1)/fselbinsize;
    //  if ((i - styr +1)/fselbinsize - k < 1e-3 ) 
    //    binindex(i)=k;
    //  else binindex(i)=k+1;
    // }
    }


     

    


  // set phases to -1, and later turn on the ones we need; 
  //   selectivity parameters
  phase_f_sel_ascend = -1;
  phase_f_sel_descend = -1;
  phase_f_sel_par = -1;
  //   time-varying parameters for logistic and double logistic
  phase_a50_devs = -1;
  phase_aslope_devs = -1;
  phase_d50_devs = -1;
  phase_dslope_devs = -1;
  if (sel_option==1)  // logistic fishery selectivity 
  {
    phase_f_sel_ascend = phase_f_sel_param;
    if(nbins > 1)
      {
        phase_a50_devs = 4;
        phase_aslope_devs = 4;
      }

  }
  else if (sel_option==2)  // double logistic fishery selectivity 
  {
    phase_f_sel_ascend = phase_f_sel_param;
    phase_f_sel_descend = phase_f_sel_param;
    if(nbins > 1)
      {
        phase_a50_devs = 4;
        phase_aslope_devs = 4;
        phase_d50_devs = 4;
        phase_dslope_devs = 4;
      }

  }
  else if (sel_option==3)  // bicubic fishery selectivity 
  {
    phase_f_sel_par = phase_f_sel_param;
    scal_age_nodes.fill_seqadd(0,1./(n_age_nodes-1));
    scal_yr_nodes.fill_seqadd(0,1./(n_yr_nodes-1));
    isel_npar = n_yr_nodes;
    jsel_npar = n_age_nodes; 
  }

  else if (sel_option==4)  // time-varying cubic fishery selectivity 
  {
    phase_f_sel_par = phase_f_sel_param;
    isel_npar = n_age_nodes;
    jsel_npar = nbins; 
  }

  if (lambda(9)==0 & lambda(10)==0 )  // turn off ebs srv sel paramaters if the age and length comps are not being fit 
    phase_s_sel_param_ebs = -1;
  if (lambda(15)==0)  // turn off ebs srv catchability if it is not being fit
    phase_q_srv_ebs = -1;

    cout <<" phase_s_sel_param_ebs is "<<phase_s_sel_param_ebs<< endl; 
   cout <<" phase_q_srv_ebs is "<<phase_q_srv_ebs<< endl;    
  //if (lambda(3)==0) phase_cpue = -1;
  //   else phase_cpue = 1;

 END_CALCS

 

  

INITIALIZATION_SECTION //-------------------------------------------------------------------------------------
  logM    -2.80836266400
  //logq_srv3  0.5 
  logq_srv3 0.153074980052
  logq_srv_ebs 0.153074980052
  //F40 0.044
  //F35 0.053
  //F30 0.063
  mat_beta1 -6.61184
  mat_beta2 0.727016
  //mean_log_rec 3.97
  mean_log_rec  4.25988340015  
  log_rinit  3.89901668107
  log_rzero  4.169
  log_avg_fmort -3.03702074878
  sel_aslope_domfish 1.214338
  sel_dslope_domfish 0.0
  sel_aslope_forfish 1.63656109222
  sel_dslope_forfish -1.63 //0.0
  sel_aslope_srv3 1.01826089348
  sel_aslope_srv_ebs 1.01826089348
  sel_a50_domfish 6.574314
  sel_d50_domfish 24.903941
  //sel_a50_forfish 6.9 
  sel_a50_forfish 7.32882926631
  sel_d50_forfish 24.903941
  //sel_a50_srv3 5.4
  sel_a50_srv3 6.05353089120
  sel_a50_srv_ebs 6.05353089120
  
  steepness 1.0
  //sigr 0.8013
  historic_F 0.000
  //logq_cpue   -4.58 
  //logq_cpue -4.67827481754
  //cv_cpue 0.3
  sel_par 0
   
 
  

PARAMETER_SECTION //-----------------------------------------------------------------------------------------

 // offset parameters
 vector offset(1,6)
 // selectivity parameters 
 //  First, the logistic curve parameters for the domestic and foreign fisheries(slope and 50% parameters)
 init_number sel_aslope_domfish(-1)
 init_number sel_dslope_domfish(-1)
 //init_number sel_aslope_forfish(phase_logist_param)
 init_number sel_aslope_forfish(phase_f_sel_ascend)
 init_number sel_dslope_forfish(phase_f_sel_descend)
 init_number sel_a50_domfish(-1) 
 init_number sel_d50_domfish(-1)
 //init_number sel_a50_forfish(phase_logist_param)
 init_number sel_a50_forfish(phase_f_sel_ascend)
 init_number sel_d50_forfish(phase_f_sel_descend)
 // fishery selectivity matrix for bicubic and time-varying cubic selectivity 
 init_matrix sel_par(1,jsel_npar,1,isel_npar,phase_f_sel_par);
 
 //init_number sel_d50_forfish(-1)
 // now the logistic curve parameters for the AI trawl survey
 //init_bounded_number sel_aslope_srv3(0.1,3.0,3)
 //init_bounded_number sel_a50_srv3(0.1,20.0,3)
 init_bounded_number sel_aslope_srv3(0.1,3.0,phase_s_sel_param_ai)
 init_bounded_number sel_a50_srv3(0.1,30.0,phase_s_sel_param_ai)

 // now the logistic curve parameters for the EBS trawl survey
 init_bounded_number sel_aslope_srv_ebs(0.1,3.0,phase_s_sel_param_ebs)
 init_bounded_number sel_a50_srv_ebs(0.1,30.0,phase_s_sel_param_ebs)

 
 // deviations on fishery selectivity
 //init_bounded_dev_vector a50_devs(1,nbins,-10.,10.,4)
 //init_bounded_dev_vector aslope_devs(1,nbins,-10.,10.,4)
 //number slptmp
 //number a50tmp
 init_bounded_dev_vector a50_devs(1,nbins,-10.,10.,phase_a50_devs)
 init_bounded_dev_vector aslope_devs(1,nbins,-10.,10.,phase_aslope_devs)
 init_bounded_dev_vector d50_devs(1,nbins,-10.,10.,phase_d50_devs)
 init_bounded_dev_vector dslope_devs(1,nbins,-10.,10.,phase_dslope_devs)
 
 number a_slptmp
 number a50tmp
 number d_slptmp
 number d50tmp
 
 // Or we could estimate the selectivity coefficients directly (not used so far)
 //init_vector log_selcoffs_domfish(1,nselages,phase_selcoff)
 //init_vector log_selcoffs_forfish(1,nselages,phase_selcoff)
 //init_vector log_selcoffs_srv2(1,nselages,phase_selcoff)
 //number avgsel_fish;  // (averages are used in the likelihood function (if used))
 //number avgsel_srv2;

 // Next, the selectivity values (logged and unlogged) 
 matrix log_sel_fish(sel_styr,sel_endyr,1,nselages)
 vector log_sel_srv3(1,nselages)
 vector log_sel_srv_ebs(1,nselages)    // EBS slope survey
 matrix sel_fish(styr,endyr_r,1,nages)
 vector sel_srv3(1,nages)
 vector sel_srv_ebs(1,nages)
 

 // The survival parameters
 init_number logM(4)
 number M
 number surv
 init_number log_avg_fmort(1)
 init_bounded_vector fmort_dev(styr_fish,endyr_r,-10,10,2)
 number avg_fmort_dev
 matrix F(styr_rec,endyr_r,1,nages)
 matrix Z(styr_rec,endyr_r,1,nages)
 matrix S(styr_rec,endyr_r,1,nages)
 matrix mort(styr_rec,endyr_r,1,nages);  // The multiplier for the mean population size
 matrix spawner_S(styr_rec,endyr_r,1,nages)
 


 // The numbers at age parameters
 init_bounded_dev_vector rec_dev(styr_rec_dev,lastyr_rec,-10,10,2) // recruitment deviations (styr to endyr)
 init_number mean_log_rec(1)			// mean recruitment
 init_bounded_number log_rinit(0,10,2)                       // initial recruitment to start model
 matrix natage(styr_rec,endyr_r+1,1,nages)		// numbers at age
 matrix natage_bin(styr_rec,endyr_r+1,1,nages_dat)    // numbers at age
 vector natagetmp(1,nages)   // temporary numbers at age (prior to start year)
 init_bounded_vector fydev(2,nages_dat,-10,10,phase_fydev)    // deviations around eq numbers in first year (stage 3)
 
 //vector spawners(styr,endyr) 			// estimated biomass of female spawners
 //vector expbiom(styr,endyr) 			// estimated exploitable biomass
 sdreport_vector totbiom(styr_rec-rec_age,endyr_r+1)		// total biomass of population
 //sdreport_number depletion			// change in totbiom over time
 //sdreport_number endbiom			// totbiom in final year
 init_number historic_F(phase_historic_F)	// historic F for computing first year age comps
   
// The parameters for evaluating the objective function
 objective_function_value obj_fun
 vector rec_like(1,3)
 vector surv_like(1,2)
 number catch_like
 number fpen
 number hf_pen
 vector age_like(1,6)
 matrix effn(1,6,1,40)
 vector rmse(1,3)  // residual mean square error for the surveys and recruitment
 vector sel_like(1,9)
 number sprpen
 vector priormq(1,2)
 vector sdnr(1,7)
 vector al_rmse(1,6)
 number mat_like

// The parameters for getting the predicted values for the CPUE index (srv1)
  //init_number logq_cpue(phase_cpue)			// q for cpue
  //number q_cpue
  //vector pred_cpue(1,nyrs_cpue)			// predicted cpue (survey 1)
  //init_number cv_cpue(-1)			// cv for cpue index  

// The parameters for getting the predicted values for the AI and EBS trawl survey index (srv3)
  init_number logq_srv3(4)		// trawl survey  catchability
  number   q_srv3
  vector pred_srv3(styr_fish,endyr_r)
  init_number logq_srv_ebs(phase_q_srv_ebs)    // trawl survey  catchability
  number   q_srv_ebs
  vector pred_srv_ebs(styr_fish,endyr_r)


// The parameters for getting the predicted values for the catch total
  vector pred_catch(styr_fish,endyr_r)  
  number ehc		   			// the estimated historic catch
  matrix catage(styr_fish,endyr_r,1,nages)

// The parameters for getting the predicted age comps
  matrix eac_fish_unbiased_mod(styr_fish,endyr_r,1,nages)
  matrix eac_surv3_unbiased_mod(styr,endyr_r,1,nages)
  matrix eac_surv_ebs_unbiased_mod(styr,endyr_r,1,nages)
  matrix eac_fish_unbiased_dat(styr_fish,endyr_r,1,nages_dat)
  matrix eac_surv3_unbiased_dat(styr,endyr_r,1,nages_dat)
  matrix eac_surv_ebs_unbiased_dat(styr,endyr_r,1,nages_dat)
 // matrix eac_fish_biased(styr_fish,endyr,1,nages)
  matrix elc_fish(styr_fish,endyr_r,1,nlen)
  matrix elc_fish_tempages(styr_fish,endyr_r,1,nages)  // the true ages for the year in which we have length comps
  matrix elc_surv3(styr,endyr_r,1,nlen) 
  matrix elc_surv3_tempages(styr,endyr_r,1,nages)  // AI survey, the true ages for the year in which we have length comps
  matrix elc_surv_ebs(styr,endyr_r,1,nlen) 
  matrix elc_surv_ebs_tempages(styr,endyr_r,1,nages)  // EBS survey, the true ages for the year in which we have length comps



// The parameters for getting the SPR values
   number F40
   number F35
   number F30
   number SB0
   number SBF40
   number SBF35
   number SBF30
   //matrix Nspr(1,4,1,nages)
   //number phizero

// numbers for getting the normalized residuals
 vector  survey_nr(1,nyrs_srv3_r)  // trawl survey normalized residuals
 //vector  cpue_nr(1,nyrs_cpue)  //  cpue normalized residuals
 //vector  fbac_nr(1,nyrs_fish_biased_ac*nages)   // fishery biased age comps
      // the normalized residuals from the multinomial, by bin (method TA1.2 in Francis 2011)
 vector  fac_nr(1,nyrs_fish_unbiased_ac_r*nages)   // fishery age comps
 vector  flc_nr(1,nyrs_fish_lc_r*nlen) 		   // fishery length comps
 vector  sac_nr(1,nyrs_surv3_unbiased_ac_r*nages)  // survey age comps
 vector  slc_nr(1,nyrs_surv3_lc_r*nlen) 	   	// survey length comps
 vector  ebs_sac_nr(1,nyrs_surv_ebs_unbiased_ac_r*nages)       // EBS survey age comps
 vector  ebs_slc_nr(1,nyrs_surv_ebs_lc_r*nlen)       // ebs survey length comps

 vector  fac_mcian_wgt(1,nyrs_fish_unbiased_ac_r)   // McAllister-Ianelli weights (method TA1.1 in Francis 2011)
 vector  flc_mcian_wgt(1,nyrs_fish_lc_r)   
 vector  sac_mcian_wgt(1,nyrs_surv3_unbiased_ac_r)
 vector  slc_mcian_wgt(1,nyrs_surv3_lc_r)
 vector  ebs_sac_mcian_wgt(1,nyrs_surv_ebs_unbiased_ac_r)
 vector  ebs_slc_mcian_wgt(1,nyrs_surv_ebs_lc_r)

 vector  fac_mcian_wgt_inv(1,nyrs_fish_unbiased_ac_r)   // inverse of McAllister-Ianelli weights for harmonic mean
 vector  flc_mcian_wgt_inv(1,nyrs_fish_lc_r)   
 vector  sac_mcian_wgt_inv(1,nyrs_surv3_unbiased_ac_r)
 vector  slc_mcian_wgt_inv(1,nyrs_surv3_lc_r)
 vector  ebs_sac_mcian_wgt_inv(1,nyrs_surv_ebs_unbiased_ac_r)
 vector  ebs_slc_mcian_wgt_inv(1,nyrs_surv_ebs_lc_r)

 vector fac_nr_fran(1,nyrs_fish_unbiased_ac_r)       // normalized residuals from the Francis method (method TA1.8 In Francis 2011)
 vector flc_nr_fran(1,nyrs_fish_lc_r)       
 vector sac_nr_fran(1,nyrs_surv3_unbiased_ac_r)       
 vector slc_nr_fran(1,nyrs_surv3_lc_r)
 vector ebs_sac_nr_fran(1,nyrs_surv_ebs_unbiased_ac_r)
 vector ebs_slc_nr_fran(1,nyrs_surv_ebs_lc_r)

 //vector  fbac_nr(1,nyrs_fish_biased_ac*nages)   // fishery biased age comps
 vector  fac_resid(1,nyrs_fish_unbiased_ac_r*nages)   // fishery age comps
 vector  flc_resid(1,nyrs_fish_lc_r*nlen)       // fishery length comps
 vector  sac_resid(1,nyrs_surv3_unbiased_ac_r*nages)  // survey age comps
 vector  slc_resid(1,nyrs_surv3_lc_r*nlen)       // survey length comps
 vector  ebs_sac_resid(1,nyrs_surv_ebs_unbiased_ac_r*nages)  // ebs survey age comps
 vector  ebs_slc_resid(1,nyrs_surv_ebs_lc_r*nlen)       // EBS survey length comps

// Parameters for doing the future projections
 //matrix nage_future(styr_fut,endyr_fut,1,nages)
 //number ftmp
 //number mean_recent_fs
 //matrix Z_future(styr_fut,endyr_fut,1,nages)
 //matrix F_future(styr_fut,endyr_fut,1,nages)
 //matrix S_future(styr_fut,endyr_fut,1,nages)
 //init_vector rec_dev_future(styr_fut,endyr_fut,phase_proj)
 //matrix catage_future(styr_fut,endyr_fut,1,nages)
 //sdreport_matrix catch_future(1,num_proj_Fs-1,styr_fut,endyr_fut)
 //sdreport_matrix biomass_future(1,num_proj_Fs,styr_fut,endyr_fut)
 //sdreport_matrix ssb_future(1,num_proj_Fs,styr_fut,endyr_fut)
 
// Stock recruitment params 
   init_number log_rzero(sr_phase)
   number bzero
   number rzero
   sdreport_vector sp_biom(styr_rec-rec_age,endyr_r)
   vector eggs(styr_rec-rec_age,endyr_r)
   vector larvae(styr_rec-rec_age,endyr_r)
   sdreport_vector est_rec(styr_rec_dev,lastyr_rec)
   number alpha   
   number beta
   init_bounded_number steepness(0.2001,1.0,sr_phase)
   vector pred_rec(styr_rec_dev,lastyr_rec)
   vector est_spb(styr_rec_dev,lastyr_rec)
   vector chi(styr_rec_dev,lastyr_rec)   // the squared difference between est and pred recruits
   //init_bounded_number sigr(0.1,2.0,2)  // the sigma parameter for the recruit likelihood
   number sumrecdev // the sum of the lognormal deviations for the recruitments
   vector SRec_spawn(1,20) // data for estimated SR curve
   vector SRec_rec(1,20) //  data for estimated SR curve

   vector xdum2(styr,endyr_r)  // for the cxx file

  // vector for projection data file
   vector Mvec(1,nages)  // natural mortality, repeated for each age

// maturity estimation
   init_bounded_number mat_beta1(-10,2,5)   // beta0 and beta1 parameters for logistic regression of maturity curve
   init_bounded_number mat_beta2(0,2,5)
   vector mat_theta(1,nages_mat)   // theta parameter estimates for logistic regression of maturity curve (uses all  ages for which we have data)
   vector maturity(1,nages) // the maturity for the ages used in the population model
   sdreport_vector maturity_bin(1,nages_dat) // maturity ogive, binned to match the data plus group
       
// updated compweights
   vector compweightsnew_ta12(1,6)     // weight by inverse of variance of normalized resid
   vector compweightsnew_ta11(1,6)     //  McAllister-Ianelli method
   vector compweightsnew_ta18(1,6)     // the Francis method

PRELIMINARY_CALCS_SECTION //-------------------------------------------------------------------------------
  // Compute the offsets for the multinomial distributions

// for ( i=1;i<=nyrs_fish_biased_ac; i++)  // fishery biased age comps
// {
// 	oac_fish_biased(i) = oac_fish_biased(i)/sum(oac_fish_biased(i)); // make sure age comps add to 1.0 for each year
// 	offset(1)-=fish_biased_ac_samp(i)*(oac_fish_biased(i))*log(1.e-13+oac_fish_biased(i)); //get the negative log like
// }
  for ( i=1;i<=nyrs_fish_unbiased_ac_r; i++)  // fishery unbiased age comps
 {
 	oac_fish_unbiased_r(i) = oac_fish_unbiased_r(i)/sum(oac_fish_unbiased_r(i)); // make sure age comps add to 1.0 for each year
 	offset(1)-=fish_unbiased_ac_samp_r(i)*(oac_fish_unbiased_r(i))*log(1.e-13+oac_fish_unbiased_r(i)); //get the negative log like
 }
 for ( i=1;i<=nyrs_fish_lc_r; i++)  // fishery length comps
 {
 	olc_fish_r(i) = olc_fish_r(i)/sum(olc_fish_r(i)); // make sure age comps add to 1.0 for each year
 	offset(2)-=fish_lc_samp_r(i)*(olc_fish_r(i))*log(1.e-13+olc_fish_r(i)); //get the negative log like
 }
  
 for ( i=1;i<=nyrs_surv3_unbiased_ac_r; i++)  // AI survey age comps
 {
 	oac_surv3_unbiased_r(i) = oac_surv3_unbiased_r(i)/sum(oac_surv3_unbiased_r(i)); // make sure age comps add to 1.0 for each year
 	offset(3)-=surv3_unbiased_ac_samp_r(i)*(oac_surv3_unbiased_r(i))*log(1.e-13+oac_surv3_unbiased_r(i)); //get the negative log like
 }

 for ( i=1;i<=nyrs_surv3_lc_r; i++)  // AI survey length comps
 {
 	olc_surv3(i) = olc_surv3_r(i)/sum(olc_surv3_r(i)); // make sure age comps add to 1.0 for each year
 	offset(4)-=surv3_lc_samp_r(i)*(olc_surv3_r(i))*log(1.e-13+olc_surv3_r(i)); //get the negative log like
 } 
 
 for ( i=1;i<=nyrs_surv_ebs_unbiased_ac_r; i++)  // EBS survey age comps
 {
  oac_surv_ebs_unbiased_r(i) = oac_surv_ebs_unbiased_r(i)/sum(oac_surv_ebs_unbiased_r(i)); // make sure age comps add to 1.0 for each year
  offset(5)-=surv_ebs_unbiased_ac_samp_r(i)*(oac_surv_ebs_unbiased_r(i))*log(1.e-13+oac_surv_ebs_unbiased_r(i)); //get the negative log like
 } 

 for ( i=1;i<=nyrs_surv_ebs_lc_r; i++)  // EBS survey length comps
 {
  olc_surv_ebs_r(i) = olc_surv_ebs_r(i)/sum(olc_surv_ebs_r(i)); // make sure age comps add to 1.0 for each year
  offset(6)-=surv_ebs_lc_samp_r(i)*(olc_surv_ebs_r(i))*log(1.e-13+olc_surv_ebs_r(i)); //get the negative log like
 } 


PROCEDURE_SECTION //-----------------------------------------------------------------------------------------
// example of using FUNCTION to structure the procedure section
    get_maturity();
    //cout << " got mat" << endl;
    get_selectivity();
    //cout << " got sel" << endl;
    get_mortality();
    //cout << "got mort"  << endl;
    get_first_year();
    //cout <<" got f year" << endl;
    get_numbers_at_age();
    //cout <<" got natage" << endl;
    get_expected_values();
    //cout <<" got expected values" << endl;
    get_sr_inputs();
    //cout <<" got sr inputs" << endl;
    get_catch_at_age();
    //cout <<" got catage" << endl;
    get_age_comps();
    //cout <<" got age comps" << endl;
    get_binned();
    //cout << " get binned"  << endl;
    
    //if (active(F40))
    //  compute_spr_rates();
//    if (current_phase()>6)
//      future_projections();

 //    if (current_phase()>=phase_f40){
 //      cout << " in procedure " << endl;
 //     get_Fspr();
 //   }
    evaluate_the_objective_function();
    //cout << " get obj"  << endl;

    get_age10();
    //cout << " get age10"  << endl;
    update_compweights();
    //cout << " get upfate comp things"  << endl;

   
     if (mceval_phase())
  {
    ofstream evalout("evalout_pop.prj", ios::app);
    evalout <<obj_fun<<" "<<M<<" "<<q_srv3<<" "<<mean_log_rec<<" "<<totbiom<<" "<<sp_biom<<" "<<est_rec<<" "<<maturity<<endl;
  }

  //cout <<" obj_fun is " <<obj_fun <<  endl;
  //cout <<" mat_beta  is "<< mat_beta << endl;
  
 
FUNCTION get_maturity   // get the maturity ogive
 mat_theta = elem_div(mfexp(mat_beta1 + mat_beta2*matages),(1. + mfexp(mat_beta1 + mat_beta2*matages)));
 maturity = mat_theta(1,nages);
 maturity(nages) = 0.5*(maturity(nages) +  mfexp(mat_beta1 + mat_beta2*100)/(1. + mfexp(mat_beta1 + mat_beta2*100)));

FUNCTION get_selectivity //---------------------------------------------------------------------------------
 // Calculate the logistic selectivity (only if being used)
 if (current_phase()>=phase_logist_sel)
 {
 
  for (j=1;j<=nselages;j++)  // Get the AI and EBS survey selectivty
   {
     log_sel_srv3(j) = -1.*log((1.0+mfexp(-1.0*sel_aslope_srv3*(ages(j)-sel_a50_srv3))));
     log_sel_srv_ebs(j) = -1.*log((1.0+mfexp(-1.0*sel_aslope_srv_ebs*(ages(j)-sel_a50_srv_ebs))));
   }
  

  if (sel_option==1)
 { 
 for (i=sel_styr; i<=sel_endyr; i++)  // Use the selectivity values for the foreign fishery through 1988
   {
     if (nbins==1) 
     {
      for (j=1;j<=nselages;j++)
        {
          log_sel_fish(i,j) = -1.*log((1.0+mfexp(-1.0*sel_aslope_forfish*(ages(j)-sel_a50_forfish))));
        }
     }
    else if (nbins>1)
    {
     a_slptmp = sel_aslope_forfish*exp(aslope_devs(binindex(i)));
     a50tmp = sel_a50_forfish*exp(a50_devs(binindex(i)));
     for (j=1;j<=nselages;j++)
       {
        log_sel_fish(i,j) = -1.*log((1.0+mfexp(-1.0*a_slptmp*(ages(j)-a50tmp))));
       }

    }
   }
 
 }
 
 else if (sel_option==2)
 {
  for (i=sel_styr; i<=sel_endyr; i++)  // Use the selectivity values for the foreign fishery through 1988 
  {
    if (nbins==1) 
    {
     for (j=1;j<=nselages;j++)
       {
        log_sel_fish(i,j) = -1.*log(1.0+mfexp(-1.0*sel_aslope_forfish*(ages(j)-sel_a50_forfish)))
                 -1.*log(1.0+mfexp(-1.0*sel_dslope_forfish*(ages(j)-sel_d50_forfish)));
       }
    }
    else if (nbins>1)
    {
     a_slptmp = sel_aslope_forfish*exp(aslope_devs(binindex(i)));
     a50tmp = sel_a50_forfish*exp(a50_devs(binindex(i)));
     d_slptmp = sel_dslope_forfish*exp(dslope_devs(binindex(i)));
     d50tmp = sel_d50_forfish*exp(d50_devs(binindex(i)));

     for (j=1;j<=nselages;j++)
       {
         log_sel_fish(i,j) = -1.*log(1.0+mfexp(-1.0*a_slptmp*(ages(j)-a50tmp)))
                   -1.*log(1.0+mfexp(-1.0*d_slptmp*(ages(j)-d50tmp)));
       }
    }  
  }
 }
 
 else if (sel_option==3)
 {
    bicubic_spline(scal_yr_nodes,scal_age_nodes,sel_par,log_sel_fish);
 }

 else if (sel_option==4)
 {
      bincount = 1;
      for (i=sel_styr; i<=sel_endyr; i++)
          {
            if (i == binstart(bincount))
            {
             log_sel_fish(i)=cubic_spline(sel_par(bincount) );
             if( bincount < nbins ) bincount++;   
            }
            if (i < endyr_r) log_sel_fish(i+1) = log_sel_fish(i);
          }
  }


  
 for (i=styr; i<=endyr_r; i++)
  {     
   if (i < sel_styr) 
      { 
      sel_fish(i)(1,nselages) = mfexp(log_sel_fish(sel_styr));
      }
   else if (i > sel_endyr)  
    {
    sel_fish(i)(1,nselages) = mfexp(log_sel_fish(sel_endyr));
    }
   else   { sel_fish(i)(1,nselages) = mfexp(log_sel_fish(i));}
  
   if (nselages<nages)  
      sel_fish(i)(nselages+1,nages)=sel_fish(i,nselages);
  }
 

   
 

   

  

 //sel_fish=mfexp(log_sel_fish);    // convert from log values

 sel_srv3(1,nselages)=mfexp(log_sel_srv3);
 sel_srv_ebs(1,nselages)=mfexp(log_sel_srv_ebs);
 
 if (nselages<nages)
  {
   sel_srv3(nselages+1,nages) = sel_srv3(nselages);
   sel_srv_ebs(nselages+1,nages) = sel_srv_ebs(nselages);
  }

 }
 
 
  // rescale so max is one

   //sel_srv3 /= max(sel_srv3);

  //for (i=styr; i<=endyr;i++)
  //{
  //  sel_fish(i) /= max(sel_fish(i));
  //}
  
 

FUNCTION get_mortality //----------------------------------------------------------------------------------
 // Calulate the values of F and Z for each age and year
   M = exp(logM);
   
   avg_fmort_dev=mean(fmort_dev);
   for (i=styr_rec; i<=styr_fish-1; i++)  // set F to zero for backup years and years with no fishery
   {
   	for (j=1;j<=nages;j++)
 	{
 	  F(i,j) = 0.0;
 	}
   }
   
   for (i=styr_fish; i<=endyr_r; i++)  // use the selectivity values for the foreign fishery through 1988
   {
   	for (j=1;j<=nages;j++)
   	{
   	   F(i,j)=sel_fish(i,j)*mfexp(log_avg_fmort + fmort_dev(i));       
  	}
   }
   
   Z=F+M;
   S=mfexp(-1.0*Z);
   spawner_S = mfexp(-spmo_frac*Z);
   
FUNCTION get_first_year  //---------------------------------------------------------------------------------------------
          // get the first year of the natage matrix. Note that in both cases below, the 'styr_rec' 
          // refers to the first year for which there is an estimated value of recruitment.  For the 
          // first year age comp option 2, this is simply the styr of the model minus the age of recruitment.  
          // For option 1, this is the styr-(nages-1)-rec_age.  The first year of the spawning biomass vector is styr_rec-rec_age
  
 surv = mfexp(-1.0*M);
 if (fyear_ac_option == 1)  // First year stochastic recruitment devs are coupled with regular recruitment deviations
 {
  natage(styr_rec,1) = mfexp(mean_log_rec)*exp((sigr*sigr)/2);
  for (j=2; j<=nages; j++)
    natage(styr_rec,j) = natage(styr_rec,j-1)*surv;
  natage(styr_rec,nages) /= (1.-surv);

  for (j=styr_rec; j<=styr_rec_dev; j++)   // deviations in the cohorts that make up the plus group are shared
    natage(j,1) = mfexp(mean_log_rec+rec_dev(styr_rec_dev));

  for (j=styr_rec_dev+1; j<styr; j++)
    natage(j,1) = mfexp(mean_log_rec+rec_dev(j));  

   for (j=styr_rec; j<styr; j++)   // get sp_biom and natage for years prior to styr
    {
      natage(j+1)(2,nages) = ++elem_prod(natage(j)(1,nages-1),S(j)(1,nages-1));
      natage(j+1,nages) += natage(j,nages)*S(j,nages);
    }  
 }

 else if (fyear_ac_option == 2) {  // Initial age comps are in equilibrium with historic catch
            			
 natage(styr,1) = mfexp(log_rinit)*exp((sigr*sigr)/2);		// first, write the first age, first year as rzero
  					
 for (j=2; j<=nages;j++)			// next, get the first year ages (2,nages)
    natage(styr,j)=natage(styr,j-1)*mfexp(-(historic_F*sel_fish(styr)(j-1)+M));
 natage(styr,nages) /= (1-mfexp(-(historic_F*sel_fish(styr)(nages)+M)));  	// Plus group for first year

    if (historic_catch > 0.) {  // estimate the historical catch
      ehc = 0;
      for (j=1;j<=nages;j++)
          {
          ehc += natage(styr,j)*wt_fsh(j)*(historic_F*sel_fish(styr,j))*
 	      (1.0-mfexp(-(historic_F*sel_fish(styr,j)+M)))/(historic_F*sel_fish(styr,j)+M);
          }
        } 
 }

  else if (fyear_ac_option == 3) {  // Initial age comps are stochastic, but have a different mean than 
  				   //  the other recruitments. fydev is noise around equilibrium decay.					
  for (j=2; j<=nages_dat;j++)			
    natage(styr,j)=mfexp((log_rinit)  -M*double(j-1) + fydev(j));
  if (nages>nages_dat)
    {                                  // for the 'extra' ages needed to account for aging error, 
     for (j=nages_dat+1; j<=nages;j++) // use a single fydev for all cohorts in the plus group                                        
       natage(styr,j)=mfexp((log_rinit)  -M*double(j-1)+fydev(nages_dat));
    }
   else   
     natage(styr,nages) = mfexp((log_rinit) - M*double(nages-1) + fydev(nages_dat))/(1-mfexp(-M));  // plus group
 }


FUNCTION get_numbers_at_age //--------------------------------------------------------------------------------
   // Get numbers for the age of recruitment for all years, and fill out the natage matrix

 // get the recruits  
 for (i=styr;i<=lastyr_rec;i++)  natage(i,1) = mfexp(mean_log_rec+rec_dev(i));   
 for (i=lastyr_rec+1;i<=endyr_r+1;i++) natage(i,1) = mfexp(mean_log_rec +sigr*sigr*0.5);   // ages where we fix the recruits

 for (i=styr;i<=endyr_r;i++)   	// get natage matrix 
  {
   natage(i+1)(2,nages)=++elem_prod(natage(i)(1,nages-1),S(i)(1,nages-1));
   natage(i+1,nages)+=natage(i,nages)*S(i,nages);  // survival of plus group
  }

FUNCTION dvar_matrix get_projection_numbers_at_age() //------------------------------------------------------------------------------------------
  // Get the numbers at age for the projection model, which uses the mean recruitment for year classes which
  // have not exceeded the criteria for the survey and/or fishery selectivity
  RETURN_ARRAYS_INCREMENT();

  // the last year of recruitment to meet the survey a10 condition 
  double lastyr_rec_a10;
  lastyr_rec_a10 = endyr_r - max(fixedrec_yrs,excludeage - (rec_age - 1));

  dvar_matrix natage_mean_tmp(lastyr_rec_a10+1,endyr_r+1,1,nages);    // numbers at age

  // get the first year
  natage_mean_tmp(lastyr_rec_a10+1) = natage(lastyr_rec_a10+1);
  
  // get the recruits  
  for (i=lastyr_rec_a10+1;i<=endyr_r+1;i++)  natage_mean_tmp(i,1) = mfexp(mean_log_rec +sigr*sigr*0.5);   
  
  for (i=lastyr_rec_a10+1;i<=endyr_r;i++)    // get natage matrix 
    {
   natage_mean_tmp(i+1)(2,nages)=++elem_prod(natage_mean_tmp(i)(1,nages-1),S(i)(1,nages-1));
   natage_mean_tmp(i+1,nages)+=natage_mean_tmp(i,nages)*S(i,nages);  // survival of plus group
    }

 RETURN_ARRAYS_DECREMENT();
 return natage_mean_tmp;

FUNCTION get_expected_values  //-----------------------------------------------------------------------------
   // get reproductive outputs, total biomass, and survey biomass   
   sp_biom.initialize();
   eggs.initialize();
   larvae.initialize();
   sp_biom(styr_rec-rec_age,styr_rec-1) = elem_prod(wt_pop,maturity)*elem_prod(natage(styr_rec)/2.,spawner_S(styr_rec));
 //  eggs(styr_rec-rec_age,styr_rec-1) = elem_prod(maturity,fecund)*elem_prod(natage(styr_rec)/2.,spawner_S(styr_rec));   
 //  larvae(styr_rec-rec_age,styr_rec-1) = elem_prod(elem_prod(maturity,fecund),lvsrv)*elem_prod(natage(styr_rec)/2.,spawner_S(styr_rec));
   totbiom(styr_rec-rec_age,styr_rec-1) = natage(styr_rec)*wt_pop;

    for (j=styr_rec; j<=endyr_r; j++)   // get sp_biom and natagetmp for years prior to styr
    {
      sp_biom(j) = elem_prod(wt_pop,maturity)*elem_prod(natage(j)/2.,spawner_S(j));
  //    eggs(j) = elem_prod(maturity,fecund)*elem_prod(natage(j)/2.,spawner_S(j));   
  //    larvae(j) = elem_prod(elem_prod(maturity,fecund),lvsrv)*elem_prod(natage(j)/2.,spawner_S(j));
      totbiom(j) = natage(j)*wt_pop;
    }
   totbiom(endyr_r+1) = natage(endyr_r+1)*wt_pop;

  // compute the predicted values for the surveys
 //q_cpue = exp(logq_cpue);
 q_srv3 = exp(logq_srv3);
 q_srv_ebs = exp(logq_srv_ebs);

 for (i=styr;i<=endyr_r;i++)
  mort(i) = elem_div((1.-mfexp(-Z(i))),Z(i));
 
 //for (i=1;i<=nyrs_cpue;i++){   // survey 1 is the cpue index
 //   pred_cpue(i)=q_cpue*elem_prod(natage(yrs_cpue(i)),mort(yrs_cpue(i)))*
 //      elem_prod(sel_fish(yrs_cpue(i)),wt_pop);}
 
 for (i=styr_fish;i<=endyr_r;i++){   // survey 3 is the AI  trawl survey 
   pred_srv3(i)=prop_bio_ai(i)*q_srv3*elem_prod(natage(i),mort(i))*
      elem_prod(sel_srv3,wt_pop);
   pred_srv_ebs(i) = (1.0-prop_bio_ai(i))*q_srv_ebs*elem_prod(natage(i),mort(i))*
      elem_prod(sel_srv_ebs,wt_pop);
    }

FUNCTION get_sr_inputs //--------------------------------------------------------------------------------------------------
  // get the inputs for the SR curve 
 // first, define rzero and set bzero, if needed
 if (active(log_rzero))
 { 
 rzero = mfexp(log_rzero);
 natagetmp = 0.0;
 natagetmp(1) = rzero;
 for (j=2; j<=nages; j++)
  natagetmp(j) = natagetmp(j-1)*surv;
 natagetmp(nages) /= (1.-surv);
 bzero = elem_prod(wt_pop,maturity)*natagetmp*0.5;
 alpha = 0.8*rzero*steepness/(steepness-0.2);
 beta = 0.2*bzero*((1.-steepness)/(steepness-0.2));  
 }

 est_rec(styr_rec_dev,lastyr_rec) = column(natage,1)(styr_rec_dev,lastyr_rec);  	// get the estimated recruitment
 dvar_vector Stmp(styr_rec_dev, lastyr_rec); 		// temporary S
 Stmp = sp_biom(styr_rec_dev-rec_age, lastyr_rec-rec_age).shift(styr_rec_dev);  //assign the ssb
 est_spb = Stmp;				// save the spawning biomass
 pred_rec = SRecruit(Stmp);			// get the predicted recruits
 
 // get the data for the fitted recruitment curve
 dvariable tmpsp=1.1*max(est_spb);
 for (i=1;i<=20;i++)
 {
   SRec_spawn(i)=tmpsp*double(i)/20.;
   SRec_rec(i)=SRecruit(SRec_spawn(i));
 }

FUNCTION get_catch_at_age // by using Baranov's Catch Equation//-----------------------------------------------
 for (i=styr_fish; i<=endyr_r; i++)
 {
    pred_catch(i) = 0.0;
    for (j=1;j<= nages;j++)
    {
 	catage(i,j) = natage(i,j)*F(i,j)*(1.0-S(i,j))/Z(i,j);
 	pred_catch(i)+=catage(i,j)*wt_fsh(j);
    }
  }
 
FUNCTION get_age_comps  // need to apply age error matrices  
  for (i=1;i<=nyrs_fish_unbiased_ac_r;i++)
   {
    eac_fish_unbiased_mod(yrs_fish_unbiased_ac_r(i))=catage(yrs_fish_unbiased_ac_r(i))/sum(catage(yrs_fish_unbiased_ac_r(i)));
    eac_fish_unbiased_mod(yrs_fish_unbiased_ac_r(i)) = unbiasedages*eac_fish_unbiased_mod(yrs_fish_unbiased_ac_r(i));
    eac_fish_unbiased_dat(yrs_fish_unbiased_ac_r(i))(1,nages_dat-1) = eac_fish_unbiased_mod(yrs_fish_unbiased_ac_r(i))(1,nages_dat-1);
    eac_fish_unbiased_dat(yrs_fish_unbiased_ac_r(i))(nages_dat) = sum(eac_fish_unbiased_mod(yrs_fish_unbiased_ac_r(i))(nages_dat,nages));
   }
 
//  for (i=1;i<=nyrs_fish_biased_ac;i++)
//   {
//    eac_fish_biased(yrs_fish_biased_ac(i))=catage(yrs_fish_biased_ac(i))/sum(catage(yrs_fish_biased_ac(i)));
//    eac_fish_biased(yrs_fish_biased_ac(i)) = biasedages*eac_fish_biased(yrs_fish_biased_ac(i));
//   }
 
  for (i=1;i<=nyrs_fish_lc_r;i++)
   {
    elc_fish_tempages(yrs_fish_lc_r(i))=catage(yrs_fish_lc_r(i))/sum(catage(yrs_fish_lc_r(i)));
    elc_fish(yrs_fish_lc_r(i))=translen*elc_fish_tempages(yrs_fish_lc_r(i));
   }
 
  for (i=1;i<=nyrs_surv3_unbiased_ac_r;i++)
   {    
    eac_surv3_unbiased_mod(yrs_surv3_unbiased_ac_r(i))=elem_prod(sel_srv3,elem_prod(  natage(yrs_surv3_unbiased_ac_r(i)),  mort(yrs_surv3_unbiased_ac_r(i))  ))/
       (sel_srv3*elem_prod(  natage(yrs_surv3_unbiased_ac_r(i)),  mort(yrs_surv3_unbiased_ac_r(i))));
    eac_surv3_unbiased_mod(yrs_surv3_unbiased_ac_r(i)) = unbiasedages*eac_surv3_unbiased_mod(yrs_surv3_unbiased_ac_r(i));
    eac_surv3_unbiased_dat(yrs_surv3_unbiased_ac_r(i))(1,nages_dat-1) = eac_surv3_unbiased_mod(yrs_surv3_unbiased_ac_r(i))(1,nages_dat-1);
    eac_surv3_unbiased_dat(yrs_surv3_unbiased_ac_r(i))(nages_dat) = sum(eac_surv3_unbiased_mod(yrs_surv3_unbiased_ac_r(i))(nages_dat,nages));
    }

  for (i=1;i<=nyrs_surv3_lc_r;i++)
   {
    elc_surv3_tempages(yrs_surv3_lc_r(i))=elem_prod(sel_srv3,natage(yrs_surv3_lc_r(i)))/
       (sel_srv3*natage(yrs_surv3_lc_r(i)));
    elc_surv3(yrs_surv3_lc_r(i))=translen*elc_surv3_tempages(yrs_surv3_lc_r(i));
   }

   for (i=1;i<=nyrs_surv_ebs_unbiased_ac_r;i++)
   {    
    eac_surv_ebs_unbiased_mod(yrs_surv_ebs_unbiased_ac_r(i))=elem_prod(sel_srv_ebs,elem_prod(  natage(yrs_surv_ebs_unbiased_ac_r(i)),  mort(yrs_surv_ebs_unbiased_ac_r(i))  ))/
       (sel_srv_ebs*elem_prod(  natage(yrs_surv_ebs_unbiased_ac_r(i)),  mort(yrs_surv_ebs_unbiased_ac_r(i))));
    eac_surv_ebs_unbiased_mod(yrs_surv_ebs_unbiased_ac_r(i)) = unbiasedages*eac_surv_ebs_unbiased_mod(yrs_surv_ebs_unbiased_ac_r(i));
    eac_surv_ebs_unbiased_dat(yrs_surv_ebs_unbiased_ac_r(i))(1,nages_dat-1) = eac_surv_ebs_unbiased_mod(yrs_surv_ebs_unbiased_ac_r(i))(1,nages_dat-1);
    eac_surv_ebs_unbiased_dat(yrs_surv_ebs_unbiased_ac_r(i))(nages_dat) = sum(eac_surv_ebs_unbiased_mod(yrs_surv_ebs_unbiased_ac_r(i))(nages_dat,nages));
    }

    for (i=1;i<=nyrs_surv_ebs_lc_r;i++)
   {
    elc_surv_ebs_tempages(yrs_surv_ebs_lc_r(i))=elem_prod(sel_srv_ebs,natage(yrs_surv_ebs_lc_r(i)))/
       (sel_srv_ebs*natage(yrs_surv_ebs_lc_r(i)));
    elc_surv_ebs(yrs_surv_ebs_lc_r(i))=translen*elc_surv_ebs_tempages(yrs_surv_ebs_lc_r(i));
   }



  
 //if (sd_phase())
 //{
 //  depletion = totbiom(endyr)/totbiom(styr);
 //  endbiom=totbiom(endyr);
 //}

FUNCTION get_binned
  // bin the natage matrix to match the plus group for the data

  for (i=styr_rec;i<=endyr_r+1;i++)
    {
       natage_bin(i)(1,nages_dat-1) = natage(i)(1,nages_dat-1);
       natage_bin(i)(nages_dat) = sum(natage(i)(nages_dat,nages)); 
    }
  maturity_bin(1,nages_dat-1) = maturity(1,nages_dat-1);
  maturity_bin(nages_dat) = mean(maturity(nages_dat,nages));

  
FUNCTION dvariable SRecruit(const dvariable& Stmp)
 RETURN_ARRAYS_INCREMENT();
 dvariable RecTmp;
 switch (sr_type)
 {
    case 1:
      RecTmp = mfexp(mean_log_rec +sigr*sigr*0.5);  // average recruitment
      break;
    case 2:
      RecTmp = alpha*Stmp*(1./(beta + Stmp)); // Beverton-Holt form
      break;
 }
 RETURN_ARRAYS_DECREMENT();
 return RecTmp;

FUNCTION dvar_vector SRecruit(const dvar_vector& Stmp)
 RETURN_ARRAYS_INCREMENT();
 dvar_vector RecTmp(Stmp.indexmin(),Stmp.indexmax());
 switch (sr_type)
 {
    case 1:
      RecTmp = mfexp(mean_log_rec +sigr*sigr*0.5);  // average recruitment
      break;
    case 2:
      RecTmp = elem_prod(alpha*Stmp,1./(beta + Stmp)); // Beverton-Holt form
      break;
 }
 RETURN_ARRAYS_DECREMENT();
 return RecTmp;

FUNCTION evaluate_the_objective_function //--------------------------------------------------------------------
 if (fyear_ac_option == 2 && historic_catch > 0.)  // computes the ssq for the historical F 
      histFpen();
 mat_likelihood();
 //cout <<" got mat_like"  << endl;
 rec_likelihood();
 //cout <<" ggt rec_like"  << endl;
 srv_likelihood();
 //cout <<" got srv_like"  << endl;
 cat_likelihood();
 //cout <<" got cat_like"  << endl;
 Fmort_pen();
 //cout <<" get Fmort_pen"  << endl;
 age_likelihood();
 //cout <<" get age_like"  << endl;
 prior();
 //cout <<" get prior"  << endl;
 sel_likelihood();
 //cout <<" sel_like"  << endl;

  
FUNCTION histFpen  // fit the historical catches if neccessary
 if (active(historic_F))
  {
  hf_pen = 500.*square(ehc - historic_catch);
  obj_fun += hf_pen;
  }
  

FUNCTION rec_likelihood   // fit the recruitment deviations
 rec_like.initialize();
 chi = log(est_rec+1e-8) - log(pred_rec+1e-8) + (sigr*sigr*0.5);  // correct for bias;
 //cout << chi << endl;

 rmse(3) = sqrt(norm2(log(est_rec+1e-8) - log(pred_rec+1e-8))/size_count(chi)+1e-13);
 sumrecdev = sum(chi); 
 
 rec_like(1) = norm2(chi)/(2.*sigr*sigr) + size_count(chi)*log(sigr);
 //rec_like(2) = square((rmse - sigr)/(sigr*sigr/size_count(chi)))/2. + log(sigr*sigr/size_count(chi)) +
 //   square(sumrecdev/(sigr/size_count(chi)))/2. + log(sigr/size_count(chi));
 
 if (fyear_ac_option == 3)
  rec_like(2) = norm2(fydev)/(2.*sigr*sigr) + size_count(fydev)*log(sigr);

 obj_fun +=lambda(1)*sum(rec_like);
 
 
 //  rec_like+=1.0*norm2(rec_dev_future); // the deviations for the future recruitments 
  
FUNCTION age_likelihood  // fit the age comps--------------------------------------------------------------------
 dvariable tmp1;
 dvariable tmp2;
 dvariable tmp3;
 dvariable tmp4;
 dvariable tmp5;
 dvariable tmp6;
 dvariable tmp7;

 age_like=0.;
 effn=0.;
 int ii;
 // offset = 1;
 // for (i=1;i<=nyrs_fish_biased_ac;i++)
 // {
 //   ii=yrs_fish_biased_ac(i);
 //   age_like(1)-=fish_biased_ac_samp(i)*oac_fish_biased(i)*log(eac_fish_biased(ii)+1.e-13);
 //   effn(1,i) = eac_fish_biased(ii)*(1.-eac_fish_biased(ii))/(norm2(eac_fish_biased(ii)-oac_fish_biased(i)));
 // }
 //age_like(1)-=offset(1);
 // age_like(1)=0.0; // do not fit the biased age comps   

 for (i=1;i<=nyrs_fish_unbiased_ac_r;i++)
 {
   ii=yrs_fish_unbiased_ac_r(i);
   age_like(1)-=fish_unbiased_ac_samp_r(i)*oac_fish_unbiased_r(i)*log(eac_fish_unbiased_dat(ii)+1.e-13);
   effn(1,i) = eac_fish_unbiased_dat(ii)*(1.-eac_fish_unbiased_dat(ii))/(norm2(eac_fish_unbiased_dat(ii)-oac_fish_unbiased_r(i))); 
 }
 age_like(1)-=offset(1);
  
 for (i=1;i<=nyrs_fish_lc_r;i++)
 {
   ii=yrs_fish_lc_r(i);
   age_like(2)-=fish_lc_samp_r(i)*olc_fish_r(i)*log(elc_fish(ii)+1.e-13);
   effn(2,i) = elc_fish(ii)*(1.-elc_fish(ii))/(norm2(elc_fish(ii)-olc_fish_r(i)));
 }
 age_like(2)-=offset(2);
 
 for (i=1;i<=nyrs_surv3_unbiased_ac_r;i++)
 {
   ii=yrs_surv3_unbiased_ac_r(i);
   age_like(3)-=surv3_unbiased_ac_samp_r(i)*oac_surv3_unbiased_r(i)*log(eac_surv3_unbiased_dat(ii)+1.e-13);
   effn(3,i) = eac_surv3_unbiased_dat(ii)*(1.-eac_surv3_unbiased_dat(ii))/(norm2(eac_surv3_unbiased_dat(ii)-oac_surv3_unbiased_r(i)));
 }
 age_like(3)-=offset(3);

 for (i=1;i<=nyrs_surv3_lc_r;i++)
 {
   ii=yrs_surv3_lc_r(i);
   age_like(4)-=surv3_lc_samp_r(i)*olc_surv3_r(i)*log(elc_surv3(ii)+1.e-13);
   effn(4,i) = elc_surv3(ii)*(1.-elc_surv3(ii))/(norm2(elc_surv3(ii)-olc_surv3_r(i)));
 }
 age_like(4)-=offset(4);

 for (i=1;i<=nyrs_surv_ebs_unbiased_ac_r;i++)
 {
   ii=yrs_surv_ebs_unbiased_ac_r(i);
   age_like(5)-=surv_ebs_unbiased_ac_samp_r(i)*oac_surv_ebs_unbiased_r(i)*log(eac_surv_ebs_unbiased_dat(ii)+1.e-13);
   effn(5,i) = eac_surv_ebs_unbiased_dat(ii)*(1.-eac_surv_ebs_unbiased_dat(ii))/(norm2(eac_surv_ebs_unbiased_dat(ii)-oac_surv_ebs_unbiased_r(i)));
 }
 age_like(5)-=offset(5);

 for (i=1;i<=nyrs_surv_ebs_lc_r;i++)
 {
   ii=yrs_surv_ebs_lc_r(i);
   age_like(6)-=surv_ebs_lc_samp_r(i)*olc_surv_ebs_r(i)*log(elc_surv_ebs(ii)+1.e-13);
   effn(6,i) = elc_surv_ebs(ii)*(1.-elc_surv_ebs(ii))/(norm2(elc_surv_ebs(ii)-olc_surv_ebs_r(i)));
 }
 age_like(6)-=offset(6); 



 //obj_fun = obj_fun + (age_like(1)*lambda(12) + age_like(2)*lambda(5) +age_like(3)*lambda(6) +age_like(4)*lambda(7));
 obj_fun = obj_fun + (age_like(1)*lambda(5) +age_like(2)*lambda(6) +age_like(3)*lambda(7) +age_like(4)*lambda(8) +age_like(5)*lambda(9) +age_like(6)*lambda(10));
 
 // now get the SDNR for the age and length comps
 // k = 0;
 // tmp1 = 0.0;
 // tmp2 = 0.0;
 // for (i=1;i<=nyrs_fish_biased_ac;i++)
 // {
 //   ii=yrs_fish_biased_ac(i);
 //   for (j=1;j<=nages;j++)
 //   {
 //     k = k+1;
 //     tmp1 = (oac_fish_biased(i,j)+0.00001) - (eac_fish_biased(ii,j) + 0.00001);
 //     tmp2 = (eac_fish_biased(ii,j)+0.00001)*(1.-(eac_fish_biased(ii,j)+0.00001)  );
 //     fbac_nr(k) = tmp1/sqrt(tmp2/fish_biased_ac_samp(i));
 //   } 
 // }

 k = 0;
 tmp1 = 0.0;
 tmp2 = 0.0;
 tmp3 = 0.0;
 tmp4 = 0.0;
 tmp5 = 0.0;
 tmp6 = 0.0;
 tmp7 = 0.0;
 for (i=1;i<=nyrs_fish_unbiased_ac_r;i++)
 {
   ii=yrs_fish_unbiased_ac_r(i);
   tmp3 = ages_dat*oac_fish_unbiased_r(i);   // the mean of the observations 
   tmp4 = ages_dat*eac_fish_unbiased_dat(ii);  // the mean of the predictions
   tmp6 = (eac_fish_unbiased_dat(ii)+0.00001)*(1.-(eac_fish_unbiased_dat(ii)+0.00001));
   tmp7 = ((oac_fish_unbiased_r(i)+0.00001) - (eac_fish_unbiased_dat(ii) + 0.00001))*
             ((oac_fish_unbiased_r(i)+0.00001) - (eac_fish_unbiased_dat(ii) + 0.00001));
   fac_mcian_wgt(i) = (tmp6/tmp7)/raw_fish_unbiased_ac_samp(i);
   fac_mcian_wgt_inv(i) = 1.0/fac_mcian_wgt(i);
   tmp5 = elem_prod(ages_dat,eac_fish_unbiased_dat(ii))*ages_dat - square(tmp4);  // the v term in Francis method
   fac_nr_fran(i) = (tmp3-tmp4)/sqrt(tmp5/raw_fish_unbiased_ac_samp(i));           


   for (j=1;j<=nages_dat;j++)
   {
     k = k+1;
     tmp1 = (oac_fish_unbiased_r(i,j)+0.00001) - (eac_fish_unbiased_dat(ii,j) + 0.00001);
     tmp2 = (eac_fish_unbiased_dat(ii,j)+0.00001)*(1.-(eac_fish_unbiased_dat(ii,j)+0.00001)  );
     fac_nr(k) = tmp1/sqrt(tmp2/raw_fish_unbiased_ac_samp(i));
     fac_resid(k) = tmp1; 
   } 
 }

 k = 0;
 tmp1 = 0.0;
 tmp2 = 0.0;
 tmp3 = 0.0;
 tmp4 = 0.0;
 tmp5 = 0.0;
 tmp6 = 0.0;
 tmp7 = 0.0;
 
 for (i=1;i<=nyrs_fish_lc_r;i++)
 {
   ii=yrs_fish_lc_r(i);
   tmp3 = lengths*olc_fish_r(i);   // the mean of the observations 
   tmp4 = lengths*elc_fish(ii);  // the mean of the predictions
   tmp6 = (elc_fish(ii)+0.00001)*(1.-(elc_fish(ii)+0.00001));
   tmp7 = ((olc_fish_r(i)+0.00001) - (elc_fish(ii) + 0.00001))*
             ((olc_fish_r(i)+0.00001) - (elc_fish(ii) + 0.00001));
   flc_mcian_wgt(i) = (tmp6/tmp7)/raw_fish_lc_samp(i);
   flc_mcian_wgt_inv(i) =  1.0/flc_mcian_wgt(i);
   tmp5 = elem_prod(lengths,elc_fish(ii))*lengths - square(tmp4);  // the v term in Francis method
   flc_nr_fran(i) = (tmp3-tmp4)/sqrt(tmp5/raw_fish_lc_samp(i));
  
   
   for (j=1;j<=nlen;j++)
   {
     k = k+1;
     tmp1 = (olc_fish_r(i,j)+0.00001) - (elc_fish(ii,j) + 0.00001);
     tmp2 = (elc_fish(ii,j)+0.00001)*(1.-(elc_fish(ii,j)+0.00001)  );
     flc_nr(k) = tmp1/sqrt(tmp2/raw_fish_lc_samp(i));
     flc_resid(k) = tmp1;

   } 
 }

 k = 0;
 tmp1 = 0.0;
 tmp2 = 0.0;
 tmp3 = 0.0;
 tmp4 = 0.0;
 tmp5 = 0.0;
 tmp6 = 0.0;
 tmp7 = 0.0;

 for (i=1;i<=nyrs_surv3_unbiased_ac_r;i++)
 {
   ii=yrs_surv3_unbiased_ac_r(i);
   tmp3 = ages_dat*oac_surv3_unbiased_r(i);  // the mean of the observations 
   tmp4 = ages_dat*eac_surv3_unbiased_dat(ii); // the mean of the predictions 
   tmp6 = (eac_surv3_unbiased_dat(ii)+0.00001)*(1.-(eac_surv3_unbiased_dat(ii)+0.00001));
   tmp7 = ((oac_surv3_unbiased_r(i)+0.00001) - (eac_surv3_unbiased_dat(ii) + 0.00001))*
             ((oac_surv3_unbiased_r(i)+0.00001) - (eac_surv3_unbiased_dat(ii) + 0.00001)); 
   sac_mcian_wgt(i) = (tmp6/tmp7)/raw_surv3_unbiased_ac_samp(i);
   sac_mcian_wgt_inv(i) = 1.0/sac_mcian_wgt(i);
   tmp5 = elem_prod(ages_dat,eac_surv3_unbiased_dat(ii))*ages_dat - square(tmp4);  // the v term in Francis method
   sac_nr_fran(i) = (tmp3-tmp4)/sqrt(tmp5/raw_surv3_unbiased_ac_samp(i));
 
   
   for (j=1;j<=nages_dat;j++)
   {
     k = k+1;
     tmp1 = (oac_surv3_unbiased_r(i,j)+0.00001) - (eac_surv3_unbiased_dat(ii,j) + 0.00001);
     tmp2 = (eac_surv3_unbiased_dat(ii,j)+0.00001)*(1.-(eac_surv3_unbiased_dat(ii,j)+0.00001)  );
     sac_nr(k) = tmp1/sqrt(tmp2/raw_surv3_unbiased_ac_samp(i));
     sac_resid(k) = tmp1; 
   } 
 }

 k = 0;
 tmp1 = 0.0;
 tmp2 = 0.0;
 tmp3 = 0.0;
 tmp4 = 0.0;
 tmp5 = 0.0;
 tmp6 = 0.0;
 tmp7 = 0.0;
 for (i=1;i<=nyrs_surv3_lc_r;i++)
 {
   ii=yrs_surv3_lc_r(i);
   tmp3 = lengths*olc_surv3_r(i);
   tmp4 = lengths*elc_surv3(ii);
   tmp6 = (elc_surv3(ii)+0.00001)*(1.-(elc_surv3(ii)+0.00001));
   tmp7 = ((olc_surv3_r(i)+0.00001) - (elc_surv3(ii) + 0.00001))*
             ((olc_surv3_r(i)+0.00001) - (elc_surv3(ii) + 0.00001)); 
   slc_mcian_wgt(i) = (tmp6/tmp7)/raw_surv3_lc_samp(i);
   slc_mcian_wgt_inv(i) = 1.0/slc_mcian_wgt(i);
   //tmp5 = sqrt((elem_prod(lengths-tmp3,lengths-tmp3)*(olc_surv3_r(i)))/surv3_lc_samp_r(i));
   tmp5 = elem_prod(lengths,elc_surv3(ii))*lengths - square(tmp4);  // the v term in Francis method
   slc_nr_fran(i) = (tmp3-tmp4)/sqrt(tmp5/raw_surv3_lc_samp(i)); 
   

   for (j=1;j<=nlen;j++)
   {
     k = k+1;
     tmp1 = (olc_surv3_r(i,j)+0.00001) - (elc_surv3(ii,j) + 0.00001);
     tmp2 = (elc_surv3(ii,j)+0.00001)*(1.-(elc_surv3(ii,j)+0.00001)  );
     slc_nr(k) = tmp1/sqrt(tmp2/raw_surv3_lc_samp(i)); 
     slc_resid(k) = tmp1;
   } 
 }

 k = 0;
 tmp1 = 0.0;
 tmp2 = 0.0;
 tmp3 = 0.0;
 tmp4 = 0.0;
 tmp5 = 0.0;
 tmp6 = 0.0;
 tmp7 = 0.0;
 for (i=1;i<=nyrs_surv_ebs_unbiased_ac_r;i++)
 {
   ii=yrs_surv_ebs_unbiased_ac_r(i);
   tmp3 = ages_dat*oac_surv_ebs_unbiased_r(i);  // the mean of the observations 
   tmp4 = ages_dat*eac_surv_ebs_unbiased_dat(ii); // the mean of the predictions 
   tmp6 = (eac_surv_ebs_unbiased_dat(ii)+0.00001)*(1.-(eac_surv_ebs_unbiased_dat(ii)+0.00001));
   tmp7 = ((oac_surv_ebs_unbiased_r(i)+0.00001) - (eac_surv_ebs_unbiased_dat(ii) + 0.00001))*
             ((oac_surv_ebs_unbiased_r(i)+0.00001) - (eac_surv_ebs_unbiased_dat(ii) + 0.00001)); 
   ebs_sac_mcian_wgt(i) = (tmp6/tmp7)/raw_surv_ebs_unbiased_ac_samp(i);
   ebs_sac_mcian_wgt_inv(i) = 1.0/ebs_sac_mcian_wgt(i);
   tmp5 = elem_prod(ages_dat,eac_surv_ebs_unbiased_dat(ii))*ages_dat - square(tmp4);  // the v term in Francis method
   ebs_sac_nr_fran(i) = (tmp3-tmp4)/sqrt(tmp5/raw_surv_ebs_unbiased_ac_samp(i));
 
   for (j=1;j<=nages_dat;j++)
   {
     k = k+1;
     tmp1 = (oac_surv_ebs_unbiased_r(i,j)+0.00001) - (eac_surv_ebs_unbiased_dat(ii,j) + 0.00001);
     tmp2 = (eac_surv_ebs_unbiased_dat(ii,j)+0.00001)*(1.-(eac_surv_ebs_unbiased_dat(ii,j)+0.00001)  );
     ebs_sac_nr(k) = tmp1/sqrt(tmp2/raw_surv_ebs_unbiased_ac_samp(i));
     ebs_sac_resid(k) = tmp1; 
   } 
 }

 k = 0;
 tmp1 = 0.0;
 tmp2 = 0.0;
 tmp3 = 0.0;
 tmp4 = 0.0;
 tmp5 = 0.0;
 tmp6 = 0.0;
 tmp7 = 0.0;
 for (i=1;i<=nyrs_surv_ebs_lc_r;i++)
 {
   ii=yrs_surv_ebs_lc_r(i);
   tmp3 = lengths*olc_surv_ebs_r(i);
   tmp4 = lengths*elc_surv_ebs(ii);
   tmp6 = (elc_surv_ebs(ii)+0.00001)*(1.-(elc_surv_ebs(ii)+0.00001));
   tmp7 = ((olc_surv_ebs_r(i)+0.00001) - (elc_surv_ebs(ii) + 0.00001))*
             ((olc_surv_ebs_r(i)+0.00001) - (elc_surv_ebs(ii) + 0.00001)); 
   ebs_slc_mcian_wgt(i) = (tmp6/tmp7)/raw_surv_ebs_lc_samp(i);
   ebs_slc_mcian_wgt_inv(i) = 1.0/ebs_slc_mcian_wgt(i); 
   //tmp5 = sqrt((elem_prod(lengths-tmp3,lengths-tmp3)*(olc_surv_ebs_r(i)))/surv_ebs_lc_samp_r(i));
   tmp5 = elem_prod(lengths,elc_surv_ebs(ii))*lengths - square(tmp4);  // the v term in Francis method
   ebs_slc_nr_fran(i) = (tmp3-tmp4)/sqrt(tmp5/raw_surv_ebs_lc_samp(i)); 
   

   for (j=1;j<=nlen;j++)
   {
     k = k+1;
     tmp1 = (olc_surv_ebs_r(i,j)+0.00001) - (elc_surv_ebs(ii,j) + 0.00001);
     tmp2 = (elc_surv_ebs(ii,j)+0.00001)*(1.-(elc_surv_ebs(ii,j)+0.00001)  );
     ebs_slc_nr(k) = tmp1/sqrt(tmp2/raw_surv_ebs_lc_samp(i)); 
     ebs_slc_resid(k) = tmp1;
   } 
 }



 // sdnr(1) = std_dev(fbac_nr);
 //sdnr(1) = std_dev(fac_nr);
 //sdnr(2) = std_dev(flc_nr);
 //sdnr(3) = std_dev(sac_nr);
 //sdnr(4) = std_dev(slc_nr);
 //sdnr(5) = std_dev(ebs_sac_nr);
 //sdnr(6) = std_dev(ebs_slc_nr);

  al_rmse = 0;
  sdnr(1,6) = 0;

 if (nyrs_fish_unbiased_ac_r>0)
 {
    al_rmse(1) = sqrt(mean(elem_prod(fac_resid,fac_resid)));
    sdnr(1) = std_dev(fac_nr);
 }
 

 if (nyrs_fish_lc_r>0)
 {
   al_rmse(2) = sqrt(mean(elem_prod(flc_resid,flc_resid)));
   sdnr(2) = std_dev(flc_nr);
 }
 
 
 if (nyrs_surv3_unbiased_ac_r>0)
 {
   al_rmse(3) = sqrt(mean(elem_prod(sac_resid,sac_resid)));
   sdnr(3) = std_dev(sac_nr);
 }
 
 
 if (nyrs_surv3_lc_r>0)
 {
   al_rmse(4) = sqrt(mean(elem_prod(slc_resid,slc_resid)));
   sdnr(4) = std_dev(slc_nr);
 }
 

 if (nyrs_surv_ebs_unbiased_ac_r>0)
 { 
   al_rmse(5) = sqrt(mean(elem_prod(ebs_sac_resid,ebs_sac_resid)));
   sdnr(5) = std_dev(ebs_sac_nr);
 }
 
 
 if (nyrs_surv_ebs_lc_r>0)
 {
   al_rmse(6) = sqrt(mean(elem_prod(ebs_slc_resid,ebs_slc_resid)));
   sdnr(6) = std_dev(ebs_slc_nr);
 }
 

FUNCTION prior  // compute the prior parts for q_surv3 and M--------------------------------------------------------------------
 priormq(1) = square(logM - log(priormean_M) + square(priorcv_M)/2.0)/(2.*square(priorcv_M)); 
 priormq(2) = square(log(q_srv3) - log(priormeansrv3q) + square(priorcvsrv3q)/2.0)/(2.*square(priorcvsrv3q)); 
 
 obj_fun += sum(priormq);
 


FUNCTION dvariable get_spr(dvariable Ftemp)    // calculation of equilibrium SPR for equilibrium recruitment
  dvariable phi;
  dvar_vector Ntmp(1,nages_dat);
  

  Ntmp(1)=1.;
  for (j=2;j<=nages_dat;j++)
     Ntmp(j)=Ntmp(j-1)*exp(-(M+Ftemp*recent_fish_sel(j-1)));  // fills in matrix for ages 2 through nages_dat-1           
   Ntmp(nages_dat)=Ntmp(nages_dat-1)*exp(-(M+Ftemp*recent_fish_sel(nages_dat-1)))/(1.-exp(-(M+Ftemp*recent_fish_sel(nages_dat))));
  
   // Kill them off until they spawn
  for (j=1;j<=nages_dat;j++) 
   phi = 0.5*elem_prod(Ntmp,maturity_bin)*elem_prod(wt_pop_bin,exp(-spmo_frac*(M+Ftemp*recent_fish_sel)));

   return(phi);

FUNCTION dvariable get_spr_rates(double spr_percent)
  dvariable df=1.e-8;
  dvariable F1;
  dvariable dd;
  F1 = 0.0;
  if (M<0.2)  
    F1 = 1.2*M*(1-spr_percent);
  else
    F1 = 0.5*M*(1-spr_percent);
  dvariable F2;
  dvariable F3;
  dvariable yld1;
  dvariable yld2;
  dvariable yld3;
  dvariable dyld;
  dvariable dyldp;
  // Newton Raphson stuff to go here
  //for (int ii=1;ii<=12;ii++)
  dd = 1.0;
  while (dd > 1.e-10)   
 {
    F2     = F1 + df;
    F3     = F1 - df;
    yld1   = -1000*square(log(spr_percent/(get_spr(F1)/SB0)));
    yld2   = -1000*square(log(spr_percent/(get_spr(F2)/SB0)));
    yld3   = -1000*square(log(spr_percent/(get_spr(F3)/SB0)));
    dyld   = (yld2 - yld3)/(2*df);                          // First derivative (to find the root of this)
    //dyldp  = (yld3-(2*yld1)+yld2)/(df*df);  // Newton-Raphson approximation for second derivitive
    //F1    -= dyld/dyldp;
    F1 -= yld1/dyld;
    dd = fabs(log(spr_percent/(get_spr(F1)/SB0)));  
  }
  return(F1); 


FUNCTION srv_likelihood  // fit to indices (lognormal) ---------------------------------------------------------
 surv_like=0;
 rmse(1,2)=0;
  int ii;
 //for (i=1;i<=nyrs_cpue;i++)
 //{
 //  surv_like(1) += (square(log(obs_cpue(i)+1e-13) - log(pred_cpue(i)+1e-13))/(2.*cv_cpue*cv_cpue) + log(cv_cpue)); 
 //  rmse(1) += square(log(obs_cpue(i)+1e-13) - log(pred_cpue(i)+1e-13) ); 
 //
 //  cpue_nr(i) = (log(obs_cpue(i)+1e-13) - log(pred_cpue(i)+1e-13) )/cv_cpue;
 //}
 //rmse(1) = sqrt(rmse(1)/nyrs_cpue);

 for (i=1;i<=nyrs_srv3_r;i++)  // AI survey likelihood
 {
   ii=yrs_srv3(i);
   surv_like(1) += square(log(obs_srv3(i)+1e-13) - log(pred_srv3(ii)+1e-13) )/(2.*cv_srv3(i)*cv_srv3(i));
   rmse(1) += square(log(obs_srv3(i)+1e-13) - log(pred_srv3(ii)+1e-13));

   survey_nr(i) = (log(obs_srv3(i)+1e-13) - log(pred_srv3(ii)+1e-13) )/cv_srv3(i);
 }
  rmse(1) = sqrt(rmse(1)/nyrs_srv3_r);
  sdnr(7) = std_dev(survey_nr);
  //sdnr(8) = std_dev(cpue_nr);

  for (i=1;i<=nyrs_srv_ebs_r;i++)  // EBS survey likelihood
 {
   ii=yrs_srv_ebs(i);
   surv_like(2) += square(log(obs_srv_ebs(i)+1e-13) - log(pred_srv_ebs(ii)+1e-13) )/(2.*cv_srv_ebs(i)*cv_srv_ebs(i));
   rmse(2) += square(log(obs_srv_ebs(i)+1e-13) - log(pred_srv_ebs(ii)+1e-13));
 }  
   rmse(2) = sqrt(rmse(2)/nyrs_srv_ebs_r);
 
   



 obj_fun+= lambda(2)*surv_like(1) + lambda(15)*surv_like(2);
 
FUNCTION cat_likelihood  // fit the catches -------------------------------------------------------------
 catch_like=norm2(log(catch_bio(styr_fish,endyr_r)+0.0001) - log(pred_catch+0.00001));
 obj_fun+=lambda(4)*catch_like;
 
 
 
 

FUNCTION Fmort_pen  // Phases less than 2, penalize high F's ---------------------------------------------
 fpen = 0.0;
 //if (current_phase()<2)
 //  fpen=10.*norm2(mfexp(fmort_dev+log_avg_fmort)-1.0);
 //else
 //  fpen=.01*norm2(mfexp(fmort_dev+log_avg_fmort)-.2);

 //if (active(fmort_dev))
   fpen+= 0.1*norm2(fmort_dev);
 
 //fpen+=100*square(avg_fmort_dev);
 obj_fun+=fpen;
 

FUNCTION sel_likelihood // penalty for smoothness on time-varying sel parameters
 sel_like = 0.0;
 dvariable s = 0.0;
 dvar_matrix trans_log_sel_fish = trans(log_sel_fish);  

 // for logistic and double logistic curves with time-varying parameters, penalize the param devs
 if (active(a50_devs))
 {
     sel_like(1) = norm2(a50_devs)/(2.*sigma_a50*sigma_a50);   
 }

 if (active(aslope_devs))
 {
     sel_like(2) = norm2(aslope_devs)/(2.*sigma_aslp*sigma_aslp);
 }
 if (active(d50_devs))
 {
     sel_like(3) = norm2(d50_devs)/(2.*sigma_d50*sigma_d50);   
 }

 if (active(dslope_devs))
 {
     sel_like(4) = norm2(dslope_devs)/(2.*sigma_dslp*sigma_dslp);
 }
 // for double logistic  and bicubic spline, penalize the dome-shape, 

  if (sel_option==2 || sel_option==3 || sel_option==4)
   {
    for (i=sel_styr;i<=sel_endyr;i++)
     {
       for (j=1;j<=nselages-1;j++)
        {
         if (log_sel_fish(i,j)>log_sel_fish(i,j+1))
           {
             sel_like(5) += lambda(11)*square(log_sel_fish(i,j)-log_sel_fish(i,j+1));  // penalize the dome shape
           }
        }
      }
    }
  // for bicubic spline and cubic spline, penalize the smoothness across ages and years, and interannual differences 
 if(sel_option==3 || sel_option==4 )
  {                   // the smoothness penalty (across ages)
   for (i=sel_styr;i<=sel_endyr;i++)
     {
       s = mean(log_sel_fish(i));
       sel_like(9) += 10000*s*s;   
       dvar_vector df2 = first_difference(first_difference(log_sel_fish(i)));
       sel_like(6) += lambda(12)/nselages*df2*df2;    
     }
 
   for (j=1;j<=nselages;j++)
     {
       dvar_vector df1 = first_difference(trans_log_sel_fish(j));
       sel_like(7) += lambda(13)/(sel_endyr-sel_styr+1)*df1*df1;  // the penalty for interannual variation

       dvar_vector df2 = first_difference(df1);           // the penalty for smoothness over time
       sel_like(8) += lambda(14)/(sel_endyr-sel_styr+1)*df2*df2; 
     }
  }    

 obj_fun += sum(sel_like);
 

FUNCTION mat_likelihood  // fit the maturity curve
 mat_like = 0.0;
 int ii;

 for (i=1;i<=nages_T;i++)
 {
  ii = ages_T(i) - rec_age +1;
  mat_like += -0.01*mat_lambda(ii)*(T_y(i)*log(mat_theta(ii)) + (T_n(i) - T_y(i))*log(1.-mat_theta(ii)+1e-15)); // -ln like, Tenbrink data
 }

 for (i=1;i<=nages_S;i++)
 {
 ii = ages_S(i) - rec_age +1;
 mat_like += -0.01*mat_lambda(ii)*(S_y(i)*log(mat_theta(ii)) + (S_n(i) - S_y(i))*log(1.-mat_theta(ii)+1e-15)); // -ln like, Shawdata
 }
 

 obj_fun += mat_like;
 

FUNCTION dvar_vector cubic_spline(const dvar_vector& spline_coffs)
  {
  RETURN_ARRAYS_INCREMENT();
  int nodes=size_count(spline_coffs);
  dvector ia(1,nodes);
  dvector fa(1,nselages);
  ia.fill_seqadd(0,1./(nodes-1));
  fa.fill_seqadd(0,1./(nselages-1));
  vcubic_spline_function ffa(ia,spline_coffs);
  RETURN_ARRAYS_DECREMENT();
  
  //some testing here
  /*dvar_vector spline_nodes(1,nodes);
    spline_nodes.fill_seqadd(-0.5,1./(nodes-1));
    cout<<spline_nodes<<endl;
    vcubic_spline_function test_ffa(ia,spline_nodes);
    cout<<test_ffa(fa)<<endl;
    exit(1);*/
  return(ffa(fa));
  }
 



 //FUNCTION future_projections
 // Start calculation for first year of numbers at age matrix in projection
 //nage_future(styr_fut)(2,nages)=++elem_prod(natage(endyr)(1,nages-1),S(endyr)(1,nages-1));
 //nage_future(styr_fut,nages)+=natage(endyr,nages)*S(endyr,nages);
 // Set future total catch biomass to zero
 //catch_future = 0.;
 //biomass_future = 0.;
 //ssb_future = 0.;
 // Compute recent F levels
 //mean_recent_fs = mfexp(sum(fmort_dev(endyr-4,endyr))/5 + log_avg_fmort);

 
// Loop to cycle different fishing mortality values through the projections
 //for (int l=1;l<=num_proj_Fs;l++)
 // {
 //    switch(l)
 //      {
 //        case 1:
 //         ftmp = F40;
 //         break;
 //        case 2:
 //         ftmp = F40/2.0;
 //         break;
 //        case 3:
 //         ftmp = F30;
 //         break;
 //        case 4:
 //         ftmp = mean_recent_fs;
 //         break; 
 // 	 case 5:
 //         ftmp = 0.0;
 //         break;
 	
 //	}
   
  // Calculation of future F's, Z and survival (S)
  //Z_future = M;
  //for (i=styr_fut;i<=endyr_fut;i++)
  //{
  //  F_future(i) = sel_fish*ftmp;
  //  Z_future(i) +=F_future(i);
  //  S_future(i) = exp(-Z_future(i));
  //}
  // Calculation of future recruitment and spawners
  //Mean average recruitment of the time-series is used for projection
  //  dvariable Rectmp=mfexp(mean_log_rec);
  //  for (i=styr_fut;i<endyr_fut;i++)
  //  {
  //   nage_future(i,1) = Rectmp*mfexp(rec_dev_future(i));
    
    // Now graduate for the next year
   //   nage_future(i+1)(2,nages) = ++elem_prod(nage_future(i)(1,nages-1),S_future(i)(1,nages-1));
   //   nage_future(i+1,nages) += nage_future(i,nages)*S_future(i,nages);
   // }
   // nage_future(endyr_fut,1)= Rectmp*mfexp(rec_dev_future(i));
  
// Calculation of catch at predicted future age composition
    //for (i=styr_fut; i<=endyr_fut;i++)
    //{
      //catage_future(i) = 0.;
      //catage_future(i) += elem_prod(nage_future(i),
        //                  elem_prod(F_future(i),
        //                  elem_div ((1.-S_future(i)),Z_future(i))));
      //if (l!=num_proj_Fs) catch_future(l,i) += catage_future(i)*wt_pop;
      //biomass_future(l,i) += nage_future(i)*wt_pop;
      //ssb_future(l,i) += (nage_future(i)/2)*elem_prod(wt_pop,maturity);
    //}
 //}
FUNCTION get_age10
 // get the age at which the survey selectivity exceeds 10% (as potentially modified by the natural mortality rate)
  dvector tmp;
  
  tmp = value(sel_srv3) - 0.10;
  for (j=1;j<=nages;j++)
    {
      if(tmp(j)<=0) tmp(j) = 100;
    }

  for (j=1;j<=nages;j++)
    {
     if(abs(tmp(j) - min(tmp)<= 0.001))  firstage = ages(j); 
    }

    firstage += round(0.05/value(M));  // modify by the natural morality
    excludeage = firstage -1;          // exclude ages at and below excludeage 
     
FUNCTION update_compweights   // update if the comp is estmated, otherwise carry over previous comp weight
  // *_ta11 -- McAllister-Ianelli weights (method TA11 in Francis 2011)
  // *_ta12 -- weight by inverse of variance of normalized resids (Method TA1.2 in Francis 2011)
  // *_ta18 -- The weights for the Francis method

  if (lambda(5)>0 && nyrs_fish_unbiased_ac_r > 0 ) 
     { 
      if (har_flag==1) compweightsnew_ta11(1) = 1.0/mean(fac_mcian_wgt_inv);
      else   compweightsnew_ta11(1) = mean(fac_mcian_wgt);    
      compweightsnew_ta12(1) = 1./var(fac_nr);
        if (nyrs_fish_unbiased_ac_r >1) compweightsnew_ta18(1) = 1/var(fac_nr_fran);
        else compweightsnew_ta18(1) = 1/var(flc_nr_fran);   // if only one data point, pair with flc  
     }
   else 
   {
     compweightsnew_ta11(1) = compweights(1);
     compweightsnew_ta12(1) = compweights(1);
     compweightsnew_ta18(1) = compweights(1);
   }  
  
  if (lambda(6)>0 && nyrs_fish_lc_r > 0 ) 
     { 
      if (har_flag==1) compweightsnew_ta11(2) = 1.0/mean(flc_mcian_wgt_inv);
      else compweightsnew_ta11(2) = mean(flc_mcian_wgt);
      compweightsnew_ta12(2) = 1./var(flc_nr);
        if (nyrs_fish_lc_r >1) compweightsnew_ta18(2) = 1/var(flc_nr_fran);
        else compweightsnew_ta18(2) = 1/var(fac_nr_fran);  // if only one data point, pair with fac 
     }
   else 
   {
     compweightsnew_ta11(2) = compweights(2);
     compweightsnew_ta12(2) = compweights(2);
     compweightsnew_ta18(2) = compweights(2);
   }
   
   if (lambda(7)>0 && nyrs_surv3_unbiased_ac_r >0 ) 
     { 
      if (har_flag==1) compweightsnew_ta11(3) = 1.0/mean(sac_mcian_wgt_inv);
      else compweightsnew_ta11(3) = mean(sac_mcian_wgt);
      compweightsnew_ta12(3) = 1./var(sac_nr);
        if (nyrs_surv3_unbiased_ac_r >1) compweightsnew_ta18(3) = 1/var(sac_nr_fran);
        else compweightsnew_ta18(3) = 1/var(slc_nr_fran);  // if only one data point, pair with slc   
     }
   else 
   {
     compweightsnew_ta11(3) = compweights(3);
     compweightsnew_ta12(3) = compweights(3);
     compweightsnew_ta18(3) = compweights(3);
   }
  
  if (lambda(8)>0 && nyrs_surv3_lc_r > 0) 
     {
      if (har_flag==1) compweightsnew_ta11(4) = 1.0/mean(slc_mcian_wgt_inv);
      else compweightsnew_ta11(4) = mean(slc_mcian_wgt);  
      compweightsnew_ta12(4) = 1./var(slc_nr);
        if (nyrs_surv3_lc_r >1) compweightsnew_ta18(4) = 1/var(slc_nr_fran);
        else compweightsnew_ta18(4) = 1/var(sac_nr_fran);  // if only one data point, pair with sac     
     }
   else 
   {
     compweightsnew_ta11(4) = compweights(4);
     compweightsnew_ta12(4) = compweights(4);
     compweightsnew_ta18(4) = compweights(4);
   }
  
   if (lambda(9)>0 && nyrs_surv_ebs_unbiased_ac_r>0) 
     {
      if (har_flag==1) compweightsnew_ta11(5) = 1.0/mean(ebs_sac_mcian_wgt_inv);
      else compweightsnew_ta11(5) = mean(ebs_sac_mcian_wgt);  
      compweightsnew_ta12(5) = 1./var(ebs_sac_nr);
        if (nyrs_surv_ebs_unbiased_ac_r >1) compweightsnew_ta18(5) = 1/var(ebs_sac_nr_fran);
        else compweightsnew_ta18(5) = 1/var(sac_nr_fran);   // if only one data point, pair with ai sac   
     }
   else 
   {
     compweightsnew_ta11(5) = compweights(5);
     compweightsnew_ta12(5) = compweights(5);
     compweightsnew_ta18(5) = compweights(5);
   }

   if (lambda(10)>0 && nyrs_surv_ebs_lc_r > 0 ) 
     {
      if (har_flag==1) compweightsnew_ta11(6) = 1.0/mean(ebs_slc_mcian_wgt_inv);
      else compweightsnew_ta11(6) = mean(ebs_slc_mcian_wgt);  
      compweightsnew_ta12(6) = 1./var(ebs_slc_nr);
        if (nyrs_surv_ebs_lc_r >1) compweightsnew_ta18(6) = 1/var(ebs_slc_nr_fran);
        else compweightsnew_ta18(6) = 1/var(ebs_sac_nr_fran);   // if only one data point, pair with ebs sac   
     }
   else 
   {
     compweightsnew_ta11(6) = compweights(6);
     compweightsnew_ta12(6) = compweights(6);
     compweightsnew_ta18(6) = compweights(6);
   }
       

REPORT_SECTION //-------------------------------------------------------------------------------------------
 rescaled_F = value(mfexp(log_avg_fmort + fmort_dev));

  for (i=styr;i<=endyr_r;i++)
      {
        rescaled_sel_fish(i) = value(sel_fish(i)(1,nages_dat));
        rescaled_F(i) = rescaled_F(i)*max(rescaled_sel_fish(i));
        rescaled_sel_fish(i) = rescaled_sel_fish(i)/max(rescaled_sel_fish(i));
      }

  for (j=1;j<=nages_dat;j++)
       recent_fish_sel(j) = mean(column(rescaled_sel_fish,j)(endyr_r-4,endyr_r));
           

  // the last year of recruitment to meet the survey a10 condition 
  double lastyr_rec_a10;
  lastyr_rec_a10 = endyr_r - max(fixedrec_yrs,excludeage - (rec_age - 1));

  // numbers at age with mean recruitments for the survey age10 year classes
   dvar_matrix natage_mean(lastyr_rec_a10+1,endyr_r+1,1,nages);    
  
   //natage_mean = get_projection_numbers_at_age();

 SB0 = get_spr(0.0);
 F40 = get_spr_rates(0.4);
 F35 = get_spr_rates(0.35);
 F30 = get_spr_rates(0.30);
 SBF40 = get_spr(F40);
 SBF35 = get_spr(F35);
 SBF30 = get_spr(F30);

 report << "the firstage is "<< firstage <<endl;
 report << "the excludeage is "<< excludeage <<endl;
 report << "the last year of recruitment is "  <<lastyr_rec << endl;
 report << "the last year of recruitment for a10 is "  << lastyr_rec_a10 << endl;
 report << "the first year where we do the new projection is " << lastyr_rec_a10+1 << endl; 
 report << "Total number of fish: years " <<styr<<" to " << endyr_r+1<< endl;
 	report << rowsum(natage) << endl;
 	report << "Numbers of fish: ages "  <<ages_dat(1)<<" to " << ages_dat(nages_dat)<< endl;
 	for (i=styr;i<= endyr_r+1;i++)
      report << i <<" "<<natage_bin(i) << endl;
	report << "Number in first year: ages "<<ages_dat(1)<<" to " << ages_dat(nages_dat)<< endl;
  report << natage_bin(styr) << endl;
  report << "Number in end year: ages " <<ages_dat(1)<<" to " << ages_dat(nages_dat)<< endl;
 	report << natage_bin(endyr_r) << endl;
  report << "Recruitments (age "<<rec_age<<"): years "<<styr<<" to " << lastyr_rec << endl;
 	report << est_rec(styr,lastyr_rec) << endl;
 	report << "Spawner biomass: years "<<styr<<" to " << lastyr_rec << endl;
 	report << est_spb(styr,lastyr_rec) << endl;
 	report << "time series spawner biomass: years "<<styr<<" to " << endyr_r << endl;
 	report << sp_biom(styr,endyr_r) << endl;
  	//report << "time series eggs: seq(1960,2004)" << endl;
 	//report << eggs(styr,endyr) << endl;
        //report << "time series larvae: seq(1960,2004)" << endl;
 	//report << larvae(styr,endyr) << endl;

	report << "SR curve SSB: seq(1,20)" << endl;
 	report << SRec_spawn << endl;
 	report << "SR curve recs: seq(1,20)"  << endl;
 	report << SRec_rec << endl;
  report << "AI Survey selectivity: "<<ages_dat(1)<<" to " << ages_dat(nages_dat)<< endl;
  report << sel_srv3(1,nselages) << endl;
  report << "ebs Survey selectivity: "<<ages_dat(1)<<" to " << ages_dat(nages_dat)<< endl;
  report << sel_srv_ebs(1,nselages) << endl;
  report << "Fishing mortality: years "<<styr<<" to " << endyr_r << endl;
 	report << rescaled_F << endl;
 	//report << "Observed CPUE: years " <<yrs_cpue << endl;
 	//report << obs_cpue << endl;
 	//report << "Predicted CPUE: years " <<yrs_cpue << endl;
  //	report << pred_cpue << endl;
	report << "Observed srv3: years " <<yrs_srv3(1,nyrs_srv3_r) << endl;
 	report << obs_srv3(1,nyrs_srv3_r) << endl;
 	report << "Observed srv3 lower: years " <<yrs_srv3(1,nyrs_srv3_r) << endl;
 	report << obs_srv3_lower(1,nyrs_srv3_r) << endl;
	report << "Observed srv3 upper: years " <<yrs_srv3(1,nyrs_srv3_r) << endl;
 	report << obs_srv3_upper(1,nyrs_srv3_r) << endl;
	report << "Predicted srv3: years "<<styr_fish<<" to " << endyr_r << endl;
 	report << pred_srv3 << endl;
 	report << "Observed srv_ebs: years " <<yrs_srv_ebs(1,nyrs_srv_ebs_r) << endl;
  report << obs_srv_ebs(1,nyrs_srv_ebs_r) << endl;
  report << "Observed srv_ebs lower: years " <<yrs_srv_ebs(1,nyrs_srv_ebs_r) << endl;
  report << obs_srv_ebs_lower(1,nyrs_srv_ebs_r) << endl;
  report << "Observed srv_ebs upper: years " <<yrs_srv_ebs(1,nyrs_srv_ebs_r) << endl;
  report << obs_srv_ebs_upper(1,nyrs_srv_ebs_r) << endl;
  report << "Predicted srv_ebs: years "<<styr_fish<<" to " << endyr_r << endl;
  report << pred_srv_ebs << endl;
  report << "Total biomass: years "<<styr_rec-rec_age<<" to " << endyr_r+1 << endl;
 	report << totbiom << endl;
 	report << "Observed catch Biomass: years "<<styr_fish<<" to " << endyr_r << endl;
	report << catch_bio(styr_fish,endyr_r) << endl;
 	report << "Predicted catch Biomass: years "<<styr_fish<<" to " << endyr_r << endl;
	report << pred_catch << endl;
 	report << "Estimated historical catch: 'est_hist_catch'"  << endl;
 	report << ehc << endl;
 	
 //	report << "Observed Prop(fishery biased): 'year', 'effn', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25'" << endl;
 //       for (i=1;i<=nyrs_fish_biased_ac; i++)
 //	{
 //         if (fish_biased_ac_samp(i)>1)
 //           {
 //	      report << yrs_fish_biased_ac(i)<<" "<<effn(1,i)<<"  "<<oac_fish_biased(i)<< endl;
 //	    }
 //	}
 //	report << "Predicted Prop(fishery biased): 'year', 'effn','3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25'" << endl;
 //	for (i=1;i<=nyrs_fish_biased_ac; i++)
 //	{
 //	  report << yrs_fish_biased_ac(i)  <<" "<<effn(1,i)<<" "<<eac_fish_biased(yrs_fish_biased_ac(i))<< endl;
 //	}

 	report << "Observed Prop(fishery lengths): year, effn, lengths " <<lengths<< endl;
        for (i=1;i<=nyrs_fish_lc_r; i++)
 	{          
 	  if (fish_lc_samp_r(i)>1)
            {
            if (yrs_fish_lc_r(i) == 1978 && yrs_fish_lc_r(i+1) == 1978)
                  report << "1978.1" <<" "<<effn(2,i)<<" "<<olc_fish(i)<< endl;
 	    else if (yrs_fish_lc_r(i) == 1978 && yrs_fish_lc_r(i-1) == 1978)
                  report << "1978.2" <<" "<<effn(2,i)<<" "<<olc_fish_r(i)<< endl;
            else
                  report << yrs_fish_lc_r(i)<<" "<<effn(2,i)<<"  "<<olc_fish_r(i)<< endl;
            }
 	}
 	report << "Predicted Prop(fishery lengths): year, effn, lengths " <<lengths<< endl;
  	for (i=1;i<=nyrs_fish_lc_r; i++)
 	{
            if (yrs_fish_lc_r(i) == 1978 && yrs_fish_lc_r(i+1) == 1978)
                  report << "1978.1" <<" "<<effn(2,i)<<"  "<<elc_fish(yrs_fish_lc(i))<< endl;
 	    else if (yrs_fish_lc_r(i) == 1978 && yrs_fish_lc_r(i-1) == 1978)
                  report << "1978.2" <<" "<<effn(2,i)<<"  "<<elc_fish(yrs_fish_lc(i))<< endl;
            else
                  report << yrs_fish_lc_r(i)<<" "<<effn(2,i)<<"  "<<elc_fish(yrs_fish_lc_r(i))<< endl;
 	}

 	report << "Observed Prop(fishery unbiased): year, effn, ages " <<ages_dat<< endl;
        for (i=1;i<=nyrs_fish_unbiased_ac_r; i++)
 	{
          if (fish_unbiased_ac_samp_r(i)>1)
            {
 	      report << yrs_fish_unbiased_ac_r(i)<<" "<<effn(1,i)<<" "<<oac_fish_unbiased_r(i)<< endl;
 	    }
 	}
 	report << "Predicted Prop(fishery unbiased): year, effn, ages " <<ages_dat<< endl;
  	for (i=1;i<=nyrs_fish_unbiased_ac_r; i++)
 	{
 	  report << yrs_fish_unbiased_ac_r(i)  <<" "<<effn(1,i)<<" "<<eac_fish_unbiased_dat(yrs_fish_unbiased_ac_r(i))<< endl;
 	}

	if (nyrs_surv3_unbiased_ac_r > 0) 
        {
        report << "Observed Prop(AI survey unbiased): year, effn, ages " <<ages_dat<< endl;
        for (i=1;i<=nyrs_surv3_unbiased_ac_r; i++)
 	{
          if (surv3_unbiased_ac_samp_r(i)>1)
            {
 	      report << yrs_surv3_unbiased_ac_r(i)<<" "<<effn(3,i)<<" "<<oac_surv3_unbiased_r(i)<< endl;
 	    }
 	}
 	report << "Predicted Prop(AI survey unbiased): year, effn, ages " <<ages_dat<< endl;
  	for (i=1;i<=nyrs_surv3_unbiased_ac_r; i++)
 	{
 	  if (surv3_unbiased_ac_samp_r(i)>1)
            {
             report << yrs_surv3_unbiased_ac_r(i)  <<" "<<effn(3,i)<<" "<<eac_surv3_unbiased_dat(yrs_surv3_unbiased_ac_r(i))<< endl;
 	    }
        }
        }
	  report << "Observed Prop(AI survey lengths): year, effn, lengths " <<lengths<< endl;
          for (i=1;i<=nyrs_surv3_lc_r; i++)
 	  {          
 	        report << yrs_surv3_lc_r(i)<<" "<<effn(4,i)<<"  "<<olc_surv3_r(i)<< endl;
          }
 	  report << "Predicted Prop(AI survey lengths): year, effn, lengths " <<lengths<< endl;
  	  for (i=1;i<=nyrs_surv3_lc_r; i++)
 	  {
              report << yrs_surv3_lc_r(i)<<" "<<effn(4,i)<<"  "<<elc_surv3(yrs_surv3_lc_r(i))<< endl;
 	  }

    if (nyrs_surv_ebs_unbiased_ac_r > 0) 
        {
        report << "Observed Prop(EBS survey unbiased): year, effn, ages " <<ages_dat<< endl;
        for (i=1;i<=nyrs_surv_ebs_unbiased_ac_r; i++)
  {
          if (surv_ebs_unbiased_ac_samp_r(i)>1)
            {
        report << yrs_surv_ebs_unbiased_ac_r(i)<<" "<<effn(5,i)<<" "<<oac_surv_ebs_unbiased_r(i)<< endl;
      }
  }


  report << "Predicted Prop(EBS survey unbiased): year, effn, ages " <<ages_dat<< endl;
    for (i=1;i<=nyrs_surv_ebs_unbiased_ac_r; i++)
  {
    if (surv_ebs_unbiased_ac_samp_r(i)>1)
            {
             report << yrs_surv_ebs_unbiased_ac_r(i)  <<" "<<effn(5,i)<<" "<<eac_surv_ebs_unbiased_dat(yrs_surv_ebs_unbiased_ac_r(i))<< endl;
      }
        }
        }

         report << "Observed Prop(EBS survey lengths): year, effn, lengths " <<lengths<< endl;
          for (i=1;i<=nyrs_surv_ebs_lc_r; i++)
    {          
          report << yrs_surv_ebs_lc_r(i)<<" "<<effn(6,i)<<"  "<<olc_surv_ebs_r(i)<< endl;
          }
    report << "Predicted Prop(EBS survey lengths): year, effn, lengths " <<lengths<< endl;
      for (i=1;i<=nyrs_surv_ebs_lc_r; i++)
    {
              report << yrs_surv_ebs_lc_r(i)<<" "<<effn(6,i)<<"  "<<elc_surv_ebs(yrs_surv_ebs_lc_r(i))<< endl;
    }
 		
	report << "the mean effective n by data type: 'fish_ac', 'fish_lc', 'AI_surv_ac', 'AI_surv_lc','ebs_surv_ac', 'ebs_surv_lc' " << endl; 
  if(nyrs_fish_unbiased_ac_r>0)
            report <<" "<<(sum(effn(1)))/nyrs_fish_unbiased_ac_r;
        else 
            report <<" "<<0; 

        if(nyrs_fish_lc_r>0)
            report <<" "<<(sum(effn(2)))/nyrs_fish_lc_r;
        else 
            report <<" "<<0;     

        if(nyrs_surv3_unbiased_ac_r>0)
            report <<" "<<(sum(effn(3)))/nyrs_surv3_unbiased_ac_r;
        else 
            report <<" "<<0;

        if(nyrs_surv3_lc_r>0)
            report <<" "<<(sum(effn(4)))/nyrs_surv3_lc_r;
        else 
            report <<" "<<0;            

        if(nyrs_surv_ebs_unbiased_ac_r>0)
            report <<" "<<(sum(effn(5)))/nyrs_surv_ebs_unbiased_ac_r;
        else 
            report <<" "<<0;

        if(nyrs_surv_ebs_lc_r>0)
            report <<" "<<(sum(effn(6)))/nyrs_surv_ebs_lc_r << endl;
        else 
            report <<" "<<0 << endl; 
  
 	report << "the sample weights by data type: 'fish_ac', 'fish_lc', 'AI_surv_ac', 'AI_surv_lc','ebs_surv_ac', 'ebs_surv_lc'" << endl; 
  if(nyrs_fish_unbiased_ac_r>0)
            report <<" "<<sum(fish_unbiased_ac_samp_r)/nyrs_fish_unbiased_ac_r;
        else 
            report <<" "<<0; 

        if(nyrs_fish_lc_r>0)
            report <<" "<<sum(fish_lc_samp_r)/nyrs_fish_lc_r;
        else 
            report <<" "<<0;     

        if(nyrs_surv3_unbiased_ac_r>0)
            report <<" "<<sum(surv3_unbiased_ac_samp_r)/nyrs_surv3_unbiased_ac_r;
        else 
            report <<" "<<0;

        if(nyrs_surv3_lc_r>0)
            report <<" "<<sum(surv3_lc_samp_r)/nyrs_surv3_lc_r;
        else 
            report <<" "<<0;            

        if(nyrs_surv_ebs_unbiased_ac_r>0)
            report <<" "<<sum(surv_ebs_unbiased_ac_samp_r)/nyrs_surv_ebs_unbiased_ac_r;
        else 
            report <<" "<<0;

        if(nyrs_surv_ebs_lc_r>0)
            report <<" "<<sum(surv_ebs_lc_samp_r)/nyrs_surv_ebs_lc_r << endl;
        else 
            report <<" "<<0 << endl;       
  


  report << "the rmse for the surveys and rec: AI','ebs','rec'"<<endl;;
  report << rmse << endl;

	report << "The likelhood components: 'histFpen', 'selpen'', 'sprpen', 'rec_likelihood', 'AI_srv', 'ebs_surv','cat_likelihood', 'Fmort_pen','fish_unbiased_ac','fish_lc', 'AI_surv_ac', 'AI_surv_lc','AI_surv_ac',";
  report <<" 'AI_surv_lc','ebs_surv_ac','ebs_surv_lc', 'mat_like'"<<endl;
 	report << hf_pen<<" "<<sel_like<<"  "<<sprpen<<" "<<sum(rec_like)<<"  "<<surv_like(1)<<" "<<surv_like(2)<<" "<<catch_like<<" "<<fpen<<" "<<age_like<<" "<<mat_like<<endl;
  report << "the prior components of the like: 'prior_m', 'prior_qsrv' " << endl;
  report << priormq(1)<<" "<<priormq(2)<<endl;
  report << " the weights for the fac, flc, sac, slc,ebs_sac, ebs_slc are "<< endl;
  report << compweightsnew_ta12<< endl; 
  report << "estimates of M ,AI q, and EBS q " << endl;
  report << M <<" "<<q_srv3<<" "<<q_srv_ebs <<   endl;
  report << "the sdnr for the fac, flc, sac, slc, ebs_sac, ebs_slc, ai_survey are   "<< endl;
  report << sdnr(1)<<" "<<sdnr(2)<<" "<<sdnr(3) << " "<<sdnr(4)<<" "<<sdnr(5)<<" "<<sdnr(6)<<" "<<sdnr(7)<< endl;
  report << " the McAllister- Ianelli weights for the fac, flc, sac, slc, ebs_ac, ebs_lc are "<< endl;
  report <<compweightsnew_ta11 <<  endl; 
  report << " the Francis weights for the fac, flc, sac, slc, ebs_ac, ebs_lc  are "<< endl;
  report <<compweightsnew_ta18 << endl;  

  report <<" the root mean square error for the fac, flc, sac, slc, ebs_sac, ebs_slc are" << endl;
  report << al_rmse(1) <<" "<<al_rmse(2)<<" "<<al_rmse(3)<<" "<<al_rmse(4)<<" "<<al_rmse(5)<<" "<<al_rmse(6) <<endl;
  report << "the fishery selectivity by year and age "<< endl;   
  report << rescaled_sel_fish << endl;
  report << "the recent selectivity (last five years)" << endl;
  report << recent_fish_sel << endl;
  report << "F spr rates " << endl;
  report << F40 <<" "<< F35 << endl;
  report << " the objective function is " << endl;
  report << obj_fun << endl;
  report << " the cubic spline parameter matrix is "<< endl;
  report << sel_par << endl;
  # include "pop-s-report_exp_r.cxx"   // ADMB code to write the S-compatible report
 	
  ofstream projfile("popproj.dat");
  projfile << "BSAI_pop"  << endl; 
  projfile << 0 <<" # SSL Species???" << endl;
  projfile << 0 <<" # Constant buffer of Dorn" << endl;
  projfile << 1 <<" # Number of fsheries" << endl;
  projfile << 1 <<" # Number of sexes??" << endl;
  projfile << mean(mfexp(log_avg_fmort + fmort_dev(endyr_r-5,endyr_r-1))) << " #  average 5 yr f " << endl;
  projfile << 1 <<" # author f" << endl;
  projfile << 0.4 <<" # ABC SPR" << endl;
  projfile << 0.35 <<" # MSY SPR" << endl;
  projfile << spawn_mo <<" # Spawnmo" << endl;
  projfile << nages_dat <<" # Number of ages" << endl;
  projfile << 1 <<" # Fratio" << endl;
  Mvec = M;
  projfile << " # Natural mortality " << endl; 
  projfile << Mvec(1,nages_dat) << endl;
  projfile << " # Maturity " << endl; 
  projfile << maturity_bin << endl;       
  projfile << " # Wt Spawn " << endl; 
  projfile << wt_pop_bin << endl;	
  projfile << " # Wt Fish " << endl; 
  projfile << wt_fsh_bin << endl;
  projfile << " # selectivity " << endl; 
  projfile << recent_fish_sel << endl;
  projfile << " # natage " << endl; 
  projfile << natage_bin(endyr_r) << endl; 
  projfile << " # Nrec " << endl; 
  projfile << lastyr_rec - (max(1977,styr) +rec_age -1)  << endl;	
  projfile << " # rec " << endl; 
  projfile << column(natage_bin,1)(max(1977,styr) +rec_age,lastyr_rec) << endl;
  projfile << " # ssb " << endl; 
  projfile << sp_biom(max(1977,styr),lastyr_rec-rec_age) << endl;

  ofstream projfile2("popproj_age10.dat");
  projfile2 << "BSAI_pop"  << endl; 
  projfile2 << 0 <<" # SSL Species???" << endl;
  projfile2 << 0 <<" # Constant buffer of Dorn" << endl;
  projfile2 << 1 <<" # Number of fsheries" << endl;
  projfile2 << 1 <<" # Number of sexes??" << endl;
  projfile2 << mean(mfexp(log_avg_fmort + fmort_dev(endyr_r-5,endyr_r-1))) << " #  average 5 yr f " << endl;
  projfile2 << 1 <<" # author f" << endl;
  projfile2 << 0.4 <<" # ABC SPR" << endl;
  projfile2 << 0.35 <<" # MSY SPR" << endl;
  projfile2 << spawn_mo <<" # Spawnmo" << endl;
  projfile2 << nages_dat <<" # Number of ages" << endl;
  projfile2 << 1 <<" # Fratio" << endl;
  Mvec = M;
  projfile2 << " # Natural mortality " << endl; 
  projfile2 << Mvec(1,nages_dat) << endl;
  projfile2 << " # Maturity " << endl; 
  projfile2 << maturity_bin << endl;       
  projfile2 << " # Wt Spawn " << endl; 
  projfile2 << wt_pop_bin << endl; 
  projfile2 << " # Wt Fish " << endl; 
  projfile2 << wt_fsh_bin << endl;
  projfile2 << " # selectivity " << endl; 
  projfile2 << recent_fish_sel << endl;
  projfile2 << " # natage " << endl; 
  projfile2 << natage_bin(endyr_r) << endl; 
  projfile2 << " # Nrec " << endl; 
  projfile2 << lastyr_rec_a10 - (max(1977,styr) +rec_age -1)  << endl; 
  projfile2 << " # rec " << endl; 
  projfile2 << column(natage_bin,1)(max(1977,styr) +rec_age,lastyr_rec_a10) << endl;
  projfile2 << " # ssb " << endl; 
  projfile2 << sp_biom(max(1977,styr),lastyr_rec_a10-rec_age) << endl;
 	
  ofstream compweightsnew_ta12_file("compweights_new_ta12.ctl");    // new comp weights based on method TA1.2 in Francis 2011
  compweightsnew_ta12_file << compweightsnew_ta12 << endl; 
  
  ofstream compweightsnew_ta18_file("compweights_new_ta18.ctl");    // new comp weights based on Francis method (TA1.8)
  compweightsnew_ta18_file << compweightsnew_ta18 << endl; 

  ofstream compweightsnew_ta11_file("compweights_new_ta11.ctl");    // new comp weights McAllister-Ianelli method (TA1.1)
  compweightsnew_ta11_file << compweightsnew_ta11 << endl; 



 	
RUNTIME_SECTION //------------------------------------------------------------------------------------------
    convergence_criteria 1.e-4 1.e-4 1.e-4 1.e-7 1.e-7 1.e-7
    maximum_function_evaluations 1000, 1000, 1000, 10000, 20000, 20000

TOP_OF_MAIN_SECTION
  arrmblsize = 1000000;

GLOBALS_SECTION
 # include "admodel.h"          // Include AD class definitions
 # include "mhp-s-funcs.cpp"    // Include S-compatible output functions (needs preceding)

 void function_minimizer::mcmc_eval(void)
        {
                // |---------------------------------------------------------------------------|
                // | Added DIC calculation.  Martell, Jan 29, 2013                             |
                // |---------------------------------------------------------------------------|
                // | DIC = pd + dbar
                // | pd  = dbar - dtheta  (Effective number of parameters)
                // | dbar   = expectation of the likelihood function (average f)
                // | dtheta = expectation of the parameter sample (average y)

          gradient_structure::set_NO_DERIVATIVES();
          initial_params::current_phase=initial_params::max_number_phases;
          uistream * pifs_psave = NULL;

        #if defined(USE_LAPLACE)
        #endif

        #if defined(USE_LAPLACE)
            initial_params::set_active_random_effects();
            int nvar1=initial_params::nvarcalc();
        #else
          int nvar1=initial_params::nvarcalc(); // get the number of active parameters
        #endif
          int nvar;

          pifs_psave= new
            uistream((char*)(ad_comm::adprogram_name + adstring(".psv")));
          if (!pifs_psave || !(*pifs_psave))
          {
            cerr << "Error opening file "
                    << (char*)(ad_comm::adprogram_name + adstring(".psv"))
               << endl;
            if (pifs_psave)
            {
              delete pifs_psave;
              pifs_psave=NULL;
              return;
            }
          }
          else
          {
            (*pifs_psave) >> nvar;
            if (nvar!=nvar1)
            {
              cout << "Incorrect value for nvar in file "
                   << "should be " << nvar1 << " but read " << nvar << endl;
              if (pifs_psave)
              {
                delete pifs_psave;
                pifs_psave=NULL;
              }
              return;
            }
          }

          int nsamp = 0;
          double sumll = 0;
          independent_variables y(1,nvar);
          independent_variables sumy(1,nvar);

          do
          {
            if (pifs_psave->eof())
            {
              break;
            }
            else
            {
              (*pifs_psave) >> y;
              sumy = sumy + y;
              if (pifs_psave->eof())
              {
                double dbar = sumll/nsamp;
                int ii=1;
                y = sumy/nsamp;
                initial_params::restore_all_values(y,ii);
                initial_params::xinit(y);
                double dtheta = 2.0 * get_monte_carlo_value(nvar,y);
                double pd     = dbar - dtheta;
                double dic    = pd + dbar;
                double dicValue      = dic;
                double dicNoPar      = pd;

                cout<<"Number of posterior samples    = "<<nsamp    <<endl;
                cout<<"Expectation of log-likelihood  = "<<dbar     <<endl;
                cout<<"Expectation of theta           = "<<dtheta   <<endl;
                cout<<"Number of estimated parameters = "<<nvar1    <<endl;
                    cout<<"Effective number of parameters = "<<dicNoPar <<endl;
                    cout<<"DIC                            = "<<dicValue <<endl;
                break;
              }
              int ii=1;
              initial_params::restore_all_values(y,ii);
              initial_params::xinit(y);
              double ll = 2.0 * get_monte_carlo_value(nvar,y);
              sumll    += ll;
              nsamp++;
              // cout<<sumy(1,3)/nsamp<<" "<<get_monte_carlo_value(nvar,y)<<endl;
            }
          }
          while(1);
          if (pifs_psave)
          {
            delete pifs_psave;
            pifs_psave=NULL;
          }
          return;
        }



   
 
 

       
