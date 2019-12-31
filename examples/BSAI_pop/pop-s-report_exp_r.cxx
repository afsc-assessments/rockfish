  //------- pop<nn>-s-report.cxx-----------------------------------------------
  // ------ Write results into S-compatible file -- PDS July 2003   -------
  // ------ Based on work from Michael Prager
  // ------ The file will be named <name>.rdat, where <name> is the TPL file name
  // ------ The file can be read by S with the statement: x <- dget("<name>.rdat")

  const char* cc = ", ";    // To reduce clutter later on
  int i, il, y;                // Vars used as indices
  
  // Open the S output file
  ofstream sfile ((char*)(adprogram_name + ".rdat"));

  // ------------------ START OVERALL S STRUCTURE (LIST) -------------------------

  sfile << "structure(list(";

   
  // ------------------ TIME SERIES DATA, 1960-2002 -------------------------

  // Write start of data frame for continous time series from 1960-2002
      sfile << "t.series = structure(list(" << endl;

  // Write vector of years:
      for(y=styr; y<=endyr_r; y++) {xdum2[y] = double(y);}
      write_s_yrcol(sfile, value(xdum2), "year", styr, endyr_r, 0);
      
  // Write time series data (1960-2002)
      write_s_yrcol(sfile, value(rowsum(natage)), "totnum", styr, endyr_r, 0);
      write_s_yrcol(sfile, value(column(natage,1)(styr,endyr_r)), "a3recs", styr, endyr_r, 0);
      write_s_yrcol(sfile, value(sp_biom), "spbiom", styr, endyr_r, 0);
      write_s_yrcol(sfile, rescaled_F, "fmort", styr, endyr_r, 0);
      write_s_yrcol(sfile, value(pred_srv3), "aipredsrv", styr, endyr_r, 0);
      write_s_yrcol(sfile, value(pred_srv_ebs), "ebspredsrv", styr, endyr_r, 0);
      write_s_yrcol(sfile, value(totbiom), "totbiom", styr, endyr_r, 0);
      write_s_yrcol(sfile, biomass2016(1960,2016), "biomass2016", 1960, 2016, 0);
      write_s_yrcol(sfile, catch_bio, "catchbio", styr, endyr_r, 1);
      write_s_rownames_int(sfile, styr, endyr_r, "df", 0);

  // ------------------ START N-AT-AGE MATRIX -------------------------
  // Matrix of numbers at age, with ages as columns & years as rows

  write_s_matrix(sfile, value(natage_bin), "natage", styr, endyr_r, 1, nages_dat);
  write_s_rownames_vec(sfile, ages_dat, 1, nages_dat, "ma", 0);  // this writes COLUMN names

 // ---------------------------PREDICTED RECRUITMENT CURVE -------------------- -------------
 // Write start of data frame for predicted recruitment curve
      
      sfile << "reccurve = structure(list(" << endl;     
      write_s_yrcol(sfile, value(SRec_spawn), "Srec.spawn", 1, 20, 0);
      write_s_yrcol(sfile, value(SRec_rec), "Srec.rec", 1, 20, 1);
      write_s_rownames_int(sfile, 1, 20, "df", 0);

// ------------------ FISHERY SELECTIVITY MATRIX -------------------------
// Matrix of fishery selectivity by year, with ages as columns & years as rows

  write_s_matrix(sfile, rescaled_sel_fish, "selfish", styr, endyr_r, 1, nages_dat);
  write_s_rownames_vec(sfile, ages_dat, 1, nages_dat, "ma", 0);  // this writes COLUMN names

 // ------------------------------SURVEY SELECTIVITY CURVES ----------------------------------------
 // write start of survey age-based selectivity curves

    sfile << "ai_srvsel = structure(list(" << endl;
    write_s_yrcol(sfile, value(sel_srv3(1,nages_dat)), "ai_srvsel", 1, nages_dat, 1);
    write_s_rownames_vec(sfile, ages_dat, 1, nages_dat, "df", 0);

    sfile << "ebs_srvsel = structure(list(" << endl;
    write_s_yrcol(sfile, value(sel_srv_ebs(1,nages_dat)), "ebs_srvsel", 1, nages_dat, 1);
    write_s_rownames_vec(sfile, ages_dat, 1, nages_dat, "df", 0);

 // ------------------ CPUE TIME SERIES DATA, 1968-1979 -------------------------
 //     sfile << "obscpue = structure(list(" << endl;
 //     write_s_yrcol(sfile, obs_cpue, "obscpue", 1, nyrs_cpue, 0);
 //     write_s_yrcol(sfile, value(pred_cpue), "predcpue", 1, nyrs_cpue, 1);
 //     write_s_rownames_vec(sfile, yrs_cpue, 1, nyrs_cpue, "df", 0);


 // ------------------ OBSERVED TRAWL SURVEY TIME SERIES DATA, 1982-2002 -------------------------
       sfile << "ai_obssrv = structure(list(" << endl;
        write_s_yrcol(sfile, obs_srv3(1,nyrs_srv3_r), "obssrv", 1, nyrs_srv3_r, 0);
        write_s_yrcol(sfile, obs_srv3_lower(1,nyrs_srv3_r), "obssrv.lower", 1, nyrs_srv3_r, 0);
        write_s_yrcol(sfile, obs_srv3_upper(1,nyrs_srv3_r), "obssrv.upper", 1, nyrs_srv3_r, 1);
        write_s_rownames_vec(sfile, yrs_srv3(1,nyrs_srv3_r), 1, nyrs_srv3_r, "df", 0);

        sfile << "ebs_obssrv = structure(list(" << endl;
        write_s_yrcol(sfile, obs_srv_ebs(1,nyrs_srv_ebs_r), "obssrv", 1, nyrs_srv_ebs_r, 0);
        write_s_yrcol(sfile, obs_srv_ebs_lower(1,nyrs_srv_ebs_r), "obssrv.lower", 1, nyrs_srv_ebs_r, 0);
        write_s_yrcol(sfile, obs_srv_ebs_upper(1,nyrs_srv_ebs_r), "obssrv.upper", 1, nyrs_srv_ebs_r, 1);
        write_s_rownames_vec(sfile, yrs_srv_ebs(1,nyrs_srv_ebs_r), 1, nyrs_srv_ebs_r, "df", 0);

 // ----------------------------AGE AND LENGTH COMPOSITIONS ----------------------------------------
 // write start of age composition matrices

 // fishery biased ages
 //   write_s_imatrix(sfile, value(oac_fish_biased), yrs_fish_biased_ac, "oac.f.biased", 1,nyrs_fish_biased_ac, 1, nages);
 //   write_s_rownames_vec(sfile, ages, 1, nages, "ma", 0);  // this writes COLUMN names 

 //   write_s_matrix(sfile, value(eac_fish_biased), "eac.f.biased", styr_fish,endyr, 1, nages);
 //   write_s_rownames_vec(sfile, ages, 1, nages, "ma", 0);  // this writes COLUMN names

// fishery unbiased ages
    if(nyrs_fish_unbiased_ac_r>0)
    {
      write_s_imatrix(sfile, value(oac_fish_unbiased_r), yrs_fish_unbiased_ac_r, "oac.f.unbiased", 1,nyrs_fish_unbiased_ac_r, 1, nages_dat);
      write_s_rownames_vec(sfile, ages_dat, 1, nages_dat, "ma", 0);  // this writes COLUMN names 

      write_s_matrix(sfile, value(eac_fish_unbiased_dat), "eac.f.unbiased", styr_fish,endyr_r, 1, nages_dat);
      write_s_rownames_vec(sfile, ages_dat, 1, nages_dat, "ma", 0);  // this writes COLUMN names
    } 
    
// ai survey ages
    if(nyrs_surv3_unbiased_ac_r>0)
    {
      write_s_imatrix(sfile, value(oac_surv3_unbiased_r), yrs_surv3_unbiased_ac_r, "oac.srv", 1,nyrs_surv3_unbiased_ac_r, 1, nages_dat);
      write_s_rownames_vec(sfile, ages_dat, 1, nages_dat, "ma", 0);  // this writes COLUMN names

      write_s_matrix(sfile, value(eac_surv3_unbiased_dat), "eac.srv", styr ,endyr_r, 1, nages_dat);
      write_s_rownames_vec(sfile, ages_dat, 1, nages_dat, "ma", 0);  // this writes COLUMN names
    }

//  fishery lengths
    if(nyrs_fish_lc_r>0)
    {
      write_s_imatrix(sfile, value(olc_fish_r), yrs_fish_lc, "olc.fish", 1,nyrs_fish_lc_r, 1, nlen);
      write_s_rownames_vec(sfile, lengths, 1, nlen, "ma", 0);  // this writes COLUMN names

      write_s_matrix(sfile, value(elc_fish), "elc.fish", styr_fish, endyr_r, 1, nlen);
      write_s_rownames_vec(sfile, lengths, 1, nlen, "ma", 0);  // this writes COLUMN names
    }

//  ai survey lengths
    if(nyrs_surv3_lc_r>0)
    {
      write_s_imatrix(sfile, value(olc_surv3_r), yrs_surv3_lc_r, "olc.srv", 1,nyrs_surv3_lc_r, 1, nlen);
      write_s_rownames_vec(sfile, lengths, 1, nlen, "ma", 0);  // this writes COLUMN names

      write_s_matrix(sfile, value(elc_surv3), "elc.srv", styr, endyr_r, 1, nlen);
      write_s_rownames_vec(sfile, lengths, 1, nlen, "ma", 0);  // this writes COLUMN names
    }

// ebs survey ages
    if(nyrs_surv_ebs_unbiased_ac_r>0)
    {
      write_s_imatrix(sfile, value(oac_surv_ebs_unbiased_r), yrs_surv_ebs_unbiased_ac_r, "oac.srv_ebs", 1,nyrs_surv_ebs_unbiased_ac_r, 1, nages_dat);
      write_s_rownames_vec(sfile, ages_dat, 1, nages_dat, "ma", 0);  // this writes COLUMN names

      write_s_matrix(sfile, value(eac_surv_ebs_unbiased_dat), "eac.srv_ebs", styr ,endyr_r, 1, nages_dat);
      write_s_rownames_vec(sfile, ages_dat, 1, nages_dat, "ma", 0);  // this writes COLUMN names
    }

    //  ebs survey lengths
    if(nyrs_surv_ebs_lc_r>0)
    {
      write_s_imatrix(sfile, value(olc_surv_ebs_r), yrs_surv_ebs_lc_r, "olc.srv_ebs", 1,nyrs_surv_ebs_lc_r, 1, nlen);
      write_s_rownames_vec(sfile, lengths, 1, nlen, "ma", 0);  // this writes COLUMN names

      write_s_matrix(sfile, value(elc_surv_ebs), "elc.srv_ebs", styr, endyr_r, 1, nlen);
      write_s_rownames_vec(sfile, lengths, 1, nlen, "ma", 0);  // this writes COLUMN names
    }  


//-----------EFFECTIVE N BY YEAR ---------------------------------------------------------------
//	sfile << "effn = structure(list(" << endl;
//	write_s_yrcol(sfile, value(effn(1)), "ac.f.biased", 1, 40, 0);
//	write_s_yrcol(sfile, value(effn(1)), "ac.f.unbiased", 1, 40, 0);
// 	write_s_yrcol(sfile, value(effn(2)), "lc.f", 1, 40, 0);
// 	write_s_yrcol(sfile, value(effn(3)), "ac.s.", 1, 40, 0);
// 	write_s_yrcol(sfile, value(effn(4)), "lc.s", 1, 40, 0);
//        write_s_yrcol(sfile, value(effn(5)), "ac.s_ebs", 1, 40, 1);
//	write_s_rownames_int(sfile, 1,40, "df", 0);

  // ------------------ START DATA LIKELIHOOD COMPONENTS VECTOR -------------------------
  // The vector of likelihood components will be stored as an S numeric vector.

  // Write start of likelihood components vector:
      sfile << "datalikecomp = structure(c(";

  //Write likelihood component values:
      sfile <<surv_like(1) << cc << surv_like(2) << cc;  
      sfile << catch_like << cc << age_like(1) << cc;
      sfile << age_like(2) << cc << age_like(3) << cc;
      sfile << age_like(4) << cc << age_like(5) << cc << age_like(6) << cc << mat_like;
      sfile << ")," << endl; 

  //Write names of likelihood components:
      sfile << ".Names = c('aisurvlike','ebssurvlike',";
      sfile << "'catch.like',";
      sfile << "'fish.unbiased.ac','fish.lc','aisrv.ac','aisrv.lc','ebssrv.ac','ebssrv.lc','mat_like'"<< endl;
      sfile << "))," << endl;   // Close vector of parameter values:

  // ------------------ START PENALTY LIKELIHOOD COMPONENTS VECTOR -------------------------
  // The vector of likelihood components will be stored as an S numeric vector.
  // This only contains the likelihood for the priors and penalties (not the data)

  // Write start of likelihood components vector:
      sfile << "pen_likecomp = structure(c(";

  //Write likelihood component values:
      sfile << hf_pen  << cc;
      sfile << sum(rec_like) << cc <<priormq(1)<< cc << priormq(2) << cc;  
      sfile << fpen << cc << sel_like(1) << cc;
      sfile << sel_like(2) << cc << sel_like(3) << cc << sel_like(4) << cc;
      sfile << sel_like(5) << cc << sel_like(6) << cc << sel_like(7) << cc;
      sfile << sel_like(8) << cc << sel_like(9);
      sfile << ")," << endl; 

  //Write names of likelihood components:
      sfile << ".Names = c('histFpen','reclike','prior_q','prior_m',";
     sfile << "'Fmortpen',";
      sfile << "'sel_like(1)','sel_like(2)','sel_like(3)', 'sel_like(4)',"<< endl;
      sfile << "'sel_like(5)','sel_like(6)','sel_like(7)', 'sel_like(8)','sel_like(9)'"<< endl;
      sfile << "))," << endl;   // Close vector of parameter values:

// ------------------ START MEAN SAMPLE WEIGHTS VECTOR -------------------------
  // The vector of mean sample weights for the age and length comps will be stored as an S numeric vector.

  // Write start of vector:
      sfile << "avgsampwts = structure(c(";

  //Write the values:
      if(nyrs_fish_unbiased_ac_r>0)  sfile << sum(fish_unbiased_ac_samp_r)/nyrs_fish_unbiased_ac_r << cc;
      else  sfile << 0 << cc; 
      
      if(nyrs_fish_lc_r>0) sfile << sum(fish_lc_samp_r)/nyrs_fish_lc_r << cc;
      else  sfile << 0 << cc; 

      if(nyrs_surv3_unbiased_ac_r>0) sfile << sum(surv3_unbiased_ac_samp_r)/nyrs_surv3_unbiased_ac_r << cc;
      else  sfile << 0 << cc; 

      if(nyrs_surv3_lc_r>0) sfile << sum(surv3_lc_samp_r)/nyrs_surv3_lc_r << cc;
      else  sfile << 0 << cc; 

      if(nyrs_surv_ebs_unbiased_ac_r>0) sfile << sum(surv_ebs_unbiased_ac_samp_r)/nyrs_surv_ebs_unbiased_ac_r << cc;
      else  sfile << 0 << cc; 

      if(nyrs_surv_ebs_lc_r>0) sfile << sum(surv_ebs_lc_samp_r)/nyrs_surv_ebs_lc_r;
      else  sfile << 0; 

      sfile << ")," << endl; 

  //Write names:
      sfile << ".Names = c('samp_ac.fish.unbiased','samp_lc.f','samp_ac.s','samp_lc.s','samp_ac.s_ebs','samp_lc.s_ebs'"<< endl;
      sfile << "))," << endl;   // Close vector of parameter values:

  // ------------------ START MEAN EFFECTIVE N VECTOR -------------------------
  // The vector of mean sample weights for the age and length comps will be stored as an S numeric vector.

  // Write start of vector:
      sfile << "meaneffn = structure(c(";

  //Write the values:
      if(nyrs_fish_unbiased_ac_r>0) sfile << (sum(effn(1)))/nyrs_fish_unbiased_ac_r << cc;
      else  sfile << 0 << cc; 

      if(nyrs_fish_lc_r>0) sfile << (sum(effn(2)))/nyrs_fish_lc_r << cc;
      else  sfile << 0 << cc; 

      if(nyrs_surv3_unbiased_ac_r>0) sfile << (sum(effn(3)))/nyrs_surv3_unbiased_ac_r << cc;
      else  sfile << 0 << cc; 

      if(nyrs_surv3_lc_r>0) sfile << (sum(effn(4)))/nyrs_surv3_lc_r << cc;
      else  sfile << 0 << cc; 

      if(nyrs_surv_ebs_unbiased_ac_r>0) sfile << (sum(effn(5)))/nyrs_surv_ebs_unbiased_ac_r << cc;
      else  sfile << 0 << cc; 
      
      if(nyrs_surv_ebs_lc_r>0) sfile << (sum(effn(6)))/nyrs_surv_ebs_lc_r;
      else  sfile << 0;   

      sfile << ")," << endl; 

  //Write names:
      sfile << ".Names = c('effn_ac.fish.unbiased','effn_lc.f','effn_ac.s','effn_lc.s','effn_ac.s_ebs','effn_lc.s_ebs'"<< endl;
//    sfile << ".Names = c('ac.fish.unbiased','lc.f','ac.s'"<< endl;
      sfile << "))," << endl;   // Close vector of parameter values:  

// ------------------ START RMSE VECTOR -------------------------
  // The vector of root mean square error for the reccruitment and survey data will be stored as an S numeric vector.

  // Write start of vector:
      sfile << "rmse = structure(c(";

  //Write the values:
      sfile << rmse(1) << cc;
      sfile << rmse(2) << cc;
      sfile << rmse(3) << cc;
      sfile << al_rmse(1) << cc;
      sfile << al_rmse(2) << cc;
      sfile << al_rmse(3) << cc;
      sfile << al_rmse(4) << cc;
      sfile << al_rmse(5) << cc;
      sfile << al_rmse(6) ;

      sfile << ")," << endl; 

  //Write names:
      sfile << ".Names = c('rmse_aisrv','rmse_ebssrv','rmse_rec','rmse_fac','rmse_flc','rmse_sac','rmse_slc','rmse_ac.s_ebs','rmse_lc.s_ebs'"<< endl;
      sfile << "))," << endl;   // Close vector of parameter values:

// ------------------ START SDNR VECTOR -------------------------
  // The standard deviation of normalized residuals of the survey and the age and length comps 

  // Write start of vector:
      sfile << "sdnr = structure(c(";

  //Write the values:
      sfile << sdnr(7) << cc;
      sfile << sdnr(1) << cc;
      sfile << sdnr(2) << cc;
      sfile << sdnr(3) << cc;
      sfile << sdnr(4) << cc;
      sfile << sdnr(5) << cc;
      sfile << sdnr(6);

      sfile << ")," << endl; 

  //Write names:
      sfile << ".Names = c('sdnr_aisrv','sdnr_ac.fish.unbiased','sdnr_lc.f','sdnr_ac.s','sdnr_lc.s','sdnr_ac.s_ebs','sdnr_lc.s_ebs'"<< endl;
      sfile << "))," << endl;   // Close vector of parameter values:


// ------------------ START SELECTIVITY BINS -------------------------
  // The number of time bins for the fishery selectivity

  // Write start of vector:
      sfile << "fishselbins = structure(c(";

  //Write the values:
      sfile << nbins;
      sfile << ")," << endl;

  //Write names:
      sfile << ".Names = c('fishselbins'"<< endl; 
      sfile << "))," << endl;   // Close 

  // ------------------ START  CONTROL RULE VECTOR -------------------------
  // The vector of info needed for mapping the control rule will be stored as an S numeric vector.

  // Write start of likelihood components vector:
      sfile << "controlrule = structure(c(";

  //Write likelihood component values:
      sfile << F40 << cc;
      sfile << F35 << cc;
      sfile << SBF40 << cc;
      sfile << SBF35 << cc;
      sfile << max(1977,styr)+rec_age << cc;
      sfile << lastyr_rec << cc;
      sfile << lastyr_rec_a10; 
      sfile << ")," << endl; 

  //Write names:
      sfile << ".Names = c('F40','F35','SBF40','SBF35','styr_rec','endyr_rec','endyr_rec_a10'"<< endl;
      sfile << "))" << endl;   // Close vector of control rule values:

  //-------- Write close of overall S list structure:--------------------------------------------
  // This is the closing punctuation plus the names of the constituent elements of the list

    sfile << "))" << endl;   // Close 
  
  //sfile << "), .Names = c('t.series'";
  //sfile << "))" << endl;

  //-------- END OF CXX FILE --------------------------------------------
