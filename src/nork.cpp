#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
 # include "admodel.h"                      // Include AD class definitions
  adstring model_name;
  adstring data_file;
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <nork.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  pad_evalout = new ofstream("evalout.prj");;
 ad_comm::change_datafile_name("goa_nr_2018.ctl");    //Read in phases, penalties and priors from "tem.ctl"
 *(ad_comm::global_datafile) >>  model_name; 
 *(ad_comm::global_datafile) >>  data_file; 
  SrType.allocate("SrType");
  styr_rec_est.allocate("styr_rec_est");
  endyr_rec_est.allocate("endyr_rec_est");
 nrecs_est = endyr_rec_est-styr_rec_est+1;
  rec_like_type.allocate("rec_like_type");
  fsh_sel_opt.allocate("fsh_sel_opt");
  srv1_sel_opt.allocate("srv1_sel_opt");
  srv2_sel_opt.allocate("srv2_sel_opt");
  n_fish_sel_ages.allocate("n_fish_sel_ages");
  n_srv1_sel_ages.allocate("n_srv1_sel_ages");
  n_srv2_sel_ages.allocate("n_srv2_sel_ages");
  ph_Fdev.allocate("ph_Fdev");
  ph_avg_F.allocate("ph_avg_F");
  ph_recdev.allocate("ph_recdev");
  ph_fish_sel.allocate("ph_fish_sel");
  ph_srv1_sel.allocate("ph_srv1_sel");
  ph_srv2_sel.allocate("ph_srv2_sel");
  mprior.allocate("mprior");
 log_mprior = log(mprior);
  cvmprior.allocate("cvmprior");
  ph_m.allocate("ph_m");
  steep_prior.allocate("steep_prior");
  cv_steep_prior.allocate("cv_steep_prior");
  ph_steepness.allocate("ph_steepness");
 if (SrType==3) ph_Rzero=-6; else ph_Rzero=2;
  sigrprior.allocate("sigrprior");
  cvsigrprior.allocate("cvsigrprior");
  ph_sigr.allocate("ph_sigr");
  q_srv1prior.allocate("q_srv1prior");
 log_q_srv1prior=log(q_srv1prior);
  cvq_srv1prior.allocate("cvq_srv1prior");
  ph_q_srv1.allocate("ph_q_srv1");
  q_srv2prior.allocate("q_srv2prior");
 log_q_srv2prior=log(q_srv2prior);
  cvq_srv2prior.allocate("cvq_srv2prior");
  ph_q_srv2.allocate("ph_q_srv2");
  ph_ESS.allocate("ph_ESS");
  yrs_ESS.allocate("yrs_ESS");
  pen_ESS.allocate("pen_ESS");
  C_like_type.allocate("C_like_type");
  S_like_type.allocate("S_like_type");
  yr_catchwt.allocate("yr_catchwt");
  wt_ssqcatch.allocate("wt_ssqcatch");
  wt_ssqcatch2.allocate("wt_ssqcatch2");
  wt_cpue.allocate("wt_cpue");
  wt_srv1.allocate("wt_srv1");
  wt_srv2.allocate("wt_srv2");
  wt_fish_age.allocate("wt_fish_age");
  wt_srv1_age.allocate("wt_srv1_age");
  wt_fish_size.allocate("wt_fish_size");
  wt_srv1_size.allocate("wt_srv1_size");
  wt_srv2_size.allocate("wt_srv2_size");
  wt_rec_var.allocate("wt_rec_var");
  wt_sel_reg_fish.allocate("wt_sel_reg_fish");
  wt_sel_reg_srv1.allocate("wt_sel_reg_srv1");
  wt_sel_reg_srv2.allocate("wt_sel_reg_srv2");
  wt_sel_dome_fish.allocate("wt_sel_dome_fish");
  wt_sel_dome_srv1.allocate("wt_sel_dome_srv1");
  wt_sel_dome_srv2.allocate("wt_sel_dome_srv2");
  wt_fmort_reg.allocate("wt_fmort_reg");
  wt_avg_sel.allocate("wt_avg_sel");
  ph_logsurv.allocate("ph_logsurv");
  initial_LMR.allocate("initial_LMR");
  wt_Rzero.allocate("wt_Rzero");
  agesamplestyle.allocate("agesamplestyle");
  lensamplestyle.allocate("lensamplestyle");
  yieldratio.allocate("yieldratio");
 ad_comm::change_datafile_name(data_file);    // Read data from the data file
  styr.allocate("styr");
  endyr.allocate("endyr");
  recage.allocate("recage");
  nages_D.allocate("nages_D");
  nages_M.allocate("nages_M");
  nlenbins.allocate("nlenbins");
  n_ageage_mat.allocate("n_ageage_mat");
  n_sizeage_mat.allocate("n_sizeage_mat");
  len_bin_labels.allocate(1,nlenbins,"len_bin_labels");
  nyrs = endyr - styr + 1;
 styr_rec = (styr - nages_M) + 1;     // First year of recruitment
 styr_sp  = styr_rec - recage ;     // First year of spawning biomass  
 endyr_sp = endyr   - recage - 1;   // endyr year of (main) spawning biomass
  yy.allocate(styr,endyr);
 yy.fill_seqadd(styr,1) ;
  aa.allocate(1,nages_M);
 aa.fill_seqadd(recage,1) ;
 ph_F50 = 4;
  spawn_fract.allocate("spawn_fract");
 spawn_fract = (spawn_fract - 1) / 12;
  wt.allocate(1,nages_M,"wt");
  wt_mature.allocate(1,nages_M);
  obs_catch_early.allocate(styr,yr_catchwt,"obs_catch_early");
  obs_catch_later.allocate(yr_catchwt+1,endyr,"obs_catch_later");
  nyrs_cpue.allocate("nyrs_cpue");
  yrs_cpue.allocate(1,nyrs_cpue,"yrs_cpue");
  obs_cpue.allocate(1,nyrs_cpue,"obs_cpue");
 if (nyrs_cpue>0) mean_obs_cpue = exp(mean(log(obs_cpue))); 
  nyrs_srv1.allocate("nyrs_srv1");
  yrs_srv1.allocate(1,nyrs_srv1,"yrs_srv1");
  obs_srv1_biom.allocate(1,nyrs_srv1,"obs_srv1_biom");
  obs_srv1_se.allocate(1,nyrs_srv1,"obs_srv1_se");
  obs_srv1_lci.allocate(1,nyrs_srv1,"obs_srv1_lci");
  obs_srv1_uci.allocate(1,nyrs_srv1,"obs_srv1_uci");
  nyrs_srv2.allocate("nyrs_srv2");
  yrs_srv2.allocate(1,nyrs_srv2,"yrs_srv2");
  obs_srv2_biom.allocate(1,nyrs_srv2,"obs_srv2_biom");
  obs_srv2_se.allocate(1,nyrs_srv2,"obs_srv2_se");
  obs_srv2_lci.allocate(1,nyrs_srv2,"obs_srv2_lci");
  obs_srv2_uci.allocate(1,nyrs_srv2,"obs_srv2_uci");
  nyrs_fish_age.allocate("nyrs_fish_age");
  yrs_fish_age.allocate(1,nyrs_fish_age,"yrs_fish_age");
  nsamples_fish_age.allocate(1,nyrs_fish_age,"nsamples_fish_age");
  nhauls_fish_age.allocate(1,nyrs_fish_age,"nhauls_fish_age");
  age_age_ind_fsh.allocate(1,nyrs_fish_age,"age_age_ind_fsh");
  oac_fish.allocate(1,nyrs_fish_age,1,nages_D,"oac_fish");
  nmulti_fish_age.allocate(1,nyrs_fish_age);
  nyrs_srv1_age.allocate("nyrs_srv1_age");
  yrs_srv1_age.allocate(1,nyrs_srv1_age,"yrs_srv1_age");
  nsamples_srv1_age.allocate(1,nyrs_srv1_age,"nsamples_srv1_age");
  nhauls_srv1_age.allocate(1,nyrs_srv1_age,"nhauls_srv1_age");
  age_age_ind_srv.allocate(1,nyrs_srv1_age,"age_age_ind_srv");
  oac_srv1.allocate(1,nyrs_srv1_age,1,nages_D,"oac_srv1");
  nmulti_srv1_age.allocate(1,nyrs_srv1_age);
  nyrs_fish_size.allocate("nyrs_fish_size");
  yrs_fish_size.allocate(1,nyrs_fish_size,"yrs_fish_size");
  nsamples_fish_size.allocate(1,nyrs_fish_size,"nsamples_fish_size");
  nhauls_fish_size.allocate(1,nyrs_fish_size,"nhauls_fish_size");
  siz_age_ind_fsh.allocate(1,nyrs_fish_size,"siz_age_ind_fsh");
  osc_fish.allocate(1,nyrs_fish_size,1,nlenbins,"osc_fish");
  nmulti_fish_size.allocate(1,nyrs_fish_size);
  nyrs_srv1_size.allocate("nyrs_srv1_size");
  yrs_srv1_size.allocate(1,nyrs_srv1_size,"yrs_srv1_size");
  nsamples_srv1_size.allocate(1,nyrs_srv1_size,"nsamples_srv1_size");
  nhauls_srv1_size.allocate(1,nyrs_srv1_size,"nhauls_srv1_size");
  siz_age_ind_srv1.allocate(1,nyrs_srv1_size,"siz_age_ind_srv1");
  osc_srv1.allocate(1,nyrs_srv1_size,1,nlenbins,"osc_srv1");
  nmulti_srv1_size.allocate(1,nyrs_srv1_size);
  nyrs_srv2_size.allocate("nyrs_srv2_size");
  yrs_srv2_size.allocate(1,nyrs_srv2_size,"yrs_srv2_size");
  nsamples_srv2_size.allocate(1,nyrs_srv2_size,"nsamples_srv2_size");
  nhauls_srv2_size.allocate(1,nyrs_srv2_size,"nhauls_srv2_size");
  siz_age_ind_srv2.allocate(1,nyrs_srv2_size,"siz_age_ind_srv2");
  osc_srv2.allocate(1,nyrs_srv2_size,1,nlenbins,"osc_srv2");
  nmulti_srv2_size.allocate(1,nyrs_srv2_size);
  nsamples_fish_age_ts.allocate(styr,endyr);
  nsamples_srv1_age_ts.allocate(styr,endyr);
  nsamples_fish_size_ts.allocate(styr,endyr);
  nsamples_srv1_size_ts.allocate(styr,endyr);
  nsamples_srv2_size_ts.allocate(styr,endyr);
  sizeage.allocate(1,n_sizeage_mat,1,nages_M,1,nlenbins,"sizeage");
  ageage.allocate(1,n_ageage_mat,1,nages_M,1,nages_D,"ageage");
  offset.allocate(1,6);
  eof.allocate("eof");
 ad_comm::change_datafile_name("mat.dat");    // Read data from the data file
  nages_mat.allocate("nages_mat");
  L_tot_na.allocate(1,nages_mat,"L_tot_na");
  L_mat_na.allocate(1,nages_mat,"L_mat_na");
  L_wt.allocate(1,nages_mat,"L_wt");
  C_tot_na.allocate(1,nages_mat,"C_tot_na");
  C_mat_na.allocate(1,nages_mat,"C_mat_na");
  C_wt.allocate(1,nages_mat,"C_wt");
  offset.initialize();
   if(eof==42) cout<<"The data has been read correctly!";
   else { cout <<"You f'ed up your data file!"<<endl;exit(1); }
   if(wt_rec_var==0)
   {
     if (ph_sigr>0)
     {
       cout << "Warning, wt_rec_var is zero, so can't estimate sigr!@"<<endl;
       cout << "turning sigr off "<<endl;
       ph_sigr =-4;
       cout << "hit any key, then enter to continue"<<endl;
       char  xxx; cin >> xxx;
     }
   }
       phase_selcoff_fsh = -1; 
       phase_logist_fsh  = ph_fish_sel;
       phase_selcoff_srv1 = -1; 
       phase_logist_srv1  = ph_srv1_sel;
       phase_selcoff_srv2 = -1; 
       phase_logist_srv2  = ph_srv2_sel;
       nmulti_fish_age = sqrt(elem_prod(nhauls_fish_age,nsamples_fish_age));
       nmulti_fish_age =  nmulti_fish_age/max( nmulti_fish_age)*100;
       nmulti_srv1_age =sqrt(elem_prod(nhauls_srv1_age,nsamples_srv1_age));
       nmulti_srv1_age = nmulti_srv1_age/max(nmulti_srv1_age)*100;
       nmulti_fish_size = sqrt(elem_prod(nhauls_fish_size,nsamples_fish_size));
       nmulti_fish_size =  nmulti_fish_size/max( nmulti_fish_size)*100;
       nmulti_srv1_size =sqrt(elem_prod(nhauls_srv1_size,nsamples_srv1_size));
       nmulti_srv1_size = nmulti_srv1_size/max(nmulti_srv1_size)*100;
       nmulti_srv2_size =sqrt(elem_prod(nhauls_srv2_size,nsamples_srv2_size));
       nmulti_srv2_size = nmulti_srv2_size/max(nmulti_srv2_size)*100;
  for (i=1; i<=nyrs_fish_age; i++)
  {
   oac_fish(i)/=sum(oac_fish(i));
   offset(1) -= nmulti_fish_age(i) *((oac_fish(i) + 0.00001)*log(oac_fish(i) + 0.00001)); 
  }
  for (i=1; i<=nyrs_srv1_age; i++)
  {
   oac_srv1(i)/=sum(oac_srv1(i));
   offset(2) -= nmulti_srv1_age(i)*((oac_srv1(i) + 0.00001)*log(oac_srv1(i) + 0.00001));
  }
  for (i=1; i<=nyrs_fish_size; i++)
  {
   osc_fish(i)/=sum(osc_fish(i));
   offset(3) -= nmulti_fish_size(i)*((osc_fish(i) + 0.00001)*log(osc_fish(i) + 0.00001));
  }
  for (i=1; i<=nyrs_srv1_size; i++)
  {
   osc_srv1(i)/=sum(osc_srv1(i));
   offset(4) -= nmulti_srv1_size(i)*((osc_srv1(i) + 0.00001)*log(osc_srv1(i) + 0.00001));
  }
  for (i=1; i<=nyrs_srv2_size; i++)
  {
   osc_srv2(i)/=sum(osc_srv2(i));
   offset(5) -= nmulti_srv2_size(i)*((osc_srv2(i) + 0.00001)*log(osc_srv2(i) + 0.00001));
  }
}

void model_parameters::initializationfunction(void)
{
  logm.set_initial_value(log_mprior);
  log_mean_rec.set_initial_value(initial_LMR);
  log_Rzero.set_initial_value(initial_LMR);
  sigr.set_initial_value(sigrprior);
  a50.set_initial_value(7.5);
  delta.set_initial_value(3.);
  a50_srv1.set_initial_value(7.3);
  delta_srv1.set_initial_value(3.8);
  a50_srv2.set_initial_value(7.3);
  delta_srv2.set_initial_value(3.8);
  Est_nsamples_fish_age.set_initial_value(50);
  Est_nsamples_srv1_age.set_initial_value(50);
  Est_nsamples_fish_size.set_initial_value(50);
  mat_a50.set_initial_value(9);
  mat_delta.set_initial_value(1);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  steepness.allocate(0.2001,0.999,ph_steepness,"steepness");
  log_Rzero.allocate(ph_Rzero,"log_Rzero");
  sam_rec.allocate(styr_rec,endyr,"sam_rec");
  #ifndef NO_AD_INITIALIZE
    sam_rec.initialize();
  #endif
  srm_rec.allocate(styr_rec,endyr,"srm_rec");
  #ifndef NO_AD_INITIALIZE
    srm_rec.initialize();
  #endif
  Sp_Biom.allocate(styr_sp,endyr,"Sp_Biom");
  #ifndef NO_AD_INITIALIZE
    Sp_Biom.initialize();
  #endif
  log_mean_rec.allocate(1,"log_mean_rec");
  sigr.allocate(0.3,10,ph_sigr,"sigr");
  sigrsq.allocate("sigrsq");
  #ifndef NO_AD_INITIALIZE
  sigrsq.initialize();
  #endif
  alpha.allocate("alpha");
  #ifndef NO_AD_INITIALIZE
  alpha.initialize();
  #endif
  beta.allocate("beta");
  #ifndef NO_AD_INITIALIZE
  beta.initialize();
  #endif
  Bzero.allocate("Bzero");
  #ifndef NO_AD_INITIALIZE
  Bzero.initialize();
  #endif
  Rzero.allocate("Rzero");
  #ifndef NO_AD_INITIALIZE
  Rzero.initialize();
  #endif
  phizero.allocate("phizero");
  #ifndef NO_AD_INITIALIZE
  phizero.initialize();
  #endif
  log_fish_sel_coffs.allocate(1,n_fish_sel_ages,phase_selcoff_fsh,"log_fish_sel_coffs");
  a50.allocate(phase_logist_fsh,"a50");
  delta.allocate(phase_logist_fsh,"delta");
  log_fish_sel.allocate(1,nages_M,"log_fish_sel");
  #ifndef NO_AD_INITIALIZE
    log_fish_sel.initialize();
  #endif
  fish_sel.allocate(1,nages_M,"fish_sel");
  #ifndef NO_AD_INITIALIZE
    fish_sel.initialize();
  #endif
  log_avgfishsel.allocate("log_avgfishsel");
  #ifndef NO_AD_INITIALIZE
  log_avgfishsel.initialize();
  #endif
  log_srv1_sel_coffs.allocate(1,n_srv1_sel_ages,phase_selcoff_srv1,"log_srv1_sel_coffs");
  a50_srv1.allocate(phase_logist_srv1,"a50_srv1");
  delta_srv1.allocate(phase_logist_srv1,"delta_srv1");
  log_srv1_sel.allocate(1,nages_M,"log_srv1_sel");
  #ifndef NO_AD_INITIALIZE
    log_srv1_sel.initialize();
  #endif
  srv1_sel.allocate(1,nages_M,"srv1_sel");
  #ifndef NO_AD_INITIALIZE
    srv1_sel.initialize();
  #endif
  log_avgsrv1sel.allocate("log_avgsrv1sel");
  #ifndef NO_AD_INITIALIZE
  log_avgsrv1sel.initialize();
  #endif
  log_srv2_sel_coffs.allocate(1,n_srv2_sel_ages,phase_selcoff_srv2,"log_srv2_sel_coffs");
  a50_srv2.allocate(phase_logist_srv2,"a50_srv2");
  delta_srv2.allocate(phase_logist_srv2,"delta_srv2");
  log_srv2_sel.allocate(1,nages_M,"log_srv2_sel");
  #ifndef NO_AD_INITIALIZE
    log_srv2_sel.initialize();
  #endif
  srv2_sel.allocate(1,nages_M,"srv2_sel");
  #ifndef NO_AD_INITIALIZE
    srv2_sel.initialize();
  #endif
  log_avgsrv2sel.allocate("log_avgsrv2sel");
  #ifndef NO_AD_INITIALIZE
  log_avgsrv2sel.initialize();
  #endif
  log_avg_F.allocate(ph_avg_F,"log_avg_F");
  log_F_devs.allocate(styr,endyr,-15.,15.,ph_Fdev,"log_F_devs");
  Fmort.allocate(styr,endyr,"Fmort");
  #ifndef NO_AD_INITIALIZE
    Fmort.initialize();
  #endif
  Z.allocate(styr,endyr,1,nages_M,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  F.allocate(styr,endyr,1,nages_M,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  S.allocate(styr,endyr,1,nages_M,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  Est_nsamples_fish_age.allocate(10,250,ph_ESS,"Est_nsamples_fish_age");
  Est_nsamples_srv1_age.allocate(10,250,ph_ESS,"Est_nsamples_srv1_age");
  Est_nsamples_fish_size.allocate(10,250,ph_ESS,"Est_nsamples_fish_size");
  mat_a50.allocate(2,"mat_a50");
  mat_delta.allocate(2,"mat_delta");
  L_pmat.allocate(1,nages_mat,"L_pmat");
  #ifndef NO_AD_INITIALIZE
    L_pmat.initialize();
  #endif
  C_pmat.allocate(1,nages_mat,"C_pmat");
  #ifndef NO_AD_INITIALIZE
    C_pmat.initialize();
  #endif
  Pred_pmat.allocate(1,nages_mat,"Pred_pmat");
  #ifndef NO_AD_INITIALIZE
    Pred_pmat.initialize();
  #endif
  p_mature.allocate(1,nages_M,"p_mature");
  #ifndef NO_AD_INITIALIZE
    p_mature.initialize();
  #endif
  wt_mature.allocate(1,nages_M,"wt_mature");
  #ifndef NO_AD_INITIALIZE
    wt_mature.initialize();
  #endif
  Like_L_vec.allocate(1,nages_mat,"Like_L_vec");
  #ifndef NO_AD_INITIALIZE
    Like_L_vec.initialize();
  #endif
  Like_C_vec.allocate(1,nages_mat,"Like_C_vec");
  #ifndef NO_AD_INITIALIZE
    Like_C_vec.initialize();
  #endif
  Like_L.allocate("Like_L");
  #ifndef NO_AD_INITIALIZE
  Like_L.initialize();
  #endif
  Like_C.allocate("Like_C");
  #ifndef NO_AD_INITIALIZE
  Like_C.initialize();
  #endif
  k.allocate("k");
  #ifndef NO_AD_INITIALIZE
  k.initialize();
  #endif
 k=0.00001;
  natmortv.allocate(1,nages_M,"natmortv");
  #ifndef NO_AD_INITIALIZE
    natmortv.initialize();
  #endif
  log_rec_dev.allocate(styr_rec+1,endyr,-10.,10.,ph_recdev,"log_rec_dev");
  natage.allocate(styr,endyr,1,nages_M,"natage");
  #ifndef NO_AD_INITIALIZE
    natage.initialize();
  #endif
  batage.allocate(styr,endyr,1,nages_M,"batage");
  #ifndef NO_AD_INITIALIZE
    batage.initialize();
  #endif
  catage.allocate(styr,endyr,1,nages_M,"catage");
  #ifndef NO_AD_INITIALIZE
    catage.initialize();
  #endif
  pred_catch_early.allocate(styr,yr_catchwt,"pred_catch_early");
  #ifndef NO_AD_INITIALIZE
    pred_catch_early.initialize();
  #endif
  pred_catch_later.allocate(yr_catchwt+1,endyr,"pred_catch_later");
  #ifndef NO_AD_INITIALIZE
    pred_catch_later.initialize();
  #endif
  log_q_srv1.allocate(ph_q_srv1,"log_q_srv1");
  log_q_srv2.allocate(ph_q_srv2,"log_q_srv2");
  cv_cpue.allocate(-1,"cv_cpue");
  logm.allocate(ph_m,"logm");
  q_cpue.allocate("q_cpue");
  #ifndef NO_AD_INITIALIZE
  q_cpue.initialize();
  #endif
  pred_cpue.allocate(1,nyrs_cpue,"pred_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_cpue.initialize();
  #endif
  pred_srv1.allocate(1,nyrs_srv1,"pred_srv1");
  #ifndef NO_AD_INITIALIZE
    pred_srv1.initialize();
  #endif
  pred_srv2.allocate(1,nyrs_srv2,"pred_srv2");
  #ifndef NO_AD_INITIALIZE
    pred_srv2.initialize();
  #endif
  eac_fish.allocate(1,nyrs_fish_age,1,nages_D,"eac_fish");
  #ifndef NO_AD_INITIALIZE
    eac_fish.initialize();
  #endif
  eac_srv1.allocate(1,nyrs_srv1_age,1,nages_D,"eac_srv1");
  #ifndef NO_AD_INITIALIZE
    eac_srv1.initialize();
  #endif
  esc_fish.allocate(1,nyrs_fish_size,1,nlenbins,"esc_fish");
  #ifndef NO_AD_INITIALIZE
    esc_fish.initialize();
  #endif
  esc_srv1.allocate(1,nyrs_srv1_size,1,nlenbins,"esc_srv1");
  #ifndef NO_AD_INITIALIZE
    esc_srv1.initialize();
  #endif
  esc_srv2.allocate(1,nyrs_srv2_size,1,nlenbins,"esc_srv2");
  #ifndef NO_AD_INITIALIZE
    esc_srv2.initialize();
  #endif
  effn_fish_age.allocate(1,nyrs_fish_age,"effn_fish_age");
  #ifndef NO_AD_INITIALIZE
    effn_fish_age.initialize();
  #endif
  effn_fish_size.allocate(1,nyrs_fish_size,"effn_fish_size");
  #ifndef NO_AD_INITIALIZE
    effn_fish_size.initialize();
  #endif
  sdnr_fish_age.allocate(1,nyrs_fish_age,"sdnr_fish_age");
  #ifndef NO_AD_INITIALIZE
    sdnr_fish_age.initialize();
  #endif
  sdnr_fish_size.allocate(1,nyrs_fish_size,"sdnr_fish_size");
  #ifndef NO_AD_INITIALIZE
    sdnr_fish_size.initialize();
  #endif
  effn_srv1_age.allocate(1,nyrs_srv1_age,"effn_srv1_age");
  #ifndef NO_AD_INITIALIZE
    effn_srv1_age.initialize();
  #endif
  effn_srv1_size.allocate(1,nyrs_srv1_size,"effn_srv1_size");
  #ifndef NO_AD_INITIALIZE
    effn_srv1_size.initialize();
  #endif
  sdnr_srv1_age.allocate(1,nyrs_srv1_age,"sdnr_srv1_age");
  #ifndef NO_AD_INITIALIZE
    sdnr_srv1_age.initialize();
  #endif
  sdnr_srv1_size.allocate(1,nyrs_srv1_size,"sdnr_srv1_size");
  #ifndef NO_AD_INITIALIZE
    sdnr_srv1_size.initialize();
  #endif
  effn_srv2_size.allocate(1,nyrs_srv2_size,"effn_srv2_size");
  #ifndef NO_AD_INITIALIZE
    effn_srv2_size.initialize();
  #endif
  sdnr_srv2_size.allocate(1,nyrs_srv2_size,"sdnr_srv2_size");
  #ifndef NO_AD_INITIALIZE
    sdnr_srv2_size.initialize();
  #endif
  effn_fish_age_ts.allocate(styr,endyr,"effn_fish_age_ts");
  #ifndef NO_AD_INITIALIZE
    effn_fish_age_ts.initialize();
  #endif
  effn_fish_size_ts.allocate(styr,endyr,"effn_fish_size_ts");
  #ifndef NO_AD_INITIALIZE
    effn_fish_size_ts.initialize();
  #endif
  sdnr_fish_age_ts.allocate(styr,endyr,"sdnr_fish_age_ts");
  #ifndef NO_AD_INITIALIZE
    sdnr_fish_age_ts.initialize();
  #endif
  sdnr_fish_size_ts.allocate(styr,endyr,"sdnr_fish_size_ts");
  #ifndef NO_AD_INITIALIZE
    sdnr_fish_size_ts.initialize();
  #endif
  effn_srv1_age_ts.allocate(styr,endyr,"effn_srv1_age_ts");
  #ifndef NO_AD_INITIALIZE
    effn_srv1_age_ts.initialize();
  #endif
  effn_srv1_size_ts.allocate(styr,endyr,"effn_srv1_size_ts");
  #ifndef NO_AD_INITIALIZE
    effn_srv1_size_ts.initialize();
  #endif
  sdnr_srv1_age_ts.allocate(styr,endyr,"sdnr_srv1_age_ts");
  #ifndef NO_AD_INITIALIZE
    sdnr_srv1_age_ts.initialize();
  #endif
  sdnr_srv1_size_ts.allocate(styr,endyr,"sdnr_srv1_size_ts");
  #ifndef NO_AD_INITIALIZE
    sdnr_srv1_size_ts.initialize();
  #endif
  effn_srv2_size_ts.allocate(styr,endyr,"effn_srv2_size_ts");
  #ifndef NO_AD_INITIALIZE
    effn_srv2_size_ts.initialize();
  #endif
  sdnr_srv2_size_ts.allocate(styr,endyr,"sdnr_srv2_size_ts");
  #ifndef NO_AD_INITIALIZE
    sdnr_srv2_size_ts.initialize();
  #endif
  tot_biom.allocate(styr,endyr,"tot_biom");
  q_srv1.allocate("q_srv1");
  q_srv2.allocate("q_srv2");
  #ifndef NO_AD_INITIALIZE
  q_srv2.initialize();
  #endif
  pred_rec.allocate(styr,endyr,"pred_rec");
  expl_rate.allocate(styr,endyr,"expl_rate");
  avg_rec.allocate("avg_rec");
  spbiom_trend.allocate("spbiom_trend");
  Depletion.allocate("Depletion");
  spawn_biom.allocate(styr,endyr,"spawn_biom");
  natmort.allocate("natmort");
  #ifndef NO_AD_INITIALIZE
  natmort.initialize();
  #endif
  LMR.allocate("LMR");
  cigar.allocate("cigar");
  q2.allocate("q2");
  nattymort.allocate("nattymort");
  mF50.allocate(0.01,1.,ph_F50,"mF50");
  mF40.allocate(0.01,1.,ph_F50,"mF40");
  mF35.allocate(0.01,1.,ph_F50,"mF35");
  F50.allocate("F50");
  F40.allocate("F40");
  F35.allocate("F35");
  SB0.allocate("SB0");
  #ifndef NO_AD_INITIALIZE
  SB0.initialize();
  #endif
  SBF50.allocate("SBF50");
  #ifndef NO_AD_INITIALIZE
  SBF50.initialize();
  #endif
  SBF40.allocate("SBF40");
  #ifndef NO_AD_INITIALIZE
  SBF40.initialize();
  #endif
  SBF35.allocate("SBF35");
  #ifndef NO_AD_INITIALIZE
  SBF35.initialize();
  #endif
  sprpen.allocate("sprpen");
  #ifndef NO_AD_INITIALIZE
  sprpen.initialize();
  #endif
  Nspr.allocate(1,4,1,nages_M,"Nspr");
  #ifndef NO_AD_INITIALIZE
    Nspr.initialize();
  #endif
  surv_like.allocate(1,3,"surv_like");
  #ifndef NO_AD_INITIALIZE
    surv_like.initialize();
  #endif
  cpue_like.allocate("cpue_like");
  #ifndef NO_AD_INITIALIZE
  cpue_like.initialize();
  #endif
  age_like.allocate(1,6,"age_like");
  #ifndef NO_AD_INITIALIZE
    age_like.initialize();
  #endif
  sel_like.allocate(1,6,"sel_like");
  #ifndef NO_AD_INITIALIZE
    sel_like.initialize();
  #endif
  rec_like.allocate("rec_like");
  #ifndef NO_AD_INITIALIZE
  rec_like.initialize();
  #endif
  ssqcatch.allocate("ssqcatch");
  #ifndef NO_AD_INITIALIZE
  ssqcatch.initialize();
  #endif
  F_mort_regularity.allocate("F_mort_regularity");
  #ifndef NO_AD_INITIALIZE
  F_mort_regularity.initialize();
  #endif
  avg_sel_penalty.allocate("avg_sel_penalty");
  #ifndef NO_AD_INITIALIZE
  avg_sel_penalty.initialize();
  #endif
  dirich_fish_age.allocate(1,nyrs_fish_age,1,nages_D,"dirich_fish_age");
  #ifndef NO_AD_INITIALIZE
    dirich_fish_age.initialize();
  #endif
  dirich_srv1_age.allocate(1,nyrs_srv1_age,1,nages_D,"dirich_srv1_age");
  #ifndef NO_AD_INITIALIZE
    dirich_srv1_age.initialize();
  #endif
  dirich_fish_size.allocate(1,nyrs_fish_size,1,nlenbins,"dirich_fish_size");
  #ifndef NO_AD_INITIALIZE
    dirich_fish_size.initialize();
  #endif
  dir_fish_age_pen.allocate(1,nyrs_fish_age,1,nages_D,"dir_fish_age_pen");
  #ifndef NO_AD_INITIALIZE
    dir_fish_age_pen.initialize();
  #endif
  dir_fa_pen.allocate(1,nyrs_fish_age,"dir_fa_pen");
  #ifndef NO_AD_INITIALIZE
    dir_fa_pen.initialize();
  #endif
  dir_srv1_age_pen.allocate(1,nyrs_srv1_age,1,nages_D,"dir_srv1_age_pen");
  #ifndef NO_AD_INITIALIZE
    dir_srv1_age_pen.initialize();
  #endif
  dir_sa_pen.allocate(1,nyrs_srv1_age,"dir_sa_pen");
  #ifndef NO_AD_INITIALIZE
    dir_sa_pen.initialize();
  #endif
  dir_fish_size_pen.allocate(1,nyrs_fish_size,1,nlenbins,"dir_fish_size_pen");
  #ifndef NO_AD_INITIALIZE
    dir_fish_size_pen.initialize();
  #endif
  dir_fs_pen.allocate(1,nyrs_fish_size,"dir_fs_pen");
  #ifndef NO_AD_INITIALIZE
    dir_fs_pen.initialize();
  #endif
  digamma_fish_age_n.allocate(1,nyrs_fish_age,1,nages_D,"digamma_fish_age_n");
  #ifndef NO_AD_INITIALIZE
    digamma_fish_age_n.initialize();
  #endif
  digamma_fish_eac.allocate(1,nyrs_fish_age,1,nages_D,"digamma_fish_eac");
  #ifndef NO_AD_INITIALIZE
    digamma_fish_eac.initialize();
  #endif
  digamma_srv1_age_n.allocate(1,nyrs_srv1_age,1,nages_D,"digamma_srv1_age_n");
  #ifndef NO_AD_INITIALIZE
    digamma_srv1_age_n.initialize();
  #endif
  digamma_srv1_eac.allocate(1,nyrs_srv1_age,1,nages_D,"digamma_srv1_eac");
  #ifndef NO_AD_INITIALIZE
    digamma_srv1_eac.initialize();
  #endif
  digamma_fish_size_n.allocate(1,nyrs_fish_size,1,nlenbins,"digamma_fish_size_n");
  #ifndef NO_AD_INITIALIZE
    digamma_fish_size_n.initialize();
  #endif
  digamma_fish_esc.allocate(1,nyrs_fish_size,1,nlenbins,"digamma_fish_esc");
  #ifndef NO_AD_INITIALIZE
    digamma_fish_esc.initialize();
  #endif
  priors.allocate(1,5,"priors");
  #ifndef NO_AD_INITIALIZE
    priors.initialize();
  #endif
  Like.allocate("Like");
  #ifndef NO_AD_INITIALIZE
  Like.initialize();
  #endif
  obj_fun.allocate("obj_fun");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  xdum2.allocate(styr,endyr,"xdum2");
  #ifndef NO_AD_INITIALIZE
    xdum2.initialize();
  #endif
  pred_catch.allocate(styr,endyr,"pred_catch");
  #ifndef NO_AD_INITIALIZE
    pred_catch.initialize();
  #endif
  obs_catch.allocate(styr,endyr,"obs_catch");
  #ifndef NO_AD_INITIALIZE
    obs_catch.initialize();
  #endif
  N_proj.allocate(endyr+1,endyr+15,1,nages_M,"N_proj");
  #ifndef NO_AD_INITIALIZE
    N_proj.initialize();
  #endif
  FABC_proj.allocate("FABC_proj");
  #ifndef NO_AD_INITIALIZE
  FABC_proj.initialize();
  #endif
  FABC_tot_proj.allocate(1,nages_M,"FABC_tot_proj");
  #ifndef NO_AD_INITIALIZE
    FABC_tot_proj.initialize();
  #endif
  FOFL_proj.allocate("FOFL_proj");
  #ifndef NO_AD_INITIALIZE
  FOFL_proj.initialize();
  #endif
  FOFL_tot_proj.allocate(1,nages_M,"FOFL_tot_proj");
  #ifndef NO_AD_INITIALIZE
    FOFL_tot_proj.initialize();
  #endif
  ABC.allocate("ABC");
  B40.allocate("B40");
  OFL.allocate("OFL");
  #ifndef NO_AD_INITIALIZE
  OFL.initialize();
  #endif
  Z_proj.allocate(1,nages_M,"Z_proj");
  #ifndef NO_AD_INITIALIZE
    Z_proj.initialize();
  #endif
  ZOFL_proj.allocate(1,nages_M,"ZOFL_proj");
  #ifndef NO_AD_INITIALIZE
    ZOFL_proj.initialize();
  #endif
  S_proj.allocate(1,nages_M,"S_proj");
  #ifndef NO_AD_INITIALIZE
    S_proj.initialize();
  #endif
  catage_proj.allocate(endyr+1,endyr+15,1,nages_M,"catage_proj");
  #ifndef NO_AD_INITIALIZE
    catage_proj.initialize();
  #endif
  catage_proj_OFL.allocate(endyr+1,endyr+15,1,nages_M,"catage_proj_OFL");
  #ifndef NO_AD_INITIALIZE
    catage_proj_OFL.initialize();
  #endif
  pred_catch_proj.allocate(endyr+1,endyr+15,"pred_catch_proj");
  #ifndef NO_AD_INITIALIZE
    pred_catch_proj.initialize();
  #endif
  pred_catch_proj_OFL.allocate(endyr+1,endyr+15,"pred_catch_proj_OFL");
  #ifndef NO_AD_INITIALIZE
    pred_catch_proj_OFL.initialize();
  #endif
  spawn_biom_proj.allocate(endyr+1,endyr+15,"spawn_biom_proj");
  tot_biom_proj.allocate(endyr+1,endyr+15,"tot_biom_proj");
  stdev_rec.allocate("stdev_rec");
  #ifndef NO_AD_INITIALIZE
  stdev_rec.initialize();
  #endif
  FOFL.allocate("FOFL");
  #ifndef NO_AD_INITIALIZE
  FOFL.initialize();
  #endif
  FABC.allocate("FABC");
  #ifndef NO_AD_INITIALIZE
  FABC.initialize();
  #endif
  FOFL2.allocate("FOFL2");
  #ifndef NO_AD_INITIALIZE
  FOFL2.initialize();
  #endif
  FABC2.allocate("FABC2");
  #ifndef NO_AD_INITIALIZE
  FABC2.initialize();
  #endif
}

void model_parameters::userfunction(void)
{
  obj_fun =0.0;
  ofstream& evalout= *pad_evalout;
  l=l+1;
  Get_Selectivity();                // Call function to get selectivities
  get_prop_mat();
  Get_Mortality_Rates();            // Call function to get fishing and natural mortality
  Get_Bzero();                      // OjO
  Get_Numbers_At_Age();             // Call function to get numbers at age per year
  Get_Catch_at_Age();               // Call function to get catch at age per year
  Get_Predicted_Values();           // Get predicted values for catch, survbio, age and size comps
  if (last_phase())
  {
   Get_Dependent_Vars();            // Solve for dependent variables like total bio, recruitment etc.
   compute_spr_rates();	            // Compute f40 etc.
   Get_Population_Projection();
  }
  Evaluate_Objective_Function();    // Minimize objective function value
 if (mceval_phase())											// For outputting MCMC simulations in text format
  {
     evalout<<sigr<<" "<<q_srv1<<" "<<q_srv2<<" "<<F40<<" "<<natmort<<" "<<spawn_biom_proj(endyr+1)<<" "<<ABC<<" "<<obj_fun<<" "<<tot_biom<<" "<<log_rec_dev<<" "<<spawn_biom<<" "<<log_mean_rec<<" "<<spawn_biom_proj<<" "<<pred_catch_proj<<" "<<N_proj(endyr+1,1)<<" "<<N_proj(endyr+2,1)<<" "<<N_proj(endyr+3,1)<<" "<<N_proj(endyr+4,1)<<" "<<N_proj(endyr+5,1)<<" "<<N_proj(endyr+6,1)<<" "<<N_proj(endyr+7,1)<<" "<<N_proj(endyr+8,1)<<" "<<N_proj(endyr+9,1)<<" "<<N_proj(endyr+10,1)<<" "<<tot_biom_proj<<endl;
     //tot_biom_proj(endyr+1)<<endl;
  }
}

void model_parameters::Get_Bzero(void)
{
  ofstream& evalout= *pad_evalout;
  if(wt_Rzero>0) sigr = sigrprior;   // Fixes sigr for northern basecase
  Bzero.initialize();
  Rzero    =  mfexp(log_Rzero); 
  sigrsq   = sigr*sigr;
  dvariable survtmp = exp(-natmort);
  dvariable spawn_adj=pow(survtmp,spawn_fract) ;
  dvar_matrix natagetmp(styr_rec,styr,1,nages_M);
  // First "year"--truly at equilibrium....in styr-nages
  natagetmp(styr_rec,1)      = Rzero;
  for (j=2; j<=nages_M; j++)
    natagetmp(styr_rec,j)    = natagetmp(styr_rec,j-1) * survtmp;
  natagetmp(styr_rec,nages_M) /= (1.-survtmp); 
  Bzero   = wt_mature * spawn_adj *natagetmp(styr_rec) ;
  phizero = Bzero/Rzero;
  // Subsequent pre-history years
  Sp_Biom.initialize();
  Sp_Biom(styr_sp,styr_rec-1) = Bzero;
  for (i=styr_rec;i<styr;i++)
  {
    if(styr_rec!=i) 
      natagetmp(i,1)        = mfexp(log_rec_dev(i) + log_Rzero); // note that this sets up initial age composition to be from Rzero, not MeanRec
    Sp_Biom(i)              = natagetmp(i)*spawn_adj * wt_mature; 
    natagetmp(i+1)(2,nages_M) = ++(natagetmp(i)(1,nages_M-1)*mfexp(-natmort ));
    natagetmp(i+1,nages_M)   += natagetmp(i,nages_M)*mfexp(-natmort);
  }
  // First year of "model"
  natagetmp(styr,1)   = mfexp(log_rec_dev(styr) + log_Rzero);
  natage(styr)        = natagetmp(styr); 
  sam_rec(styr_rec,styr) = column(natagetmp,1);
  Sp_Biom(styr)          = natagetmp(styr)*spawn_adj * wt_mature; 
  // Set Sr params 
  switch (SrType)
  {
    case 1:
      alpha  = log(-4.*steepness/(steepness-1.));
      break;
    case 2:
      alpha  =  Bzero * (1. - (steepness - 0.2) / (0.8*steepness) ) / Rzero;
      beta   = (5. * steepness - 1.) / (4. * steepness * Rzero);
      break;
    case 4:
      beta   = log(5.*steepness)/(0.8*Bzero) ;
      alpha  = log(Rzero/Bzero)+beta*Bzero;
      break;
  }
}

void model_parameters::Get_Selectivity(void)
{
  ofstream& evalout= *pad_evalout;
    for (j=1;j<=nages_M;j++)
      fish_sel(j) = 1./(1. + mfexp(-2.944438979*(double(j+1)-a50)/delta));
      for (j=1;j<=nages_M;j++)
        srv1_sel(j) = 1./(1. + mfexp(-2.944438979*(double(j+1)-a50_srv1)/delta_srv1));
    srv2_sel       = srv1_sel;
}

void model_parameters::get_prop_mat(void)
{
  ofstream& evalout= *pad_evalout;
  for (int i=1;i<=nages_mat;i++){
   if(L_tot_na(i)>0){L_pmat(i) = L_mat_na(i)/L_tot_na(i);}
   if(C_tot_na(i)>0){C_pmat(i) = C_mat_na(i)/C_tot_na(i);}
   Pred_pmat(i) = 1/(1+exp(-1.0*mat_delta*(i-mat_a50)));}
  for (int i=1;i<=nages_M;i++){
  //p_mature(i) = value(Pred_pmat(i));
  p_mature(i) = 1/(1+exp(-1.0*mat_delta*((i+1)-mat_a50)));
  wt_mature(i) = wt(i)*p_mature(i)/2;}
}

void model_parameters::Get_Mortality_Rates(void)
{
  ofstream& evalout= *pad_evalout;
  natmort = mfexp(logm);			// setting natural mortality to arithmetic scale
   if(ph_m>0) nattymort=natmort;
   else nattymort=log_mean_rec;
  Fmort = mfexp(log_avg_F +  log_F_devs);	//setting fishing mortaltiy to arithmetic scale
  for (iyr=styr; iyr<=endyr; iyr++)
    F(iyr) = Fmort(iyr) * fish_sel;		// Getting fully selected fishing mortality
  Z = F + natmort;				// Fully selected total mortality
  S = mfexp(-1.0*Z);
}

void model_parameters::Get_Numbers_At_Age(void)
{
  ofstream& evalout= *pad_evalout;
  if (SrType==3)
  {
    int itmp;
    for (j=2;j<nages_M;j++)
    {
      itmp = styr+1-j;
      natage(styr,j) = mfexp(log_mean_rec  - natmort * double(j-1)+ log_rec_dev(itmp)); 
    }
    natage(styr,nages_M) = mfexp(log_mean_rec - natmort * (nages_M-1)) / (1. - exp(-natmort) );
  }
  for ( i=styr;i < endyr;i++)
  {
    natage(i,1)           = mfexp(log_rec_dev(i) + log_mean_rec );
    natage(i+1)(2,nages_M)  = ++elem_prod(natage(i)(1,nages_M-1),S(i)(1,nages_M-1));       // Following year
    natage(i+1,nages_M)    += natage(i,nages_M)*S(i,nages_M);
    sam_rec(i)            = natage(i,1); 
     Sp_Biom(i)           = natage(i) * wt_mature;  // Old way
 }
  natage(endyr,1)         = mfexp(log_rec_dev(endyr) + log_mean_rec ); 
  sam_rec(endyr)          = natage(endyr,1); 
  Sp_Biom(endyr)          = natage(endyr)* wt_mature; //Old way
}

void model_parameters::Get_Catch_at_Age(void)
{
  ofstream& evalout= *pad_evalout;
  pred_catch_early.initialize();
  pred_catch_later.initialize();
  for (iyr=styr; iyr<=yr_catchwt; iyr++)
  {
    catage(iyr) = elem_div(elem_prod(elem_prod(natage(iyr),F(iyr)),(1.-S(iyr))),Z(iyr));
    pred_catch_early(iyr) = catage(iyr)*wt;
  }
  for (iyr=yr_catchwt+1; iyr<=endyr; iyr++)
  {
    catage(iyr) = elem_div(elem_prod(elem_prod(natage(iyr),F(iyr)),(1.-S(iyr))),Z(iyr));
    pred_catch_later(iyr) = catage(iyr)*wt;
  }
}

void model_parameters::Get_Dependent_Vars(void)
{
  ofstream& evalout= *pad_evalout;
  for (i=styr;i<=endyr;i++)
  {
    pred_rec(i)   = natage(i,1);  					    // Setting up results based on estimated paramters
    tot_biom(i)   = wt * natage(i);				     	// Total biomass results
    expl_rate(i)  = pred_catch(i)/tot_biom(i);  // Setting up results based on estimated paramters
    spawn_biom(i) = Sp_Biom(i) ;		            // Spawning biomass result
  }
  avg_rec        = mean(pred_rec);
  Depletion      = spawn_biom(endyr)/spawn_biom(styr);      // 1-Depletion
  spbiom_trend   = spawn_biom(endyr)/spawn_biom(endyr-1);
}

void model_parameters::Get_Predicted_Values(void)
{
  ofstream& evalout= *pad_evalout;
  nsamples_fish_age_ts.initialize();
  nsamples_srv1_age_ts.initialize();
  nsamples_fish_size_ts.initialize();
  nsamples_srv1_size_ts.initialize();
  nsamples_srv2_size_ts.initialize();
  q_srv1         = exp(log_q_srv1);                         // Survey catchability at arithmetic scale
  q_srv2         = exp(log_q_srv2);                         // Survey catchability at arithmetic scale
  for (i=1;i<=nyrs_srv1;i++)
    pred_srv1(i) = q_srv1 * (natage(yrs_srv1(i))*elem_prod(srv1_sel,wt));  	// Predicted Survey biomass
 if(nyrs_srv2>0) {  for (i=1;i<=nyrs_srv2;i++)
    pred_srv2(i) = q_srv2 * (natage(yrs_srv2(i))*elem_prod(srv2_sel,wt));  	// Predicted Survey biomass
 }
  for (i=1;i<=nyrs_fish_age;i++) {
   eac_fish(i)  = catage(yrs_fish_age(i))/sum(catage(yrs_fish_age(i)))
                   * ageage(age_age_ind_fsh(i));            // Predicted Fishery age comps
   effn_fish_age(i) = (1-eac_fish(i))*eac_fish(i)/norm2(oac_fish(i)-eac_fish(i));	//effective n fish age comps
   sdnr_fish_age(i) = sdnr(eac_fish(i),oac_fish(i),double(nmulti_fish_age(i)));	//sdnr fish age comps	
   nsamples_fish_age_ts(yrs_fish_age(i))=nmulti_fish_age(i);  
   effn_fish_age_ts(yrs_fish_age(i))=effn_fish_age(i);  
   sdnr_fish_age_ts(yrs_fish_age(i))=sdnr_fish_age(i);}
  for (i=1;i<=nyrs_srv1_age;i++) {
   eac_srv1(i)  = elem_prod(srv1_sel,natage(yrs_srv1_age(i)))/(natage(yrs_srv1_age(i)) 
                   * srv1_sel)* ageage(age_age_ind_srv(i)); // Predicted Survey age comps
   effn_srv1_age(i) = (1-eac_srv1(i))*eac_srv1(i)/norm2(oac_srv1(i)-eac_srv1(i));
   sdnr_srv1_age(i) = sdnr(eac_srv1(i),oac_srv1(i),double(nmulti_srv1_age(i)));
   nsamples_srv1_age_ts(yrs_srv1_age(i))=nmulti_srv1_age(i);  
   effn_srv1_age_ts(yrs_srv1_age(i))=effn_srv1_age(i);  
   sdnr_srv1_age_ts(yrs_srv1_age(i))=sdnr_srv1_age(i);}
  for (i=1;i<=nyrs_fish_size;i++) {											           // Lets you use a second matrix for part of it
   esc_fish(i)  = catage(yrs_fish_size(i))/sum(catage(yrs_fish_size(i)))
                  * sizeage(siz_age_ind_fsh(i));                                       // Second Predicted Fishery size comps for 80s and 90s
   effn_fish_size(i) = (1-esc_fish(i))*esc_fish(i)/norm2(osc_fish(i)-esc_fish(i));	   // effective n fish size comps			
   sdnr_fish_size(i) = sdnr(esc_fish(i),osc_fish(i),double(nmulti_fish_size(i)));    // sdnr fish size comps
   nsamples_fish_size_ts(yrs_fish_size(i))=nmulti_fish_size(i);  
   effn_fish_size_ts(yrs_fish_size(i))=effn_fish_size(i);  
   sdnr_fish_size_ts(yrs_fish_size(i))=sdnr_fish_size(i);}
  for ( i=1;i<=nyrs_srv1_size;i++) {
   esc_srv1(i)  = elem_prod(srv1_sel,natage(yrs_srv1_size(i)))
                   /(natage(yrs_srv1_size(i)) * srv1_sel)* sizeage(siz_age_ind_srv1(i));   // Predicted Survey size comps (not used in POP model)
   effn_srv1_size(i) = (1-esc_srv1(i))*esc_srv1(i)/norm2(osc_srv1(i)-esc_srv1(i)); 		   // effective n survey 1 size comps		
   sdnr_srv1_size(i) = sdnr(esc_srv1(i),osc_srv1(i),double(nmulti_srv1_size(i)));        // sdnr survey 1 size comps
   nsamples_srv1_size_ts(yrs_srv1_size(i))=nmulti_srv1_size(i);  
   effn_srv1_size_ts(yrs_srv1_size(i))=effn_srv1_size(i);  
   sdnr_srv1_size_ts(yrs_srv1_size(i))=sdnr_srv1_size(i);}
  if(nyrs_srv2>0) {
   for ( i=1;i<=nyrs_srv2_size;i++) {
    esc_srv2(i)  = elem_prod(srv2_sel,natage(yrs_srv2_size(i)))
                   /(natage(yrs_srv2_size(i)) * srv2_sel)* sizeage(siz_age_ind_srv2(i));   // Predicted Survey size comps (not used in POP model)
    effn_srv2_size(i) = (1-esc_srv2(i))*esc_srv2(i)/norm2(osc_srv2(i)-esc_srv2(i));        // effective n survey 2 size comps
	sdnr_srv2_size(i) = sdnr(esc_srv2(i),osc_srv2(i),double(nmulti_srv2_size(i)));       // sdnr survey 2 size comps
   nsamples_srv2_size_ts(yrs_srv2_size(i))=nmulti_srv2_size(i);  
   effn_srv2_size_ts(yrs_srv2_size(i))=effn_srv2_size(i);  
   sdnr_srv2_size_ts(yrs_srv2_size(i))=sdnr_srv2_size(i);}}	
  if (nyrs_cpue>0)
  {
    int yy;
    for (i=1;i<=nyrs_cpue;i++) 
    {
      yy           = yrs_cpue(i);
      pred_cpue(i) = wt*elem_div(elem_prod(natage(yy),fish_sel),Z(yy)); 
    } 
    q_cpue = mean_obs_cpue/mfexp(mean(log(pred_cpue)));
    pred_cpue *=  q_cpue;
  } 
  pred_catch(styr,yr_catchwt)    = pred_catch_early;
  pred_catch(yr_catchwt+1,endyr) = pred_catch_later;
  obs_catch(styr,yr_catchwt)    = obs_catch_early;
  obs_catch(yr_catchwt+1,endyr) = obs_catch_later;
  if(ph_q_srv2>0) q2=mfexp(log_q_srv2); else q2=mfexp(log_q_srv1);
  cigar= sigr;
  LMR = log_mean_rec;
}

double model_parameters::round(double r)
{
  ofstream& evalout= *pad_evalout;
}

void model_parameters::Get_Population_Projection(void)
{
  ofstream& evalout= *pad_evalout;
  int k;
  if(mceval_phase()) {
  stdev_rec = sqrt(norm2(value(log_rec_dev(1979,endyr-recage))-mean(value(log_rec_dev(1979,endyr-recage))))/(size_count(value(log_rec_dev(1979,endyr-recage)))-1));
   k=round(value(stdev_rec)*10000);
     N_proj(endyr+1,1)= mfexp(value(log(mean(value(pred_rec(1979,endyr-recage))))-square(stdev_rec)/2+stdev_rec*randn(k+l)));
      cout<<stdev_rec<<" "<<k<<" "<<l<<" "<<endl;
   }
  else {      	N_proj(endyr+1,1)= value(mean(pred_rec(1979,endyr-recage))); }
    for (j=1; j<nages_M-1;j++) {
           N_proj(endyr+1,j+1)=natage(endyr,j)*S(endyr,j); }
    N_proj(endyr+1,nages_M) = natage(endyr,nages_M-1)*S(endyr,nages_M-1)+ natage(endyr,nages_M)*S(endyr,nages_M);
   tot_biom_proj(endyr+1)=N_proj(endyr+1)*wt;
   spawn_biom_proj(endyr+1) = elem_prod(N_proj(endyr+1),pow(mfexp(-yieldratio*FABC_tot_proj-natmort),spawn_fract)) * wt_mature;
  for (i=endyr+1;i<=endyr+15;i++)
  {
    if (spawn_biom_proj(i)/B40 > 1.) {
      FABC_proj =F40;
      FOFL_proj = F35; }
    else {
      FABC_proj = F40 * (spawn_biom_proj(i)/B40 - 0.05)/(1 - 0.05); 
	  FOFL_proj = F35*(spawn_biom_proj(i)/B40 - 0.05)/(1 - 0.05);  }
    for (j=1;j<=nages_M;j++)
    {  
      FOFL_tot_proj(j) = fish_sel(j)*FOFL_proj;
      FABC_tot_proj(j) = fish_sel(j)* FABC_proj ;
      Z_proj(j)   =FABC_tot_proj(j)+ natmort;
      ZOFL_proj(j)   = FOFL_tot_proj(j)+ natmort;
      S_proj(j)   = mfexp(-1.0* Z_proj(j));
      }
    for (j=1;j<=nages_M;j++)
     { 
      catage_proj(i,j) = yieldratio*N_proj(i,j)* FABC_tot_proj(j)/Z_proj(j)*(1.-S_proj(j));
      catage_proj_OFL(i,j) = yieldratio*N_proj(i,j)* FOFL_tot_proj(j)/ZOFL_proj(j)*(1.-mfexp(-ZOFL_proj(j)));
     }
    pred_catch_proj(i)     = catage_proj(i)*wt/yieldratio;
    pred_catch_proj_OFL(i)     = catage_proj_OFL(i)*wt/yieldratio;
    if (i < endyr+15)
    {
 if(mceval_phase()) {
  stdev_rec = sqrt(norm2(value(log_rec_dev(1979,endyr-recage))-mean(value(log_rec_dev(1979,endyr-recage))))/(size_count(value(log_rec_dev(1979,endyr-recage)))-1));
     k=round(value(spawn_biom(endyr)*10000))+i;
  k=k+i;
     N_proj(i+1,1)= mfexp((log(mean(value(pred_rec(1979,endyr-recage))))-square(stdev_rec)/2+stdev_rec*randn(k+l)));  }
     else { 	N_proj(i+1,1)= value(mean(pred_rec(1979,endyr-recage))); }
      for (j=1; j<nages_M-1;j++) {
        N_proj(i+1,j+1) = N_proj(i,j)       * mfexp(-yieldratio*FABC_tot_proj(j)-natmort); }
      N_proj(i+1,nages_M) = N_proj(i,nages_M-1) * mfexp(-yieldratio*FABC_tot_proj(nages_M-1)-natmort)+ N_proj(i,nages_M)   * mfexp(-yieldratio*FABC_tot_proj(nages_M)-natmort);
     //  spawn_biom_proj(i+1)    = N_proj(i+1)*wt_mature;  //Old way
       spawn_biom_proj(i+1)        = elem_prod(N_proj(i+1),pow(mfexp(-yieldratio*FABC_tot_proj-natmort),spawn_fract)) * wt_mature;  // Right way
       tot_biom_proj(i+1)=N_proj(i+1)*wt;
   }
    }
     if (spawn_biom_proj(endyr+1)/B40 > 1.) {
      FABC = F40;
      FOFL = F35; 
      FABC2 = F40;
      FOFL2 = F35; }
    else {
      FABC = F40 * (spawn_biom_proj(endyr+1)/B40 - 0.05)/(1 - 0.05); 
	  FOFL = F35*(spawn_biom_proj(endyr+1)/B40 - 0.05)/(1 - 0.05);  
      FABC2 = F40 * (spawn_biom_proj(endyr+2)/B40 - 0.05)/(1 - 0.05); 
	  FOFL2 = F35*(spawn_biom_proj(endyr+2)/B40 - 0.05)/(1 - 0.05);  }
     OFL=pred_catch_proj_OFL(endyr+1);
     ABC=pred_catch_proj(endyr+1);
}

void model_parameters::Evaluate_Objective_Function(void)
{
  ofstream& evalout= *pad_evalout;
  ssqcatch.initialize();
  Like.initialize();
  rec_like.initialize();
  cpue_like.initialize();
  surv_like.initialize();		// Likelihood values for survey biomasses, allowance for up to 3 surveys
  age_like.initialize();			// Likelihood values for age and size compositions allowance for up 6 comps
  sel_like.initialize();			// LIkelihood values for selectivities with alowance for up to 6 selectivities
  F_mort_regularity.initialize();
  avg_sel_penalty.initialize();
  Catch_Likelihood();
  Surv_Likelihood();                            // Likelihood function for survey biomass
  Size_Age_Like();                              // Multinomial likelihood
  Calc_priors();                                // Priors
  Sel_Like();                                   // Penalty function for selectivity
  Rec_Like();                                   // Penalty function for selectivity
  if(active(log_F_devs))                        // Penalty function for fishing mortality deviations
    F_mort_regularity  = wt_fmort_reg * norm2(log_F_devs);
  if (active(log_fish_sel_coffs))
    avg_sel_penalty   = square(log_avgfishsel); // Average fishery selectivity penalty
  if (active(log_srv1_sel_coffs))
    avg_sel_penalty += square(log_avgsrv1sel);  // Average survey selectivity penalty
  if (active(log_srv2_sel_coffs))
    avg_sel_penalty += square(log_avgsrv2sel);  // Average survey selectivity penalty
  Like              += ssqcatch ;
  Like              += cpue_like;
  Like              += sum(surv_like);
  Like              += sum(age_like);
  obj_fun           += Like;                      // Put here to capture the data likelihood
  obj_fun           += wt_rec_var *rec_like;
  obj_fun           += sum(sel_like);
  if(active(log_fish_sel_coffs))
    obj_fun         += wt_avg_sel*avg_sel_penalty;       // Conditions the model to have mean selectivity =1
  if(active(log_F_devs))                          // Penalty function for fishing mortality deviations
    obj_fun         += F_mort_regularity;
  obj_fun           += sum(priors);               //Add priors
  if (active(mF50)&&last_phase())
    obj_fun         += sprpen;                    // To solve for the F40 etc.     
  if (current_phase()<3)
      obj_fun         += norm2(F);        
    obj_fun         += wt_Rzero*square(log_Rzero-log_mean_rec);   // Penalty early on to scale population...                
  // reassiging variables for northern model outputs and hessian
   if(ph_sigr>0) cigar=sigr; else cigar=log_mean_rec;
  for (int i=1;i<=nages_mat;i++){
   if(L_tot_na(i)>0){
    Like_L_vec(i) = L_wt(i)*(log(L_mat_na(i)+k)+gammln(L_mat_na(i)+k)+log((L_tot_na(i)-L_mat_na(i))+k)+gammln((L_tot_na(i)-L_mat_na(i))+k)-log(L_tot_na(i)+k)-gammln(L_tot_na(i)+k)-(L_mat_na(i)+k)*log(Pred_pmat(i)+k)-((L_tot_na(i)-L_mat_na(i))+k)*log((1-Pred_pmat(i))+k));}
   if(C_tot_na(i)>0){
    Like_C_vec(i) = C_wt(i)*(log(C_mat_na(i)+k)+gammln(C_mat_na(i)+k)+log((C_tot_na(i)-C_mat_na(i))+k)+gammln((C_tot_na(i)-C_mat_na(i))+k)-log(C_tot_na(i)+k)-gammln(C_tot_na(i)+k)-(C_mat_na(i)+k)*log(Pred_pmat(i)+k)-((C_tot_na(i)-C_mat_na(i))+k)*log((1-Pred_pmat(i))+k));}}
  Like_L = sum(Like_L_vec);
  Like_C = sum(Like_C_vec);
  obj_fun   += Like_L;
  obj_fun   += Like_C;
}

void model_parameters::Catch_Likelihood(void)
{
  ofstream& evalout= *pad_evalout;
  ssqcatch  +=  wt_ssqcatch *norm2(log(obs_catch_early+.00001)-log(pred_catch_early+.00001));
  ssqcatch  +=  wt_ssqcatch2 *norm2(log(obs_catch_later+.00001)-log(pred_catch_later+.00001));
}

void model_parameters::Surv_Likelihood(void)
{
  ofstream& evalout= *pad_evalout;
  for (i=1; i<=nyrs_srv1; i++) {
      if(ph_logsurv>1) 
	  surv_like(1) += square((log(obs_srv1_biom(i))-log(pred_srv1(i)) ))/ (2.*square(obs_srv1_se(i)/obs_srv1_biom(i))); // likelihood for survey biomass
	  else 
	  	  surv_like(1) += square(obs_srv1_biom(i)-pred_srv1(i) )/ (2.*square(obs_srv1_se(i))); }  // likelihood for survey biomass
 if(nyrs_srv2>0) {  for (i=1; i<=nyrs_srv2; i++) {
      if(ph_logsurv>1) 
	  surv_like(2) += square((log(obs_srv2_biom(i))-log(pred_srv2(i)) ))/ (2.*square(obs_srv2_se(i)/obs_srv2_biom(i))); // likelihood for survey biomass
	  else 
	  	  surv_like(2) += square(obs_srv2_biom(i)-pred_srv2(i) )/ (2.*square(obs_srv2_se(i))); } // likelihood for survey biomass
 }
  surv_like(1) *= wt_srv1 ;  
  surv_like(2) *= wt_srv2 ;  
  if (nyrs_cpue>0)
    cpue_like  = norm2(log(obs_cpue)-log(pred_cpue)) / (2.*cv_cpue*cv_cpue); // likelihood for fishery cpue
}

void model_parameters::Size_Age_Like(void)
{
  ofstream& evalout= *pad_evalout;
  for (i=1; i <= nyrs_fish_age; i++)
    age_like(1) -= nmulti_fish_age(i)*((oac_fish(i) + 0.00001) * log(eac_fish(i) + 0.00001)) ;
  age_like(1)   -= offset(1);                       // Subract offsets
  for (i=1; i <= nyrs_srv1_age; i++)
    age_like(2) -= nmulti_srv1_age(i)*((oac_srv1(i) + 0.00001) * log(eac_srv1(i) + 0.00001)) ;
  age_like(2)   -= offset(2);                       // Subract offsets
  for (i=1; i <= nyrs_fish_size; i++)
    age_like(3) -= nmulti_fish_size(i)*((osc_fish(i) + 0.00001) * log(esc_fish(i) + 0.00001)) ;
  age_like(3)   -= offset(3);                       // Subract offsets
  for (i=1; i <= nyrs_srv1_size; i++)
    age_like(4) -= nmulti_srv1_size(i)*((osc_srv1(i) + 0.00001) * log(esc_srv1(i) + 0.00001)) ;
  age_like(4)   -= offset(4);                       // Subract offsets
 if(nyrs_srv2>0) {
  for (i=1; i <= nyrs_srv2_size; i++)
    age_like(5) -= nmulti_srv2_size(i)*((osc_srv2(i) + 0.00001) * log(esc_srv2(i) + 0.00001)) ;
  age_like(5)   -= offset(5);   }                    // Subract offsets
  age_like(1) *= wt_fish_age;                   // Multiple each likelihood by their weights from .ctl file
  age_like(2) *= wt_srv1_age;
  age_like(3) *= wt_fish_size;
  age_like(4) *= wt_srv1_size;
  age_like(5) *= wt_srv2_size;
}

void model_parameters::Calc_priors(void)
{
  ofstream& evalout= *pad_evalout;
    priors.initialize();
    if (active(sigr))
      priors(1)    = square(log(sigr/sigrprior))/(2.*square(cvsigrprior));
    if (active(log_q_srv1))
      priors(2)    = square(log_q_srv1-log_q_srv1prior)/(2.*square(cvq_srv1prior));
    if (active(log_q_srv2))
      priors(5)    = square(log_q_srv2-log_q_srv2prior)/(2.*square(cvq_srv2prior));
    if (active(steepness))
      priors(3)    = square(log(steepness/steep_prior))/(2.*cv_steep_prior); // not used in POP model
    if (active(logm))
      priors(4)    = square(logm-log(mprior))/(2.*square(cvmprior));
}

void model_parameters::Sel_Like(void)
{
  ofstream& evalout= *pad_evalout;
  sel_like(1)   +=wt_sel_reg_fish * norm2(first_difference(first_difference(log_fish_sel)));  // Constrains selectivities to be smooth
  if (active(log_srv1_sel_coffs) )
    sel_like(2) +=wt_sel_reg_srv1 * norm2(first_difference(first_difference(log_srv1_sel)));  // Constrains selectivities to be smooth
  if (active(log_srv2_sel_coffs) )
    sel_like(3) +=wt_sel_reg_srv2 * norm2(first_difference(first_difference(log_srv2_sel)));  // Constrains selectivities to be smooth
  for (j=1;j<nages_M;j++)
    if (log_fish_sel(j)>log_fish_sel(j+1))
      sel_like(4) += wt_sel_dome_fish *square(log_fish_sel(j)-log_fish_sel(j+1));  //Prevents dome-shapedness
  if (active(log_srv1_sel_coffs) ) {
    for (j=1;j<nages_M;j++)
      if (log_srv1_sel(j)>log_srv1_sel(j+1))
        sel_like(5) +=wt_sel_dome_srv1 *square(log_srv1_sel(j)-log_srv1_sel(j+1)); }
  if (active(log_srv2_sel_coffs) ) {
    for (j=1;j<nages_M;j++)
      if (log_srv2_sel(j)>log_srv2_sel(j+1))
        sel_like(6) +=wt_sel_dome_srv2 *square(log_srv2_sel(j)-log_srv2_sel(j+1)); }
}

void model_parameters::Rec_Like(void)
{
  ofstream& evalout= *pad_evalout;
  if(SrType== 3)
  {
    if (rec_like_type==1)
      rec_like      = norm2(log_rec_dev)/(2*square(sigr)) + (size_count(log_rec_dev)*log(sigr));
    else
      if (active(sigr))
        rec_like      = norm2(log_rec_dev+sigr*sigr/2.)/(2.*square(sigr)) + size_count(log_rec_dev)*log(sigr);
      else 
        rec_like      = norm2(log_rec_dev+sigr*sigr/2.)/(2.*square(sigr)) ;
  }
  else  
  {
    dvar_vector stmp(styr_rec,endyr);
    for (i=styr_rec;i<=endyr;i++)
      stmp(i) = Sp_Biom(i-recage);
    srm_rec   = SRecruit(stmp);
    dvar_vector   chi(styr_rec_est,endyr_rec_est);
    chi         = log(elem_div(sam_rec(styr_rec_est,endyr_rec_est) , srm_rec(styr_rec_est,endyr_rec_est)));
    dvariable SSQRec = norm2( chi + sigrsq/2.) ;
    rec_like    = .5*SSQRec/sigrsq + nrecs_est*log(sigr); 
    rec_like   += .5*norm2( log_rec_dev(styr_rec,styr_rec_est) )/sigrsq + (styr_rec_est-styr_rec)*log(sigr) ; 
    if (endyr>endyr_rec_est)
      rec_like += .5*norm2( log_rec_dev(endyr_rec_est,endyr  ) )/sigrsq + (endyr-endyr_rec_est)  *log(sigr) ; 
  }
}

dvar_vector model_parameters::SRecruit(const dvar_vector& Stmp)
{
  ofstream& evalout= *pad_evalout;
  RETURN_ARRAYS_INCREMENT();
  dvar_vector RecTmp(Stmp.indexmin(),Stmp.indexmax());
      RecTmp = mfexp(log_mean_rec);                    //Avg recruitment
  RETURN_ARRAYS_DECREMENT();
  return RecTmp;
}

void model_parameters::compute_spr_rates(void)
{
  ofstream& evalout= *pad_evalout;
  //Compute SPR Rates and add them to the likelihood for Females 
  SB0=0.;
  SBF50=0.;
  SBF40=0.;
  SBF35=0.;
  // Scale F-spr rates to be on full-selected values
  F50  = mF50*max(fish_sel);
  F40  = mF40*max(fish_sel);
  F35  = mF35*max(fish_sel);
  for (i=1;i<=4;i++)
    Nspr(i,1)=1.;
  for (j=2;j<nages_M;j++)
  {
    Nspr(1,j)=Nspr(1,j-1)*mfexp(-1.*natmort);
    Nspr(2,j)=Nspr(2,j-1)*mfexp(-1.*(natmort+mF50*fish_sel(j-1)));
    Nspr(3,j)=Nspr(3,j-1)*mfexp(-1.*(natmort+mF40*fish_sel(j-1)));
    Nspr(4,j)=Nspr(4,j-1)*mfexp(-1.*(natmort+mF35*fish_sel(j-1)));
  }
  Nspr(1,nages_M)=Nspr(1,nages_M-1)*mfexp(-1.*natmort)/(1.-mfexp(-1.*natmort));
  Nspr(2,nages_M)=Nspr(2,nages_M-1)*mfexp(-1.* (natmort+mF50*fish_sel(nages_M-1)))/(1.-mfexp(-1.*(natmort+mF50*fish_sel(nages_M))));
  Nspr(3,nages_M)=Nspr(3,nages_M-1)*mfexp(-1.* (natmort+mF40*fish_sel(nages_M-1)))/ (1.-mfexp(-1.*(natmort+mF40*fish_sel(nages_M))));
  Nspr(4,nages_M)=Nspr(4,nages_M-1)*mfexp(-1.* (natmort+mF35*fish_sel(nages_M-1)))/ (1.-mfexp(-1.*(natmort+mF35*fish_sel(nages_M))));
  for (j=1;j<=nages_M;j++)
  {
   // Kill them off till (spawn_fract)
    SB0    += Nspr(1,j)*wt_mature(j)*mfexp(-spawn_fract*natmort);
    SBF50  += Nspr(2,j)*wt_mature(j)*mfexp(-spawn_fract*(natmort+mF50*fish_sel(j)));
    SBF40  += Nspr(3,j)*wt_mature(j)*mfexp(-spawn_fract*(natmort+mF40*fish_sel(j)));
    SBF35  += Nspr(4,j)*wt_mature(j)*mfexp(-spawn_fract*(natmort+mF35*fish_sel(j)));
   }
  sprpen    = 100.*square(SBF50/SB0-0.5);
  sprpen   += 100.*square(SBF40/SB0-0.4);
  sprpen   += 100.*square(SBF35/SB0-0.35);
  B40= SBF40*mean(pred_rec(1979,endyr-recage));
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  cout<<"-------------Finished: "<<current_phase()<<" "<<Like<<" "<<age_like<<endl;
  if (last_phase())
    write_proj();
  report<<"****Executive Summary Material*****"<<endl;
  report<<"     Model name"     <<endl;
  report<<model_name<<endl;
  report<<"     .dat file"     <<endl;
  report<<data_file<<endl;
  report<<"     Number parameters estimated"     <<endl;
  report<<initial_params::nvarcalc()<<endl;
  report<<"     TotalBiomass for "<<endyr+1<<endl;
  report<<N_proj(endyr+1)*wt<<endl;
  report<<"     TotalBiomass for "<<endyr+2     <<endl;
  report<<N_proj(endyr+2)*wt<<endl;
  report<<"     Female_Spawning Biomass for "<<endyr+1     <<endl;
  report<<spawn_biom_proj(endyr+1)<<endl;
  report<<"     Female_Spawning_Biomass for "<<endyr+2     <<endl;
  report<<spawn_biom_proj(endyr+2)<<endl;
  report<<"     B_zero"     <<endl;
  report<<SB0*mean(pred_rec(1979,endyr-recage))<<endl;
  report<<"     B_40"     <<endl;
  report<<B40<<endl;
  report<<"     B_35"     <<endl;
  report<<SBF35*mean(pred_rec(1979,endyr-recage))<<endl;
  report<<"     F_40"     <<endl;
  report<<F40<<endl;
  report<<"     F_35"     <<endl;
  report<<F35<<endl;
  report<<"     F_ABC for "<<endyr+1     <<endl;
  report<<FABC<<endl;
  report<<"     F_ABC for "<<endyr+2     <<endl;
  report<<FABC2<<endl;
  report<<"     ABC for "<<endyr+1     <<endl;
  report<<pred_catch_proj(endyr+1)<<endl;
  report<<"     ABC for "<<endyr+2     <<endl;
  report<<pred_catch_proj(endyr+2)<<endl;
  report<<"     F_OFL for "<<endyr+1     <<endl;
  report<<FOFL<<endl;
  report<<"     F_OFL for "<<endyr+2     <<endl;
  report<<FOFL2<<endl;
  report<<"     OFL for "<<endyr+1     <<endl;
  report<<OFL<<endl; 
  report<<"     OFL for "<<endyr+2     <<endl;
  report<<pred_catch_proj_OFL(endyr+2)<<endl; 
  report<<"     Total likelihood"     <<endl;
  report<<obj_fun<<endl;
  report<<"     Data likelihood"     <<endl;
  report<<Like<<endl<<endl;
  report<<" ************   Some more parameter estimates and their SDs ************"<<endl;
  if(last_phase()) {
    // add standard deviation data types    
  report<<"   q_trawl   "<<endl;
  report<<q_srv1<<" "<<q_srv1.sd<<endl;
  report<<"   q_longline  "<<endl;
  report<<q_srv2<<" "<<q2.sd<<endl;
  report<<"   nat_mort  "<<endl;
  report<<natmort<<" "<<nattymort.sd<<endl;
  report<<"  sigr   "<<endl;  
  report<<sigr<<" "<<cigar.sd<<endl;  
  report<<"   log_mean_rec"<<endl;
  report<<log_mean_rec<<" "<<LMR.sd<<endl;
  report<<"   F_40"<<endl;
  report<<F40<<" "<<F40.sd<<endl;
  report<<"    tot_biom"<<endl;
  report<<tot_biom_proj(endyr+1)<<" "<<tot_biom_proj.sd(endyr+1)<<endl;
  report<<"   spawn_biom"<<endl;
  report<<spawn_biom_proj(endyr+1)<<" "<<spawn_biom_proj.sd(endyr+1)<<endl;
  report<<"    B40"<<endl;
  report<<B40<<" "<<B40.sd<<endl;
  report<<"   ABC"<<endl;
  report<<ABC<<" "<<ABC.sd<<endl<<endl;
 }
  // end exec
  report<<"************Rest of data output **********"<<endl<<endl;
  report << "Year "<< yy <<endl;
  report << "Pred_Catch "<< pred_catch_early<<" Pred_catch_later "<<pred_catch_later <<endl;
  report << "Obs_Catch "<< obs_catch_early<<" Obs_Catch_Later "<<obs_catch_later <<endl;
  report << "Catch_at_age "<<aa <<endl;
  for (i=styr;i<=endyr;i++) report << i<<" "<<catage(i) <<endl; report<<endl;
  report << "Numbers "<<aa <<endl;
  for (i=styr;i<=endyr;i++) report << i<<" "<<natage(i) <<endl; report<<endl;
  if (nyrs_cpue>0) 
  {
    report <<"Years_CPUE: "     <<yrs_cpue <<endl; 
    report <<"Predicted_CPUE: " <<pred_cpue<<endl; 
    report <<"observed_CPUE: "  <<obs_cpue <<endl; 
    report <<"q_CPUE: "         <<q_cpue   <<endl; 
  }
  report << "Obs_P_fish_age"<<aa <<endl;
  for (i=1;i<=nyrs_fish_age;i++) report << yrs_fish_age(i)<<" "<<oac_fish(i) 
      <<" eff_N "<<(1-eac_fish(i))*eac_fish(i)/norm2(oac_fish(i)-eac_fish(i))  <<" N "<<nmulti_fish_age(i)
      <<" SDNR "<< sdnr(eac_fish(i),oac_fish(i),double(nmulti_fish_age(i)))<<endl; report<<endl;
  report << "Pred_P_fish_age"<<aa <<endl;
  for (i=1;i<=nyrs_fish_age;i++) report << yrs_fish_age(i)<<" "<<eac_fish(i) <<endl; report<<endl;
  report << "Obs_P_fish_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish_size;i++) report << yrs_fish_size(i)<<" "<<osc_fish(i) 
      <<" eff_N "<<(1-esc_fish(i))*esc_fish(i)/norm2(osc_fish(i)-esc_fish(i))  <<" N "<<nmulti_fish_size(i)
      <<" SDNR "<< sdnr(esc_fish(i),osc_fish(i),double(nmulti_fish_size(i)))<<endl; report<<endl;
  report << "Pred_P_fish_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish_size;i++) report << yrs_fish_size(i)<<" "<<esc_fish(i) <<endl; report<<endl;
  report << "Obs_P_srv1_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv1_age;i++) report << yrs_srv1_age(i)<<" "<<oac_srv1(i) 
      <<" eff_N "<<(1-eac_srv1(i))*eac_srv1(i)/norm2(oac_srv1(i)-eac_srv1(i)) <<" N "<<nmulti_srv1_age(i)
      <<" SDNR "<< sdnr(eac_srv1(i),oac_srv1(i),double(nmulti_srv1_age(i)))<<endl; report<<endl;
  report << "Pred_P_srv1_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv1_age;i++) report << yrs_srv1_age(i)<<" "<<eac_srv1(i) <<endl; report<<endl;
  report << "Obs_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) report << yrs_srv1_size(i)<<" "<<osc_srv1(i) 
      <<" eff_N "<<(1-esc_srv1(i))*esc_srv1(i)/norm2(osc_srv1(i)-esc_srv1(i)) <<" N "<<nmulti_srv1_size(i)
      <<" SDNR "<< sdnr(esc_srv1(i),osc_srv1(i),double(nmulti_srv1_size(i)))<<endl; report<<endl;
  report << "Pred_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) report << yrs_srv1_size(i)<<" "<<esc_srv1(i) <<endl; report<<endl;
  report << "Obs_P_srv2_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) report << yrs_srv2_size(i)<<" "<<osc_srv2(i) 
      <<" eff_N "<<(1-esc_srv2(i))*esc_srv2(i)/norm2(osc_srv2(i)-esc_srv2(i)) <<" N "<<nmulti_srv2_size(i)
      <<" SDNR "<< sdnr(esc_srv2(i),osc_srv2(i),double(nmulti_srv2_size(i)))<<endl; report<<endl;
  report << "Pred_P_srv2_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) report << yrs_srv2_size(i)<<" "<<esc_srv2(i) <<endl; report<<endl;
    report << "Survey Biomass " <<endl;
  report << "Year:     " << yrs_srv1  <<endl;
  report << "Predicted:   " << pred_srv1  <<endl;
  report << "Observed:   " << obs_srv1_biom  <<endl<< endl;
  report << "Observed_SE:   " << obs_srv1_se  <<endl<< endl;
  report << "Survey Biomass " <<endl;
  report << "Year:     " << yrs_srv2  <<endl;
  report << "Predicted:   " << pred_srv2  <<endl;
  report << "Observed:   " << obs_srv2_biom  <<endl<< endl;
  report << "Observed_SE:   " << obs_srv2_se  <<endl<< endl;
  report << "Age  "<<aa<< endl;
  report << "Weight "<< wt << endl;
  report << "Maturity "<<p_mature<<endl<<endl;; 
  report << "Year " << yy<< endl;
  report << "Fully_selected_F "<< Fmort*max(fish_sel) <<endl<<endl;;
  report << "Year " << yy<< endl;
  report << "SpBiom "<< spawn_biom <<endl;
  report << "Tot_biom "<< tot_biom   <<endl;
    // report << tot_biom.sd <<endl;
  report << "Age  "<<aa<< endl;
  report << "Fishery_Selectivity " << fish_sel / max(fish_sel) <<endl;
  report << "TWL Survey_Selectivity " << srv1_sel / max(srv1_sel) <<endl<<endl;
  report << "LL Survey_Selectivity " << srv2_sel / max(srv2_sel) <<endl<<endl;
  report << "F35 F40 F50 "<<endl;
  report <<  F35 << " "<< F40 <<" "<<  F50 <<endl <<endl <<endl <<endl <<endl <<endl <<endl <<endl;
  report << "Wts_and_Likelihoods  (Data-Like: " <<Like<< endl;
  report << wt_ssqcatch <<" "<<ssqcatch     <<" " ; report << "SSQ_Catch_Likelihood"                << endl;
  report << wt_cpue     <<" "<<cpue_like    <<" " ; report << "Fishery_CPUE_Likelihood"             << endl;
  report << wt_srv1     <<" "<<surv_like(1) <<" " ; report << "TWL Survey_Abundance_Index_Likelihood"   << endl;
  report << wt_srv2     <<" "<<surv_like(2) <<" " ; report << "LL Survey_Abundance_Index_Likelihood"   << endl;
  report << wt_fish_age <<" "<<age_like(1)  <<" " ; report << "Fishery_Age_Composition_Likelihood"  << endl;
  report << wt_srv1_age <<" "<<age_like(2)  <<" " ; report << "Survey_Age_Composition_Likelihood"   << endl;
  report << wt_fish_size<<" "<<age_like(3)  <<" " ; report << "Fishery_Size_Composition_Likelihood" << endl;
  report << wt_srv1_size<<" "<<age_like(4)  <<" " ; report << "TWL Survey_Size_Composition_Likelihood"  << endl;
  report << wt_srv2_size<<" "<<age_like(5)  <<" " ; report << "LL Survey_Size_Composition_Likelihood"  << endl;
  report << wt_rec_var  <<" "<<rec_like     <<" " ; report << "Recruitment_Deviations_Likelihood"   << endl;
  report << wt_sel_reg_fish <<" "<<sel_like(1)      <<" " ; report << "Fish_sel_Regularity_Penalty "<< endl  ;
  report << wt_sel_reg_srv1 <<" "<<sel_like(2)      <<" " ; report << "Surv_sel_Regularity_Penalty "<< endl  ;
  report << wt_sel_reg_srv2 <<" "<<sel_like(3)      <<" " ; report << "LL Surv_sel_Regularity_Penalty "<< endl  ;
  report << wt_sel_dome_fish<<" "<<sel_like(4)      <<" " ; report << "Fish_Sel_Domeshapedness_Penalty "<<endl  ;
  report << wt_sel_dome_srv1<<" "<<sel_like(5)      <<" " ; report << "Surv_Sel_Domeshapedness_Penalty "<<endl  ;
  report << wt_sel_dome_srv2<<" "<<sel_like(6)      <<" " ; report << "LL Surv_Sel_Domeshapedness_Penalty "<<endl  ;
  report << wt_avg_sel   <<" "<<avg_sel_penalty  <<" " ; report << "Average_Selectivity_Condition_Penalty "<<endl  ;
  report << wt_fmort_reg     <<" "<<F_mort_regularity<<" " ; report << "Fishing_Mortality_Regularity_Penalty" << endl;
  report << " "<<priors(1)  <<" " ; report << "priors sigr"     <<endl;
  report << " "<<priors(2)  <<" " ; report << "priors q TWL survey (1)" <<endl;
  report << " "<<priors(5)  <<" " ; report << "priors q LL survey (2)" <<endl;
  report << " "<<priors(4)  <<" " ; report << "priors M"<<endl;
  report << " "<<Like_L  <<" " ; report << "L mat like"<<endl;
  report << " "<<obj_fun    <<" " ; report << "obj_fun"         <<endl;
  report << " "<<Like       <<" " ; report << "data likelihood" <<endl;//(2*square(sigr))+ size_count(log_rec_dev)*log(sigr)<<endl;
  report << "SigmaR: "<<sigr<< " Nat_Mort: "<<natmort<<" Survey_q: "<<q_srv1<<" LL Survey_q: "<<q_srv2<<" Spawning Per Recruit "<< " "<<SBF40<<" "<<SB0<<" Virgin SPR "<<endl;
  report << "Stock-recruitment, type: "<<SrType<<" 1=Ricker, 2=B-Holt, 3=Mean"<<endl;
  report << "YearClass SSB SR_Pred R_Est "<<endl;
  report<< styr_rec<<" "<<Sp_Biom(styr_rec-recage)<<" "<<srm_rec(styr_rec)<<" "<<sam_rec(styr_rec)<<" 0 "<< endl;
  for (i=styr_rec+1;i<=endyr;i++)
    report<< i-recage <<" "<<Sp_Biom(i-recage)<<" "<<srm_rec(i)<<" "<<sam_rec(i)<<" "<< wt_rec_var * (log_rec_dev(i)*log_rec_dev(i)/(2.*square(sigr)) + log(sigr))<<endl;
  report<<"Projection outputs"<<endl;
  report << "N_at_age projected "<<endl<<N_proj<<endl<<" spawn_bio projected"<<endl<<spawn_biom_proj<<endl;
}

void model_parameters::write_sumreport(void)
{
  ofstream& evalout= *pad_evalout;
    ofstream sumreport("report.rep");
  sumreport<<"****Executive Summary Material*****"<<endl;
  sumreport<<"     Model name"     <<endl;
  sumreport<<model_name<<endl;
  sumreport<<"     .dat file"     <<endl;
  sumreport<<data_file<<endl;
  sumreport<<"     Number parameters estimated"     <<endl;
  sumreport<<initial_params::nvarcalc()<<endl;
  sumreport<<"     TotalBiomass for "<<endyr+1<<endl;
  sumreport<<N_proj(endyr+1)*wt<<endl;
  sumreport<<"     TotalBiomass for "<<endyr+2     <<endl;
  sumreport<<N_proj(endyr+2)*wt<<endl;
  sumreport<<"     Female_Spawning Biomass for "<<endyr+1     <<endl;
  sumreport<<spawn_biom_proj(endyr+1)<<endl;
  sumreport<<"     Female_Spawning_Biomass for "<<endyr+2     <<endl;
  sumreport<<spawn_biom_proj(endyr+2)<<endl;
  sumreport<<"     B_zero"     <<endl;
  sumreport<<SB0*mean(pred_rec(1979,endyr-recage))<<endl;
  sumreport<<"     B_40"     <<endl;
  sumreport<<B40<<endl;
  sumreport<<"     B_35"     <<endl;
  sumreport<<SBF35*mean(pred_rec(1979,endyr-recage))<<endl;
  sumreport<<"     F_40"     <<endl;
  sumreport<<F40<<endl;
  sumreport<<"     F_35"     <<endl;
  sumreport<<F35<<endl;
  sumreport<<"     F_ABC for "<<endyr+1     <<endl;
  sumreport<<FABC<<endl;
  sumreport<<"     F_ABC for "<<endyr+2     <<endl;
  sumreport<<FABC2<<endl;
  sumreport<<"     ABC for "<<endyr+1     <<endl;
  sumreport<<pred_catch_proj(endyr+1)<<endl;
  sumreport<<"     ABC for "<<endyr+2     <<endl;
  sumreport<<pred_catch_proj(endyr+2)<<endl;
  sumreport<<"     F_OFL for "<<endyr+1     <<endl;
  sumreport<<FOFL<<endl;
  sumreport<<"     F_OFL for "<<endyr+2     <<endl;
  sumreport<<FOFL2<<endl;
  sumreport<<"     OFL for "<<endyr+1     <<endl;
  sumreport<<OFL<<endl; 
  sumreport<<"     OFL for "<<endyr+2     <<endl;
  sumreport<<pred_catch_proj_OFL(endyr+2)<<endl; 
  sumreport<<"     Total likelihood"     <<endl;
  sumreport<<obj_fun<<endl;
  sumreport<<"     Data likelihood"     <<endl;
  sumreport<<Like<<endl<<endl;
  sumreport<<" ************   Some more parameter estimates and their SDs ************"<<endl;
  if(last_phase()) {
    // add standard deviation data types    
  sumreport<<"   q_trawl   "<<endl;
  sumreport<<q_srv1<<" "<<q_srv1.sd<<endl;
  sumreport<<"   q_longline  "<<endl;
  sumreport<<q_srv2<<" "<<q2.sd<<endl;
  sumreport<<"   nat_mort  "<<endl;
  sumreport<<natmort<<" "<<nattymort.sd<<endl;
  sumreport<<"  sigr   "<<endl;  
  sumreport<<sigr<<" "<<cigar.sd<<endl;  
  sumreport<<"   log_mean_rec"<<endl;
  sumreport<<log_mean_rec<<" "<<LMR.sd<<endl;
  sumreport<<"   F_40"<<endl;
  sumreport<<F40<<" "<<F40.sd<<endl;
  sumreport<<"    tot_biom"<<endl;
  sumreport<<tot_biom_proj(endyr+1)<<" "<<tot_biom_proj.sd(endyr+1)<<endl;
  sumreport<<"   spawn_biom"<<endl;
  sumreport<<spawn_biom_proj(endyr+1)<<" "<<spawn_biom_proj.sd(endyr+1)<<endl;
  sumreport<<"    B40"<<endl;
  sumreport<<B40<<" "<<B40.sd<<endl;
  sumreport<<"   ABC"<<endl;
  sumreport<<ABC<<" "<<ABC.sd<<endl<<endl;
 }
  // end exec
  sumreport<<"************Rest of data output **********"<<endl<<endl;
  sumreport << "Year "<< yy <<endl;
  sumreport << "Pred_Catch "<< pred_catch_early<<" Pred_catch_later "<<pred_catch_later <<endl;
  sumreport << "Obs_Catch "<< obs_catch_early<<" Obs_Catch_Later "<<obs_catch_later <<endl;
  sumreport << "Catch_at_age "<<aa <<endl;
  for (i=styr;i<=endyr;i++) sumreport << i<<" "<<catage(i) <<endl; sumreport<<endl;
  sumreport << "Numbers "<<aa <<endl;
  for (i=styr;i<=endyr;i++) sumreport << i<<" "<<natage(i) <<endl; sumreport<<endl;
  if (nyrs_cpue>0) 
  {
    sumreport <<"Years_CPUE: "     <<yrs_cpue <<endl; 
    sumreport <<"Predicted_CPUE: " <<pred_cpue<<endl; 
    sumreport <<"observed_CPUE: "  <<obs_cpue <<endl; 
    sumreport <<"q_CPUE: "         <<q_cpue   <<endl; 
  }
  sumreport << "Obs_P_fish_age"<<aa <<endl;
  for (i=1;i<=nyrs_fish_age;i++) sumreport << yrs_fish_age(i)<<" "<<oac_fish(i) 
      <<" eff_N "<<(1-eac_fish(i))*eac_fish(i)/norm2(oac_fish(i)-eac_fish(i))  <<" N "<<nmulti_fish_age(i)
      <<" SDNR "<< sdnr(eac_fish(i),oac_fish(i),double(nmulti_fish_age(i)))<<endl; sumreport<<endl;
  sumreport << "Pred_P_fish_age"<<aa <<endl;
  for (i=1;i<=nyrs_fish_age;i++) sumreport << yrs_fish_age(i)<<" "<<eac_fish(i) <<endl; sumreport<<endl;
  sumreport << "Obs_P_fish_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish_size;i++) sumreport << yrs_fish_size(i)<<" "<<osc_fish(i) 
      <<" eff_N "<<(1-esc_fish(i))*esc_fish(i)/norm2(osc_fish(i)-esc_fish(i))  <<" N "<<nmulti_fish_size(i)
      <<" SDNR "<< sdnr(esc_fish(i),osc_fish(i),double(nmulti_fish_size(i)))<<endl; sumreport<<endl;
  sumreport << "Pred_P_fish_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish_size;i++) sumreport << yrs_fish_size(i)<<" "<<esc_fish(i) <<endl; sumreport<<endl;
  sumreport << "Obs_P_srv1_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv1_age;i++) sumreport << yrs_srv1_age(i)<<" "<<oac_srv1(i) 
      <<" eff_N "<<(1-eac_srv1(i))*eac_srv1(i)/norm2(oac_srv1(i)-eac_srv1(i)) <<" N "<<nmulti_srv1_age(i)
      <<" SDNR "<< sdnr(eac_srv1(i),oac_srv1(i),double(nmulti_srv1_age(i)))<<endl; sumreport<<endl;
  sumreport << "Pred_P_srv1_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv1_age;i++) sumreport << yrs_srv1_age(i)<<" "<<eac_srv1(i) <<endl; sumreport<<endl;
  sumreport << "Obs_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) sumreport << yrs_srv1_size(i)<<" "<<osc_srv1(i) 
      <<" eff_N "<<(1-esc_srv1(i))*esc_srv1(i)/norm2(osc_srv1(i)-esc_srv1(i)) <<" N "<<nmulti_srv1_size(i)
      <<" SDNR "<< sdnr(esc_srv1(i),osc_srv1(i),double(nmulti_srv1_size(i)))<<endl; sumreport<<endl;
  sumreport << "Pred_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) sumreport << yrs_srv1_size(i)<<" "<<esc_srv1(i) <<endl; sumreport<<endl;
  sumreport << "Obs_P_srv2_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) sumreport << yrs_srv2_size(i)<<" "<<osc_srv2(i) 
      <<" eff_N "<<(1-esc_srv2(i))*esc_srv2(i)/norm2(osc_srv2(i)-esc_srv2(i)) <<" N "<<nmulti_srv2_size(i)
      <<" SDNR "<< sdnr(esc_srv2(i),osc_srv2(i),double(nmulti_srv2_size(i)))<<endl; sumreport<<endl;
  sumreport << "Pred_P_srv2_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) sumreport << yrs_srv2_size(i)<<" "<<esc_srv2(i) <<endl; sumreport<<endl;
    sumreport << "Survey Biomass " <<endl;
  sumreport << "Year:     " << yrs_srv1  <<endl;
  sumreport << "Predicted:   " << pred_srv1  <<endl;
  sumreport << "Observed:   " << obs_srv1_biom  <<endl<< endl;
  sumreport << "Observed_SE:   " << obs_srv1_se  <<endl<< endl;
  sumreport << "Survey Biomass " <<endl;
  sumreport << "Year:     " << yrs_srv2  <<endl;
  sumreport << "Predicted:   " << pred_srv2  <<endl;
  sumreport << "Observed:   " << obs_srv2_biom  <<endl<< endl;
  sumreport << "Observed_SE:   " << obs_srv2_se  <<endl<< endl;
  sumreport << "Age  "<<aa<< endl;
  sumreport << "Weight "<< wt << endl;
  sumreport << "Maturity "<<p_mature<<endl<<endl;; 
  sumreport << "Year " << yy<< endl;
  sumreport << "Fully_selected_F "<< Fmort*max(fish_sel) <<endl<<endl;;
  sumreport << "Year " << yy<< endl;
  sumreport << "SpBiom "<< spawn_biom <<endl;
  sumreport << "Tot_biom "<< tot_biom   <<endl;
    // sumreport << tot_biom.sd <<endl;
  sumreport << "Age  "<<aa<< endl;
  sumreport << "Fishery_Selectivity " << fish_sel / max(fish_sel) <<endl;
  sumreport << "TWL Survey_Selectivity " << srv1_sel / max(srv1_sel) <<endl<<endl;
  sumreport << "LL Survey_Selectivity " << srv2_sel / max(srv2_sel) <<endl<<endl;
  sumreport << "F35 F40 F50 "<<endl;
  sumreport <<  F35 << " "<< F40 <<" "<<  F50 <<endl <<endl <<endl <<endl <<endl <<endl <<endl <<endl;
  sumreport << "Wts_and_Likelihoods  (Data-Like: " <<Like<< endl;
  sumreport << wt_ssqcatch <<" "<<ssqcatch     <<" " ; sumreport << "SSQ_Catch_Likelihood"                << endl;
  sumreport << wt_cpue     <<" "<<cpue_like    <<" " ; sumreport << "Fishery_CPUE_Likelihood"             << endl;
  sumreport << wt_srv1     <<" "<<surv_like(1) <<" " ; sumreport << "TWL Survey_Abundance_Index_Likelihood"   << endl;
  sumreport << wt_srv2     <<" "<<surv_like(2) <<" " ; sumreport << "LL Survey_Abundance_Index_Likelihood"   << endl;
  sumreport << wt_fish_age <<" "<<age_like(1)  <<" " ; sumreport << "Fishery_Age_Composition_Likelihood"  << endl;
  sumreport << wt_srv1_age <<" "<<age_like(2)  <<" " ; sumreport << "Survey_Age_Composition_Likelihood"   << endl;
  sumreport << wt_fish_size<<" "<<age_like(3)  <<" " ; sumreport << "Fishery_Size_Composition_Likelihood" << endl;
  sumreport << wt_srv1_size<<" "<<age_like(4)  <<" " ; sumreport << "TWL Survey_Size_Composition_Likelihood"  << endl;
  sumreport << wt_srv2_size<<" "<<age_like(5)  <<" " ; sumreport << "LL Survey_Size_Composition_Likelihood"  << endl;
  sumreport << wt_rec_var  <<" "<<rec_like     <<" " ; sumreport << "Recruitment_Deviations_Likelihood"   << endl;
  sumreport << wt_sel_reg_fish <<" "<<sel_like(1)      <<" " ; sumreport << "Fish_sel_Regularity_Penalty "<< endl  ;
  sumreport << wt_sel_reg_srv1 <<" "<<sel_like(2)      <<" " ; sumreport << "Surv_sel_Regularity_Penalty "<< endl  ;
  sumreport << wt_sel_reg_srv2 <<" "<<sel_like(3)      <<" " ; sumreport << "LL Surv_sel_Regularity_Penalty "<< endl  ;
  sumreport << wt_sel_dome_fish<<" "<<sel_like(4)      <<" " ; sumreport << "Fish_Sel_Domeshapedness_Penalty "<<endl  ;
  sumreport << wt_sel_dome_srv1<<" "<<sel_like(5)      <<" " ; sumreport << "Surv_Sel_Domeshapedness_Penalty "<<endl  ;
  sumreport << wt_sel_dome_srv2<<" "<<sel_like(6)      <<" " ; sumreport << "LL Surv_Sel_Domeshapedness_Penalty "<<endl  ;
  sumreport << wt_avg_sel   <<" "<<avg_sel_penalty  <<" " ; sumreport << "Average_Selectivity_Condition_Penalty "<<endl  ;
  sumreport << wt_fmort_reg     <<" "<<F_mort_regularity<<" " ; sumreport << "Fishing_Mortality_Regularity_Penalty" << endl;
  sumreport << " "<<priors(1)  <<" " ; sumreport << "priors sigr"     <<endl;
  sumreport << " "<<priors(2)  <<" " ; sumreport << "priors q TWL survey (1)" <<endl;
  sumreport << " "<<priors(5)  <<" " ; sumreport << "priors q LL survey (2)" <<endl;
  sumreport << " "<<priors(4)  <<" " ; sumreport << "priors M"<<endl;
  sumreport << " "<<obj_fun    <<" " ; sumreport << "obj_fun"         <<endl;
  sumreport << " "<<Like       <<" " ; sumreport << "data likelihood" <<endl;//(2*square(sigr))+ size_count(log_rec_dev)*log(sigr)<<endl;
  sumreport << "SigmaR: "<<sigr<< " Nat_Mort: "<<natmort<<" Survey_q: "<<q_srv1<<" LL Survey_q: "<<q_srv2<<" Spawning Per Recruit "<< " "<<SBF40<<" "<<SB0<<" Virgin SPR "<<endl;
  sumreport << "Stock-recruitment, type: "<<SrType<<" 1=Ricker, 2=B-Holt, 3=Mean"<<endl;
  sumreport << "YearClass SSB SR_Pred R_Est "<<endl;
  sumreport<< styr_rec<<" "<<Sp_Biom(styr_rec-recage)<<" "<<srm_rec(styr_rec)<<" "<<sam_rec(styr_rec)<<" 0 "<< endl;
  for (i=styr_rec+1;i<=endyr;i++)
    sumreport<< i-recage <<" "<<Sp_Biom(i-recage)<<" "<<srm_rec(i)<<" "<<sam_rec(i)<<" "<< wt_rec_var * (log_rec_dev(i)*log_rec_dev(i)/(2.*square(sigr)) + log(sigr))<<endl;
  sumreport<<"Projection outputs"<<endl;
  sumreport << "N_at_age projected "<<endl<<N_proj<<endl<<" spawn_bio projected"<<endl<<spawn_biom_proj<<endl;
}

double model_parameters::sdnr(const dvar_vector& pred,const dvector& obs,double m)
{
  ofstream& evalout= *pad_evalout;
  RETURN_ARRAYS_INCREMENT();
  double sdnr;
  dvector pp = value(pred);
  int ntmp = -obs.indexmin()+obs.indexmax();
  sdnr = std_dev(elem_div(obs-pp,sqrt(elem_prod(pp,(1.-pp))/m)));
  RETURN_ARRAYS_DECREMENT();
  return sdnr;
}

void model_parameters::write_proj(void)
{
  ofstream& evalout= *pad_evalout;
 ofstream newproj("proj.dat");
 newproj <<"#Species name here:"<<endl;
 newproj <<model_name+"_"+data_file<<endl;
 newproj <<"#SSL Species?"<<endl;
 newproj <<"0"<<endl;
 newproj <<"#Constant buffer of Dorn?"<<endl;
 newproj <<"0"<<endl;
 newproj <<"#Number of fisheries?"<<endl;
 newproj <<"1"<<endl;
 newproj <<"#Number of sexes?"<<endl;
 newproj <<"1"<<endl;
 newproj <<"#5year_Average_F(endyr-4,endyr_as_estimated_by_ADmodel)"<<endl;
 newproj << mean(Fmort(endyr-4,endyr))<<endl;
 newproj <<"#_Author_F_as_fraction_F_40%"<<endl;
 newproj <<"1"<<endl;
 newproj <<"#ABC SPR" <<endl;
 newproj <<"0.4"<<endl;
 newproj <<"#MSY SPR" <<endl;
 newproj <<"0.35"<<endl;
 newproj <<"#_Spawn_month"<<endl;
 newproj << spawn_fract*12+1<<endl;
 newproj <<"#_Number_of_ages"<<endl;
 newproj <<nages_M<<endl;
 newproj <<"#_F_ratio(must_sum_to_one_only_one_fishery)"<<endl;
 newproj <<"1"<<endl;
 for (j=1;j<=nages_M;j++) natmortv = natmort; 
 newproj <<"#_Natural_Mortality" << aa << endl;
 newproj <<natmortv<<endl;
 newproj <<"#_Maturity_divided_by_2(projection_program_uses_to_get_female_spawning_biomass_if_divide_by_2"<<aa<<endl<<p_mature<< endl;
 newproj <<"#_Wt_at_age_spawners"<<aa<<endl<<wt<< endl;
 newproj <<"#_Wt_at_age_fishery" <<aa<<endl<<wt<< endl;
 newproj <<"#_Selectivity_fishery_scaled_to_max_at_one"<<aa<<endl<<fish_sel/max(fish_sel)<< endl;
 newproj <<"#_Numbers_at_age_end_year"<<aa<<endl<<natage(endyr)<< endl;
 newproj <<"#_N_recruitment_years"<<endl<<endyr-recage-1979+1<< endl;
 newproj <<"#_Recruitment_start_at_1977_yearclass=1979_for_age_2_recruits"<<yy(1979,endyr-recage)<<endl<<pred_rec(1979,endyr-recage)<< endl;
 newproj <<"#_Spawners per recruitment (starting at 1977)"<<endl<<spawn_biom(1977,endyr-recage)/1000<< endl;
 newproj.close();
}

void model_parameters::write_r_report(void)
{
  ofstream& evalout= *pad_evalout;
   ofstream rreport("rtem.rep");
  rreport<<"$model_rownames"<<endl;
  rreport<<"model_name"<<endl;
  rreport<<"dat_file"<<endl;
  rreport<<"number_parameters_estimated"<<endl;
  rreport<<"SR_type_1Ricker_2BH_3Mean"<<endl;
  rreport<<"$model"<<endl;
  rreport<<model_name<<endl;
  rreport<<data_file<<endl;
  rreport<<initial_params::nvarcalc()<<endl;
  rreport<<SrType<<endl;
  rreport<<"$summary_rownames"<<endl;
  rreport<<"Total_Biomass_one_year"<<endl;
  rreport<<"Total_Biomass_two_year"<<endl;
  rreport<<"Female_Spawning_Biomass_one_year"<<endl;
  rreport<<"Female_Spawning_Biomass_two_year"<<endl;
  rreport<<"B_zero"<<endl;
  rreport<<"B_40"<<endl;
  rreport<<"B_35"<<endl;
  rreport<<"M"<<endl;
  rreport<<"F_40"<<endl;
  rreport<<"F_35"<<endl;
  rreport<<"F_ABC_one_year"<<endl;
  rreport<<"F_ABC_two_year"<<endl;
  rreport<<"ABC_one_year"<<endl;
  rreport<<"ABC_two_year"<<endl;
  rreport<<"F_OFL_one_year"<<endl;
  rreport<<"F_OFL_two_year"<<endl;
  rreport<<"OFL_one_year"<<endl;
  rreport<<"OFL_two_year"<<endl;
  rreport<<"$summary"<<endl;
  rreport<<N_proj(endyr+1)*wt<<endl;
  rreport<<N_proj(endyr+2)*wt<<endl;
  rreport<<spawn_biom_proj(endyr+1)<<endl;
  rreport<<spawn_biom_proj(endyr+2)<<endl;
  rreport<<SB0*mean(pred_rec(1979,endyr-recage))<<endl;
  rreport<<B40<<endl;
  rreport<<SBF35*mean(pred_rec(1979,endyr-recage))<<endl;
  rreport<<natmort<<endl;
  rreport<<F40<<endl;
  rreport<<F35<<endl;
  rreport<<FABC<<endl;
  rreport<<FABC2<<endl;
  rreport<<pred_catch_proj(endyr+1)<<endl;
  rreport<<pred_catch_proj(endyr+2)<<endl;
  rreport<<FOFL<<endl;
  rreport<<FOFL2<<endl;
  rreport<<OFL<<endl;
  rreport<<pred_catch_proj_OFL(endyr+2)<<endl;
  rreport<<"$year"<<endl;
  rreport<<yy<<endl;
  rreport<<"$age"<<endl;
  rreport<<aa<<endl;
  rreport<<"$size"<<endl;
  rreport<<len_bin_labels<<endl;
  rreport<<"$weight"<<endl;
  rreport<<wt<< endl;
  rreport<<"$maturity"<<endl;
  rreport<<p_mature<<endl; 
  rreport<<"$tseries_colnames"<<endl;
  rreport<<"year"<<endl;
  rreport<<"obs_catch"<<endl;
  rreport<<"pred_catch"<<endl;
  rreport<<"ful_sel_F"<<endl;
  rreport<<"tot_biom"<<endl;
  rreport<<"spawn_biom"<<endl;
  rreport<<"$tseries"<<endl;
  rreport<<yy<<endl;
  rreport<<obs_catch_early<<" "<<obs_catch_later<<endl;
  rreport<<pred_catch_early<<" "<<pred_catch_later<<endl;
  rreport<<Fmort*max(fish_sel)<<endl;;
  rreport<<tot_biom<<endl;
  rreport<<spawn_biom<<endl;
  rreport<<"$catage"<<endl;
  for (i=styr;i<=endyr;i++) rreport << i<<" "<<catage(i) <<endl; rreport<<endl;
  rreport<<"$natage"<<endl;
  for (i=styr;i<=endyr;i++) rreport << i<<" "<<natage(i) <<endl; rreport<<endl;
  rreport<<"$twl_rownames"<<endl;
  rreport<<"Year"<<endl;
  rreport<<"Observed"<<endl;
  rreport<<"Observed_SE"<<endl;
  rreport<<"Observed_LCI"<<endl;
  rreport<<"Observed_UCI"<<endl;
  rreport<<"Predicted"<<endl;
  rreport<<"$twl"<<endl;
  rreport<<yrs_srv1<<endl;
  rreport<<obs_srv1_biom<<endl;
  rreport<<obs_srv1_se<<endl;
  rreport<<obs_srv1_lci<<endl;
  rreport<<obs_srv1_uci<<endl;
  rreport<<pred_srv1<<endl;
  rreport<<"$ll_rownames"<<endl;
  rreport<<"Year"<<endl;
  rreport<<"Observed"<<endl;
  rreport<<"Observed_SE"<<endl;
  rreport<<"Observed_LCI"<<endl;
  rreport<<"Observed_UCI"<<endl;
  rreport<<"Predicted"<<endl;
  rreport<<"$ll"<<endl;
  rreport<<yrs_srv2<<endl;
  rreport<<obs_srv2_biom<<endl;
  rreport<<obs_srv2_se<<endl;
  rreport<<obs_srv2_lci<<endl;
  rreport<<obs_srv2_uci<<endl;
  rreport<<pred_srv2<<endl;
  rreport<<"$selectivity_rownames"<<endl;
  rreport<<"fish_sel"<<endl;
  rreport<<"twl_sel"<<endl;
  rreport<<"ll_sel"<<endl;
  rreport<<"$selectivity"<<endl;
  rreport<<fish_sel/max(fish_sel)<<endl;
  rreport<<srv1_sel/max(srv1_sel)<<endl;
  rreport<<srv2_sel/max(srv2_sel)<<endl;
  if (nyrs_cpue>0) 
  {
    rreport<<"$cpue_fish_rownames"<<endl;
	rreport<<"Year"<<endl; 
    rreport<<"Observed"<<endl; 
    rreport<<"Predicted"<<endl; 
    rreport<<"q_cpue"<<endl; 
    rreport<<"$cpue_fish"<<endl;
    rreport<<yrs_cpue<<endl; 
    rreport<<obs_cpue<<endl; 
    rreport<<pred_cpue<<endl; 
    rreport<<q_cpue<<endl; 
  }
  rreport<<"$yrs_fish_age"<<endl;
  for (i=1;i<=nyrs_fish_age;i++) rreport << yrs_fish_age(i)<<endl; rreport<<endl;
  rreport<<"$oac_fish_sample_colnames"<<endl;
  rreport<<"N"<<endl;
  rreport<<"eff_n"<<endl;
  rreport<<"sdnr"<<endl;
  rreport<<"$oac_fish_sample"<<endl;
  for (i=1;i<=nyrs_fish_age;i++) rreport << nmulti_fish_age(i)
      <<" "<<(1-eac_fish(i))*eac_fish(i)/norm2(oac_fish(i)-eac_fish(i))
  	  <<" "<<sdnr(eac_fish(i),oac_fish(i),double(nmulti_fish_age(i)))<<endl; rreport<<endl;
  rreport<<"$oac_fish"<<endl;
  for (i=1;i<=nyrs_fish_age;i++) rreport << oac_fish(i) <<endl; rreport<<endl;
  rreport<<"$eac_fish"<<endl;
  for (i=1;i<=nyrs_fish_age;i++) rreport << eac_fish(i) <<endl; rreport<<endl;
  rreport<<"$yrs_fish_size"<<endl;
  for (i=1;i<=nyrs_fish_size;i++) rreport << yrs_fish_size(i)<<endl; rreport<<endl;
  rreport<<"$osc_fish_sample_colnames"<<endl;
  rreport<<"N"<<endl;
  rreport<<"eff_n"<<endl;
  rreport<<"sdnr"<<endl;
  rreport<<"$osc_fish_sample"<<endl;
  for (i=1;i<=nyrs_fish_size;i++) rreport << nmulti_fish_size(i)
      <<" "<<(1-esc_fish(i))*esc_fish(i)/norm2(osc_fish(i)-esc_fish(i))
      <<" "<< sdnr(esc_fish(i),osc_fish(i),double(nmulti_fish_size(i)))<<endl; rreport<<endl;
  rreport<<"$osc_fish"<<endl;
  for (i=1;i<=nyrs_fish_size;i++) rreport << osc_fish(i) <<endl; rreport<<endl;
  rreport<<"$esc_fish"<<endl;
  for (i=1;i<=nyrs_fish_size;i++) rreport << esc_fish(i) <<endl; rreport<<endl;
  rreport<<"$yrs_srv1_age"<<endl;
  for (i=1;i<=nyrs_srv1_age;i++) rreport << yrs_srv1_age(i)<<endl; rreport<<endl;
  rreport<<"$oac_srv1_sample_colnames"<<endl;
  rreport<<"N"<<endl;
  rreport<<"eff_n"<<endl;
  rreport<<"sdnr"<<endl;
  rreport<<"$oac_srv1_sample"<<endl;
  for (i=1;i<=nyrs_srv1_age;i++) rreport << nmulti_srv1_age(i)
      <<" "<<(1-eac_srv1(i))*eac_srv1(i)/norm2(oac_srv1(i)-eac_srv1(i))
      <<" "<< sdnr(eac_srv1(i),oac_srv1(i),double(nmulti_srv1_age(i)))<<endl; rreport<<endl;
  rreport<<"$oac_srv1"<<endl;
  for (i=1;i<=nyrs_srv1_age;i++) rreport << oac_srv1(i) <<endl; rreport<<endl;
  rreport<<"$eac_srv1"<<endl;
  for (i=1;i<=nyrs_srv1_age;i++) rreport << eac_srv1(i) <<endl; rreport<<endl;
  rreport<<"$yrs_srv1_size"<<endl;
  for (i=1;i<=nyrs_srv1_size;i++) rreport << yrs_srv1_size(i)<<endl; rreport<<endl;
  rreport<<"$osc_srv1_sample_colnames"<<endl;
  rreport<<"N"<<endl;
  rreport<<"eff_n"<<endl;
  rreport<<"sdnr"<<endl;
  rreport<<"$osc_srv1_sample"<<endl;  
  for (i=1;i<=nyrs_srv1_size;i++) rreport << nmulti_srv1_size(i)
      <<" "<<(1-esc_srv1(i))*esc_srv1(i)/norm2(osc_srv1(i)-esc_srv1(i)) 
      <<" "<< sdnr(esc_srv1(i),osc_srv1(i),double(nmulti_srv1_size(i)))<<endl; rreport<<endl;
  rreport<<"$osc_srv1"<<endl;
  for (i=1;i<=nyrs_srv1_size;i++) rreport << osc_srv1(i) <<endl; rreport<<endl;
  rreport<<"$esc_srv1"<<endl;
  for (i=1;i<=nyrs_srv1_size;i++) rreport << esc_srv1(i) <<endl; rreport<<endl;
  rreport<<"$yrs_srv2_size"<<endl;
  for (i=1;i<=nyrs_srv2_size;i++) rreport << yrs_srv2_size(i)<<endl; rreport<<endl;
  rreport<<"$osc_srv2_sample_colnames"<<endl;
  rreport<<"N"<<endl;
  rreport<<"eff_n"<<endl;
  rreport<<"sdnr"<<endl;
  rreport<<"$osc_srv2_sample"<<endl;  
  for (i=1;i<=nyrs_srv2_size;i++) rreport << nmulti_srv2_size(i)
      <<" "<<(1-esc_srv2(i))*esc_srv2(i)/norm2(osc_srv2(i)-esc_srv2(i)) 
      <<" "<< sdnr(esc_srv2(i),osc_srv2(i),double(nmulti_srv2_size(i)))<<endl; rreport<<endl;
  rreport<<"$osc_srv2"<<endl;
  for (i=1;i<=nyrs_srv2_size;i++) rreport << osc_srv2(i) <<endl; rreport<<endl;
  rreport<<"$esc_srv2"<<endl;
  for (i=1;i<=nyrs_srv2_size;i++) rreport << esc_srv2(i) <<endl; rreport<<endl;
  rreport<<"$like_rownames"<<endl;
  rreport<<"SSQ_Catch_Likelihood"<<endl;
  rreport<<"Fishery_CPUE_Likelihood"<<endl;
  rreport<<"TWL_Survey_Abundance_Index_Likelihood"<<endl;
  rreport<<"LL_Survey_Abundance_Index_Likelihood"<<endl;
  rreport<<"Fishery_Age_Composition_Likelihood"<<endl;
  rreport<<"Survey_Age_Composition_Likelihood"<<endl;
  rreport<<"Fishery_Size_Composition_Likelihood"<<endl;
  rreport<<"TWL_Survey_Size_Composition_Likelihood"<<endl;
  rreport<<"LL_Survey_Size_Composition_Likelihood"<<endl;
  rreport<<"$like"<<endl;
  rreport<<ssqcatch<<" "<<wt_ssqcatch<<endl; 
  rreport<<cpue_like<<" "<<wt_cpue<<endl; 
  rreport<<surv_like(1)<<" "<<wt_srv1<<endl; 
  rreport<<surv_like(2)<<" "<<wt_srv2<<endl; 
  rreport<<age_like(1)<<" "<<wt_fish_age<<endl; 
  rreport<<age_like(2)<<" "<<wt_srv1_age<<endl; 
  rreport<<age_like(3)<<" "<<wt_fish_size<<endl; 
  rreport<<age_like(4)<<" "<<wt_srv1_size<<endl; 
  rreport<<age_like(5)<<" "<<wt_srv2_size<<endl; 
  rreport<<"$pen_rownames"<<endl;
  rreport<<"Recruitment_Deviations"<<endl;
  rreport<<"Fish_sel_Regularity_Penalty"<<endl;
  rreport<<"Surv_sel_Regularity_Penalty"<<endl;
  rreport<<"LL_Surv_sel_Regularity_Penalty"<<endl;
  rreport<<"Fish_Sel_Domeshapedness_Penalty"<<endl;
  rreport<<"Surv_Sel_Domeshapedness_Penalty"<<endl;
  rreport<<"LL_Surv_Sel_Domeshapedness_Penalty"<<endl;
  rreport<<"Average_Selectivity_Condition_Penalty"<<endl;
  rreport<<"Fishing_Mortality_Regularity_Penalty"<<endl;
  rreport<<"$pen"<<endl;
  rreport<<rec_like<<" "<<wt_rec_var<<endl; 
  rreport<<sel_like(1)<<" "<<wt_sel_reg_fish<<endl;
  rreport<<sel_like(2)<<" "<<wt_sel_reg_srv1<<endl;
  rreport<<sel_like(3)<<" "<<wt_sel_reg_srv2<<endl;
  rreport<<sel_like(4)<<" "<<wt_sel_dome_fish<<endl;
  rreport<<sel_like(5)<<" "<<wt_sel_dome_srv1<<endl;
  rreport<<sel_like(6)<<" "<<wt_sel_dome_srv2<<endl;
  rreport<<avg_sel_penalty<<" "<<wt_avg_sel<<endl;
  rreport<<F_mort_regularity<<" "<<wt_fmort_reg<<endl;
  rreport<<"$prior_rownames"<<endl;
  rreport<<"priors_sigr"<<endl;
  rreport<<"priors_q_TWL_survey"<<endl;
  rreport<<"priors_q_LL_survey"<<endl;
  rreport<<"priors_M"<<endl;
  rreport<<"$prior"<<endl;
  rreport<<priors(1)<<endl; 
  rreport<<priors(2)<<endl; 
  rreport<<priors(5)<<endl; 
  rreport<<priors(4)<<endl; 
  rreport<<"$obj_rownames"<<endl;
  rreport<<"Data_Likelihood"<<endl;
  rreport<<"Penalty_Total"<<endl;
  rreport<<"Objective_Function"<<endl;
  rreport<<"$obj"<<endl;
  rreport<<Like<<endl; 
  rreport<<obj_fun - Like<<endl; 
  rreport<<obj_fun<<endl; 
  rreport<<"$parm_est_rownames"<<endl;
  rreport<<"q_trawl"<<endl;
  rreport<<"q_longline"<<endl;
  rreport<<"nat_mort"<<endl;
  rreport<<"sigr"<<endl;  
  rreport<<"log_mean_rec"<<endl;
  rreport<<"F_40"<<endl;
  rreport<<"tot_biom"<<endl;
  rreport<<"spawn_biom"<<endl;
  rreport<<"B40"<<endl;
  rreport<<"ABC"<<endl;
  // add standard deviation data types    
  rreport<<"$parm_est"<<endl;
  rreport<<q_srv1<<" "<<q_srv1.sd<<endl;
  rreport<<q_srv2<<" "<<q2.sd<<endl;
  rreport<<natmort<<" "<<nattymort.sd<<endl;
  rreport<<sigr<<" "<<cigar.sd<<endl;  
  rreport<<log_mean_rec<<" "<<LMR.sd<<endl;
  rreport<<F40<<" "<<F40.sd<<endl;
  rreport<<tot_biom_proj(endyr+1)<<" "<<tot_biom_proj.sd(endyr+1)<<endl;
  rreport<<spawn_biom_proj(endyr+1)<<" "<<spawn_biom_proj.sd(endyr+1)<<endl;
  rreport<<B40<<" "<<B40.sd<<endl;
  rreport<<ABC<<" "<<ABC.sd<<endl;
  // add size_age
  for (i=1;i<=n_sizeage_mat;i++) {
	   rreport<<"$sizeage"<<i<<endl;
	   rreport<<sizeage(i)<<endl;			//size comp #1
   }
  for (i=1;i<=n_ageage_mat;i++) {
	   rreport<<"$ageage"<<i<<endl;
	   rreport<<ageage(i)<<endl;			//size comp #1
   }
   rreport.close();
}

void model_parameters::final_calcs()
{
  write_r_report();
  write_sumreport();
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{
  delete pad_evalout;
  pad_evalout = NULL;
}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  gradient_structure::set_MAX_NVAR_OFFSET(1000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(100000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(500);
  arrmblsize=10000000;
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
