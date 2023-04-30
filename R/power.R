#' Power calculation for specific design, using simulated data derived from fixed parameters
#'
#' @param method Which model to use,can be Bliss, HSA or RA
#' @param type Which metric to use, can be AUC or TGI
#' @param n number of animals per group
#' @param tv0 initial average tumor volume
#' @param cv_tv0 standard deviation of initial tumor volume, also assuming standard deviation of tumor volume is sd0
#' @param vehicle_double_time days to take for double initial tumor size
#' @param t length of study
#' @param tgi_a TGI of drug a at day 21
#' @param tgi_b TGI of durg b at day 21
#' @param tgi_c TGI of combo at day 21
#' @param n_core number of cores used for parallel computing
#' @param cv_r ratio between standard deviation of tumor growth rate and tumor growth rate for vehicle, reflects heterogeneity of tumor growth
#' @param eu_tv mouse is killed if TV exceeds eu_tv
#' @param kit_a kick in time for drug a, default 0
#' @param kit_b kick in time for drug b, default 0
#' @param kit_c kick in time for drug combo, default 0
#' @param rp_a resistant proportion of drug a
#' @param rp_b resistant proportion of drug b
#' @param rp_c resistant proportion of drug combo
#' @param sigma random error, default to 1
#' @param n_sim number of simulations
#' @param irt_a day of resistance for drug a (if not 0, induced resistance)
#' @param irt_b day of resistance for drug b
#' @param irt_c day of resistance for drug c
#'
#' @return power of the experiment
#'
#' @export
#'
#' @examples
#' power_calc()
#'
power_calc=function(method="Bliss",type='TGI',n_sim=100,n=5,tv0=50,cv_tv0=0.2,vehicle_double_time=3,t=21,cv_r=0.05,eu_tv=3000,tgi_a=0.5,tgi_b=0.4,tgi_c=0.8,kit_a=0,kit_b=0,kit_c=0,rp_a=0,rp_b=0,rp_c=0,irt_a=0,irt_b=0,irt_c=0,sigma=1,n_core=1){
  #Power calculation using simulated dataset, generate 100 simulated dataset, calculate p-value using different methods
  out=list()
  n_sim=n_sim #generate 100 simulated datasets

  max_core=parallel::detectCores()-1
  n_core=min(max_core,n_core)
  future::plan(future::multisession, workers = n_core)
  out <- suppressMessages(furrr::future_map(1:n_sim, .f = function(x){
    tv=simu_TV(n,tv0,cv_tv0,vehicle_double_time,t,cv_r,eu_tv,tgi_a,tgi_b,tgi_c,kit_a,kit_b,kit_c,rp_a,rp_b,rp_c,sigma)
    t_vec=round(seq(0,t,3.5)) #define vector of days to collect TV
    if(t>max(t_vec)) t_vec=c(t_vec,t)
    sel_day=tv %>% filter(Group=='Group 1') %>% count(Day) %>% mutate(prop=n/max(n)) %>% filter(prop>=0.8) %>% pull(Day) %>% max()

    result=switch (type,
      'TGI' = tv %>% getTGI(sel_day=sel_day) %>% TGI_synergy(method = method,display=F,save = F),
      'CombPDX_CI' = tv %>% CombPDX_CI(sel_day=sel_day,method=method),
      'AUC' = tv %>% get_mAUCr(sel_day=t) %>% AUC_synergy(t = t,method=method,boot_n=10,display=F,save = F) %>% slice(1),
      'CombPDX_GlobalCI' = tv %>% global_CI_synergy(method=method,boot_n=100,display = F,save = F),
      'lmm' = tv %>% lmm_synergy(sel_day = t) %>% select(estimate,p.value) %>% slice(n()) %>% rename(p.val=p.value)
    )
    result
  }))
  #######################Deprecated code###############################################
  #cl=parallel::makeCluster(n_core)
  #parallel::clusterEvalQ(cl,library(SynergyInVivo))
  #parallel::clusterExport(cl,list('expand_tv'))
  #sim_tv_lst=parallel::clusterApply(cl,1:n_sim,simu_TV,n=n,tv0=tv0,sd0=sd0,vehicle_double_time=vehicle_double_time,
  #                               t=t,sigma_r=sigma_r,tgi_a=tgi_a,tgi_b=tgi_b,tgi_c=tgi_c)
  #parallel::stopCluster(cl)
  #getS=NULL
  #if(type=='TGI'){
   # getS=function(tv,t){
  #     tv %>% getTGI(sel_day=t) %>% TGI_synergy(save = F)
  #  }
  #  cl=parallel::makeCluster(n_core)
  #  parallel::clusterEvalQ(cl,library(SynergyInVivo))
   # parallel::clusterExport(cl,list('getS'),envir = environment())
   # sim_synergy=parallel::parLapply(cl,sim_tv_lst,getS,t)
  #}
  #doParallel::registerDoParallel(cl)
  #sim_tv_lst=foreach::foreach(1:n_sim,.packages = c('SynergyInVivo')) %dopar% {
  #  simu_TV(n,tv0,sd0,vehicle_double_time,t,sigma_r,tgi_a,tgi_b,tgi_c)
  #}
  #parallel::stopCluster(cl)
  bind_rows(out)
}


#' Calculate power curve for specific paramter combination)
#'
#' @param ... named parameters, same as power_calc
#'
#' @return a dataframe of powers for specific paramter combination
#' @export
#'
#' @examples
#' power_stats=sim_power(method=c('Bliss','HSA'),type=c('TGI','AUC','lmm'),n_sim=10,n=5:6)
sim_power=function(...){
  params=list(...)
  param_grid=expand.grid(params,stringsAsFactors = F) #Be careful with stringAsFactor, chars converted to factors will give erroreous results
  results=param_grid %>% purrr::pmap(power_calc,n_core=parallel::detectCores()-2)
  #results=purrr::pmap(param_grid,c)
  power=sapply(results,function(df) sum(df[['p.val']]<0.05,na.rm=T)/nrow(df))
  bind_cols(param_grid,"Power"=power)
}


