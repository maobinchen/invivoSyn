#' Generate simulated tumor volume data
#'
#' @param n number of animals per group
#' @param tv0 initial average tumor volume
#' @param vehicle_double_time days to take for double initial tumor size
#' @param t length of study
#' @param cv_r ratio between standard deviation of tumor growth rate and tumor growth rate for vehicle, reflects heterogeneity of tumor growth
#' @param tgi_a TGI of drug a at day 21, tgi is defined by 1-RTV(tr)/RTV(c)
#' @param tgi_b TGI of durg b at day 21
#' @param tgi_c TGI of combo at day 21
#' @param eu_tv mouse is killed if TV exceeds eu_tv
#' @param cv_tv0 coefficient of variation for TV0
#' @param kit_a kick in time for drug a, default 0
#' @param kit_b kick in time for drug b, default 0
#' @param kit_c kick in time for drug combo, default 0
#' @param rp_a resistant proportion of drug a
#' @param rp_b resistant proportion of drug b
#' @param rp_c resistant proportion of drug combo
#' @param sigma random error, default to 1
#' @param irt_a day of resistance for drug a (if not 0, induced resistance)
#' @param irt_b day of resistance for drug b
#' @param irt_c day of resistance for drug c
#' @param ...
#'
#' @return a tumor growth data data frame
#' @export
#'
#' @examples
#' simu_TV(n=5,tv0=50,cv_tv0=0.2,vehicle_double_time=3,t=21,cv_r=0.1,eu_tv=3000,tgi_a=0.5,tgi_b=0.4,tgi_c=0.8)
simu_TV=function(n=5,tv0=50,cv_tv0=0.2,vehicle_double_time=3,t=21,cv_r=0.1,eu_tv=3000,tgi_a=0.5,tgi_b=0.4,tgi_c=0.8,kit_a=0,kit_b=0,kit_c=0,rp_a=0,rp_b=0,rp_c=0,irt_a=0,irt_b=0,irt_c=0,sigma=1,...){
  #simulation of tumor volume data, given TV0 (assuming equal for each group), tumor doubling time for vehicle group(In(2)/T),
  #length of study(t=14 or 21), tumor growth inhibition of group a, group b and combo,variance of tumor growth rate at mouse level
  #t_vec=round(seq(0,t,3.5)) #define vector of days to collect TV
  #if(t>max(t_vec)) t_vec=c(t_vec,t)
  trs=c('Vehicle control (C)','Drug A (A)','Drug B (B)','Drug A+B (A+B)')
  treatments=factor(trs,levels=trs)
  sd0=cv_tv0*tv0
  tv0s=as.vector((sapply(1:4,function(i) rnorm(n,mean=tv0,sd=sd0))))
  k_v=log(2)/vehicle_double_time #calcul-1ate growth rate for vehicle group
  sd_k=cv_r*k_v #standard deviation of tumor growth rate
  drug_eff_vec=rep(c(0,log(1-tgi_a)/21,log(1-tgi_b)/21,log(1-tgi_c)/21),each=n)#drug effect
  #k_a=log(1-tgi_a)/t+k_v #using formula TGI=1-RTV(a)/RTV(c)
  #k_a=log((1-tgi_a)*(exp(k_v*t)-1)+1)/t #using formula TGI=1-(v0*e(ka*t)-v0)/(v0*e(kv*t)-v0) to calculate ka, assuming v0 is equal between treatment and vehicle
  #k_b=log(1-tgi_b)/t+k_v
  #k_b=log((1-tgi_b)*(exp(k_v*t)-1)+1)/t
  #k_c=log(1-tgi_c)/t+k_v
  #k_c=log((1-tgi_c)*(exp(k_v*t)-1)+1)/t
  kit_vec=rep(c(0,kit_a,kit_b,kit_c),each=n) #kickin time vector
  rp_vec=rep(c(0,rp_a,rp_b,rp_c),each=n)
  irt_vec=rep(c(0,irt_a,irt_b,irt_c),each=n) #induced resistance time vector
  k_vec=as.vector(sapply(1:4,function(i) rnorm(n,mean=k_v,sd=sd_k)))#vector of initial tumor growth rate
  #k_vec=as.vector(sapply(c(k_v,k_a,k_b,k_c),function(k) rnorm(n,mean=k,sd=sigma_k)))
  tv=data.frame(Treatment=rep(treatments,each=n),Mouse=1:(4*n),Day=0,TV=tv0s,kit=kit_vec,rp=rp_vec,irt=irt_vec,k_c=k_vec,drug_eff=drug_eff_vec) %>%
    mutate(k=ifelse(kit==0,k_c+drug_eff,k_c),rpv=ifelse(irt==0,rp_vec*TV,0))
  tvs=list()
  tvs[[1]]=tv
  for(i in 1:t){
    #prev_t=t_vec[i-1]
    #cur_t=t_vec[i]
    #time_span=cur_t-prev_t
    #using exponential growth model to model tumor growth(can also try gompertz model)
    tvs[[i+1]]=tvs[[i]] %>% mutate(Day=i) %>% mutate(k=ifelse(Day>=kit,k_c+drug_eff,k_c),rpv=ifelse(Day>=irt & rpv==0,rp*TV,rpv)) %>% mutate(TV=(TV-rpv)*exp(k)+rpv*exp(k_c),rpv=rpv*exp(k_c)) %>%
      mutate(TV=TV+rnorm(n*4,sd=sigma)) %>% mutate(TV=ifelse(TV>=eu_tv,NA,TV))
    #tv_t=tvs[[1]] %>% mutate(Day=cur_t,TV=TV*exp(k*cur_t)) %>% mutate(TV=TV+rnorm(n*4,sd=sd0))
    #tvs[[i]]=tv_t
  }
  tv_long=do.call(rbind,tvs) %>% as.data.frame()
  tv_long=tv_long %>% mutate(Group=factor(paste0("Group ",as.numeric(Treatment)))) %>% select(Group,Treatment,Mouse,Day,TV)
  tv_long=expand_tv(tv_long)
  #list(tv_long=tv_long,tv0=tv)
  tv_long
}
