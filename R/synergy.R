#' Synergy calculation based on TGI (TGI can either be defined by 1-delta(t)/delta(c) or 1-RTV(t)/RTV(c))
#'
#' @param TGI_lst An list object returned by getTGI function
#' @param ci confidence intervaln for TGI based synergy score, can be 0.95,0.9,0.8, etc
#' @param ci_type type of bootstrap confidence interval, can be stud,perc or bca
#' @param save save image
#' @param method method for synergy calculation, can be Bliss,HSA or RA
#' @param display whether or not display figure
#'
#' @return Result of synergy calculation
#' @export
#'
#' @examples
#' TGI_lst=getTGI(LS_1034,17)
#' bliss_synergy_TGI=TGI_synergy(TGI_lst)
TGI_synergy=function(TGI_lst,method='Bliss',ci=0.95,ci_type='perc',display=TRUE,save=TRUE){
  data=TGI_lst$bsTGI_r$data
  bsTGI_df=TGI_lst$bsTGI_df
  d1_TGI=bsTGI_df[1,'TGI']
  d2_TGI=bsTGI_df[2,'TGI']
  expected_TGI=switch(method,"Bliss"=100*(d1_TGI/100+d2_TGI/100-d1_TGI*d2_TGI/10000),
                      "HSA"=max(d1_TGI,d2_TGI),
                      "RA"=d1_TGI+d2_TGI)
  #synergy type,can either be Combination Index(CI:log(Observed effect/Expeced effect)) or synergy score (Observed effect-expected effect)
  synergy_score=bsTGI_df[3,'TGI']-expected_TGI
  CI=log(bsTGI_df[3,'TGI']/expected_TGI)

  p1=ggplot(bsTGI_df, aes(x=Treatment, y=TGI)) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_errorbar(aes(ymin=TGI, ymax=TGI+std.err), width=.2,size=1,
                  position=position_dodge(.9))+xlab("Treatment")+ylab("Inhibition(%)")+
    geom_hline(yintercept = expected_TGI,color='red',linetype=2)+
    geom_text(aes(x=-Inf,y=expected_TGI+2),label=paste0('Expeted inhibition by ',method),vjust=0,hjust=0,size=rel(1.2))+
    theme_Publication()


  bsTGI_all=TGI_lst$bsTGI_r$t
  colnames(bsTGI_all)=c('D1','D2','Combo')
  bsTGI_all = bsTGI_all %>% as.data.frame() %>% mutate('expected_TGI'=switch(method,"Bliss"=D1+D2-D1*D2,"HSA"=pmax(D1,D2),"RA"=D1+D2))
  bsTGI_all = bsTGI_all %>% mutate(Synergy_score=Combo-expected_TGI)
  bsTGI_all = bsTGI_all*100 #in percentage
  pval=round(1-sum(bsTGI_all$Combo>bsTGI_all$expected_TGI)/nrow(bsTGI_all),4)
  p2=bsTGI_all %>% ggplot()+aes(expected_TGI,Combo)+geom_point()+xlab("Expected TGI")+ylab("Observed Combo TGI")+
    geom_abline(slope = 1,color='red',linetype=2)+annotate('text',x=-Inf,y=Inf,hjust=0,vjust=1,label=paste0("Bootstrap Pvalue=",pval),size=7)+
    theme_Publication()
  figure=ggpubr::ggarrange(p1,p2,labels=c("A","B"),ncol=2)
  if(display) print(figure)
  if(save) ggsave(paste0("TGI synergy plot by ",method,'.png'),width=16,height=8,dpi=300)

  bsTGI_r=TGI_lst$bsTGI_r
  bsTGI_r$t=bsTGI_all
  bsTGI_r$t0=c(TGI_lst$bsTGI_df$TGI,expected_TGI,synergy_score)
  ss_ci=getCI(bsTGI_r,c(5,5),conf=ci,ci_type=ci_type)
  c('Expected TGI'=expected_TGI,'Observed TGI'=bsTGI_df[3,'TGI'],'p.val'=pval,'CI'=CI,
    "Synergy score"=synergy_score,'ss_lb'=ss_ci[1],'ss_ub'=ss_ci[2])
}


#' AUC based synergy calculation for bootstrap
#'
#' @param auc_mouse auc dataframe for each mouse
#' @param t date to estimate survival (default 14 days)
#' @param method method for synergy calculation, can be Bliss,HSA or RA
#' @param idx bootstrap index
#'
#' @return combination index (CI) and synergy score
bs_AUC_synergy=function(auc_mouse,t=21,method='Bliss',idx){
  #here normalized auc is the estimate of tumor growth rate(k),survival is calculated using exp(k*t)
  auc_mouse=auc_mouse[idx,]
  #auc_v=auc_mouse %>% filter(Group=="Group 1") %>% pull(AUC)
  #auc_a=auc_mouse %>% filter(Group=="Group 2") %>% pull(AUC)
  #auc_b=auc_mouse %>% filter(Group=="Group 3") %>% pull(AUC)
  #auc_c=auc_mouse %>% filter(Group=="Group 4") %>% pull(AUC)
  #aucs=expand.grid('AUC_v'=auc_v,'AUC_a'=auc_a,'AUC_b'=auc_b,'AUC_c'=auc_c)
  #aucs=aucs[complete.cases(aucs),]
  #aucs=aucs %>% mutate(s_a=exp((AUC_a-AUC_v)*t),s_b=exp((AUC_b-AUC_v)*t),s_c=exp((AUC_c-AUC_v)*t))
  #not making sense for AUC based synergy, considering expected RA survivalcan be less than 0
  #aucs = aucs %>% mutate('s_e'=switch(method,"Bliss"=s_a*s_b,"HSA"=pmin(s_a,s_b),"RA"=s_a+s_b-1))
  #aucs = aucs %>% mutate(CI=s_c/s_e,Synergy_score=100*(s_e-s_c))
  #c(median_CI=median(aucs$CI,na.rm = T),median_ss=median(aucs$Synergy_score,na.rm=T))
  #c(median_CI=median(aucs$CI,na.rm = T),median_ss=median(aucs$Synergy_score,na.rm=T))
  #auc_mean=auc_mouse %>% filter(!is.na(AUC)) %>% group_by(Group) %>% summarise(auc=mean(AUC))
  auc_mean=auc_mouse %>% filter(!is.na(AUC)) %>% group_by(Group) %>% summarise(auc=mean(AUC))
  auc_mean=auc_mean %>% mutate(auc_v=auc_mean[['auc']][1]) %>% mutate(s=exp((auc-auc_v)*t))
  s_vec=auc_mean %>% pull(s)
  s_e=switch(method,"Bliss"=s_vec[2]*s_vec[3],"HSA"=min(s_vec[2],s_vec[3]),"RA"=s_vec[2]+s_vec[3]-1)
  c(CI=s_vec[4]/s_e,Synergy_score=100*(s_e-s_vec[4]))
}


#' Calculatr AUC based synergy scores and its bootstrap confidence interval
#'
#' @param auc_lst results from AUC calculation
#' @param t date to estimate survival (default 14 days)
#' @param method method for synergy calculation, can be Bliss,HSA or RA
#' @param ci confidence intervaln for AUC based synergy score, can be 0.95,0.9,0.8, etc
#' @param ci_type type of bootstrap confidence interval, can be stud,perc or bca
#' @param save save image
#' @param display print image
#' @param boot_n number of bootstrap resample
#'
#' @return a dataframe of synergy scores and its bootstrap confidence interval
#' @export
#'
#' @examples
#' auc_lst=get_mAUCr(LS_1034)
#' bliss_synergy_AUC=AUC_synergy(auc_lst)
AUC_synergy=function(auc_lst,t=21,method='Bliss',boot_n=1000,ci=0.95,ci_type='perc',display=TRUE,save=TRUE,kw='Test'){
  auc_mouse=as.data.frame(auc_lst$auc_mouse)
  bsAUCci_r=boot::boot(data=auc_mouse,statistic=bs_AUC_synergy,t=t,method=method,strata=auc_mouse$Group,R=boot_n)
  #bsAUCci_r=readRDS('SW837_boot_paper.Rdata')
  n=length(bsAUCci_r$t0)
  cis=do.call(rbind,lapply(1:n,function(i) getCI(bsAUCci_r,i,conf=ci,ci_type=ci_type)))
  out_df=cbind(broom::tidy(bsAUCci_r),cis)
  out_df=out_df[,-3]
  names(out_df)=c('Metric','Value','std.err','lb','ub')
  bs_df=bsAUCci_r$t %>% as.data.frame()
  #colnames(bs_df)=c("mAUC-CI","mAUC-Synergy Score")
  #define name of synergy scores
  ss_names=paste0(method,c(" CI"," Synergy Score"),'(invivoSyn)')
  colnames(bs_df)=ss_names
  pval_CI=mean(bs_df[[ss_names[1]]]>=1,na.rm=T)
  pval_SS=mean(bs_df[[ss_names[2]]]<=0,na.rm=T)
  p1=plot_density(bs_df,ss_names[1],1,pval_CI,pe=out_df[1,'Value'],lb = out_df[1,'lb'],ub = out_df[1,'ub'])
  p2=plot_density(bs_df,ss_names[2],0,pval_SS,pe=out_df[2,'Value'],lb = out_df[2,'lb'],ub = out_df[2,'ub'])
  figure=ggpubr::ggarrange(p1,p2,ncol=2)#labels=c("A","B"),
  if(display) print(figure)
  if(save) ggsave(paste0(kw," AUC synergy plot by ",method,'.png'),width=17,height=8,dpi=300)
  out_df=bind_cols(out_df,data.frame(p.val=c(pval_CI,pval_SS)))
  out_df
}


#' Calculate synergy based on linear mixed model
#'
#' @param tv tumor growth data
#' @param sel_day the last day selected for tumor growth data, if not defined, use all data
#'
#' @return a data frame, linear mixed model results
#' @export
#'
#' @examples
#' lmm_synergy(LS_1034)
lmm_synergy=function(tv,sel_day=NA){
  if(!is.na(sel_day)) tv=tv %>% filter(Day <= sel_day)
  tv=tv %>% mutate(Treatment=make.names(Treatment))
  group_info=tv %>% select(Group,Treatment) %>% distinct()
  trts=as.character(group_info$Treatment[2:3])
  for(tr in trts) tv[[tr]]=ifelse(grepl(tr,tv$Treatment),1,0)
  lmm_f=paste0("logTV~log(TV0)+Day+Day:(",trts[1],'*',trts[2],')')
  lmm1<-nlme::lme(as.formula(lmm_f), random=list(Mouse=~Day-1), method='REML', data=tv)
  lmm_results=broom.mixed::tidy(lmm1)
  lmm_results=lmm_results %>% filter(effect=='fixed') %>% select(3:8)
  lmm_results
}

#' Calculate Day-specifc CI based on CombPDX's method
#'
#' @param tv tumor growth data
#' @param sel_day the specific day for CI calcuation
#' @param method method for synergy calculation, can be Bliss,HSA or RA
#' @param ci confidence interval for CombPDX global CI, can be 0.95,0.9,0.8, etc;standard error based on delta method
#'
#' @return a vector of CI and it's confidence interval
#' @export
#'
#' @examples
CombPDX_CI=function(tv,sel_day=NA,method='Bliss',ci=0.95){
  if(!is.na(sel_day)) tv=tv %>% filter(Day == sel_day)
  RTV=tv %>% group_by(Group) %>% summarise(mean_RTV=mean(RTV),n=n(),variance=var(RTV))
  RTV_vec=RTV %>% pull(mean_RTV)
  group_idx=ifelse(RTV_vec[2]<RTV_vec[3],2,3) #deterimine which group have better efficacy
  n_vec=RTV %>% pull(n)
  var_vec=RTV %>% pull(variance)
  CI=switch(method,"Bliss"=log(RTV_vec[2])+log(RTV_vec[3])-log(RTV_vec[1])-log(RTV_vec[4]),
            "HSA"=log(RTV_vec[group_idx])-log(RTV_vec[4]),
            "RA"=log(RTV_vec[2]+RTV_vec[3])-log(RTV_vec[1]+RTV_vec[4]))
  var_CI=switch(method,"Bliss"=var_vec[1]/(n_vec[1]*(RTV_vec[1]^2))+var_vec[2]/(n_vec[2]*(RTV_vec[2]^2))+
                  var_vec[3]/(n_vec[3]*(RTV_vec[3]^2))+var_vec[4]/(n_vec[4]*(RTV_vec[4]^2)),
            "HSA"=var_vec[group_idx]/(n_vec[group_idx]*(RTV_vec[group_idx]^2))+var_vec[4]/(n_vec[4]*(RTV_vec[4]^2)),
            "RA"=var_vec[2]/(n_vec[2]*((RTV_vec[2]+RTV_vec[3])^2))+var_vec[3]/(n_vec[3]*((RTV_vec[2]+RTV_vec[3])^2))+
              var_vec[1]/(n_vec[1]*((RTV_vec[1]+RTV_vec[4])^2))+var_vec[4]/(n_vec[4]*((RTV_vec[1]+RTV_vec[4])^2)))
  q1=qnorm(ci)
  c('CI'=CI,'std.err'=sqrt(var_CI),'lb'=CI-q1*sqrt(var_CI),'ub'=CI+q1*sqrt(var_CI),'p.val'=1-pnorm(CI,sd=sqrt(var_CI)))#p.val'=2*(1-pnorm(abs(CI),sd=sqrt(var_CI))))
}


#' Bootstrap function to Calculate global CI based on TGI definition
#'
#' @param tv_wide tumor growth data in wide format (can be eiter deltaTV or RTV (mouse,Treatment,Group,Day...)), for bootstrap mouse within group
#' @param method method for synergy calculation, can be Bliss,HSA or RA
#' @param idx bootstrap index
#'
#' @return a dataframe of synergy scores and its bootstrap confidence interval
#' @export
#'
#' @examples
bs_global_CI=function(tv_wide,method='Bliss',idx){
  tv_wide=tv_wide[idx,]
  #if(!is.na(sel_day)) tv=tv %>% filter(Day <= sel_day)
  rtv_long=tv_wide %>% pivot_longer(c(-1,-2,-3),names_to='Day',values_to = 'RTV')
  rtv_long=rtv_long[complete.cases(rtv_long),]
  mean_rtv=suppressMessages(rtv_long %>% group_by(Day,Group) %>% summarise(mean_RTV=mean(RTV)) %>%
    pivot_wider(id_cols = Day,names_from = Group,values_from = mean_RTV) %>% mutate(Day=as.numeric(Day)) %>%
    arrange(Day) %>% na.omit() %>% filter(Day>0))
  #When a lot of RTV equals to 0, could be problematic
  CI=switch(method,"Bliss"=mean_rtv %>% mutate(CI=log(`Group 2`)+log(`Group 3`)-log(`Group 1`)-log(`Group 4`)) %>% pull(CI) %>% mean,
            "HSA"=mean_rtv %>% mutate(lowm=min(`Group 2`,`Group 3`)) %>% mutate(CI=log(lowm)-log(`Group 4`)) %>% pull(CI) %>% mean,
            "RA"=mean_rtv %>% mutate(CI=log(`Group 2`+`Group 3`)-log(`Group 1`+`Group 4`)) %>% pull(CI) %>% mean)
  CI
}

#' Calculatr global CI from CombPDX and its bootstrap confidence interval
#'
#' @param tv tumor growth data
#' @param method method for synergy calculation, can be Bliss,HSA or RA
#' @param ci confidence interval for CombPDX global CI, can be 0.95,0.9,0.8, etc
#' @param ci_type type of bootstrap confidence interval, can be stud,perc or bca
#' @param save save image
#' @param display print image
#' @param boot_n number of bootstrap resample
#'
#' @return a dataframe of synergy scores and its bootstrap confidence interval
#' @export
#'
#' @examples
global_CI_synergy=function(tv,method='Bliss',boot_n=1000,ci=0.95,ci_type='perc',display=TRUE,save=TRUE){
  RTV_wide=tv %>% pivot_wider(1:3,names_from = Day,values_from = RTV)
  bsAUCci_r=boot::boot(data=RTV_wide,statistic=bs_global_CI,method=method,strata=RTV_wide$Group,R=boot_n)
  n=length(bsAUCci_r$t0)
  cis=do.call(rbind,lapply(1:n,function(i) getCI(bsAUCci_r,i,conf=ci,ci_type=ci_type)))
  out_df=cbind(broom::tidy(bsAUCci_r),cis)
  names(out_df)=c('Global CI','Bias','std.err','lb','ub')
  bs_df=bsAUCci_r$t %>% as.data.frame()
  colnames(bs_df)=c("CombPDX-gCI")
  pval_CI=mean(bs_df[['CombPDX-gCI']]<=0,na.rm=T)
  p1=plot_density(bs_df,'CombPDX-gCI',0,pval_CI,pe=out_df[1,'Global CI'],lb = out_df[1,'lb'],ub = out_df[1,'ub'])
  if(display) print(p1)
  if(save) ggsave(paste0("Global CI synergy plot by ",method,'.png'),width=8,height=8,dpi=300)
  out_df=bind_cols(out_df,data.frame(p.val=pval_CI))
  out_df
}


