#' Synergy calculation based on TGI (TGI can either be defined by 1-delta(t)/delta(c) or 1-RTV(t)/RTV(c))
#'
#' Generalized to N-drug combinations: handles any number of single-agent
#' arms via the `roles` list carried in `TGI_lst$roles`.
#'
#' @param TGI_lst An list object returned by getTGI function
#' @param ci confidence intervaln for TGI based synergy score, can be 0.95,0.9,0.8, etc
#' @param ci_type type of bootstrap confidence interval, can be stud,perc or bca
#' @param save save image
#' @param method method for synergy calculation, can be Bliss,HSA or RA
#' @param display whether or not display figure
#' @param file if save is TRUE, the name of output file
#' @param width width (in inches) of the saved figure
#' @param height height (in inches) of the saved figure
#'
#' @return a data frame of synergy scores (CI and Synergy Score) and their
#'   bootstrap confidence intervals, in the same format as [AUC_synergy()]
#' @export
#'
#' @examples
#' TGI_lst=getTGI(LS_1034,17)
#' bliss_synergy_TGI=TGI_synergy(TGI_lst)
TGI_synergy=function(TGI_lst,method='Bliss',ci=0.95,ci_type='perc',display=TRUE,save=TRUE,file="TGI_synergy",width=16,height=8){
  bsTGI_df=TGI_lst$bsTGI_df
  roles=TGI_lst$roles
  if(is.null(roles)) stop("TGI_lst missing `roles`; re-run getTGI() on a tv produced by read_tv().")
  combo_row=which(bsTGI_df$Group==roles$combo_grp)
  if(length(combo_row)!=1L) stop("Combo group not found in bsTGI_df.")
  single_rows=which(bsTGI_df$Group %in% roles$single_grps)
  single_tgi=bsTGI_df$TGI[single_rows]
  single_surv=1-single_tgi/100
  combo_tgi=bsTGI_df$TGI[combo_row]
  expected_TGI=switch(method,"Bliss"=100*(1-prod(single_surv)),
                      "HSA"=max(single_tgi),
                      "RA"=sum(single_tgi))
  synergy_score=combo_tgi-expected_TGI
  CI=log(combo_tgi/expected_TGI)

  p1=ggplot(bsTGI_df, aes(x=Treatment, y=TGI)) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_errorbar(aes(ymin=TGI, ymax=TGI+std.err), width=.2,size=1,
                  position=position_dodge(.9))+xlab("Treatment")+ylab("Inhibition(%)")+
    geom_hline(yintercept = expected_TGI,color='red',linetype=2)+
    geom_text(aes(x=-Inf,y=expected_TGI+2),label=paste0('Expected inhibition by ',method),vjust=0,hjust=0,size=7)+
    theme_Publication()+
    theme(axis.text.x=element_text(angle=15,hjust=1))

  bsTGI_all=TGI_lst$bsTGI_r$t
  bs_cols=c(as.character(roles$single_grps),'Combo')
  colnames(bsTGI_all)=bs_cols
  bsTGI_all=as.data.frame(bsTGI_all)
  single_mat=as.matrix(bsTGI_all[,seq_along(roles$single_grps),drop=FALSE])
  expected_TGI_bs=switch(method,
                         "Bliss"=1-apply(1-single_mat,1,prod),
                         "HSA"=apply(single_mat,1,max),
                         "RA"=rowSums(single_mat))
  bsTGI_all$expected_TGI=expected_TGI_bs
  bsTGI_all$Synergy_score=bsTGI_all$Combo-bsTGI_all$expected_TGI
  bsTGI_all=bsTGI_all*100
  # CI (log ratio) is scale-invariant, so it is computed after the *100 rescale.
  bsTGI_all$CI=log(bsTGI_all$Combo/bsTGI_all$expected_TGI)
  pval_syn=round(1-sum(bsTGI_all$Combo>bsTGI_all$expected_TGI)/nrow(bsTGI_all),4)
  pval_anta=round(1-sum(bsTGI_all$Combo<bsTGI_all$expected_TGI)/nrow(bsTGI_all),4)
  p2=bsTGI_all %>% ggplot()+aes(expected_TGI,Combo)+geom_point()+xlab("Expected TGI")+ylab("Observed Combo TGI")+
    geom_abline(slope = 1,color='red',linetype=2)+annotate('text',x=-Inf,y=Inf,hjust=0,vjust=1,label=paste0("Bootstrap Pvalue=",pval_syn),size=7)+
    theme_Publication()
  figure=ggpubr::ggarrange(p1,p2,labels=c("A","B"),ncol=2)
  if(display) print(figure)
  if(save) ggsave(paste0(file,'.png'),width=width,height=height,dpi=300)

  # Bootstrap confidence intervals for both the synergy score and the CI.
  bsTGI_r=TGI_lst$bsTGI_r
  bsTGI_r$t=as.matrix(bsTGI_all)
  bsTGI_r$t0=c(TGI_lst$bsTGI_df$TGI,expected_TGI,synergy_score,CI)
  ncol_t=ncol(bsTGI_r$t)
  ss_idx=ncol_t-1L  # Synergy_score column
  ci_idx=ncol_t     # CI column
  ss_ci=getCI(bsTGI_r,ss_idx,conf=ci,ci_type=ci_type)
  ci_ci=getCI(bsTGI_r,ci_idx,conf=ci,ci_type=ci_type)
  # Mirror AUC_synergy()'s output: one row per metric (CI, Synergy Score).
  out_df=data.frame(
    Metric=paste0(method,c(" CI"," Synergy Score"),'(invivoSyn)'),
    Value=c(CI,synergy_score),
    std.err=c(stats::sd(bsTGI_all$CI,na.rm=TRUE),stats::sd(bsTGI_all$Synergy_score,na.rm=TRUE)),
    lb=c(ci_ci[1],ss_ci[1]),
    ub=c(ci_ci[2],ss_ci[2]),
    stringsAsFactors=FALSE
  )
  out_df=bind_cols(out_df,data.frame('p.val.Synergy'=c(pval_syn,pval_syn),
                                     'p.val.Antagonism'=c(pval_anta,pval_anta)))
  out_df
}


#' AUC based synergy calculation for bootstrap
#'
#' @param auc_mouse auc dataframe for each mouse
#' @param t date to estimate survival (default 21 days)
#' @param method method for synergy calculation, can be Bliss,HSA or RA
#' @param idx bootstrap index
#' @param roles roles list (output of [get_roles()])
#'
#' @return combination index (CI) and synergy score
bs_AUC_synergy=function(auc_mouse,t=21,method='Bliss',roles=NULL,idx){
  auc_mouse=auc_mouse[idx,]
  if(is.null(roles)){
    grps=levels(droplevels(auc_mouse$Group))
    roles=list(vehicle_grp=grps[1],
               single_grps=grps[2:(length(grps)-1)],
               combo_grp=grps[length(grps)])
  }
  auc_mean=auc_mouse %>% filter(!is.na(AUC)) %>% group_by(Group) %>% summarise(auc=mean(AUC),.groups='drop')
  auc_v=auc_mean$auc[auc_mean$Group==roles$vehicle_grp]
  s_single=exp((auc_mean$auc[match(roles$single_grps,auc_mean$Group)]-auc_v)*t)
  s_combo=exp((auc_mean$auc[auc_mean$Group==roles$combo_grp]-auc_v)*t)
  s_e=switch(method,"Bliss"=prod(s_single),"HSA"=min(s_single),"RA"=sum(1-s_single))
  c(CI=s_combo/s_e,Synergy_score=100*(s_e-s_combo))
}


#' Calculatr AUC based synergy scores and its bootstrap confidence interval
#'
#' @param auc_lst results from AUC calculation
#' @param t date to estimate survival (default 21 days)
#' @param method method for synergy calculation, can be Bliss,HSA or RA
#' @param ci confidence intervaln for AUC based synergy score, can be 0.95,0.9,0.8, etc
#' @param ci_type type of bootstrap confidence interval, can be stud,perc or bca
#' @param save save image
#' @param display print image
#' @param boot_n number of bootstrap resample
#' @param kw legend keyword for plotting
#' @param file if save is TRUE, the name of output file
#' @param parallel whether or not usint parallel computing
#'
#' @return a dataframe of synergy scores and its bootstrap confidence interval
#' @import parallel
#' @export
#'
#' @examples
#' auc_lst=get_mAUCr(LS_1034)
#' bliss_synergy_AUC=AUC_synergy(auc_lst)
AUC_synergy=function(auc_lst,t=21,method='Bliss',boot_n=1000,ci=0.95,
                     ci_type='perc',kw='Test',file="eGR_synergy",save=T,display=F,
                     parallel='no'){
  auc_mouse=as.data.frame(auc_lst$auc_mouse)
  roles=auc_lst$roles
  if(is.null(roles)) stop("auc_lst missing `roles`; re-run get_mAUCr() on a tv produced by read_tv().")
  set.seed(123456)
  bsAUCci_r=boot::boot(data=auc_mouse,statistic=bs_AUC_synergy,t=t,method=method,roles=roles,strata=auc_mouse$Group,
                       R=boot_n,parallel = parallel,ncpus = parallel::detectCores()-1)
  n=length(bsAUCci_r$t0)
  cis=do.call(rbind,lapply(1:n,function(i) getCI(bsAUCci_r,i,conf=ci,ci_type=ci_type)))
  # Build name-aligned rather than via broom::tidy(), which drops rows for NA
  # statistics and would then mis-align against `cis`.
  out_df=data.frame(
    Metric=names(bsAUCci_r$t0),
    Value=as.numeric(bsAUCci_r$t0),
    std.err=apply(bsAUCci_r$t,2,stats::sd,na.rm=TRUE),
    lb=cis[,1],
    ub=cis[,2],
    stringsAsFactors=FALSE
  )
  bs_df=bsAUCci_r$t %>% as.data.frame()
  ss_names=paste0(method,c(" CI"," Synergy Score"),'(invivoSyn)')
  colnames(bs_df)=ss_names
  pval_CI_syn=mean(bs_df[[ss_names[1]]]>1,na.rm=T)
  pval_SS_syn=mean(bs_df[[ss_names[2]]]<0,na.rm=T)
  pval_CI_anta=mean(bs_df[[ss_names[1]]]<1,na.rm=T)
  pval_SS_anta=mean(bs_df[[ss_names[2]]]>0,na.rm=T)
  p1=plot_density(bs_df,ss_names[1],1,pval_CI_syn,pe=out_df[1,'Value'],lb = out_df[1,'lb'],ub = out_df[1,'ub'])
  p2=plot_density(bs_df,ss_names[2],0,pval_SS_syn,pe=out_df[2,'Value'],lb = out_df[2,'lb'],ub = out_df[2,'ub'])
  figure=ggpubr::ggarrange(p1,p2,ncol=2)
  if(display) print(figure)
  if(save){
    png(paste0(file,'.png'),width = 17,height = 8,units = 'in',res=300)
    print(figure)
    dev.off()
  }
  out_df=bind_cols(out_df,data.frame('p.val.Synergy'=c(pval_CI_syn,pval_SS_syn),
                                     'p.val.Antagonism'=c(pval_CI_anta,pval_SS_anta)))
  out_df
}


#' Calculate synergy based on linear mixed model
#'
#' Generalized to N-drug combinations. Builds indicator variables for each
#' single-agent treatment in `roles$singles` and fits
#' `logTV ~ log(TV0) + Day + Day:(d1*d2*...*dN)` so that the N-way
#' interaction term tests for synergy. For N>3 the number of interaction
#' terms grows as 2^N, so a warning is emitted.
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
  roles=get_roles(tv)
  if(!is.na(sel_day)) tv=tv %>% filter(Day <= sel_day)
  participating=c(roles$vehicle_grp,roles$single_grps,roles$combo_grp)
  tv=tv %>% filter(Group %in% participating)
  tv=tv %>% mutate(Treatment=make.names(as.character(Treatment)))
  singles_nm=make.names(roles$singles)
  combo_nm=make.names(roles$combo)
  if(length(singles_nm)>3L)
    warning("lmm_synergy: ", length(singles_nm), " single agents requested; the full ",
            length(singles_nm), "-way interaction has ", 2^length(singles_nm)-1,
            " terms and may not be identifiable.", call.=FALSE)
  for(tr in singles_nm){
    tv[[tr]]=as.integer(tv$Treatment == tr | tv$Treatment == combo_nm)
  }
  first_order=paste0("Day:",singles_nm,collapse='+')
  highest=paste0("Day:",paste(singles_nm,collapse=':'))
  lmm_f=paste0("logTV~log(TV0)+Day+",first_order,"+",highest)
  lmm1<-nlme::lme(as.formula(lmm_f), random=list(Mouse=~Day-1), method='REML', data=tv)
  broom.mixed::tidy(lmm1) %>% filter(effect=='fixed') %>% select(3:8)
}

#' Calculate Day-specifc CI based on CombPDX's method
#'
#' Generalized to N-drug combinations via the closed-form Bliss/HSA/RA
#' formulas in RTV space. Variance is computed via the delta method.
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
  roles=get_roles(tv)
  if(!is.na(sel_day)) tv=tv %>% filter(Day == sel_day)
  RTV=tv %>% group_by(Group) %>% summarise(mean_RTV=mean(RTV),n=n(),variance=var(RTV),.groups='drop')
  pick=function(g) which(RTV$Group==g)
  i_v=pick(roles$vehicle_grp); i_c=pick(roles$combo_grp)
  i_s=vapply(roles$single_grps,pick,integer(1L))
  N=length(i_s)
  rtv_v=RTV$mean_RTV[i_v]; rtv_c=RTV$mean_RTV[i_c]; rtv_s=RTV$mean_RTV[i_s]
  var_v=RTV$variance[i_v]; var_c=RTV$variance[i_c]; var_s=RTV$variance[i_s]
  n_v=RTV$n[i_v];           n_c=RTV$n[i_c];           n_s=RTV$n[i_s]
  if(method=='Bliss'){
    CI=sum(log(rtv_s))-(N-1)*log(rtv_v)-log(rtv_c)
    var_CI=sum(var_s/(n_s*rtv_s^2))+(N-1)^2*var_v/(n_v*rtv_v^2)+var_c/(n_c*rtv_c^2)
  } else if(method=='HSA'){
    j=which.min(rtv_s)
    CI=log(rtv_s[j])-log(rtv_c)
    var_CI=var_s[j]/(n_s[j]*rtv_s[j]^2)+var_c/(n_c*rtv_c^2)
  } else if(method=='RA'){
    sum_s=sum(rtv_s); denom=(N-1)*rtv_v+rtv_c
    CI=log(sum_s)-log(denom)
    var_CI=sum(var_s/(n_s*sum_s^2))+(N-1)^2*var_v/(n_v*denom^2)+var_c/(n_c*denom^2)
  } else stop("Unknown method: ",method)
  q1=qnorm(ci)
  c('CI'=CI,'std.err'=sqrt(var_CI),'lb'=CI-q1*sqrt(var_CI),'ub'=CI+q1*sqrt(var_CI),
    'p.val'=1-pnorm(CI,sd=sqrt(var_CI)))
}


#' Bootstrap function to Calculate global CI based on TGI definition
#'
#' @param tv_wide tumor growth data in wide format (mouse, Treatment, Group, day columns); bootstrap mouse within group
#' @param method method for synergy calculation, can be Bliss,HSA or RA
#' @param roles roles list (output of [get_roles()])
#' @param idx bootstrap index
#'
#' @return mean CI across days (a scalar)
#' @export
bs_global_CI=function(tv_wide,method='Bliss',roles=NULL,idx){
  tv_wide=tv_wide[idx,]
  if(is.null(roles)){
    grps=levels(droplevels(tv_wide$Group))
    roles=list(vehicle_grp=grps[1],
               single_grps=grps[2:(length(grps)-1)],
               combo_grp=grps[length(grps)])
  }
  rtv_long=tv_wide %>% pivot_longer(c(-1,-2,-3),names_to='Day',values_to = 'RTV')
  rtv_long=rtv_long[complete.cases(rtv_long),]
  mean_rtv=suppressMessages(rtv_long %>% group_by(Day,Group) %>% summarise(mean_RTV=mean(RTV),.groups='drop') %>%
    pivot_wider(id_cols = Day,names_from = Group,values_from = mean_RTV) %>% mutate(Day=as.numeric(Day)) %>%
    arrange(Day) %>% na.omit() %>% filter(Day>0))
  N=length(roles$single_grps)
  v=mean_rtv[[roles$vehicle_grp]]; cc=mean_rtv[[roles$combo_grp]]
  s_mat=as.matrix(mean_rtv[,as.character(roles$single_grps),drop=FALSE])
  if(method=='Bliss'){
    ci_per_day=rowSums(log(s_mat))-(N-1)*log(v)-log(cc)
  } else if(method=='HSA'){
    ci_per_day=log(apply(s_mat,1,min))-log(cc)
  } else if(method=='RA'){
    ci_per_day=log(rowSums(s_mat))-log((N-1)*v+cc)
  } else stop("Unknown method: ",method)
  mean(ci_per_day)
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
  roles=get_roles(tv)
  participating=c(roles$vehicle_grp,roles$single_grps,roles$combo_grp)
  tv=tv %>% filter(Group %in% participating) %>% mutate(Group=droplevels(Group))
  RTV_wide=tv %>% pivot_wider(1:3,names_from = Day,values_from = RTV)
  bsAUCci_r=boot::boot(data=RTV_wide,statistic=bs_global_CI,method=method,roles=roles,
                       strata=RTV_wide$Group,R=boot_n)
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
