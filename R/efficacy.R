#' Calculate TGI statistic for bootstrap
#'
#' @param tvs tumor growth data of one specific day
#' @param ref_group reference group,by default the first group (vehicle)
#' @param idx bootstrap index
#' @param tv_var variable for TGI calculation, can be either DeltaTV or RTV (RTV: definiton of CombPDX for TGI)
#'
#' @return TGI statistics for each group
bsTGI=function(tvs,tv_var='DeltaTV',ref_group="Group 1",idx){
  tvs=tvs[idx,]
  ms=tapply(tvs[[tv_var]],tvs$Group,mean)
  ms1=ms[-which(names(ms)==ref_group)]/ms[ref_group]
  ms1=1-ms1
  ms1
}

#' Calculate TGI using stratified bootstrap,with confidence interval
#'
#' @param tv object return from read_tv
#' @param sel_day timepoint for TGI calculation
#' @param ref_group name of reference group
#' @param ci confidence interval, can be 0.95,0.9,0.8, etc
#' @param ci_type type of bootstrap confidence interval, can be stud,perc or bca
#' @param tv_var variable for TGI calculation, can be either DeltaTV or RTV (RTV: definiton of CombPDX for TGI)
#'
#' @return a list of bootstrap result and TGI statistic (TGI and its confidence interval for each group)
#' @export
#'
#' @examples
#' getTGI(LS_1034,17,ci=0.9,ci_type='bca')
getTGI=function(tv,sel_day,tv_var='DeltaTV',ref_group='Group 1',ci=0.95,ci_type='perc',n_rep=1000){
  dayTV=subset(tv,Day==sel_day)
  rownames(dayTV)=dayTV$Mouse
  bsTGI_r=boot::boot(data=dayTV,statistic=bsTGI,tv_var=tv_var,ref_group=ref_group,strata=dayTV$Group,R=n_rep)
  group_n=length(bsTGI_r$t0)
  cis=do.call(rbind,lapply(1:group_n,function(i) getCI(bsTGI_r,i,conf=ci,ci_type=ci_type)))
  bsTGI_df=cbind(broom::tidy(bsTGI_r),cis)
  bsTGI_df=bsTGI_df[,-3]
  names(bsTGI_df)=c('Group','TGI','std.err','lb','ub')
  bsTGI_df[,-1]=100*bsTGI_df[,-1]
  group_info=tv %>% select(Group,Treatment) %>% distinct()
  bsTGI_df=bsTGI_df %>% left_join(group_info) %>% select(Group,Treatment,everything())
  list(bsTGI_r=bsTGI_r,bsTGI_df=bsTGI_df)
}

#' Calculate normalized AUC for tumor growth curve
#'
#' @param df tumor growth data for a single mouse
#'
#' @return an AUC dataframe
#' @export
#'
#' @examples get_AUC(mouse_tv)
get_AUC=function(df){
  df=as.data.frame(df)
  df=df[order(df$Day),]
  day_v=df$Day
  TV_v=log(df$TV+1)
  n1=length(day_v)
  AUC=NA
  if(n1>1){
    day_span=day_v[n1]-day_v[1]
    AUC=pracma::trapz(x=day_v,y=TV_v) - pracma::trapz(x=day_v,y=rep(TV_v[1],length(day_v)))
    AUC = 2*AUC/(day_span*day_span) #tumor growth rate under exponential growth model
  }
  auc_df=data.frame('Mouse'=unique(df$Mouse),'AUC'=AUC,'Day'=day_v[n1],'Group'=unique(as.character(df$Group)))
  auc_df
}


#' Calculate medianAUCRatio for two groups
#'
#' @param df AUC data frame
#' @param grp grp id
#' @param ref_grp reference group id
#'
#' @return A number
mAUCr=function(df,grp,ref_grp){
  auc_v=df[df$Group==ref_grp,'AUC']
  auc_t=df[df$Group==grp,'AUC']
  aucs=expand.grid('AUC_tr'=auc_t,'AUC_v'=auc_v)
  aucs=aucs[complete.cases(aucs),]
  medianAUCRatio=with(aucs,median(AUC_tr/AUC_v))
  medianAUCRatio
}

#' Calculate medianAUCRatio statistic for bootstrap
#'
#' @param auc_df AUC data frame for each mouse
#' @param ref_group reference group
#' @param idx bootstrap index
#'
#' @return A vector of medianAUCRatio statistic for each group
bs_mAUCr=function(auc_df,ref_group="Group 1",idx){
  auc_df=auc_df[idx,]
  grps=unique(auc_df$Group)
  grps=grps[grps != ref_group]
  mAUCrs=sapply(grps,function(g) mAUCr(auc_df,g,ref_group))
  names(mAUCrs)=grps
  mAUCrs
}



#' Calculate median AUC ratio using bootstrap
#'
#' @param tv tumor growth data
#' @param sel_day the last day selected for tumor growth data, if not defined, use all data
#' @param ref_group reference group
#' @param ci confidence interval
#' @param ci_type type of bootstrap confidence interval, can be stud,perc or bca
#' @param nrep number of bootstrap replcates
#'
#' @return A list of boostrap mAUCr and mAUCr dataframe
#' @export
#'
#' @examples get_mAUCr(LS_1034,ci=0.9,ci_type='bca')
get_mAUCr=function(tv,sel_day=NA,ref_group='Group 1',ci=0.95,ci_type='perc',nrep=1000){
  if(!is.na(sel_day)) tv=tv %>% filter(Day <= sel_day)
  mouses=split(tv,tv$Mouse)
  auc_mouse=do.call(rbind,lapply(mouses,get_AUC)) %>% as.data.frame()
  auc_mouse=auc_mouse %>% mutate(Group = factor(Group))
  mAUCr_r=boot::boot(data=auc_mouse,statistic=bs_mAUCr,ref_group=ref_group,strata=auc_mouse$Group,R=nrep)
  group_n=length(mAUCr_r$t0)
  cis=do.call(rbind,lapply(1:group_n,function(i) getCI(mAUCr_r,i,conf=ci,ci_type=ci_type)))
  bsAUC_df=cbind(broom::tidy(mAUCr_r),cis)
  bsAUC_df=bsAUC_df[,-3]
  names(bsAUC_df)=c('Group','mAUCr','std.err','lb','ub')
  group_info=tv %>% select(Group,Treatment) %>% distinct()
  bsAUC_df=bsAUC_df %>% left_join(group_info) %>% select(Group,Treatment,everything())
  list(mAUCr_r=mAUCr_r,bsAUC_df=bsAUC_df,auc_mouse=auc_mouse)
}

