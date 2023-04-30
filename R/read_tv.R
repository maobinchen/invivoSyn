#' read_tv
#'
#' @param file the path of input file in a wide format, should have Treatment, Mouse and days in numeric number,
#'        treatment should be in the order of vehicle,group 1, group 2, and combo
#'
#' @return a tumor volume data frame, including TV,RTV and deltaTV
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @import tidyr
#' @import dplyr
#'
#' @examples
#' library(SynergyInVivo)
#' read_tv(system.file("extdata","test.csv",package="SynergyInVivo"))
read_tv=function(file){
  dat=utils::read.csv(file,check.names = F)
  #colnames(dat)[1:2]=c("Treatment","Mouse")
  tv=dat %>% mutate(Treatment=factor(Treatment,levels = unique(Treatment))) %>% mutate(Group=factor(paste0("Group ",as.numeric(Treatment))))
  tv=tv %>% select(Group,Treatment,Mouse,everything())
  tv_long=tv %>% pivot_longer(-1:-3,names_to = 'Day',values_to = 'TV') %>% filter(!is.na(TV))
  tv_long=expand_tv(tv_long)
  tv_long
}

#' Expand tumor volume data frame, add TV0, RTV and logTV
#'
#' @param tv_long tumor volume data frame, including Treatment, Mouse, Day and TV
#'
#' @return a tumor volume data frame, including TV,RTV and deltaTV
#' @export
#'
#' @examples
#' tv_long=expand_tv(tv_long)
expand_tv=function(tv_long){
  tv_long=tv_long %>% mutate(Day=as.numeric(Day))
  tv_long = tv_long %>% left_join(tv_long %>% group_by(Mouse) %>% summarise(Min_day=min(Day))) %>%
    mutate(Day=Day-min(Day))
  tv_long=tv_long %>% select(-Min_day) %>% filter(!is.na(TV))
  sel_day=max(tv_long %>% filter(Group=="Group 1") %>% pull(Day))
  tv_long=tv_long  %>% left_join(tv_long %>% filter(Day==0) %>% select(Mouse,TV) %>% rename(TV0=TV))
  tv_long=tv_long %>% mutate(logTV=log(TV+1),RTV=ifelse(TV0==0,NA,TV/TV0),DeltaTV=TV-TV0)
  tv_long
}
