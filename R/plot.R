#' Draw tumor growth curve
#'
#' @param tv object return from read_tv
#' @param y  variable to plot,can be TV, RTV, DeltaTV or logTV
#' @param save whether or not save the plot
#'
#' @return a ggplot object
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' plot_tumor_growth_curve(LS_1034)
plot_tumor_growth_curve=function(tv,y='TV',save=TRUE,kw='Test'){
  stopifnot("y must be TV,RTV,logTV or DeltaTV"= y %in% c('TV','RTV','logTV','DeltaTV'))
  y=sym(y)
  p=tv %>% ggplot()+aes(Day,!!y,group=Mouse,color=Group)+geom_point()+geom_line()+
    theme_Publication()+facet_wrap(~Treatment)+
    ggsci::scale_color_jco() + theme(axis.title.x=element_text(face = "bold",size = rel(1.6)),axis.text=element_text(size=rel(1.8)),axis.ticks=element_line(size=rel(5)))+
    theme(legend.position = 'none')
  if(save) ggsave(paste0(kw,'_',y,"_tumor_growth_curve.png"),width=10,height=8,dpi=300)
  p
}

#' Draw boxplot to compare group means for each day
#'
#' @param tv object return from read_tv
#' @param y  variable to plot,can be TV, RTV, DeltaTV or logTV
#' @param save whether or not save the plot
#'
#' @import ggplot2
#' @return a ggplot object
#' @export
#'
#' @examples
#' plot_group_by_day(LS_1034,y='RTV')
plot_group_by_day=function(tv,y="TV",save=TRUE){
  stopifnot("y must be TV,RTV,logTV or DeltaTV"= y %in% c('TV','RTV','logTV','DeltaTV'))
  cmps=as.list(as.data.frame(combn(as.character(levels(tv$Group)),2)))
  n_day=length(unique(tv$Day))
  col_n=ceiling(sqrt(n_day))
  row_n=ceiling(n_day/col_n)
  p=ggpubr::ggboxplot(tv,x='Group',y=y,add=c('jitter'),color='Treatment',add.params = list(shape=17),xlab='Day',repel = TRUE)+facet_wrap(~Day,scales = 'free_y')+
    ggpubr::stat_compare_means(method='t.test',comparisons = cmps,aes(label = paste0("p=", ..p.format..)),label.x.npc = 'left',label.y.npc = 'top')+
    theme_Publication()+ggsci::scale_color_npg()
  if(save) ggsave(paste0(y,"~Group_by_day_boxplot.png"),width=col_n*5,height = row_n*4,dpi=300)
  p
}

#' Density plot for selected variable in a data frame
#'
#' @param df data frame
#' @param sel_var variable to plot
#' @param pval p value (optional)
#' @param cutoff
#' @param pe point estimate
#' @param lb lower confidence interval
#' @param ub higher confidence interval
#'
#' @return a ggplot object
#' @export
#'
#' @examples
plot_density=function(df,sel_var,cutoff=NA,pval=NA,pe=NA,lb=NA,ub=NA){
  x1=ifelse(is.na(lb),quantile(df[[sel_var]],0.025),lb)
  x2=ifelse(is.na(ub),quantile(df[[sel_var]],0.975),ub)
  xrange=range(df[[sel_var]])
  sel_var=sym(sel_var)
  pal=ggpubr::get_palette("jco",4)

  p=df %>% ggplot()+aes(!!sel_var) + geom_histogram(aes(y=..density..),size=1.2,color=pal[1],fill='lightgray',alpha=0.5)+
    theme_Publication()+geom_density(fill=pal[3],alpha=0.2,size=1.2,color=pal[4])+ylab("Density")+
    theme(axis.text=element_text(size=rel(1.4)),axis.ticks=element_line(size=rel(2.5)),axis.title = element_text(size=rel(1.2)))
  if(!is.na(pval)) p=p+annotation_custom(grid::textGrob(label = paste0("Bootstrap pval=",pval),
                                                           x = unit(0.83, "npc"), y = unit(0.95, "npc"),
                                                           gp = grid::gpar(cex = 1.5)))
  #if(!is.na(cutoff)) p=p+geom_vline(xintercept = cutoff,size=1.2,linetype=2,color=pal[4])
  if(!is.na(cutoff)) p=p+geom_segment(x = cutoff,xend=cutoff,y=0,yend=10000,size=1.5,linetype=2,color=pal[4])
  #if(!is.na(pe)) p=p+geom_vline(xintercept = pe,size=1.2,linetype=1,color='green')
  if(!is.na(pe)) p=p+geom_point(aes(x=pe,y=0),shape=17,size=rel(6),color='black',fill='black')
  densitys=ggplot_build(p)$data[[1]]$density
  yrange=range(densitys)
  p=p+geom_segment(x=x1,xend=x2,y=0,yend=0,color="blue",size=2,lineend="round")+
    annotation_custom(grid::segmentsGrob(x0 = unit(0.5,"npc"),y0 = unit(0.9,"npc"),x1=unit(0.6,"npc"),y1=unit(0.9,"npc"),gp=grid::gpar(col="blue",lwd=3)))+
    annotation_custom(grid::textGrob(label = "95% Confidence Interval",
                                     x = unit(0.8, "npc"), y = unit(0.9, "npc"),
                                     gp = grid::gpar(cex = 1.5)))
  p
}




