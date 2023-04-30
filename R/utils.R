#' Theme for publication
#'
#' @param base_size text base size for ggplot
#' @param base_family text base family
#'
#' @return a theme for publication
#' @import ggplot2
#' @export
#'
#' @examples
#' theme_Publication()
theme_Publication <- function(base_size=14, base_family="helvetica") {
  (ggthemes::theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(2.4), hjust = 0),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1.6)),
            axis.title.y = element_text(angle=90,vjust =2,size = rel(1.4)),
            axis.title.x = element_text(vjust = -0.2,size=rel(1.4)),
            axis.text = element_text(face="bold", size=rel(1.2)),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "top",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold",size=rel(1.4))
    ))

}

#' Scale fill for publication
#'
#' @param ...
#'
#' @return Scale fill for publication
#' @export
#'
#'
#' @examples
#' scale_fill_Publication()
scale_fill_Publication <- function(...){
  discrete_scale("fill","Publication",scales::manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

#' Scale color for publication
#'
#' @param ...
#'
#' @return Scale color for publication
#' @export
#'
#' @examples
#' scale_colour_Publication()
scale_colour_Publication <- function(...){
  discrete_scale("colour","Publication",scales::manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

#' Get bootstrap confidence interval
#'
#' @param res bootstrap results
#' @param conf confidence interval
#' @param ci_type type of bootstrap confidence interval, can be perc,stud or bca
#' @param i index of bootstrap statistic
#'
#' @return boostrap confidence interval
getCI=function(res,i,conf=0.95,ci_type='perc'){
  x=tryCatch({boot::boot.ci(res,conf=conf,index=c(i,i),type=ci_type)},error=function(e){})
  #x=boot::boot.ci(res,conf=conf,index=c(i,i),type=ci_type)
  if(is.null(x)){
    return(rep(NA,2))
  }else{
    return(x[[4]][1,4:5])
  }
}
