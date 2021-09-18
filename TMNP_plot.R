# getwd()

library(ggplot2)
library(ggpubr)
library(cowplot)
library("RColorBrewer")
library("genefilter")
library(dplyr)
library(tidyr)
library(tibble)
# library(VennDiagram)
library(pheatmap)


# huangqin_module_res <- list(huangqin_target,huangqin_kegg,huangqin_go,huangqin_pp,huangqin_cell,huangqin_tissue)

# res_all <- huangqin_module_res
# names(res_all) <- c("Target","KEGG","GO","PP","Cell","Tissue")#all the result of module correlation
# res_name <- colnames(huangqin_all_LFC)
# res_name <- "HQ"

#the function for the module 
plot_fun <- function(res_all,Module=1,res_name){
  if (length(res_name) >1) {
    hist_plot <- plot_hist(res_all,Module,res_name)
    corr_plot <- plot_corr(res_all,Module,res_name)
    fq_plot <- plot_fq(res_all,Module,res_name)
    #venn_plot <- plot_venn(res_all,Module)
    heat_plot <- plot_heatmap(res_all,Module,res_name)
    # plot_grid(hist_plot,corr_plot,fq_plot,heat_plot, nrow = 4)
     return(list(fq_plot,heat_plot,hist_plot,corr_plot))
  }else if (length(res_name) == 1){
    heat_plot <- plot_heatmap(res_all,Module,res_name)
    fq_plot <- plot_fq(res_all,Module,res_name)
    # plot_grid(heat_plot,fq_plot, nrow = 2)
    return(list(fq_plot,heat_plot))
  }
}

# plot_res <- plot_fun(res_all = res_all,Module = "Tissue",res_name = res_name)

#correlation
plot_corr <- function(res_all,Module,res_name){
  res_all <- res_all[[Module]]
  res_ncs <- apply(data.frame(res_all[,endsWith(colnames(res_all),"NCS")]),2,as.numeric)
  res_ncs <- data.frame(res_ncs)
  # rownames(res_ncs) <- rownames(res_ncs)
  colnames(res_ncs) <- res_name
  num_p <- combn(res_name,2)
  lis_p <- paste0("p",seq(1:ncol(num_p)))
  for (i in 1:ncol(num_p)) {
    aes_xy <- num_p[,i]
    cor_value <- cor(res_ncs[,c(aes_xy)])
    cor_value <- as.character(round(cor_value[2],2))
    assign(lis_p[i],ggplot(res_ncs, aes_string(x = aes_xy[1],y= aes_xy[2]))+
             geom_point(size = 1) +
             geom_smooth(method=lm , color="red", fill="#69b3a2", se=FALSE,formula = y~x) +
             theme_bw()+
             # +theme()
             theme(axis.text = element_text(size = 18),legend.title = element_text(size = 25),title = element_text(size = 25),
                   panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
             # annotate("text",x=0,y=4,label= cor_value,color="black",size=8)
             stat_cor(data=res_ncs, method = "pearson")#theme_ipsum(grid = TRUE)
    )
  }
  plot_list <- eval(parse(text = paste0("cowplot::plot_grid(",noquote(paste0(noquote(lis_p),",",collapse = "")),"nrow = 1)")))
  return(plot_list)
}


#frequence

plot_fq <- function(res_all,Module,res_name) {
  res_all <- res_all[[Module]]
  res_ncs <- apply(data.frame(res_all[,endsWith(colnames(res_all),"NCS")]),2,as.numeric)
  res_ncs <- data.frame(res_ncs)
  # rownames(res_ncs) <- rownames(res_ncs)
  colnames(res_ncs) <- res_name
  ns <- length(res_name)
  lis_p <- paste0("p",seq(1:ns))
  for (i in 1:ns) {
    aes_x <- colnames(res_ncs)[i] 
    assign(lis_p[i],
           ggplot(res_ncs,aes_string(x= aes_x))+geom_histogram(aes(y=..density..),colour="white",fill="lightblue")+
             theme_bw()+theme(title = element_text(size = 17),axis.text = element_text(size = 18),axis.title = element_text(size = 25))+labs(x= paste(res_name[i],"_NCS"))  
    )
  }
  plot_list <- eval(parse(text = paste0("cowplot::plot_grid(",noquote(paste0(noquote(lis_p),",",collapse = "")),"nrow = 1)")))
  return(plot_list)
}


# plot_fq(res_all = res_all,Module = "Tissue")

#pheatmap

plot_heatmap <- function(res_all,Module,res_name){
  res_all <- res_all[[Module]]
  res_ncs <- apply(data.frame(res_all[,endsWith(colnames(res_all),"NCS")]),2,as.numeric)
  rownames(res_ncs) <- rownames(res_all)
  colnames(res_ncs) <- res_name
  mat <- data.frame(res_ncs)
  p1 <- pheatmap(mat, 
                 cluster_col=FALSE,cluster_rows = F, show_rownames = F,cellwidth = 80,
                 fontsize = 20,fontsize_row = 15,fontsize_col = 30)
  
  return(p1)
}



#hist 

plot_hist <- function(res_all,Module,res_name){
  res_all <- res_all[[Module]]
  res_all <- data.frame(res_all)
  res_ncs_fdr <- cbind(apply(data.frame(res_all[,endsWith(colnames(res_all),"NCS")]),2,as.numeric),
                       apply(data.frame(res_all[,endsWith(colnames(res_all),"FDR")]),2,as.numeric))
  rownames(res_ncs_fdr) <- rownames(res_all) 
  ns <- length(res_name)
  res_p_count <- matrix(ncol = ns,nrow = 2)
  for (i in 1:ns) {
    res_p_count[1,i] <- length(rownames(res_ncs_fdr)[which(res_ncs_fdr[,ns+i] < 0.05 & res_ncs_fdr[,i]>=0)])
    res_p_count[2,i] <- length(rownames(res_ncs_fdr)[which(res_ncs_fdr[,ns+i] < 0.05 & res_ncs_fdr[,i]<0)])
  }
  colnames(res_p_count) <- res_name
  rownames(res_p_count) <- c("up","down")
  
  dat_hist <- data.frame(res_p_count)%>%
    rownames_to_column(var = "regulate")%>%
    gather(key = "group",value = "count",-1)
  factor <- factor(dat_hist$regulate,levels = c("up","down"))
  p1 <- ggplot(dat_hist,aes(x = group,y=count,fill = factor))+
    geom_col(show.legend = T,width = 0.5) +
    xlab('histgram of result')+
    theme_bw()+
    theme(axis.text = element_text(size = 25),
          legend.text = element_text(size = 25),legend.title = element_text(size = 25),
          strip.background = element_rect(fill = NULL),
          title = element_text(size = 25))
  return(p1)
}

text_desc_fun <- function(res_all,Module,res_name){
  res_module <- res_all
  res_ncs_fdr <- cbind(apply(data.frame(res_module[,endsWith(colnames(res_module),"NCS")]),2,as.numeric),
                       # apply(data.frame(res_module[,endsWith(colnames(res_module),"Pvalue")]),2,as.numeric),
                       apply(data.frame(res_module[,endsWith(colnames(res_module),"FDR")]),2,as.numeric))
  # rownames(res_ncs_fdr) <- rownames(res_module) 
  ns <- length(res_name)
  # res_p_count <- matrix(ncol = ns,nrow = 2)
  text_desc <- list()
  module_lis <- c("targets","pathways","biological Processes","pathological processes","cells","tissues")
  for (i in 1:ns) {
    num1 <- sum(res_ncs_fdr[,ns+i] < 0.05 & res_ncs_fdr[,i]>=0)
    num2 <- sum(res_ncs_fdr[,ns+i] < 0.05 & res_ncs_fdr[,i]<0)
    text_desc[[i]] <- paste(res_name[i],"is postively and negatively realted to",num1,"and",num2,
                            module_lis[Module],"at FDR <0.05, respectively",sep = " ")
  }
  
  return(text_desc)
}









