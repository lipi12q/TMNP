#go kegg geneset ###########
load("KEGG_geneset.Rdata")
load("GO_geneset.Rdata")
#cell pp ##########
load("cell_marker_human.Rdata")
load("data/pp_markers.Rdata")
#tissue module######
load("tissue_tar_module_dn.RData")
load("tissue_tar_module_up.RData")
#target module######
load("target_module_dn_symbol_3275.RData")
load("data/target_module_up_symbol_3275.RData")


fun_all <- function(querydata){
  result_target<- fun_target(querydata)
  result_pathway<-fun_Pathway(querydata)
  result_go<-fun_GO(querydata)
  result_pp<- fun_pp(querydata)
  # result_sym<-fun_sym(querydata)
  result_cell<-fun_cell(querydata)
  result_tissue<-fun_tissue(querydata)
  return(list(result_target[[1]],result_pathway[[1]],result_go[[1]],result_pp[[1]],result_cell[[1]],result_tissue[[1]]))
}



###########  
fun_target <- function(querydata){
  nsets<- length(target_module_dn_symbol[,1])
  genenames<- rownames(querydata)
  n_sample <- ncol(querydata)
  A <-apply(querydata, 2,  function(x){order(x,decreasing = TRUE)})  
  TES <- matrix(ncol = n_sample,nrow = nsets)
  
  for(q in 1:nsets){
    up_gene_order <- match(target_module_up_symbol[q,], genenames)
    up_gene_order <- up_gene_order[!is.na(up_gene_order)]
    dw_gene_order <- match(target_module_dn_symbol[q,], genenames)
    dw_gene_order <- dw_gene_order[!is.na(dw_gene_order)]
    for (p in 1:n_sample) {
      gene.list2 <- A[,p]
      upES <- unlist(GSEA.EnrichmentScore2(gene.list = gene.list2, gene.set = up_gene_order, weighted.score.type = 0))
      downES <- unlist(GSEA.EnrichmentScore2(gene.list = gene.list2, gene.set = dw_gene_order, weighted.score.type = 0))
      TES[q,p] <- upES - downES
    }
  }
  query.rowname <- rownames(target_module_dn_symbol)[!is.na(TES[,1])]
  TES <- matrix(TES[!is.na(TES[,1]),])
  
  ###2、random TESvalue
  rows <- length(genenames)
  nperm <- 10
  phi <- matrix(ncol = 10,nrow = nsets)
  for(q in 1:nsets){
    up_gene_order <- match(target_module_up_symbol[q,], genenames)
    up_gene_order <- up_gene_order[!is.na(up_gene_order)]
    dw_gene_order <- match(target_module_dn_symbol[q,], genenames)
    dw_gene_order <- dw_gene_order[!is.na(dw_gene_order)]
    for (p in 1:10) {
      reshuffled.gene.labels <- sample(1:rows)#order
      
      upES <- unlist(GSEA.EnrichmentScore2(gene.list = reshuffled.gene.labels, gene.set = up_gene_order, weighted.score.type = 0))
      downES <- unlist(GSEA.EnrichmentScore2(gene.list = reshuffled.gene.labels, gene.set = dw_gene_order, weighted.score.type = 0))
      phi[q,p] <- upES - downES
    }
  }  
  phi <- phi[!is.na(phi[,1]),]
  nsets1 <- length(phi[,1])
  ###3. p value  
  pvalue <- matrix(nrow = nsets1,ncol = n_sample)
  for(i in 1:n_sample){
    for (j in 1:nsets1) {
      pvalue[j,i] <- ifelse(TES[j,i] >= 0, sum(phi[j,] >= TES[j,i])/sum(phi[j,]>=0), 
                            sum(phi[j,] < TES[j,i])/sum(phi[j,]<0))
    } 
  }
  ###3、NESvalue and fdr value
  # 
  phi.pos.mean <- apply(phi,1,function(x){
    mean(x[x>=0])
  })
  phi.neg.mean <- apply(phi,1,function(x){
    mean(abs(x[x<0]))
  })
  
  phi.norm.pos <- phi
  phi.norm.pos[phi<0]<-0
  phi.norm.pos <-phi.norm.pos/phi.pos.mean
  phi.norm.neg <- phi
  phi.norm.neg[phi>=0]<-0
  phi.norm.neg <-phi.norm.neg/phi.neg.mean
  phi.norm <- phi
  phi.norm[phi<0] <- phi.norm.neg[phi<0]
  phi.norm[phi>=0] <- phi.norm.pos[phi>=0]
  #r#r#r#r#r#r#r#r#r#r#r#r#r#r
  ES_query <- TES
  ES_query.pos <- ES_query
  ES_query.pos[ES_query<0]<-0
  ES_query.pos <-ES_query.pos/phi.pos.mean
  ES_query.neg <- ES_query
  ES_query.neg[ES_query>=0]<-0
  ES_query.neg <-ES_query.neg/phi.neg.mean
  
  ES_query.norm <- ES_query
  ES_query.norm[ES_query<0] <- ES_query.neg[ES_query<0]
  ES_query.norm[ES_query>=0] <- ES_query.pos[ES_query>=0]
  
  FDRvalue <- matrix(nrow = nsets1,ncol = n_sample)
  
  for(i in 1:n_sample){
    for (j in 1:nsets1) {
      if(ES_query.norm[j,i] >= 0){
        A <- sum(phi.norm >= ES_query.norm[j,i])/sum(phi.norm>=0)
        B <- sum(ES_query.norm >= ES_query.norm[j,i])/sum(ES_query.norm>=0)
        FDRvalue[j,i] <- A/B 
      }
      else{
        A <- sum(phi.norm <= ES_query.norm[j,i])/sum(phi.norm<0)
        B <- sum(ES_query.norm <= ES_query.norm[j,i])/sum(ES_query.norm<0)
        FDRvalue[j,i] <- A/B
      } 
    }                   
  }
  
  query_target_gsea <- cbind(query.rowname,TES,ES_query.norm, pvalue, FDRvalue)
  rownames(query_target_gsea) <- query.rowname
  colnames(query_target_gsea) <- c("Target name",paste(rep(colnames(querydata),times=4),rep(c("_CS","_NCS","_Pvalue","_FDR"),each=n_sample),sep = ""))
  gc()
  return(list(query_target_gsea=query_target_gsea))
}
############

#pathway GSEA analysis####
fun_Pathway <- function(querydata){
  gene.sybmol<- rownames(querydata) #r
  max.Ng <- length(KEGG_geneset$description)
  temp.size.G <- lengths(KEGG_geneset$hsa)
  max.size.G <- max(temp.size.G)      
  kegg_geneset <- matrix(rep("null", max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
  temp.names <- vector(length = max.Ng, mode = "character")
  temp.desc <- vector(length = max.Ng, mode = "character")
  gs.count <- 1
  gene.set.name <- names(KEGG_geneset$hsa)
  gene.set.desc <- KEGG_geneset$description
  for (i in 1:max.Ng) {
    gene.set.size <- length(KEGG_geneset$hsa[[i]])
    
    existing.set <- is.element(KEGG_geneset$hsa[[i]], gene.sybmol)
    set.size <- length(existing.set[existing.set == T])
    if ((set.size < 15) || (set.size > 2000)) next
    temp.size.G[gs.count] <- set.size
    kegg_geneset[gs.count,] <- c(KEGG_geneset$hsa[[i]][existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
    temp.names[gs.count] <- gene.set.name[i]
    temp.desc[gs.count] <- gene.set.desc[i]
    gs.count <- gs.count + 1
  } 
  Ng <- gs.count - 1#set length
  gs.names <- vector(length = Ng, mode = "character")#ES_query_rowname
  gs.desc <- vector(length = Ng, mode = "character")#description
  gs.names <- temp.names[1:Ng]
  gs.desc <- temp.desc[1:Ng] 
  
  n_sample <- ncol(querydata) #r
  nsets <- Ng
  A_query <- apply(querydata, 2, function(x){order(x, decreasing = TRUE)}) #r
  ES_query <- matrix(nrow = nsets,ncol = n_sample) 
  
  for(q in 1:nsets){
    gene_set <- kegg_geneset[q,kegg_geneset[q,] != "null"]
    gene_order <- vector(length = length(gene_set),mode = "numeric")
    gene_order <- match(gene_set, gene.sybmol)
    for(p in 1:n_sample) {
      gene.list2 <- A_query[, p]
      ES_query[q, p] <- unlist(GSEA.EnrichmentScore2(gene.list = gene.list2, gene.set = gene_order, weighted.score.type = 0))
    }
  }
  
  # nperm=10
  
  rows=length(gene.sybmol)
  nperm=10
  phi <- matrix(nrow = nsets,ncol = 10)
  
  for(q in 1:nsets){
    gene_set <- kegg_geneset[q,kegg_geneset[q,] != "null"]
    gene_order <- vector(length = length(gene_set),mode = "numeric")
    gene_order <- match(gene_set, gene.sybmol)
    for (r in 1:nperm) {
      reshuffled.gene.labels <- sample(1:rows)
      gsea.result <- GSEA.EnrichmentScore2(gene.list = reshuffled.gene.labels, gene.set = gene_order, weighted.score.type = 0)
      phi[q,r] <- gsea.result$ES
    }
  }
  
 
  pvalue <- matrix(nrow = nsets,ncol = n_sample)
  for(i in 1:n_sample){
    for (j in 1:nsets) {
      
      pvalue[j,i] <- ifelse(ES_query[j,i] >= 0, sum(phi[j,] >= ES_query[j,i])/sum(phi[j,]>=0), 
                            sum(phi[j,] < ES_query[j,i])/sum(phi[j,]<0))
    } 
  }  
  
  
  # 
  phi.pos.mean <- apply(phi,1,function(x){
    mean(x[x>=0])
  })
  phi.neg.mean <- apply(phi,1,function(x){
    mean(abs(x[x<0]))
  })
  
  phi.norm.pos <- phi
  phi.norm.pos[phi<0]<-0
  phi.norm.pos <-phi.norm.pos/phi.pos.mean
  phi.norm.neg <- phi
  phi.norm.neg[phi>=0]<-0
  phi.norm.neg <-phi.norm.neg/phi.neg.mean
  
  phi.norm <- phi
  phi.norm[phi<0] <- phi.norm.neg[phi<0]
  phi.norm[phi>=0] <- phi.norm.pos[phi>=0]
  #r#r#r#r#r#r#r#r#r#r#r#r#r#r
  ES_query.pos <- ES_query
  ES_query.pos[ES_query<0]<-0
  ES_query.pos <-ES_query.pos/phi.pos.mean
  ES_query.neg <- ES_query
  ES_query.neg[ES_query>=0]<-0
  ES_query.neg <-ES_query.neg/phi.neg.mean
  
  ES_query.norm <- ES_query
  ES_query.norm[ES_query<0] <- ES_query.neg[ES_query<0]
  ES_query.norm[ES_query>=0] <- ES_query.pos[ES_query>=0]
  
  FDRvalue <- matrix(nrow = nsets,ncol = n_sample)
  
  for(i in 1:n_sample){
    for (j in 1:nsets) {
      if(ES_query.norm[j,i] >= 0){
        A <- sum(phi.norm >= ES_query.norm[j,i])/sum(phi.norm>=0)
        B <- sum(ES_query.norm >= ES_query.norm[j,i])/sum(ES_query.norm>=0)
        FDRvalue[j,i] <- A/B 
      }
      else{
        A <- sum(phi.norm <= ES_query.norm[j,i])/sum(phi.norm<0)
        B <- sum(ES_query.norm <= ES_query.norm[j,i])/sum(ES_query.norm<0)
        FDRvalue[j,i] <- A/B
      } 
    }                   
  }
  
  query_kegg_gsea <- cbind(gs.desc,ES_query,ES_query.norm, pvalue, FDRvalue)
  rownames(query_kegg_gsea) <- gs.names
  colnames(query_kegg_gsea) <- c("description",paste(rep(colnames(querydata),times=4),rep(c("_CS","_NCS", "_Pvalue","_FDR"),each=n_sample),sep = ""))
  return(list(query_kegg_gsea=query_kegg_gsea))
  gc()
  
}
###########

#go GSEA analysis####
fun_GO <- function(querydata){
  gene.sybmol<- rownames(querydata) #r
  max.Ng <- length(GO_geneset$description)
  temp.size.G <- lengths(GO_geneset$go_bp)
  max.size.G <- max(temp.size.G)      
  go_geneset <- matrix(rep("null", max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
  temp.names <- vector(length = max.Ng, mode = "character")
  temp.desc <- vector(length = max.Ng, mode = "character")
  gs.count <- 1
  gene.set.name <- names(GO_geneset$go_bp)
  gene.set.desc <- GO_geneset$description
  for (i in 1:max.Ng) {
    gene.set.size <- length(GO_geneset$go_bp[[i]])
    existing.set <- is.element(GO_geneset$go_bp[[i]], gene.sybmol)
    set.size <- length(existing.set[existing.set == T])
    if ((set.size < 15) || (set.size > 2000)) next
    temp.size.G[gs.count] <- set.size
    go_geneset[gs.count,] <- c(GO_geneset$go_bp[[i]][existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
    temp.names[gs.count] <- gene.set.name[i]
    temp.desc[gs.count] <- gene.set.desc[i]
    gs.count <- gs.count + 1
  } 
  Ng <- gs.count - 1#set length
  gs.names <- vector(length = Ng, mode = "character")#ES_query_rowname
  gs.desc <- vector(length = Ng, mode = "character")#description
  gs.names <- temp.names[1:Ng]
  gs.desc <- temp.desc[1:Ng] 
  
  n_sample <- ncol(querydata) #r
  nsets <- Ng
  A_query <- apply(querydata, 2, function(x){order(x, decreasing = TRUE)}) #r
  ES_query <- matrix(nrow = nsets,ncol = n_sample) 
  
  for(q in 1:nsets){
    gene_set <- go_geneset[q,go_geneset[q,] != "null"]
    gene_order <- vector(length = length(gene_set),mode = "numeric")
    gene_order <- match(gene_set, gene.sybmol)
    for(p in 1:n_sample) {
      gene.list2 <- A_query[, p]
      ES_query[q, p] <- unlist(GSEA.EnrichmentScore2(gene.list = gene.list2, gene.set = gene_order, weighted.score.type = 0))
    }
  }
  
  # nperm=10
  
  rows=length(gene.sybmol)
  nperm=10
  phi <- matrix(nrow = nsets,ncol = 10)
  
  for(q in 1:nsets){
    gene_set <- go_geneset[q,go_geneset[q,] != "null"]
    gene_order <- vector(length = length(gene_set),mode = "numeric")
    gene_order <- match(gene_set, gene.sybmol)
    for (r in 1:nperm) {
      reshuffled.gene.labels <- sample(1:rows)
      gsea.result <- GSEA.EnrichmentScore2(gene.list = reshuffled.gene.labels, gene.set = gene_order, weighted.score.type = 0)
      phi[q,r] <- gsea.result$ES
    }
  }
  

  
  pvalue <- matrix(nrow = nsets,ncol = n_sample)
  for(i in 1:n_sample){
    for (j in 1:nsets) {
      pvalue[j,i] <- ifelse(ES_query[j,i] >= 0, sum(phi[j,] >= ES_query[j,i])/sum(phi[j,]>=0), 
                            sum(phi[j,] < ES_query[j,i])/sum(phi[j,]<0))
    } 
  }  
  
  
  #
  phi.pos.mean <- apply(phi,1,function(x){
    mean(x[x>=0])
  })
  phi.neg.mean <- apply(phi,1,function(x){
    mean(abs(x[x<0]))
  })
  
  phi.norm.pos <- phi
  phi.norm.pos[phi<0]<-0
  phi.norm.pos <-phi.norm.pos/phi.pos.mean
  phi.norm.neg <- phi
  phi.norm.neg[phi>=0]<-0
  phi.norm.neg <-phi.norm.neg/phi.neg.mean
  
  phi.norm <- phi
  phi.norm[phi<0] <- phi.norm.neg[phi<0]
  phi.norm[phi>=0] <- phi.norm.pos[phi>=0]
  #r#r#r#r#r#r#r#r#r#r#r#r#r#r
  ES_query.pos <- ES_query
  ES_query.pos[ES_query<0]<-0
  ES_query.pos <-ES_query.pos/phi.pos.mean
  ES_query.neg <- ES_query
  ES_query.neg[ES_query>=0]<-0
  ES_query.neg <-ES_query.neg/phi.neg.mean
  
  ES_query.norm <- ES_query
  ES_query.norm[ES_query<0] <- ES_query.neg[ES_query<0]
  ES_query.norm[ES_query>=0] <- ES_query.pos[ES_query>=0]
  
  FDRvalue <- matrix(nrow = nsets,ncol = n_sample)
  
  for(i in 1:n_sample){
    for (j in 1:nsets) {
      if(ES_query.norm[j,i] >= 0){
        A <- sum(phi.norm >= ES_query.norm[j,i])/sum(phi.norm>=0)
        B <- sum(ES_query.norm >= ES_query.norm[j,i])/sum(ES_query.norm>=0)
        FDRvalue[j,i] <- A/B 
      }
      else{
        A <- sum(phi.norm <= ES_query.norm[j,i])/sum(phi.norm<0)
        B <- sum(ES_query.norm <= ES_query.norm[j,i])/sum(ES_query.norm<0)
        FDRvalue[j,i] <- A/B
      } 
    }                   
  }
  
  query_go_gsea <- cbind(gs.desc,ES_query,ES_query.norm, pvalue, FDRvalue)
  rownames(query_go_gsea) <- gs.names
  colnames(query_go_gsea) <- c("description",paste(rep(colnames(querydata),times=4),rep(c("_CS","_NCS", "_Pvalue","_FDR"),each=n_sample),sep = ""))
  return(list(query_go_gsea=query_go_gsea))
  gc()
  
}
###########

##pathological process GSEA analysis####
fun_pp <- function(querydata) {
  nsets<-length(pp_marker_use[,1])
  gene.sybmol<- rownames(querydata)
  n_sample <- ncol(querydata)
  
  ES_query <- matrix(nrow = nsets,ncol = n_sample) 
  A_query <- apply(querydata, 2, function(x){order(x, decreasing = TRUE)})
  
  
  for(q in 1:nsets){
    up_gene_order <- match(pp_marker_use[q,], gene.sybmol)
    up_gene_order <- up_gene_order[!is.na(up_gene_order)]
    for(p in 1:n_sample) {
      gene.list2 <- A_query[, p]
      
      ES_query[q, p] <- unlist(GSEA.EnrichmentScore2(gene.list = gene.list2, gene.set = up_gene_order, weighted.score.type = 0))
    }
  }
  query.rowname <- rownames(pp_marker_use)[!is.na(ES_query[,1])]
  ES_query <- matrix(ES_query[!is.na(ES_query[,1]),])
  # compute enrichment for random permutations####
  rows=length(gene.sybmol)
  nperm=10
  phi <- matrix(nrow = nsets,ncol = 10)
  
  for(q in 1:nsets){
    up_gene_order <- match(pp_marker_use[q,], gene.sybmol)
    up_gene_order <- up_gene_order[!is.na(up_gene_order)]
    for (r in 1:nperm) {
      reshuffled.gene.labels <- sample(1:rows)
      gsea.result <- GSEA.EnrichmentScore2(gene.list = reshuffled.gene.labels, gene.set = up_gene_order, weighted.score.type = 0)
      phi[q,r] <- gsea.result$ES
    }
  }
  phi <- phi[!is.na(phi[,1]),]
  nsets1 <- length(phi[,1])
  
 
  
  pvalue <- matrix(nrow = nsets1,ncol = n_sample)
  for(i in 1:n_sample){
    for (j in 1:nsets1) {
      
      pvalue[j,i] <- ifelse(ES_query[j,i] >= 0, sum(phi[j,] >= ES_query[j,i])/sum(phi[j,]>=0), 
                            sum(phi[j,] < ES_query[j,i])/sum(phi[j,]<0))
    } 
  }
  
  # 
  phi.pos.mean <- apply(phi,1,function(x){
    mean(x[x>=0])
  })
  phi.neg.mean <- apply(phi,1,function(x){
    mean(abs(x[x<0]))
  })
  
  phi.norm.pos <- phi
  phi.norm.pos[phi<0]<-0
  phi.norm.pos <-phi.norm.pos/phi.pos.mean
  phi.norm.neg <- phi
  phi.norm.neg[phi>=0]<-0
  phi.norm.neg <-phi.norm.neg/phi.neg.mean
  
  phi.norm <- phi
  phi.norm[phi<0] <- phi.norm.neg[phi<0]
  phi.norm[phi>=0] <- phi.norm.pos[phi>=0]
  
  ES_query.pos <- ES_query
  ES_query.pos[ES_query<0]<-0
  ES_query.pos <-ES_query.pos/phi.pos.mean
  ES_query.neg <- ES_query
  ES_query.neg[ES_query>=0]<-0
  ES_query.neg <-ES_query.neg/phi.neg.mean
  
  ES_query.norm <- ES_query
  ES_query.norm[ES_query<0] <- ES_query.neg[ES_query<0]
  ES_query.norm[ES_query>=0] <- ES_query.pos[ES_query>=0]
  
  FDRvalue <- matrix(nrow = nsets1,ncol = n_sample)
  
  for(i in 1:n_sample){
    for (j in 1:nsets1) {
      if(ES_query.norm[j,i] >= 0){
        A <- sum(phi.norm >= ES_query.norm[j,i])/sum(phi.norm>=0)
        B <- sum(ES_query.norm >= ES_query.norm[j,i])/sum(ES_query.norm>=0)
        FDRvalue[j,i] <- A/B 
      }
      else{
        A <- sum(phi.norm <= ES_query.norm[j,i])/sum(phi.norm<0)
        B <- sum(ES_query.norm <= ES_query.norm[j,i])/sum(ES_query.norm<0)
        FDRvalue[j,i] <- A/B
      } 
    }                   
  }
  
  query_pp_gsea <- cbind(pp_info[match(query.rowname,rownames(pp_marker_use))],ES_query,ES_query.norm, pvalue, FDRvalue)
  rownames(query_pp_gsea) <- query.rowname
  colnames(query_pp_gsea) <- c("Pathological process",paste(rep(colnames(querydata),times=4),rep(c("_CS","_NCS", "_Pvalue","_FDR"),each=n_sample),sep = ""))
  gc()
  return(list(query_pp_gsea=query_pp_gsea))
}
########

##cell GSEA analysis####
fun_cell <- function(querydata) {
  gene.sybmol<- rownames(querydata)
  n_sample <- ncol(querydata)
  nsets <- length(cell_marker_human_use[,1])
  ES_query <- matrix(nrow = nsets, ncol = n_sample) 
  A_query <- apply(querydata, 2, function(x){order(x, decreasing = TRUE)})#
  
  
  for(q in 1:nsets){
    up_gene_order <- match(cell_marker_human_use[q,], gene.sybmol)
    up_gene_order <- up_gene_order[!is.na(up_gene_order)]
    for(p in 1:n_sample) {
      gene.list2 <- A_query[, p]
      ES_query[q, p] <- unlist(GSEA.EnrichmentScore2(gene.list = gene.list2, gene.set = up_gene_order, weighted.score.type = 0))
    }
  }
  
  query.rowname <- rownames(cell_marker_human_use)[!is.na(ES_query[,1])]
  ES_query <- matrix(ES_query[!is.na(ES_query[,1]),])
  
  # compute enrichment for random permutations####
  rows=length(gene.sybmol)
  nperm=10
  phi <- matrix(nrow = length(cell_marker_human_use[,1]),ncol = 10)
  
  for(q in 1:nsets){
    up_gene_order <- match(cell_marker_human_use[q,], gene.sybmol)
    up_gene_order <- up_gene_order[!is.na(up_gene_order)]
    for (r in 1:nperm) {
      reshuffled.gene.labels <- sample(1:rows)
      gsea.result <- GSEA.EnrichmentScore2(gene.list = reshuffled.gene.labels, gene.set = up_gene_order, weighted.score.type = 0)
      phi[q,r] <- gsea.result$ES
    }
  }
  
  phi <- phi[!is.na(phi[,1]),]
  nsets1 <- length(phi[,1])
 
  
  pvalue <- matrix(nrow = nsets1,ncol = n_sample)
  for(i in 1:n_sample){
    for (j in 1:nsets1) {
      
      pvalue[j,i] <- ifelse(ES_query[j,i] >= 0, sum(phi[j,] >= ES_query[j,i])/sum(phi[j,]>=0), 
                            sum(phi[j,] < ES_query[j,i])/sum(phi[j,]<0))
    } 
  }
  
  # 
  phi.pos.mean <- apply(phi,1,function(x){
    mean(x[x>=0])
  })
  phi.neg.mean <- apply(phi,1,function(x){
    mean(abs(x[x<0]))
  })
  
  phi.norm.pos <- phi
  phi.norm.pos[phi<0]<-0
  phi.norm.pos <- phi.norm.pos/phi.pos.mean
  phi.norm.neg <- phi
  phi.norm.neg[phi>=0]<-0
  phi.norm.neg <-phi.norm.neg/phi.neg.mean
  
  phi.norm <- phi
  phi.norm[phi<0] <- phi.norm.neg[phi<0]
  phi.norm[phi>=0] <- phi.norm.pos[phi>=0]
  
  ES_query.pos <- ES_query
  ES_query.pos[ES_query<0]<-0
  ES_query.pos <-ES_query.pos/phi.pos.mean
  ES_query.neg <- ES_query
  ES_query.neg[ES_query>=0]<-0
  ES_query.neg <-ES_query.neg/phi.neg.mean
  
  ES_query.norm <- ES_query
  ES_query.norm[ES_query<0] <- ES_query.neg[ES_query<0]
  ES_query.norm[ES_query>=0] <- ES_query.pos[ES_query>=0]
  
  FDRvalue <- matrix(nrow = nsets1,ncol = n_sample)
  
  for(i in 1:n_sample){
    for (j in 1:nsets1) {
      if(ES_query.norm[j,i] >= 0){
        A <- sum(phi.norm >= ES_query.norm[j,i])/sum(phi.norm>=0)
        B <- sum(ES_query.norm >= ES_query.norm[j,i])/sum(ES_query.norm>=0)
        FDRvalue[j,i] <- A/B 
      }
      else{
        A <- sum(phi.norm <= ES_query.norm[j,i])/sum(phi.norm<0)
        B <- sum(ES_query.norm <= ES_query.norm[j,i])/sum(ES_query.norm<0)
        FDRvalue[j,i] <- A/B
      } 
    }  
  }
  query_cell_gsea <- cbind(cell_info_hum[match(query.rowname,rownames(cell_marker_human_use)),],ES_query,ES_query.norm, pvalue, FDRvalue)
  rownames(query_cell_gsea) <- query.rowname
  colnames(query_cell_gsea) <- c(colnames(cell_info_hum),paste(rep(colnames(querydata),times=4),rep(c("_CS","_NCS", "_Pvalue","_FDR"),each=n_sample),sep = ""))
  gc()
  return(list(query_cell_gsea=query_cell_gsea))
  
}
###############


#tissue GSEA analysis####
fun_tissue <- function(querydata){
  nsets <- length(tissue_tar_module_dn[,1])
  genenames<- rownames(querydata)
  n_sample <- ncol(querydata)
  A <-apply(querydata, 2,  function(x){order(x,decreasing = TRUE)})  
  # upES <- matrix(ncol = ncol(querydata),nrow = length(tissue_tar_module_dn[,1]))
  # downES <- matrix(ncol = ncol(querydata),nrow = length(tissue_tar_module_dn[,1]))
  TES <- matrix(ncol = n_sample,nrow = nsets)
  
  for(q in 1:nsets){
    up_gene_order <- match(tissue_tar_module_up[q,], genenames)
    up_gene_order <- up_gene_order[!is.na(up_gene_order)]
    dw_gene_order <- match(tissue_tar_module_dn[q,], genenames)
    dw_gene_order <- dw_gene_order[!is.na(dw_gene_order)]
    for (p in 1:n_sample) {
      gene.list2 <- A[,p]
      # upES[q,p] <- unlist(GSEA.EnrichmentScore2(gene.list = gene.list2, gene.set = up_gene_order, weighted.score.type = 0))
      # downES[q,p] <- unlist(GSEA.EnrichmentScore2(gene.list = gene.list2, gene.set = dw_gene_order, weighted.score.type = 0))
      upES <- unlist(GSEA.EnrichmentScore2(gene.list = gene.list2, gene.set = up_gene_order, weighted.score.type = 0))
      downES <- unlist(GSEA.EnrichmentScore2(gene.list = gene.list2, gene.set = dw_gene_order, weighted.score.type = 0))
      TES[q,p] <- upES - downES
    }
  }
  
  ###2、random TESvalue
  rows <- length(genenames)
  nperm <- 10
  phi <- matrix(ncol = 10,nrow = nsets)
  for(q in 1:nsets){
    up_gene_order <- match(tissue_tar_module_up[q,], genenames)
    up_gene_order <- up_gene_order[!is.na(up_gene_order)]
    dw_gene_order <- match(tissue_tar_module_dn[q,], genenames)
    dw_gene_order <- dw_gene_order[!is.na(dw_gene_order)]
    for (p in 1:10) {
      reshuffled.gene.labels <- sample(1:rows)#
      
      upES <- unlist(GSEA.EnrichmentScore2(gene.list = reshuffled.gene.labels, gene.set = up_gene_order, weighted.score.type = 0))
      downES <- unlist(GSEA.EnrichmentScore2(gene.list = reshuffled.gene.labels, gene.set = dw_gene_order, weighted.score.type = 0))
      phi[q,p] <- upES - downES
    }
  }  
  
  ###3、 p value  
  pvalue <- matrix(nrow = nsets,ncol = n_sample)
  for(i in 1:n_sample){
    for (j in 1:nsets) {
      pvalue[j,i] <- ifelse(TES[j,i] >= 0, sum(phi[j,] >= TES[j,i])/sum(phi[j,]>=0), 
                            sum(phi[j,] < TES[j,i])/sum(phi[j,]<0))
    } 
  }
  ###3、NESvalue and fdr value
  # 
  phi.pos.mean <- apply(phi,1,function(x){
    mean(x[x>=0])
  })
  phi.neg.mean <- apply(phi,1,function(x){
    mean(abs(x[x<0]))
  })
  
  phi.norm.pos <- phi
  phi.norm.pos[phi<0]<-0
  phi.norm.pos <-phi.norm.pos/phi.pos.mean
  phi.norm.neg <- phi
  phi.norm.neg[phi>=0]<-0
  phi.norm.neg <-phi.norm.neg/phi.neg.mean
  phi.norm <- phi
  phi.norm[phi<0] <- phi.norm.neg[phi<0]
  phi.norm[phi>=0] <- phi.norm.pos[phi>=0]
  #r#r#r#r#r#r#r#r#r#r#r#r#r#r
  ES_query <- TES
  ES_query.pos <- ES_query
  ES_query.pos[ES_query<0]<-0
  ES_query.pos <-ES_query.pos/phi.pos.mean
  ES_query.neg <- ES_query
  ES_query.neg[ES_query>=0]<-0
  ES_query.neg <-ES_query.neg/phi.neg.mean
  
  ES_query.norm <- ES_query
  ES_query.norm[ES_query<0] <- ES_query.neg[ES_query<0]
  ES_query.norm[ES_query>=0] <- ES_query.pos[ES_query>=0]
  
  FDRvalue <- matrix(nrow = nsets,ncol = n_sample)
  
  for(i in 1:n_sample){
    for (j in 1:nsets) {
      if(ES_query.norm[j,i] >= 0){
        A <- sum(phi.norm >= ES_query.norm[j,i])/sum(phi.norm>=0)
        B <- sum(ES_query.norm >= ES_query.norm[j,i])/sum(ES_query.norm>=0)
        FDRvalue[j,i] <- A/B 
      }
      else{
        A <- sum(phi.norm <= ES_query.norm[j,i])/sum(phi.norm<0)
        B <- sum(ES_query.norm <= ES_query.norm[j,i])/sum(ES_query.norm<0)
        FDRvalue[j,i] <- A/B
      } 
    }                   
  }
  
  query_tissue_gsea <- cbind(rownames(tissue_tar_module_dn),TES,ES_query.norm, pvalue, FDRvalue)
  rownames(query_tissue_gsea) <- rownames(tissue_tar_module_dn)
  colnames(query_tissue_gsea) <- c("Tissue",paste(rep(colnames(querydata),times=4),rep(c("_CS","_NCS", "_Pvalue","_FDR"),each=n_sample),sep = ""))
  gc()
  return(list(query_tissue_gsea=query_tissue_gsea))
}
############



#GSEA.EnrichmentScore2 ####

GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. 
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain") 
  #   RES: Numerical vector containing the running enrichment score for all locations in the gene list 
  #   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
  no.tag.indicator <- 1 - tag.indicator 
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  if (weighted.score.type == 0) {
    correl.vector <- rep(1, N)
  }
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector**alpha)
  sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > - min.ES) {
    #      ES <- max.ES
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } else {
    #      ES <- min.ES
    ES <- signif(min.ES, digits=5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))  
}

GSEA.EnrichmentScore2 <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. It is the same calculation as in 
  # GSEA.EnrichmentScore but faster (x8) without producing the RES, arg.RES and tag.indicator outputs.
  # This call is intended to be used to asses the enrichment of random permutations rather than the 
  # observed one.
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  
  loc.vector <- vector(length=N, mode="numeric")
  peak.res.vector <- vector(length=Nh, mode="numeric")
  valley.res.vector <- vector(length=Nh, mode="numeric")
  tag.correl.vector <- vector(length=Nh, mode="numeric")
  tag.diff.vector <- vector(length=Nh, mode="numeric")
  tag.loc.vector <- vector(length=Nh, mode="numeric")
  
  loc.vector[gene.list] <- seq(1, N)
  tag.loc.vector <- loc.vector[gene.set]
  
  tag.loc.vector <- sort(tag.loc.vector, decreasing = F)
  
  if (weighted.score.type == 0) {
    tag.correl.vector <- rep(1, Nh)
  } else if (weighted.score.type == 1) {
    tag.correl.vector <- correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else if (weighted.score.type == 2) {
    tag.correl.vector <- correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else {
    tag.correl.vector <- correl.vector[tag.loc.vector]**weighted.score.type
    tag.correl.vector <- abs(tag.correl.vector)
  }
  
  norm.tag <- 1.0/sum(tag.correl.vector)
  tag.correl.vector <- tag.correl.vector * norm.tag
  norm.no.tag <- 1.0/Nm
  tag.diff.vector[1] <- (tag.loc.vector[1] - 1) 
  tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
  tag.diff.vector <- tag.diff.vector * norm.no.tag
  peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
  valley.res.vector <- peak.res.vector - tag.correl.vector
  max.ES <- max(peak.res.vector)
  min.ES <- min(valley.res.vector)
  ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
  
  return(list(ES = ES))
}


