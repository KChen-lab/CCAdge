library(irlba)
library(umap)
library(rdist)
library(EnhancedVolcano)
library(gtools)
library(Seurat)

# color code for plot 

vega_20 = c(
    '#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD',
    '#8C564B', '#E377C2', '#BCBD22', '#17BECF', '#AEC7E8',
    '#FFBB78', '#98DF8A', '#FF9896', '#C5B0D5', '#C49C94',
    '#F7B6D2', '#DBDB8D', '#9EDAE5', '#AD494A', '#8C6D31')

zeileis_28 = c(
    "#023fa5", "#7d87b9", "#bec1d4", "#d6bcc0", "#bb7784", "#8e063b", "#4a6fe3",
    "#8595e1", "#b5bbe3", "#e6afb9", "#e07b91", "#d33f6a", "#11c638", "#8dd593",
    "#c6dec7", "#ead3c6", "#f0b98d", "#ef9708", "#0fcfc0", "#9cded6", "#d5eae7",
    "#f3e1eb", "#f6c4e1", "#f79cd4",
    '#7f7f7f', "#c7c7c7", "#1CE6FF", "#336600")


# comparsion with two clusters


irlCCA <- function(X=NULL,Y=NULL,K=NULL){
  out <- irlba(t(X)%*%Y,nv=K,nu=K)
  x <- X%*%out$u
  y <- Y%*%out$v
  out<-list()
  out$x_dim <- x
  out$y_dim <-y
  out$diff <- (x-y)
  return(out)
}



### permutation to derive null model 

perMu_base <- function(x=NULL){
  mat <- matrix(0, nrow(x),ncol(x))
  vec <- seq(1,ncol(x),1)
  for(i in seq(1,nrow(mat),1)){
    mat[i,] <-  x[i,permute(vec)]
  }
  return(mat)
}


# High-order gene expression change test 

HighOrderScore_singleSample <- function(obj=NULL, 
                       factor=NULL, 
                       reduction=NULL, 
                       sample = NULL, 
                       assay=NULL, 
                       K = 20, 
                       genelst=NULL,
                       T_latent=10^(-5),
                       min.frac= 0.05){
  if(sample!=factor){
    stop("The sample and factor should be the same if in single-sample mode!")
  }
  sample<-factor
  #lowDim <- obj[[reduction]]@cell.embeddings
  compare <- names(table(obj@meta.data[,sample]))
  pos1 <- which(obj@meta.data[,sample]==compare[1])
  pos2 <- which(obj@meta.data[,sample]==compare[2])
  X <- as.matrix(obj[[assay]][,pos1])
  Y <- as.matrix(obj[[assay]][,pos2])
  fct_x <- geneFrac(x=X)
  fct_y <- geneFrac(x=Y)
  gene_lst <- names(fct_x)[which(fct_x>min.frac | fct_y>min.frac)]
  
  message("Single sample mode started!")
  message(paste0("A total of ", length(gene_lst), " genes are used for testing"))
  message("Perform high-order gene change detection...")
  res <- irlCCA(X=X[gene_lst,], Y=Y[gene_lst,], K=K)
  
  # perform permuation to derive the null model of distribution 
  message("Perform permutation to generate the null model...")
  X_permu <- perMu_base(X[gene_lst,])
  Y_permu <- perMu_base(Y[gene_lst,])
  res_permu <- irlCCA(X=X_permu, Y=Y_permu, K=K)
 
  
  # determine effective latent space detection 
  null_mod <- as.vector(as.vector(res_permu$diff[,-c(1)]))
  p_ks <- rep(0,K)
  latent_ks <- rep(0,K)
  for(i in seq(1,K,1)){
    a <- ks.test(res$diff[,i], null_mod)
    if(a$p.value==0){a$p.value<-2.2*10^(-16)}
    p_ks[i] <-a$p.value
    if(a$p.value<=T_latent){
      latent_ks[i] <-1
    }
  }
  
  # using latent variable vectors with p value <10^(-5)
  n_latent <- sum(latent_ks)
  message(paste0(n_latent," informative latent factors were detected "))
  highorder_diff <- res$diff[,which(latent_ks==1)]
  
  # squre of z-score matrix 
  out <- list()
  sd_null_mod <- sd(null_mod)
  mean_null_mod <- mean(null_mod)
  out$diff <- (res$diff - mean_null_mod)/sd_null_mod
  x2_diff <- (highorder_diff-mean_null_mod)^2/(sd_null_mod^2)
  x2_score <- rowSums(x2_diff[,-c(1)])
  #  x2_score follow x2(n) distrubiton 
  p_vec <- pchisq(x2_score, df=n_latent, lower.tail = FALSE)
  out$sum_score<-data.frame("X2_score"=x2_score,"X2_pval"=p_vec, "N"=n_latent,"Diff_depth"=res$diff[,1])
  out$latent_score <- x2_diff
  out$latent_pval <- p_ks
  out$null_mod <- null_mod
  out$x_dim <- res$x_dim
  out$y_dim <- res$y_dim
 
  message("Finished!")
  return(out)
}

# High-order gene expression change test multiple sample 

HighOrderScore_multipleSample <- function(obj=NULL, 
                                          factor=NULL, 
                                          sample=NULL,
                                          reduction=NULL, 
                                          assay=NULL, 
                                          K = 20, 
                                          genelst=NULL,
                                          T_latent=10^(-5),
                                          min.frac=min.frac){
  if(sample==factor){
    stop("The sample and factor should be the different if in multiple-sample mode!")
  }
  
  if(!factor%in%colnames(obj@meta.data)){
    stop("Please make sure the test factor is in meta data")
  }
  factor_lst <- as.vector(unique(obj@meta.data[,factor]))
  sampleTable <- unique(obj@meta.data[,c(factor, sample)])
  colnames(sampleTable)<-c("factor","sample")
  factor_sample <- list()
  for(i in seq(1,length(factor_lst),1)){
    tp <- unique(obj@meta.data[obj@meta.data[,factor]==factor_lst[i],sample])
    n <- length(tp)
    message(paste0("factor ", factor_lst[i], " includes ",n, " samples"))
  }
  
  # perform within group test 
  sampleLst <- unique(obj@meta.data[,sample])
  out <- list()
  for(k in seq(1,length(sampleLst),1)){
    for(l in seq(1, length(sampleLst),1)){
      if(l>k){
          message(paste0("Performing high-order detection between  ",sampleLst[k], " and ",sampleLst[l]))
          id_pair <- paste0(sampleLst[k], "-",sampleLst[l])
          a1 <- sampleTable$factor[sampleTable$sample==sampleLst[k]]
          a2 <- sampleTable$factor[sampleTable$sample==sampleLst[l]] 
          test_pair <- paste0(a1,"-",a2)
          obj_tp <- subset(obj, cells=colnames(obj)[obj@meta.data[,sample]%in%c(sampleLst[k], sampleLst[l])])
          res <- HighOrderScore_singleSample(obj=obj_tp, 
                                             factor = sample, 
                                             sample = sample, 
                                             reduction = NULL,
                                             assay = "RNA", 
                                             K =20,
                                             genelst = NULL,
                                             T_latent=10^(-5),
                                             min.frac = min.frac)
          res$pair <- test_pair
          out[[id_pair]] <- res
          
      }
    }
  }
  message("Performing multiple sample level testing")
  sta_table <- list()
  gene_all <-c()
  test_pair_all <-c()
  id_pair_all <-c()
  for(k in seq(1,length(sampleLst),1)){
    for(l in seq(1, length(sampleLst),1)){
      if(l>k){
        #message(paste0("Performing high-order detection between  ",sampleLst[k], " and ",sampleLst[l]))
        id_pair <- paste0(sampleLst[k], "-",sampleLst[l])
        a1 <- sampleTable$factor[sampleTable$sample==sampleLst[k]]
        a2 <- sampleTable$factor[sampleTable$sample==sampleLst[l]] 
        test_pair <- paste0(a1,"-",a2)
        test_pair_all <-c(test_pair_all, test_pair)
        id_pair_all <-c(id_pair_all, id_pair)
        gene_all <- unique(c(gene_all, rownames(out[[id_pair]]$sum_score)))
        print(test_pair)
      }
    }
  }
  score_mat <- matrix(-1, length(gene_all), length(test_pair_all))
  rownames(score_mat) <- gene_all
  colnames(score_mat) <- test_pair_all
  for(i in seq(1,ncol(score_mat),1)){
    tp <- out[[id_pair_all[i]]]$sum_score
    tp$X2_score <- tp$X2_score/tp$N
    score_mat[rownames(tp),i] <- tp$X2_score
  }
  
  # add multiple-sample level p value calculation 
  pval_mat <- matrix(0, length(gene_all),5)
  colnames(pval_mat) <-c("A1","A2","miss_frq1","miss_frq2","miss_frq3")
  for(i in seq(1,nrow(score_mat),1)){
    pval_mat[i,] <-T_test(x=score_mat[i,], flag=unique(test_pair_all))
  }
  rownames(pval_mat) <- rownames(score_mat)
  pval_mat <- as.data.frame(pval_mat)
  message("Finished")
  return(pval_mat)
}

T_test <- function(x=NULL, flag=NULL){
  #x <- score_mat[148,]
  #flag<-c("sn-sn","sn-sc","sc-sc")
  #x <- x[!x==(-1)]
  grp1 <- x[names(x)==flag[1]]
  grp2 <- x[names(x)==flag[2]]
  grp3 <- x[names(x)==flag[3]]
  p12 <- t.test(grp1, grp2, alternative = c("less"))
  p31 <- t.test(grp3, grp2, alternative = c("less"))
  res <- c(p12$p.value, p31$p.value, sum(grp1==(-1))/length(grp1),sum(grp2==(-1))/length(grp2), sum(grp3==(-1))/length(grp3))
  return(res)
  
}

# Main function 

HighOrder_Test <- function(obj=NULL, 
                           factor=NULL, 
                           sample=NULL,
                           reduction=NULL, 
                           assay=NULL, 
                           K = 20, 
                           genelst=NULL,
                           T_latent=10^(-5),
                           mode="single",
                           min.frac = 0.05){
  if(!factor%in%colnames(obj@meta.data)){
    stop("Please make sure the test factor is listed in meta data!")
  }
  if(!sample%in%colnames(obj@meta.data)){
    stop("Please make sure the sample factor is listed in meta data!")
  }
  if(!mode %in%c("single","multiple")){
    stop("The mode can be only single or multiple!")
  }
  if(is.null(genelst)){genelst <- rownames(obj)}else{
    if(!sum(genelst%in%rownames(obj))==length(genelst)){
      stop("Some genes are not listed in seurat object!")
    }
  }
  if(mode=="single"){
    if(factor!=sample){
      stop("Please make sure the sample and factor be the same in single sample mode!")
    }
    out <- HighOrderScore_singleSample(
            obj=obj, 
            factor=factor, 
            reduction=NULL, 
            sample = sample, 
            assay=assay, 
            K = 20, 
            genelst=genelst,
            T_latent=T_latent,
            min.frac = min.frac)
  }
  if(mode=="multiple"){
    if(factor==sample){
      stop("Please make sure the sample and factor be different in multiple sample mode!")
    }
    out <- HighOrderScore_multipleSample(
      obj=obj, 
      factor=factor, 
      reduction=NULL, 
      sample = sample, 
      assay=assay, 
      K = 20, 
      genelst=genelst,
      T_latent=T_latent,
      min.frac = min.frac)
  }
  return(out)
}

geneFrac <- function(x=NULL){
  cnt <-abs(x)
  cnt[which(cnt>0)]<-1
  frq<-rowSums(cnt)/ncol(cnt)
  names(frq) <- rownames(x)
  return(frq)
}

getKNN <- function(x=NULL, k=NULL){
  a<-sort(x)[k]
  pos<-which(x<=a)
  return(pos)
}


DEGs_plot<-function(dt=NULL,T1=NULL,T2=NULL,T3=NULL,size=4){
  library(ggrepel)
  # colum name should be X2_score, X2_pval, avg_logFC, p_adj
  if(sum(c("X2_score","X2_pval","avg_logFC","p_adj")%in%colnames(dt))<4) {
    stop("make sure the following 4 column inside: X2_score, X2_pval, avg_logFC, p_adj")
  }
  
  dt$label <-"NA"
  dt$label[abs(dt$avg_logFC)>T1 & dt$X2_pval>10^(0-T2)]<-"DE genes only"
  dt$label[abs(dt$avg_logFC)<T1 & dt$X2_pval<10^(0-T3)]<-"High-order genes only"
  dt$label[abs(dt$avg_logFC)>T1 & dt$X2_pval<10^(0-T3)]<-"Both"
  dt$gene <- rownames(dt)
  dt$gene[is.na(dt$label)]<-""
  
  dt$label <-factor(dt$label,levels=c("NA","DE genes only","High-order genes only", "Both"))
  dt$p_mismatch <-0-(log10(dt$X2_pval))
  max_lim <- max(abs(dt$avg_logFC))
  p <- ggplot(dt,aes(x=avg_logFC,y=p_mismatch,color=label)) + 
    geom_point() +  theme_bw() +  xlab("avg_logFC") + 
    ylab("-log10(p)[Cell population mismatch]") + 
    geom_text_repel(data=dt[dt$label%in%c("DE genes only","High-order genes only","Both"),], 
                    aes(avg_logFC,p_mismatch,label=gene),nudge_x = 0,size=size) + 
    geom_vline(xintercept = 0, linetype="dotted",color="black",size=0.5) + 
    scale_colour_manual(values=vega_20) + xlim(c(0-max_lim, max_lim))
  
  return(p)
}

marker_plot<-function(gene=NULL, obj1=NULL, obj2=NULL, merge=NULL,lst=c("X","Y")){
  
  #lst<-"CD69"
  #lst<-"TCR-V-2"
  # fold change 
  p1 <- VlnPlot(merge,features=gene)
  p2 <- FeaturePlot(obj1,features=gene)+ggtitle(lst[1])
  p2 <- customize_Seurat_FeaturePlot(p2)
  p3 <- FeaturePlot(obj2,features=gene)+ggtitle(lst[2])
  p3 <-customize_Seurat_FeaturePlot(p3)
  p <- ggarrange(p1,p2,p3,ncol=3)
  return(p)
}


Diff_Plot <- function(dt=NULL, nullmod=NULL){
  library(reshape2)
  dt <- out$diff
  dt <- melt(dt)
  dt$gene <- as.character(dt$Var1)
  dt$gene[abs(dt$value)<4.1]<-""
  dt$rank <- rep(seq(1,nrow(out$diff),1),20)
  #p <- 
  for(i in seq(1,ncol(dt),1)){
    a <- dt[,i]
  }
  colnames(dt)<-paste0("CCA",seq(1,20,1))
  
  
  
  
  nullmod <- out$null_mod
  dt <- out$diff[,3]
  nullmod <- nullmod[abs(nullmod)<max(abs(dt))]
  dt <- data.frame("X2"=c(dt, nullmod[seq(1,length(dt),1)]),
                   "method"=c(rep("Observed", length(dt)), rep("NullModel", length(dt))),
                   "gene"=c(names(dt),names(dt)))
  dt <- dt[order(dt$X2),]
  dt$label <- as.character(dt$gene)
  
  dt$rank <- seq(1,nrow(dt),1)
  T <- quantile(nullmod,probs=c(0.995))
  dt$label[abs(dt$X2)<=T[[1]] | dt$method=="NullModel"] <-"filter"
  write.csv(dt[dt$method=="Observed",],file="~/project/tmp/test.3.csv")
  
  
  dt$value[abs(dt$value)>20] <- 20
  p <- ggplot(dt,aes(x=rank,y=value)) + 
    geom_point() +  theme_bw() +  xlab("Rank") + ylab("Difference") + 
    facet_wrap(~Var2,ncol=5) + 
    geom_hline(yintercept=-4.1,linetype="dashed",color="red") + 
    geom_hline(yintercept=4.1,linetype="dashed",color="red")
  
  
  
  
  tiff("~/project/tmp/highOrder_gene.tiff", width=20, height =20, res =300, units = "in", compression = "lzw")
  print(p)
  dev.off()
  
    
  
    geom_text_repel(data=dt[!dt$label=="filter",], max.overlaps=Inf,
                    aes(rank,X2,label=label),nudge_x =0,size=1)
  
  
    geom_text_repel(data=dt[dt$label%in%c("DE genes only","High-order genes only","Both"),], 
                    aes(avg_logFC,p_mismatch,label=gene),nudge_x = 0,size=size) + 
    geom_vline(xintercept = 0, linetype="dotted",color="black",size=0.5) + 
    scale_colour_manual(values=vega_20) + xlim(c(0-max_lim, max_lim))

}

score_plot<-function(gene=NULL,merge=NULL,value=NULL){
  
  #lst<-"CD69"
  #lst<-"TCR-V-2"
  # fold change 
  value<-as.data.frame(value)
  merge@meta.data$gene_diff <- value$Diff_level
  merge@meta.data$MismatchScore <- value$Mismatch_score
  p1 <- VlnPlot(merge,features=gene,pt.size=0)
  p2 <- FeaturePlot(merge,features="gene_diff")+ggtitle("Difference level")
  p2 <- customize_Seurat_FeaturePlot(p2, expression.threshold = -100,gradient.use = c("red", "grey",  "blue"))
  p3 <- FeaturePlot(merge,features="MismatchScore")+ggtitle("Mismatch score")
  p3 <- customize_Seurat_FeaturePlot(p3, expression.threshold = -100,gradient.use = c( "red", "blue"))
  p <- ggarrange(p1,p2,ncol=2)
  return(p)
}

#### For FeaturePlot() ####
# Function to pseudo customize FeaturePlots
customize_Seurat_FeaturePlot <- function(p, alpha.use = 1, 
                                         gradient.use = c("gray", "red"), 
                                         expression.threshold = 0, 
                                         is.log1p.transformed = F) {
  
  #print(p$data)
  
  # Order data by gene expresion level
  p$data <- p$data[order(p$data[,4]),]
  
  # define 
  if(p$data[,4])
  
  # Define lower limit of gene expression level
  if (isTRUE(is.log1p.transformed)) {
    expression.threshold <- expression.threshold
  } else {
    expression.threshold <- log1p(expression.threshold)
  }
  
  # Compute maximum value in gene expression
  max.exp <- max(p$data[,4])
  
  # Fill points using the gene expression levels
  p$layers[[1]]$mapping$fill <- p$layers[[1]]$mapping$colour
  
  # Define transparency of points
  p$layers[[1]]$mapping$alpha <- alpha.use
  
  # Change fill and colour gradient values
  p <- p + scale_colour_gradientn(colours = gradient.use, guide = F, limits = c(expression.threshold, max.exp), na.value = "grey") +
    scale_fill_gradientn(colours = gradient.use, name = expression(atop(Expression, (log))), limits = c(expression.threshold, max.exp), na.value = "grey") +
    scale_alpha_continuous(range = alpha.use, guide = F)
  
}


VolcanoPlot <- function(res=NULL, features=NULL, pCutoff=NULL, FCcutoff=NULL, pointSize=2, labSize=3){

  p <- EnhancedVolcano(res,
                  lab = rownames(res),
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  selectLab = features,
                  xlab = bquote(~Log[2]~ 'fold change'),
                  pCutoff = pCutoff,
                  FCcutoff = FCcutoff,
                  pointSize = pointSize,
                  labSize = labSize,
                  labCol = 'black',
                  labFace = 'bold',
                  boxedLabels = TRUE,
                  title="",
                  subtitle="",
                  colAlpha = 4/5,
                  legendPosition = 'right',
                  legendLabSize = 1,
                  legendIconSize =1,
                  drawConnectors = TRUE,
                  widthConnectors = 1.0,
                  colConnectors = 'black')
  return(p)
}


KNN_smooth <- function(mat=NULL, reduction=NULL, K=10){
  
  dis <- cdist(reduction,reduction)
  knn_mat <- matrix(0, nrow(mat),ncol(mat))
  for(i in seq(1,ncol(mat),1)){
    pos <- getKNN(dis[i,],k=K)
    tp <- rowMeans(as.matrix(mat[,pos]))
    knn_mat[,i] <- tp
  }
  return(knn_mat)
}

