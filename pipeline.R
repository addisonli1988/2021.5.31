###Part1 TCGA RNAseq数据下载###
library(GDCRNATools)

setwd("F:\\DATA")
project <- 'TCGA-UCEC'
rnadir <- paste(project, 'RNAseq', sep='/')
mirdir <- paste(project, 'miRNAs', sep='/')
####### Download RNAseq data #######
gdcRNADownload(project.id     = 'TCGA-UCEC', 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = rnadir)
####### Parse RNAseq metadata #######
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-UCEC',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)
####### Filter duplicated samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)
####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)
####### Filter Solid Tissue Normal samples in RNAseq metadata and marker m6A mutant #######
metaMatrix.RNA <-  subset(metaMatrix.RNA,sample_type =="PrimaryTumor")
write.csv(x = metaMatrix.RNA,file = "metaMatrix_RNA.csv" ,row.names = TRUE)
metaMatrix.RNA <- read.table("F:\\DATA\\UCEC_analysis_2020.9.15\\res.csv",
                             header = T,sep=",")
####### Merge RNAseq data #######
rnaCounts <- gdcRNAMerge(metadata  = metaMatrix.RNA, 
                         path      = rnadir, # the folder in which the data stored
                         organized = FALSE, # if the data are in separate folders
                         data.type = 'RNAseq')
####### Merge clinical data #######
clinicalDa <- gdcClinicalMerge(path = clinicaldir, key.info = TRUE)
####### Normalization of RNAseq data #######
rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)

###Part2 差异基因筛选及Heatmap绘制###
####### Differential gene expression analysis
DEall_mut <- gdcDEAnalysis(counts   = rnaCounts, 
                           group      = metaMatrix.RNA$Mutant, 
                           comparison = 'Mut-WT', 
                           method     = 'DESeq2')

### All DEGs
DEall_mut_sig <- gdcDEReport(deg = DEall_mut, gene.type = 'all')
gdcBarPlot(deg = DEall_mut_sig, angle = 45, data.type = 'RNAseq')
write.csv(x = DEall_mut_sig,file = "DEall_mut_sig.csv" ,row.names = TRUE)

##### DEG Heatmap #####
degALLName = rownames(DEall_mut_sig)
deg.id = degALLName
metadata = metaMatrix.RNA
rna.expr = rnaExpr
degDa <- rna.expr[deg.id,]

sampleCol = data.frame(SampleType = factor(metadata$Mutant))
rownames(sampleCol) = metadata$sample
ann_colors = list(SampleType = c(WT = "#436eee", Mut = "#EE0000"))

degDa_scale = apply(degDa,1,scale)
row.names(degDa_scale) <- colnames(degDa)
write.csv(x = degDa_scale,file = "degDa_scale.csv" ,row.names = TRUE)
Mut <- read.table("Mut.csv",header = T, sep=",")
Mut <- t(Mut)
WT <- read.table("WT.csv",header = T, sep=",")
WT <- t(WT)
library(ComplexHeatmap)
library(circlize)
col_fun <- colorRamp2(c(-4,0,4),c("navy", "white", "firebrick3"))
p1 <- Heatmap(Mut,name = "m6A altered", column_title = "m6A altered",
              col = col_fun,clustering_distance_columns = "euclidean",
              clustering_method_columns = "ward.D2",width = unit(8,"cm"),show_row_names = FALSE)

p2 <- Heatmap(WT,name = "m6A not altered",column_title = "m6A not altered",
              col = col_fun,clustering_distance_columns = "euclidean",
              clustering_method_columns = "ward.D2",width = unit(8,"cm"),show_row_names = FALSE)

p1+p2

####Part3 生存分析####
####单因素cox分析####
survOutput_RNA_sorted <- gdcSurvivalAnalysis(gene  = rownames(DEall_mut_sig), 
                                  method   = 'coxph', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA)
write.csv(x = survOutput_RNA_sorted,file = "survOutput_RNA_mut.csv" ,row.names = TRUE)
####Lasso回归分析####
library("survival")
library("survminer")
options(stringsAsFactors = F)
exprSet = na.omit(rnaExpr)
meta = metaMatrix.RNA

head(meta)
colnames(meta)
meta[,12][is.na(meta[,12])]=0
meta[,13][is.na(meta[,13])]=0
meta$days=as.numeric(meta[,12])+as.numeric(meta[,13])

meta=meta[,c(4,14,9,8,17)] 
colnames(meta)=c("ID","event","age","gender","days")
meta$event=ifelse(meta$event=="Alive",0,1)
meta$age=as.numeric(meta$age)
meta$days <- ifelse(meta$days <=0,1,meta$days)
meta$time=meta$days/30

phe = meta
phe$ID = toupper(phe$ID)
phe=phe[match(colnames(exprSet),phe$ID),]
colnames(exprSet)==phe$ID

surV_sorted_pvalue <-  subset(survOutput_RNA_sorted,pValue<0.05)
cox_genes = row.names(surV_sorted_pvalue)
cox_expr = exprSet[cox_genes,]
cox_expr = t(cox_expr)
colnames(cox_expr)= cox_genes
dat = cbind(phe,cox_expr)
dat$gender=factor(dat$gender)

library(glmnet)
x <- as.matrix(cox_expr)
y <- phe[,c("time","event")]

names(y) <- c("time","status")
y$time <- as.double(y$time)
y$status <- as.double(y$status)
y <- as.matrix(survival::Surv(y$time,y$status))
fit <- glmnet(x,y,family = "cox")
plot(fit,label=TRUE)
cvfit = cv.glmnet(x,y, family="cox")
plot(cvfit)
coef.min =coef(cvfit,s=cvfit$lambda.min)
Active.Index <- which(as.numeric(coef.min)!=0)
Active.coef <- as.numeric(coef.min)[Active.Index]
sig_gene_coef <- rownames(coef.min)[Active.Index]

####多因素cox回归分析建模####
multi_var <- paste0(sig_gene_coef,collapse = "+")
fml <- as.formula(paste0("Surv(time,event)~",multi_var))
multi_var_cox <- coxph(fml,data = dat)
ph_hypo_multi <- cox.zph(multi_var_cox)
ph_hypo_table <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]

correlation <- cor(dat[,rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05]],method = "pearson")
library("GGally")
ggpairs(dat[,rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05]],
        axisLabels = "show")+
  theme_bw()+
  theme(panel.background = element_rect(colour = "black",size = 1,fill = "white"),
        panel.grid = element_blank())
library("rms")
vif<- rms::vif(multi_var_cox)
sqrt(vif) < 2

ggforest(model = multi_var_cox,data = dat,main = "Hazard ratio of candidate genes",fontsize = 1)

####risk score计算####
riskscore <- function(survival_cancer_df,candidate_genes_for_cox, cox_report) {
      library('dplyr')
      risk_score_table <- survival_cancer_df[,candidate_genes_for_cox]
      for(each_sig_gene in colnames(risk_score_table)){
            risk_score_table$each_sig_gene <- risk_score_table[,each_sig_gene]*(summary(cox_report)$coefficients[each_sig_gene,1])
      }
  risk_score_table <- cbind(risk_score_table,'total_risk_score'=exp(rowSums(risk_score_table))) %>%
    cbind(survival_cancer_df[,c('ID','days','event')])
  risk_score_table <- risk_score_table[,c('ID','days','event', candidate_genes_for_cox, 'total_risk_score')]
  risk_score_table
}
candidate_genes_for_cox <- c(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05])
risk_score_table_multi_cox <- riskscore(dat_train,candidate_genes_for_cox,multi_var_cox)

####ROC曲线绘制option 1####
ROC_marginal <-timeROC(T=risk_score_table_multi_cox$days,
                  delta=risk_score_table_multi_cox$event,marker=risk_score_table_multi_cox$total_risk_score,
                  cause=1,weighting="cox",
                  times=quantile(pbc$time,probs=seq(0.2,0.8,0.1)))
ROC_marginal

plot(ROC_marginal,time=365,lwd=2)        
plot(ROC_marginal,time=1095,add=TRUE,col="blue",lwd=2) 
plot(ROC_marginal,time=1825,add=TRUE,col="grey50",lwd=2)
legend("bottomright",c("Y-1","Y-3","Y-5"),col=c("red","blue","grey50"),lty=1,lwd=2)

####ROC曲线绘制option 2####
multi_ROC <- function(time_vector, risk_score_table){
  library('survivalROC')
  single_ROC <- function(single_time){
  for_ROC <- survivalROC(Stime = risk_score_table$days,
                         status = risk_score_table$event,
                         marker = risk_score_table$total_risk_score,
                         predict.time = single_time, method = 'KM')
  data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP, 
             'Cut_values'=for_ROC$cut.values, 'Time_point'=rep(single_time, length(for_ROC$TP)),
             'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}
for_multi_ROC <- multi_ROC(time_vector = c(365*seq(3,5,0.2)), risk_score_table = risk_score_table_multi_cox)

#maybe AUCs are identical in different time points. So select the last time point indicating longer survival.
AUC_max_time <- for_multi_ROC$Time_point[which(for_multi_ROC$AUC == AUC_max)]
AUC_max_time <- AUC_max_time[!duplicated(AUC_max_time)]
AUC_max_time <- AUC_max_time[length(AUC_max_time)]
for_multi_ROC$Time_point <- as.factor(for_multi_ROC$Time_point)

#visualization of the ROC curves of multiple time points.
pROC <- ggplot(for_multi_ROC, aes(x = False_positive, y = True_positive, label = Cut_values, color = Time_point)) + 
  geom_roc(labels = F, stat = 'identity', n.cuts = 0) + 
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 2)+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'), 
        panel.grid = element_blank())+
  annotate("text",x = 0.75, y = 0.15,
           label = paste("AUC max = ", round(AUC_max, 2), '\n', 'AUC max time = ', AUC_max_time, ' days', sep = ''))
pROC

AUC_max <- max(for_multi_ROC$AUC)

#find the optimal cutoff value within the ROC curve of the optimal time point.
optimal_time_ROC_df <- for_multi_ROC[which(for_multi_ROC$Time_point == AUC_max_time),]
cut.off <- optimal_time_ROC_df$Cut_values[which.max(optimal_time_ROC_df$True_positive-optimal_time_ROC_df$False_positive)]
high_low <- (risk_score_table_multi_cox$total_risk_score > cut.off)
high_low[high_low == TRUE] <- 'high'
high_low[high_low == FALSE] <- 'low'
risk_score_table_multi_cox <- cbind(risk_score_table_multi_cox, high_low)

#KM_plot generation.
library('survminer')
#first edit the status of patients with OS > AUC max time. (censoring status=0 (Alive), OS=365*5 days)
risk_score_table_multi_cox$event[which(risk_score_table_multi_cox$days > AUC_max_time)] <- 0
risk_score_table_multi_cox$days[which(risk_score_table_multi_cox$days > AUC_max_time)] <- AUC_max_time
fit_km <- survfit(Surv(days, event) ~high_low, data = risk_score_table_multi_cox)     
ggsurvplot(fit_km, conf.int = F,pval = T,legend.title="total risk score",
           legend.labs=c(paste0('>',as.character(round(cut.off,2))), paste0('<=',as.character(round(cut.off,2)))), risk.table = T, 
           palette = c('red','blue'), surv.median.line = 'hv')


