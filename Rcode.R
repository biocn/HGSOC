setwd('/P199_OV_All/OV_GSE')
source('codes/mg_base.R')

cox_result=readMatrix('cox/GSE_OV_cox.results')
#head(cox_result)
#row.names(cox_result)[cox_result$p.value<=0.01]
cox_genes=row.names(cox_result)[cox_result$p.value<=0.01]
cox_result_enrich=mg_clusterProfiler(genes = cox_genes)

g1=dotplot_batch(cox_result_enrich$Enrich_tab,FDR = F,dbs = c('pathway_KEGG','geneontology_Biological_Process'))

summary(cox_result_enrich$KEGG)
summary(cox_result_enrich$GO_BP)

cox_result_enrich$Enrich_tab

gene_exp=readMatrix('data_work/OV_counts_normalizted_symbol_id_clinic.txt')
gene_exp=t(gene_exp)
head(gene_exp)
colnames()
gene_exp1=gene_exp[-c(1:5),]
head(gene_exp1)
cluster=readMatrix('cox/cluster.txt')

head(cluster)
paste0('Cluster',cluster$Cluster[match(colnames(gene_exp1),row.names(cluster))])
gene_exp1=crbind2DataFrame(gene_exp1)
head(gene_exp)
sum(as.numeric(gene_exp[1,])<50)
table(gene_exp[4,])

dim(gene_exp1)
gene_exp1=gene_exp1[-18739,]

gene_exp1[18738,]

risk=ifelse(mod1$Score>median(mod1$Score),'High risk','Low risk')
risk=risk[order(risk)]
head(gene_exp1)
gsea_risk=mg_RunGSEA_exp(exp_data = gene_exp1[,match(names(risk),colnames(gene_exp1))]
                    ,sample_group =risk
                    ,outFolder = 'GSEA')

gsea_kegg_risk=parseGSEAResult('GSEA/my_analysis.Gsea.1591586310175')
tail(gsea_kegg_risk$EnrichTable)
gsea_kegg_risk$EnrichTable[gsea_kegg_risk$EnrichTable$FDR<0.05,1]

gsea1=plot_GSEA_By_node(gsea_kegg_risk,TermName = gsea_kegg_risk$EnrichTable[gsea_kegg_risk$EnrichTable$FDR<0.05,1][1])
gsea2=plot_GSEA_By_node(gsea_kegg_risk,TermName = gsea_kegg_risk$EnrichTable[gsea_kegg_risk$EnrichTable$FDR<0.05,1][2])

gsea.gal=ggpubr::ggarrange(gsea1,gsea2, ncol = 2, nrow = 1,labels = c('A','B'))

savePDF(filename = 'Fig8.pdf',width = 8,height = 4,plot = gsea.gal)

gsea=mg_RunGSEA_exp(exp_data = gene_exp1[,match(row.names(cluster),colnames(gene_exp1))]
                    ,sample_group =paste0('Cluster',cluster$Cluster)
                    ,outFolder = 'GSEA')

gsea_kegg=parseGSEAResult('GSEA/my_analysis.Gsea.1591260147979')


gsea_kegg$EnrichTable[gsea_kegg$EnrichTable$NP<0.01,]

g2=plot_GSEA_By_nodes(gsea_kegg,TermNames = gsea_kegg$EnrichTable[gsea_kegg$EnrichTable$NP<0.001&gsea_kegg$EnrichTable$FDR<0.2,]$Term)
gal=ggpubr::ggarrange(g1,g2, ncol = 1, nrow = 2,labels = c('','C'))
gal
savePDF('cox_gene_enrich.pdf',gal,12,10)

gsea=mg_RunGSEA_exp(exp_data = gene_exp1[,match(row.names(cluster),colnames(gene_exp1))]
                    ,sample_group =paste0('Cluster',cluster$Cluster),gmt_Path='HALLMARK'
                    ,outFolder = 'GSEA')
gsea_HALLMARK=parseGSEAResult('GSEA/my_analysis.Gsea.1591260536762')

gsea_HALLMARK$EnrichTable

sig.genes=c('ATHL1','PMFBP1','LSAMP','NLRP7','PLTP')
sig.genes=c('AKR1B10','ANGPT4')

os=as.numeric(gene_exp[3,])
ev=as.numeric(gene_exp[4,])
sig.exp=crbind2DataFrame(gene_exp[match(sig.genes,row.names(gene_exp)),])
head(sig.exp)
i=1
plts=list()
cx_al=rbind()
for(i in 1:2){
x=as.numeric(sig.exp[i,])
cx_al=rbind(cx_al,coxRun(data.frame(os,ev,x)))
g=ggplotKMCox(data.frame(os,ev,ifelse(x>median(x),'High','Low')),title = row.names(sig.exp)[i])
plts=c(plts,list(g))
}

i=2
x=as.numeric(sig.exp[i,])
#roc1=mg_plot_cox(os,ev,as.numeric(sig.exp[1,]))
#roc2=mg_plot_cox(os,ev,as.numeric(sig.exp[2,]))

match(names(mod1$Score),colnames(gene_exp))

cli.mod=createCoxModel(crbind2DataFrame(cbind(t(gene_exp1[1:2,]),risk=mod1$Score)),time = os,event =ev )
cox_batch(dat = t(crbind2DataFrame(cbind(t(gene_exp[1:2,]),risk=mod1$Score))),time = os,event = ev)
#kmg=ggpubr::ggarrange(plotlist = plts, ncol = 2, nrow = 1,labels =LETTERS[1:2])
summary(cli.mod$Cox)
roc1=ggplotTimeROC(os,ev,as.numeric(sig.exp[1,]))
roc2=ggplotTimeROC(os,ev,as.numeric(sig.exp[2,]))
dim(sig.exp)
kmg=ggpubr::ggarrange(roc1,plts[[1]],roc2,plts[[2]], ncol = 2, nrow = 2,labels =LETTERS[3:6])

rg=mg_ridges_plot(cbind(paste0('Cluster',cluster$Cluster)[match(colnames(sig.exp),row.names(cluster))],t(sig.exp)))
mg_ridges_plot()

sig.exp.clust=gene_exp1[match(sig.genes,row.names(gene_exp1)),match(row.names(cluster),colnames(gene_exp1))]

mg_muti_cor_plot(t(sig.exp.clust))

mg_quick_cor_plot(dat1 = t(sig.exp))
rg1=mg_ridges_plot(cbind(Repeats=glaso$repeats),s_height = 100,xlab = 'Number of repeats')
rg.al=ggpubr::ggarrange(rg1,rg, ncol = 2, nrow = 1,labels =c('A','B'))

fig4=ggpubr::ggarrange(rg.al,kmg, ncol = 1, nrow = 2,labels =c('','E'),heights = c(0.3,1))

savePDF('Fig4.pdf',fig4,width = 8,height = 10)

g_st=groupViolin(crbind2DataFrame(t(sig.exp.clust)),group = paste0('Cluster',cluster$Cluster),melt = F
            ,ylab = 'Gene expression')
savePDF('Fig4B.pdf',g_st,width = 6,height = 6)
stg=paste0('Stage ',gene_exp[2,match(row.names(cluster),colnames(gene_exp1))])

groupViolin(crbind2DataFrame(t(sig.exp.clust))[which(stg!='Stage 1'&stg!='Stage 2'),]
            ,group = stg[which(stg!='Stage 1'&stg!='Stage 2')],melt = F
            ,ylab = 'Gene expression')



head(gene_exp)

table(gene_exp[2,match(row.names(cluster),colnames(gene_exp1))]
      ,paste0('Cluster',cluster$Cluster))
table(gene_exp[4,match(row.names(cluster),colnames(gene_exp1))]
      ,paste0('Cluster',cluster$Cluster))


mod1=createCoxModel(t(sig.exp),time = os,event =ev )
mod1.gal=plotCoxModel_Batch(mod1$Score,dat = t(sig.exp),time = os,event = ev
                   ,cutoff = median(mod1$Score))
savePDF('Fig5.pdf',mod1.gal,width = 10,height = 8)

tcga_exp=readMatrix('TCGA/TCGA_OV_counts_normalized_clinic.txt')
tcga_samples=readMatrix('TCGA/samples.txt')
head(tcga_exp)
tcga_exp_risk=readMatrix('TCGA/riskScore.txt',header = F)
tcga_exp_risk=tcga_exp_risk[match(row.names(tcga_exp),row.names(tcga_exp_risk)),]
tcga_exp[,1:2]
dim(tcga_exp)
library(pROC)
library(ggplot2)
data3=data.frame(tcga_exp[,2],tcga_exp[,1],tcga_exp_risk)
colnames(data3)=c('OS.day','VitalStatus','Risk')
fit_validation <- survivalROC::survivalROC(Stime = data3$OS.day, status = data3$VitalStatus, marker = data3$Risk, predict.time = 100, method = "KM")
optimalCutoff_training <- median(fit_validation$cut.values)
median(data3$Risk)
median(tcga_exp_risk)
dim(data3)
getwd()
writeMatrix(gsub('\\.','-',row.names(tcga_exp)),outpath = 'Supplementary Table1.txt',row=F)

data3=data3[data3[,1]>30,]

mg_plot_cox(data3[,1],data3[,2],data3[,3])
risk_training=ifelse(data3$Risk<median(data3$Risk),'Low','High')
dat=data.frame(data3[,1],data3[,2],risk_training)
colnames(dat)=c('OS.day','VitalStatus','Risk')
risk_training <- survfit( Surv(OS.day, VitalStatus) ~ Risk,data = dat )
par(mar = c(4,4,2,1))
gt12=ggpubr::ggarrange(gt1$plot,gt1$table, ncol = 1, nrow = 2)
gt2=ggplotKMCox(dat)

ggpubr::ggarrange(gt12,gt2, ncol = 2, nrow = 1,labels =LETTERS[1:2])
dat[which(data3[,1]>365*3),2]=0
dat[which(data3[,1]>365*3),1]=365*3
#dat=dat[dat[,1]>30,]
#sf<-survival::survfit(Surv(OS.day,VitalStatus) ~ Risk,data=dat)
risk_training <- survfit( Surv(OS.day, VitalStatus) ~ Risk,data = dat )
ggsurvplot(risk_training,data = dat,pval = T,conf.int =T,conf.int.style ="step"
           ,risk.table="abs_pct",palette=c('red','blue')
           ,risk.table.y.text = FALSE,ncensor.plot = TRUE,fontsize=3)  

survminer::ggsurvplot(sf, data = dat
                      #, palette = palette, #jco palette 
                          , pval = TRUE
                          ,surv.median.line='hv'
                           #,conf.int = T
                           ,conf.int.style ='step'
                           #, pval.coord=c(0, 0.2), #Add p-value 
                           ,risk.table = "abs_pct"
                      ,risk.table.y.text = FALSE
                           #legend.title = title
                           #,legend.labs = labs
                           #,conf.int=show_confint
)

head(dat)
head(data3)
coxRun(data3)

plot(risk_training, mark.time = TRUE,col=c('red','blue')
     #,xlab=paste("Survival time in day","\np=",round(p,5))
     ,ylab = "Survival probabilities"
     #,main=title
     )
legend('topright',paste0(gsub('groups=','',names(sf$strata)),'(N=',sf$n,')'), col = colKm,
       lty = c(1,1, 1, 1),lwd=c(1,1,1,1),merge = TRUE)

plotKMCox(dat)

tcga_exp[,2]
inds=which(tcga_exp[,2]>30&tcga_exp[,2]<365*3)
p.os=tcga_exp[inds,2]
p.ev=tcga_exp[inds,1]
p.ev[which(p.os>365*3)]=0
p.os[which(p.os>365*3)]=365*3
plotCoxModel_Batch(tcga_exp_risk[inds],dat = tcga_exp[inds,match(sig.genes,colnames(tcga_exp))]
                   ,time = p.os
                   ,event = p.ev
                   ,cutoff = median(tcga_exp_risk[inds]))

plotRiskScoreModel(tcga_exp_risk[inds],dat = tcga_exp[inds,match(sig.genes,colnames(tcga_exp))]
                   ,time = p.os
                   ,event = p.ev
                   ,cutoff=median(tcga_exp_risk[inds]))
ggplotKMCox(data.frame(p.os,p.ev
                       ,ifelse(tcga_exp_risk[inds]>median(tcga_exp_risk[inds]),'H','L')))

tcga.cmp=intersect(gsub('\\.','-',row.names(tcga_exp)),row.names(tcga_samples))
tcga_exp=tcga_exp[match(tcga.cmp,gsub('\\.','-',row.names(tcga_exp))),]
tcga_exp=tcga_exp[which(tcga_exp[,2]>30,),]
sig.tcga.exp=tcga_exp[,c(1:2,match(sig.genes,colnames(tcga_exp)))]

cmp.smaples=tcga_samples[match(tcga.cmp,row.names(tcga_samples)),]
sum(cmp.smaples$age_at_initial_pathologic_diagnosis<50)
table(cmp.smaples$A6_Stage)
table(tcga_exp[,1])
tcga_samples

#head(sig.tcga.exp)
all.tcga=cox_batch(dat = t(tcga_exp[,3:ncol(tcga_exp)]),time = tcga_exp[,2],event = tcga_exp[,1])
#all.tcga$
dim(tcga_exp)
sig.cox=cox_result[cox_result$p.value<0.01,]
sig.vd1=all.tcga[which(all.tcga$p.value<0.05),]
sig.vd2=all.vd1.cox[which(all.vd1.cox$p.value<0.05),]
sig.vd3=all.vd2.cox[which(all.vd2.cox$p.value<0.05),]
sig.vd4=all.vd3.cox[which(all.vd3.cox$p.value<0.05),]

glaso=readMatrix('lasoo/lasso_gene_index.txt')
glaso[match(lst.cmp,glaso$Gene),]
sum(glaso$repeats>30)

plot(density(glaso$repeats))

lst.cmp=intersect(intersect(row.names(sig.cox),row.names(sig.vd3)),row.names(sig.vd4))

lst.cmp=intersect(row.names(sig.cox),row.names(sig.vd1))
#mg_venn_plot(list(A=row.names(sig.cox),B=row.names(sig.vd1),C=row.names(sig.vd2),D=row.names(sig.vd3),E=row.names(sig.vd4)))

cbind(glaso[match(lst.cmp,glaso$Gene),2],sig.cox[match(lst.cmp,row.names(sig.cox)),c(5,2)],
sig.vd1[match(lst.cmp,row.names(sig.vd1)),c(1,2)],
all.vd1.cox[match(lst.cmp,row.names(all.vd1.cox)),c(1,2)],
all.vd2.cox[match(lst.cmp,row.names(all.vd2.cox)),c(1,2)],
all.vd3.cox[match(lst.cmp,row.names(all.vd3.cox)),c(1,2)])[glaso[match(lst.cmp,glaso$Gene),2]>0,]

tcga.risk.score=as.matrix(sig.tcga.exp[,3:7])%*%c(-0.05470921,2.61980185,1.05829238,0.42534627,1.01587661)
dim(sig.tcga.exp)
cx_mod=createCoxModel(sig.tcga.exp[,3:4],time = sig.tcga.exp[,2],event = sig.tcga.exp[,1])
tcga.vd.roc=plotCoxModel_Batch(riskScore = cx_mod$Score,time = sig.tcga.exp[,2],event = sig.tcga.exp[,1]
                   ,cutoff = median(cx_mod$Score),dat = sig.tcga.exp[,3:4])

#run_c
sig.tcga.exp[which(sig.tcga.exp[,2]>365*3),1]=0
sig.tcga.exp[which(sig.tcga.exp[,2]>365*3),2]=365*3
sig.tcga.exp=sig.tcga.exp[which(sig.tcga.exp[,2]>30),]
mg_plot_cox(os = sig.tcga.exp[,2]
            ,event = sig.tcga.exp[,1]
            ,rickscore = as.numeric(cx_mod$Score))
ggplotKMCox(data.frame(sig.tcga.exp[,2],sig.tcga.exp[,1],ifelse(cx_mod$Score>median(cx_mod$Score),'H','L')))


sig.genes1=c('PGGHG',sig.genes[-1])

GSE53963=getGEOExpData('GSE53963')
GSE53963.exp=exp_probe2symbol_v2(GSE53963$Exp$GPL6480_41000_Data_col1,GPL = 'GPL6480')

head(GSE53963.exp)
GSE53963.exp.sig=GSE53963.exp[match(sig.genes1,row.names(GSE53963.exp)),]
GSE53963.os=as.numeric(as.character(GSE53963$Sample$`channel2:time_fu_months`))
GSE53963.ev=ifelse(as.character(GSE53963$Sample$`channel2:vital_status`)=='Dead',1,0)

all.vd1.cox=cox_batch(dat = GSE53963.exp,time = GSE53963.os,event = GSE53963.ev)

vd1_mod=createCoxModel(t(GSE53963.exp.sig),time = GSE53963.os,event = GSE53963.ev)
mg_plot_cox(os = GSE53963.os*30,event = GSE53963.ev,rickscore = vd1_mod$Score)
ggplotKMCox(data.frame(GSE53963.os,GSE53963.ev,ifelse(vd1_mod$Score>median(vd1_mod$Score),'H','L')))

plotCoxModel_Batch(vd1_mod$Score,dat = t(GSE53963.exp.sig),time = GSE53963.os
                   ,event = GSE53963.ev,cutoff = median(vd1_mod$Score))




GSE49997=getGEOExpData('GSE49997')
GSE49997.os=as.numeric(as.character(GSE49997$Sample$`os month`))
GSE49997.ev=as.numeric(as.character(GSE49997$Sample$`os event`))
GSE49997.exp=exp_probe2symbol_v2(GSE49997$Exp$GPL2986_32878_Data_col1,GPL='GPL2986')
GSE49997.exp.sig=GSE49997.exp[match(sig.genes1,row.names(GSE49997.exp)),]
all.vd2.cox=cox_batch(dat = GSE49997.exp,time = GSE49997.os,event = GSE49997.ev)

vd2_mod=createCoxModel(t(GSE49997.exp.sig),time = GSE49997.os,event = GSE49997.ev)
mg_plot_cox(os = GSE49997.os,event = GSE49997.ev,rickscore = vd2_mod$Score)
ggplotKMCox(data.frame(GSE49997.os,GSE49997.ev,ifelse(vd2_mod$Score>median(vd2_mod$Score),'H','L')))
dim(GSE49997.exp.sig)
sum(is.na( GSE49997.os))
#GSE51088=getGEOExpData('GSE51088')
#GSE51088.ev=as.numeric(gsub('Alive','1',gsub('Dead','1',as.character(GSE51088$Sample$`channel2:patient status`))))
#GSE51088.os=as.numeric(as.character(GSE51088$Sample$`channel2:follow up months`))
#GSE51088.exp=exp_probe2symbol_v2(GSE51088$Exp$GPL7264_20173_Data_col1,GPL='GPL7264')
#dim(GSE51088.exp)
#GSE51088.exp.sig=GSE51088.exp[match(sig.genes1,row.names(GSE51088.exp)),]
#vd2_mod=createCoxModel(t(GSE49997.exp.sig),time = GSE49997.os,event = GSE49997.ev)
#mg_plot_cox(os = GSE49997.os,event = GSE49997.ev,rickscore = vd2_mod$Score)
#ggplotKMCox(data.frame(GSE49997.os,GSE49997.ev,ifelse(vd2_mod$Score>median(vd2_mod$Score),'H','L')))


GSE26712=getGEOExpData('GSE26712')
GSE26712.exp=exp_probe2symbol_v2(GSE26712$Exp$GPL96_22283_Data_col1,GPL='GPL96')
GSE26712.exp.sig=GSE26712.exp[match(intersect(sig.genes1,row.names(GSE26712.exp)),row.names(GSE26712.exp)),]
GSE26712.ev=as.numeric(gsub('AWD \\(alive with disease\\)','0',gsub('DOD \\(dead of disease\\)','1',as.character(GSE26712$Sample$status))))
GSE26712.os=as.numeric(as.character(GSE26712$Sample$`survival years`))
vd3_mod=createCoxModel(t(GSE26712.exp.sig),time = GSE26712.os,event = GSE26712.ev)
mg_plot_cox(os = GSE26712.os,event = GSE26712.ev,rickscore = vd3_mod$Score)
ggplotKMCox(data.frame(GSE26712.os,GSE26712.ev,ifelse(vd3_mod$Score>median(vd3_mod$Score),'H','L')))

all.vd3.cox=cox_batch(dat = GSE26712.exp,time = GSE26712.os,event = GSE26712.ev)
GSE26712$Anno$GPL96

dim(GSE26712.exp)
vd4.sig.exp=t(GSE26712.exp[match(intersect(sig.genes,row.names(GSE26712.exp)),row.names(GSE26712.exp)),])
vd4.inds=which(!is.na(GSE26712.os))
vd4.cx_mod=createCoxModel(vd4.sig.exp[vd4.inds,],time = GSE26712.os[vd4.inds],event = GSE26712.ev[vd4.inds])
dim(vd4.sig.exp)
dim(vd4.sig.exp[vd4.inds,])
table(GSE26712.ev[vd4.inds])
head(GSE26712$Sample)

gse.roc=plotCoxModel_Batch(riskScore = vd4.cx_mod$Score,time = GSE26712.os[vd4.inds],event = GSE26712.ev[vd4.inds]
                   ,cutoff = median(vd4.cx_mod$Score),dat = vd4.sig.exp[vd4.inds,])
tcga.vd.roc

fig4=ggpubr::ggarrange(rg.al,kmg, ncol = 1, nrow = 2,labels =c('','E'),heights = c(0.3,1))
savePDF('Fig6.pdf',tcga.vd.roc,width = 10,height = 8)
savePDF('Fig7.pdf',gse.roc,width = 10,height = 8)

GSE135886=getGEOExpData('GSE135886')
GSE135886.exp=exp_probe2symbol_v2(GSE135886$Exp$GPL15314_25848_Data_col1
                                  ,anno = cbind(as.character(GPL15314$match_uniq$Probe),as.character(GPL15314$match_uniq$vals2)))
head(GSE135886.exp)
match(sig.genes,row.names(GSE135886.exp))
#GSE135886$Anno

#GPL15314=re_annotation_by_probe(GSE135886$Anno$GPL15314[,c(1,8)])
sum(row.names(GSE135886)=='ASHG19A3A017244')

as.character(GPL15314$match_uniq$Probe)[grep(sig.genes[2],as.character(GPL15314$match_uniq$vals2))]

GSE135886$Sample$`disease state`

dim(GSE135886$Sample)
ov.dat=mg_data_hub_OV()

ov.dat$ICGC_OV_US$CLINI
sig.genes=c('PMFBP1','LSAMP','NLRP7','PLTP')
icgc.xc=cox_batch(ov.dat$ICGC_OV_US$EXP,ov.dat$ICGC_OV_US$OS[,1],ov.dat$ICGC_OV_US$OS[,2])
icgc.xc[icgc.xc[,1]>0.0001&icgc.xc[,1]<0.01,]
icgc.xc[order(icgc.xc[,1]),]
sig.genes=c('AKR1B10','ANGPT4','COL8A1')
icgc.sig.dat=ov.dat$ICGC_OV_US$EXP[match(sig.genes,row.names(ov.dat$ICGC_OV_US$EXP)),]
icgc.md=createCoxModel(t(icgc.sig.dat),time = ov.dat$ICGC_OV_US$OS[,1],event = ov.dat$ICGC_OV_US$OS[,2])
icgc.md$Cox
icgc.md$Score
figS1=mg_plot_cox(ov.dat$ICGC_OV_US$OS[,1], ov.dat$ICGC_OV_US$OS[,2],icgc.md$Score)
figS1
dim(icgc.sig.dat)
savePDF('FigS1.pdf',figS1,width = 10,height = 4)


ov.dat$GSE26712$OS
ov.dat$GSE26712$CLINI
ov.dat$GSE17260$INFO
geo.sig.dat=ov.dat$GSE17260$EXP[match(sig.genes[-3],row.names(ov.dat$GSE17260$EXP)),]
geo.md=createCoxModel(t(geo.sig.dat),time = ov.dat$GSE17260$OS[,1],event = ov.dat$GSE17260$OS[,2])
geo.md$Cox
geo.md$Score
figS2=mg_plot_cox(ov.dat$GSE17260$OS[,1], ov.dat$GSE17260$OS[,2],geo.md$Score)

figS12=mg_merge_plot(figS1,figS2,ncol = 1,nrow=2)
savePDF('FigS1.pdf',figS12,width = 10,height = 8)

PMID32096871=c('PRKG1','SDF2L1','PPP1R12A')
PMID33193589=c('TGFBI','SFRP1','COL16A1','THY1','PPIB','BGN')
PMID30569721=c('IGF2','PEG3','DCN','LYPD1','RARRES1')

cal_cindex=function(data,genes,os,event){
   geo.sig.dat=data[match(genes,row.names(data)),]
   geo.md=createCoxModel(t(geo.sig.dat),time = os,event = event)
   cdx=mg_cal_cox_auc_cindex(geo.md$Score,event,os)
   return(cdx)
}
crbind2DataFrame(t(tcga_exp))->tcga.t

all.cindx=rbind(cal_cindex(ov.dat$GSE17260$EXP,sig.genes[-3],ov.dat$GSE17260$OS[,1],ov.dat$GSE17260$OS[,2])
,cal_cindex(ov.dat$GSE26712$EXP,sig.genes[-3],ov.dat$GSE26712$OS[,1],ov.dat$GSE26712$OS[,2])
,cal_cindex(ov.dat$GSE102073$EXP,sig.genes[-3],ov.dat$GSE102073$OS[,1],ov.dat$GSE102073$OS[,2])
,cal_cindex(ov.dat$ICGC_OV_US$EXP,sig.genes,ov.dat$ICGC_OV_US$OS[,1],ov.dat$ICGC_OV_US$OS[,2])
,cal_cindex(tcga.t,sig.genes[-3],as.numeric(tcga_exp[,2]),tcga_exp[,1]))
row.names(all.cindx)=c('GSE17260','GSE26712','GSE102073','ICGC','TCGA')

all_get_cindex=function(genes){
all.cindx=rbind(cal_cindex(ov.dat$GSE17260$EXP,genes,ov.dat$GSE17260$OS[,1],ov.dat$GSE17260$OS[,2])
                ,cal_cindex(ov.dat$GSE26712$EXP,genes,ov.dat$GSE26712$OS[,1],ov.dat$GSE26712$OS[,2])
                ,cal_cindex(ov.dat$GSE102073$EXP,genes,ov.dat$GSE102073$OS[,1],ov.dat$GSE102073$OS[,2])
                ,cal_cindex(ov.dat$ICGC_OV_US$EXP,genes,ov.dat$ICGC_OV_US$OS[,1],ov.dat$ICGC_OV_US$OS[,2])
                ,cal_cindex(tcga.t,genes,as.numeric(tcga_exp[,2]),tcga_exp[,1]))
row.names(all.cindx)=c('GSE17260','GSE26712','GSE102073','ICGC','TCGA')
return(all.cindx)
}

all_get_cindex(PMID32096871)->PMID32096871.cx
all_get_cindex(PMID33193589)->PMID33193589.cx
PMID30569721.cx=all_get_cindex(PMID30569721)

cbind(all.cindx[,6],PMID32096871.cx[,6],PMID33193589.cx[,6],PMID30569721.cx[,6])

figS3=mg_barplot(cbind(row.names(all.cindx),RiskScore=all.cindx[,6]
                 ,ThreeGeneModel=PMID32096871.cx[,6],SixGeneModel=PMID33193589.cx[,6]
                 ,FiveGeneModel=PMID30569721.cx[,6]),bar_width = 0.8,show_txt = F
           ,ylab = 'C-index')


figS12=mg_merge_plot(figS1,figS2,figS3,ncol = 1,nrow=3,labels = c('','','C'))
savePDF('FigS1.pdf',figS12,width = 10,height = 12)

