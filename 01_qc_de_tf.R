# Usage (from archive root): Rscript scripts/01_qc_de_tf.R
suppressPackageStartupMessages({library(edgeR);library(limma);library(DESeq2)})
w=function(d,f) write.csv(d,file.path("results",f),row.names=FALSE)
x=as.matrix(read.delim(gzfile("inputs/GSE274558_raw_counts.txt.gz"),row.names=1,check.names=FALSE))
stopifnot(all(x>=0),all(x==round(x)),!anyNA(x),!anyDuplicated(rownames(x)))
six=c("NM1","NM2","NM3","NC1","NC2","NC3"); x=x[,six]
g=factor(sub("[0-9]$","",six),levels=c("NC","NM")); des=model.matrix(~g)
keep=filterByExpr(DGEList(x),des); y=x[keep,]
d=calcNormFactors(DGEList(y))
w(data.frame(sample=six,group=g,library_size=colSums(x),detected_genes=colSums(x>0),
  library_size_filtered=colSums(y),TMM_factor=d$samples$norm.factors),"qc/sample_qc.csv")
lc=cpm(d,log=TRUE,prior.count=2)
w(data.frame(gene=rownames(lc),lc),"qc/log2cpm.csv")
p=prcomp(t(lc)); v=p$sdev^2/sum(p$sdev^2)
w(data.frame(sample=six,PC1=p$x[,1],PC2=p$x[,2],PC1_var=v[1],PC2_var=v[2]),"qc/pca.csv")
cr=cor(lc,method="spearman"); w(data.frame(sample=six,cr),"qc/spearman.csv")
h=hclust(as.dist(1-cr),method="average"); w(data.frame(h$merge,height=h$height),"qc/hclust_merge.csv")

fit_both=function(y){
  d=calcNormFactors(DGEList(y)); d=estimateDisp(d,des)
  e=topTags(glmQLFTest(glmQLFit(d,des),coef=2),Inf,sort.by="none")$table
  dds=DESeq(DESeqDataSetFromMatrix(y,data.frame(g=g,row.names=colnames(y)),~g),quiet=TRUE)
  r=results(dds,contrast=c("g","NM","NC"),alpha=0.05)
  data.frame(gene=rownames(y),edgeR_logFC=e$logFC,edgeR_logCPM=e$logCPM,edgeR_F=e$F,edgeR_p=e$PValue,edgeR_FDR=e$FDR,
    DESeq2_baseMean=r$baseMean,DESeq2_log2FC=r$log2FoldChange,DESeq2_stat=r$stat,DESeq2_p=r$pvalue,DESeq2_FDR=r$padj)
}
flag=function(t){
  s1=!is.na(t$edgeR_FDR)&t$edgeR_FDR<.05; s2=!is.na(t$DESeq2_FDR)&t$DESeq2_FDR<.05
  t$concordant=s1&s2&sign(t$edgeR_logFC)==sign(t$DESeq2_log2FC)
  t$concordant_lfc1=t$concordant&abs(t$edgeR_logFC)>=1&abs(t$DESeq2_log2FC)>=1
  t
}
summ=function(t,label){
  one=function(m,fc,fdr) data.frame(filter=label,method=m,genes_tested=nrow(t),
    sig=sum(fdr<.05,na.rm=TRUE),up_NM=sum(fdr<.05&fc>0,na.rm=TRUE),up_NC=sum(fdr<.05&fc<0,na.rm=TRUE),
    sig_lfc1=sum(fdr<.05&abs(fc)>=1,na.rm=TRUE),up_NM_lfc1=sum(fdr<.05&fc>=1,na.rm=TRUE),up_NC_lfc1=sum(fdr<.05&fc<=-1,na.rm=TRUE))
  both=function(m,col) {k=t[[col]]; data.frame(filter=label,method=m,genes_tested=nrow(t),sig=sum(k),
    up_NM=sum(k&t$edgeR_logFC>0),up_NC=sum(k&t$edgeR_logFC<0),sig_lfc1=sum(t$concordant_lfc1),
    up_NM_lfc1=sum(t$concordant_lfc1&t$edgeR_logFC>0),up_NC_lfc1=sum(t$concordant_lfc1&t$edgeR_logFC<0))}
  rbind(one("edgeR",t$edgeR_logFC,t$edgeR_FDR),one("DESeq2",t$DESeq2_log2FC,t$DESeq2_FDR),both("concordant","concordant"))
}
main=flag(fit_both(y)); w(main,"de/NM_vs_NC_all_genes.csv")
w(main[main$concordant,],"de/concordant_genes.csv"); w(main[main$concordant_lfc1,],"de/concordant_lfc1_genes.csv")
S=summ(main,"filterByExpr (primary)"); sens=list()
for(k in c(0.5,1,2)){
  t=flag(fit_both(x[rowSums(cpm(x)>=k)>=3,])); lab=paste0("CPM>=",k," in >=3")
  w(t,paste0("de/sensitivity_CPM",k,".csv")); S=rbind(S,summ(t,lab))
  a=main$gene[main$concordant]; b=t$gene[t$concordant]
  sens[[lab]]=data.frame(filter=lab,primary_concordant=length(a),this=length(b),shared=length(intersect(a,b)),
    jaccard=length(intersect(a,b))/length(union(a,b)))
}
w(S,"de/summary_counts.csv"); w(do.call(rbind,sens),"de/filter_sensitivity_overlap.csv"); print(S,row.names=FALSE)

# Orthologs: unique in both directions across all pairs, then one2one and confidence 1
o=read.delim("inputs/Mouse_Human_Orthologs_Ensembl116.tsv",quote="")
names(o)=c("mid","msym","hid","hsym","type","conf")
pr=unique(o[nzchar(o$hid),]); pr=pr[!duplicated(paste(pr$mid,pr$hid)),]
u=pr[!(pr$mid%in%pr$mid[duplicated(pr$mid)])&!(pr$hid%in%pr$hid[duplicated(pr$hid)]),]
u=u[u$type=="ortholog_one2one"&u$conf==1,]
sym=unique(o[,c("mid","msym")]); main$mouse_symbol=sym$msym[match(main$gene,sym$mid)]
j=match(main$gene,u$mid); main$human_id=u$hid[j]; main$human_symbol=u$hsym[j]
w(main,"de/NM_vs_NC_all_genes_annotated.csv")
w(data.frame(step=c("one2one conf1 pairs","genes passing primary filter","of these mapped (background)",
  "concordant","concordant mapped","concordant lfc1","concordant lfc1 mapped"),
  n=c(nrow(u),nrow(main),sum(!is.na(main$human_id)),sum(main$concordant),sum(main$concordant&!is.na(main$human_id)),
  sum(main$concordant_lfc1),sum(main$concordant_lfc1&!is.na(main$human_id)))),"mapping/mapping_audit.csv")

# TF activity: same scoring rule as the original run.R
z=t(scale(t(lc))); z=z[complete.cases(z),]
m=sym[sym$mid%in%rownames(z)&nzchar(sym$msym),]; m=m[!duplicated(m$msym)&!duplicated(m$msym,fromLast=TRUE),]
rownames(z)[match(m$mid,rownames(z))]=m$msym
n=read.csv("inputs/mapped_regulons.csv"); total=table(n$TF); n=n[n$target%in%rownames(z),]; rr=split(n,n$TF)
for(k in c(10,20)){
  r=rr[sapply(rr,nrow)>=k]
  sc=t(sapply(r,function(a) colMeans(z[a$target,,drop=FALSE]*a$sign)))
  f=eBayes(lmFit(sc,des))
  res=data.frame(TF=rownames(sc),targets_used=sapply(r,nrow),targets_in_file=as.integer(total[rownames(sc)]),
    NM_minus_NC=f$coefficients[,2],p=f$p.value[,2],FDR=p.adjust(f$p.value[,2],"BH"),
    name_is_mouse_symbol=rownames(sc)%in%o$msym)
  res=res[order(res$p),]; w(res,paste0("tf/TF_NM_vs_NC_min",k,".csv")); w(data.frame(TF=rownames(sc),sc),paste0("tf/scores_min",k,".csv"))
  cat("TF min",k,": tested",nrow(res),"FDR<0.05",sum(res$FDR<.05),"of which mouse symbols",sum(res$FDR<.05&res$name_is_mouse_symbol),"\n")
}
writeLines(capture.output(sessionInfo()),"results/sessionInfo_R.txt")
