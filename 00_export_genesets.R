# Usage: Rscript 00_export_genesets.R MSIGDBR_7.5.1_LIBRARY_DIR OUT.gmt
a=commandArgs(TRUE); e=new.env(); lazyLoad(file.path(a[1],"R/sysdata"),e)
s=as.data.frame(e$msigdbr_genesets); s=s[s$gs_subcat%in%c("GO:BP","GO:MF","GO:CC","CP:REACTOME"),]
m=merge(e$msigdbr_geneset_genes,unique(as.data.frame(e$msigdbr_genes)[,c("gene_id","human_ensembl_gene")]),by="gene_id")
m=m[m$gs_id%in%s$gs_id & !is.na(m$human_ensembl_gene),]; g=split(m$human_ensembl_gene,m$gs_id)
lab=c("GO:BP"="GOBP","GO:MF"="GOMF","GO:CC"="GOCC","CP:REACTOME"="REACTOME")
writeLines(sapply(seq_len(nrow(s)),function(i) paste(c(s$gs_name[i],lab[s$gs_subcat[i]],sort(unique(g[[s$gs_id[i]]]))),collapse="\t")),a[2])
