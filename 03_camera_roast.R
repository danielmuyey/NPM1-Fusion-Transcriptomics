# Usage (archive root): Rscript scripts/03_camera_roast.R
# Correlation-aware pathway tests on all GO/Reactome sets mapped to mouse through one-to-one orthologs.
suppressMessages({library(edgeR); library(limma)})
raw <- read.delim(gzfile("inputs/GSE274558_raw_counts.txt.gz"), row.names = 1, check.names = FALSE)[, c("NC1","NC2","NC3","NM1","NM2","NM3")]
group <- factor(rep(c("NC","NM"), each = 3), levels = c("NC","NM")); design <- model.matrix(~group)
de <- read.csv("inputs/NM_vs_NC_all_genes_annotated.csv")
y <- calcNormFactors(DGEList(raw[filterByExpr(DGEList(raw, group = group), design), ], group = group))
v <- voom(y, design)
hid <- de$human_id[match(rownames(v), de$gene)]
excl <- read.csv("results/marker_enrichment/excluded_lineage_markers.csv")$gene
run <- function(drop) {
  x <- if (drop) v[!rownames(v) %in% excl, ] else v
  h <- de$human_id[match(rownames(x), de$gene)]
  G <- strsplit(readLines("inputs/msigdb_7.5.1_human_GO_Reactome_ensembl.gmt"), "\t")
  idx <- lapply(G, function(f) which(h %in% f[-(1:2)])); names(idx) <- sapply(G, `[`, 1)
  coll <- setNames(sapply(G, `[`, 2), names(idx))
  idx <- idx[lengths(idx) >= 10 & lengths(idx) <= 500]
  cm <- camera(x, idx, design, contrast = 2, inter.gene.cor = NA, sort = FALSE)
  c1 <- camera(x, idx, design, contrast = 2, sort = FALSE)
  fr <- fry(x, idx, design, contrast = 2, sort = "none")
  r <- data.frame(collection = coll[names(idx)], term = names(idx), n = cm$NGenes, camera_direction = cm$Direction,
    inter_gene_correlation = cm$Correlation, camera_p = cm$PValue, camera_preset_p = c1$PValue,
    fry_direction = fr$Direction, fry_p = fr$PValue, fry_mixed_p = fr$PValue.Mixed)
  for (m in c("camera","camera_preset","fry","fry_mixed")) r[[paste0(m,"_FDR")]] <- ave(r[[paste0(m,"_p")]], r$collection, FUN = function(p) p.adjust(p, "BH"))
  r
}
for (d in c(FALSE, TRUE)) {
  r <- run(d); tag <- if (d) "marker_excluded" else "all_genes"
  write.csv(r, sprintf("results/camera_roast/GO_Reactome_camera_fry_%s.csv", tag), row.names = FALSE)
  cat("==", tag, nrow(r), "\n")
  print(aggregate(cbind(camera = camera_FDR < .05, camera_preset = camera_preset_FDR < .05, fry = fry_FDR < .05, fry_mixed = fry_mixed_FDR < .05) ~ collection, r, sum))
  cat("min FDR camera/fry:", min(r$camera_FDR), min(r$fry_FDR), "\n")
  print(head(r[order(r$camera_p), c("term","n","camera_direction","inter_gene_correlation","camera_p","camera_FDR","fry_FDR")], 12), row.names = FALSE)
}
