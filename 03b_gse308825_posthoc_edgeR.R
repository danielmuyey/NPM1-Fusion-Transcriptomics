# Usage (archive root): Rscript scripts/03b_gse308825_posthoc_edgeR.R
# POST HOC robustness check, added after the planned limma-voom analysis returned no DE genes. Does not replace it.
suppressMessages(library(edgeR))
cnt <- read.csv("inputs/external/GSE308825_counts.csv.gz", row.names = 1, check.names = FALSE, fileEncoding = "UTF-8-BOM")
s <- read.csv("results/GSE308825/sample_qc.csv")
out <- do.call(rbind, lapply(list(list("AML903","Quiz"), list("AML001", c("Quiz","Palb","Combo"))), function(a) {
  si <- s[s$model == a[[1]] & s$treatment %in% c("DMSO", a[[2]]), ]; grp <- factor(si$treatment, levels = c("DMSO", a[[2]]))
  d <- model.matrix(~grp); y <- DGEList(cnt[, si$library], group = grp); y <- y[filterByExpr(y, d), , keep.lib.sizes = FALSE]
  y <- estimateDisp(calcNormFactors(y), d, robust = TRUE); f <- glmQLFit(y, d, robust = TRUE)
  do.call(rbind, lapply(seq_along(a[[2]]), function(i) { t <- topTags(glmQLFTest(f, coef = i + 1), Inf)$table
    data.frame(contrast = paste(a[[1]], a[[2]][i], sep = "_"), genes_tested = nrow(t), up_FDR05_lfc1 = sum(t$FDR < .05 & t$logFC >= 1),
      down_FDR05_lfc1 = sum(t$FDR < .05 & t$logFC <= -1), FDR05_any = sum(t$FDR < .05), min_FDR = min(t$FDR),
      FOXM1_logFC = t["FOXM1", "logFC"], FOXM1_FDR = t["FOXM1", "FDR"],
      genes_FDR05 = paste(head(rownames(t)[t$FDR < .05], 40), collapse = ";")) })) }))
write.csv(out, "results/GSE308825/posthoc_edgeR_QL_summary.csv", row.names = FALSE); print(out)
