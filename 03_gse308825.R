# Usage (archive root): Rscript scripts/03_gse308825.R
# Drug-EXPOSURE-associated expression in two paediatric AML xenograft models. Not patient response.
suppressMessages({library(edgeR); library(limma)})
dir.create("results/GSE308825", recursive = TRUE, showWarnings = FALSE)
cnt <- read.csv("inputs/external/GSE308825_counts.csv.gz", row.names = 1, check.names = FALSE, fileEncoding = "UTF-8-BOM")
s <- read.csv("metadata/GSE308825_samples.csv")
s$model <- sub("_.*", "", s$title); s$treatment <- sub("^[^_]+_([^_]+)_.*", "\\1", s$title); s$replicate <- sub(".*_", "", s$title)
stopifnot(setequal(s$library, colnames(cnt)))
cnt <- as.matrix(cnt[, s$library]); stopifnot(all(cnt == round(cnt)), all(cnt >= 0))
s$total_counts <- colSums(cnt); s$genes_detected <- colSums(cnt > 0)
write.csv(s, "results/GSE308825/sample_qc.csv", row.names = FALSE)
allcpm <- cpm(calcNormFactors(DGEList(cnt)), log = TRUE, prior.count = 2); colnames(allcpm) <- s$title
write.csv(data.frame(gene = rownames(cnt), allcpm, check.names = FALSE), "results/GSE308825/logCPM_all_samples.csv", row.names = FALSE)
pc <- do.call(rbind, lapply(unique(s$model), function(m) {
  x <- allcpm[, s$model == m]; keep <- rowMeans(x) > 1; v <- apply(x[keep, ], 1, var); x <- x[names(sort(v, decreasing = TRUE))[1:1000], ]
  p <- prcomp(t(x)); data.frame(model = m, sample = colnames(x), PC1 = p$x[, 1], PC2 = p$x[, 2],
    PC1_pct = 100 * p$sdev[1]^2 / sum(p$sdev^2), PC2_pct = 100 * p$sdev[2]^2 / sum(p$sdev^2)) }))
write.csv(pc, "results/GSE308825/pca.csv", row.names = FALSE)

fusion <- list(S398 = read.csv("signatures/fusion_398_human.csv"), S336 = read.csv("signatures/marker_excluded_human.csv"))
fx <- read.csv("signatures/foxm1_targets_human.csv"); fx <- fx[!is.na(fx$human_symbol) & nzchar(fx$human_symbol), ]
prolif <- c("MKI67","TOP2A","PCNA","MCM2","CCNB1","CDK1","BUB1","AURKB")

fit_contrasts <- function(model, treatments, block = FALSE, label) {
  k <- s$model == model & s$treatment %in% c("DMSO", treatments); si <- s[k, ]
  grp <- factor(si$treatment, levels = c("DMSO", treatments))
  design <- if (block) model.matrix(~ si$replicate + grp) else model.matrix(~ grp)
  colnames(design) <- sub("^grp", "", colnames(design))
  y <- DGEList(cnt[, si$library], group = grp); y <- y[filterByExpr(y, design), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  v <- voom(y, design); f <- eBayes(lmFit(v, design))
  lapply(setNames(treatments, treatments), function(t) {
    tt <- topTable(f, coef = t, number = Inf, sort.by = "none")
    list(id = paste(model, t, sep = "_"), analysis = label, model = model, treatment = t, n_drug = sum(grp == t), n_dmso = sum(grp == "DMSO"),
         tt = tt, v = v, design = design, coef = which(colnames(design) == t)) }) }

primary <- c(fit_contrasts("AML903", "Quiz", label = "primary"), fit_contrasts("AML001", c("Quiz","Palb","Combo"), label = "primary"))
sens <- fit_contrasts("AML001", c("Quiz","Palb","Combo"), block = TRUE, label = "replicate_blocked")
run <- c(primary, sens)

summ <- do.call(rbind, lapply(run, function(r) { t <- r$tt; data.frame(analysis = r$analysis, contrast = r$id, n_drug = r$n_drug, n_dmso = r$n_dmso,
  genes_tested = nrow(t), up_FDR05_lfc1 = sum(t$adj.P.Val < .05 & t$logFC >= 1), down_FDR05_lfc1 = sum(t$adj.P.Val < .05 & t$logFC <= -1),
  up_FDR05 = sum(t$adj.P.Val < .05 & t$logFC > 0), down_FDR05 = sum(t$adj.P.Val < .05 & t$logFC < 0)) }))
write.csv(summ, "results/GSE308825/DE_summary.csv", row.names = FALSE)
for (r in run) write.csv(data.frame(gene = rownames(r$tt), r$tt), gzfile(sprintf("results/GSE308825/DE_%s_%s.csv.gz", r$analysis, r$id)), row.names = FALSE)

settest <- function(r, genes, w = NULL) {
  j <- match(genes, rownames(r$v)); ok <- !is.na(j); j <- j[ok]; w <- if (is.null(w)) rep(1, length(j)) else w[ok]
  v <- r$v; v$E[j[w < 0], ] <- -v$E[j[w < 0], ]    # signed sets: flip repressed targets
  ca <- camera(v, index = list(s = j), design = r$design, contrast = r$coef, inter.gene.cor = NA)
  fr <- fry(v, index = list(s = j), design = r$design, contrast = r$coef)
  data.frame(genes_in_set = length(genes), genes_tested = length(j), mean_signed_logFC = mean(r$tt[j, "logFC"] * w),
    camera_direction = ca$Direction, camera_correlation = ca$Correlation, camera_p = ca$PValue, fry_direction = fr$Direction, fry_p = fr$PValue) }

cam <- do.call(rbind, lapply(run, function(r) do.call(rbind, lapply(names(fusion), function(sn) do.call(rbind, lapply(c("higher_in_NM","higher_in_NC"), function(d)
  cbind(analysis = r$analysis, contrast = r$id, signature = sn, signature_direction = d, settest(r, fusion[[sn]]$human_symbol[fusion[[sn]]$direction == d]))))))))
fox <- do.call(rbind, lapply(run, function(r) rbind(
  cbind(analysis = r$analysis, contrast = r$id, set = "FOXM1_targets_signed", settest(r, fx$human_symbol, fx$sign)),
  cbind(analysis = r$analysis, contrast = r$id, set = "FOXM1_targets_without_proliferation_set", settest(r, fx$human_symbol[!fx$human_symbol %in% prolif], fx$sign[!fx$human_symbol %in% prolif])),
  cbind(analysis = r$analysis, contrast = r$id, set = "proliferation_set", settest(r, prolif)))))
for (a in c("primary", "replicate_blocked")) {     # separate FDR families per analysis
  i <- cam$analysis == a; cam$camera_FDR[i] <- p.adjust(cam$camera_p[i], "BH"); cam$fry_FDR[i] <- p.adjust(cam$fry_p[i], "BH")
  i <- fox$analysis == a; fox$camera_FDR[i] <- p.adjust(fox$camera_p[i], "BH"); fox$fry_FDR[i] <- p.adjust(fox$fry_p[i], "BH") }
write.csv(cam, "results/GSE308825/signature_camera_fry.csv", row.names = FALSE)
write.csv(fox, "results/GSE308825/foxm1_set_tests.csv", row.names = FALSE)

mrna <- do.call(rbind, lapply(run, function(r) { t <- r$tt["FOXM1", ]; data.frame(analysis = r$analysis, contrast = r$id,
  tested = "FOXM1" %in% rownames(r$tt), logFC = t$logFC, p = t$P.Value, FDR_genomewide = t$adj.P.Val) }))
write.csv(mrna, "results/GSE308825/foxm1_mRNA_DE.csv", row.names = FALSE)
# per-sample values, z-scored within model (includes the single AML903 palbociclib sample: descriptive only)
val <- do.call(rbind, lapply(unique(s$model), function(m) {
  x <- allcpm[, s$model == m]; x <- x[rowMeans(cnt[, s$library[s$model == m]]) >= 10, ]; z <- t(scale(t(x))); z <- z[complete.cases(z), ]
  g <- fx$human_symbol %in% rownames(z); p <- intersect(prolif, rownames(z))
  data.frame(sample = colnames(x), model = m, treatment = s$treatment[s$model == m], FOXM1_logCPM = allcpm["FOXM1", s$model == m],
    FOXM1_target_score = colMeans(z[fx$human_symbol[g], ] * fx$sign[g]), targets_used = sum(g), proliferation_score = colMeans(z[p, ])) }))
write.csv(val, "results/GSE308825/foxm1_sample_values.csv", row.names = FALSE)
cov <- do.call(rbind, lapply(names(fusion), function(sn) data.frame(signature = sn, signature_genes = nrow(fusion[[sn]]),
  present_in_matrix = sum(fusion[[sn]]$human_symbol %in% rownames(cnt)),
  absent = paste(setdiff(fusion[[sn]]$human_symbol, rownames(cnt)), collapse = ";"))))
write.csv(cov, "results/GSE308825/signature_coverage.csv", row.names = FALSE)
print(s[, c("title","library","total_counts","genes_detected")]); print(summ); print(cam[, c(1:4,6:8,10:11,13:14)], digits = 2)
print(mrna, digits = 3); print(fox[, c(1:3,5,6,7,9,10,12,13)], digits = 3); print(cov[, 1:3]); print(val, digits = 3)
