# Usage (archive root): Rscript scripts/02_foxm1_gse274558.R
# Exploratory FOXM1 summary in GSE274558 (NM vs NC, pre-leukaemic marrow). Preprocessing as in Phase 1/2.
suppressMessages({library(edgeR); library(limma)})
dir.create("results", showWarnings = FALSE)
o <- read.delim("inputs/Mouse_Human_Orthologs_Ensembl116.tsv", quote = ""); names(o) <- c("mid","msym","hid","hsym","type","conf")
sym <- unique(o[, c("mid","msym")])
raw <- read.delim(gzfile("inputs/GSE274558_raw_counts.txt.gz"), row.names = 1, check.names = FALSE)[, c("NC1","NC2","NC3","NM1","NM2","NM3")]
group <- factor(rep(c("NC","NM"), each = 3), levels = c("NC","NM")); design <- model.matrix(~group)
y <- calcNormFactors(DGEList(raw[filterByExpr(DGEList(raw, group = group), design), ], group = group))
lc <- cpm(y, log = TRUE, prior.count = 2); v <- voom(y, design)
z <- t(scale(t(lc))); z <- z[complete.cases(z), ]; stopifnot(identical(rownames(z), rownames(v)))
m <- sym[sym$mid %in% rownames(z) & nzchar(sym$msym), ]; m <- m[!duplicated(m$msym) & !duplicated(m$msym, fromLast = TRUE), ]
rownames(z)[match(m$mid, rownames(z))] <- m$msym
n <- read.csv("inputs/mapped_regulons.csv"); map <- read.csv("inputs/regulator_identifier_mapping.csv"); map <- map[nzchar(map$mouse_symbol), ]
i <- match(n$TF, map$original); n$TF <- ifelse(is.na(i), n$TF, map$mouse_symbol[i])
fx <- unique(n[n$TF == "Foxm1", c("target","sign")]); stopifnot(!any(duplicated(fx$target)))
prolif <- c("Mki67","Top2a","Pcna","Mcm2","Ccnb1","Cdk1","Bub1","Aurkb")
used <- fx[fx$target %in% rownames(z), ]
sig <- read.csv("signatures/fusion_522_mouse.csv")
used$in_proliferation_set <- used$target %in% prolif
used$in_fusion_522 <- used$target %in% sig$mouse_symbol
used$NM_minus_NC_mean_z <- rowMeans(z[used$target, 4:6]) - rowMeans(z[used$target, 1:3])
used$limma_logFC <- topTable(eBayes(lmFit(v, design)), coef = 2, number = Inf, sort.by = "none")[match(used$target, rownames(z)), "logFC"]
write.csv(data.frame(fx, tested = fx$target %in% rownames(z)), "signatures/foxm1_targets_mouse.csv", row.names = FALSE)
write.csv(used, "results/foxm1_GSE274558_targets.csv", row.names = FALSE)
score <- function(g, s = rep(1, length(g))) colMeans(z[g, , drop = FALSE] * s)
sc <- rbind(Foxm1_targets_signed = score(used$target, used$sign),
  Foxm1_targets_without_proliferation_set = score(used$target[!used$in_proliferation_set], used$sign[!used$in_proliferation_set]),
  proliferation_set = score(intersect(prolif, rownames(z))))
fid <- sym$mid[sym$msym == "Foxm1"]; stopifnot(length(fid) == 1, fid %in% rownames(lc))
sc <- rbind(Foxm1_mRNA_log2CPM = lc[fid, ], sc)
write.csv(data.frame(measure = rownames(sc), sc, check.names = FALSE), "results/foxm1_GSE274558_sample_values.csv", row.names = FALSE)
tt <- t(apply(sc, 1, function(x) c(mean_NC = mean(x[1:3]), mean_NM = mean(x[4:6]), NM_minus_NC = mean(x[4:6]) - mean(x[1:3]),
  welch_p = t.test(x[4:6], x[1:3])$p.value, all_NM_above_all_NC = min(x[4:6]) > max(x[1:3]))))
write.csv(data.frame(measure = rownames(tt), tt), "results/foxm1_GSE274558_summary.csv", row.names = FALSE)
# gene-set tests: self-contained (fry), competitive (CAMERA, estimated correlation), expression-matched random sets
sets <- list(targets = used, targets_no_prolif = used[!used$in_proliferation_set, ])
ave <- rowMeans(lc); bin <- cut(rank(ave, ties.method = "first"), 20, labels = FALSE); names(bin) <- rownames(z)
delta <- rowMeans(z[, 4:6]) - rowMeans(z[, 1:3])
set.seed(20261126)
res <- do.call(rbind, lapply(names(sets), function(nm) {
  a <- sets[[nm]]; j <- match(a$target, rownames(z)); w <- numeric(nrow(v)); w[j] <- a$sign
  fr <- fry(v, index = list(s = j), design = design, contrast = 2, gene.weights = w)
  # competitive test on signed values: flip repressed targets before CAMERA
  v2 <- v; neg <- j[a$sign < 0]; v2$E[neg, ] <- -v2$E[neg, ]
  ca <- camera(v2, index = list(s = j), design = design, contrast = 2, inter.gene.cor = NA)
  obs <- mean(delta[a$target] * a$sign); pool <- split(seq_along(delta), bin)
  null <- replicate(10000, { k <- unlist(lapply(split(seq_along(j), bin[j]), function(q) sample(pool[[as.character(bin[j[q[1]]])]], length(q))))
    mean(delta[k] * a$sign[unlist(split(seq_along(j), bin[j]))]) })
  data.frame(set = nm, genes = nrow(a), activating = sum(a$sign > 0), repressing = sum(a$sign < 0), fry_direction = fr$Direction, fry_p = fr$PValue,
    camera_direction = ca$Direction, camera_correlation = ca$Correlation, camera_p = ca$PValue,
    camera_p_preset_cor_0.01 = camera(v2, index = list(s = j), design = design, contrast = 2, inter.gene.cor = 0.01)$PValue,
    observed_mean_delta = obs, random_median = median(null), random_p_one_sided = (1 + sum(null >= obs)) / 10001) }))
write.csv(res, "results/foxm1_GSE274558_set_tests.csv", row.names = FALSE)
d <- read.csv("inputs/NM_vs_NC_all_genes_annotated.csv")
write.csv(d[d$gene == fid, ], "results/foxm1_GSE274558_mRNA_DE.csv", row.names = FALSE)
cat("targets in file:", nrow(fx), "tested:", nrow(used), "in 522:", sum(used$in_fusion_522), "in proliferation set:", paste(used$target[used$in_proliferation_set], collapse = ","), "\n")
print(tt, digits = 3); print(res, digits = 3); print(round(sc, 2))
