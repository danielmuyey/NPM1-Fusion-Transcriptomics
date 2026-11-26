# Usage (archive root): Rscript scripts/04_tf_resolved.R
# Repeat the Phase 1 TF scoring after translating regulator names to mouse symbols.
suppressMessages({library(edgeR); library(limma)})
o <- read.delim("inputs/Mouse_Human_Orthologs_Ensembl116.tsv", quote = ""); names(o) <- c("mid","msym","hid","hsym","type","conf")
sym <- unique(o[, c("mid","msym")])
raw <- read.delim(gzfile("inputs/GSE274558_raw_counts.txt.gz"), row.names = 1, check.names = FALSE)[, c("NC1","NC2","NC3","NM1","NM2","NM3")]
group <- factor(rep(c("NC","NM"), each = 3), levels = c("NC","NM")); design <- model.matrix(~group)
y <- calcNormFactors(DGEList(raw[filterByExpr(DGEList(raw, group = group), design), ], group = group))
lc <- cpm(y, log = TRUE, prior.count = 2)   # as in Phase 1
z <- t(scale(t(lc))); z <- z[complete.cases(z), ]
m <- sym[sym$mid %in% rownames(z) & nzchar(sym$msym), ]; m <- m[!duplicated(m$msym) & !duplicated(m$msym, fromLast = TRUE), ]
rownames(z)[match(m$mid, rownames(z))] <- m$msym
v <- voom(y, design); stopifnot(identical(nrow(v), nrow(z)))
n <- read.csv("inputs/mapped_regulons.csv")
map <- read.csv("results/tf/regulator_identifier_mapping.csv"); map <- map[nzchar(map$mouse_symbol), ]
n$TF_original <- n$TF; i <- match(n$TF, map$original); n$TF <- ifelse(is.na(i), n$TF, map$mouse_symbol[i])
key <- paste(n$TF, n$target); conflict <- names(which(tapply(n$sign, key, function(s) length(unique(s)) > 1)))
cat("edges:", nrow(n), " duplicated after renaming:", sum(duplicated(key)), " sign conflicts removed:", length(conflict), "\n")
cat("target names not mouse symbols:", length(setdiff(n$target, o$msym)), "of", length(unique(n$target)), "\n")
n <- n[!key %in% conflict & !duplicated(key), ]
total <- table(n$TF); n <- n[n$target %in% rownames(z), ]; rr <- split(n, n$TF)
for (k in c(10, 20)) {
  r <- rr[sapply(rr, nrow) >= k]
  sc <- t(sapply(r, function(a) colMeans(z[a$target, , drop = FALSE] * a$sign)))
  f <- eBayes(lmFit(sc, design))
  fy <- do.call(rbind, lapply(r, function(a) {
    w <- numeric(nrow(v)); j <- match(a$target, rownames(z)); w[j] <- a$sign
    fry(v, index = list(s = j), design = design, contrast = 2, gene.weights = w)[, c("Direction","PValue")] }))
  res <- data.frame(TF = rownames(sc), targets_used = sapply(r, nrow), targets_in_file = as.integer(total[rownames(sc)]),
    renamed = rownames(sc) %in% map$mouse_symbol, NM_minus_NC = f$coefficients[, 2], p = f$p.value[, 2], FDR = p.adjust(f$p.value[, 2], "BH"),
    fry_direction = fy$Direction, fry_p = fy$PValue, fry_FDR = p.adjust(fy$PValue, "BH"))
  res <- res[order(res$p), ]; write.csv(res, sprintf("results/tf/TF_resolved_min%d.csv", k), row.names = FALSE)
  write.csv(data.frame(TF = rownames(sc), sc), sprintf("results/tf/TF_resolved_scores_min%d.csv", k), row.names = FALSE)
  cat("min", k, "tested", nrow(res), "limma FDR<.05:", sum(res$FDR < .05), " fry FDR<.05:", sum(res$fry_FDR < .05), " median delta", median(res$NM_minus_NC), " frac>0", mean(res$NM_minus_NC > 0), "\n")
  print(head(res, 25), row.names = FALSE, digits = 3)
}
a <- read.csv("results/tf/TF_resolved_min10.csv"); b <- read.csv("results/tf/TF_resolved_min20.csv")
cat("significant in both thresholds:", intersect(a$TF[a$FDR < .05], b$TF[b$FDR < .05]), "\n")
print(a[a$TF %in% c("Sp1","Stat1","Stat3","Runx1","Trp53","Foxm1","Gata1","Spi1","Cebpa","Hoxa9","Myc","Pax5","Meis1"), ], row.names = FALSE, digits = 3)
# checks: sign weights flip direction; regulators whose names did not change reproduce Phase 1
j <- 1:30; d1 <- fry(v, index = list(s = j), design = design, contrast = 2, gene.weights = rep(1, nrow(v)))$Direction
d2 <- fry(v, index = list(s = j), design = design, contrast = 2, gene.weights = rep(-1, nrow(v)))$Direction
stopifnot(d1 != d2)
for (k in c(10, 20)) {
  old <- read.csv(sprintf("inputs/TF_NM_vs_NC_min%d.csv", k)); new <- read.csv(sprintf("results/tf/TF_resolved_min%d.csv", k))
  x <- merge(old, new, by = "TF"); cat("min", k, "unchanged names:", nrow(x), " max |delta difference|:", max(abs(x$NM_minus_NC.x - x$NM_minus_NC.y)),
    " tested old/new:", nrow(old), nrow(new), "\n")
}
