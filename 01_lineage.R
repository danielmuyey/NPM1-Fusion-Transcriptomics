# Usage (from archive root): Rscript scripts/01_lineage.R
suppressMessages({library(edgeR); library(limma)})
w <- function(x, f) write.csv(x, file.path("results", f), row.names = FALSE)
load1 <- function(f) { e <- new.env(); load(f, e); get(ls(e)[1], e) }
ig <- load1("ref/singler/data/immgen.rda"); mr <- load1("ref/singler/data/mouse.rnaseq.rda")

## ---- signatures -----------------------------------------------------------
grp_immgen <- function(types, main) {
  s <- sub(".*\\((.*)\\)$", "\\1", types); g <- rep(NA_character_, length(types))
  g[s %in% c("LTHSC","SC.LT34F","SC.ST34F","SC.STSL","SC.MPP34F","SC.CD150-CD48-")] <- "HSC_MPP"
  g[s %in% c("GMP","SC.CMP.DR","SC.MDP")] <- "Myeloid_progenitor"
  g[main == "Neutrophils"] <- "Granulocyte_neutrophil"
  g[main == "Monocytes"] <- "Monocyte"
  g[main %in% c("B cells","B cells, pro") & !grepl("CLP", s)] <- "B_cell"
  g[main == "T cells" & !grepl("^T\\.(DN|DP|SP|ISP|ETP|4SP|8SP|4\\+8int|4int8)", s)] <- "T_cell"
  g[main == "Mast cells"] <- "Mast_cell"; g[main == "Basophils"] <- "Basophil"
  o <- is.na(g) & (main %in% c("DC","Macrophages","ILC","NK cells","NKT","Tgd","Eosinophils") |
                   s %in% c("SC.MEP","SC.CDP","MLP") | grepl("CLP", s))
  g[o] <- paste0("other:", ifelse(main[o] %in% c("Stem cells","B cells"), s[o], main[o]))
  g
}
grp_rnaseq <- function(main) {
  m <- c(Erythrocytes="Erythroid", Granulocytes="Granulocyte_neutrophil", Monocytes="Monocyte", `B cells`="B_cell", `T cells`="T_cell",
         `NK cells`="other:NK", `Dendritic cells`="other:DC", Macrophages="other:Macrophages")
  unname(m[main])
}
derive <- function(ref, g, targets, source, drop = character()) {
  keep <- !is.na(g); x <- ref$data[, keep]; g <- g[keep]; ty <- ref$types[keep]
  key <- paste(g, ty, sep = "||")
  tm <- sapply(split(seq_along(key), key), function(i) rowMeans(x[, i, drop = FALSE]))
  tg <- sub("\\|\\|.*", "", colnames(tm))
  gm <- sapply(split(seq_along(tg), tg), function(i) rowMeans(tm[, i, drop = FALSE]))
  design <- do.call(rbind, lapply(unique(g), function(l) data.frame(source=source, group=l, n_samples=sum(g==l), cell_types=paste(sort(unique(ty[g==l])), collapse="; "))))
  sig <- do.call(rbind, lapply(targets, function(l) {
    comp <- setdiff(colnames(gm), c(l, drop))
    if (l == "T_cell") comp <- setdiff(comp, c("other:NKT","other:Tgd"))
    margin <- gm[, l] - apply(gm[, comp, drop = FALSE], 1, max)
    s <- sort(margin[margin >= 1], decreasing = TRUE); s <- head(s, 50)
    if (!length(s)) return(NULL)
    data.frame(source = source, lineage = l, symbol = names(s), log2_margin = round(unname(s), 3))
  }))
  list(sig = sig, design = design, means = gm)
}
lin <- c("HSC_MPP","Myeloid_progenitor","Granulocyte_neutrophil","Monocyte","B_cell","T_cell","Mast_cell","Basophil")
A <- derive(ig, grp_immgen(ig$types, ig$main_types), lin, "ImmGen")
B <- derive(mr, grp_rnaseq(mr$main_types), c("Erythroid","Granulocyte_neutrophil","Monocyte","B_cell","T_cell"), "MouseRNAseq")
# DEVIATION from the plan, added after the first run: the main rule gave only 1 myeloid-progenitor gene.
# This version lets progenitors compete with mature lineages only (not with HSC/MPP or other progenitors).
prog <- c("HSC_MPP","other:SC.MEP","other:SC.CDP","other:MLP","other:SC.CLP","other:proB.CLP")
P <- derive(ig, grp_immgen(ig$types, ig$main_types), "Myeloid_progenitor", "ImmGen vs mature only (post hoc)", drop = prog)
P$sig$lineage <- "Myeloid_progenitor_vs_mature"
w(rbind(A$design, B$design), "signatures/reference_samples_by_group.csv")
primary <- rbind(A$sig, B$sig[B$sig$lineage == "Erythroid", ])
secondary <- B$sig[B$sig$lineage != "Erythroid", ]
secondary <- rbind(secondary, P$sig)
e <- new.env(); load("inputs/mMCPcounter_signatures_GCRm38.RData", e); mm <- e$mMCPcounter_signatures_GCRm38
hemat <- c("T cells","CD8 T cells","NK cells","B derived","Memory B cells","Monocytes / macrophages","Monocytes","Granulocytes","Mast cells","Eosinophils","Neutrophils","Basophils")
mm <- unique(data.frame(source="mMCPcounter", lineage=mm$Denomination, symbol=mm$Gene.Symbol, ensembl=mm$ENSEMBL.ID))
mm <- mm[mm$lineage %in% hemat, ]

## ---- expression -----------------------------------------------------------
raw <- read.delim(gzfile("inputs/GSE274558_raw_counts.txt.gz"), row.names = 1, check.names = FALSE)
raw <- raw[, c("NC1","NC2","NC3","NM1","NM2","NM3")]
group <- factor(rep(c("NC","NM"), each = 3), levels = c("NC","NM"))
de <- read.csv("inputs/NM_vs_NC_all_genes_annotated.csv")
y <- DGEList(raw[de$gene, ], group = group, lib.size = colSums(raw))
design <- model.matrix(~group)
keep <- filterByExpr(DGEList(raw, group = group), design)
stopifnot(setequal(rownames(raw)[keep], de$gene))
y <- DGEList(raw[keep, ], group = group); y <- calcNormFactors(y)
L <- cpm(y, log = TRUE, prior.count = 1)
sym <- de$mouse_symbol[match(rownames(L), de$gene)]
ids <- function(s) rownames(L)[!is.na(sym) & sym %in% s]
mm$gene <- ifelse(mm$ensembl %in% rownames(L), mm$ensembl, NA)
primary$gene <- rownames(L)[match(primary$symbol, sym)]; secondary$gene <- rownames(L)[match(secondary$symbol, sym)]
allsig <- rbind(primary, secondary); allsig$tested <- !is.na(allsig$gene)
allsig$edgeR_logFC_NM_vs_NC <- de$edgeR_logFC[match(allsig$gene, de$gene)]
allsig$concordant_DE <- de$concordant[match(allsig$gene, de$gene)]
w(allsig, "signatures/derived_lineage_markers.csv")
mm$tested <- !is.na(mm$gene); w(mm, "signatures/mMCPcounter_haematopoietic_markers.csv")
sets <- lapply(split(primary$gene, primary$lineage), function(v) unique(v[!is.na(v)]))
sets2 <- lapply(split(secondary$gene, secondary$lineage), function(v) unique(v[!is.na(v)]))
names(sets2) <- paste0(names(sets2), ifelse(names(sets2) == "Myeloid_progenitor_vs_mature", " [post hoc]", " [RNAseq ref]"))
setsM <- lapply(split(mm$gene, mm$lineage), function(v) unique(v[!is.na(v)])); setsM <- setsM[lengths(setsM) > 0]
names(setsM) <- paste0(names(setsM), " [mMCP]")

## ---- scores ---------------------------------------------------------------
Z <- t(scale(t(L))); R <- apply(L, 2, function(v) rank(v) / length(v)); rownames(R) <- rownames(L)
score <- function(S, M, f = colMeans) t(sapply(S, function(g) f(M[g, , drop = FALSE])))
all_sets <- c(sets, sets2)
sc <- list(mean_z = score(all_sets, Z), rank_percentile = score(all_sets, R),
           mMCP_median_log2CPM = score(setsM, L, function(m) apply(m, 2, median)))
stat <- do.call(rbind, lapply(names(sc), function(m) {
  s <- sc[[m]]; S <- c(all_sets, setsM)[rownames(s)]
  r <- data.frame(method = m, signature = rownames(s), genes_tested = lengths(S),
    mean_NM = rowMeans(s[, group == "NM"]), mean_NC = rowMeans(s[, group == "NC"]))
  r$diff_NM_minus_NC <- r$mean_NM - r$mean_NC
  r$welch_p <- apply(s, 1, function(v) t.test(v[group == "NM"], v[group == "NC"])$p.value)
  r$welch_FDR <- p.adjust(r$welch_p, "BH")
  r$samples_separate <- apply(s, 1, function(v) min(v[group=="NM"]) > max(v[group=="NC"]) | max(v[group=="NM"]) < min(v[group=="NC"]))
  r$flag <- ifelse(r$genes_tested < 10, "fewer than 10 genes: unreliable", "")
  cbind(r, s)
}))
w(stat, "scores/lineage_scores_by_sample.csv")
w(data.frame(gene = rownames(L), symbol = sym, L), "scores/log2CPM_tested_genes.csv")

## ---- CAMERA and fry on lineage sets ---------------------------------------
v <- voom(y, design)
test_sets <- Filter(function(g) length(g) >= 5, c(all_sets, setsM))
idx <- ids2indices(test_sets, rownames(v))
cm <- camera(v, idx, design, contrast = 2, inter.gene.cor = NA, sort = FALSE)
cf <- camera(v, idx, design, contrast = 2, sort = FALSE)   # preset correlation 0.01
fr <- fry(v, idx, design, contrast = 2, sort = "none")
w(data.frame(signature = names(idx), n = cm$NGenes, direction = cm$Direction,
  camera_estimated_correlation = cm$Correlation, camera_p = cm$PValue, camera_FDR = p.adjust(cm$PValue, "BH"),
  camera_preset0.01_p = cf$PValue, camera_preset0.01_FDR = p.adjust(cf$PValue, "BH"),
  fry_direction = fr$Direction, fry_p = fr$PValue, fry_FDR = p.adjust(fr$PValue, "BH"), fry_mixed_p = fr$PValue.Mixed),
  "camera_roast/lineage_sets_camera_fry.csv")

## ---- are DE genes enriched for lineage markers? ----------------------------
enr <- function(universe, q, S, level) do.call(rbind, lapply(names(S), function(n) {
  s <- intersect(S[[n]], universe); k <- length(intersect(q, s))
  data.frame(level = level, signature = n, background = length(universe), query = length(q), markers_in_background = length(s),
    overlap = k, fold = if (length(s) && length(q)) (k/length(q))/(length(s)/length(universe)) else NA,
    p = phyper(k - 1, length(s), length(universe) - length(s), length(q), lower.tail = FALSE),
    overlap_genes = paste(sort(sym[match(intersect(q, s), rownames(L))]), collapse = ";")) }))
E <- list(); S <- c(all_sets, setsM)
S$`ANY primary or mMCP marker` <- unique(unlist(c(sets, setsM)))
hm <- de[!is.na(de$human_id), ]
Sh <- lapply(S, function(g) unique(hm$human_id[match(g, hm$gene, nomatch = 0)]))
for (d in c("both","up_in_NM","up_in_NC")) {
  pick <- function(t) t$concordant & (d == "both" | (d == "up_in_NM" & t$edgeR_logFC > 0) | (d == "up_in_NC" & t$edgeR_logFC < 0))
  a <- enr(de$gene, de$gene[pick(de)], S, "mouse genes"); a$direction <- d
  b <- enr(hm$human_id, hm$human_id[pick(hm)], Sh, "mapped human orthologs"); b$direction <- d
  b$overlap_genes <- sapply(b$signature, function(n) paste(sort(hm$mouse_symbol[hm$human_id %in% intersect(Sh[[n]], hm$human_id[pick(hm)])]), collapse=";"))
  E <- c(E, list(a, b))
}
E <- do.call(rbind, E); E$FDR <- ave(E$p, E$level, E$direction, FUN = function(p) p.adjust(p, "BH"))
w(E, "marker_enrichment/DE_lineage_marker_enrichment.csv")
excl <- unique(unlist(c(sets, setsM)))
w(data.frame(gene = excl, mouse_symbol = sym[match(excl, rownames(L))], human_id = de$human_id[match(excl, de$gene)]), "marker_enrichment/excluded_lineage_markers.csv")
writeLines(capture.output(sessionInfo()), "results/sessionInfo_R.txt")

cat("genes derived / tested\n"); print(table(allsig$source, allsig$lineage)); print(table(allsig$source, allsig$lineage, allsig$tested)[,, "TRUE"])
print(stat[, 1:10], digits = 3, row.names = FALSE)
print(read.csv("results/camera_roast/lineage_sets_camera_fry.csv"), digits = 3)
print(E[E$direction != "both" & E$level == "mouse genes", c("direction","signature","markers_in_background","overlap","fold","FDR")], digits = 3, row.names = FALSE)
cat("excluded markers:", length(excl), " concordant among them:", sum(de$gene[de$concordant] %in% excl), " human-mapped concordant:", sum(hm$gene[hm$concordant] %in% excl), "\n")
