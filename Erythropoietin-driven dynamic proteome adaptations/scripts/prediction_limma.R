library(org.Mm.eg.db)
library(limma)
library(xlsx)
library(ggplot2)

###########
# LOAD DATA
setwd(file.path("../prediction"))
pred_int <- read.delim("predicted_protein_intensity.txt")

#######
# LIMMA
mymat <- pred_int
design <- model.matrix(~ 0+factor(c(0,rep(1, 5), rep(2, 3), rep(3, 3))))
colnames(design) <- c("t0", "early", "switch", "late")
fit <- lmFit(mymat, design)

# late Vs. early
contrast.matrix <- makeContrasts(late-early, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

res <- topTable(fit2, adjust="BH", number=nrow(mymat), sort.by="P")

up <- rownames(res)[res$logFC > 0 & res$adj.P.Val <= 0.05]
down <- rownames(res)[res$logFC < 0 & res$adj.P.Val <= 0.05]

# early Vs. t0
contrast.matrix <- makeContrasts(early-t0, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

res <- topTable(fit2, adjust="BH", number=nrow(mymat), sort.by="P")

up <- rownames(res)[res$logFC > 0 & res$adj.P.Val <= 0.05]
down <- rownames(res)[res$logFC < 0 & res$adj.P.Val <= 0.05]

# late Vs. t0
contrast.matrix <- makeContrasts(late-t0, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

res <- topTable(fit2, adjust="BH", number=nrow(mymat), sort.by="P")

up <- rownames(res)[res$logFC > 0 & res$adj.P.Val <= 0.05]
down <- rownames(res)[res$logFC < 0 & res$adj.P.Val <= 0.05]

# early Vs. switch
contrast.matrix <- makeContrasts(early-switch, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

res <- topTable(fit2, adjust="BH", number=nrow(mymat), sort.by="P")

up <- rownames(res)[res$logFC > 0 & res$adj.P.Val <= 0.05]
down <- rownames(res)[res$logFC < 0 & res$adj.P.Val <= 0.05]

# late Vs. switch
contrast.matrix <- makeContrasts(late-switch, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

res <- topTable(fit2, adjust="BH", number=nrow(mymat), sort.by="P")

up <- rownames(res)[res$logFC > 0 & res$adj.P.Val <= 0.05]
down <- rownames(res)[res$logFC < 0 & res$adj.P.Val <= 0.05]
