## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----installation, eval=FALSE-------------------------------------------------
# install.packages(pubchem.bio, dependencies = TRUE)
# library(pubchem.bio)

## ----get.pubchem, eval=FALSE--------------------------------------------------
# pubchem.bio::get.pubchem.ftp(pc.directory ="C:/Temp/20250703")
# 

## ----build.cid.lca, eval=FALSE------------------------------------------------
# build.cid.lca(pc.directory = "C:/Temp/20250703", tax.sources =  "LOTUS - the natural products occurrence database")

## ----lca.demo, eval = TRUE----------------------------------------------------
data(sub.taxid.hierarchy, package = 'pubchem.bio')
sub.taxid.hierarchy

## ----build.pubchem.bio, eval=FALSE--------------------------------------------
# pc.bio <- build.pubchem.bio(pc.directory = "C:/Temp/20250703")

## ----build.taxon.metabolome, eval=FALSE---------------------------------------
# pc.bio.sub <- build.taxon.metabolome(taxid = 1710960, pc.directory = "C:/Temp/20250703", get.properties = FALSE)

## ----build.salsa.metabolome, eval=FALSE---------------------------------------
# pc.bio.sub <- build.taxon.metabolome(taxid = c(4072, 4107, 4047), pc.directory = "C:/Temp/20250703", get.properties = FALSE)

## ----lca.demo2, eval = TRUE---------------------------------------------------
data(pc.bio.subset, package = 'pubchem.bio')
pc.bio.subset[,c("cid", "name", "taxonomy.lca.similarity.4072", "taxonomy.lca.similarity.4107", "taxonomy.lca.similarity.4047", "taxonomy.lca.similarity.aggregate")]

## ----export.data, eval=FALSE--------------------------------------------------
# export.msfinder(pc.bio.object = pc.bio.sub, export.file.name = "C:/Temp/20250703/pc.bio.msfinder.tsv")
# export.sirius(pc.bio.object = pc.bio.sub, export.file.name = "C:/Temp/20250703/pc.bio.sirius.tsv")
# export.pubchem.bio(pc.bio.object = pc.bio.sub, export.file.name = "C:/Temp/20250703/pc.bio.full.pubchem.bio.tsv")

## ----all.steps, eval=FALSE----------------------------------------------------
# 
# library(pubchem.bio)
# local.pubchem.directory <- "C:/Temp/20250828"  ## i suggest this naming scheme, but feel free to chose your own.
# 
# ## run code, with timings recorded, b-a will give you the time to run the 'get.pubchem.ftp' function, etc.
# a <- Sys.time()
# pubchem.bio::get.pubchem.ftp(pc.directory = local.pubchem.directory, threads = 8)
# b <- Sys.time()
# 
# cid.lca <- pubchem.bio::build.cid.lca(pc.directory = local.pubchem.directory, tax.sources =  "LOTUS - the natural products occurrence database", threads = 8, min.taxid.table.length = 3)
# cid.lca[cid.lca$cid == 5793,]
# cid.lca[cid.lca$cid == 1548943,]
# cid.lca[cid.lca$cid == 44254980,]
# c <- Sys.time()
# 
# pc.bio <- pubchem.bio::build.pubchem.bio(pc.directory = local.pubchem.directory, get.properties = TRUE, threads = 8)
# d <- Sys.time()
# 
# pc.bio.tomato <- build.taxon.metabolome(taxid = c(4081), pc.directory = local.pubchem.directory, get.properties = FALSE)
# pc.bio.tomato[which(pc.bio.tomato == 5793),]
# pc.bio.tomato[which(pc.bio.tomato == 1548943),]
# pc.bio.tomato[which(pc.bio.tomato == 44254980),]
# e <- Sys.time()
# 
# pc.bio.salsa <- build.taxon.metabolome(taxid = c(4081, 4072, 4047, 4679, 4682), pc.directory = local.pubchem.directory, get.properties = FALSE)
# pc.bio.salsa[which(pc.bio.salsa == 5793),]
# pc.bio.salsa[which(pc.bio.salsa == 1548943),]
# pc.bio.salsa[which(pc.bio.salsa == 44254980),]
# f <- Sys.time()
# 
# load(paste0(local.pubchem.directory, "/cid.lca.Rdata"))
# load(paste0(local.pubchem.directory, "/cid.taxid.Rdata"))
# load(paste0(local.pubchem.directory, "/cid.monoisotopic.mass.Rdata"))
# load(paste0(local.pubchem.directory, "/pc.bio.Rdata"))
# load(paste0(local.pubchem.directory, "/cid.pwid.Rdata"))
# 
# primary.metabolome <- build.primary.metabolome(
#   pc.directory = local.pubchem.directory, get.properties = FALSE,
#   min.tax.ct = 2^5, keep.primary.only = TRUE)
# g <- Sys.time()
# 
# pubchem.bio::export.msfinder(pc.bio.object = pc.bio.salsa, export.file.name = paste0(local.pubchem.directory, "/pc.bio.msfinder.tsv"))
# pubchem.bio::export.sirius(pc.bio.object = pc.bio.salsa, export.file.name = paste0(local.pubchem.directory, "/pc.bio.sirius.tsv"))
# pubchem.bio::export.pubchem.bio(pc.bio.object = pc.bio.salsa, export.file.name = paste0(local.pubchem.directory, "/pc.bio.full.pubchem.bio.tsv"))
# h <- Sys.time()
# 
# out.dir <- "C:/Users/cbroeckl/OneDrive - Colostate/ARC/manuscripts/20250721_pubchem.bio/"
# 
# ggpl <- ggplot2::ggplot(pc.bio.salsa, ggplot2::aes(taxonomy.lca.similarity.aggregate)) +
#   ggplot2::geom_histogram(color = "#000000", fill = "#1E4D2B") +
#   ggplot2::theme_classic() +
#   ggplot2::annotate("text", x=0.7, y=75000, label= paste0("n = ", length(which(!is.na(pc.bio.salsa$taxonomy.lca.similarity.aggregate)))))
# 
# ggplot2::ggsave(ggpl, file = paste0(out.dir, "figure.2.eps"), scale = 0.75)
# 
# length(which(pc.bio.salsa$taxonomy.lca.similarity.4081 == 1))
# length(which(pc.bio.salsa$taxonomy.lca.similarity.4072 == 1))
# length(which(pc.bio.salsa$taxonomy.lca.similarity.4047 == 1))
# length(which(pc.bio.salsa$taxonomy.lca.similarity.4679 == 1))
# length(which(pc.bio.salsa$taxonomy.lca.similarity.4682 == 1))
# 
# ## percentage of tomato metabolites also in pepper
# tom <- which(pc.bio.salsa$taxonomy.lca.similarity.4081 == 1)
# pep <- which(pc.bio.salsa$taxonomy.lca.similarity.4072 == 1)
# length(which(pc.bio.salsa$cid[tom] %in% pc.bio.salsa$cid[pep]))/length(pc.bio.salsa$cid[tom])
# 
# ## percentage of tomato metabolites also in pepper
# garlic <- which(pc.bio.salsa$taxonomy.lca.similarity.4682 == 1)
# onion <- which(pc.bio.salsa$taxonomy.lca.similarity.4679 == 1)
# length(which(pc.bio.salsa$cid[garlic] %in% pc.bio.salsa$cid[onion]))/length(pc.bio.salsa$cid[garlic])
# 
# ## percentage of tomato metabolites also in garlic
# length(which(pc.bio.salsa$cid[tom] %in% pc.bio.salsa$cid[garlic]))/length(pc.bio.salsa$cid[tom])
# 
# 
# 
# 
# cids <- c(174174, 1548943) #atropine, capsaicin
# cid.lca[cid.lca$cid %in% cids,]
# 
# pc.bio.datura <- pubchem.bio::build.taxon.metabolome(taxid = c(4074), pc.directory = local.pubchem.directory, get.properties = FALSE)
# pc.bio.datura[pc.bio.datura$cid %in% cids,]
# 
# 
# n.tax <- 2^(seq(0, 10, 0.2))
# n.mets <- rep(NA, length(n.tax))
# for(i in 1:length(n.tax)) {
#   tmp <- build.primary.metabolome(
#     pc.directory = local.pubchem.directory, get.properties = FALSE,
#     min.tax.ct = n.tax[i], keep.primary.only = TRUE)
#   n.mets[i] <- nrow(tmp)
#   rm(tmp); gc()
# }
# 
# pc.mw <- ggplot2::ggplot(cid.monoisotopic.mass, ggplot2::aes(monoisotopic.mass)) +
#   ggplot2::geom_density(color = "#000000", fill = "#1E4D2B", bw = 50) +
#   ggplot2::theme_classic() +
#   ggplot2::scale_x_continuous(limits = c(0, 2000)) +
#   ggplot2::annotate("text", x=1600, y=0.0025, label= paste0("n = ", nrow(cid.monoisotopic.mass)))
# 
# pc.mw.bio <- ggplot2::ggplot(pc.bio, ggplot2::aes(monoisotopic.mass)) +
#   ggplot2::geom_density(color = "#000000", fill = "#1E4D2B", bw = 50) +
#   ggplot2::theme_classic() +
#   ggplot2::scale_x_continuous(limits = c(0, 2000)) +
#   ggplot2::annotate("text", x=1600, y=0.0011, label= paste0("n = ", nrow(pc.bio)))
# 
# pc.mw.prim <- ggplot2::ggplot(primary.metabolome, ggplot2::aes(monoisotopic.mass)) +
#   ggplot2::geom_density(color = "#000000", fill = "#1E4D2B", bw = 50) +
#   ggplot2::theme_classic() +
#   ggplot2::scale_x_continuous(limits = c(0, 2000)) +
#   ggplot2::xlab("monoisotopic mass") +
#   ggplot2::annotate("text", x=1600, y=0.002, label= paste0("n = ", nrow(primary.metabolome)))
# 
# 
# rarefaction <- data.frame(
#   taxa.count = n.tax,
#   primary.metabolites = n.mets
# )
# 
# pc.prim.curve <- ggplot2::ggplot(rarefaction, ggplot2::aes(x=taxa.count, y=primary.metabolites)) +
#   ggplot2::geom_line() +
#   ggplot2::scale_x_continuous(
#     trans = "log2",
#     labels = scales::math_format(2^.x, format = log2)
#   ) +
#   ggplot2::theme_classic()
# 
# 
# p1 <- pc.mw
# p2 <- pc.mw.bio
# p3 <- pc.prim.curve
# p4 <- pc.mw.prim
# 
# 
# fig.1 <- gridExtra::grid.arrange(
#   gridExtra::arrangeGrob(p1, left = grid::textGrob("a)", x = grid::unit(1, "npc"),
#                                        y = grid::unit(.95, "npc"))),
#   gridExtra::arrangeGrob(p2, left = grid::textGrob("b)", x = grid::unit(1, "npc"),
#                                   y = grid::unit(.95, "npc"))),
#   gridExtra::arrangeGrob(p3, left = grid::textGrob("c)", x = grid::unit(1, "npc"),
#                                   y = grid::unit(.95, "npc"))),
#   gridExtra::arrangeGrob(p4, left = grid::textGrob("d)", x = grid::unit(1, "npc"),
#                                   y = grid::unit(.95, "npc"))))
# ggplot2::ggsave(fig.1, file = paste0(out.dir, "figure.1.eps"), scale = 1.25)
# 
# pubchem.bio::export.pubchem.bio(pc.bio.object = primary.metabolome, export.file.name = "primary.metabolome.tsv")
# write.csv(primary.metabolome, file = "C:/Users/cbroeckl/OneDrive - Colostate/ARC/manuscripts/20250721_pubchem.bio/primary.metabolome.csv", row.names = FALSE)
# 

