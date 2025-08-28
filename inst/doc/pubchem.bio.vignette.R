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
# install.packages(pubchem.bio, dependencies = TRUE)
# library(pubchem.bio)
# local.pubchem.directory <- "C:/Temp/20250801"  ## i suggest this naming scheme, but feel free to chose your own.
# pubchem.bio::get.pubchem.ftp(pc.directory = local.pubchem.directory)
# cid.lca <- pubchem.bio::build.cid.lca(pc.directory = local.pubchem.directory, tax.sources =  "LOTUS - the natural products occurrence database")
# pc.bio <- pubchem.bio::build.pubchem.bio(pc.directory = local.pubchem.directory)
# pc.bio.salsa <- pubchem.bio::build.taxon.metabolome(taxid = c(4081, 4072, 4047, 4679, 4682), pc.directory = local.pubchem.directory, get.properties = FALSE)
# export.msfinder(pc.bio.object = pc.bio.salsa, export.file.name = paste0(local.pubchem.directory, "/pc.bio.msfinder.tsv"))
# export.sirius(pc.bio.object = pc.bio.salsa, export.file.name = paste0(local.pubchem.directory, "/pc.bio.sirius.tsv"))
# export.pubchem.bio(pc.bio.object = pc.bio.salsa, export.file.name = paste0(local.pubchem.directory, "/pc.bio.full.pubchem.bio.tsv"))

