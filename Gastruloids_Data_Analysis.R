################################################# ~ Good Luck ~ #########################################################
###### >>>>>>                          TITLE: Gastruloids project: Wang Bin                                 <<<<<< ######



### ========================
### 1st part: Global setting ----
### ========================

### >>> 1. data saving
setwd("/home/yhw/bioinfo/project-wb")
load("CodeData/WB_downstream_v1.RData")
#save.image("scripts/downstream_v1.RData")
.libPaths("/home/laborer/yhw_shared/library")


### >>> 2. packages loading
library("dplyr")
library("tidyr")
library("ggplot2")
library("stringr")
library("Seurat")
library("tibble")
library("SeuratDisk")
library("harmony")
library("moonBook")
library("webr")


### >>> 3. define output directory
sr.out <- "/home/yhw/bioinfo/project-wb/hs_gastruloids/analysis/graphs/seurat_output/"


### >>> 4. color setting
# - for human embryo
# raw
pd.col <- c("#449945", "#ea7827", "#c22f2f", "#1f70a9", "#83639f",
            "#e9963e", "#f23b27", "#65a9d7", "#304f9e",
            "#7d3c98", "#F08080")
c("#58539f","#bbbbd6", "#eebabb", "#d86967")
# with name
pd.col.fix <- c("#54beca", "#3b77b3", "#b0c6e4", "#cb7eaf",
                "#845a4d", "#eb8838", "#b7be70", "#4f996e",
                "#f2adad", "#6d51a0", "#c33b30")
names(pd.col.fix) <- c("ECT", "EPI", "PS", "NM",
                       "AxM", "End", "EM", "AM",
                       "ExM", "HEP", "Ery")
names(pd.col.fix) <- c("Ectoderm", "Epiblast", "Primitive streak", "Nascent mesoderm",
                       "Axial Mesoderm", "Endoderm", "Emergent Mesoderm", "Advanced Mesoderm",
                       "Extraembryonic mesoderm", "Haemato-endothelial progenitors", "Erythroblasts")
# with name from wb
pd.col.fix <- c("#54beca", "#b0c6e4", "#cb7eaf", "#b7be70", "#845a4d", "#eb8838",
                "#4f996e", "#3b77b3", "#f2adad", "#6d51a0", "#c33b30")
names(pd.col.fix) <- c("ECT", "EPI", "PS", "NM","AxM", "End",
                       "EM", "AM","ExM", "HEP", "Ery")
# - for monkey embryo
mk.col <- c("#5e5343", "#d3ba9b", "#906257", "#edc34d", "#ba9b70",
            "#20467f", "#de8fb7", "#ce6595", "#b48eb7", "#449089",
            "#c9e4f5", "#7e6ca8", "#bf7b30", "#ecb690", "#f4ddce",
            "#c2a893", "#3a4b30", "#d9987b", "#d2543d", "#c2c096",
            "#6a7a56", "#9c352c", "#db873e", "#f0eeae", "#e8bdca",
            "#7b6473", "#686868", "#969696", "#1a1a1a")
names(mk.col) <- c("EPI", "PS", "PGC", "APS", "Node",
                   "DE", "Gut", "Nas.Meso", "LP.Meso", "Caud.Meso",
                   "Pharyg.Meso", "Al", "ExE.Meso", "ys.Meso1", "ys.Meso2",
                   "Mes", "EC", "BP", "Mac", "Ery1",
                   "Ery2", "ECT", "NC", "SE1", "SE2",
                   "AM", "VE", "ys.Endo1", "ys.Endo2")


### >>> 5. marker gene list
# - human CS7 markers from Nature paper
mk.ge.hs <- list()
mk.ge.hs$EAE <- c("DLX5","TFAP2A","GATA3")
mk.ge.hs$EPI <- c("SOX2","OTX2")
mk.ge.hs$PS <- c("TBXT","CDH1","FST")
mk.ge.hs$NM <- c("TBXT","MESP1","PDGFRA","TBX6")
mk.ge.hs$AxM <- c("TBXT","CHRD","NOTO")
mk.ge.hs$Endo <- c("SOX17","GATA6","FOXA2","TTR")
mk.ge.hs$EM <- c("MESP1","LHX1","LEFTY2")
mk.ge.hs$AM <- c("PDGFRA","GATA6","HAND1","BMP4","FOXF1","SNAI2")
mk.ge.hs$ExM <- c("POSTN","ANXA1")
mk.ge.hs$HEP <- c("PECAM1","MEF2C","RUNX1")
mk.ge.hs$Ery <- c("HBZ","HBE1","GATA1")
# - human embryo markers from wb
mk.ge.hs$wb <- data.frame(Lineage=c(rep("EAE",4), rep("EPI",4), rep("PS", 4), rep("NM", 4), rep("AxM", 3), rep("Endo", 5), rep("EM", 5),
                                    rep("AM", 6), rep("ExM",6), rep("HEP", 6), rep("Ery", 2)),
                          Marker=c("TFAP2A","GATA3","DLX5","TPM1",
                                   "SOX2","POU5F1","OTX2","CDH1",
                                   "TBXT","SP5","FST","CDX1",
                                   "MESP1","TBX6","PDGFRA","SFRP2",
                                   "TBXT","CHRD","NOTO",
                                   "SOX17","FOXA2","CST1","TTR","AFP",
                                   "LHX1","GATA6","IRX3","FOXC1","LEFTY2",
                                   "MSX1","FOXF1","HAND1","BMP4","TNNI1","MYL7",
                                   "POSTN","ANXA1","LUM","NID2","DCN","MAB21L2",
                                   "MEF2C","PECAM1","RUNX1","ITGA2B","CDH5","TEK",
                                   "HBZ","GATA1"))
# - monkey CS8 makers from Nature paper
mk.ge.mk <- list()
mk.ge.mk$EPI <- c("POU5F1", "SOX2", "DNMT3B", "POU3F1", "UCHL1")
mk.ge.mk$PS <- c("TBXT", "POU5FA", "WNT5B", "CDX2", "CDX1")
mk.ge.mk$PGC <- c("DPPA3", "NANOG", "NANOS3", "SOX17", "TFAP2C", "PRDM14")
mk.ge.mk$APS <- c("CER1", "FOXA2", "HHEX")
mk.ge.mk$Node <- c("FOXA2", "CHRD", "NOTO", "TBXT")
mk.ge.mk$DE <- c("SOX17", "FOXA2", "OTX2", "HHEX", "CER1", "EOMES")
mk.ge.mk$Gut <- c("CDX2", "ISL1", "CDH2", "OSR1", "SOX17", "FOXA2")
mk.ge.mk$Nas.Meso <- c("POU5FA", "LEFTY2", "TBXT", "MESP1", "TBX6", "FOXC1", "FOXC2")
mk.ge.mk$LP.Meso <- c("GATA6", "IRX3", "FOXC1", "FOXF1", "GATA4")
mk.ge.mk$Caud.Meso <- c("TBXT", "CDX1", "WNT3A", "HES7", "EOMES", "CDX2", "MIXL1")
mk.ge.mk$Pharyg.Meso <- c("GATA6", "TCF21", "HAND1", "MYL7")
mk.ge.mk$Al <- c("FOXF1", "HAND1", "COL3A1", "COL1A1", "ISL1", "TCF21")
mk.ge.mk$ExE.Meso <- c("HAND1", "CDX2", "FOXF1", "IGFBP5", "IRX3", "ISL1", "GATA6")
mk.ge.mk$ys.Meso1 <- c("ANXA1", "ANXA8", "COL3A1", "COL1A1", "HAND1", "APOE")
mk.ge.mk$ys.Meso2 <- c("ANXA1", "ANXA8", "COL3A1", "COL1A1", "HAND1", "APOE")
mk.ge.mk$Mes <- c("COL3A1", "COL1A1", "ANXA8")
mk.ge.mk$EC <- c("PECAM1", "KDR", "VWF", "HHEX")
mk.ge.mk$BP <- c("PECAM1", "PPBP", "GP1BB", "MPO", "HBZ", "HBM")
mk.ge.mk$Mac <- c("FOLR2", "MPO")
mk.ge.mk$Ery1 <- c("HBZ", "HBM")
mk.ge.mk$Ery2 <- c("HBZ", "HBM")
mk.ge.mk$ECT <- c("POU5FA", "SOX2", "DNMT3B", "OTX2")
mk.ge.mk$NC <- c("SOX10", "SOX9", "PAX3")
mk.ge.mk$SE1 <- c("SOX9", "TFAP2A", "KRT18", "ISL1", "KRT8")
mk.ge.mk$SE2 <- c("TFAP2A", "KRT18", "TFAP2C")
mk.ge.mk$AM <- c("TFAP2A", "SIX1", "TFAP2C", "POU5F1", "SOX2")
mk.ge.mk$VE <- c("GATA4", "GATA6", "FOXA2", "SOX17", "APOE", "TTR")
mk.ge.mk$ys.Endo1 <- c("APOE", "AFP", "TTR", "KRT18", "KRT8")
mk.ge.mk$ys.Endo2 <- c("APOE", "AFP", "TTR")
marker.gene <- list(Organizer = c("FST", "OTX2", "GSC", "EOMES", "FOXH1", "CER1", "MIXL1", "DKK1", "NODAL",
                                  "FOXA2", "LEFTY1", "LEFTY2", "T", "SP5", "NOG", "CHRD"),
                    Buchong = c("NANOG", "HHEX", "CDH1", "POU3F1", "CDX2", "CDX1", "SOX17",
                                "GATA6", "GATA4", "LUM", "NID2", "POU5F1", "SOX2", "MYC",
                                "LIN28", "GDF3", "TDGF1", "SNAIL", "SNAIL2", "TBX6",
                                "MESP1", "LHX1", "IRX3", "HAND1", "FOXF1"),
                    FMH.marker  = c("LHX5","ZNF616","GBX2","RREB1","TEAD4","EN2","HAND1",
                                    "HESX1","FOXG1","NFATC4","PAX7","DBX1","PAX5","SMARCA5",
                                    "PAX6","IRX3","EBF1","HOXA3","MEIS1","WRNIP1","MAFB"),
                    FMH = c("PAX6","OTX2","EN1","PAX8","GBX2","EGR2","PAX3","DBX2","DBX1","IRX3","OLIG2","NKX2.2"),
                    Neural = c("SOX1","SOX3","PAX6","TUBB3","ELAVL3","OLIG2","NEUROG1","NEUROD1","NIKX2.2"),
                    NMP = unique(c("LIMD1","MARK3","CIT","MAPK14","FAT4","LATS2","STK4","WWC1","SHANK2","SOX11","TEAD1","PJA2",
                                   "TEAD2","VGLL5","SAV1","MAP2K3","TIAL1","NF2","STK3","LATS1","THBS1","TEAD3","AMOT","YAP1",
                                   "AMOTL1","AMOTL2","NEK8","WTIP","DLG5","TEAD4",
                                   "NKX1-2","SP5","FGF17","WNT3A","CDX1","CDX2","FGF8","HES7","CDH2",
                                   "TBX6","FGF15","WNT5B","WNT10A","DUSP6","AXIN2","HOXA5","HOXB8","HOXA13",
                                   "HOXC13","POU5F1","SALL4","TET1","TET2","TET3","MED12","WNT8C","WNT8A",
                                   "CTNNB1","SP8","AXIN2","TCF1","LEF1","WNT3","VANGL2","WNT5A","WNT11","FGF4",
                                   "FGFR1","ALDH1A2","RARB","DELTA1","GFD8","GDF11","BMP4","POU3F1","UTF1","FOXA2",
                                   "CER1","LEFTY1","LEFTY2","MESP1","TBX3","EOMES","GATA6","MEIS1","OTX2","RSPO3",
                                   "HOXA1","IDL1","HOXA1","MIXL1","CER1","IDI1","HMGCS1","UCHL1","PLP1","KISS1","HES6",
                                   "LPAR6","CTGF","LMO2","CCNB2","EBPL","TBXT","SOX2","SALL1","BASP1","VCAN","FASN",
                                   "ID3","AURKA","MGST1","ALPL","BCAT1","SFRP2","MSX1","HOXA9","DIT4","PENK","EGR1",
                                   "MLLT3","FOSB","GFP")),
                    Amnion = c("GABRP","TGFB1","VTCN1","S100A10","GADD45G","ACTC1",
                               "IGFBP3","MCM5","CLDN10","TUBG1","UCP2","BIN3","GGA2","FGFR1","GALM"),
                    EPI = c("CDH1", "POU5F1", "SOX2", "NANOG", "LIN28AP1", "MYC", "KLF4", "TDGF1", "GDF3", "GFP"),
                    Embryonic = c("TBXT", "MESP1", "PDGFRA", "LHX1", "IRX3", "GATA6", "MYL7", "MSX1",
                                  "MSX2", "FOXF1", "HAND1", "HAND2", "BMP4", "MAB21L2", "CST1", "MEF2C",
                                  "PECAM1", "CDH5", "TEK", "TFAP2A", "TFAP2C", "GATA3", "DLX5", "LHX5",
                                  "PAX6", "SOX3", "PAX3", "PAX7"),
                    Extraembryonic = c("LUM", "NID2", "ANXA1", "POSTN", "GATA2", "KRT7", "VIM", "GATA4",
                                       "DCN", "TP63", "KRT18", "NR2F2", "ISL1", "HEY1", "CDH10", "CTSV",
                                       "ARID5B", "PLAGL1", "CREB3L1", "HOXA10", "HOXA11", "HOXA9",
                                       "HOXA13", "WNT6", "GABRP", "HEY1", "BST2"),
                    anterior = c("POU5F1","NANOG","TDGF1","NES","ID1","CDH1","SOX2","POU3F1","SOX3","IRX2",
                                 "CRABP1","PBX1","BTG1","TFAP2A","SFRP1","SFRP2","OTX2","FGFR1","HESX1","LHX5"),
                    medial = c("TBXT","TBX6","MESP1","EOMES","IRX3","LHX1","FOXC1","HOXB1","HOXA1","MEIS2","HOXB2",
                               "HOXB3","HOXA3","DLL1","GATA6","GATA4","FST","MIXL1","SP5","CER1","GFP"),
                    posterior = c("CDX2","CDX4","HOXB9","HOXC4","HOXC8","HOXC6","HOXB7","HOXB8","HOXB6","HOXB5",
                                  "VIM","SHH","FOXA2","CDH2","DUSP6","HOXA10","HOXA13","HAND1","HAND2","MSX1","GFP"),
                    megakaryocyte = c("GP1BB", "ITGA2B", "NFE2"),
                    erythroid = c("GATA1", "KLF1", "GYPB", "HBE1"),
                    myeloid.progenitor = c("CD36", "CSF1R", "LYVE1", "KIT", "CSF1R", "MYB", "SPI1", "CD34", "CD45", "CD52", "NFE2"),
                    wanghongmei = c("AFP","ANXA1","ANXA8","APOE","APOE ","CDH1","CDH2","CDX1","CDX2","CER1","CHRD",
                                    "COL1A1","COL3A1","DNMT3B","DPPA3","EOMES","FGF8","FGF8 HES7","FOLR2","FOXA2",
                                    "FOXC1","FOXC2","FOXF1","GATA1","GATA4","GATA6","GP1BB","HAND1","HBE1","HBM",
                                    "HBZ","HES7","HHEX","IGFBP5","IRX3","ISL1","KDR","KRT18","KRT8","LEFTY2","LFNG",
                                    "MEF2C","MIXL1","MPO","MYL7","NANOG","NANOS3","NKX2-5","NOTO","OSR1","OTX2","PAX3",
                                    "PAX3 ","PAX6","PAX7","PAX8","PECAM1","POU3F1","POU5F1","POU5FA","PPBP","PRDM14",
                                    "RUNX1","SIX1","SOX10","SOX17","SOX2","SOX2 ","SOX9","TBXT","TBX6","TCF15","TCF21",
                                    "TFAP2A","TFAP2C","TFC15","MESP1","TTR","UCHL1","VWF","WNT3A","WNT5B","WT1"),
                    neuron_NMP = c("CYSTM1","EM6","GRHL3","PAX6","PAX3","PAX2"),
                    meso_NMP = c("IRX3","HOXB9","CDX4","EPHA5","HES3","HBB-BH1","LMO2","ANXA5",
                                 "PMP22","GBX2","HOXB1","APELA","DNMT3B","SNAIL","SNAIL2","TTR"),
                    GUT = c("SHH","ISL1","IHH","ZG16","CDX1","URAD","TACSTD2","CDX2","TNNC1","COL9A2","MSX2","TBX3",
                            "CDH2","ADD3","EPCAM","CTHRC1","EBPL","CTNNBL1","FN1","PENK","FDX1","CXCL12","RBP1",
                            "CTSV","FOXJ1","C15H9ORF116","EGFL6","BMP7","S100A16","S100A13","TTR","AFP","GJB1",
                            "HHEX","SHISA2","CXCR4","NKX2.1","PPY","IRX1","APOA2","AKR1D1","APOB","VTN","DPYS",
                            "VEPH1","SERPINE2","DENND2C","S100A11","COL4A1","SPOCK3","AC117945.2","SIX3","CENPU",
                            "TK1","HMGN2","PCLAF","TMSB15A","SFRP5","APOM","APOC1","MEST","TEKT1","S100A1","VEGFA",
                            "LAMA1","ZIM2","LARP7","GJA1","FOS","UBE2C","COL18A1","IGDCC3","CLDN4","SLC39A2",
                            "LPAR6","BAMBI","CDKN1C","NKX2-5","HHEX","IRX2","OSR1"),
                    ExE.meso = c("PENK","APLNR","HAS2","HAPLN1","TMEM88","PMP22","BAMBI","PITX1","CDH11","KDR"),
                    Cardiac.meso = c("NKX2.5","TBX5","TNNI1","TNNT2","MYH10","MAB21I2","HOPX","TGFBI"),
                    NSCs = c("HOXB1", "HOXB2", "HOXB5", "HOXB6", "HOXB8", "HOXB9", "SOX2", "SOX10", "SOX9", "FOXD3", "TJP1", "NESTIN"),
                    Neuron.epi = c("NES", "NOTCH1", "HES1", "BMI1", "TBRG1", "FABP7", "NEUROD1"),
                    Neuron.crest = c("NGFR", "HNK1", "CD49D", "ZIC1", "SNAI1"), 
                    buchong = c("PAX6","SOX2","SOX3","HES1","HESX1","IRX2","OTX2","LHX5","NES","TJP1","SIX1","SIX2","SIX3","SIX5",
                                "SIX6","TBRG1，PAX3","TDGF1","CDX2","CDX4","APLNR","HAND1","MSX1","MSX2","TFAP2A","ZIC3","HAND2",
                                "PMP22","SNAI1","HOXA3","HOXA10","HOXA13","HOXB2","HOXB6","HOXB7","HOXB8","HOXB9","HOXC4","HOXC6",
                                "HOXC8","DMBX1","BRN3C","SOX5","SOX6","LMX1B","SOX8","TWIST","TFAP2B","PAX7","PITX2","MYF5","MYF6",
                                "MYOD1","MYOG","NEB","MYH3","SOX9","SOX10","ETS1"))



### =============================
### 2nd part: Data pre-processing ----
### =============================

### >>> 1. load count table
ge.co <- list()
# data1
ge.co$d1 <- read.csv("data1/output/Combined_Sample_RSEC_ReadsPerCell_ForR.csv", stringsAsFactors = F, row.names = 1)
ge.co$d1 <- t(ge.co$d1)
ge.co$d1 <- subset(ge.co$d1, rowSums(ge.co$d1)>0)
# data2
ge.co$d2 <- read.csv("data2/output/Combined_data2Sample_RSEC_ReadsPerCell_ForR.csv", stringsAsFactors = F, row.names = 1)
ge.co$d2 <- t(ge.co$d2)
ge.co$d2 <- subset(ge.co$d2, rowSums(ge.co$d2)>0)
dim(ge.co$d2)
dim(ge.co$d1)

### >>> 2. load cell metadata
cell.meta <- list()
cell.meta$d1 <- read.csv("data1/output/Sample_Sample_Tag_Calls_ForR.csv", stringsAsFactors = F, row.names = 1)
cell.meta$d2 <- read.csv("data2/output/data2Sample_Sample_Tag_Calls_ForR.csv", stringsAsFactors = F, row.names = 1)

### >>> 3. create Seurat object
# - raw
sr <- list()
# data1
sr$d1 <- CreateSeuratObject(counts = ge.co$d1, project = "Data1", meta.data = cell.meta$d1) %>%
  PercentageFeatureSet(col.name = "Pect.mt", pattern = "^MT\\.")
pdf(paste0(sr.out,"data1_quality_raw_vlnplot.pdf"), height = 4, width = 8)
VlnPlot(sr$d1, features = c("nFeature_RNA", "nCount_RNA", "Pect.mt"), ncol = 3, group.by = "Sample_Name")
dev.off()
# data2
sr$d2 <- CreateSeuratObject(counts = ge.co$d2, project = "Data2", meta.data = cell.meta$d2) %>%
  PercentageFeatureSet(col.name = "Pect.mt", pattern = "^MT\\.")
pdf(paste0(sr.out,"data2_quality_raw_vlnplot.pdf"), height = 4, width = 8)
VlnPlot(sr$d2, features = c("nFeature_RNA", "nCount_RNA", "Pect.mt"), ncol = 3, group.by = "Sample_Name")
dev.off()

# - filtering cells by quality
sr <- merge(sr$d1, sr$d2)
sr <- subset(sr, nFeature_RNA >= 500 & nCount_RNA >= 5000 & Pect.mt <= 25)
sr <- subset(sr, Sample_Name%in%c("EPI_D0","GAS_D1", "GAS_D2", "GAS_D3", "GAS_D5", "ORG_D1", "ORG_D2", "ORG_D3"))
# plotting
pdf(paste0(sr.out,"all_sample_quality_after_filtering_vlnplot.pdf"), height = 4, width = 8)
VlnPlot(sr, features = c("nFeature_RNA", "nCount_RNA", "Pect.mt"), ncol = 3, group.by = "Sample_Name")
dev.off()
# add GFP count to metadata
tmp <- GetAssayData(sr, slot = "count")
tmp <- tmp["GFP",]
if (identical(names(tmp),colnames(sr))) {
  sr$GFP_count <- as.vector(tmp)
}; rm(tmp)
rm(ge.co)

# - visualization of GFP expression level
sr.list <- list()
# cell line
sr.list$cl <- subset(sr, Sample_Name%in%c("EPI_D0","ORG_D1", "ORG_D2", "ORG_D3"))
sr.list$cl <- NormalizeData(sr.list$cl) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sr.list$cl))
sr.list$cl <- RunPCA(sr.list$cl, features = VariableFeatures(object = sr.list$cl))
ElbowPlot(sr.list$cl)
sr.list$cl <- RunUMAP(sr.list$cl, dims = 1:10) %>% RunTSNE(dims = 1:10)
pdf(paste0(sr.out,"cell_line_umap_tsne_plot.pdf"), height = 4, width = 10)
p1 <- DimPlot(sr.list$cl, reduction = "umap", group.by = c("Sample_Name"), cols = pd.col)
p2 <- DimPlot(sr.list$cl, reduction = "tsne", group.by = c("Sample_Name"), cols = pd.col)
p1 + p2
dev.off(); rm(p1, p2)
pdf(paste0(sr.out,"GFP_expression_in_cell_line_featureplot.pdf"), height = 4, width = 4.5)
FeaturePlot(sr.list$cl, reduction = "umap", features = c("GFP"), cols = c("#ebebeb", "#7d3c98"))
dev.off()
pdf(paste0(sr.out,"GFP_expression_in_cell_line_vlnplot.pdf"), height = 3, width = 4)
VlnPlot(sr.list$cl, features = "GFP", group.by = "Sample_Name", cols = pd.col)
dev.off()
tmp <- subset(sr.list$cl, Sample_Name%in%c("EPI_D0", "ORG_D2"))
tmp@meta.data %>% mutate(Type=case_when(GFP_count > 0 ~ "GFP_pos", TRUE ~ "GFP_neg")) -> tmp@meta.data
table(tmp@meta.data[,c("Sample_Name","Type")])
tmp <- NormalizeData(tmp) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(tmp))
tmp <- RunPCA(tmp, features = VariableFeatures(object = tmp))
ElbowPlot(tmp)
tmp <- RunUMAP(tmp, dims = 1:5) %>% RunTSNE(dims = 1:5)
sr.out <- "/home/yhw/bioinfo/project-wb/hs_gastruloids/analysis/graphs/seurat_output/GASD0_marker_gene"
dir.create(sr.out, recursive = T)
pdf(file.path(sr.out, "cell_line_umap_tsne_plot.pdf"), height = 4, width = 10)
p1 <- DimPlot(tmp, reduction = "umap", group.by = c("Sample_Name"), cols = c("#B91122","#4CAF49"))
p2 <- DimPlot(tmp, reduction = "tsne", group.by = c("Sample_Name"), cols = c("#B91122","#4CAF49"))
p1 + p2
dev.off(); rm(p1, p2)
marker.gene <- list(EPI = c("CDH1", "POU5F1", "SOX2", "NANOG", "LIN28AP1", "MYC", "KLF4", "TDGF1", "GDF3", "GFP"),
                    Organizer = c("FST", "OTX2", "GSC", "EOMES", "FOXH1", "CER1", "MIXL1", "DKK1", "NODAL",
                                  "FOXA2", "LEFTY1", "LEFTY2", "TBXT", "SOX17", "SP5", "NOG", "CHRD"),
                    Embryonic = c("TBXT", "MESP1", "PDGFRA", "LHX1", "IRX3", "GATA6", "MYL7", "MSX1",
                                  "MSX2", "FOXF1", "HAND1", "HAND2", "BMP4", "MAB21L2", "CST1", "MEF2C",
                                  "PECAM1", "CDH5", "TEK", "TFAP2A", "TFAP2C", "GATA3", "DLX5", "LHX5",
                                  "PAX6", "SOX3", "PAX3", "PAX7"),
                    Extraembryonic = c("LUM", "NID2", "ANXA1", "POSTN", "GATA2", "KRT7", "VIM", "GATA4",
                                       "DCN", "TP63", "KRT18", "NR2F2", "ISL1", "HEY1", "CDH10", "CTSV",
                                       "ARID5B", "PLAGL1", "CREB3L1", "HOXA10", "HOXA11", "HOXA9",
                                       "HOXA13", "WNT6", "GABRP", "HEY1", "BST2"))
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("GASD0_marker_gene_expression_of_", marker.gene[[type]][g], "_featureplot.pdf")),
        height = 3.5, width = 4)
    print(FeaturePlot(tmp, features = marker.gene[[type]][g],
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 0.5, slot = "data"))
    dev.off()
  }; rm(g)
}
# GAS
sr.list$gas <- subset(sr, Sample_Name%in%c("GAS_D1", "GAS_D2", "GAS_D3", "GAS_D5"))
sr.list$gas <- NormalizeData(sr.list$gas) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sr.list$gas))
sr.list$gas <- RunPCA(sr.list$gas, features = VariableFeatures(object = sr.list$gas))
ElbowPlot(sr.list$gas)
sr.list$gas <- RunUMAP(sr.list$gas, dims = 1:15) %>% RunTSNE(dims = 1:15)
pdf(paste0(sr.out,"GAS_umap_tsne_plot.pdf"), height = 4, width = 10)
p1 <- DimPlot(sr.list$gas, reduction = "umap", group.by = c("Sample_Name"), cols = pd.col)
p2 <- DimPlot(sr.list$gas, reduction = "tsne", group.by = c("Sample_Name"), cols = pd.col)
p1 + p2
dev.off(); rm(p1, p2)
pdf(paste0(sr.out,"GFP_expression_in_GAS_featureplot.pdf"), height = 4, width = 4.5)
FeaturePlot(sr.list$gas, reduction = "umap", features = c("GFP"), cols = c("#ebebeb", "#7d3c98"))
dev.off()
pdf(paste0(sr.out,"GFP_expression_in_GAS_vlnplot.pdf"), height = 3, width = 4)
VlnPlot(sr.list$gas, features = "GFP", group.by = "Sample_Name", cols = pd.col)
dev.off()
# cell line2
sr.list$cl2 <- subset(sr, Sample_Name%in%c("EPI_D0","ORG_D2"))
sr.list$cl2 <- NormalizeData(sr.list$cl2) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sr.list$cl2))
sr.list$cl2 <- RunPCA(sr.list$cl2, features = VariableFeatures(object = sr.list$cl2))
ElbowPlot(sr.list$cl2)
sr.list$cl2 <- RunUMAP(sr.list$cl2, dims = 1:10) %>% RunTSNE(dims = 1:10)
pdf(paste0(sr.out,"EPI_D0_and_ORG_D2_umap_tsne_plot.pdf"), height = 4, width = 10)
p1 <- DimPlot(sr.list$cl2, reduction = "umap", group.by = c("Sample_Name"), cols = c("#C9DDFB","#8F5090"), pt.size = 1.25)
p2 <- DimPlot(sr.list$cl2, reduction = "tsne", group.by = c("Sample_Name"), cols = c("#C9DDFB","#8F5090"), pt.size = 1.25)
p1 + p2
dev.off(); rm(p1, p2)
pdf(paste0(sr.out,"GFP_expression_in_EPI_D0_and_ORG_D2_featureplot.pdf"), height = 4, width = 4.5)
FeaturePlot(sr.list$cl2, reduction = "umap", features = c("GFP"), cols = c("#ebebeb", "#7d3c98"), pt.size = 1.5)
dev.off()
pdf(paste0(sr.out,"GFP_expression_in_EPI_D0_and_ORG_D2_vlnplot.pdf"), height = 3, width = 4)
VlnPlot(sr.list$cl2, features = "GFP", group.by = "Sample_Name", cols = pd.col)
dev.off()
saveRDS(sr.list$cl2, "/home/cmq/bioinfo/project-tmp/wb/analysis/rds/epi.d0.org.d2.sr.rds")


### >>> 4. pre-processing of GAS-D5 data
sr.gas <- list()
sr.gas$d5 <- subset(sr.list$gas, Sample_Name=="GAS_D5") %>% NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
sr.gas$d5 <- ScaleData(sr.gas$d5, features = rownames(sr.gas$d5)) %>%
  RunPCA(features = VariableFeatures(sr.gas$d5))
ElbowPlot(sr.gas$d5)
sr.gas$d5 <- RunUMAP(sr.gas$d5, dims = 1:10) %>% RunTSNE(dims = 1:10) %>% FindNeighbors(dims = 1:10)
sr.gas$d5@meta.data %>%
  mutate(Type = case_when(GFP_count > 0 ~ "GFP_pos", TRUE ~ "GFP_neg")) -> sr.gas$d5@meta.data
pdf(paste0(sr.out,"GFP_expression_in_GAS_D5_featureplot.pdf"), height = 4, width = 4.5)
FeaturePlot(sr.gas$d5, reduction = "umap", features = c("GFP"), cols = c("#ebebeb", "#7d3c98"))
dev.off()
pdf(paste0(sr.out,"GAS_D5_umap_tsne_plot.pdf"), height = 4, width = 10)
p1 <- DimPlot(sr.gas$d5, reduction = "umap", group.by = c("Type"), cols = c("#1eaf90","#fa3f89"))
p2 <- DimPlot(sr.gas$d5, reduction = "tsne", group.by = c("Type"), cols = c("#1eaf90","#fa3f89"))
p1 + p2
dev.off(); rm(p1, p2)
sr.gas$d5$DataSet <- "GAS"


### >>> 5. load public dataset
# - human CS7 data
hs.cs7 <- readRDS("/home2/data/publicdata/scRNAseq/E-MTAB-9388/analysis2/results/r/rds/sr.rds")
DefaultAssay(hs.cs7) <- 'RNA'
hs.cs7@assays <- hs.cs7@assays[1]
hs.cs7 <- hs.cs7 %>% FindVariableFeatures(selection.method = "vst", nfeatures = 4000)
hs.cs7 <- hs.cs7 %>% RunPCA(features = VariableFeatures(hs.cs7))
ElbowPlot(hs.cs7, ndims = 50)
hs.cs7 <- RunUMAP(hs.cs7, dims = 1:18, return.model = T) %>% RunTSNE(dims = 1:18) %>% FindNeighbors(dims = 1:18)
hs.cs7@meta.data <- hs.cs7@meta.data[,-c(11:16)]
hs.cs7@meta.data <- hs.cs7@meta.data %>%
  mutate(AuthorLabel=case_when(AuthorLabel=="ectodermal cell" ~ "ECT",
                               AuthorLabel=="epiblast cell" ~ "EPI",
                               AuthorLabel=="primitive streak" ~ "PS",
                               AuthorLabel=="nascent mesoderm" ~ "NM",
                               AuthorLabel=="axial mesoderm" ~ "AxM",
                               AuthorLabel=="endodermal cell" ~ "End",
                               AuthorLabel=="emergent mesoderm" ~ "EM",
                               AuthorLabel=="advanced mesoderm" ~ "AM",
                               AuthorLabel=="yolk sac mesoderm" ~ "ExM",
                               AuthorLabel=="hemogenic endothelial progenitor" ~ "HEP",
                               AuthorLabel=="erythrocyte" ~ "Ery"))
hs.cs7$CellType <- factor(hs.cs7$AuthorLabel,
                          levels = c("ECT", "EPI", "PS", "NM",
                                     "AxM", "End", "EM", "AM",
                                     "ExM", "HEP", "Ery"))
pdf(paste0(sr.out, "Human_CS7_umap_tsne_plot.pdf"), height = 4, width = 9.5)
p1 <- DimPlot(hs.cs7, reduction = "umap", group.by = "CellType", cols = pd.col.fix)
p2 <- DimPlot(hs.cs7, reduction = "tsne", group.by = "CellType", cols = pd.col.fix)
p1 + p2
dev.off(); rm(p1, p2)
pdf(paste0(sr.out,"Human_CS7_marker_genes_dotplot.pdf"), height = 3.25, width = 13)
DotPlot(hs.cs7, features = unique(mk.ge.hs$wb$Marker), group.by = "CellType",
        col.min = -2, col.max = 2.5) +
  scale_color_gradient2(high = "#cd120c", mid = "#ffffff", low = "#04579c", limits = c(-2, 2.5)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        panel.grid = element_blank(),
        panel.background= element_rect(color = "black", fill = NA, linewidth = 1.25),
        axis.text.y = element_text(size = 12),
        legend.box = "horizontal")
dev.off()
hs.cs7$DataSet <- "NE"
# plotting
Vis.Annotation.Ratio(sr.meta = hs.cs7@meta.data, annotation = c("CellType", "Tissue"),
                     pd.title = "Stat by Cell Type", pd.col = NULL,
                     seq.x = c("caudal", "rostral", "yolk sac"),
                     seq.y = c("AxM", "NM", "PS", "EM", "AM", "EPI", "ECT", "End", "ExM", "HEP", "Ery"),
                     pd.height = 15, pd.width = 15,
                     res.out = file.path(sr.out, "composition"))

# - monkey CS8 data
mk.cs8 <- readRDS("/home/yhw/bioinfo/project-mine/MultiOmics/R/CodeData/Monkey_embryo.rds")
tmp1 <- GetAssayData(mk.cs8, assay = "RNA", slot = "count")
rownames(tmp1)[grep("^T$", rownames(tmp1))] <- "TBXT"
tmp2 <- mk.cs8@meta.data
mk.cs8 <- CreateSeuratObject(counts = tmp1, meta.data = tmp2); rm(tmp1, tmp2)
# subset data
mk.cs8 <- SampleSeurat(mk.cs8, "cell_type", 0.4, 50)
table(mk.cs8$cell_type)
mk.cs8 <- mk.cs8 %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(mk.cs8))
mk.cs8 <- mk.cs8 %>% RunPCA(features = VariableFeatures(mk.cs8))
ElbowPlot(mk.cs8, ndims = 30)
mk.cs8 <- RunUMAP(mk.cs8, dims = 1:16, return.model = T) %>% RunTSNE(dims = 1:16) %>% FindNeighbors(dims = 1:16)
mk.cs8@meta.data <- mk.cs8@meta.data %>% mutate(cell_type=case_when(cell_type=="AI" ~ "Al", TRUE ~ cell_type))
mk.cs8$CellType <- factor(mk.cs8$cell_type, levels = names(mk.col))
pdf(paste0(sr.out,"Monkey_CS8_umap_tsne_plot.pdf"), height = 4, width = 12.5)
p1 <- DimPlot(mk.cs8, reduction = "umap", group.by = "CellType", cols = mk.col)
p2 <- DimPlot(mk.cs8, reduction = "tsne", group.by = "CellType", cols = mk.col)
p1 + p2
dev.off(); rm(p1, p2)
pdf(paste0(sr.out,"Monkey_CS8_marker_genes_dotplot.pdf"), height = 6.5, width = 16)
DotPlot(mk.cs8, features = intersect(unique(as.vector(unlist(mk.ge.mk))), rownames(mk.cs8)), group.by = "CellType",
        col.min = -1.5, col.max = 2) +
  scale_color_gradient2(high = "#cd120c", mid = "#ffffff", low = "#04579c", limits = c(-1.5, 2)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        panel.grid = element_blank(),
        panel.background= element_rect(color = "black", fill = NA, linewidth = 1.25),
        axis.text.y = element_text(size = 12),
        legend.box = "horizontal")
dev.off()
mk.cs8$DataSet <- "MK"



### ========================================
### 3rd part: Annotation of GAS-D5 by Seurat ----
### ========================================

### >>> 1. taking human CS7 data as reference
i <- 25
for (i in c(25)) {
  # - gene to be plotted
  pd.ge <- c("TFAP2A","DLX5","GATA3","CDH1","SOX2","POU5F1","OTX2","TBXT","SP5",
             "FST","MESP1","TBX6","PDGFRA","LHX1","IRX3","GATA6","MYL7","MSX1",
             "FOXF1","HAND1","BMP4","LUM","NID2","ANXA1","MAB21L2","SOX17","FOXA2",
             "CST1","MEF2C","PECAM1","CDH5","TEK")
  # - mapping and visualization
  sr.pdt <- list()
  anchors <- FindTransferAnchors(reference = hs.cs7, query = sr.gas$d5, dims = 1:i, reference.reduction = "pca")
  # cell type
  sr.pdt$gasd5.to.hs <- MapQuery(anchorset = anchors, reference = hs.cs7, query = sr.gas$d5,
                                 refdata = list(celltype = "CellType"), reference.reduction = "pca", reduction.model = "umap")
  sr.pdt$gasd5.to.hs@meta.data <- sr.pdt$gasd5.to.hs@meta.data %>%
    mutate(CellType=factor(predicted.celltype,
                           levels = names(pd.col.fix)))
  # umap and tsne plot
  p1 <- DimPlot(hs.cs7, reduction = "umap", group.by = "CellType", label = F, pt.size = 1) +
    ggtitle("Human embryos") + scale_color_manual(values = pd.col.fix) + theme(aspect.ratio=4/7) + ylim(c(-6,12))
  p2 <- DimPlot(sr.pdt$gasd5.to.hs, reduction = "ref.umap", group.by = "CellType", label = F, pt.size = 1) +
    ggtitle("GAS") + scale_color_manual(values = pd.col.fix) + theme(aspect.ratio=4/7) + ylim(c(-6,12))
  pdf(paste0(sr.out, "annotation_by_humanCS7_final/dim", i, "/GAS_D5_and_human_CS7_embryo_visualized_separately_umap.pdf"),
      height = 3.5, width = 10)
  print(p1+p2)
  dev.off()
  # GAS_D5 umap plot (GFP pos and neg separately)
  p1 <- DimPlot(sr.pdt$gasd5.to.hs, reduction = "ref.umap", group.by = "CellType", label = F, pt.size = 2,
                cells = colnames(subset(sr.pdt$gasd5.to.hs,Type=="GFP_pos"))) +
    ggtitle("GFP_pos") + scale_color_manual(values = pd.col.fix) + theme(aspect.ratio=4/7)
  p2 <- DimPlot(sr.pdt$gasd5.to.hs, reduction = "ref.umap", group.by = "CellType", label = F, pt.size = 2,
                cells = colnames(subset(sr.pdt$gasd5.to.hs,Type=="GFP_neg"))) +
    ggtitle("GFP_neg") + scale_color_manual(values = pd.col.fix) + theme(aspect.ratio=4/7)
  pdf(paste0(sr.out, "annotation_by_humanCS7_final/dim", i, "/GAS_D5_GFP_pos_and_neg_visualized_separately_umap.pdf"),
      height = 3.5, width = 10)
  print(p1+p2)
  dev.off()
  pdf(paste0(sr.out, "annotation_by_humanCS7_final/dim", i, "/GAS_D5_group_by_GFP_pos_and_neg_umap.pdf"),
      height = 3.5, width = 5)
  DimPlot(sr.pdt$gasd5.to.hs, reduction = "ref.umap", group.by = "Type", label = F, pt.size = 2, order = c("GFP_pos","GFP_neg")) +
    ggtitle("GAS") + scale_color_manual(values = c("#C9DDFB","#8F5090")) + theme(aspect.ratio=4/7)
  dev.off()
  # featureplot hs CS7
  dir.create(paste0(sr.out, "annotation_by_humanCS7_final/dim", i, "/featureplot/"))
  for (g in seq(1, length(pd.ge))) {
    pdf(paste0(sr.out, "annotation_by_humanCS7_final/dim", i, "/featureplot/humanCS7_expression_of_", pd.ge[g], "_featureplot.pdf"),
        height = 3.5, width = 4)
    print(FeaturePlot(hs.cs7, features = pd.ge[g], cols = c("#ebebeb", "#7d3c98"), pt.size = 1.5) + theme(aspect.ratio=4/7) + ylim(c(-6,12)))
    dev.off()
  }; rm(g)
  # featureplot GAS-D5
  for (g in seq(1, length(pd.ge))) {
    pdf(paste0(sr.out, "annotation_by_humanCS7_final/dim", i, "/featureplot/GAS_D5_expression_of_", pd.ge[g], "_featureplot.pdf"),
        height = 3.5, width = 4)
    print(FeaturePlot(sr.pdt$gasd5.to.hs, features = pd.ge[g], cols = c("#ebebeb", "#7d3c98"), reduction = "ref.umap", pt.size = 1.5) +
            theme(aspect.ratio=4/7) + ylim(c(-6,12)))
    dev.off()
  }; rm(g)
  # - integrated human CS7 and GAS-D5 visualization
  sr.pdt$integrated.hs <- merge(sr.pdt$gasd5.to.hs, hs.cs7)
  sr.pdt$integrated.hs@meta.data <- sr.pdt$integrated.hs@meta.data[,-c(13:19)]
  sr.pdt$integrated.hs@meta.data <- sr.pdt$integrated.hs@meta.data %>%
    mutate(CellType.sub=paste(DataSet, CellType, sep = "-"))
  tmp <- rbind(Embeddings(sr.pdt$gasd5.to.hs, reduction = "ref.umap"),
               Embeddings(hs.cs7, reduction = "umap"))
  sr.pdt$integrated.hs[['umap']] <- CreateDimReducObject(embeddings = tmp, key = "UMAP_")
  sr.pdt$integrated.hs <- sr.pdt$integrated.hs %>% NormalizeData() %>% ScaleData(rownames(sr.pdt$integrated.hs))
  sr.pdt$integrated.hs@meta.data <- sr.pdt$integrated.hs@meta.data %>% mutate(DataSet=factor(DataSet, levels = c("NE","GAS")))
  # umap + tsne plot
  pdf(paste0(sr.out, "annotation_by_humanCS7_final/dim", i, "/GAS_D5_and_humanCS7_integrated_umap_tsne.pdf"),
       height = 4, width = 10)
  p1 <- DimPlot(sr.pdt$integrated.hs, reduction = "umap", group.by = c("CellType"), cols = pd.col.fix, pt.size = 1.5) +
    theme(aspect.ratio=4/7)
  p2 <- DimPlot(sr.pdt$integrated.hs, reduction = "umap", group.by = c("DataSet"), cols = c("#C9DDFB","#8F5090"),
                pt.size = 1.5, order = c("GAS","NE")) + theme(aspect.ratio=4/7)
  print(p1+p2)
  dev.off()
  sr.pdt$integrated.hs@meta.data <- sr.pdt$integrated.hs@meta.data %>%
    mutate(Type=case_when(DataSet=="NE" ~ "GFP_neg", TRUE ~ Type))
  pdf(paste0(sr.out, "annotation_by_humanCS7_final/dim", i, "/GAS_D5_and_humanCS7_integrated_umap_group_by_GFP_pos_and_neg.pdf"),
      height = 3.5, width = 5)
  DimPlot(sr.pdt$integrated.hs, reduction = "umap", group.by = "Type", label = F, pt.size = 1.5, order = c("GFP_pos","GFP_neg")) +
    ggtitle("NE and GAS integrated") + scale_color_manual(values = c("#C9DDFB","#8F5090")) + theme(aspect.ratio=4/7)
  dev.off()
  # featureplot
  for (g in seq(1, length(pd.ge))) {
    pdf(paste0(sr.out, "annotation_by_humanCS7_final/dim", i, "/featureplot/GAS_D5_and_humanCS7_integrated_expression_of_",
               pd.ge[g], "_featureplot.pdf"),
        height = 3.5, width = 4)
    print(FeaturePlot(sr.pdt$integrated.hs, features = pd.ge[g], cols = c("#ebebeb", "#7d3c98"), pt.size = 1.5) +
            theme(aspect.ratio=4/7))
    dev.off()
  }; rm(g)
  # cell type proportion
  p <- sr.pdt$integrated.hs@meta.data[,c("CellType","DataSet")]
  p1 <- table(subset(p, DataSet=="GAS")$CellType) %>% as.data.frame() %>% mutate(Group="GAS")
  p2 <- table(subset(p, DataSet=="NE")$CellType) %>% as.data.frame() %>% mutate(Group="NE")
  p <- rbind(p1,p2)
  colnames(p)[1] <- "CellType"
  ct <- setdiff(p[p$Group=="NE",]$CellType, p[p$Group=="GAS",]$CellType)
  p <- rbind(p, data.frame(CellType=ct,
                           Freq=rep(0, length(ct)),
                           Group=rep("GAS", length(ct))))
  # barplot
  pdf(paste0(sr.out, "annotation_by_humanCS7_final/dim", i, "/cell_type_proportion_bar_plot.pdf"), height = 3, width = 3)
  print(p %>% mutate(CellType=factor(CellType, levels = c("ECT","EPI","PS","NM","AxM","EM","AM","ExM","End","HEP","Ery")),
                     Group=factor(Group, levels = c("NE", "GAS"))) %>%
          ggplot(aes(fill = CellType, y = Freq, x = Group)) +
          geom_bar(position = "fill", stat = "identity") +
          scale_fill_manual(values = pd.col.fix) + theme_bw() + ylab("Ratio") +
          theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 12),
                axis.text.x = element_text(size = 8, hjust = 1, angle = 45), axis.text.y = element_text(size = 12)))
  dev.off()
  table(sr.pdt$integrated.hs@meta.data[,c("CellType","DataSet")])
  # dotplot
  p %>% group_by(Group) %>% mutate(total = sum(Freq)) %>%
    mutate(Ratio = (Freq/total)*100) %>%
    as.data.frame() -> p
  p <- p[,c("Group", "CellType", "Ratio")]
  p <- spread(p, Group, Ratio)
  pdf(paste0(sr.out, "annotation_by_humanCS7_final/dim", i, "/cell_type_proportion_dot_plot.pdf"), height = 3.75, width = 4.5)
  print(p %>% mutate(CellType=factor(CellType, levels = names(pd.col.fix))) %>%
          ggplot(aes(x=`GAS`, y=`NE`, group=CellType)) +
          geom_point(aes(color = CellType), size = 2.5) +
          scale_color_manual(values = pd.col.fix) +
          geom_abline(intercept = 0) +
          xlim(0,40) + ylim(0,40) + xlab("GAS (Cell Type Proportion)") + ylab("NE (Cell Type Proportion)") +
          theme_bw() + theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
                             axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)))
  dev.off(); rm(p, p1, p2, ct)
  # cell type and GFP type
  p <- table(sr.pdt$gasd5.to.hs@meta.data[,c("CellType","Type")]) %>% as.data.frame() %>% filter(!CellType%in%c("Ery","AxM"))
  pdf(paste0(sr.out, "annotation_by_humanCS7_final/dim", i, "/cell_type_and_GFP_positive_proportion_bar_plot.pdf"),
      height = 2.25, width = 3.25)
  print(p %>% mutate(CellType=factor(CellType, levels = c("ECT","EPI","PS","NM","EM","AM","ExM","End","HEP"))) %>%
          ggplot(aes(fill = Type, y = Freq, x = CellType)) +
          geom_bar(position = "fill", stat = "identity") +
          scale_fill_manual(values = c("#D2C8E6","#1eaf90")) + theme_bw() + ylab("Ratio") +
          theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 12),
                axis.text.x = element_text(size = 8, hjust = 1, angle = 45), axis.text.y = element_text(size = 12)))
  dev.off()

  # saving seurat object
  saveRDS(sr.pdt, paste0("/home/cmq/bioinfo/project-tmp/wb/analysis/rds/annotation_by_humanCS7_final/sr.pdt", ".dim", i, ".rds"))
  # dotplot
  predicted.ct <- c(paste0("NE-", unique(as.vector(sr.pdt$gasd5.to.hs$CellType))),
                    paste0("GAS-", unique(as.vector(sr.pdt$gasd5.to.hs$CellType))))
  sr.pdt$integrated.hs <- subset(sr.pdt$integrated.hs, CellType.sub%in%predicted.ct) %>%
    NormalizeData() %>% ScaleData(rownames(sr.pdt$integrated.hs))
  sr.pdt$integrated.hs@meta.data <- sr.pdt$integrated.hs@meta.data %>%
    mutate(CellType.sub=factor(CellType.sub,
                               levels = c("GAS-HEP", "NE-HEP", "GAS-End", "NE-End", "GAS-ExM", "NE-ExM",
                                          "GAS-AM", "NE-AM", "GAS-EM", "NE-EM", "GAS-NM", "NE-NM",
                                          "GAS-PS", "NE-PS", "GAS-EPI", "NE-EPI", "GAS-ECT", "NE-ECT")))
  saveRDS(sr.pdt$integrated.hs, paste0("/home/cmq/bioinfo/project-tmp/wb/analysis/rds/annotation_by_humanCS7_final/sr.pdt.filtered", ".dim", i, ".rds"))
  # marker genes dotplot
  pdf(paste0(sr.out, "annotation_by_humanCS7_final/dim", i, "/GAS_D5_and_humanCS7_integrated_expression_of_humanCS7_marker_genes_dotplot1.pdf"),
       height = 5, width = 12)
  print(DotPlot(sr.pdt$integrated.hs, features = pd.ge, group.by = "CellType.sub",
                col.min = -2, col.max = 2.5) +
          scale_color_gradient2(high = "#830923", mid = "#C3F1F9", low = "#16518F", limits = c(-2, 2.5)) +
          theme_bw() +
          theme(axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                panel.grid = element_blank(),
                panel.background= element_rect(color = "black", fill = NA, linewidth = 1.25),
                axis.text.y = element_text(size = 12),
                legend.box = "horizontal"))
  dev.off()
  sr.pdt$integrated.hs@meta.data <- sr.pdt$integrated.hs@meta.data %>%
    mutate(CellType.sub=factor(CellType.sub,
                               levels = c("GAS-HEP", "GAS-End", "GAS-ExM", "GAS-AM", "GAS-EM", "GAS-NM", "GAS-PS", "GAS-EPI", "GAS-ECT",
                                          "NE-HEP", "NE-End", "NE-ExM", "NE-AM", "NE-EM", "NE-NM", "NE-PS", "NE-EPI", "NE-ECT" )))
  pdf(paste0(sr.out, "annotation_by_humanCS7_final/dim", i, "/GAS_D5_and_humanCS7_integrated_expression_of_humanCS7_marker_genes_dotplot2.pdf"),
      height = 5, width = 12)
  print(DotPlot(sr.pdt$integrated.hs, features = pd.ge, group.by = "CellType.sub",
                col.min = -2, col.max = 2.5) +
          scale_color_gradient2(high = "#830923", mid = "#C3F1F9", low = "#16518F", limits = c(-2, 2.5)) +
          theme_bw() +
          theme(axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                panel.grid = element_blank(),
                panel.background= element_rect(color = "black", fill = NA, linewidth = 1.25),
                axis.text.y = element_text(size = 12),
                legend.box = "horizontal"))
  dev.off()
  # - delete useless variables
  rm(p1, p2, p3, anchors, predicted.ct, tmp)
  # - cell type correlation analysis
  expr1 = AverageExpression(sr.pdt$gasd5.to.hs, assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
  p <- colnames(expr1)
  colnames(expr1) = paste("GAS", gsub("RNA.", "", colnames(expr1)), sep = "-")
  expr2 = AverageExpression(hs.cs7, assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
  expr2 <- expr2[,p]
  colnames(expr2) = paste("NE", gsub("RNA.", "", colnames(expr2)), sep = "-")
  deg1 = VariableFeatures(sr.pdt$gasd5.to.hs)
  deg2 = VariableFeatures(hs.cs7)
  spe.score.seurat <- CorrComparePlot(ExpressionTableSpecies1 = expr1, DEgenesSpecies1 = deg1,
                                      ExpressionTableSpecies2 = expr2, DEgenesSpecies2 = deg2,
                                      Species1 = "human", Species2 = "human",
                                      filename = file.path(paste0(sr.out, "annotation_by_humanCS7_final/dim", i),
                                                           paste0("Specificity_score_correlation_plot_of_All_cell_types_in_GAS_D5_and_humanCS7")))
  rm(expr1, expr2, deg1, deg2, spe.score.seurat)
}; rm(i)
# - reverse the umap plot
sr.pdt <- readRDS(paste0("/home/cmq/bioinfo/project-tmp/wb/analysis/rds/annotation_by_humanCS7_final/sr.pdt", ".dim", i, ".rds"))
saveRDS(sr.pdt$integrated.hs, "/home/cmq/bioinfo/project-tmp/wb/analysis/rds/sr.pdt.integrated.hs.rds")
# human CS7
tmp <- Embeddings(hs.cs7, reduction = "umap")
tmp[,2] <- -tmp[,2]
hs.cs7[['umap_rv']] <- CreateDimReducObject(tmp)
# GAS-D5
tmp <- Embeddings(sr.pdt$gasd5.to.hs, reduction = "ref.umap")
tmp[,2] <- -tmp[,2]
sr.pdt$gasd5.to.hs[['umap_rv']] <- CreateDimReducObject(tmp)
# Integrated
tmp <- Embeddings(sr.pdt$integrated.hs, reduction = "umap")
tmp[,2] <- -tmp[,2]
sr.pdt$integrated.hs[['umap_rv']] <- CreateDimReducObject(tmp)
# umap seperatly
p1 <- DimPlot(hs.cs7, reduction = "umap_rv", group.by = "CellType", label = F, pt.size = 0.5) +
  ggtitle("Natural embryos") + scale_color_manual(values = pd.col.fix)
p2 <- DimPlot(sr.pdt$gasd5.to.hs, reduction = "umap_rv", group.by = "CellType", label = F, pt.size = 0.5) +
  ggtitle("GAS") + scale_color_manual(values = pd.col.fix)
pdf(paste0(sr.out, "annotation_by_humanCS7_final/dim", i, "/GAS_D5_and_human_CS7_embryo_visualized_separately_umap_reverse.pdf"),
    height = 4, width = 9)
print(p1+p2)
dev.off()
# featureplot hs CS7
dir.create(paste0(sr.out, "annotation_by_humanCS7_final/dim25/featureplot_reverse/"))
sr.out <- "/home/yhw/bioinfo/project-wb/hs_gastruloids/analysis/graphs/seurat_output/"
pd.ge <- c("AMOT","AMOTL1","AMOTL2","CIT","DBX1","DBX2","DLG5","EBF1","EGR2",
           "ELAVL3","EN1","EN2","FAT4","FOXG1","GBX2","HAND1","HESX1","HOXA3",
           "IRX3","LATS1","LATS2","LHX5","LIMD1","MAFB","MAP2K3","MAPK14","MARK3",
           "MEIS1","NEK8","NEUROD1","NEUROG1","NF2","NFATC4",
           "OLIG2","OTX2","PAX3","PAX5","PAX6","PAX7","PAX8","PJA2","RREB1","SAV1",
           "SHANK2","SMARCA5","SOX1","SOX11","SOX3","STK3","STK4","TEAD1","TEAD2",
           "TEAD3","TEAD4","THBS1","TIAL1","TUBB3","WRNIP1","WTIP","WWC1","YAP1","ZNF616")
hox.gene <- read.table("/home/yhw/refgenome/embl/rs107/homo_sapiens/Homo_sapiens.GRCh38.107.chr_CellRanger_gene_pro1k.bed")
hox.gene <- hox.gene[grep("HOX", hox.gene$V4),]
hox.gene <- subset(hox.gene, V4 %in% intersect(rownames(hs.cs7), intersect(rownames(sr.pdt$integrated.hs), rownames(sr.pdt$gasd5.to.hs))))
pd.ge <- hox.gene$V4
for (g in seq(1, length(pd.ge))) {
  pdf(paste0(sr.out, "annotation_by_humanCS7_final/dim25/featureplot_reverse/humanCS7_expression_of_", pd.ge[g], "_featureplot.pdf"),
      height = 3.5, width = 4)
  print(FeaturePlot(hs.cs7, features = pd.ge[g], cols = c("#ebebeb", "#7d3c98"), pt.size = 1.35, reduction = "umap"))
  dev.off()
}; rm(g)
sr.out <- "/home/yhw/bioinfo/project-wb/hs_gastruloids/analysis/graphs/seurat_output/annotation_by_humanCS7_final/dim25/CS7_marker_gene"
dir.create(sr.out, recursive = T)
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("CS7_marker_gene_expression_of_",
                                  marker.gene[[type]][g], "_featureplot.pdf")),
        height = 3.5, width = 4)
    print(FeaturePlot(hs.cs7, features = marker.gene[[type]][g],
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 1.35, reduction = "umap") +
            xlim(c(-10, 12.5)) + ylim(c(-8, 12.5)) +
            theme(aspect.ratio=4/7))
    dev.off()
  }
}
# featureplot GAS-D5
for (g in seq(1, length(pd.ge))) {
  pdf(paste0(sr.out, "annotation_by_humanCS7_final/dim25/featureplot_reverse/GAS_D5_expression_of_", pd.ge[g], "_featureplot.pdf"),
      height = 3.5, width = 4)
  print(FeaturePlot(sr.pdt$gasd5.to.hs, features = pd.ge[g], cols = c("#ebebeb", "#7d3c98"), reduction = "ref.umap", pt.size = 1.35))
  dev.off()
}; rm(g)
sr.out <- "/home/yhw/bioinfo/project-wb/hs_gastruloids/analysis/graphs/seurat_output/annotation_by_humanCS7_final/dim25/GASD5_marker_gene"
dir.create(sr.out, recursive = T)
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("GAS_D5_marker_gene_expression_of_",
                                  marker.gene[[type]][g], "_featureplot.pdf")),
        height = 3.5, width = 4)
    print(FeaturePlot(sr.pdt$gasd5.to.hs, features = marker.gene[[type]][g],
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 1.35, reduction = "ref.umap") +
            xlim(c(-10, 12.5)) + ylim(c(-8, 12.5)) +
            theme(aspect.ratio=4/7))
    dev.off()
  }
}
# integrated umap
pdf(paste0(sr.out, "annotation_by_humanCS7_final/dim25/GAS_D5_and_humanCS7_integrated_umap_tsne_reverse.pdf"),
    height = 4, width = 14)
p1 <- DimPlot(sr.pdt$integrated.hs, reduction = "umap_rv", group.by = c("CellType"), cols = pd.col.fix, pt.size = 0.5)
p2 <- DimPlot(sr.pdt$integrated.hs, reduction = "umap_rv", group.by = c("DataSet"), cols = c("#1eaf90","#fd79a8"), pt.size = 0.5)
p3 <- DimPlot(sr.pdt$integrated.hs, reduction = "umap_rv", group.by = c("DataSet"), cols = c("#fd79a8", "#1eaf90"), pt.size = 0.5)
print(p1+p2+p3)
dev.off()
# integrated featureplot
for (g in seq(1, length(pd.ge))) {
  pdf(paste0(sr.out, "annotation_by_humanCS7_final/dim25/featureplot_reverse/GAS_D5_and_humanCS7_integrated_expression_of_",
             pd.ge[g], "_featureplot.pdf"),
      height = 3.5, width = 4)
  print(FeaturePlot(sr.pdt$integrated.hs, features = pd.ge[g], cols = c("#ebebeb", "#7d3c98"), pt.size = 1.35, reduction = "umap"))
  dev.off()
}; rm(g)
# plot expression
sr.out <- "/home/yhw/bioinfo/project-wb/hs_gastruloids/analysis/graphs/seurat_output/annotation_by_humanCS7_final/dim25/integrated_marker_gene"
dir.create(sr.out, recursive = T)
marker.gene <- list(FMH.marker  = c("LHX5","ZNF616","GBX2","RREB1","TEAD4","EN2","HAND1",
                                    "HESX1","FOXG1","NFATC4","PAX7","DBX1","PAX5","SMARCA5",
                                    "PAX6","IRX3","EBF1","HOXA3","MEIS1","WRNIP1","MAFB"),
                    FMH = c("PAX6","OTX2","EN1","PAX8","GBX2","EGR2",
                            "PAX3","DBX2","DBX1","IRX3","OLIG2","NKX2.2"),
                    Neural = c("SOX1","SOX3","PAX6","TUBB3","ELAVL3",
                               "OLIG2","NEUROG1","NEUROD1","NIKX2.2"),
                    NMP = c("LIMD1","MARK3","CIT","MAPK14","FAT4","LATS2",
                            "STK4","WWC1","SHANK2","SOX11","TEAD1","PJA2",
                            "TEAD2","VGLL5","SAV1","MAP2K3","TIAL1","NF2",
                            "STK3","LATS1","THBS1","TEAD3","AMOT","YAP1",
                            "AMOTL1","AMOTL2","NEK8","WTIP","DLG5","TEAD4"),
                    Amnion = c("GABRP","TGFB1","VTCN1","S100A10","GADD45G","ACTC1",
                               "IGFBP3","MCM5","CLDN10","TUBG1","UCP2","BIN3",
                               "GGA2","FGFR1","GALM"),
                    Anterior.marker = c("LHX5","POU5F1","SOX2","POU3F1","EGR2","MAFB","OTX2",
                                        "HESX1","ID1","NR2F1","SOX3","IRX2","SIX1","SIX3",
                                        "SIX6","HOXA2","POU3F2","HOXB-AS1","CRABP1","TUBB3",
                                        "PBX1","BTG1","NES","MEIS2","SFRP1","HES5","PAX6","CDH1",
                                        "FOXC1","HOXA1","HOXA3","HOXB1","HOXB2","HOXB3","HOXD1",
                                        "HOXD3","IRX1","MEIS2","ISL1","FGFR1"),
                    Posterior.marker = c("CDX4","HOXA-AS3","CDX2","PRICKLE1","DUSP6","ID3",
                                         "SFRP2","TBXT","GSC","HOXA5","HOXA7","HOXA9","HOXA10",
                                         "HOXA11","HOXA13","HOXB4","HOXB5","HOXB6","HOXB7",
                                         "HOXB8","HOXB9","HOXC4","HOXC6","HOXC8","VIM","FGF8",
                                         "FGF17","FOXA2","SHH","EOMES","DLL1","CDH2","BMP4",
                                         "WNT3A","MESP1","FOXB1","TBX6"),
                    EPI.marker = c("CDH1","POU5F1","SOX2","NANOG","LIN28",
                                   "MYC","KLF4","TDGF1","GDF3"),
                    Organizer = c("FST","OTX2","GSC","EOMES","FOXH1","CER1","MIXL1",
                                  "DKK1","NODAL","FOXA2","LEFTY1","LEFTY2","TBXT",
                                  "SOX17","SP5","NOGGIN","CHRD"),
                    Embryonic.marker = c("TBXT6","MESP1","PDGFRA","LHX1","IRX3","GATA6",
                                         "MYL7","MSX1","MSX2","FOXF1","HAND1","HAND2",
                                         "BMP4,","MAB21L2","CST1","MEF2C","PECAM1","CDH5",
                                         "TEK","TFAP2A","TFAP2C","GATA3","DLX5","LHX5",
                                         "PAX6","SOX3","PAX3","PAX7","HOXB9"),
                    Extraembryonic = c("LUM","NID2","ANXA1","POSTN","GATA2","KRT7","VIM",
                                       "GATA4","DCN","TP63","KRT18","NR2F2","ISL1","HEY1",
                                       "CDH10","CTSV","ARID5B","PLAGL1","CREB3L1","HOXA10",
                                       "HOXA11","HOXA9","HOXA13","WNT6","GABRP","HEY1","BST2"))
for (type in names(marker.gene)) {
  #marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(sr.pdt$integrated.hs))
  #marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(sr.pdt$gasd5.to.hs))
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(hs.cs7))
}
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("GAS_D5_and_humanCS7_integrated_marker_gene_expression_of_",
                                  marker.gene[[type]][g], "_featureplot.pdf")),
        height = 3.5, width = 4)
    print(FeaturePlot(sr.pdt$integrated.hs, features = marker.gene[[type]][g],
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 1.35, reduction = "umap") +
            xlim(c(-10, 12.5)) + ylim(c(-8, 12.5)) +
            theme(aspect.ratio=4/7))
    dev.off()
  }; rm(g)
}
# removing useless variables
rm(tmp, p, p1, p2, p3, i)

# - GO analysis of cell in GAS D5
# find marker genes
Idents(sr.pdt$gasd5.to.hs) <- sr.pdt$gasd5.to.hs$CellType
Idents(hs.cs7) <- hs.cs7$CellType
gasd5.to.hs.mk <- FindAllMarkers(sr.pdt$gasd5.to.hs, assay = "RNA", only.pos = T, logfc.threshold = log2(1.5))
cs7.mk <- FindAllMarkers(hs.cs7, assay = "RNA", only.pos = T, logfc.threshold = log2(1.5))
# GO
sr.out <- "/home/yhw/bioinfo/project-wb/hs_gastruloids/analysis/graphs/seurat_output/annotation_by_humanCS7_final/dim25/GO_markers"
dir.create(sr.out, recursive = T)
gas.markers.go <- list()
for (type in unique(as.character(gasd5.to.hs.mk$cluster))) {
  pd.gene <- gasd5.to.hs.mk %>%
    mutate(cluster = as.character(cluster)) %>%
    filter(cluster == type)
  gas.markers.go[[type]] <- Pipe.GO(species = "human", genelist = pd.gene$gene, basename = type,
                                    genetype = "SYMBOL", res.out = file.path(sr.out, "GASD5"))
}
cs7.markers.go <- list()
for (type in unique(as.character(cs7.mk$cluster))) {
  pd.gene <- cs7.mk %>%
    mutate(cluster = as.character(cluster)) %>%
    filter(cluster == type)
  cs7.markers.go[[type]] <- Pipe.GO(species = "human", genelist = pd.gene$gene, basename = type,
                                    genetype = "SYMBOL", res.out = file.path(sr.out, "CS7"))
}
go.terms <- list(ECT = c("GO:0016331","GO:0072009","GO:0001892","GO:0043588",
                         "GO:0060706","GO:0001890","GO:0060706","GO:0060669",
                         "GO:0001655","GO:0001822","GO:0048562","GO:0045669",
                         "GO:0001657","GO:0072164","GO:0072163","GO:0001823",
                         "GO:0003007","GO:0071772","GO:0071773","GO:0001838"),
                 EPI = c("GO:0060795","GO:0001708","GO:0035019","GO:0060026",
                         "GO:0042249","GO:0090175","GO:0060071","GO:0045664",
                         "GO:0007272","GO:0045664","GO:0090178","GO:0090175",
                         "GO:1904338","GO:0090179","GO:0008630","GO:1902230",
                         "GO:0019827","GO:2000035"),
                 PS = c("GO:0060795","GO:0001704","GO:0007492","GO:0007369",
                        "GO:0042552","GO:0035019","GO:0001708"),
                 NM = c("GO:0009952","GO:0003002","GO:0007389","GO:0007369",
                        "GO:0001704","GO:0007498","GO:0060485","GO:0003007",
                        "GO:0001756","GO:0035282","GO:0061053","GO:0061371",
                        "GO:0001707","GO:0048332","GO:0061311","GO:0035050",
                        "GO:0001837","GO:0061053"),
                 EM = c("GO:0001704","GO:0007369","GO:0007498","GO:0009952",
                        "GO:0061311","GO:0007389","GO:0035282","GO:0003002",
                        "GO:0001756","GO:0061053","GO:0009948","GO:0009948",
                        "GO:0061371","GO:0032525","GO:0048732","GO:0001947",
                        "GO:0003231","GO:0003206","GO:0007368","GO:0001822"),
                 AM = c("GO:0048568","GO:0003007","GO:0060485","GO:0030198",
                        "GO:0043062","GO:0045229","GO:0048562","GO:0048562",
                        "GO:0014706","GO:0048762","GO:0048705","GO:0001704",
                        "GO:0007369","GO:0003206","GO:0060537","GO:0048738",
                        "GO:0003205","GO:0007389","GO:0051216","GO:0001503",
                        "GO:0042476","GO:0022612","GO:0060411"),
                 ExM = c("GO:0030198","GO:0043062","GO:0045229","GO:0001667",
                         "GO:0061448","GO:0003206","GO:0060485","GO:0003007",
                         "GO:0045765","GO:0050673","GO:1904018","GO:0045766",
                         "GO:0050678","GO:1901342","GO:0002040","GO:0032963",
                         "GO:0070482","GO:0001666","GO:0036293","GO:0022604"),
                 End = c("GO:0030198","GO:0007492","GO:0001706","GO:0043062",
                         "GO:0045229","GO:0048608","GO:0061458","GO:0060485",
                         "GO:0060562","GO:0001655","GO:0048762","GO:0014706",
                         "GO:0072001","GO:0048568","GO:0007440","GO:0009948"),
                 HEP = c("GO:0043542","GO:0010631","GO:0090132","GO:0003158",
                         "GO:0010594","GO:0090050","GO:0002040","GO:0045766",
                         "GO:0010595","GO:1901342","GO:0001945","GO:0001938",
                         "GO:0001935","GO:0043534","GO:0010632","GO:0007599",
                         "GO:0007596","GO:0050817","GO:1903706","GO:0060216"))
sr.out <- "/home/yhw/bioinfo/project-wb/hs_gastruloids/analysis/graphs/seurat_output/annotation_by_humanCS7_final/dim25/GO_markers"
for (i in names(go.terms)) {
  pdf(file.path(sr.out, paste0("GAS_GO_", i, ".pdf")), height = 10, width = 10)
  print(subset(gas.markers.go[[i]]$GO.BP, ID %in% go.terms[[i]]) %>%
          mutate(Description = as.factor(Description), `-Log10(pvalue)` = -log10(pvalue)) %>%
          mutate(Description = fct_reorder(Description, `-Log10(pvalue)`)) %>%
          ggplot(aes(x = `-Log10(pvalue)`, y = Description)) +
          geom_bar(stat = "identity", fill = "#ffffff", linewidth = 1, color = pd.col.fix[i]) +
          geom_text(aes(x = 0.05, label = Description), hjust = 0, size = 3.25, stat = "identity") +
          xlab("-log10(pvalue)") +
          theme_bw() +
          ggtitle(i) +
          theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank()
          ))
  dev.off()
  pdf(file.path(sr.out, paste0("CS7_GO_", i, ".pdf")), height = 10, width = 10)
  print(subset(cs7.markers.go[[i]]$GO.BP, ID %in% go.terms[[i]]) %>%
          mutate(Description = as.factor(Description), `-Log10(pvalue)` = -log10(pvalue)) %>%
          mutate(Description = fct_reorder(Description, `-Log10(pvalue)`)) %>%
          ggplot(aes(x = `-Log10(pvalue)`, y = Description)) +
          geom_bar(stat = "identity", fill = "#ffffff", linewidth = 1, color = pd.col.fix[i]) +
          geom_text(aes(x = 0.05, label = Description), hjust = 0, size = 3.25, stat = "identity") +
          xlab("-log10(pvalue)") +
          theme_bw() +
          ggtitle(i) +
          theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank()
          ))
  dev.off()
}


### >>> 2. taking monkey CS8 as reference
# - load human-monkey homologous gene list
rs109.hs.mk.homo <- read.table(gzfile("/home/cmq/document/ensembl/homologous/human_monkey/release109_GRCh38_macFas5_mart_export.txt.gz"),
                               header = T, stringsAsFactors = F, sep = "\t")
rs109.hs.mk.homo <- rs109.hs.mk.homo[,c(1,3,10,11,16)] %>% unique()
colnames(rs109.hs.mk.homo) <- c("Hs_Gene_ID", "Hs_Gene_SYMBOL", "MK_Gene_ID", "MK_Gene_SYMBOL", "Type")
# - filtering non-homologous genes
tmp <- rs109.hs.mk.homo %>% filter(MK_Gene_SYMBOL!="" & Type=="ortholog_one2one")
tmp <- tmp[tmp$Hs_Gene_SYMBOL==tmp$MK_Gene_SYMBOL,]
tmp <- tmp %>% filter(Hs_Gene_SYMBOL%in%rownames(sr.gas$d5) & MK_Gene_SYMBOL%in%rownames(mk.cs8))
sr.pdt.mk <- list()
sr.pdt.mk$mk.cs8 <- mk.cs8[rownames(mk.cs8)%in%c(tmp$MK_Gene_SYMBOL, intersect(intersect(setdiff(as.vector(unlist(mk.ge.mk)), tmp$Hs_Gene_SYMBOL),
                                                                                         rownames(sr.gas$d5)), rownames(mk.cs8))),]
sr.pdt.mk$gas.d5 <- sr.gas$d5[rownames(sr.gas$d5)%in%c(tmp$Hs_Gene_SYMBOL, intersect(intersect(setdiff(as.vector(unlist(mk.ge.mk)), tmp$Hs_Gene_SYMBOL),
                                                                                               rownames(sr.gas$d5)), rownames(mk.cs8))),]
rm(tmp)
# - running seurat mapping
i <- 25
sr.pdt.mk$gasd5.to.mk <- readRDS(paste0("/home/cmq/bioinfo/project-tmp/wb/analysis/rds/annotation_by_monkeyCS8/sr.pdt.mk.gasd5.to.mk",
                                        ".dim", i, ".rds"))
sr.pdt.mk$integrated.mk <- readRDS(paste0("/home/cmq/bioinfo/project-tmp/wb/analysis/rds/annotation_by_monkeyCS8/sr.pdt.mk.integrated.mk",
                                          ".dim", i, ".rds"))

for (i in c(15,20,30,35)) {
  dir.create(paste0(sr.out, "annotation_by_monkeyCS8/dim", i))
  pd.ge <- intersect(unique(as.vector(unlist(mk.ge.mk))), rownames(sr.pdt.mk$mk.cs8))
  # - mapping and visualization
  anchors <- FindTransferAnchors(reference = sr.pdt.mk$mk.cs8, query = sr.pdt.mk$gas.d5, dims = 1:i, reference.reduction = "pca")
  # cell type
  sr.pdt.mk$gasd5.to.mk <- MapQuery(anchorset = anchors, reference = sr.pdt.mk$mk.cs8, query = sr.pdt.mk$gas.d5,
                                    refdata = list(celltype = "CellType"), reference.reduction = "pca", reduction.model = "umap")
  sr.pdt.mk$gasd5.to.mk@meta.data <- sr.pdt.mk$gasd5.to.mk@meta.data %>%
    mutate(CellType=factor(predicted.celltype, levels = names(mk.col)))
  # umap and tsne plot
  p1 <- DimPlot(sr.pdt.mk$mk.cs8, reduction = "umap", group.by = "CellType", label = F, pt.size = 0.5) +
    ggtitle("Natural embryos") + scale_color_manual(values = mk.col)
  p2 <- DimPlot(sr.pdt.mk$gasd5.to.mk, reduction = "ref.umap", group.by = "CellType", label = F, pt.size = 0.5) +
    ggtitle("GAS") + scale_color_manual(values = mk.col)
  pdf(paste0(sr.out, "annotation_by_monkeyCS8/dim", i, "/GAS_D5_and_monkey_CS8_embryo_visualized_separately_umap.pdf"),
      height = 4, width = 12)
  print(p1+p2)
  dev.off()
  # umap and tsne plot named by cell type from human
  if (all(rownames(sr.pdt$gasd5.to.hs@meta.data)==rownames(sr.pdt.mk$gasd5.to.mk@meta.data))) {
    sr.pdt.mk$gasd5.to.mk@meta.data$CellType.hs <- sr.pdt$gasd5.to.hs@meta.data$CellType
  }
  p1 <- DimPlot(sr.pdt.mk$mk.cs8, reduction = "umap", group.by = "CellType", label = F, pt.size = 0.5) +
    ggtitle("Monkey embryos") + scale_color_manual(values = mk.col)
  p2 <- DimPlot(sr.pdt.mk$gasd5.to.mk, reduction = "ref.umap", group.by = "CellType.hs", label = F, pt.size = 0.5) +
    ggtitle("GAS") + scale_color_manual(values =pd.col.fix)
  pdf(paste0(sr.out, "annotation_by_monkeyCS8/dim", i, "/GAS_D5_and_monkey_CS8_embryo_visualized_separately_umap_named_by_CellType_from_HumanCS7.pdf"),
      height = 4, width = 12)
  print(p1+p2)
  dev.off()
  # featureplot mk CS8
  dir.create(paste0(sr.out, "annotation_by_monkeyCS8/dim", i, "/featureplot/"))
  for (g in seq(1, length(pd.ge))) {
    pdf(paste0(sr.out, "annotation_by_monkeyCS8/dim", i, "/featureplot/monkeyCS8_expression_of_", pd.ge[g], "_featureplot.pdf"),
        height = 3.5, width = 4)
    print(FeaturePlot(sr.pdt.mk$mk.cs8, features = pd.ge[g], cols = c("#ebebeb", "#7d3c98"), pt.size = 1.35))
    dev.off()
  }; rm(g)
  # featureplot GAS-D5
  for (g in seq(1, length(pd.ge))) {
    pdf(paste0(sr.out, "annotation_by_monkeyCS8/dim", i, "/featureplot/GAS_D5_expression_of_", pd.ge[g], "_featureplot.pdf"),
        height = 3.5, width = 4)
    print(FeaturePlot(sr.pdt.mk$gasd5.to.mk, features = pd.ge[g], cols = c("#ebebeb", "#7d3c98"), reduction = "ref.umap", pt.size = 1.35))
    dev.off()
  }; rm(g)
  # - integrated monkey CS8 and GAS-D5 visualization
  sr.pdt.mk$integrated.mk <- merge(sr.pdt.mk$mk.cs8, sr.pdt.mk$gasd5.to.mk)
  sr.pdt.mk$integrated.mk@meta.data <- sr.pdt.mk$integrated.mk@meta.data %>%
    mutate(CellType.sub=paste(DataSet, CellType, sep = "-"))
  tmp <- rbind(Embeddings(sr.pdt.mk$mk.cs8, reduction = "umap"),
               Embeddings(sr.pdt.mk$gasd5.to.mk, reduction = "ref.umap"))
  sr.pdt.mk$integrated.mk[['umap']] <- CreateDimReducObject(embeddings = tmp, key = "UMAP_")
  sr.pdt.mk$integrated.mk <- sr.pdt.mk$integrated.mk %>% NormalizeData() %>% ScaleData(rownames(sr.pdt.mk$integrated.mk))
  # umap + tsne plot
  pdf(paste0(sr.out, "annotation_by_monkeyCS8/dim", i, "/GAS_D5_and_monkeyCS8_integrated_umap_tsne.pdf"),
      height = 4, width = 15.5)
  p1 <- DimPlot(sr.pdt.mk$integrated.mk, reduction = "umap", group.by = c("CellType"), cols = mk.col, pt.size = 0.5)
  p2 <- DimPlot(sr.pdt.mk$integrated.mk, reduction = "umap", group.by = c("DataSet"), cols = c("#1eaf90","#fd79a8"), pt.size = 0.5)
  p3 <- DimPlot(sr.pdt.mk$integrated.mk, reduction = "umap", group.by = c("DataSet"), cols = c("#fd79a8", "#1eaf90"), pt.size = 0.5)
  print(p1+p2+p3)
  dev.off()
  # featureplot
  for (g in seq(1, length(pd.ge))) {
    pdf(paste0(sr.out, "annotation_by_monkeyCS8/dim", i, "/featureplot/GAS_D5_and_monkeyCS8_integrated_expression_of_", pd.ge[g], "_featureplot.pdf"),
        height = 3.5, width = 4)
    print(FeaturePlot(sr.pdt.mk$integrated.mk, features = pd.ge[g], cols = c("#ebebeb", "#7d3c98"), pt.size = 1.35))
    dev.off()
  }; rm(g)
  # cell type and GFP type
  p <- table(sr.pdt.mk$gasd5.to.mk@meta.data[,c("CellType","Type")]) %>% as.data.frame()
  pdf(paste0(sr.out, "annotation_by_monkeyCS8/dim", i, "/cell_type_and_GFP_positive_proportion_bar_plot.pdf"), height = 2.25, width = 6)
  print(p %>%
          ggplot(aes(fill = Type, y = Freq, x = CellType)) +
          geom_bar(position = "fill", stat = "identity") +
          scale_fill_manual(values = c("#FF7F50","#1eaf90")) + theme_bw() + ylab("Ratio") +
          theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 12),
                axis.text.x = element_text(size = 8, hjust = 1, angle = 45), axis.text.y = element_text(size = 12)))
  dev.off()
  # saving seurat object
  saveRDS(sr.pdt.mk$gasd5.to.mk, paste0("/home/cmq/bioinfo/project-tmp/wb/analysis/rds/annotation_by_monkeyCS8/sr.pdt.mk.gasd5.to.mk",
                                        ".dim", i, ".rds"))
  saveRDS(sr.pdt.mk$integrated.mk, paste0("/home/cmq/bioinfo/project-tmp/wb/analysis/rds/annotation_by_monkeyCS8/sr.pdt.mk.integrated.mk",
                                          ".dim", i, ".rds"))
  # dotplot
  predicted.ct <- c(paste0("MK-", unique(as.vector(sr.pdt.mk$gasd5.to.mk$CellType))),
                    paste0("GAS-", unique(as.vector(sr.pdt.mk$gasd5.to.mk$CellType))))
  sr.pdt.mk$integrated.mk <- subset(sr.pdt.mk$integrated.mk, CellType.sub%in%predicted.ct) %>%
    NormalizeData() %>% ScaleData(rownames(sr.pdt.mk$integrated.mk))
  all.ct <- c("MK-EPI", "GAS-EPI", "MK-PS", "GAS-PS", "MK-PGC", "GAS-PGC", "MK-APS", "GAS-APS", "MK-Node", "GAS-Node",
              "MK-DE", "GAS-DE", "MK-Gut", "GAS-Gut", "MK-Nas.Meso", "GAS-Nas.Meso", "MK-LP.Meso", "GAS-LP.Meso",
              "MK-Caud.Meso", "GAS-Caud.Meso",
              "MK-Pharyg.Meso", "GAS-Pharyg.Meso", "MK-Al", "GAS-Al", "MK-ExE.Meso", "GAS-ExE.Meso", "MK-ys.Meso1",
              "GAS-ys.Meso1", "MK-ys.Meso2", "GAS-ys.Meso2",
              "MK-Mes", "GAS-Mes", "MK-EC", "GAS-EC", "MK-BP", "GAS-BP", "MK-Mac", "GAS-Mac", "MK-Ery1", "GAS-Ery1",
              "MK-Ery2", "GAS-Ery2", "MK-ECT", "GAS-ECT", "MK-NC", "GAS-NC", "MK-SE1", "GAS-SE1", "MK-SE2", "GAS-SE2",
              "MK-AM", "GAS-AM", "MK-VE", "GAS-VE", "MK-ys.Endo1", "GAS-ys.Endo1", "MK-ys.Endo2", "GAS-ys.Endo2")
  sr.pdt.mk$integrated.mk@meta.data <- sr.pdt.mk$integrated.mk@meta.data %>%
    mutate(CellType.sub=factor(CellType.sub,
                               levels = all.ct[all.ct%in%predicted.ct]))
  saveRDS(sr.pdt.mk$integrated.mk, paste0("/home/cmq/bioinfo/project-tmp/wb/analysis/rds/annotation_by_monkeyCS8/sr.pdt.mk.integrated.filtered",
                                          ".dim", i, ".rds"))
  # marker genes dotplot
  pdf(paste0(sr.out, "annotation_by_monkeyCS8/dim", i, "/GAS_D5_and_monkeyCS8_integrated_expression_of_monkeyCS8_marker_genes_dotplot.pdf"),
      height = 7, width = 16)
  print(DotPlot(sr.pdt.mk$integrated.mk, features = pd.ge, group.by = "CellType.sub",
                col.min = -2, col.max = 2.5) +
          scale_color_gradient2(high = "#cd120c", mid = "#ffffff", low = "#04579c", limits = c(-2, 2.5)) +
          theme_bw() +
          theme(axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                panel.grid = element_blank(),
                panel.background= element_rect(color = "black", fill = NA, linewidth = 1.25),
                axis.text.y = element_text(size = 12),
                legend.box = "horizontal"))
  dev.off()
  # - delete useless variables
  rm(p1, p2, p3, anchors, predicted.ct, tmp)
  # - cell type correlation analysis
  expr1 = AverageExpression(sr.pdt.mk$gasd5.to.mk, assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
  p <- colnames(expr1)
  colnames(expr1) = paste("GAS", gsub("RNA.", "", colnames(expr1)), sep = "-")
  expr2 = AverageExpression(sr.pdt.mk$mk.cs8, assays = "RNA", group.by = "CellType", slot = "data") %>% as.data.frame()
  expr2 <- expr2[,p]
  colnames(expr2) = paste("MK", gsub("RNA.", "", colnames(expr2)), sep = "-")
  deg1 = VariableFeatures(sr.pdt.mk$gasd5.to.mk)
  deg2 = VariableFeatures(sr.pdt.mk$mk.cs8)
  spe.score.seurat <- CorrComparePlot(ExpressionTableSpecies1 = expr1, DEgenesSpecies1 = deg1,
                                      ExpressionTableSpecies2 = expr2, DEgenesSpecies2 = deg2,
                                      Species1 = "human", Species2 = "human",
                                      filename = file.path(paste0(sr.out, "annotation_by_monkeyCS8/dim", i),
                                                           paste0("Specificity_score_correlation_plot_of_All_cell_types_in_GAS_D5_and_monkeyCS8")))
  rm(expr1, expr2, deg1, deg2, spe.score.seurat, p)
}; rm(i)

# - integrating huamn and monkey
# genes filtering
tmp <- rs109.hs.mk.homo %>% filter(MK_Gene_SYMBOL!="" & Type=="ortholog_one2one")
tmp <- tmp[tmp$Hs_Gene_SYMBOL==tmp$MK_Gene_SYMBOL,]
tmp <- tmp %>% filter(Hs_Gene_SYMBOL%in%rownames(hs.cs7) & MK_Gene_SYMBOL%in%rownames(mk.cs8))
sr.pdt.hs.mk <- list()
sr.pdt.hs.mk$mk.cs8 <- mk.cs8[rownames(mk.cs8)%in%c(tmp$MK_Gene_SYMBOL, intersect(intersect(setdiff(as.vector(unlist(mk.ge.mk)), tmp$Hs_Gene_SYMBOL),
                                                                                            rownames(hs.cs7)), rownames(mk.cs8))),]
sr.pdt.hs.mk$hs.cs7 <- hs.cs7[rownames(hs.cs7)%in%c(tmp$Hs_Gene_SYMBOL, intersect(intersect(setdiff(as.vector(unlist(mk.ge.mk)), tmp$Hs_Gene_SYMBOL),
                                                                                            rownames(hs.cs7)), rownames(mk.cs8))),]

sr.pdt.hs.mk$gas.d7 <- sr.pdt$gasd5.to.hs[rownames(sr.pdt$gasd5.to.hs)%in%c(tmp$Hs_Gene_SYMBOL, intersect(intersect(setdiff(as.vector(unlist(mk.ge.mk)), tmp$Hs_Gene_SYMBOL),
                                                                                            rownames(sr.pdt$gasd5.to.hs)), rownames(mk.cs8))),]

# integrating
anchors <- FindTransferAnchors(reference = sr.pdt.hs.mk$mk.cs8, query = sr.pdt.hs.mk$hs.cs7, dims = 1:25, reference.reduction = "pca")
sr.pdt.hs.mk$hs.to.mk <- MapQuery(anchorset = anchors, reference = sr.pdt.hs.mk$mk.cs8, query = sr.pdt.hs.mk$hs.cs7,
                                  refdata = list(celltype = "CellType"), reference.reduction = "pca", reduction.model = "umap")
anchors <- FindTransferAnchors(reference = sr.pdt.hs.mk$mk.cs8, query = sr.pdt.hs.mk$gas.d7, dims = 1:25, reference.reduction = "pca")
sr.pdt.hs.mk$gas.to.mk <- MapQuery(anchorset = anchors, reference = sr.pdt.hs.mk$mk.cs8, query = sr.pdt.hs.mk$gas.d7,
                                   refdata = list(celltype = "CellType"), reference.reduction = "pca", reduction.model = "umap")
mk.col2 <- c("#607d8b","#795548","#ff5722","#ffc107","#cddc39","#4caf50","#009688",
             "#00bcd4","#2196f3","#3f51b5","#673ab7","#9c27b0","#e91e63","#f44336",
             "#b0bec5","#bcaaa4","#ffab91","#ffe082","#e6ee9c","#a5d6a7","#80cbc4",
             "#80deea","#90caf9","#9fa8da","#b39ddb","#ce93d8","#f48fb1","#ef9a9a",
             "#37474f","#4e342e","#d84315","#ff8f00","#9e9d24","#2e7d32","#00695c",
             "#00838f","#1565c0","#283593","#4527a0","#6a1b9a","#ad1457","#c62828")
mk.col2.fix <- mk.col2[1:length(unique(sr.pdt.hs.mk$mk.cs8$CellType))]
names(mk.col2.fix) <- sort(unique(sr.pdt.hs.mk$mk.cs8$CellType))
sr.pdt.hs.mk$mk.cs8$CellType <- as.character(sr.pdt.hs.mk$mk.cs8$CellType)
p1 <- DimPlot(sr.pdt.hs.mk$mk.cs8, reduction = "umap", group.by = "CellType", label = F, pt.size = 0.75) +
  ggtitle("Monkey embryos") +
  scale_color_manual(values = mk.col2.fix) +
  xlim(c(-15, 20)) +
  ylim(c(-10, 15))
p2 <- DimPlot(sr.pdt.hs.mk$hs.to.mk, reduction = "ref.umap", group.by = "predicted.celltype", label = F, pt.size = 0.75) +
  ggtitle("Human embryos") +
  scale_color_manual(values = mk.col2.fix) +
  xlim(c(-15, 20)) +
  ylim(c(-10, 15))
p3 <- DimPlot(sr.pdt.hs.mk$gas.to.mk, reduction = "ref.umap", group.by = "predicted.celltype", label = F, pt.size = 0.75) +
  ggtitle("Human embryos") +
  scale_color_manual(values = mk.col2.fix) +
  xlim(c(-15, 20)) +
  ylim(c(-10, 15))
pdf(file.path(sr.out, "Integration_of_Hs_CS7_GAS_D5_and_monkey_CS8_embryo_visualized_separately_umap.pdf"),
    height = 4.5, width = 20)
print(p1 + p2 + p3)
dev.off()
# with raw cell types
p1 <- DimPlot(sr.pdt.hs.mk$mk.cs8, reduction = "umap", group.by = "CellType", label = F, pt.size = 0.75) +
  ggtitle("Monkey embryos") +
  scale_color_manual(values = mk.col2.fix) +
  xlim(c(-15, 20)) +
  ylim(c(-10, 15))
p2 <- DimPlot(sr.pdt.hs.mk$hs.to.mk, reduction = "ref.umap", group.by = "CellType", label = F, pt.size = 0.75) +
  ggtitle("Human embryos") +
  scale_color_manual(values = pd.col.fix) +
  xlim(c(-15, 20)) +
  ylim(c(-10, 15))
p3 <- DimPlot(sr.pdt.hs.mk$gas.to.mk, reduction = "ref.umap", group.by = "CellType", label = F, pt.size = 0.75) +
  ggtitle("Human embryos") +
  scale_color_manual(values = pd.col.fix) +
  xlim(c(-15, 20)) +
  ylim(c(-10, 15))
pdf(file.path(sr.out, "Integration_of_Hs_CS7_GAS_D5_with_raw_cell_types_and_monkey_CS8_embryo_visualized_separately_umap.pdf"),
    height = 4.5, width = 18)
print(p1 + p2 + p3)
dev.off()
p1 <- Embeddings(sr.pdt.hs.mk$mk.cs8, reduction = "umap") %>% as.data.frame() %>% mutate(Datasets = "Monkey")
colnames(p1) <- c("UMAP_1", "UMAP_2", "Datasets")
p2 <- Embeddings(sr.pdt.hs.mk$hs.to.mk, reduction = "ref.umap") %>% as.data.frame() %>% mutate(Datasets = "Human")
colnames(p2) <- c("UMAP_1", "UMAP_2", "Datasets")
p3 <- Embeddings(sr.pdt.hs.mk$gas.to.mk, reduction = "ref.umap") %>% as.data.frame() %>% mutate(Datasets = "GAS_D5")
colnames(p3) <- c("UMAP_1", "UMAP_2", "Datasets")
pdf(file.path(sr.out, "Integration_of_Hs_CS7_GAS_D5_with_raw_cell_types_and_monkey_CS8_embryo_visualized_merged_umap_pt1.pdf"),
    height = 5, width = 6.5)
rbind(p1, p2, p3) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Datasets), size = 1) +
  theme_classic() +
  xlim(c(-15, 20)) +
  ylim(c(-10, 15))
dev.off()
pdf(file.path(sr.out, "Integration_of_Hs_CS7_GAS_D5_with_raw_cell_types_and_monkey_CS8_embryo_visualized_merged_umap_pt0.75.pdf"),
    height = 5, width = 6.5)
rbind(p1, p2, p3) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Datasets), size = 0.75) +
  theme_classic() +
  xlim(c(-15, 20)) +
  ylim(c(-10, 15))
dev.off()
pdf(file.path(sr.out, "Integration_of_Hs_CS7_GAS_D5_with_raw_cell_types_and_monkey_CS8_embryo_visualized_merged_umap_pt0.5.pdf"),
    height = 5, width = 6.5)
rbind(p1, p2, p3) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Datasets), size = 0.5) +
  theme_classic() +
  xlim(c(-15, 20)) +
  ylim(c(-10, 15))
dev.off()
pdf(file.path(sr.out, "Integration_of_Hs_CS7_GAS_D5_with_raw_cell_types_and_monkey_CS8_embryo_visualized_merged_umap_pt0.25.pdf"),
    height = 5, width = 6.5)
rbind(p1, p2, p3) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Datasets), size = 0.25) +
  theme_classic() +
  xlim(c(-15, 20)) +
  ylim(c(-10, 15))
dev.off()
# plot expression
sr.out <- file.path(sr.out, "Marker_gene_plots")
dir.create(sr.out, recursive = T)
marker.gene <- list(FMH.marker  = c("LHX5","ZNF616","GBX2","RREB1","TEAD4","EN2","HAND1",
                                    "HESX1","FOXG1","NFATC4","PAX7","DBX1","PAX5","SMARCA5",
                                    "PAX6","IRX3","EBF1","HOXA3","MEIS1","WRNIP1","MAFB"),
                    FMH = c("PAX6","OTX2","EN1","PAX8","GBX2","EGR2","PAX3","DBX2","DBX1","IRX3","OLIG2","NKX2.2"),
                    Neural = c("SOX1","SOX3","PAX6","TUBB3","ELAVL3","OLIG2","NEUROG1","NEUROD1","NIKX2.2"),
                    NMP = c("LIMD1","MARK3","CIT","MAPK14","FAT4","LATS2","STK4","WWC1","SHANK2","SOX11","TEAD1","PJA2",
                            "TEAD2","VGLL5","SAV1","MAP2K3","TIAL1","NF2","STK3","LATS1","THBS1","TEAD3","AMOT","YAP1",
                            "AMOTL1","AMOTL2","NEK8","WTIP","DLG5","TEAD4"),
                    Amnion = c("GABRP","TGFB1","VTCN1","S100A10","GADD45G","ACTC1",
                               "IGFBP3","MCM5","CLDN10","TUBG1","UCP2","BIN3","GGA2","FGFR1","GALM"),
                    EPI = c("CDH1", "POU5F1", "SOX2", "NANOG", "LIN28AP1", "MYC", "KLF4", "TDGF1", "GDF3", "GFP"),
                    Organizer = c("FST", "OTX2", "GSC", "EOMES", "FOXH1", "CER1", "MIXL1", "DKK1", "NODAL",
                                  "FOXA2", "LEFTY1", "LEFTY2", "TBXT", "SOX17", "SP5", "NOG", "CHRD"),
                    Embryonic = c("TBXT", "MESP1", "PDGFRA", "LHX1", "IRX3", "GATA6", "MYL7", "MSX1",
                                  "MSX2", "FOXF1", "HAND1", "HAND2", "BMP4", "MAB21L2", "CST1", "MEF2C",
                                  "PECAM1", "CDH5", "TEK", "TFAP2A", "TFAP2C", "GATA3", "DLX5", "LHX5",
                                  "PAX6", "SOX3", "PAX3", "PAX7"),
                    Extraembryonic = c("LUM", "NID2", "ANXA1", "POSTN", "GATA2", "KRT7", "VIM", "GATA4",
                                       "DCN", "TP63", "KRT18", "NR2F2", "ISL1", "HEY1", "CDH10", "CTSV",
                                       "ARID5B", "PLAGL1", "CREB3L1", "HOXA10", "HOXA11", "HOXA9",
                                       "HOXA13", "WNT6", "GABRP", "HEY1", "BST2"),
                    anterior = c("POU5F1","NANOG","TDGF1","NES","ID1","CDH1","SOX2","POU3F1","SOX3","IRX2",
                                 "CRABP1","PBX1","BTG1","TFAP2A","SFRP1","SFRP2","OTX2","FGFR1","HESX1","LHX5"),
                    medial = c("TBXT","TBX6","MESP1","EOMES","IRX3","LHX1","FOXC1","HOXB1","HOXA1","MEIS2","HOXB2",
                               "HOXB3","HOXA3","DLL1","GATA6","GATA4","FST","MIXL1","SP5","CER1","GFP"),
                    posterior = c("CDX2","CDX4","HOXB9","HOXC4","HOXC8","HOXC6","HOXB7","HOXB8","HOXB6","HOXB5",
                                  "VIM","SHH","FOXA2","CDH2","DUSP6","HOXA10","HOXA13","HAND1","HAND2","MSX1","GFP"),
                    megakaryocyte = c("GP1BB", "ITGA2B", "NFE2"),
                    erythroid = c("GATA1", "KLF1", "GYPB", "HBE1"),
                    myeloid.progenitor = c("CD36", "CSF1R", "LYVE1", "KIT", "CSF1R", "MYB", "SPI1", "CD34", "CD45", "CD52", "NFE2"),
                    wanghongmei = c("AFP","ANXA1","ANXA8","APOE","APOE ","CDH1","CDH2","CDX1","CDX2","CER1","CHRD",
                                    "COL1A1","COL3A1","DNMT3B","DPPA3","EOMES","FGF8","FGF8 HES7","FOLR2","FOXA2",
                                    "FOXC1","FOXC2","FOXF1","GATA1","GATA4","GATA6","GP1BB","HAND1","HBE1","HBM",
                                    "HBZ","HES7","HHEX","IGFBP5","IRX3","ISL1","KDR","KRT18","KRT8","LEFTY2","LFNG",
                                    "MEF2C","MIXL1","MPO","MYL7","NANOG","NANOS3","NKX2-5","NOTO","OSR1","OTX2","PAX3",
                                    "PAX3 ","PAX6","PAX7","PAX8","PECAM1","POU3F1","POU5F1","POU5FA","PPBP","PRDM14",
                                    "RUNX1","SIX1","SOX10","SOX17","SOX2","SOX2 ","SOX9","TBXT","TBX6","TCF15","TCF21",
                                    "TFAP2A","TFAP2C","TFC15","MESP1","TTR","UCHL1","VWF","WNT3A","WNT5B","WT1"))
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], Reduce(intersect, list(rownames(sr.pdt.hs.mk$mk.cs8),
                                                                               rownames(sr.pdt.hs.mk$hs.cs7),
                                                                               rownames(sr.pdt.hs.mk$gas.d7))))
}
for (type in "wanghongmei") {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 12.5, width = 5)
    p1 <- FeaturePlot(sr.pdt.hs.mk$mk.cs8, features = marker.gene[[type]][g], order = T, reduction = "umap",
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 0.75, slot = "data") +
      xlim(c(-15, 20)) +
      ylim(c(-10, 15))
    p2 <- FeaturePlot(sr.pdt.hs.mk$hs.to.mk, features = marker.gene[[type]][g], order = T, reduction = "ref.umap",
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 0.75, slot = "data") +
      xlim(c(-15, 20)) +
      ylim(c(-10, 15))
    p3 <- FeaturePlot(sr.pdt.hs.mk$gas.to.mk, features = marker.gene[[type]][g], order = T, reduction = "ref.umap",
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 0.75, slot = "data") +
      xlim(c(-15, 20)) +
      ylim(c(-10, 15))

    print(p1 + p2 + p3)
    dev.off()
  }; rm(g)
}
# plot relationships
Vis.Annotation.Relationship(sr.meta = sr.pdt.hs.mk$gas.to.mk@meta.data,
                            annotation = c("Sample_Name", "CellType", "predicted.celltype"),
                            pd.title = "Relationship between various annotations",
                            pd.col = NULL, split.meta = TRUE, split.by = "CellType",
                            pd.height = 15, pd.width = 15,
                            res.out = file.path(sr.out, "Sankey_plot_GAS_D5_sample"))
Vis.Annotation.Relationship(sr.meta = sr.pdt.hs.mk$gas.to.mk@meta.data,
                            annotation = c("Sample_Name", "CellType", "predicted.celltype"),
                            pd.title = "Relationship between various annotations",
                            pd.col = NULL, split.meta = FALSE, split.by = "CellType",
                            pd.height = 15, pd.width = 15,
                            res.out = file.path(sr.out, "Sankey_plot_GAS_D5_sample"))
Vis.Annotation.Relationship(sr.meta = sr.pdt.hs.mk$hs.to.mk@meta.data,
                            annotation = c("Tissue", "CellType", "predicted.celltype"),
                            pd.title = "Relationship between various annotations",
                            pd.col = NULL, split.meta = TRUE, split.by = "CellType",
                            pd.height = 15, pd.width = 15,
                            res.out = file.path(sr.out, "Sankey_plot_Hs_CS7_sample"))
Vis.Annotation.Relationship(sr.meta = sr.pdt.hs.mk$hs.to.mk@meta.data,
                            annotation = c("Tissue", "CellType", "predicted.celltype"),
                            pd.title = "Relationship between various annotations",
                            pd.col = NULL, split.meta = FALSE, split.by = "CellType",
                            pd.height = 15, pd.width = 15,
                            res.out = file.path(sr.out, "Sankey_plot_Hs_CS7_sample"))
# correlation analysis with monkey scRNAseq data (GAS vs Monkey)
expr.monkey <- AverageExpression(sr.pdt.hs.mk$mk.cs8, assays = "RNA", slot = "data", group.by = "CellType")
expr.gas <- AverageExpression(sr.pdt.hs.mk$gas.to.mk, assays = "RNA", slot = "data", group.by = "predicted.celltype")
hvg.monkey <- VariableFeatures(sr.pdt.hs.mk$mk.cs8)
hvg.gas <- VariableFeatures(sr.pdt.hs.mk$gas.to.mk)
homo.gene <- data.frame(species.1 = intersect(rownames(sr.pdt.hs.mk$mk.cs8), rownames(sr.pdt.hs.mk$gas.to.mk)),
                        species.2 = intersect(rownames(sr.pdt.hs.mk$mk.cs8), rownames(sr.pdt.hs.mk$gas.to.mk)))
homo.gene[nrow(homo.gene)+1,] <- c("TBXT", "T")
CorrComparePlot(ExpressionTableSpecies1 = expr.monkey$RNA, DEgenesSpecies1 = hvg.monkey,
                ExpressionTableSpecies2 = expr.gas$RNA, DEgenesSpecies2 = hvg.gas,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(sr.out, "Correlation_between_monkey_embryo_and_GAS_D5_with_monkey_cell_types"),
                Height = 15, Width = 15)
# correlation analysis with monkey scRNAseq data (GAS vs Human)
expr.hs <- AverageExpression(sr.pdt.hs.mk$hs.to.mk, assays = "RNA", slot = "data", group.by = "predicted.celltype")
expr.gas <- AverageExpression(sr.pdt.hs.mk$gas.to.mk, assays = "RNA", slot = "data", group.by = "predicted.celltype")
hvg.hs <- VariableFeatures(sr.pdt.hs.mk$hs.to.mk)
hvg.gas <- VariableFeatures(sr.pdt.hs.mk$gas.to.mk)
homo.gene <- data.frame(species.1 = intersect(rownames(sr.pdt.hs.mk$hs.to.mk), rownames(sr.pdt.hs.mk$gas.to.mk)),
                        species.2 = intersect(rownames(sr.pdt.hs.mk$hs.to.mk), rownames(sr.pdt.hs.mk$gas.to.mk)))
CorrComparePlot(ExpressionTableSpecies1 = expr.hs$RNA, DEgenesSpecies1 = hvg.hs,
                ExpressionTableSpecies2 = expr.gas$RNA, DEgenesSpecies2 = hvg.gas,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(sr.out, "Correlation_between_human_embryo_and_GAS_D5_with_monkey_cell_types"),
                Height = 15, Width = 15)
# correlation analysis with monkey scRNAseq data (GAS vs Monkey)
expr.monkey <- AverageExpression(sr.pdt.hs.mk$mk.cs8, assays = "RNA", slot = "data", group.by = "CellType")
expr.gas <- AverageExpression(sr.pdt.hs.mk$gas.to.mk, assays = "RNA", slot = "data", group.by = "CellType")
hvg.monkey <- VariableFeatures(sr.pdt.hs.mk$mk.cs8)
hvg.gas <- VariableFeatures(sr.pdt.hs.mk$gas.to.mk)
homo.gene <- data.frame(species.1 = intersect(rownames(sr.pdt.hs.mk$mk.cs8), rownames(sr.pdt.hs.mk$gas.to.mk)),
                        species.2 = intersect(rownames(sr.pdt.hs.mk$mk.cs8), rownames(sr.pdt.hs.mk$gas.to.mk)))
homo.gene[nrow(homo.gene)+1,] <- c("TBXT", "T")
CorrComparePlot(ExpressionTableSpecies1 = expr.monkey$RNA, DEgenesSpecies1 = hvg.monkey,
                ExpressionTableSpecies2 = expr.gas$RNA, DEgenesSpecies2 = hvg.gas,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(sr.out, "Correlation_between_monkey_embryo_and_GAS_D5_with_human_cell_types"),
                Height = 15, Width = 15)
# correlation analysis with monkey scRNAseq data (GAS vs Human)
expr.hs <- AverageExpression(sr.pdt.hs.mk$hs.to.mk, assays = "RNA", slot = "data", group.by = "CellType")
expr.gas <- AverageExpression(sr.pdt.hs.mk$gas.to.mk, assays = "RNA", slot = "data", group.by = "CellType")
hvg.hs <- VariableFeatures(sr.pdt.hs.mk$hs.to.mk)
hvg.gas <- VariableFeatures(sr.pdt.hs.mk$gas.to.mk)
homo.gene <- data.frame(species.1 = intersect(rownames(sr.pdt.hs.mk$hs.to.mk), rownames(sr.pdt.hs.mk$gas.to.mk)),
                        species.2 = intersect(rownames(sr.pdt.hs.mk$hs.to.mk), rownames(sr.pdt.hs.mk$gas.to.mk)))
CorrComparePlot(ExpressionTableSpecies1 = expr.hs$RNA, DEgenesSpecies1 = hvg.hs,
                ExpressionTableSpecies2 = expr.gas$RNA, DEgenesSpecies2 = hvg.gas,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(sr.out, "Correlation_between_human_embryo_and_GAS_D5_with_human_cell_types"),
                Height = 15, Width = 15)
# correlation analysis with monkey scRNAseq data (Human vs Monkey)
expr.monkey <- AverageExpression(sr.pdt.hs.mk$mk.cs8, assays = "RNA", slot = "data", group.by = "CellType")
expr.hs <- AverageExpression(sr.pdt.hs.mk$hs.to.mk, assays = "RNA", slot = "data", group.by = "predicted.celltype")
hvg.monkey <- VariableFeatures(sr.pdt.hs.mk$mk.cs8)
hvg.hs <- VariableFeatures(sr.pdt.hs.mk$hs.to.mk)
homo.gene <- data.frame(species.1 = intersect(rownames(sr.pdt.hs.mk$mk.cs8), rownames(sr.pdt.hs.mk$hs.to.mk)),
                        species.2 = intersect(rownames(sr.pdt.hs.mk$mk.cs8), rownames(sr.pdt.hs.mk$hs.to.mk)))
homo.gene[nrow(homo.gene)+1,] <- c("TBXT", "T")
CorrComparePlot(ExpressionTableSpecies1 = expr.monkey$RNA, DEgenesSpecies1 = hvg.monkey,
                ExpressionTableSpecies2 = expr.hs$RNA, DEgenesSpecies2 = hvg.hs,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(sr.out, "Correlation_between_monkey_embryo_and_Hs_CS7_with_monkey_cell_types"),
                Height = 15, Width = 15)
# GFP expression
pdf(file.path(sr.out, "GAS_D5_GFP_expression_with_monkey_cell_types.pdf"), height = 4.5, width = 5)
FeaturePlot(sr.pdt.hs.mk$gas.to.mk, features = "GFP_count", order = T, reduction = "ref.umap",
            cols = c("#ebebeb", "#7d3c98"), pt.size = 0.75, slot = "data") +
  xlim(c(-15, 20)) +
  ylim(c(-10, 15))
dev.off()
# Cell type and GFP type
Vis.Annotation.Ratio(sr.pdt.hs.mk$gas.to.mk@meta.data, annotation = c("predicted.celltype", "Type"),
                     pd.title = "GAS D5 annotated by Monkey data", pd.col = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279"),
                     pd.height = 10, pd.width = 10, res.out = file.path(sr.out, "GAS_D5_Cell_Type_Ratio"))


### >>> 3. Find marker genes
hs.tfs <- read.table("/home/yhw/document/AnimalTFDBs/Homo_sapiens_TF.txt", header = T, sep = "\t")
sr.out <- "/home/yhw/bioinfo/project-wb/hs_gastruloids/analysis/graphs/seurat_output/annotation_by_monkeyCS8"
Idents(mk.cs8) <- mk.cs8$CellType
mk.cs8 <- ScaleData(mk.cs8, features = rownames(mk.cs8))
library(future)
plan("multiprocess", workers = 6)
mk.cs8.markers <- FindAllMarkers(mk.cs8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(mk.cs8.markers, file.path(sr.out, "monkey_CS8_embryo_all_marker_genes.csv"), col.names = T)
pd.gene <- mk.cs8.markers[-grep("^LOC", mk.cs8.markers$gene),] %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
pdf(file.path(sr.out, "monkey_CS8_embryo_top5_marker_genes.pdf"), height = 15, width = 15)
DoHeatmap(mk.cs8, features = pd.gene$gene, slot = "scale.data")
dev.off()
pd.gene <- mk.cs8.markers %>%
  filter(cluster == "Gut") %>%
  filter(gene %in% hs.tfs$Symbol)
pdf(file.path(sr.out, "monkey_CS8_embryo_Gut_top20_marker_genes.pdf"), height = 15, width = 15)
VlnPlot(mk.cs8, features = pd.gene$gene, slot = "scale.data", stack = T, flip = T)
dev.off()
pdf(file.path(sr.out, "Scatter_plot_to_show_expression_level_of_monkey_CS8_embryo_Gut_marker_genes_on_monkey_umap.pdf"),
    height = 45, width = 19)
FeaturePlot(sr.pdt.hs.mk$mk.cs8, reduction = "umap", features = pd.gene$gene,
            cols = c("#ebebeb", "#7d3c98"), slot = "scale.data")
dev.off()
pdf(file.path(sr.out, "Scatter_plot_to_show_expression_level_of_monkey_CS8_embryo_Gut_marker_genes_on_GAS_umap.pdf"),
    height = 45, width = 19)
FeaturePlot(sr.pdt.hs.mk$gas.to.mk, reduction = "ref.umap", features = pd.gene$gene,
            cols = c("#ebebeb", "#7d3c98"), slot = "scale.data")
dev.off()



### ================================================
### 4th part: Trajectory time analysis with Monocle2 ----
### ================================================

### >>> 1. function loading
#' Perform trajectory analysis with Monocle2.
#'
#' This function uses the Monocle2 to infer trajectory.
#'
#' @param sr.obj A Seurat object.
#' @param outdir Output directory.
#' @param sp.name Prefix of output pdf file.
#' @param gp.name Column name to group cells.
#' @param pd.color The named color vector to map cell types.
#' @param point.size The point size.
#'
#' @import monocle
#' @import Seurat
#' @import ggplot2
#' @import dplyr
#'
#' @return A monocle object.
#' @export
#'
TI.Monocle2 <- function(sr.obj, outdir, sp.name, gp.name, pd.color, point.size) {
  library(monocle)
  library(Seurat)
  # - 1. check parameters
  if (!(class(sr.obj)[1] %in% c("Seurat"))) {
    stop("[Error] The data type is not supported by this function !!!")
  }

  # - 2. create a new CellDataSet object
  DefaultAssay(sr.obj) <- "RNA"
  sr.obj@meta.data$group <- sr.obj@meta.data[, gp.name]
  expr.mtx <- GetAssayData(sr.obj, assay = "RNA", slot = "counts")
  sp.meta <- sr.obj@meta.data
  gene.anno <- data.frame(
    gene_short_name = rownames(sr.obj),
    row.names = row.names(sr.obj),
    gene_id = row.names(sr.obj)
  )
  pd <- new("AnnotatedDataFrame", data = sp.meta)
  fd <- new("AnnotatedDataFrame", data = gene.anno)
  mncle <- newCellDataSet(expr.mtx, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

  # - 3. estimate size factors and dispersions
  mncle <- estimateSizeFactors(mncle)
  mncle <- estimateDispersions(mncle)

  # - 4. filtering low-quality cells
  mncle <- detectGenes(mncle, min_expr = 0.1)
  ordering.genes <- VariableFeatures(sr.obj)
  mncle <- setOrderingFilter(mncle, ordering.genes)

  # - 5. reduce data dimensionality
  mncle <- reduceDimension(mncle, max_components = 2, method = "ICA")

  # - 6. order cells along the trajectory
  mncle <- orderCells(mncle)

  # - 7. visualize the trajectory in the reduced dimensional space
  pdf(paste(outdir, sp.name, "_cell_trajectory_color_by_", gp.name, ".pdf", sep = ""), height = 4.5, width = 4)
  print(plot_cell_trajectory(mncle, color_by = gp.name, cell_size = point.size, cell_name_size = 2, show_branch_points = F) +
          scale_color_manual(values = pd.color[sort(as.character(unique(sr.obj@meta.data[, gp.name])))]))
  dev.off()

  # - 8. output data
  return(mncle)
}


### >>> 2. run Monocle2
mncle.ob <- list()
# - all cells in GAS D5
mncle.ob$gasd5.all <- TI.Monocle2(sr.pdt$gasd5.to.hs, "/home/cmq/bioinfo/project-tmp/wb/analysis/graphs/mncle2/",
                                  "GAS_D5_all_cells", "CellType", pd.col.fix, 1)
plot_cell_trajectory(mncle.ob$gasd5.all, color_by = "State", cell_size = 2, cell_name_size = 2)
mncle.ob$gasd5.all <- orderCells(mncle.ob$gasd5.all, root_state = "1")
pdf(paste("/home/cmq/bioinfo/project-tmp/wb/analysis/graphs/mncle2/", "GAS_D5_all_cells", "_cell_trajectory_color_by_pseudotime.pdf", sep = ""),
    height = 4.5, width = 4)
plot_cell_trajectory(mncle.ob$gasd5.all, color_by = "Pseudotime", cell_size = 1, cell_name_size = 2, show_branch_points = F)
dev.off()
# - GFP_pos cells in GAS D5
sr.mncle <- list()
sr.mncle$gasd5.gfp.pos <- subset(sr.pdt$gasd5.to.hs, Type=="GFP_pos")
sr.mncle$gasd5.gfp.pos <- ScaleData(sr.mncle$gasd5.gfp.pos, rownames(sr.mncle$gasd5.gfp.pos)) %>%
  FindVariableFeatures()
mncle.ob$gasd5.gfp.pos <- TI.Monocle2(sr.mncle$gasd5.gfp.pos, "/home/cmq/bioinfo/project-tmp/wb/analysis/graphs/mncle2/", "GAS_D5_GFP_pos_cells",
                                      "CellType", pd.col.fix, 4)
plot_cell_trajectory(mncle.ob$gasd5.gfp.pos, color_by = "State", cell_size = 2, cell_name_size = 2)
mncle.ob$gasd5.gfp.pos <- orderCells(mncle.ob$gasd5.gfp.pos, root_state = "3")
pdf(paste("/home/cmq/bioinfo/project-tmp/wb/analysis/graphs/mncle2/", "GAS_D5_GFP_pos_cells", "_cell_trajectory_color_by_pseudotime.pdf", sep = ""),
    height = 4.5, width = 4)
plot_cell_trajectory(mncle.ob$gasd5.gfp.pos, color_by = "Pseudotime", cell_size = 4, cell_name_size = 2, show_branch_points = F)
dev.off()
# - GFP_neg cells in GAS D5
sr.mncle$gasd5.gfp.neg <- subset(sr.pdt$gasd5.to.hs, Type=="GFP_neg")
sr.mncle$gasd5.gfp.neg <- ScaleData(sr.mncle$gasd5.gfp.neg, rownames(sr.mncle$gasd5.gfp.neg)) %>%
  FindVariableFeatures()
mncle.ob$gasd5.gfp.neg <- TI.Monocle2(sr.mncle$gasd5.gfp.neg, "/home/cmq/bioinfo/project-tmp/wb/analysis/graphs/mncle2/",
                                      "GAS_D5_GFP_neg_cells", "CellType", pd.col.fix, 1)
plot_cell_trajectory(mncle.ob$gasd5.gfp.neg, color_by = "State", cell_size = 2, cell_name_size = 2)
mncle.ob$gasd5.gfp.neg <- orderCells(mncle.ob$gasd5.gfp.neg, root_state = "1")
pdf(paste("/home/cmq/bioinfo/project-tmp/wb/analysis/graphs/mncle2/", "GAS_D5_GFP_neg_cells", "_cell_trajectory_color_by_pseudotime.pdf", sep = ""),
    height = 4.5, width = 4)
plot_cell_trajectory(mncle.ob$gasd5.gfp.neg, color_by = "Pseudotime", cell_size = 1, cell_name_size = 2, show_branch_points = F)
dev.off()
# - all cells in human CS7
sr.mncle$hs.cs7 <- hs.cs7 %>% FindVariableFeatures(nfeatures = 2500)
mncle.ob$hs.cs7 <- TI.Monocle2(sr.mncle$hs.cs7, "/home/cmq/bioinfo/project-tmp/wb/analysis/graphs/mncle2/",
                               "Human_CS7", "CellType", pd.col.fix, 2.5)
plot_cell_trajectory(mncle.ob$hs.cs7, color_by = "State", cell_size = 2, cell_name_size = 2)
mncle.ob$hs.cs7 <- orderCells(mncle.ob$hs.cs7, root_state = "")
pdf(paste("/home/cmq/bioinfo/project-tmp/wb/analysis/graphs/mncle2/", "GAS_D5_GFP_neg_cells", "_cell_trajectory_color_by_pseudotime.pdf", sep = ""),
    height = 4.5, width = 4)
plot_cell_trajectory(mncle.ob$hs.cs7, color_by = "Pseudotime", cell_size = 2.5, cell_name_size = 2, show_branch_points = F)
dev.off()
# -
sr.mncle$gas.merge <- gas.merge.urd
mncle.ob$gas.merge <- TI.Monocle2(sr.obj = sr.mncle$gas.merge,
                                  outdir = "/home/yhw/bioinfo/project-wb/hs_gastruloids/analysis/graphs/mncle2/",
                                  sp.name = "Human_GAS_merge",
                                  gp.name = "CellType", pd.color = gas.merge.col, point.size = 1)
plot_cell_trajectory(mncle.ob$gas.merge, color_by = "State", cell_size = 1, cell_name_size = 2)
mncle.ob$gas.merge <- orderCells(mncle.ob$gas.merge, root_state = "4")
pdf(paste("/home/yhw/bioinfo/project-wb/hs_gastruloids/analysis/graphs/mncle2", "GAS_D5_GFP_neg_cells", "_cell_trajectory_color_by_pseudotime.pdf", sep = ""),
    height = 4.5, width = 4)
plot_cell_trajectory(mncle.ob$gas.merge, color_by = "Pseudotime", cell_size = 1,
                     cell_name_size = 2, show_branch_points = F)
plot_cell_trajectory(mncle.ob$gas.merge, color_by = "Sample_Group", cell_size = 1,
                     cell_name_size = 2, show_branch_points = F)
plot_cell_trajectory(mncle.ob$gas.merge, color_by = "CellType", cell_size = 1,
                     cell_name_size = 2, show_branch_points = F) +
  scale_color_manual(values = gas.merge.col)
dev.off()



### ===========================================
### 5th part: Trajectory time analysis with URD ----
### ===========================================

### >>> 1. Load data
marker.genes <- c("TFAP2A","TFAP2C","DLX5","GATA3", # Ectoderm (ECT)
                  "CDH1","SOX2","POU5F1","OTX2", # Epiblast (EPI)
                  "TBXT","SP5","FST","MIXL1", # Primitive streak (PS)
                  "EOMES","MESP1","TBX6", # Nascent mesoderm (NM)
                  "PDGFRA","LHX1","IRX3","GATA6", # Emergent mesoderm (EM)
                  "FOXF1","HAND1","HAND2","MYL7","MSX1", # Advanced mesoderm (AM)
                  "LUM","NID2","POSTN","BMP4","ANXA1","MAB21L2", # Extraembryonic mesoderm (ExM)
                  "SOX17","GSC","FOXA2","CST1", # Endoderm (End)
                  "MEF2C","PECAM1","CDH5","TEK" # Haemato-endothelial progenitors (HEP)
)
pd <- hs.cs7
pd <- NormalizeData(pd) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)


### >>> 2. URD pipeline before building trees
urd.gasd5 <- list()
urd.hscs7 <- list()
# End_AM
urd.gasd5[["knn56_sigma16_End_AM"]] <- URD.Btree(sr.obj = sr.pdt$gasd5.to.hs, group.by = c("CellType", "Type"),
                                                 tip.gp = c("AM", "End"), start.cell = "EPI", urd.knn = 56, urd.sigma = 16,
                                                 group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                 pt.gene = marker.genes, res.out = "graphs/urd/GAS_D5_knn56_sigma16_End_AM")
urd.gasd5[["knn100_sigma16_End_AM"]] <- URD.Btree(sr.obj = sr.pdt$gasd5.to.hs, group.by = c("CellType", "Type"),
                                                  tip.gp = c("AM", "End"), start.cell = "EPI", urd.knn = 100, urd.sigma = 16,
                                                  group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                  pt.gene = marker.genes, res.out = "graphs/urd/GAS_D5_knn100_sigma16_End_AM")
urd.hscs7[["knn35_sigma0_End_AM"]] <- URD.Btree(sr.obj = pd, group.by = c("CellType"),
                                                tip.gp = c("AM", "End"), start.cell = "EPI", urd.knn = 35, urd.sigma = NULL,
                                                group.col = list(CellType = pd.col.fix),
                                                pt.gene = marker.genes, res.out = "graphs/urd/Hs_CS7_knn35_sigma0_End_AM")
urd.gasd5[["knn56_sigma0_End_AM"]] <- URD.Btree(sr.obj = sr.pdt$gasd5.to.hs, group.by = c("CellType", "Type"),
                                                tip.gp = c("AM", "End"), start.cell = "EPI", urd.knn = 56, urd.sigma = NULL,
                                                group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                pt.gene = marker.genes, res.out = "graphs/urd/GAS_D5_knn56_sigma0_End_AM")
# End_AM_ECT
urd.gasd5[["knn56_sigma16_End_AM_ECT"]] <- URD.Btree(sr.obj = sr.pdt$gasd5.to.hs, group.by = c("CellType", "Type"),
                                                     tip.gp = c("AM", "ECT", "End"), start.cell = "EPI", urd.knn = 56, urd.sigma = 16,
                                                     group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                     pt.gene = marker.genes, res.out = "graphs/urd/GAS_D5_knn56_sigma16_End_AM_ECT")
urd.gasd5[["knn100_sigma16_End_AM_ECT"]] <- URD.Btree(sr.obj = sr.pdt$gasd5.to.hs, group.by = c("CellType", "Type"),
                                                      tip.gp = c("AM", "ECT", "End"), start.cell = "EPI", urd.knn = 56, urd.sigma = 16,
                                                      group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                      pt.gene = marker.genes, res.out = "graphs/urd/GAS_D5_knn100_sigma16_End_AM_ECT")
urd.hscs7[["knn35_sigma0_End_AM_ECT"]] <- URD.Btree(sr.obj = pd, group.by = c("CellType"),
                                                    tip.gp = c("AM", "End", "ECT"), start.cell = "EPI", urd.knn = 35, urd.sigma = NULL,
                                                    group.col = list(CellType = pd.col.fix),
                                                    pt.gene = marker.genes, res.out = "graphs/urd/Hs_CS7_knn35_sigma0_End_AM_ECT")
urd.gasd5[["knn56_sigma0_End_AM_ECT"]] <- URD.Btree(sr.obj = sr.pdt$gasd5.to.hs, group.by = c("CellType", "Type"),
                                                    tip.gp = c("AM", "ECT", "End"), start.cell = "EPI", urd.knn = 56, urd.sigma = NULL,
                                                    group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                    pt.gene = marker.genes, res.out = "graphs/urd/GAS_D5_knn56_sigma0_End_AM_ECT")
# End_AM_ECT_ExM
urd.gasd5[["knn56_sigma16_End_AM_ECT_ExM"]] <- URD.Btree(sr.obj = sr.pdt$gasd5.to.hs, group.by = c("CellType", "Type"),
                                                         tip.gp = c("AM","ECT","End","ExM"), start.cell = "EPI", urd.knn = 56, urd.sigma = 16,
                                                         group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                         pt.gene = marker.genes, res.out = "graphs/urd/GAS_D5_knn56_sigma16_End_AM_ECT_ExM")
urd.gasd5[["knn100_sigma16_End_AM_ECT_ExM"]] <- URD.Btree(sr.obj = sr.pdt$gasd5.to.hs, group.by = c("CellType", "Type"),
                                                          tip.gp = c("AM","ECT","End","ExM"), start.cell = "EPI", urd.knn = 56, urd.sigma = 16,
                                                          group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                          pt.gene = marker.genes, res.out = "graphs/urd/GAS_D5_knn100_sigma16_End_AM_ECT_ExM")
urd.hscs7[["knn35_sigma0_End_AM_ECT_ExM"]] <- URD.Btree(sr.obj = pd, group.by = c("CellType"),
                                                        tip.gp = c("AM", "End", "ECT", "ExM"), start.cell = "EPI", urd.knn = 35, urd.sigma = NULL,
                                                        group.col = list(CellType = pd.col.fix),
                                                        pt.gene = marker.genes, res.out = "graphs/urd/Hs_CS7_knn35_sigma0_End_AM_ECT_ExM")
urd.gasd5[["knn56_sigma0_End_AM_ECT_ExM"]] <- URD.Btree(sr.obj = sr.pdt$gasd5.to.hs, group.by = c("CellType", "Type"),
                                                        tip.gp = c("AM", "ECT", "End", "ExM"), start.cell = "EPI", urd.knn = 56, urd.sigma = NULL,
                                                        group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                        pt.gene = marker.genes, res.out = "graphs/urd/GAS_D5_knn56_sigma0_End_AM_ECT_ExM")


### >>> 3. Building trees
tree.gasd5 <- list()
tree.hscs7 <- list()
#
tree.gasd5[["knn56_sigma16_End_AM"]] <- URD.Atree(urd = urd.gasd5[["knn56_sigma16_End_AM"]], tip.gp = c("AM", "End"),
                                                  tip.cluster = c("3", "6"), group.by = c("CellType", "Type"),
                                                  group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                  pt.gene = marker.genes, res.out = "graphs/urd/GAS_D5_knn56_sigma16_End_AM")
tree.gasd5[["knn100_sigma16_End_AM"]] <- URD.Atree(urd = urd.gasd5[["knn100_sigma16_End_AM"]], tip.gp = c("AM", "End"),
                                                   tip.cluster = c("3", "6"), group.by = c("CellType", "Type"),
                                                   group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                   pt.gene = marker.genes, res.out = "graphs/urd/GAS_D5_knn100_sigma16_End_AM")
tree.hscs7[["knn35_sigma0_End_AM"]] <- URD.Atree(urd = urd.hscs7[["knn35_sigma0_End_AM"]], tip.gp = c("AM", "End"),
                                                 tip.cluster = c("3", "1"), group.by = c("CellType"),
                                                 group.col = list(CellType = pd.col.fix),
                                                 pt.gene = marker.genes, res.out = "graphs/urd/Hs_CS7_knn35_sigma0_End_AM")
tree.gasd5[["knn56_sigma0_End_AM"]] <- URD.Atree(urd = urd.gasd5[["knn56_sigma0_End_AM"]], tip.gp = c("AM", "End"),
                                                 tip.cluster = c("3", "6"), group.by = c("CellType", "Type"),
                                                 group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                 pt.gene = marker.genes, res.out = "graphs/urd/GAS_D5_knn56_sigma0_End_AM")
# End_AM_ECT
tree.gasd5[["knn56_sigma16_End_AM_ECT"]] <- URD.Atree(urd = urd.gasd5[["knn56_sigma16_End_AM_ECT"]], tip.gp = c("AM", "End", "ECT"),
                                                      tip.cluster = c("6", "7", "2"), group.by = c("CellType", "Type"),
                                                      group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                      pt.gene = marker.genes, res.out = "graphs/urd/GAS_D5_knn56_sigma16_End_AM_ECT")
tree.gasd5[["knn100_sigma16_End_AM_ECT"]] <- URD.Atree(urd = urd.gasd5[["knn100_sigma16_End_AM_ECT"]], tip.gp = c("AM", "End", "ECT"),
                                                       tip.cluster = c("6", "7", "2"), group.by = c("CellType", "Type"),
                                                       group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                       pt.gene = marker.genes, res.out = "graphs/urd/GAS_D5_knn100_sigma16_End_AM_ECT")
tree.hscs7[["knn35_sigma0_End_AM_ECT"]] <- URD.Atree(urd = urd.hscs7[["knn35_sigma0_End_AM_ECT"]], tip.gp = c("AM", "End", "ECT"),
                                                     tip.cluster = c("3", "1", "6"), group.by = c("CellType"),
                                                     group.col = list(CellType = pd.col.fix),
                                                     pt.gene = marker.genes, res.out = "graphs/urd/Hs_CS7_knn35_sigma0_End_AM_ECT")
tree.gasd5[["knn56_sigma0_End_AM_ECT"]] <- URD.Atree(urd = urd.gasd5[["knn56_sigma0_End_AM_ECT"]], tip.gp = c("AM", "End", "ECT"),
                                                     tip.cluster = c("8", "7", "2"), group.by = c("CellType", "Type"),
                                                     group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                     pt.gene = marker.genes, res.out = "graphs/urd/GAS_D5_knn56_sigma0_End_AM_ECT")
# End_AM_ECT_ExM
tree.gasd5[["knn56_sigma16_End_AM_ECT_ExM"]] <- URD.Atree(urd = urd.gasd5[["knn56_sigma16_End_AM_ECT_ExM"]], tip.gp = c("AM", "End", "ECT", "ExM"),
                                                          tip.cluster = c("2", "7", "9", "3"), group.by = c("CellType", "Type"),
                                                          group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                          pt.gene = marker.genes, res.out = "graphs/urd/GAS_D5_knn56_sigma16_End_AM_ECT_ExM")
tree.gasd5[["knn100_sigma16_End_AM_ECT_ExM"]] <- URD.Atree(urd = urd.gasd5[["knn100_sigma16_End_AM_ECT"]], tip.gp = c("AM", "End", "ECT", "ExM"),
                                                           tip.cluster = c("2", "7", "9", "3"), group.by = c("CellType", "Type"),
                                                           group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                           pt.gene = marker.genes, res.out = "graphs/urd/GAS_D5_knn100_sigma16_End_AM_ECT_ExM")
tree.hscs7[["knn35_sigma0_End_AM_ECT_ExM"]] <- URD.Atree(urd = urd.hscs7[["knn35_sigma0_End_AM_ECT_ExM"]], tip.gp = c("AM", "End", "ECT", "ExM"),
                                                         tip.cluster = c("5", "2", "6", "3"), group.by = c("CellType"),
                                                         group.col = list(CellType = pd.col.fix),
                                                         pt.gene = marker.genes, res.out = "graphs/urd/Hs_CS7_knn35_sigma0_End_AM_ECT_ExM")
tree.gasd5[["knn56_sigma0_End_AM_ECT_ExM"]] <- URD.Atree(urd = urd.gasd5[["knn56_sigma0_End_AM_ECT_ExM"]], tip.gp = c("AM", "End", "ECT", "ExM"),
                                                         tip.cluster = c("2", "7", "9", "3"), group.by = c("CellType", "Type"),
                                                         group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                         pt.gene = marker.genes, res.out = "graphs/urd/GAS_D5_knn56_sigma0_End_AM_ECT_ExM")


### >>> 4. Re-plot
# gastruloids
pd.gene <- c("NCAM", "SOX3", "ERNI", "SIX3", "HOXA1", "HOXA2", "HOXB1", "HOXB2", "MAFB", "EPHA4", "EPHB2", "SOX10")
pd.gene <- c(unname(unlist(marker.gene)))
pd.gene <- intersect(pd.gene, rownames(urd.gasd5$knn56_sigma16_End_AM@logupx.data))
pd.gene <- "LUM"
for (name in names(urd.gasd5)[5]) {
  URD.Replot(urd = urd.gasd5[[name]], tree = tree.gasd5[[name]], mode = c("gene", "group"),
             group.by = c("CellType", "Type"),
             group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#29C169")),
             expr.col = c("#ebebeb", "#7d3c98"), pt.gene = pd.gene, pt.size = 1.5,
             res.out = paste0("graphs/urd/Replot_GAS_D5_", name))
}
# human embryo
pd.gene <- c("NCAM", "SOX3", "ERNI", "SIX3", "HOXA1", "HOXA2", "HOXB1", "HOXB2", "MAFB", "EPHA4", "EPHB2", "SOX10")
pd.gene <- c(unname(unlist(marker.gene)))
pd.gene <- intersect(pd.gene, rownames(urd.hscs7$knn35_sigma0_End_AM@logupx.data))
if (all(colnames(urd.hscs7[[name]]@logupx.data) == colnames(hs.cs7))) {
  urd.hscs7[[name]]@logupx.data <- GetAssayData(hs.cs7, slot = "data")
  tree.hscs7[[name]]@logupx.data <- GetAssayData(hs.cs7, slot = "data")
}
urd.hscs7[[name]]@meta
for (name in names(urd.hscs7)[2]) {
  URD.Replot(urd = urd.hscs7[[name]], tree = tree.hscs7[[name]], mode = c("gene", "group"),
             group.by = c("CellType"),
             group.col = list(CellType = pd.col.fix, Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#29C169")),
             expr.col = c("#ebebeb", "#7d3c98"), pt.gene = pd.gene, pt.size = 1.5,
             res.out = paste0("graphs/urd/Replot_Hs_CS7_", name))
}
# gene expression level along pseudotime
pd.gene <- list(cell = c("EPI","PS","End","NM","EM","AM","ECT","ExM","HEP"),
                gene = c("ANXA1","BMP4","CDH1","CDH5","CST1","DLX5","FOXA2","FOXF1",
                         "FST","GATA3","GATA6","HAND1","IRX3","LHX1","LUM","MAB21L2",
                         "MEF2C","MESP1","MYL7","NID2","OTX2","PDGFRA","PECAM1",
                         "POU5F1","SOX17","SOX2","SP5","TBX6","TBXT","TEK","TFAP2A"))
URD.ExprTime(urd = urd.gasd5$knn56_sigma16_End_AM_ECT_ExM, para = pd.gene, group.by = "CellType",
             line.cols = "#D9212A", se.col = "#D9212A",
             res.out = "graphs/urd/ExprTime_GAS_D5_knn56_sigma16_End_AM_ECT_ExM_geneset1")
URD.ExprTime(urd = urd.hscs7$knn35_sigma0_End_AM_ECT, para = pd.gene, group.by = "CellType",
             line.cols = "#045EC3", se.col = "#045EC3",
             res.out = "graphs/urd/ExprTime_Hs_CS7_knn35_sigma0_End_AM_ECT_geneset1")
#
pd.gene <- list(cell = c("EPI","PS","End","NM","EM","AM","ECT","ExM","HEP"),
                gene = c("AFP","ANXA3","APLNR","AQP3","CDX2","EOMES","FGF2","FGF8",
                         "FOXC1","FURIN","GABRP","HAND2","IGFBP3","ISL1","LHX5",
                         "MIXL1","MSX1","MSX2","NANOG","PAX6","PGF","PITX2","POSTN",
                         "PTN","SNAIL2","SNAL1","TDGF1","TFAP2B","TFAP2C","TTR","VTCN1","WNT11"))
URD.ExprTime(urd = urd.gasd5$knn56_sigma16_End_AM_ECT_ExM, para = pd.gene, group.by = "CellType",
             line.cols = "#D9212A", se.col = "#D9212A",
             res.out = "graphs/urd/ExprTime_GAS_D5_knn56_sigma16_End_AM_ECT_ExM_geneset2")
URD.ExprTime(urd = urd.hscs7$knn35_sigma0_End_AM_ECT, para = pd.gene, group.by = "CellType",
             line.cols = "#045EC3", se.col = "#045EC3",
             res.out = "graphs/urd/ExprTime_Hs_CS7_knn35_sigma0_End_AM_ECT_geneset2")
#
pd.gene <- list(cell = c("EPI","PS","End","NM","EM","AM","ECT","ExM","HEP"),
                gene = c("DBX1","DBX2","DLX1","EBF1","EGRZ","ELAVL3","EN1","EN2",
                         "FGF1","GBX2","HESX1","MAFB","NEUROD1","NEUROG1","NKX2.2",
                         "OLIG2","PAX3","PAX5","PAX6","PAX7","SIX3","SIX6","SOX1","SOX3","TUBB3"))
URD.ExprTime(urd = urd.gasd5$knn56_sigma16_End_AM_ECT_ExM, para = pd.gene, group.by = "CellType",
             line.cols = "#D9212A", se.col = "#D9212A",
             res.out = "graphs/urd/ExprTime_GAS_D5_knn56_sigma16_End_AM_ECT_ExM_geneset3")
URD.ExprTime(urd = urd.hscs7$knn35_sigma0_End_AM_ECT, para = pd.gene, group.by = "CellType",
             line.cols = "#045EC3", se.col = "#045EC3",
             res.out = "graphs/urd/ExprTime_Hs_CS7_knn35_sigma0_End_AM_ECT_geneset3")





### ===============================
### 6th part: RNA velocity analysis ----
### ===============================


### >>> 1. Setting output directory
outdir <- file.path(getwd(), "graphs/velocyto")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}
pd.gene <- c("AMOT","AMOTL1","AMOTL2","CIT","DBX1","DBX2","DLG5","EBF1","EGR2",
             "ELAVL3","EN1","EN2","FAT4","FOXG1","GBX2","HAND1","HESX1","HOXA3",
             "IRX3","LATS1","LATS2","LHX5","LIMD1","MAFB","MAP2K3","MAPK14","MARK3",
             "MEIS1","NEK8","NEUROD1","NEUROG1","NF2","NFATC4","NIKX2.2","NKX2.2",
             "OLIG2","OTX2","PAX3","PAX5","PAX6","PAX7","PAX8","PJA2","RREB1","SAV1",
             "SHANK2","SMARCA5","SOX1","SOX11","SOX3","STK3","STK4","TEAD1","TEAD2",
             "TEAD3","TEAD4","THBS1","TIAL1","TUBB3","VGLL5","WRNIP1","WTIP","WWC1","YAP1","ZNF616")
pd.gene <- intersect(pd.gene, rownames(hs.cs7))
neuron.genes <- list(FMH = c("DBX1","DBX2","EBF1","EGR2","EN1","EN2","FOXG1",
                             "GBX2","HAND1","HESX1","HOXA3","IRX3","LHX5",
                             "MAFB","MEIS1","NFATC4","NKX2.2","OLIG2","OTX2",
                             "PAX3","PAX5","PAX6","PAX7","PAX8","RREB1",
                             "SMARCA5","TEAD4","WRNIP1","ZNF616"),
                     Neural = c("ELAVL3","NEUROD1","NEUROG1","NIKX2.2",
                                "OLIG2","PAX6","SOX1","SOX3","TUBB3"),
                     NMP = c("AMOT","AMOTL1","AMOTL2","CIT","DLG5","FAT4",
                             "LATS1","LATS2","LIMD1","MAP2K3","MAPK14","MARK3",
                             "NEK8","NF2","PJA2","SAV1","SHANK2","SOX11","STK3",
                             "STK4","TEAD1","TEAD2","TEAD3","TEAD4","THBS1",
                             "TIAL1","VGLL5","WTIP","WWC1","YAP1"))
ave.exp <- AverageExpression(object = hs.cs7, assays = "RNA",
                             features = pd.gene, group.by = "CellType")
keep.gene <- rownames(ave.exp$RNA)[rowSums(ave.exp$RNA)!=0]
for (name in names(neuron.genes)) {
  neuron.genes[[name]] <- intersect(neuron.genes[[name]], keep.gene)
}


### >>> 2. GAS D5
# all cells --> ["#3b77b3", "#54beca", "#4f996e", "#b0c6e4", "#eb8838", "#f2adad", "#6d51a0", "#b7be70", "#cb7eaf"]
SeuratToVelocyto(sr.ob = sr.pdt$gasd5.to.hs,
                 reduction = "tsne", out.dir = file.path(outdir, "GAS_D5"),
                 sample.name = "GAS_D5_all", str.in.barcode = "_.*")
for (gene in pd.gene) {
  pdf(file.path(outdir, paste0("/GAS_D5GAS_D5_all_tsne/gene_expression_level_of_", gene, ".pdf")),
      height = 5, width = 6)
  print(FeaturePlot(sr.pdt$gasd5.to.hs, features = gene, reduction = "tsne", slot = "data",
                    cols = c("#ebebeb", "#7d3c98"), pt.size = 1.5))
  dev.off()
}
for (type in names(neuron.genes)) {
  pdf(file.path(outdir, paste0("/GAS_D5GAS_D5_all_tsne/Violin_plot_to_show_gene_expression_level_of_", type, ".pdf")),
      height = 10, width = 10)
  print(VlnPlot(sr.pdt$gasd5.to.hs, features = neuron.genes[[type]],
                group.by = "CellType", stack = T, flip = T, slot = "data"))
  dev.off()
}
SeuratPlotExpr(seu.obj = sr.pdt$gasd5.to.hs, gene.list = neuron.genes,
               dimension = "tsne", pd.assay = "RNA", pd.data = "data",
               group.by = "CellType", pt.size = 2, res.out = "/home/yhw/hahaha")
# "PS","End","NM","EM","AM" --> ["#3b77b3", "#4f996e", "#eb8838", "#b7be70", "#cb7eaf"]
set.seed(123)
gasd5.sub <- subset(sr.pdt$gasd5.to.hs, CellType %in% c("PS","End","NM","EM","AM"))
gasd5.sub <- NormalizeData(gasd5.sub) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
gasd5.sub <- ScaleData(gasd5.sub, verbose = FALSE, features = rownames(gasd5.sub))
gasd5.sub <- RunPCA(gasd5.sub, npcs = 30, verbose = FALSE, features = VariableFeatures(gasd5.sub))
ElbowPlot(gasd5.sub, ndims = 30)
dim.n <- 15
gasd5.sub <- FindNeighbors(gasd5.sub, reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 2)
gasd5.sub <- RunUMAP(gasd5.sub, reduction = "pca", dims = 1:dim.n) %>% RunTSNE(reduction = "pca", dims = 1:dim.n)
DimPlot(gasd5.sub, reduction = "tsne", group.by = "CellType", cols = pd.col.fix)
SeuratToVelocyto(sr.ob = gasd5.sub,
                 reduction = "tsne", out.dir = "graphs/velocyto/GAS_D5",
                 sample.name = "GAS_D5_PS_End_NM_EM_AM_Dim15", str.in.barcode = "_.*")
for (gene in pd.gene) {
  pdf(file.path(outdir, paste0("/GAS_D5GAS_D5_PS_End_NM_EM_AM_Dim15_tsne/gene_expression_level_of_", gene, ".pdf")),
      height = 5, width = 6)
  print(FeaturePlot(gasd5.sub, features = gene, reduction = "tsne", slot = "data",
                    cols = c("#ebebeb", "#7d3c98"), pt.size = 1.5))
  dev.off()
}
for (type in names(neuron.genes)) {
  pdf(file.path(outdir, paste0("/GAS_D5GAS_D5_PS_End_NM_EM_AM_Dim15_tsne/Violin_plot_to_show_gene_expression_level_of_", type, ".pdf")),
      height = 10, width = 10)
  print(VlnPlot(gasd5.sub, features = neuron.genes[[type]],
                group.by = "CellType", stack = T, flip = T, slot = "data"))
  dev.off()
}
# "PS","End","NM","EM","AM","ECT" --> ["#3b77b3", "#54beca", "#4f996e", "#eb8838", "#b7be70", "#cb7eaf"]
set.seed(123)
gasd5.sub2 <- subset(sr.pdt$gasd5.to.hs, CellType %in% c("PS","End","NM","EM","AM","ECT"))
gasd5.sub2 <- NormalizeData(gasd5.sub2) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
gasd5.sub2 <- ScaleData(gasd5.sub2, verbose = FALSE, features = rownames(gasd5.sub2))
gasd5.sub2 <- RunPCA(gasd5.sub2, npcs = 30, verbose = FALSE, features = VariableFeatures(gasd5.sub2))
ElbowPlot(gasd5.sub2, ndims = 30)
dim.n <- 15
gasd5.sub2 <- FindNeighbors(gasd5.sub2, reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 2)
gasd5.sub2 <- RunUMAP(gasd5.sub2, reduction = "pca", dims = 1:dim.n) %>% RunTSNE(reduction = "pca", dims = 1:dim.n)
DimPlot(gasd5.sub2, reduction = "tsne", group.by = "CellType", cols = pd.col.fix)
SeuratToVelocyto(sr.ob = gasd5.sub2,
                 reduction = "tsne", out.dir = "graphs/velocyto/GAS_D5",
                 sample.name = "GAS_D5_PS_End_NM_EM_AM_ECT_Dim15", str.in.barcode = "_.*")
for (gene in pd.gene) {
  pdf(file.path(outdir, paste0("/GAS_D5GAS_D5_PS_End_NM_EM_AM_ECT_Dim15_tsne/gene_expression_level_of_", gene, ".pdf")),
      height = 5, width = 6)
  print(FeaturePlot(gasd5.sub2, features = gene, reduction = "tsne", slot = "data",
                    cols = c("#ebebeb", "#7d3c98"), pt.size = 1.5))
  dev.off()
}
for (type in names(neuron.genes)) {
  pdf(file.path(outdir, paste0("/GAS_D5GAS_D5_PS_End_NM_EM_AM_ECT_Dim15_tsne/Violin_plot_to_show_gene_expression_level_of_", type, ".pdf")),
      height = 10, width = 10)
  print(VlnPlot(gasd5.sub2, features = neuron.genes[[type]],
                group.by = "CellType", stack = T, flip = T, slot = "data"))
  dev.off()
}
# "PS","End","NM","EM","AM","ECT","ExM" -->  ["#3b77b3", "#54beca", "#4f996e", "#eb8838", "#f2adad", "#b7be70", "#cb7eaf"]
set.seed(123)
gasd5.sub3 <- subset(sr.pdt$gasd5.to.hs, CellType %in% c("PS","End","NM","EM","AM","ECT","ExM"))
gasd5.sub3 <- NormalizeData(gasd5.sub3) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
gasd5.sub3 <- ScaleData(gasd5.sub3, verbose = FALSE, features = rownames(gasd5.sub3))
gasd5.sub3 <- RunPCA(gasd5.sub3, npcs = 30, verbose = FALSE, features = VariableFeatures(gasd5.sub3))
ElbowPlot(gasd5.sub3, ndims = 30)
dim.n <- 15
gasd5.sub3 <- FindNeighbors(gasd5.sub3, reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 2)
gasd5.sub3 <- RunUMAP(gasd5.sub3, reduction = "pca", dims = 1:dim.n) %>% RunTSNE(reduction = "pca", dims = 1:dim.n)
DimPlot(gasd5.sub3, reduction = "tsne", group.by = "CellType", cols = pd.col.fix)
SeuratToVelocyto(sr.ob = gasd5.sub3,
                 reduction = "tsne", out.dir = "graphs/velocyto/GAS_D5",
                 sample.name = "GAS_D5_PS_End_NM_EM_AM_ECT_ExM_Dim15", str.in.barcode = "_.*")
for (gene in pd.gene) {
  pdf(file.path(outdir, paste0("/GAS_D5GAS_D5_PS_End_NM_EM_AM_ECT_ExM_Dim15_tsne/gene_expression_level_of_", gene, ".pdf")),
      height = 5, width = 6)
  print(FeaturePlot(gasd5.sub3, features = gene, reduction = "tsne", slot = "data",
                    cols = c("#ebebeb", "#7d3c98"), pt.size = 1.5))
  dev.off()
}
for (type in names(neuron.genes)) {
  pdf(file.path(outdir, paste0("/GAS_D5GAS_D5_PS_End_NM_EM_AM_ECT_ExM_Dim15_tsne/Violin_plot_to_show_gene_expression_level_of_", type, ".pdf")),
      height = 10, width = 10)
  print(VlnPlot(gasd5.sub3, features = neuron.genes[[type]],
                group.by = "CellType", stack = T, flip = T, slot = "data"))
  dev.off()
}
# "EPI","PS","NM","End" --> ["#b0c6e4", "#eb8838", "#b7be70", "#cb7eaf"]
set.seed(123)
gasd5.sub4 <- subset(sr.pdt$gasd5.to.hs, CellType %in% c("EPI","PS","End","NM"))
gasd5.sub4 <- NormalizeData(gasd5.sub4) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
gasd5.sub4 <- ScaleData(gasd5.sub4, verbose = FALSE, features = rownames(gasd5.sub4))
gasd5.sub4 <- RunPCA(gasd5.sub4, npcs = 30, verbose = FALSE, features = VariableFeatures(gasd5.sub4))
ElbowPlot(gasd5.sub4, ndims = 30)
dim.n <- 10
gasd5.sub4 <- FindNeighbors(gasd5.sub4, reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 2)
gasd5.sub4 <- RunUMAP(gasd5.sub4, reduction = "pca", dims = 1:dim.n) %>% RunTSNE(reduction = "pca", dims = 1:dim.n)
DimPlot(gasd5.sub4, reduction = "tsne", group.by = "CellType", cols = pd.col.fix)
SeuratToVelocyto(sr.ob = gasd5.sub4,
                 reduction = "tsne", out.dir = "graphs/velocyto/GAS_D5",
                 sample.name = "GAS_D5_EPI_PS_End_NM_Dim10", str.in.barcode = "_.*")
for (gene in pd.gene) {
  pdf(file.path(outdir, paste0("/GAS_D5GAS_D5_EPI_PS_End_NM_Dim10_tsne/gene_expression_level_of_", gene, ".pdf")),
      height = 5, width = 6)
  print(FeaturePlot(gasd5.sub4, features = gene, reduction = "tsne", slot = "data",
                    cols = c("#ebebeb", "#7d3c98"), pt.size = 1.5))
  dev.off()
}
for (type in names(neuron.genes)) {
  pdf(file.path(outdir, paste0("/GAS_D5GAS_D5_EPI_PS_End_NM_Dim10_tsne/Violin_plot_to_show_gene_expression_level_of_", type, ".pdf")),
      height = 10, width = 10)
  print(VlnPlot(gasd5.sub4, features = neuron.genes[[type]],
                group.by = "CellType", stack = T, flip = T, slot = "data"))
  dev.off()
}
# "PS","NM","END": --> ["#eb8838", "#b7be70", "#cb7eaf"]
["#54beca", "#b0c6e4", "#cb7eaf", "#b7be70", "#845a4d", "#eb8838", "#4f996e", "#3b77b3", "#f2adad", "#6d51a0", "#c33b30"]
set.seed(123)
gasd5.sub5 <- subset(sr.pdt$gasd5.to.hs, CellType %in% c("PS","End","NM"))
gasd5.sub5 <- NormalizeData(gasd5.sub5) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
gasd5.sub5 <- ScaleData(gasd5.sub5, verbose = FALSE, features = rownames(gasd5.sub5))
gasd5.sub5 <- RunPCA(gasd5.sub5, npcs = 30, verbose = FALSE, features = VariableFeatures(gasd5.sub5))
ElbowPlot(gasd5.sub5, ndims = 30)
dim.n <- 10
gasd5.sub5 <- FindNeighbors(gasd5.sub5, reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 2)
gasd5.sub5 <- RunUMAP(gasd5.sub5, reduction = "pca", dims = 1:dim.n) %>% RunTSNE(reduction = "pca", dims = 1:dim.n)
DimPlot(gasd5.sub5, reduction = "tsne", group.by = "CellType", cols = pd.col.fix)
SeuratToVelocyto(sr.ob = gasd5.sub5,
                 reduction = "tsne", out.dir = "graphs/velocyto/GAS_D5",
                 sample.name = "GAS_D5_PS_End_NM_Dim10", str.in.barcode = "_.*")
for (gene in pd.gene) {
  pdf(file.path(outdir, paste0("/GAS_D5GAS_D5_PS_End_NM_Dim10_tsne/gene_expression_level_of_", gene, ".pdf")),
      height = 5, width = 6)
  print(FeaturePlot(gasd5.sub5, features = gene, reduction = "tsne", slot = "data",
                    cols = c("#ebebeb", "#7d3c98"), pt.size = 1.5))
  dev.off()
}
for (type in names(neuron.genes)) {
  pdf(file.path(outdir, paste0("/GAS_D5GAS_D5_PS_End_NM_Dim10_tsne//Violin_plot_to_show_gene_expression_level_of_", type, ".pdf")),
      height = 10, width = 10)
  print(VlnPlot(gasd5.sub5, features = neuron.genes[[type]],
                group.by = "CellType", stack = T, flip = T, slot = "data"))
  dev.off()
}
# "EPI","PS","End","NM","EM","AM": ["#3b77b3","#4f996e","#b0c6e4","#eb8838","#b7be70","#cb7eaf"]
set.seed(123)
gasd5.sub6 <- subset(sr.pdt$gasd5.to.hs, CellType %in% c("EPI","PS","End","NM","EM","AM"))
gasd5.sub6 <- NormalizeData(gasd5.sub6) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
gasd5.sub6 <- ScaleData(gasd5.sub6, verbose = FALSE, features = rownames(gasd5.sub6))
gasd5.sub6 <- RunPCA(gasd5.sub6, npcs = 30, verbose = FALSE, features = VariableFeatures(gasd5.sub6))
ElbowPlot(gasd5.sub6, ndims = 30)
dim.n <- 15
gasd5.sub6 <- FindNeighbors(gasd5.sub6, reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 2)
gasd5.sub6 <- RunUMAP(gasd5.sub6, reduction = "pca", dims = 1:dim.n) %>% RunTSNE(reduction = "pca", dims = 1:dim.n)
DimPlot(gasd5.sub6, reduction = "tsne", group.by = "CellType", cols = pd.col.fix)
SeuratToVelocyto(sr.ob = gasd5.sub6,
                 reduction = "tsne", out.dir = "graphs/velocyto/GAS_D5",
                 sample.name = "GAS_D5_EPI_PS_End_NM_EM_AM_Dim15", str.in.barcode = "_.*")
for (gene in pd.gene) {
  pdf(file.path(outdir, paste0("/GAS_D5GAS_D5_EPI_PS_End_NM_EM_AM_Dim15_tsne/gene_expression_level_of_", gene, ".pdf")),
      height = 5, width = 6)
  print(FeaturePlot(gasd5.sub6, features = gene, reduction = "tsne", slot = "data",
                    cols = c("#ebebeb", "#7d3c98"), pt.size = 1.5))
  dev.off()
}
for (type in names(neuron.genes)) {
  pdf(file.path(outdir, paste0("/GAS_D5GAS_D5_EPI_PS_End_NM_EM_AM_Dim15_tsne/Violin_plot_to_show_gene_expression_level_of_", type, ".pdf")),
      height = 10, width = 10)
  print(VlnPlot(gasd5.sub6, features = neuron.genes[[type]],
                group.by = "CellType", stack = T, flip = T, slot = "data"))
  dev.off()
}


### >>> 2. Human CS7
# all cells --> ["#3b77b3", "#845a4d", "#54beca", "#4f996e", "#b0c6e4", "#eb8838", "#c33b30", "#f2adad", "#6d51a0", "#b7be70", "#cb7eaf"]
set.seed(123)
table(hs.cs7$CellType)
hs.cs7.sub0 <- subset(hs.cs7, CellType %in% c("ECT","EPI","PS","NM","AxM","End","EM","AM","ExM","HEP","Ery"))
hs.cs7.sub0 <- NormalizeData(hs.cs7.sub0) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
hs.cs7.sub0 <- ScaleData(hs.cs7.sub0, verbose = FALSE, features = rownames(hs.cs7.sub0))
hs.cs7.sub0 <- RunPCA(hs.cs7.sub0, npcs = 30, verbose = FALSE, features = VariableFeatures(hs.cs7.sub0))
ElbowPlot(hs.cs7.sub0, ndims = 30)
dim.n <- 20
hs.cs7.sub0 <- FindNeighbors(hs.cs7.sub0, reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.cs7.sub0 <- RunUMAP(hs.cs7.sub0, reduction = "pca", dims = 1:dim.n) %>% RunTSNE(reduction = "pca", dims = 1:dim.n)
DimPlot(hs.cs7.sub0, reduction = "tsne", group.by = "CellType", cols = pd.col.fix)
SeuratToVelocyto(sr.ob = hs.cs7.sub0,
                 reduction = "tsne", out.dir = "graphs/velocyto/Hs_cs7",
                 sample.name = "all", str.in.barcode = NULL)
for (gene in pd.gene) {
  pdf(file.path(outdir, paste0("/Hs_cs7all_tsne/gene_expression_level_of_", gene, ".pdf")),
      height = 5, width = 6)
  print(FeaturePlot(hs.cs7.sub0, features = gene, reduction = "tsne", slot = "data",
                    cols = c("#ebebeb", "#7d3c98"), pt.size = 1.5))
  dev.off()
}
# "PS","End","NM","EM","AM" --> ["#3b77b3", "#4f996e", "#eb8838", "#b7be70", "#cb7eaf"]
set.seed(123)
hs.cs7.sub <- subset(hs.cs7, CellType %in% c("PS","End","NM","EM","AM"))
hs.cs7.sub <- NormalizeData(hs.cs7.sub) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
hs.cs7.sub <- ScaleData(hs.cs7.sub, verbose = FALSE, features = rownames(hs.cs7.sub))
hs.cs7.sub <- RunPCA(hs.cs7.sub, npcs = 30, verbose = FALSE, features = VariableFeatures(hs.cs7.sub))
ElbowPlot(hs.cs7.sub, ndims = 30)
dim.n <- 15
hs.cs7.sub <- FindNeighbors(hs.cs7.sub, reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.cs7.sub <- RunUMAP(hs.cs7.sub, reduction = "pca", dims = 1:dim.n) %>% RunTSNE(reduction = "pca", dims = 1:dim.n)
DimPlot(hs.cs7.sub, reduction = "tsne", group.by = "CellType", cols = pd.col.fix)
SeuratToVelocyto(sr.ob = hs.cs7.sub,
                 reduction = "tsne", out.dir = "graphs/velocyto/Hs_cs7",
                 sample.name = "PS_End_NM_EM_AM_Dim15", str.in.barcode = NULL)
for (gene in pd.gene) {
  pdf(file.path(outdir, paste0("/Hs_cs7PS_End_NM_EM_AM_Dim15_tsne/gene_expression_level_of_", gene, ".pdf")),
      height = 5, width = 6)
  print(FeaturePlot(hs.cs7.sub, features = gene, reduction = "tsne", slot = "data",
                    cols = c("#ebebeb", "#7d3c98"), pt.size = 1.5))
  dev.off()
}
# "PS","End","NM","EM","AM","ECT" --> ["#3b77b3", "#54beca", "#4f996e", "#eb8838", "#b7be70", "#cb7eaf"]
set.seed(123)
hs.cs7.sub2 <- subset(hs.cs7, CellType %in% c("PS","End","NM","EM","AM","ECT"))
hs.cs7.sub2 <- NormalizeData(hs.cs7.sub2) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
hs.cs7.sub2 <- ScaleData(hs.cs7.sub2, verbose = FALSE, features = rownames(hs.cs7.sub2))
hs.cs7.sub2 <- RunPCA(hs.cs7.sub2, npcs = 30, verbose = FALSE, features = VariableFeatures(hs.cs7.sub2))
ElbowPlot(hs.cs7.sub2, ndims = 30)
dim.n <- 15
hs.cs7.sub2 <- FindNeighbors(hs.cs7.sub2, reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.cs7.sub2 <- RunUMAP(hs.cs7.sub2, reduction = "pca", dims = 1:dim.n) %>% RunTSNE(reduction = "pca", dims = 1:dim.n)
DimPlot(hs.cs7.sub2, reduction = "tsne", group.by = "CellType", cols = pd.col.fix)
SeuratToVelocyto(sr.ob = hs.cs7.sub2,
                 reduction = "tsne", out.dir = "graphs/velocyto/Hs_cs7",
                 sample.name = "PS_End_NM_EM_AM_ECT_Dim15", str.in.barcode = NULL)
for (gene in pd.gene) {
  pdf(file.path(outdir, paste0("/Hs_cs7PS_End_NM_EM_AM_ECT_Dim15_tsne/gene_expression_level_of_", gene, ".pdf")),
      height = 5, width = 6)
  print(FeaturePlot(hs.cs7.sub2, features = gene, reduction = "tsne", slot = "data",
                    cols = c("#ebebeb", "#7d3c98"), pt.size = 1.5))
  dev.off()
}
# "PS","End","NM","EM","AM","ECT","ExM" -->  ["#3b77b3", "#54beca", "#4f996e", "#eb8838", "#f2adad", "#b7be70", "#cb7eaf"]
set.seed(123)
hs.cs7.sub3 <- subset(hs.cs7, CellType %in% c("PS","End","NM","EM","AM","ECT","ExM"))
hs.cs7.sub3 <- NormalizeData(hs.cs7.sub3) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
hs.cs7.sub3 <- ScaleData(hs.cs7.sub3, verbose = FALSE, features = rownames(hs.cs7.sub3))
hs.cs7.sub3 <- RunPCA(hs.cs7.sub3, npcs = 30, verbose = FALSE, features = VariableFeatures(hs.cs7.sub3))
ElbowPlot(hs.cs7.sub3, ndims = 30)
dim.n <- 15
hs.cs7.sub3 <- FindNeighbors(hs.cs7.sub3, reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.cs7.sub3 <- RunUMAP(hs.cs7.sub3, reduction = "pca", dims = 1:dim.n) %>% RunTSNE(reduction = "pca", dims = 1:dim.n)
DimPlot(hs.cs7.sub3, reduction = "tsne", group.by = "CellType", cols = pd.col.fix)
SeuratToVelocyto(sr.ob = hs.cs7.sub3,
                 reduction = "tsne", out.dir = "graphs/velocyto/Hs_cs7",
                 sample.name = "PS_End_NM_EM_AM_ECT_ExM_Dim15", str.in.barcode = NULL)
for (gene in pd.gene) {
  pdf(file.path(outdir, paste0("/Hs_cs7PS_End_NM_EM_AM_ECT_ExM_Dim15_tsne/gene_expression_level_of_", gene, ".pdf")),
      height = 5, width = 6)
  print(FeaturePlot(hs.cs7.sub3, features = gene, reduction = "tsne", slot = "data",
                    cols = c("#ebebeb", "#7d3c98"), pt.size = 1.5))
  dev.off()
}
# "EPI","PS","NM","End" --> ["#b0c6e4", "#eb8838", "#b7be70", "#cb7eaf"]
set.seed(123)
hs.cs7.sub4 <- subset(hs.cs7, CellType %in% c("EPI","PS","End","NM"))
hs.cs7.sub4 <- NormalizeData(hs.cs7.sub4) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
hs.cs7.sub4 <- ScaleData(hs.cs7.sub4, verbose = FALSE, features = rownames(hs.cs7.sub4))
hs.cs7.sub4 <- RunPCA(hs.cs7.sub4, npcs = 30, verbose = FALSE, features = VariableFeatures(hs.cs7.sub4))
ElbowPlot(hs.cs7.sub4, ndims = 30)
dim.n <- 10
hs.cs7.sub4 <- FindNeighbors(hs.cs7.sub4, reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.cs7.sub4 <- RunUMAP(hs.cs7.sub4, reduction = "pca", dims = 1:dim.n) %>% RunTSNE(reduction = "pca", dims = 1:dim.n)
DimPlot(hs.cs7.sub4, reduction = "tsne", group.by = "CellType", cols = pd.col.fix)
SeuratToVelocyto(sr.ob = hs.cs7.sub4,
                 reduction = "tsne", out.dir = "graphs/velocyto/Hs_cs7",
                 sample.name = "EPI_PS_End_NM_Dim10", str.in.barcode = NULL)
for (gene in pd.gene) {
  pdf(file.path(outdir, paste0("/Hs_cs7EPI_PS_End_NM_Dim10_tsne/gene_expression_level_of_", gene, ".pdf")),
      height = 5, width = 6)
  print(FeaturePlot(hs.cs7.sub4, features = gene, reduction = "tsne", slot = "data",
                    cols = c("#ebebeb", "#7d3c98"), pt.size = 1.5))
  dev.off()
}
# "PS","NM","END": --> ["#eb8838", "#b7be70", "#cb7eaf"]
set.seed(123)
hs.cs7.sub5 <- subset(hs.cs7, CellType %in% c("PS","End","NM"))
hs.cs7.sub5 <- NormalizeData(hs.cs7.sub5) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
hs.cs7.sub5 <- ScaleData(hs.cs7.sub5, verbose = FALSE, features = rownames(hs.cs7.sub5))
hs.cs7.sub5 <- RunPCA(hs.cs7.sub5, npcs = 30, verbose = FALSE, features = VariableFeatures(hs.cs7.sub5))
ElbowPlot(hs.cs7.sub5, ndims = 30)
dim.n <- 10
hs.cs7.sub5 <- FindNeighbors(hs.cs7.sub5, reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.cs7.sub5 <- RunUMAP(hs.cs7.sub5, reduction = "pca", dims = 1:dim.n) %>% RunTSNE(reduction = "pca", dims = 1:dim.n)
DimPlot(hs.cs7.sub5, reduction = "tsne", group.by = "CellType", cols = pd.col.fix)
SeuratToVelocyto(sr.ob = hs.cs7.sub5,
                 reduction = "tsne", out.dir = "graphs/velocyto/Hs_cs7",
                 sample.name = "PS_End_NM_Dim10", str.in.barcode = NULL)
for (gene in pd.gene) {
  pdf(file.path(outdir, paste0("/Hs_cs7PS_End_NM_Dim10_tsne/gene_expression_level_of_", gene, ".pdf")),
      height = 5, width = 6)
  print(FeaturePlot(hs.cs7.sub5, features = gene, reduction = "tsne", slot = "data",
                    cols = c("#ebebeb", "#7d3c98"), pt.size = 1.5))
  dev.off()
}
# "EPI","PS","End","NM","EM","AM": ["#3b77b3","#845a4d","#4f996e","#b0c6e4","#eb8838","#b7be70","#cb7eaf"]
set.seed(123)
hs.cs7.sub6 <- subset(hs.cs7, CellType %in% c("EPI","PS","End","NM","EM","AM","AxM"))
hs.cs7.sub6 <- NormalizeData(hs.cs7.sub6) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
hs.cs7.sub6 <- ScaleData(hs.cs7.sub6, verbose = FALSE, features = rownames(hs.cs7.sub6))
hs.cs7.sub6 <- RunPCA(hs.cs7.sub6, npcs = 30, verbose = FALSE, features = VariableFeatures(hs.cs7.sub6))
ElbowPlot(hs.cs7.sub6, ndims = 30)
dim.n <- 15
hs.cs7.sub6 <- FindNeighbors(hs.cs7.sub6, reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.cs7.sub6 <- RunUMAP(hs.cs7.sub6, reduction = "pca", dims = 1:dim.n) %>% RunTSNE(reduction = "pca", dims = 1:dim.n)
DimPlot(hs.cs7.sub6, reduction = "tsne", group.by = "CellType", cols = pd.col.fix)
SeuratToVelocyto(sr.ob = hs.cs7.sub6,
                 reduction = "tsne", out.dir = "graphs/velocyto/Hs_cs7",
                 sample.name = "EPI_PS_End_NM_EM_AM_AxM_Dim15", str.in.barcode = NULL)
for (gene in pd.gene) {
  pdf(file.path(outdir, paste0("/Hs_cs7EPI_PS_End_NM_EM_AM_AxM_Dim15_tsne/gene_expression_level_of_", gene, ".pdf")),
      height = 5, width = 6)
  print(FeaturePlot(hs.cs7.sub6, features = gene, reduction = "tsne", slot = "data",
                    cols = c("#ebebeb", "#7d3c98"), pt.size = 1.5))
  dev.off()
}
# "PS","End","NM","EM","AM","AxM" --> ["#3b77b3", "#845a4d", "#4f996e", "#eb8838", "#b7be70", "#cb7eaf"]
set.seed(123)
hs.cs7.sub7 <- subset(hs.cs7, CellType %in% c("PS","End","NM","EM","AM","AxM"))
hs.cs7.sub7 <- NormalizeData(hs.cs7.sub7) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
hs.cs7.sub7 <- ScaleData(hs.cs7.sub7, verbose = FALSE, features = rownames(hs.cs7.sub7))
hs.cs7.sub7 <- RunPCA(hs.cs7.sub7, npcs = 30, verbose = FALSE, features = VariableFeatures(hs.cs7.sub7))
ElbowPlot(hs.cs7.sub7, ndims = 30)
dim.n <- 15
hs.cs7.sub7 <- FindNeighbors(hs.cs7.sub7, reduction = "pca", dims = 1:dim.n) %>% FindClusters(resolution = 2)
hs.cs7.sub7 <- RunUMAP(hs.cs7.sub7, reduction = "pca", dims = 1:dim.n) %>% RunTSNE(reduction = "pca", dims = 1:dim.n)
DimPlot(hs.cs7.sub7, reduction = "tsne", group.by = "CellType", cols = pd.col.fix)
SeuratToVelocyto(sr.ob = hs.cs7.sub7,
                 reduction = "tsne", out.dir = "graphs/velocyto/Hs_cs7",
                 sample.name = "PS_End_NM_EM_AM_AxM_Dim15", str.in.barcode = NULL)
for (gene in pd.gene) {
  pdf(file.path(outdir, paste0("/Hs_cs7PS_End_NM_EM_AM_AxM_Dim15_tsne/gene_expression_level_of_", gene, ".pdf")),
      height = 5, width = 6)
  print(FeaturePlot(hs.cs7.sub7, features = gene, reduction = "tsne", slot = "data",
                    cols = c("#ebebeb", "#7d3c98"), pt.size = 1.5))
  dev.off()
}



### =================================
### 7th part: Gastruloids development ----
### =================================


### >>> 1. Setting output directory
outdir <- file.path(getwd(), "graphs/gas_dev")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}


### >>> 2. Define cell types before GAS D5
gas.before <- Reduce(function(x, y) merge(x, y, add.cell.ids = c(x@project.name, y@project.name)), sr.list[2:3])
gas.before <- subset(gas.before, Sample_Name != "GAS_D5")
gas.before <- subset(gas.before, Sample_Name != "GAS_D2")
gas.before <- NormalizeData(gas.before) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(gas.before))
gas.before <- RunPCA(gas.before, features = VariableFeatures(object = gas.before))
ElbowPlot(gas.before)
dim.n <- 10
gas.before <- RunUMAP(gas.before, reduction = "pca", dims = 1:dim.n) %>%
  RunTSNE(dims = 1:dim.n)
DimPlot(gas.before, reduction = "umap", group.by = c("Sample_Name"), pt.size = 0.5)
gas.before@meta.data <- gas.before@meta.data %>%
  mutate(Batch = case_when(Sample_Name %in% c("EPI_D0", "GAS_D1", "GAS_D3") ~ "Batch1",
                           Sample_Name %in% c("GAS_D5", "ORG_D2") ~ "Batch2"))
table(gas.before$Batch)
set.seed(123)
gas.before <- RunHarmony(gas.before, group.by.vars = "Batch")
gas.before <- FindNeighbors(gas.before, reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 0.75)
gas.before <- RunUMAP(gas.before, reduction = "harmony", dims = 1:dim.n, return.model = T)
pdf(file.path(outdir, "Scatter_plot_to_show_clustering_of_gastruloid_development_before_d5.pdf"), height = 5, width = 12)
p1 <- DimPlot(gas.before, reduction = "umap", group.by = c("Sample_Name"), pt.size = 0.5, label = T) +
  scale_colour_brewer(palette = "Set1")
p2 <- DimPlot(gas.before, reduction = "umap", group.by = c("seurat_clusters"), pt.size = 0.5, label = T)
p1 + p2
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_marker_expression_of_gastruloid_development_before_d5.pdf"), height = 20, width = 25)
FeaturePlot(gas.before, features = c("SOX2", "POU5F1", "OTX2", "CDH1",
                                     "TBXT", "MIXL1", "SP5", "FST",
                                     "SOX17", "FOXA2", "GSC", "CST1",
                                     "TBX6", "MESP1", "EOMES",
                                     "LHX1", "IRX3", "GATA6",
                                     "FOXF1", "HAND1", "HAND2",
                                     "POSTN", "LUM", "NID2",
                                     "MEF2C", "PECAM1", "CDH5", "TEK",
                                     "DLX5", "TFAP2A", "TFAP2C", "GATA3","GFP"), ncol = 6, cols = c("#ebebeb", "#7d3c98"))
dev.off()
# rename cells
new.cluster.ids <- c("EPILCs","EPILCs","EndLCs","EPILCs","EPILCs","EPILCs","MXCs2",
                     "EPILCs","MXCs1","EPILCs","MXCs1","EPILCs","EndLCs","EndLCs")
names(new.cluster.ids) <- levels(gas.before)
gas.before <- RenameIdents(gas.before, new.cluster.ids)
pdf(file.path(outdir, "Scatter_plot_to_show_annotated_clustering_of_gastruloid_development_before_d5.pdf"),
    height = 5, width = 12)
pd.col.fix3 <- c(ORG_D2 = "#4CAF49", EPI_D0 = "#B91122", GAS_D1 = "#9ECAE1", GAS_D3 = "#EB8838", GAS_D5 = "#08519C")
p1 <- DimPlot(gas.before, reduction = "umap", group.by = c("Sample_Name"), pt.size = 0.5, label = T) +
  scale_colour_manual(values = pd.col.fix3)
pd.col.fix2 <- c(EPILCs = "#B91122", EndLCs = "#A65E25", MXCs1 = "#cb7eaf", MXCs2 = "#b7be70")
p2 <- DimPlot(gas.before, reduction = "umap", group.by = c("ident"), pt.size = 0.5, label = T) +
  scale_colour_manual(values = pd.col.fix2)
p1 + p2
dev.off()
gas.before$CellType <- as.character(Idents(gas.before))
# plot expression
sr.out <- "/home/yhw/bioinfo/project-wb/hs_gastruloids/analysis/graphs/gas_dev/GAS_before_D5_marker_gene"
dir.create(sr.out, recursive = T)
marker.gene <- list(EPI = c("CDH1", "POU5F1", "SOX2", "NANOG", "LIN28AP1", "MYC", "KLF4", "TDGF1", "GDF3", "GFP"),
                    Organizer = c("FST", "OTX2", "GSC", "EOMES", "FOXH1", "CER1", "MIXL1", "DKK1", "NODAL",
                                  "FOXA2", "LEFTY1", "LEFTY2", "TBXT", "SOX17", "SP5", "NOG", "CHRD"),
                    Embryonic = c("TBXT", "MESP1", "PDGFRA", "LHX1", "IRX3", "GATA6", "MYL7", "MSX1",
                                  "MSX2", "FOXF1", "HAND1", "HAND2", "BMP4", "MAB21L2", "CST1", "MEF2C",
                                  "PECAM1", "CDH5", "TEK", "TFAP2A", "TFAP2C", "GATA3", "DLX5", "LHX5",
                                  "PAX6", "SOX3", "PAX3", "PAX7"),
                    Extraembryonic = c("LUM", "NID2", "ANXA1", "POSTN", "GATA2", "KRT7", "VIM", "GATA4",
                                       "DCN", "TP63", "KRT18", "NR2F2", "ISL1", "HEY1", "CDH10", "CTSV",
                                       "ARID5B", "PLAGL1", "CREB3L1", "HOXA10", "HOXA11", "HOXA9",
                                       "HOXA13", "WNT6", "GABRP", "HEY1", "BST2"))
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("GAS_before_D5_marker_gene_expression_of_", marker.gene[[type]][g], "_featureplot.pdf")),
        height = 3.5, width = 4)
    print(FeaturePlot(gas.before, features = marker.gene[[type]][g],
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 0.5, slot = "data"))
    dev.off()
  }; rm(g)
}
# plot cell ratio
pdf(file.path(outdir, "Bar_plot_to_show_cell_type_ratio_after_annotating_clustering_gastruloid_development_before_d5.pdf"),
    height = 5, width = 7.5)
gas.before@meta.data[,c("Sample_Name", "CellType")] %>%
  group_by(Sample_Name, CellType) %>%
  summarise(CellNum = length(CellType)) %>%
  mutate(CellType = factor(CellType, levels = c("EPILCs", "EndLCs", "MXCs1", "MXCs2")),
         Sample_Name = factor(Sample_Name, levels = c("EPI_D0", "ORG_D2", "GAS_D1", "GAS_D3"))) %>%
  ggplot(aes(x = Sample_Name, y = CellNum, fill = CellType)) +
  geom_bar(position = "stack", stat = "identity", width = 0.75) +
  scale_fill_manual(values = pd.col.fix2) +
  theme_bw() +
  labs(x = "Sample Group", y = "Number of cells", fill = "Cell Type") +
  theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 14, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 14, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        panel.grid = element_blank())
dev.off()


### >>> 3. Sub-clustering
# - MXCs1
gas.before.sub1 <- subset(gas.before, CellType %in% c("MXCs1"))
gas.before.sub1 <- NormalizeData(gas.before.sub1) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(gas.before.sub1))
gas.before.sub1 <- RunPCA(gas.before.sub1, features = VariableFeatures(object = gas.before.sub1))
ElbowPlot(gas.before.sub1)
dim.n <- 7
gas.before.sub1 <- RunUMAP(gas.before.sub1, reduction = "pca", dims = 1:dim.n) %>%
  RunTSNE(dims = 1:dim.n)
DimPlot(gas.before.sub1, reduction = "umap", group.by = c("Sample_Name"), pt.size = 0.5)
table(gas.before.sub1$Batch)
set.seed(123)
gas.before.sub1 <- RunHarmony(gas.before.sub1, group.by.vars = "Batch", plot_convergence = FALSE, project.dim = F)
gas.before.sub1 <- FindNeighbors(gas.before.sub1, reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 0.75)
gas.before.sub1 <- RunUMAP(gas.before.sub1, reduction = "harmony", dims = 1:dim.n, return.model = T)
pdf(file.path(outdir, "Scatter_plot_to_show_clustering_of_gastruloid_development_before_d5_sub1.pdf"), height = 4.5, width = 15)
p1 <- DimPlot(gas.before.sub1, reduction = "umap", group.by = c("Sample_Name"), pt.size = 1, label = T) +
  scale_colour_manual(values = pd.col.fix3)
p2 <- DimPlot(gas.before.sub1, reduction = "umap", group.by = c("seurat_clusters"), pt.size = 1, label = T) +
  scale_colour_brewer(palette = "Paired")
pd.col.fix2 <- c(EPILCs = "#B91122", EndLCs = "#A65E25", MXCs1 = "#cb7eaf", MXCs2 = "#b7be70")
p3 <- DimPlot(gas.before.sub1, reduction = "umap", group.by = c("CellType"), pt.size = 1, label = T) +
  scale_colour_manual(values = pd.col.fix2)
p1 + p2 + p3
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_marker_expression_of_gastruloid_development_before_d5_sub1.pdf"), height = 20, width = 25)
FeaturePlot(gas.before.sub1, features = c("SOX2", "POU5F1", "OTX2", "CDH1",
                                     "TBXT", "MIXL1", "SP5", "FST",
                                     "SOX17", "FOXA2", "GSC", "CST1",
                                     "TBX6", "MESP1", "EOMES",
                                     "LHX1", "IRX3", "GATA6",
                                     "FOXF1", "HAND1", "HAND2",
                                     "POSTN", "LUM", "NID2",
                                     "MEF2C", "PECAM1", "CDH5", "TEK",
                                     "DLX5", "TFAP2A", "TFAP2C", "GATA3","GFP"),
            ncol = 6, cols = c("#ebebeb", "#7d3c98"))
dev.off()
# rename cells
new.cluster.ids <- c("EPILCs","EPILCs","EPILCs","ALCs","PSLCs","EPILCs","PSLCs")
names(new.cluster.ids) <- levels(gas.before.sub1)
gas.before.sub1 <- RenameIdents(gas.before.sub1, new.cluster.ids)
pdf(file.path(outdir, "Scatter_plot_to_show_annotated_clustering_of_gastruloid_development_before_d5_sub1.pdf"),
    height = 5, width = 12)
pd.col.fix3 <- c(ORG_D2 = "#4CAF49", EPI_D0 = "#B91122", GAS_D1 = "#9ECAE1", GAS_D3 = "#EB8838", GAS_D5 = "#08519C")
p1 <- DimPlot(gas.before.sub1, reduction = "umap", group.by = c("Sample_Name"), pt.size = 2, label = T) +
  scale_colour_manual(values = pd.col.fix3)
pd.col.fix2 <- c(EPILCs = "#B91122", PSLCs = "#A12F77", ALCs = "#05B5FC")
p2 <- DimPlot(gas.before.sub1, reduction = "umap", group.by = c("ident"), pt.size = 2, label = T) +
  scale_colour_manual(values = pd.col.fix2)
p1 + p2
dev.off()
gas.before.sub1$CellType.new <- as.character(Idents(gas.before.sub1))
# plot expression
sr.out <- "/home/yhw/bioinfo/project-wb/hs_gastruloids/analysis/graphs/gas_dev/GAS_before_D5_MSX1_marker_gene"
dir.create(sr.out, recursive = T)
marker.gene <- list(EPI = c("CDH1", "POU5F1", "SOX2", "NANOG", "LIN28AP1", "MYC", "KLF4", "TDGF1", "GDF3", "GFP"),
                    Organizer = c("FST", "OTX2", "GSC", "EOMES", "FOXH1", "CER1", "MIXL1", "DKK1", "NODAL",
                                  "FOXA2", "LEFTY1", "LEFTY2", "TBXT", "SOX17", "SP5", "NOG", "CHRD"),
                    Embryonic = c("TBXT", "MESP1", "PDGFRA", "LHX1", "IRX3", "GATA6", "MYL7", "MSX1",
                                  "MSX2", "FOXF1", "HAND1", "HAND2", "BMP4", "MAB21L2", "CST1", "MEF2C",
                                  "PECAM1", "CDH5", "TEK", "TFAP2A", "TFAP2C", "GATA3", "DLX5", "LHX5",
                                  "PAX6", "SOX3", "PAX3", "PAX7"),
                    Extraembryonic = c("LUM", "NID2", "ANXA1", "POSTN", "GATA2", "KRT7", "VIM", "GATA4",
                                       "DCN", "TP63", "KRT18", "NR2F2", "ISL1", "HEY1", "CDH10", "CTSV",
                                       "ARID5B", "PLAGL1", "CREB3L1", "HOXA10", "HOXA11", "HOXA9",
                                       "HOXA13", "WNT6", "GABRP", "HEY1", "BST2"))
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("MSX1_marker_gene_expression_of_", marker.gene[[type]][g], "_featureplot.pdf")),
        height = 3.5, width = 4)
    print(FeaturePlot(gas.before.sub1, features = marker.gene[[type]][g],
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 0.5, slot = "data"))
    dev.off()
  }; rm(g)
}
# plot cell ratio
pdf(file.path(outdir, "Bar_plot_to_show_cell_type_ratio_after_annotating_clustering_gastruloid_development_before_d5_sub1.pdf"),
    height = 5, width = 7.5)
gas.before.sub1@meta.data[,c("Sample_Name", "CellType.new")] %>%
  group_by(Sample_Name, CellType.new) %>%
  summarise(CellNum = length(CellType.new)) %>%
  mutate(CellType.new = factor(CellType.new, levels = c("EPILCs", "PSLCs", "ALCs")),
         Sample_Name = factor(Sample_Name, levels = c("EPI_D0", "ORG_D2", "GAS_D1", "GAS_D3"))) %>%
  ggplot(aes(x = Sample_Name, y = CellNum, fill = CellType.new)) +
  geom_bar(position = "stack", stat = "identity", width = 0.75) +
  scale_fill_manual(values = pd.col.fix2) +
  theme_bw() +
  labs(x = "Sample Group", y = "Number of cells", fill = "Cell Type") +
  theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 14, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 14, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        panel.grid = element_blank())
dev.off()
# - MXCs2
gas.before.sub2 <- subset(gas.before, CellType %in% c("MXCs2"))
gas.before.sub2 <- NormalizeData(gas.before.sub2) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(gas.before.sub2))
gas.before.sub2 <- RunPCA(gas.before.sub2, features = VariableFeatures(object = gas.before.sub2))
ElbowPlot(gas.before.sub2)
dim.n <- 6
gas.before.sub2 <- RunUMAP(gas.before.sub2, reduction = "pca", dims = 1:dim.n) %>%
  RunTSNE(dims = 1:dim.n)
DimPlot(gas.before.sub2, reduction = "umap", group.by = c("Sample_Name"), pt.size = 0.5)
table(gas.before.sub2$Batch, gas.before.sub2$Sample_Name)
set.seed(123)
gas.before.sub2 <- RunHarmony(gas.before.sub2, group.by.vars = "Batch", plot_convergence = FALSE, project.dim = F)
gas.before.sub2 <- FindNeighbors(gas.before.sub2, reduction = "harmony", dims = 1:dim.n) %>% FindClusters(resolution = 0.75)
gas.before.sub2 <- RunUMAP(gas.before.sub2, reduction = "harmony", dims = 1:dim.n, return.model = T)
pdf(file.path(outdir, "Scatter_plot_to_show_clustering_of_gastruloid_development_before_d5_sub2.pdf"), height = 4.5, width = 15)
p1 <- DimPlot(gas.before.sub2, reduction = "umap", group.by = c("Sample_Name"), pt.size = 1, label = T) +
  scale_colour_manual(values = pd.col.fix3)
p2 <- DimPlot(gas.before.sub2, reduction = "umap", group.by = c("seurat_clusters"), pt.size = 1, label = T) +
  scale_colour_brewer(palette = "Paired")
pd.col.fix2 <- c(EPILCs = "#B91122", EndLCs = "#A65E25", MXCs1 = "#cb7eaf", MXCs2 = "#b7be70")
p3 <- DimPlot(gas.before.sub2, reduction = "umap", group.by = c("CellType"), pt.size = 1, label = T) +
  scale_colour_manual(values = pd.col.fix2)
p1 + p2 + p3
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_marker_expression_of_gastruloid_development_before_d5_sub2.pdf"), height = 20, width = 25)
FeaturePlot(gas.before.sub2, features = c("SOX2", "POU5F1", "OTX2", "CDH1",
                                          "TBXT", "MIXL1", "SP5", "FST",
                                          "SOX17", "FOXA2", "GSC", "CST1",
                                          "TBX6", "MESP1", "EOMES",
                                          "LHX1", "IRX3", "GATA6",
                                          "FOXF1", "HAND1", "HAND2",
                                          "POSTN", "LUM", "NID2",
                                          "MEF2C", "PECAM1", "CDH5", "TEK",
                                          "DLX5", "TFAP2A", "TFAP2C", "GATA3","GFP"),
            ncol = 6, cols = c("#ebebeb", "#7d3c98"))
dev.off()
# rename cells
new.cluster.ids <- c("EMLCs","ExMLCs","EMLCs","EMLCs","EMLCs","EMLCs","EMLCs","EMLCs","EMLCs")
names(new.cluster.ids) <- levels(gas.before.sub2)
gas.before.sub2 <- RenameIdents(gas.before.sub2, new.cluster.ids)
pdf(file.path(outdir, "Scatter_plot_to_show_annotated_clustering_of_gastruloid_development_before_d5_sub2.pdf"),
    height = 5, width = 12)
pd.col.fix3 <- c(ORG_D2 = "#4CAF49", EPI_D0 = "#B91122", GAS_D1 = "#9ECAE1", GAS_D3 = "#EB8838", GAS_D5 = "#08519C")
p1 <- DimPlot(gas.before.sub2, reduction = "umap", group.by = c("Sample_Name"), pt.size = 2, label = T) +
  scale_colour_manual(values = pd.col.fix3)
pd.col.fix2 <- c(EMLCs = "#7D861E", ExMLCs = "#CE6666")
p2 <- DimPlot(gas.before.sub2, reduction = "umap", group.by = c("ident"), pt.size = 2, label = T) +
  scale_colour_manual(values = pd.col.fix2)
p1 + p2
dev.off()
gas.before.sub2$CellType.new <- as.character(Idents(gas.before.sub2))
# plot expression
sr.out <- "/home/yhw/bioinfo/project-wb/hs_gastruloids/analysis/graphs/gas_dev/GAS_before_D5_MSX2_marker_gene"
dir.create(sr.out, recursive = T)
marker.gene <- list(EPI = c("CDH1", "POU5F1", "SOX2", "NANOG", "LIN28AP1", "MYC", "KLF4", "TDGF1", "GDF3", "GFP"),
                    Organizer = c("FST", "OTX2", "GSC", "EOMES", "FOXH1", "CER1", "MIXL1", "DKK1", "NODAL",
                                  "FOXA2", "LEFTY1", "LEFTY2", "TBXT", "SOX17", "SP5", "NOG", "CHRD"),
                    Embryonic = c("TBXT", "MESP1", "PDGFRA", "LHX1", "IRX3", "GATA6", "MYL7", "MSX1",
                                  "MSX2", "FOXF1", "HAND1", "HAND2", "BMP4", "MAB21L2", "CST1", "MEF2C",
                                  "PECAM1", "CDH5", "TEK", "TFAP2A", "TFAP2C", "GATA3", "DLX5", "LHX5",
                                  "PAX6", "SOX3", "PAX3", "PAX7"),
                    Extraembryonic = c("LUM", "NID2", "ANXA1", "POSTN", "GATA2", "KRT7", "VIM", "GATA4",
                                       "DCN", "TP63", "KRT18", "NR2F2", "ISL1", "HEY1", "CDH10", "CTSV",
                                       "ARID5B", "PLAGL1", "CREB3L1", "HOXA10", "HOXA11", "HOXA9",
                                       "HOXA13", "WNT6", "GABRP", "HEY1", "BST2"))
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("MSX2_marker_gene_expression_of_", marker.gene[[type]][g], "_featureplot.pdf")),
        height = 3.5, width = 4)
    print(FeaturePlot(gas.before.sub2, features = marker.gene[[type]][g],
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 0.5, slot = "data"))
    dev.off()
  }; rm(g)
}
# plot cell ratio
pdf(file.path(outdir, "Bar_plot_to_show_cell_type_ratio_after_annotating_clustering_gastruloid_development_before_d5_sub2.pdf"),
    height = 5, width = 7.5)
gas.before.sub2@meta.data[,c("Sample_Name", "CellType.new")] %>%
  group_by(Sample_Name, CellType.new) %>%
  summarise(CellNum = length(CellType.new)) %>%
  mutate(CellType.new = factor(CellType.new, levels = c("EMLCs", "ExMLCs")),
         Sample_Name = factor(Sample_Name, levels = c("EPI_D0", "ORG_D2", "GAS_D1", "GAS_D3"))) %>%
  ggplot(aes(x = Sample_Name, y = CellNum, fill = CellType.new)) +
  geom_bar(position = "stack", stat = "identity", width = 0.75) +
  scale_fill_manual(values = pd.col.fix2) +
  theme_bw() +
  labs(x = "Sample Group", y = "Number of cells", fill = "Cell Type") +
  theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 14, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 14, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        panel.grid = element_blank())
dev.off()


### >>> 4. Merge all data (GAS + EPI)
# sub1
tmp1 <- gas.before.sub1
tmp1$CellType <- tmp1$CellType.new
tmp1@meta.data <- tmp1@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA",
                                     "Sample_Tag", "Sample_Name", "Pect.mt",
                                     "GFP_count", "CellType")]
# sub2
tmp2 <- gas.before.sub2
tmp2$CellType <- tmp2$CellType.new
tmp2@meta.data <- tmp2@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA",
                                     "Sample_Tag", "Sample_Name", "Pect.mt",
                                     "GFP_count", "CellType")]
# sub3
tmp3 <- gas.before[, !(gas.before$CellType %in% c("MXCs1", "MXCs2"))]
tmp3@meta.data <- tmp3@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA",
                                     "Sample_Tag", "Sample_Name", "Pect.mt",
                                     "GFP_count", "CellType")]
# sub4
tmp4 <- sr.pdt$gasd5.to.hs
tmp4@meta.data <- tmp4@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA",
                                     "Sample_Tag", "Sample_Name", "Pect.mt",
                                     "GFP_count", "CellType")]
# merge
tmp <- list(tmp1, tmp2, tmp3)
gas.d5before.merge <- Reduce(function(x, y) merge(x, y, add.cell.ids = c(x@project.name,y@project.name)), tmp)
gas.d5before.merge@meta.data %>%
  mutate(Type = case_when(GFP_count > 0 ~ "GFP_pos",
                          TRUE ~ "GFP_neg")) -> gas.d5before.merge@meta.data
table(gas.d5before.merge$Sample_Name)
table(gas.d5before.merge$Type)
table(gas.d5before.merge$CellType)
gas.d5before.merge.col <- c(EMLCs = "#7D861E", EndLCs = "#A65E25", EPILCs = "#B91122",
                            ExMLCs = "#CE6666", PSLCs = "#A12F77", ALCs = "#05B5FC")
gas.d5before.merge <- NormalizeData(gas.d5before.merge) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(gas.d5before.merge))
gas.d5before.merge <- RunPCA(gas.d5before.merge, features = VariableFeatures(object = gas.d5before.merge))
ElbowPlot(gas.d5before.merge)
dim.n <- 6
gas.d5before.merge <- RunUMAP(gas.d5before.merge, reduction = "pca", dims = 1:dim.n) %>%
  RunTSNE(dims = 1:dim.n)
DimPlot(gas.d5before.merge, reduction = "umap", group.by = c("Sample_Name","CellType"), pt.size = 0.5)
gas.d5before.merge@meta.data <- gas.d5before.merge@meta.data %>%
  mutate(Batch = case_when(Sample_Name %in% c("EPI_D0", "GAS_D1", "GAS_D3") ~ "Batch1",
                           Sample_Name %in% c("GAS_D5", "ORG_D2") ~ "Batch2"))
table(gas.d5before.merge$Batch)
set.seed(123)
gas.d5before.merge <- RunHarmony(gas.d5before.merge, group.by.vars = "Batch", plot_convergence = FALSE, project.dim = F)
gas.d5before.merge <- FindNeighbors(gas.d5before.merge, reduction = "harmony", dims = 1:dim.n) %>%
  FindClusters(resolution = 0.75)
gas.d5before.merge <- RunUMAP(gas.d5before.merge, reduction = "harmony", dims = 1:dim.n, return.model = T) %>%
  RunTSNE(reduction = "harmony", dims = 1:dim.n)
pdf(file.path(outdir, "Scatter_plot_to_show_clustering_of_gastruloid_development_before_d5_sub1.pdf"),
    height = 4.5, width = 12)
p1 <- DimPlot(gas.d5before.merge, reduction = "umap", group.by = c("Sample_Name"), pt.size = 1, label = T) +
  scale_colour_manual(values = pd.col.fix3)
p2 <- DimPlot(gas.d5before.merge, reduction = "umap", group.by = c("CellType"), pt.size = 1, label = T) +
  scale_colour_manual(values = gas.d5before.merge.col)
p3 <- DimPlot(gas.d5before.merge, reduction = "umap", group.by = c("seurat_clusters"), pt.size = 1, label = T)
p1 + p2 + p3
dev.off()
rm.cells <- gas.d5before.merge$CellType == "EPILCs" & gas.d5before.merge$seurat_clusters == 7
gas.d5before.merge.sub <- gas.d5before.merge[,!rm.cells]
gas.d5before.merge.sub <- NormalizeData(gas.d5before.merge.sub) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(gas.d5before.merge.sub))
gas.d5before.merge.sub <- RunPCA(gas.d5before.merge.sub, features = VariableFeatures(gas.d5before.merge.sub))
ElbowPlot(gas.d5before.merge.sub)
dim.n <- 6
gas.d5before.merge.sub <- RunUMAP(gas.d5before.merge.sub, reduction = "pca", dims = 1:dim.n) %>%
  RunTSNE(dims = 1:dim.n)
DimPlot(gas.d5before.merge.sub, reduction = "umap", group.by = c("Sample_Name","CellType"), pt.size = 0.5)
set.seed(123)
gas.d5before.merge.sub <- RunHarmony(gas.d5before.merge.sub, group.by.vars = "Batch", plot_convergence = FALSE, project.dim = F)
gas.d5before.merge.sub <- FindNeighbors(gas.d5before.merge.sub, reduction = "harmony", dims = 1:dim.n) %>%
  FindClusters(resolution = 0.75)
gas.d5before.merge.sub <- RunUMAP(gas.d5before.merge.sub, reduction = "harmony", dims = 1:dim.n, return.model = T) %>%
  RunTSNE(reduction = "harmony", dims = 1:dim.n)
pdf(file.path(outdir, "Scatter_plot_to_show_clustering_of_gastruloid_development_before_d5_merged.pdf"),
    height = 4.5, width = 11)
p1 <- DimPlot(gas.d5before.merge.sub, reduction = "umap", group.by = c("Sample_Name"), pt.size = 1, label = T) +
  scale_colour_manual(values = pd.col.fix3)
p2 <- DimPlot(gas.d5before.merge.sub, reduction = "umap", group.by = c("CellType"), pt.size = 1, label = T) +
  scale_colour_manual(values = gas.d5before.merge.col)
p1 + p2
dev.off()
pdf(file.path(outdir, "Bar_plot_to_show_cell_type_ratio_after_annotating_clustering_gastruloid_development_before_d5_merged.pdf"),
    height = 5, width = 7.5)
gas.d5before.merge.sub@meta.data[,c("Sample_Name", "CellType")] %>%
  group_by(Sample_Name, CellType) %>%
  summarise(CellNum = length(CellType)) %>%
  mutate(CellType = factor(CellType, levels = c("EPILCs", "PSLCs", "EMLCs", "EndLCs", "ExMLCs", "ALCs")),
         Sample_Name = factor(Sample_Name, levels = c("EPI_D0", "ORG_D2", "GAS_D1", "GAS_D3"))) %>%
  ggplot(aes(x = Sample_Name, y = CellNum, fill = CellType)) +
  geom_bar(position = "stack", stat = "identity", width = 0.75) +
  scale_fill_manual(values = gas.d5before.merge.col) +
  theme_bw() +
  labs(x = "Sample Group", y = "Number of cells", fill = "Cell Type") +
  theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 14, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 14, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        panel.grid = element_blank())
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_marker_expression_of_gastruloid_development_before_d5_merged.pdf"), 
    height = 20, width = 25)
FeaturePlot(gas.d5before.merge.sub, features = c("SOX2", "POU5F1", "OTX2", "CDH1",
                                                 "TBXT", "MIXL1", "SP5", "FST",
                                                 "SOX17", "FOXA2", "GSC", "CST1",
                                                 "TBX6", "MESP1", "EOMES",
                                                 "LHX1", "IRX3", "GATA6",
                                                 "FOXF1", "HAND1", "HAND2",
                                                 "POSTN", "LUM", "NID2",
                                                 "MEF2C", "PECAM1", "CDH5", "TEK",
                                                 "DLX5", "TFAP2A", "TFAP2C", "GATA3","GFP"),
            ncol = 6, cols = c("#ebebeb", "#7d3c98"))
dev.off()
pd.ge <- c("SOX2", "POU5F1", "OTX2", "CDH1",
           "TBXT", "MIXL1", "SP5", "FST",
           "SOX17", "FOXA2", "GSC", "CST1",
           "TBX6", "MESP1", "EOMES",
           "LHX1", "IRX3", "GATA6",
           "FOXF1", "HAND1", "HAND2",
           "POSTN", "LUM", "NID2",
           "MEF2C", "PECAM1", "CDH5", "TEK",
           "DLX5", "TFAP2A", "TFAP2C", "GATA3","GFP")
for (g in seq(1, length(pd.ge))) {
  pdf(paste0(outdir, "/Scatter_plot_to_show_", pd.ge[g], "_expression_of_gastruloid_development_before_d5_merged_featureplot.pdf"),
      height = 3.5, width = 4)
  print(FeaturePlot(gas.d5before.merge.sub, features = pd.ge[g], cols = c("#ebebeb", "#7d3c98"), pt.size = 1.35, reduction = "umap"))
  dev.off()
}; rm(g)

# merge all subs
tmp <- list(tmp1, tmp2, tmp3, tmp4)
gas.merge <- Reduce(function(x, y) merge(x, y, add.cell.ids = c(x@project.name,y@project.name)), tmp)
gas.merge@meta.data %>%
  mutate(Type = case_when(GFP_count > 0 ~ "GFP_pos",
                          TRUE ~ "GFP_neg")) -> gas.merge@meta.data
table(gas.merge$Sample_Name)
table(gas.merge$Type)
table(gas.merge$CellType)
gas.merge.col <- c(c(EMLCs = "#7D861E", EndLCs = "#A65E25", EPILCs = "#B91122",
                     ExMLCs = "#CE6666", PSLCs = "#A12F77", ALCs = "#05B5FC"),
                   pd.col.fix)
gas.merge <- NormalizeData(gas.merge) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(gas.merge))
gas.merge <- RunPCA(gas.merge, features = VariableFeatures(object = gas.merge))
ElbowPlot(gas.merge)
dim.n <- 10
gas.merge <- RunUMAP(gas.merge, reduction = "pca", dims = 1:dim.n) %>%
  RunTSNE(dims = 1:dim.n)
DimPlot(gas.merge, reduction = "umap", group.by = c("Sample_Name","CellType"), pt.size = 0.5, label = T) +
  scale_color_manual(values = gas.merge.col)
gas.merge@meta.data <- gas.merge@meta.data %>%
  mutate(Batch = case_when(Sample_Name %in% c("EPI_D0", "GAS_D1", "GAS_D3") ~ "Batch1",
                           Sample_Name %in% c("GAS_D5", "ORG_D2") ~ "Batch2"))
table(gas.merge$Batch)
set.seed(123)
gas.merge <- RunHarmony(gas.merge, group.by.vars = "Batch", plot_convergence = FALSE, project.dim = F)
gas.merge <- FindNeighbors(gas.merge, reduction = "harmony", dims = 1:dim.n) %>%
  FindClusters(resolution = 0.75)
gas.merge <- RunUMAP(gas.merge, reduction = "harmony", dims = 1:dim.n, return.model = T) %>%
  RunTSNE(reduction = "harmony", dims = 1:dim.n)
pdf(file.path(outdir, "Scatter_plot_to_show_clustering_of_gastruloid_development_merged.pdf"),
    height = 5, width = 17)
p1 <- DimPlot(gas.merge, reduction = "umap", group.by = c("Sample_Name"), pt.size = 0.5, label = T, repel = T) +
  scale_colour_manual(values = pd.col.fix3)
p2 <- DimPlot(gas.merge, reduction = "umap", group.by = c("CellType"), pt.size = 0.5, label = T, repel = T) +
  scale_colour_manual(values = gas.merge.col)
p3 <- DimPlot(gas.merge, reduction = "umap", group.by = c("seurat_clusters"), pt.size = 0.5, label = T, repel = T)
p1 + p2 + p3
dev.off()

for (i in unique(gas.merge@meta.data$Sample_Name)) {
  Vis.Annotation.Ratio(subset(gas.merge@meta.data, Sample_Name == i), 
                       annotation = c("CellType", "Type"),
                       pd.title = i, pd.col = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279"),
                       pd.height = 10, pd.width = 10, 
                       res.out = file.path(outdir, paste0("Cell_Type_Ratio_", i)))
}
pdf(file.path(outdir, "Bar_plot_to_show_cell_type_ratio_after_annotating_clustering_gastruloid_development_merged.pdf"),
    height = 5, width = 7.5)
gas.merge@meta.data[,c("Sample_Name", "CellType")] %>%
  group_by(Sample_Name, CellType) %>%
  summarise(CellNum = length(CellType)) %>%
  mutate(Sample_Name = factor(Sample_Name, levels = c("EPI_D0", "ORG_D2", "GAS_D1", "GAS_D3", "GAS_D5"))) %>%
  ggplot(aes(x = Sample_Name, y = CellNum, fill = CellType)) +
  geom_bar(position = "stack", stat = "identity", width = 0.75) +
  scale_fill_manual(values = gas.merge.col) +
  theme_bw() +
  labs(x = "Sample Group", y = "Number of cells", fill = "Cell Type") +
  theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 14, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 14, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        panel.grid = element_blank())
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_marker_expression_of_gastruloid_development_merged.pdf"), height = 20, width = 25)
FeaturePlot(gas.merge, features = c("SOX2", "POU5F1", "OTX2", "CDH1",
                                    "TBXT", "MIXL1", "SP5", "FST",
                                    "SOX17", "FOXA2", "GSC", "CST1",
                                    "TBX6", "MESP1", "EOMES",
                                    "LHX1", "IRX3", "GATA6",
                                    "FOXF1", "HAND1", "HAND2",
                                    "POSTN", "LUM", "NID2",
                                    "MEF2C", "PECAM1", "CDH5", "TEK",
                                    "DLX5", "TFAP2A", "TFAP2C", "GATA3","GFP"),
            ncol = 6, cols = c("#ebebeb", "#7d3c98"))
dev.off()
# plot expression
sr.out <- "/home/yhw/bioinfo/project-wb/hs_gastruloids/analysis/graphs/gas_dev/GAS_all_marker_gene"
dir.create(sr.out, recursive = T)
marker.gene <- list(EPI = c("CDH1", "POU5F1", "SOX2", "NANOG", "LIN28AP1", "MYC", "KLF4", "TDGF1", "GDF3", "GFP"),
                    Organizer = c("FST", "OTX2", "GSC", "EOMES", "FOXH1", "CER1", "MIXL1", "DKK1", "NODAL",
                                  "FOXA2", "LEFTY1", "LEFTY2", "TBXT", "SOX17", "SP5", "NOG", "CHRD"),
                    Embryonic = c("TBXT", "MESP1", "PDGFRA", "LHX1", "IRX3", "GATA6", "MYL7", "MSX1",
                                  "MSX2", "FOXF1", "HAND1", "HAND2", "BMP4", "MAB21L2", "CST1", "MEF2C",
                                  "PECAM1", "CDH5", "TEK", "TFAP2A", "TFAP2C", "GATA3", "DLX5", "LHX5",
                                  "PAX6", "SOX3", "PAX3", "PAX7"),
                    Extraembryonic = c("LUM", "NID2", "ANXA1", "POSTN", "GATA2", "KRT7", "VIM", "GATA4",
                                       "DCN", "TP63", "KRT18", "NR2F2", "ISL1", "HEY1", "CDH10", "CTSV",
                                       "ARID5B", "PLAGL1", "CREB3L1", "HOXA10", "HOXA11", "HOXA9",
                                       "HOXA13", "WNT6", "GABRP", "HEY1", "BST2"))
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("GASD0_marker_gene_expression_of_", marker.gene[[type]][g], "_featureplot.pdf")),
        height = 3.5, width = 4)
    print(FeaturePlot(gas.merge, features = marker.gene[[type]][g],
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 0.5, slot = "data"))
    dev.off()
  }; rm(g)
}
# delete all useless data
rm(tmp, tmp1, tmp2, tmp3, tmp4)


### >>> 4. Gastruloid data before D5
# - reduce cells D5 before
keep.cell <- c(grep("ALCs", gas.d5before.merge$CellType),
               grep("ExMLCs", gas.d5before.merge$CellType),
               grep("PSLCs", gas.d5before.merge$CellType),
               sample(grep("EMLCs", gas.d5before.merge$CellType), 300),
               sample(grep("EndLCs", gas.d5before.merge$CellType), 500),
               sample(grep("EPILCs", gas.d5before.merge$CellType), 700))
table(gas.d5before.merge$CellType[keep.cell])
gas.d5before.merge.urd <- gas.d5before.merge[,keep.cell]
table(gas.d5before.merge.urd$CellType)
gas.d5before.merge.urd <- NormalizeData(gas.d5before.merge.urd) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(gas.d5before.merge.urd))
gas.d5before.merge.urd <- RunPCA(gas.d5before.merge.urd, features = VariableFeatures(gas.d5before.merge.urd))
ElbowPlot(gas.d5before.merge.urd)
dim.n <- 6
gas.d5before.merge.urd <- RunUMAP(gas.d5before.merge.urd, reduction = "pca", dims = 1:dim.n) %>%
  RunTSNE(dims = 1:dim.n)
set.seed(123)
gas.d5before.merge.urd <- RunHarmony(gas.d5before.merge.urd, group.by.vars = "Batch", plot_convergence = FALSE, project.dim = F)
gas.d5before.merge.urd <- FindNeighbors(gas.d5before.merge.urd, reduction = "harmony", dims = 1:dim.n) %>%
  FindClusters(resolution = 0.75)
gas.d5before.merge.urd <- RunUMAP(gas.d5before.merge.urd, reduction = "harmony", dims = 1:dim.n, return.model = T) %>%
  RunTSNE(reduction = "harmony", dims = 1:dim.n)
DimPlot(gas.d5before.merge.urd, reduction = "tsne", group.by = c("Sample_Name", "CellType"), pt.size = 1, label = T)
gas.d5before.merge.urd.tsne <- gas.d5before.merge.urd@reductions$tsne@cell.embeddings %>% as.data.frame()
# - reduce cells all
table(gas.merge$CellType)
keep.cell <- c(grep("^ALCs$", gas.merge$CellType),
               sample(grep("^AM$", gas.merge$CellType), 400),
               grep("^ECT$", gas.merge$CellType),
               grep("^EM$", gas.merge$CellType),
               sample(grep("^EMLCs$", gas.merge$CellType), 300),
               grep("^End$", gas.merge$CellType),
               sample(grep("^EndLCs$", gas.merge$CellType), 300),
               sample(grep("^EPI$", gas.merge$CellType), 300),
               sample(grep("^EPILCs$", gas.merge$CellType), 300),
               grep("^ExM$", gas.merge$CellType),
               grep("^ExMLCs$", gas.merge$CellType),
               grep("^HEP$", gas.merge$CellType),
               grep("^NM$", gas.merge$CellType),
               sample(grep("^PS$", gas.merge$CellType), 700),
               grep("^PSLCs$", gas.merge$CellType))
table(gas.merge$CellType[keep.cell])
gas.merge.urd <- gas.merge[,keep.cell]
table(gas.merge.urd$CellType)
gas.merge.urd <- NormalizeData(gas.merge.urd) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(gas.merge.urd))
gas.merge.urd <- RunPCA(gas.merge.urd, features = VariableFeatures(gas.merge.urd))
ElbowPlot(gas.merge.urd)
dim.n <- 10
gas.merge.urd <- RunUMAP(gas.merge.urd, reduction = "pca", dims = 1:dim.n) %>%
  RunTSNE(dims = 1:dim.n)
table(gas.merge.urd$Batch, gas.merge.urd$Sample_Name)
set.seed(123)
gas.merge.urd <- RunHarmony(gas.merge.urd, group.by.vars = "Batch", plot_convergence = FALSE, project.dim = F)
gas.merge.urd <- FindNeighbors(gas.merge.urd, reduction = "harmony", dims = 1:dim.n) %>%
  FindClusters(resolution = 0.75)
gas.merge.urd <- RunUMAP(gas.merge.urd, reduction = "harmony", dims = 1:dim.n, return.model = T) %>%
  RunTSNE(reduction = "harmony", dims = 1:dim.n)
p1 <- DimPlot(gas.merge.urd, reduction = "tsne", group.by = c("Sample_Name"), pt.size = 1, label = T, repel = T) +
  scale_colour_manual(values = pd.col.fix3)
p2 <- DimPlot(gas.merge.urd, reduction = "tsne", group.by = c("CellType"), pt.size = 1, label = T, repel = T) +
  scale_colour_manual(values = gas.merge.col)
p1 + p2
gas.merge.urd.tsne <- gas.merge.urd@reductions$tsne@cell.embeddings %>% as.data.frame()
# - URD pipeline before building trees
marker.genes <- c("TFAP2A","TFAP2C","DLX5","GATA3", # Ectoderm (ECT)
                  "CDH1","SOX2","POU5F1","OTX2", # Epiblast (EPI)
                  "TBXT","SP5","FST","MIXL1", # Primitive streak (PS)
                  "EOMES","MESP1","TBX6", # Nascent mesoderm (NM)
                  "PDGFRA","LHX1","IRX3","GATA6", # Emergent mesoderm (EM)
                  "FOXF1","HAND1","HAND2","MYL7","MSX1", # Advanced mesoderm (AM)
                  "LUM","NID2","POSTN","BMP4","ANXA1","MAB21L2", # Extraembryonic mesoderm (ExM)
                  "SOX17","GSC","FOXA2","CST1", # Endoderm (End)
                  "MEF2C","PECAM1","CDH5","TEK" # Haemato-endothelial progenitors (HEP)
                  )
# End_Meso
urd.gasd5.before <- list()
urd.gasd5.merge <- list()
urd.gasd5.before[["knn44_sigma0_End_Meso"]] <- URD.Btree(sr.obj = gas.d5before.merge.urd, group.by = c("CellType", "Type"),
                                                         tip.gp = c("EndLCs", "EMLCs"), start.cell = "EPILCs",
                                                         urd.knn = 44, urd.sigma = NULL, embedding = gas.d5before.merge.urd.tsne,
                                                         group.col = list(CellType = gas.d5before.merge.col,
                                                                          Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                         pt.gene = marker.genes,
                                                         res.out = "graphs/urd/GAS_D5_Before_knn44_sigma0_End_Meso")
urd.gasd5.merge[["knn58_sigma0_End_AM"]] <- URD.Btree(sr.obj = gas.merge.urd, group.by = c("CellType", "Type", "Sample_Name"),
                                                      tip.gp = c("End", "AM"), start.cell = "EPI",
                                                      urd.knn = 58, urd.sigma = NULL, embedding = gas.merge.urd.tsne,
                                                      group.col = list(CellType = gas.merge.col,
                                                                       Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279"),
                                                                       Sample_Name = c(ORG_D2 = "#4DAF4A", EPI_D0 = "#FF7F00",
                                                                                       GAS_D1 = "#9ECAE1", GAS_D3 = "#4292C6",
                                                                                       GAS_D5 = "#08519C")),
                                                      pt.gene = marker.genes,
                                                      res.out = "graphs/urd/GAS_merged_knn58_sigma0_End_AM")
# End_Meso_ExM
urd.gasd5.before[["knn44_sigma0_End_Meso_ExM"]] <- URD.Btree(sr.obj = gas.d5before.merge.urd, group.by = c("CellType", "Type"),
                                                             tip.gp = c("EndLCs", "EMLCs", "ExMLCs"), start.cell = "EPILCs",
                                                             urd.knn = 44, urd.sigma = NULL, embedding = gas.d5before.merge.urd.tsne,
                                                             group.col = list(CellType = gas.d5before.merge.col,
                                                                              Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                             pt.gene = marker.genes,
                                                             res.out = "graphs/urd/GAS_D5_Before_knn44_sigma0_End_Meso_ExM")
urd.gasd5.merge[["knn58_sigma0_End_AM_ExM"]] <- URD.Btree(sr.obj = gas.merge.urd, group.by = c("CellType", "Type", "Sample_Name"),
                                                          tip.gp = c("End", "AM", "ExM"), start.cell = "EPI",
                                                          urd.knn = 58, urd.sigma = NULL, embedding = gas.merge.urd.tsne,
                                                          group.col = list(CellType = gas.merge.col,
                                                                           Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279"),
                                                                           Sample_Name = c(ORG_D2 = "#4DAF4A", EPI_D0 = "#FF7F00",
                                                                                           GAS_D1 = "#9ECAE1", GAS_D3 = "#4292C6",
                                                                                           GAS_D5 = "#08519C")),
                                                          pt.gene = marker.genes,
                                                          res.out = "graphs/urd/GAS_merged_knn58_sigma0_End_AM_ExM")
# End_Meso_ExM_ECT
urd.gasd5.merge[["knn58_sigma0_End_AM_ExM_ECT"]] <- URD.Btree(sr.obj = gas.merge.urd, group.by = c("CellType", "Type", "Sample_Name"),
                                                              tip.gp = c("End", "AM", "ExM", "ECT"), start.cell = "EPI",
                                                              urd.knn = 58, urd.sigma = NULL, embedding = gas.merge.urd.tsne,
                                                              group.col = list(CellType = gas.merge.col,
                                                                               Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279"),
                                                                               Sample_Name = c(ORG_D2 = "#4DAF4A", EPI_D0 = "#FF7F00",
                                                                                               GAS_D1 = "#9ECAE1", GAS_D3 = "#4292C6",
                                                                                               GAS_D5 = "#08519C")),
                                                              pt.gene = marker.genes,
                                                              res.out = "graphs/urd/GAS_merged_knn58_sigma0_End_AM_ExM_ECT")
urd.gasd5.merge[["knn58_sigma0_End_AM_ExM_ECT_new"]] <- URD.Btree(sr.obj = gas.merge.urd, group.by = c("CellType", "Type", "Sample_Name"),
                                                                  tip.gp = c("End", "AM", "ExM", "ECT"), start.cell = "EPILCs",
                                                                  urd.knn = 58, urd.sigma = NULL, embedding = gas.merge.urd.tsne,
                                                                  group.col = list(CellType = gas.merge.col,
                                                                                   Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279"),
                                                                                   Sample_Name = c(ORG_D2 = "#4DAF4A", EPI_D0 = "#FF7F00",
                                                                                                   GAS_D1 = "#9ECAE1", GAS_D3 = "#4292C6",
                                                                                                   GAS_D5 = "#08519C")),
                                                                  pt.gene = marker.genes,
                                                                  res.out = "graphs/urd/GAS_merged_knn58_sigma0_End_AM_ExM_ECT_new")
# End_Meso_ECT
urd.gasd5.merge[["knn58_sigma0_End_AM_ECT"]] <- URD.Btree(sr.obj = gas.merge.urd, group.by = c("CellType", "Type", "Sample_Name"),
                                                          tip.gp = c("End", "AM", "ECT"), start.cell = "EPI",
                                                          urd.knn = 58, urd.sigma = NULL, embedding = gas.merge.urd.tsne,
                                                          group.col = list(CellType = gas.merge.col,
                                                                           Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279"),
                                                                           Sample_Name = c(ORG_D2 = "#4DAF4A", EPI_D0 = "#FF7F00",
                                                                                           GAS_D1 = "#9ECAE1", GAS_D3 = "#4292C6",
                                                                                           GAS_D5 = "#08519C")),
                                                          pt.gene = marker.genes,
                                                          res.out = "graphs/urd/GAS_merged_knn58_sigma0_End_AM_ECT")
urd.gasd5.merge[["knn58_sigma0_End_AM_ECT_new"]] <- URD.Btree(sr.obj = gas.merge.urd, group.by = c("CellType", "Type", "Sample_Name"),
                                                              tip.gp = c("End", "AM", "ECT"), start.cell = "EPILCs",
                                                              urd.knn = 58, urd.sigma = NULL, embedding = gas.merge.urd.tsne,
                                                              group.col = list(CellType = gas.merge.col,
                                                                               Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279"),
                                                                               Sample_Name = c(ORG_D2 = "#4DAF4A", EPI_D0 = "#FF7F00",
                                                                                               GAS_D1 = "#9ECAE1", GAS_D3 = "#4292C6",
                                                                                               GAS_D5 = "#08519C")),
                                                              pt.gene = marker.genes,
                                                              res.out = "graphs/urd/GAS_merged_knn58_sigma0_End_AM_ECT_new")
# - building trees
tree.gasd5.before <- list()
tree.gasd5.merge <- list()
# End_Meso
tree.gasd5.before[["knn44_sigma0_End_Meso"]] <- URD.Atree(urd = urd.gasd5.before[["knn44_sigma0_End_Meso"]],
                                                          tip.gp = c("EndLCs", "EMLCs"),
                                                          tip.cluster = c("3", "2"),
                                                          group.by = c("CellType", "Type"),
                                                          group.col = list(CellType = gas.d5before.merge.col,
                                                                           Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                          pt.gene = marker.genes,
                                                          res.out = "graphs/urd/GAS_D5_Before_knn44_sigma0_End_Meso")
tree.gasd5.merge[["knn58_sigma0_End_AM"]] <- URD.Atree(urd = urd.gasd5.merge[["knn58_sigma0_End_AM"]],
                                                       tip.gp = c("End", "AM"),
                                                       tip.cluster = c("1", "2"),
                                                       group.by = c("CellType", "Type", "Sample_Name"),
                                                       group.col = list(CellType = gas.merge.col,
                                                                        Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279"),
                                                                        Sample_Name = c(ORG_D2 = "#4DAF4A", EPI_D0 = "#FF7F00",
                                                                                        GAS_D1 = "#9ECAE1", GAS_D3 = "#4292C6",
                                                                                        GAS_D5 = "#08519C")),
                                                       pt.gene = marker.genes,
                                                       res.out = "graphs/urd/GAS_merged_knn58_sigma0_End_AM")
# End_Meso_ExM
tree.gasd5.before[["knn44_sigma0_End_Meso_ExM"]] <- URD.Atree(urd = urd.gasd5.before[["knn44_sigma0_End_Meso_ExM"]],
                                                              tip.gp = c("EndLCs", "EMLCs", "ExMLCs"),
                                                              tip.cluster = c("1", "2", "7"),
                                                              group.by = c("CellType", "Type"),
                                                              group.col = list(CellType = gas.d5before.merge.col,
                                                                               Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279")),
                                                              pt.gene = marker.genes,
                                                              res.out = "graphs/urd/GAS_D5_Before_knn44_sigma0_End_Meso_ExM")
tree.gasd5.merge[["knn58_sigma0_End_AM_ExM"]] <- URD.Atree(urd = urd.gasd5.merge[["knn58_sigma0_End_AM_ExM"]],
                                                           tip.gp = c("End", "AM", "ExM"),
                                                           tip.cluster = c("3", "2", "5"),
                                                           group.by = c("CellType", "Type", "Sample_Name"),
                                                           group.col = list(CellType = gas.merge.col,
                                                                            Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279"),
                                                                            Sample_Name = c(ORG_D2 = "#4DAF4A", EPI_D0 = "#FF7F00",
                                                                                            GAS_D1 = "#9ECAE1", GAS_D3 = "#4292C6",
                                                                                            GAS_D5 = "#08519C")),
                                                           pt.gene = marker.genes,
                                                           res.out = "graphs/urd/GAS_merged_knn58_sigma0_End_AM_ExM")
# End_Meso_ExM_ECT
tree.gasd5.merge[["knn58_sigma0_End_AM_ExM_ECT"]] <- URD.Atree(urd = urd.gasd5.merge[["knn58_sigma0_End_AM_ExM_ECT"]],
                                                               tip.gp = c("End", "AM", "ExM", "ECT"),
                                                               tip.cluster = c("7", "1", "5", "2"),
                                                               group.by = c("CellType", "Type", "Sample_Name"),
                                                               group.col = list(CellType = gas.merge.col,
                                                                                Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279"),
                                                                                Sample_Name = c(ORG_D2 = "#4DAF4A", EPI_D0 = "#FF7F00",
                                                                                                GAS_D1 = "#9ECAE1", GAS_D3 = "#4292C6",
                                                                                                GAS_D5 = "#08519C")),
                                                               pt.gene = marker.genes,
                                                               res.out = "graphs/urd/GAS_merged_knn58_sigma0_End_AM_ExM_ECT")
tree.gasd5.merge[["knn58_sigma0_End_AM_ExM_ECT_new"]] <- URD.Atree(urd = urd.gasd5.merge[["knn58_sigma0_End_AM_ExM_ECT_new"]],
                                                                   tip.gp = c("End", "AM", "ExM", "ECT"),
                                                                   tip.cluster = c("7", "1", "5", "2"),
                                                                   group.by = c("CellType", "Type", "Sample_Name"),
                                                                   group.col = list(CellType = gas.merge.col,
                                                                                    Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279"),
                                                                                    Sample_Name = c(ORG_D2 = "#4DAF4A", EPI_D0 = "#FF7F00",
                                                                                                    GAS_D1 = "#9ECAE1", GAS_D3 = "#4292C6",
                                                                                                    GAS_D5 = "#08519C")),
                                                                   pt.gene = marker.genes,
                                                                   res.out = "graphs/urd/GAS_merged_knn58_sigma0_End_AM_ExM_ECT_new")
# End_Meso_ECT
tree.gasd5.merge[["knn58_sigma0_End_AM_ECT"]] <- URD.Atree(urd = urd.gasd5.merge[["knn58_sigma0_End_AM_ECT"]],
                                                           tip.gp = c("End", "AM", "ECT"),
                                                           tip.cluster = c("8", "1", "7"),
                                                           group.by = c("CellType", "Type", "Sample_Name"),
                                                           group.col = list(CellType = gas.merge.col,
                                                                            Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279"),
                                                                            Sample_Name = c(ORG_D2 = "#4DAF4A", EPI_D0 = "#FF7F00",
                                                                                            GAS_D1 = "#9ECAE1", GAS_D3 = "#4292C6",
                                                                                            GAS_D5 = "#08519C")),
                                                           pt.gene = marker.genes,
                                                           res.out = "graphs/urd/GAS_merged_knn58_sigma0_End_AM_ECT")
tree.gasd5.merge[["knn58_sigma0_End_AM_ECT_new"]] <- URD.Atree(urd = urd.gasd5.merge[["knn58_sigma0_End_AM_ECT_new"]],
                                                               tip.gp = c("End", "AM", "ECT"),
                                                               tip.cluster = c("8", "1", "7"),
                                                               group.by = c("CellType", "Type", "Sample_Name"),
                                                               group.col = list(CellType = gas.merge.col,
                                                                                Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279"),
                                                                                Sample_Name = c(ORG_D2 = "#4DAF4A", EPI_D0 = "#FF7F00",
                                                                                                GAS_D1 = "#9ECAE1", GAS_D3 = "#4292C6",
                                                                                                GAS_D5 = "#08519C")),
                                                               pt.gene = marker.genes,
                                                               res.out = "graphs/urd/GAS_merged_knn58_sigma0_End_AM_ECT_new")
# - replot
pd.gene <- c("HSPG2", "DAG1")
for (name in names(urd.gasd5.before)) {
  URD.Replot(urd = urd.gasd5.before[[name]], tree = tree.gasd5.before[[name]], mode = c("gene", "group"),
             group.by = c("CellType", "Type", "Sample_Name"),
             group.col = list(CellType = gas.d5before.merge.col,
                              Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279"),
                              Sample_Name = c(ORG_D2 = "#4CAF49", EPI_D0 = "#B91122",
                                              GAS_D1 = "#9ECAE1", GAS_D3 = "#EB8838", GAS_D5 = "#08519C")),
             expr.col = c("#ebebeb", "#7d3c98"), pt.gene = pd.gene, pt.size = 1.5,
             res.out = paste0("graphs/urd/Replot_", name))
}
for (name in names(urd.gasd5.merge)) {
  URD.Replot(urd = urd.gasd5.merge[[name]], tree = tree.gasd5.merge[[name]], mode = c("gene", "group"),
             group.by = c("CellType", "Type", "Sample_Name"),
             group.col = list(CellType = gas.merge.col,
                              Type = c(GFP_neg = "#AFC5E4", GFP_pos = "#BBC279"),
                              Sample_Name = c(ORG_D2 = "#4CAF49", EPI_D0 = "#B91122",
                                              GAS_D1 = "#9ECAE1", GAS_D3 = "#EB8838", GAS_D5 = "#08519C")),
             expr.col = c("#ebebeb", "#7d3c98"), pt.gene = pd.gene, pt.size = 1.5,
             res.out = paste0("graphs/urd/Replot_", name))
}



### ===========================
### 8th part: CellChat analysis ----
### ===========================


### >>> 1. Setting output directory
outdir <- file.path(getwd(), "graphs/cellchat")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}


### >>> 2. CellChat analysis
library("AfterChat")
library("CellChat")
db.cellchat.hs <- CellChatDB.human
write.csv(db.cellchat.hs$interaction, file.path(outdir, "pathway_database.csv"),
          col.names = T, row.names = F)
# gastruloid - between cell types
table(gas.merge.urd$Type, gas.merge.urd$CellType)
gas.merge.urd@meta.data <- gas.merge.urd@meta.data %>%
  mutate(Sample_Group = case_when(Sample_Name %in% c("EPI_D0", "ORG_D2") ~ "GAS_D0",
                                  TRUE ~ Sample_Name))
table(gas.merge.urd$Sample_Group, gas.merge.urd$CellType)
gas.ct <- list()
gas.ct <- Pipe.CellChat(sc.obj = gas.merge.urd, cell.db = db.cellchat.hs,
                        group.by = "CellType", split.by = "Sample_Group")
for (i in names(gas.ct)) {
  gas.ct[[i]] <- computeCommunProb(gas.ct[[i]], type = "truncatedMean", trim = 0.1) %>%
    computeCommunProbPathway() %>%
    aggregateNet()
}
for (i in names(gas.ct)[4]) {
  #AfterChat::PathCentrality(ct.obj = gas.ct[[i]], outdir = file.path(outdir, "Gastruloids"),
  #                          file.prefix = i)
  #AfterChat::PathInteracion(ct.obj = gas.ct[[i]], outdir = file.path(outdir, "Gastruloids"),
  #                          file.prefix = i)
  #AfterChat::LRsContribution(ct.obj = gas.ct[[i]], outdir = file.path(outdir, "Gastruloids"),
  #                           file.prefix = i)
  LRsInteraction(ct.obj = gas.ct[[i]], outdir = file.path(outdir, "Gastruloids"), file.prefix = i,
                 cell.source = c("ALCs", "AM", "ECT", "EM", "EMLCs", "End", "EndLCs", "ExM", "ExMLCs", "HEP", "NM", "PS", "PSLCs"),
                 cell.target = c("EPI", "EPILCs"))
  #AfterChat::PathClustering(ct.obj = gas.ct[[i]], outdir = file.path(outdir, "Gastruloids"),
  #                          file.prefix = i)
}
# gastruloid - GFP + EPI
gas.merge.urd@meta.data <- gas.merge.urd@meta.data %>%
  mutate(Sample_Group = case_when(Sample_Name %in% c("EPI_D0", "ORG_D2") ~ "GAS_D0",
                                  TRUE ~ Sample_Name),
         CellType.ct = case_when(CellType %in% c("EPILCs", "EPI") & Type == "GFP_neg" ~ "EPI",
                                 (Type == "GFP_pos") & !(CellType %in% c("EPILCs", "EPI")) ~ "GFP"))
gas.merge.urd.ct <- subset(gas.merge.urd, CellType.ct %in% c("EPI", "GFP"))
table(gas.merge.urd.ct$CellType.ct, gas.merge.urd.ct$Sample_Group)
gas.ct.two <- list()
gas.ct.two <- Pipe.CellChat(sc.obj = gas.merge.urd.ct, cell.db = db.cellchat.hs,
                            group.by = "CellType.ct", split.by = "Sample_Group")
for (i in names(gas.ct.two)) {
  gas.ct.two[[i]] <- computeCommunProb(gas.ct.two[[i]], type = "truncatedMean", trim = 0.1) %>%
    computeCommunProbPathway() %>%
    aggregateNet()
}
for (i in names(gas.ct.two)) {
  #AfterChat::PathCentrality(ct.obj = gas.ct.two[[i]], outdir = file.path(outdir, "Gastruloids_two"),
  #                          file.prefix = i)
  #AfterChat::PathInteracion(ct.obj = gas.ct.two[[i]], outdir = file.path(outdir, "Gastruloids_two"),
  #                          file.prefix = i)
  #AfterChat::LRsContribution(ct.obj = gas.ct.two[[i]], outdir = file.path(outdir, "Gastruloids_two"),
  #                           file.prefix = i)
  LRsInteraction(ct.obj = gas.ct.two[[i]], outdir = file.path(outdir, "Gastruloids_two"),
                 file.prefix = i,
                 cell.source = c("GFP"),
                 cell.target = c("EPI"))
  #AfterChat::PathClustering(ct.obj = gas.ct.two[[i]], outdir = file.path(outdir, "Gastruloids_two"),
  #                          file.prefix = i)
}



### ===========================
### 9th part: A-P axis analysis ----
### ===========================


### >>> 1. Setting output directory
outdir <- file.path(getwd(), "graphs/AP_axis")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}


### >>> 2. Cluter genes
# load markers
marker.gene <- list(anterior = c("POU5F1","NANOG","TDGF1","NES","ID1","CDH1","SOX2","POU3F1","SOX3","IRX2",
                                 "CRABP1","PBX1","BTG1","TFAP2A","SFRP1","SFRP2","OTX2","FGFR1","HESX1","LHX5"),
                    medial = c("TBXT","TBX6","MESP1","EOMES","IRX3","LHX1","FOXC1","HOXB1","HOXA1","MEIS2","HOXB2",
                               "HOXB3","HOXA3","DLL1","GATA6","GATA4","FST","MIXL1","SP5","CER1","GFP"),
                    posterior = c("CDX2","CDX4","HOXB9","HOXC4","HOXC8","HOXC6","HOXB7","HOXB8","HOXB6","HOXB5",
                                  "VIM","SHH","FOXA2","CDH2","DUSP6","HOXA10","HOXA13","HAND1","HAND2","MSX1","GFP"))
# clustering
mfuzz.expr <- AverageExpression(sr.pdt$gasd5.to.hs, group.by = "CellType",
                                assays = "RNA", slot = "data")
mfuzz.expr <- mfuzz.expr$RNA[,c("EPI", "ECT", "PS", "NM", "EM", "AM", "ExM", "End", "HEP")] %>%
  as.data.frame()
colnames(mfuzz.expr) <- paste0("S", 1:ncol(mfuzz.expr))
mfuzz <- list()
for (i in 9:12) {
  name <- paste0("GASD5_C", i)
  # run Mfuzz
  mfuzz[[name]] <- Pipe.Mfuzz(expr.mtx = pd, res.out = file.path(outdir, "GASD5"),
                              prefix = name, cluster.num = i, min.acore = 0.5)
  # plot expression dynamics of gene in clusters over timeline
  Plot.Mfuzz(mfuzz = mfuzz[[name]], geneset = NULL,
             mode = "line",
             prefix = name,
             res.out = file.path(outdir, "GASD5"))
}
# line plot (grouped genes)
for (i in 9:12) {
  name <- paste0("GASD5_C", i)
  # plot custom gene over timeline
  for (type in names(marker.gene)) {
    Plot.Mfuzz(expr = mfuzz.expr,
               mfuzz = mfuzz[[name]],
               geneset = marker.gene[[type]],
               mode = "line",
               prefix = paste0(name, "-", type),
               res.out = file.path(outdir, "GASD5"))
  }
}
# line plot (single genes)
for (i in 9:12) {
  name <- paste0("GASD5_C", i)
  # plot custom gene over timeline
  for (type in names(marker.gene)) {
    for (j in marker.gene[[type]]) {
      Plot.Mfuzz(expr = mfuzz.expr,
                 mfuzz = mfuzz[[name]],
                 geneset = j,
                 mode = "line",
                 prefix = paste0(name, "-", type, "-", j),
                 res.out = file.path(outdir, paste0("GASD5/", name)))
    }
  }
}
# heatmap plot
mfuzz.expr <- GetAssayData(sr.pdt$gasd5.to.hs, slot = "scale.data")
cell.seq <- order(match(sr.pdt$gasd5.to.hs$CellType,
                        c("EPI", "ECT", "PS", "NM", "EM", "AM", "ExM", "End", "HEP")))
cell.label <- as.character(sr.pdt$gasd5.to.hs@meta.data$CellType[cell.seq])
mfuzz.expr <- mfuzz.expr[, cell.seq]
for (i in 9:12) {
  name <- paste0("GASD5_C", i)
  # plot custom gene over timeline
  for (type in names(marker.gene)) {
    Plot.Mfuzz(expr = mfuzz.expr,
               mfuzz = mfuzz[[name]],
               geneset = marker.gene[[type]],
               mode = "heatmap",
               cell.group = cell.label,
               cell.col = pd.col.fix,
               prefix = paste0(name, "-", type),
               res.out = file.path(outdir, "GASD5"))
  }
}



### ====================
### 10th part: Organizer ----
### ====================


### >>> 1. Setting output directory
outdir <- file.path(getwd(), "graphs/organizer")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}
mk.col <- c("#607d8b","#795548","#ff5722","#ffc107","#cddc39","#4caf50","#009688",
            "#00bcd4","#2196f3","#3f51b5","#673ab7","#9c27b0","#e91e63","#f44336",
            "#b0bec5","#bcaaa4","#ffab91","#ffe082","#e6ee9c","#a5d6a7","#80cbc4",
            "#80deea","#90caf9","#9fa8da","#b39ddb","#ce93d8","#f48fb1","#ef9a9a",
            "#37474f","#4e342e","#d84315","#ff8f00","#9e9d24","#2e7d32","#00695c",
            "#00838f","#1565c0","#283593","#4527a0","#6a1b9a","#ad1457","#c62828")
organizer.gene <- c("GSC", "OTX2", "FOXA2", "FST", "CER1", "DKK1", "HHEX", "CHRD", "T")


### >>> 2. Process monkey data
mf.embryo <- Read10X("/home/yhw/document/public_data/GSE193007_monkey_post-gastrulation/GSE193007_MFE-filtered_feature_bc_matrix")
meta <- read.csv("/home/yhw/document/public_data/GSE193007_monkey_post-gastrulation/MFE56636-meta.csv", row.names = 1)
mf.embryo <- mf.embryo[,str_split_fixed(rownames(meta), "_", 2)[,2]]
rownames(meta) <- str_split_fixed(rownames(meta), "_", 2)[,2]
if (all(colnames(mf.embryo) == rownames(meta))) {
  mf.embryo <- CreateSeuratObject(counts = mf.embryo, meta.data = meta)
}; rm(meta)
mf.embryo <- subset(mf.embryo, sample %in% c("CS8-e1", "CS8-e2", "CS9-e1", "CS9-e2"))
mf.embryo <- NormalizeData(mf.embryo) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
mf.embryo <- ScaleData(mf.embryo, verbose = FALSE, features = rownames(mf.embryo))
mf.embryo <- RunPCA(mf.embryo, npcs = 50, verbose = FALSE, features = VariableFeatures(mf.embryo))
ElbowPlot(mf.embryo, ndims = 30)
dim.n <- 15
mf.embryo <- FindNeighbors(mf.embryo, reduction = "pca", dims = 1:dim.n) %>%
  FindClusters(resolution = 2)
mf.embryo <- RunUMAP(mf.embryo, reduction = "pca", dims = 1:dim.n, return.model = T) %>%
  RunTSNE(reduction = "pca", dims = 1:dim.n)
# plot clustering
pdf(file.path(outdir, "Scatter_plot_to_show_monkey_data_clustering.pdf"), width = 13.5, height = 5)
DimPlot(mf.embryo, reduction = "umap", group.by = c("stage", "cell_type"),
        label = TRUE, repel = TRUE, pt.size = 0.35, cols = mk.col)
dev.off()
# plot expression
sr.out <- file.path(outdir, "Marker_gene_monkey_cells")
dir.create(sr.out, recursive = T)
marker.gene <- list(Organizer = c("FST", "OTX2", "GSC", "EOMES", "FOXH1", "CER1", "MIXL1", "DKK1", "NODAL",
                                  "FOXA2", "LEFTY1", "LEFTY2", "T", "SP5", "NOG", "CHRD"),
                    Buchong = c("NANOG", "HHEX", "CDH1", "POU3F1", "CDX2", "CDX1", "SOX17",
                                "GATA6", "GATA4", "LUM", "NID2", "POU5F1", "SOX2", "MYC",
                                "LIN28", "GDF3", "TDGF1", "SNAIL", "SNAIL2", "TBX6",
                                "MESP1", "LHX1", "IRX3", "HAND1", "FOXF1"),
                    FMH.marker  = c("LHX5","ZNF616","GBX2","RREB1","TEAD4","EN2","HAND1",
                                    "HESX1","FOXG1","NFATC4","PAX7","DBX1","PAX5","SMARCA5",
                                    "PAX6","IRX3","EBF1","HOXA3","MEIS1","WRNIP1","MAFB"),
                    FMH = c("PAX6","OTX2","EN1","PAX8","GBX2","EGR2","PAX3","DBX2","DBX1","IRX3","OLIG2","NKX2.2"),
                    Neural = c("SOX1","SOX3","PAX6","TUBB3","ELAVL3","OLIG2","NEUROG1","NEUROD1","NIKX2.2"),
                    NMP = c("LIMD1","MARK3","CIT","MAPK14","FAT4","LATS2","STK4","WWC1","SHANK2","SOX11","TEAD1","PJA2",
                            "TEAD2","VGLL5","SAV1","MAP2K3","TIAL1","NF2","STK3","LATS1","THBS1","TEAD3","AMOT","YAP1",
                            "AMOTL1","AMOTL2","NEK8","WTIP","DLG5","TEAD4"),
                    Amnion = c("GABRP","TGFB1","VTCN1","S100A10","GADD45G","ACTC1",
                               "IGFBP3","MCM5","CLDN10","TUBG1","UCP2","BIN3","GGA2","FGFR1","GALM"),
                    EPI = c("CDH1", "POU5F1", "SOX2", "NANOG", "LIN28AP1", "MYC", "KLF4", "TDGF1", "GDF3", "GFP"),
                    Embryonic = c("TBXT", "MESP1", "PDGFRA", "LHX1", "IRX3", "GATA6", "MYL7", "MSX1",
                                  "MSX2", "FOXF1", "HAND1", "HAND2", "BMP4", "MAB21L2", "CST1", "MEF2C",
                                  "PECAM1", "CDH5", "TEK", "TFAP2A", "TFAP2C", "GATA3", "DLX5", "LHX5",
                                  "PAX6", "SOX3", "PAX3", "PAX7"),
                    Extraembryonic = c("LUM", "NID2", "ANXA1", "POSTN", "GATA2", "KRT7", "VIM", "GATA4",
                                       "DCN", "TP63", "KRT18", "NR2F2", "ISL1", "HEY1", "CDH10", "CTSV",
                                       "ARID5B", "PLAGL1", "CREB3L1", "HOXA10", "HOXA11", "HOXA9",
                                       "HOXA13", "WNT6", "GABRP", "HEY1", "BST2"),
                    anterior = c("POU5F1","NANOG","TDGF1","NES","ID1","CDH1","SOX2","POU3F1","SOX3","IRX2",
                                 "CRABP1","PBX1","BTG1","TFAP2A","SFRP1","SFRP2","OTX2","FGFR1","HESX1","LHX5"),
                    medial = c("TBXT","TBX6","MESP1","EOMES","IRX3","LHX1","FOXC1","HOXB1","HOXA1","MEIS2","HOXB2",
                               "HOXB3","HOXA3","DLL1","GATA6","GATA4","FST","MIXL1","SP5","CER1","GFP"),
                    posterior = c("CDX2","CDX4","HOXB9","HOXC4","HOXC8","HOXC6","HOXB7","HOXB8","HOXB6","HOXB5",
                                  "VIM","SHH","FOXA2","CDH2","DUSP6","HOXA10","HOXA13","HAND1","HAND2","MSX1","GFP"),
                    megakaryocyte = c("GP1BB", "ITGA2B", "NFE2"),
                    erythroid = c("GATA1", "KLF1", "GYPB", "HBE1"),
                    myeloid.progenitor = c("CD36", "CSF1R", "LYVE1", "KIT", "CSF1R", "MYB", "SPI1", "CD34", "CD45", "CD52", "NFE2"),
                    wanghongmei = c("AFP","ANXA1","ANXA8","APOE","APOE ","CDH1","CDH2","CDX1","CDX2","CER1","CHRD",
                                    "COL1A1","COL3A1","DNMT3B","DPPA3","EOMES","FGF8","FGF8 HES7","FOLR2","FOXA2",
                                    "FOXC1","FOXC2","FOXF1","GATA1","GATA4","GATA6","GP1BB","HAND1","HBE1","HBM",
                                    "HBZ","HES7","HHEX","IGFBP5","IRX3","ISL1","KDR","KRT18","KRT8","LEFTY2","LFNG",
                                    "MEF2C","MIXL1","MPO","MYL7","NANOG","NANOS3","NKX2-5","NOTO","OSR1","OTX2","PAX3",
                                    "PAX3 ","PAX6","PAX7","PAX8","PECAM1","POU3F1","POU5F1","POU5FA","PPBP","PRDM14",
                                    "RUNX1","SIX1","SOX10","SOX17","SOX2","SOX2 ","SOX9","TBXT","TBX6","TCF15","TCF21",
                                    "TFAP2A","TFAP2C","TFC15","MESP1","TTR","UCHL1","VWF","WNT3A","WNT5B","WT1"))
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(mf.embryo))
}
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(mf.embryo, features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 0.35, slot = "data"))
    dev.off()
  }; rm(g)
}


### >>> 3. Define organizers in monkey embryo
# filter cells according to marker gene expression
for (gene in organizer.gene) {
  mf.embryo@meta.data[, gene] <- GetAssayData(mf.embryo, slot = "data")[gene, ]
}
keep.cell <- apply(mf.embryo@meta.data[, organizer.gene[1:3]], 1, function(x){sum(x > 0)}) >= 3
table(mf.embryo$cell_type[keep.cell])
mf.embryo.organizer <- mf.embryo[, keep.cell]
mf.embryo.organizer <- subset(mf.embryo.organizer, cell_type %in% c("APS","DE","ECT","Nas.Meso","Node","VE"))
# seurat pipeline
mf.embryo.organizer <- NormalizeData(mf.embryo.organizer) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
mf.embryo.organizer <- ScaleData(mf.embryo.organizer, verbose = FALSE,
                                 features = rownames(mf.embryo.organizer))
mf.embryo.organizer <- RunPCA(mf.embryo.organizer, npcs = 50, verbose = FALSE,
                              features = VariableFeatures(mf.embryo.organizer))
ElbowPlot(mf.embryo.organizer, ndims = 30)
dim.n <- 6
mf.embryo.organizer <- FindNeighbors(mf.embryo.organizer, reduction = "pca", dims = 1:dim.n) %>%
  FindClusters(resolution = 0.75)
mf.embryo.organizer <- RunUMAP(mf.embryo.organizer, reduction = "pca", dims = 1:dim.n, return.model = T) %>%
  RunTSNE(reduction = "pca", dims = 1:dim.n)
# plot clustering
DimPlot(mf.embryo.organizer, reduction = "umap", group.by = c("stage", "cell_type"),
        label = TRUE, repel = TRUE, pt.size = 1, cols = mk.col)
# remove batch effect
mf.embryo.organizer <- RunHarmony(mf.embryo.organizer, group.by.vars = "sample")
mf.embryo.organizer <- RunUMAP(mf.embryo.organizer, reduction = "harmony", dims = 1:dim.n, return.model = T)
pdf(file.path(outdir, "Scatter_plot_to_show_monkey_organizer_clustering.pdf"), width = 11.25, height = 4.5)
DimPlot(mf.embryo.organizer, reduction = "umap", group.by = c("stage", "cell_type"),
        label = TRUE, repel = TRUE, pt.size = 2.5, cols = mk.col)
dev.off()
DimPlot(mf.embryo.organizer, reduction = "umap", group.by = c("stage", "cell_type", "seurat_clusters"),
        label = TRUE, repel = TRUE, pt.size = 1, cols = mk.col)
# calculate gene set score
library("msigdbr")
hs.msigdbr <- msigdbr(species = "Homo sapiens")
hs.msigdbr.final <- hs.msigdbr %>%
  filter(gs_subcat %in% c("GO:BP", "CP:KEGG", "CP:REACTOME"))
genesets <- unique(grep("PATHWAY", hs.msigdbr.final$gs_name, value = T))
hs.msigdbr.final <- subset(hs.msigdbr.final, gs_name %in% genesets)
table(hs.msigdbr.final$gs_name) %>% length()
hs.msigdbr.final <- with(hs.msigdbr.final, split(gene_symbol, gs_name))
names(hs.msigdbr.final) <- str_to_title(names(hs.msigdbr.final))
mf.embryo.organizer.cal <- ScoreGeneset(sr.obj = mf.embryo.organizer, gene.set = hs.msigdbr.final)
for (i in names(mf.embryo.organizer.cal@assays)[2]) {
  DefaultAssay(mf.embryo.organizer.cal) <- i
  mf.embryo.organizer.cal[[i]]@scale.data <- mf.embryo.organizer.cal[[i]]@data
}
for (i in names(mf.embryo.organizer.cal@assays)[-1:-2]) {
  DefaultAssay(mf.embryo.organizer.cal) <- i
  mf.embryo.organizer.cal <- ScaleData(mf.embryo.organizer.cal, 
                                       features = rownames(mf.embryo.organizer.cal))
}
for (i in row.names(GetAssayData(mf.embryo.organizer.cal, assay = "AUCell"))) {
  pdf(file.path(outdir, paste0("Pathway_", i, ".pdf")), height = 5, width = 10)
  print(VlnPlot(mf.embryo.organizer.cal, features = gsub("_", "-", i),
                assay = "AUCell", cols = mk.col, pt.size = 0) +
          NoLegend() + ggtitle("AUCell"))
  dev.off()
}
# plot expression
sr.out <- file.path(outdir, "Marker_gene_of_monkey_organizer")
dir.create(sr.out, recursive = T)
marker.gene <- list(Organizer = c("FST", "OTX2", "GSC", "EOMES", "FOXH1", "CER1", "MIXL1", "DKK1", "NODAL",
                                  "FOXA2", "LEFTY1", "LEFTY2", "T", "SP5", "NOG", "CHRD"),
                    Ligands = intersect(unique(CellChatDB.human$interaction$ligand), rownames(mf.embryo.organizer)),
                    Buchong = intersect(c("NANOG", "HHEX", "CDH1", "POU3F1", "CDX2", "CDX1", "SOX17",
                                          "GATA6", "GATA4", "LUM", "NID2", "POU5F1", "SOX2", "MYC",
                                          "LIN28", "GDF3", "TDGF1", "SNAIL", "SNAIL2", "TBX6",
                                          "MESP1", "LHX1", "IRX3", "HAND1", "FOXF1"), rownames(mf.embryo.organizer)))
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(mf.embryo.organizer, features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 2, slot = "data"))
    dev.off()
  }; rm(g)
}
# plot multiple markers
pdf(file.path(outdir, "Scatter_plot_to_show_monkey_organizer_gene_expression_GSC_OTX2_FOXA2_FST_CER1_DKK1_HHEX_CHRD.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = mf.embryo.organizer,
                        features = c("GSC", "OTX2", "FOXA2", "FST", "CER1", "DKK1", "HHEX", "CHRD"), pt.size = 2)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_monkey_organizer_gene_expression_GSC_OTX2_FOXA2_FST_CER1_DKK1_HHEX.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = mf.embryo.organizer,
                        features = c("GSC", "OTX2", "FOXA2", "FST", "CER1", "DKK1", "HHEX"), pt.size = 2)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_monkey_organizer_gene_expression_GSC_OTX2_FOXA2_FST_CER1_DKK1.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = mf.embryo.organizer, features = c("GSC", "OTX2", "FOXA2", "FST", "CER1", "DKK1"), pt.size = 2)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_monkey_organizer_gene_expression_GSC_OTX2_FOXA2_FST_CER1.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = mf.embryo.organizer, features = c("GSC", "OTX2", "FOXA2", "FST", "CER1"), pt.size = 2)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_monkey_organizer_gene_expression_GSC_OTX2_FOXA2_FST.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = mf.embryo.organizer, features = c("GSC", "OTX2", "FOXA2", "FST"), pt.size = 2)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_monkey_organizer_gene_expression_GSC_OTX2_FOXA2.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = mf.embryo.organizer, features = c("GSC", "OTX2", "FOXA2"), pt.size = 2)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_monkey_organizer_gene_expression_OTX2_FOXA2.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = mf.embryo.organizer, features = c("OTX2", "FOXA2"), pt.size = 2)
dev.off()
# reclustering with LRs
library("CellChat")
lrs <- unique(c(unique(CellChatDB.human$interaction$ligand),
                unique(CellChatDB.human$interaction$receptor),
                str_split_fixed(grep("_", unique(CellChatDB.human$interaction$receptor), value = T), "_", 2)[,1],
                str_split_fixed(grep("_", unique(CellChatDB.human$interaction$receptor), value = T), "_", 2)[,2]))
mf.embryo.organizer.lrs <- RunPCA(mf.embryo.organizer, npcs = 50, verbose = FALSE,
                                  features = intersect(VariableFeatures(mf.embryo.organizer), lrs))
ElbowPlot(mf.embryo.organizer.lrs, ndims = 30)
dim.n <- 8
mf.embryo.organizer.lrs <- FindNeighbors(mf.embryo.organizer.lrs, reduction = "pca", dims = 1:dim.n) %>%
  FindClusters(resolution = 2)
mf.embryo.organizer.lrs <- RunUMAP(mf.embryo.organizer.lrs, reduction = "pca", dims = 1:dim.n, return.model = T) %>%
  RunTSNE(reduction = "pca", dims = 1:dim.n)
DimPlot(mf.embryo.organizer.lrs, reduction = "umap", group.by = c("stage", "cell_type"),
        label = TRUE, repel = TRUE, pt.size = 1, cols = mk.col)
mf.embryo.organizer.lrs <- RunHarmony(mf.embryo.organizer.lrs, group.by.vars = "sample")
mf.embryo.organizer.lrs <- RunUMAP(mf.embryo.organizer.lrs, reduction = "harmony", dims = 1:dim.n, return.model = T)
#pdf(file.path(outdir, "Scatter_plot_to_show_monkey_organizer_clustering.pdf"), width = 11.5, height = 4.5)
DimPlot(mf.embryo.organizer.lrs, reduction = "umap", group.by = c("stage", "cell_type"),
        label = TRUE, repel = TRUE, pt.size = 1, cols = mk.col)
#dev.off()
Plot_Density_Joint_Only(seurat_object = mf.embryo.organizer.lrs, features = c("OTX2", "FOXA2"), pt.size = 2)


### >>> 3. Define organizers in monkey embryo
# extract cells
cell.coor <- data.frame(group = c("Organizer.1", "Organizer.2", "Organizer.3", "Organizer.4"),
                        xmin = c(-4, 0, 6, 3),
                        xmax = c(-2.5, 3, 7.25, 7),
                        ymin = c(-4.5, 2, 1.5, -2.5),
                        ymax = c(-3.5, 2.75, 2.75, -0.5))
organizer.cells <- ExtractCellByPos(object = mf.embryo.organizer, object.type = "seurat",
                                    dim.name = "umap", group.coor = cell.coor, group.name = "CellType.sub",
                                    pt.size = 2)
pdf(file.path(outdir, "Scatter_plot_to_show_monkey_organizer_4_populations_clustering.pdf"), width = 6, height = 4.5)
organizer.cells$plot
dev.off()
# add new cell types
mf.embryo.organizer@meta.data <- mf.embryo.organizer@meta.data %>%
  rownames_to_column(var = "barcode") %>%
  mutate(CellType.organizer = case_when(barcode %in% organizer.cells$id$Organizer.1 ~ "Organizer.1",
                                        barcode %in% organizer.cells$id$Organizer.2 ~ "Organizer.2",
                                        barcode %in% organizer.cells$id$Organizer.3 ~ "Organizer.3",
                                        barcode %in% organizer.cells$id$Organizer.4 ~ "Organizer.4"))
rownames(mf.embryo.organizer@meta.data) <- mf.embryo.organizer@meta.data$barcode
mf.embryo.organizer@meta.data$CellType.organizer[is.na(mf.embryo.organizer@meta.data$CellType.organizer)] <- "Others"
table(mf.embryo.organizer$CellType.organizer, mf.embryo.organizer$cell_type)
# make a stat
pdf(file.path(outdir, "Scatter_plot_to_show_monkey_organizer_4_populations_stage_ratio.pdf"), width = 5, height = 5)
PieDonut(subset(mf.embryo.organizer@meta.data, CellType.organizer %in% paste0("Organizer.", 1:4)),
         aes(pies = CellType.organizer, donuts = stage),
         r0 = 0.4, r1 = 0.6, r2 = 0.8, title = "Type vs Stage", showRatioThreshold = 0.001)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_monkey_organizer_1_populations_stage_ratio.pdf"), width = 5, height = 5)
PieDonut(mf.embryo.organizer@meta.data,
         aes(pies = CellType.organizer2, donuts = theiler_stage),
         r0 = 0.4, r1 = 0.6, r2 = 0.8, title = "Type vs Stage", showRatioThreshold = 0.001)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_monkey_organizer_1_populations_celltype_ratio.pdf"), width = 5, height = 5)
PieDonut(mf.embryo.organizer@meta.data,
         aes(pies = CellType.organizer2, donuts = cell_type),
         r0 = 0.4, r1 = 0.6, r2 = 0.8, title = "Type vs Stage", showRatioThreshold = 0.001)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_monkey_organizer_4_populations_cell_ratio.pdf"), width = 5, height = 5)
PieDonut(subset(mf.embryo.organizer@meta.data, CellType.organizer != "Others"),
         aes(pies = CellType.organizer, donuts = cell_type),
         r0 = 0.4, r1 = 0.6, r2 = 0.8, title = "Type vs Cell", showRatioThreshold = 0.001)
dev.off()
# find marker genes and GO
mf.embryo.organizer@meta.data %>%
  mutate(CellType.organizer = case_when(CellType.organizer == "Others" & cell_type == "APS" ~ "Others.APS",
                                        CellType.organizer == "Others" & cell_type == "DE" ~ "Others.DE",
                                        CellType.organizer == "Others" & cell_type == "ECT" ~ "Others.ECT",
                                        CellType.organizer == "Others" & cell_type == "Node" ~ "Others.Node",
                                        CellType.organizer == "Others" & cell_type == "VE" ~ "Others.VE",
                                        TRUE ~ CellType.organizer)) -> mf.embryo.organizer@meta.data
pdf(file.path(outdir, "Scatter_plot_to_show_monkey_organizer_4_populations_clustering_new.pdf"), width = 11, height = 4.5)
DimPlot(mf.embryo.organizer, reduction = "umap", group.by = c("stage", "CellType.organizer"),
        label = TRUE, repel = TRUE, pt.size = 2, cols = mk.col)
dev.off()
Idents(mf.embryo.organizer) <- mf.embryo.organizer$CellType.organizer
organizer.markers <- FindAllMarkers(mf.embryo.organizer, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
organizer.markers.go <- list()
for (i in as.character(unique(organizer.markers$cluster))) {
  pd.gene <- organizer.markers %>%
    mutate(cluster = as.character(cluster)) %>%
    filter(cluster == i) %>%
    top_n(n = 100, wt = avg_log2FC)
  prefix <- gsub("\\.", "_", i)
  organizer.markers.go[[prefix]] <- Pipe.GO(species = "human", genelist = pd.gene$gene,
                                            basename = prefix, genetype = "SYMBOL",
                                            res.out = file.path(outdir, "GO_Marker_gene_of_organizer"))
}
pd.gene <- organizer.markers %>%
  mutate(cluster = as.character(cluster)) %>%
  filter(cluster %in% paste0("Organizer.", 1:4)) %>%
  group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC)
organizer.markers.go$Organizer.all <- Pipe.GO(species = "human", genelist = pd.gene$gene,
                                              basename = "Organizer.all", genetype = "SYMBOL",
                                              res.out = file.path(outdir, "GO_Marker_gene_of_monkey_organizer"))
pd.gene <- organizer.markers %>%
  mutate(cluster = as.character(cluster)) %>%
  filter(cluster %in% grep("Others", as.character(unique(organizer.markers$cluster)), value = T)) %>%
  group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC)
organizer.markers.go$Others.all <- Pipe.GO(species = "human", genelist = pd.gene$gene,
                                           basename = "Others.all", genetype = "SYMBOL",
                                           res.out = file.path(outdir, "GO_Marker_gene_of_monkey_organizer"))
pd.go <- c("GO:0038092", "GO:0016055", "GO:0007369")
organizer.markers.go[[grep("Organizer", names(organizer.markers.go))[1]]]$GO.BP[pd.go,]
organizer.markers.go[[grep("Organizer", names(organizer.markers.go))[2]]]$GO.BP[pd.go,]
organizer.markers.go[[grep("Organizer", names(organizer.markers.go))[3]]]$GO.BP[pd.go,]
organizer.markers.go[[grep("Organizer", names(organizer.markers.go))[4]]]$GO.BP[pd.go,]
pdf(file.path(outdir, "GO_Marker_gene_of_monkey_organizer/final_organizer_go.pdf"), height = 4, width = 5)
organizer.markers.go$Organizer.all$GO.BP[pd.go,] %>% 
  mutate(
    GeneRatio.new = sapply(GeneRatio, function(x) eval(parse(text = as.character(x)))),
    Description = as.factor(Description), `-Log10(pvalue)` = -log10(pvalue)
  ) %>% 
  mutate(Description = fct_reorder(Description, `-Log10(pvalue)`)) %>%
  ggplot(aes(x = `-Log10(pvalue)`, y = Description)) +
  geom_bar(stat = "identity", fill = "#ffffff", linewidth = 1, color = "#B2DF8A") +
  geom_text(aes(x = 0.05, label = Description), hjust = 0, size = 3.25, stat = "identity") +
  geom_vline(xintercept = -log10(0.05), colour = "#000000", linetype = "longdash", size = 1, alpha = 0.5) +
  xlab("-log10(pvalue)") +
  theme_bw() +
  ggtitle("Monkey organizer") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )
dev.off()
# marker gene expression
pd.gene <- organizer.markers[-grep("^LOC[1-9].*", organizer.markers$gene),] %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC)
pdf(file.path(outdir, "Dot_plot_to_show_monkey_organizer_4_populations_markers_top30.pdf"), width = 35, height = 5)
cols <- colorRampPalette(c("#034896", "#52A7D5", "#FFFFFF", "#FF6262", "#B20000"))
Clustered_DotPlot(seurat_object = mf.embryo.organizer,
                  features = pd.gene$gene, k = 5, flip = T,
                  exp_color_min = -1, exp_color_max = 1,
                  colors_use_exp = cols(30))
dev.off()
pd.gene <- organizer.markers[-grep("^LOC[1-9].*", organizer.markers$gene),] %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
pdf(file.path(outdir, "Dot_plot_to_show_monkey_organizer_4_populations_markers_top20.pdf"), width = 30, height = 5)
cols <- colorRampPalette(c("#034896", "#52A7D5", "#FFFFFF", "#FF6262", "#B20000"))
Clustered_DotPlot(seurat_object = mf.embryo.organizer,
                  features = pd.gene$gene, k = 5, flip = T,
                  exp_color_min = -1, exp_color_max = 1,
                  colors_use_exp = cols(30))
dev.off()
pd.gene <- organizer.markers[-grep("^LOC[1-9].*", organizer.markers$gene),] %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
pdf(file.path(outdir, "Dot_plot_to_show_monkey_organizer_4_populations_markers_top10.pdf"), width = 20, height = 5)
cols <- colorRampPalette(c("#034896", "#52A7D5", "#FFFFFF", "#FF6262", "#B20000"))
Clustered_DotPlot(seurat_object = mf.embryo.organizer,
                  features = pd.gene$gene, k = 5, flip = T,
                  exp_color_min = -1, exp_color_max = 1,
                  colors_use_exp = cols(30))
dev.off()
# top50 ligands
pd.gene <- organizer.markers[-grep("^LOC[1-9].*", organizer.markers$gene),] %>%
  group_by(cluster) %>%
  filter(gene %in% unique(CellChatDB.human$interaction$ligand)) %>%
  top_n(n = 50, wt = avg_log2FC)
pdf(file.path(outdir, "Dot_plot_to_show_monkey_organizer_4_populations_markers_top50_ligands.pdf"), width = 15, height = 5)
cols <- colorRampPalette(c("#034896", "#52A7D5", "#FFFFFF", "#FF6262", "#B20000"))
Clustered_DotPlot(seurat_object = mf.embryo.organizer,
                  features = pd.gene$gene, k = 5, flip = T,
                  exp_color_min = -1, exp_color_max = 1,
                  colors_use_exp = cols(30))
dev.off()
# top50 LRs
pd.gene <- organizer.markers[-grep("^LOC[1-9].*", organizer.markers$gene),] %>%
  group_by(cluster) %>%
  filter(gene %in% lrs) %>%
  top_n(n = 50, wt = avg_log2FC)
pdf(file.path(outdir, "Dot_plot_to_show_monkey_organizer_4_populations_markers_top50_LRs.pdf"), width = 20, height = 5)
cols <- colorRampPalette(c("#034896", "#52A7D5", "#FFFFFF", "#FF6262", "#B20000"))
Clustered_DotPlot(seurat_object = mf.embryo.organizer,
                  features = pd.gene$gene, k = 5, flip = T,
                  exp_color_min = -1, exp_color_max = 1,
                  colors_use_exp = cols(30))
dev.off()


### >>> 4. ORG scRNA-seq data
table(sr.list$cl$Sample_Name)
sr.org <- subset(sr.list$cl, Sample_Name %in% c("ORG_D1", "ORG_D2", "ORG_D3"))
sr.org <- NormalizeData(sr.org) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
sr.org <- ScaleData(sr.org, verbose = FALSE, features = rownames(sr.org))
sr.org <- RunPCA(sr.org, npcs = 50, verbose = FALSE, features = VariableFeatures(sr.org))
ElbowPlot(sr.org, ndims = 30)
dim.n <- 10
sr.org <- RunUMAP(sr.org, reduction = "pca", dims = 1:dim.n, return.model = T) %>%
  RunTSNE(reduction = "pca", dims = 1:dim.n)
# plot clustering
DimPlot(sr.org, reduction = "umap", group.by = c("Sample_Name"),
        label = TRUE, repel = TRUE, pt.size = 0.5, cols = mk.col)
# remove batch effect
sr.org <- RunHarmony(sr.org, group.by.vars = "Sample_Name")
sr.org <- RunUMAP(sr.org, reduction = "harmony", dims = 1:dim.n, return.model = T)
sr.org <- FindNeighbors(sr.org, reduction = "pca", dims = 1:dims) %>%
  FindClusters(resolution = seq(0, 2, 0.2))
sr.org$Pect.mt <- round(sr.org$Pect.mt, 1)
pdf(file.path(outdir, "Tree_plot_to_show_ORG_clustering.pdf"), width = 10, height = 10)
clustree(sr.org, prefix = "RNA_snn_res.",
         node_label = "Pect.mt", node_label_aggr = "median") +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired")
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_ORG_clustering.pdf"), width = 11, height = 4.5)
DimPlot(sr.org, reduction = "umap", group.by = c("Sample_Name", "RNA_snn_res.1"),
        label = TRUE, repel = TRUE, pt.size = 0.5, cols = mk.col)
dev.off()
sr.org$Type <- "Negative"
sr.org$Type[GetAssayData(sr.org, slot = "count")["GFP", ] > 0] <- "Positive"
pdf(file.path(outdir, "Scatter_plot_to_show_ORG_clustering_with_GFP.pdf"), width = 11, height = 4.5)
DimPlot(sr.org, reduction = "umap", group.by = c("Sample_Name", "Type"),
        label = TRUE, repel = TRUE, pt.size = 0.5, cols = mk.col)
dev.off()
Plot_Density_Joint_Only(seurat_object = sr.org, features = c("OTX2", "FOXA2"), pt.size = 0.5)
pdf(file.path(outdir, "Scatter_plot_to_show_ORG_gene_expression_GSC_OTX2_FOXA2_FST_CER1_DKK1_HHEX_CHRD.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = sr.org,
                        features = c("GSC", "OTX2", "FOXA2", "FST", "CER1", "DKK1", "HHEX", "CHRD"), pt.size = 0.5)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_ORG_gene_expression_GSC_OTX2_FOXA2_FST_CER1_DKK1_HHEX.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = sr.org,
                        features = c("GSC", "OTX2", "FOXA2", "FST", "CER1", "DKK1", "HHEX"), pt.size = 0.5)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_ORG_gene_expression_GSC_OTX2_FOXA2_FST_CER1_DKK1.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = sr.org, features = c("GSC", "OTX2", "FOXA2", "FST", "CER1", "DKK1"), pt.size = 0.5)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_ORG_gene_expression_GSC_OTX2_FOXA2_FST_CER1.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = sr.org, features = c("GSC", "OTX2", "FOXA2", "FST", "CER1"), pt.size = 0.5)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_ORG_gene_expression_GSC_OTX2_FOXA2_FST.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = sr.org, features = c("GSC", "OTX2", "FOXA2", "FST"), pt.size = 0.5)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_ORG_gene_expression_GSC_OTX2_FOXA2.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = sr.org, features = c("GSC", "OTX2", "FOXA2"), pt.size = 0.5)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_ORG_gene_expression_OTX2_FOXA2.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = sr.org, features = c("OTX2", "FOXA2"), pt.size = 0.5)
dev.off()
# plot expression
sr.out <- file.path(outdir, "Marker_gene_of_ORG_cells")
dir.create(sr.out, recursive = T)
marker.gene <- list(Organizer = c("FST", "OTX2", "GSC", "EOMES", "FOXH1", "CER1", "MIXL1", "DKK1", "NODAL",
                                  "FOXA2", "LEFTY1", "LEFTY2", "T", "SP5", "NOG", "CHRD"),
                    #Ligands = intersect(unique(CellChatDB.human$interaction$ligand), rownames(sr.org)),
                    Buchong = c("NANOG", "HHEX", "CDH1", "POU3F1", "CDX2", "CDX1", "SOX17",
                                "GATA6", "GATA4", "LUM", "NID2", "POU5F1", "SOX2", "MYC",
                                "LIN28", "GDF3", "TDGF1", "SNAIL", "SNAIL2", "TBX6",
                                "MESP1", "LHX1", "IRX3", "HAND1", "FOXF1"),
                    FMH.marker  = c("LHX5","ZNF616","GBX2","RREB1","TEAD4","EN2","HAND1",
                                    "HESX1","FOXG1","NFATC4","PAX7","DBX1","PAX5","SMARCA5",
                                    "PAX6","IRX3","EBF1","HOXA3","MEIS1","WRNIP1","MAFB"),
                    FMH = c("PAX6","OTX2","EN1","PAX8","GBX2","EGR2","PAX3","DBX2","DBX1","IRX3","OLIG2","NKX2.2"),
                    Neural = c("SOX1","SOX3","PAX6","TUBB3","ELAVL3","OLIG2","NEUROG1","NEUROD1","NIKX2.2"),
                    NMP = c("LIMD1","MARK3","CIT","MAPK14","FAT4","LATS2","STK4","WWC1","SHANK2","SOX11","TEAD1","PJA2",
                            "TEAD2","VGLL5","SAV1","MAP2K3","TIAL1","NF2","STK3","LATS1","THBS1","TEAD3","AMOT","YAP1",
                            "AMOTL1","AMOTL2","NEK8","WTIP","DLG5","TEAD4"),
                    Amnion = c("GABRP","TGFB1","VTCN1","S100A10","GADD45G","ACTC1",
                               "IGFBP3","MCM5","CLDN10","TUBG1","UCP2","BIN3","GGA2","FGFR1","GALM"),
                    EPI = c("CDH1", "POU5F1", "SOX2", "NANOG", "LIN28AP1", "MYC", "KLF4", "TDGF1", "GDF3", "GFP"),
                    Embryonic = c("TBXT", "MESP1", "PDGFRA", "LHX1", "IRX3", "GATA6", "MYL7", "MSX1",
                                  "MSX2", "FOXF1", "HAND1", "HAND2", "BMP4", "MAB21L2", "CST1", "MEF2C",
                                  "PECAM1", "CDH5", "TEK", "TFAP2A", "TFAP2C", "GATA3", "DLX5", "LHX5",
                                  "PAX6", "SOX3", "PAX3", "PAX7"),
                    Extraembryonic = c("LUM", "NID2", "ANXA1", "POSTN", "GATA2", "KRT7", "VIM", "GATA4",
                                       "DCN", "TP63", "KRT18", "NR2F2", "ISL1", "HEY1", "CDH10", "CTSV",
                                       "ARID5B", "PLAGL1", "CREB3L1", "HOXA10", "HOXA11", "HOXA9",
                                       "HOXA13", "WNT6", "GABRP", "HEY1", "BST2"),
                    anterior = c("POU5F1","NANOG","TDGF1","NES","ID1","CDH1","SOX2","POU3F1","SOX3","IRX2",
                                 "CRABP1","PBX1","BTG1","TFAP2A","SFRP1","SFRP2","OTX2","FGFR1","HESX1","LHX5"),
                    medial = c("TBXT","TBX6","MESP1","EOMES","IRX3","LHX1","FOXC1","HOXB1","HOXA1","MEIS2","HOXB2",
                               "HOXB3","HOXA3","DLL1","GATA6","GATA4","FST","MIXL1","SP5","CER1","GFP"),
                    posterior = c("CDX2","CDX4","HOXB9","HOXC4","HOXC8","HOXC6","HOXB7","HOXB8","HOXB6","HOXB5",
                                  "VIM","SHH","FOXA2","CDH2","DUSP6","HOXA10","HOXA13","HAND1","HAND2","MSX1","GFP"),
                    megakaryocyte = c("GP1BB", "ITGA2B", "NFE2"),
                    erythroid = c("GATA1", "KLF1", "GYPB", "HBE1"),
                    myeloid.progenitor = c("CD36", "CSF1R", "LYVE1", "KIT", "CSF1R", "MYB", "SPI1", "CD34", "CD45", "CD52", "NFE2"),
                    wanghongmei = c("AFP","ANXA1","ANXA8","APOE","APOE ","CDH1","CDH2","CDX1","CDX2","CER1","CHRD",
                                    "COL1A1","COL3A1","DNMT3B","DPPA3","EOMES","FGF8","FGF8 HES7","FOLR2","FOXA2",
                                    "FOXC1","FOXC2","FOXF1","GATA1","GATA4","GATA6","GP1BB","HAND1","HBE1","HBM",
                                    "HBZ","HES7","HHEX","IGFBP5","IRX3","ISL1","KDR","KRT18","KRT8","LEFTY2","LFNG",
                                    "MEF2C","MIXL1","MPO","MYL7","NANOG","NANOS3","NKX2-5","NOTO","OSR1","OTX2","PAX3",
                                    "PAX3 ","PAX6","PAX7","PAX8","PECAM1","POU3F1","POU5F1","POU5FA","PPBP","PRDM14",
                                    "RUNX1","SIX1","SOX10","SOX17","SOX2","SOX2 ","SOX9","TBXT","TBX6","TCF15","TCF21",
                                    "TFAP2A","TFAP2C","TFC15","MESP1","TTR","UCHL1","VWF","WNT3A","WNT5B","WT1"),
                    pe.marker = c("GATA4", "GATA6", "SOX17", "SOX7", "PDGFRA", "PDGFA", "CLDN3", "FGFR4", "APOA1", "APOA1", "FN1", "S100A14",
                                  "AKR1D1","APOB","VTN","AFP","DPYS","VEPH1","SERPINE2","DENND2C","S100A11","COL4A1","HHEX","CXCR4","SHISAS2",
                                  "OTX2","SOX17","DKK1","SPOCK3","ACC117945.2","SIX3","CENPU","TK1","HMGN2","PCLAF","TMSB15A","CER1","EOMES","TTR","GJB1"))
for (type in names(marker.gene)[18]) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(sr.org))
}
for (type in names(marker.gene)[18]) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(sr.org, features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 0.5, slot = "data"))
    dev.off()
  }; rm(g)
}
pdf(file.path(tmp.out, "PE_marker_expression_level_in_ORG.pdf"), height = 20, width = 8)
VlnPlot(sr.org, group.by = "Sample_Name", features = marker.gene$pe.marker, stack = T, flip = T)
dev.off()
pdf(file.path(tmp.out, "PE_marker_expression_level_in_Hs_CS7.pdf"), height = 20, width = 15)
VlnPlot(hs.cs7, group.by = "CellType", features = marker.gene$pe.marker, stack = T, flip = T)
dev.off()
Vis.Annotation.Ratio(sr.org@meta.data, annotation = c("CellType.organizer", "Sample_Name"),
                     pd.title = "ORG",
                     pd.height = 10, pd.width = 10, res.out = file.path(outdir, "ORG_Cell_Type_Ratio"))
# correlation analysis with monkey scRNAseq data
expr.monkey <- AverageExpression(mf.embryo.organizer, assays = "RNA", slot = "data", group.by = "CellType.organizer")
expr.org <- AverageExpression(sr.org, assays = "RNA", slot = "data", group.by = "RNA_snn_res.1")
hvg.monkey <- VariableFeatures(mf.embryo.organizer)
hvg.org <- VariableFeatures(sr.org)
homo.gene <- data.frame(species.1 = intersect(rownames(mf.embryo.organizer), rownames(sr.org)),
                        species.2 = intersect(rownames(mf.embryo.organizer), rownames(sr.org)))
homo.gene[nrow(homo.gene)+1,] <- c("TBXT", "T")
CorrComparePlot(ExpressionTableSpecies1 = expr.monkey$RNA, DEgenesSpecies1 = hvg.monkey,
                ExpressionTableSpecies2 = expr.org$RNA, DEgenesSpecies2 = hvg.org,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(outdir, "Correlation_between_monkey_organizer_and_ORG_clusters"),
                Height = 10, Width = 10)
# extract cells
cell.coor <- data.frame(group = c("Organizer.1", "Organizer.2", "Organizer.3", "Organizer.4"),
                        xmin = c(-3.5, -3.5, -3, -2),
                        xmax = c(-2.5, -2.5, -2, -1),
                        ymin = c(8.5, 6, -8, -6),
                        ymax = c(9.5, 7, -7, -5))
organizer.cells <- ExtractCellByPos(object = sr.org, object.type = "seurat",
                                    dim.name = "umap", group.coor = cell.coor, group.name = "CellType.sub",
                                    pt.size = 0.5)
pdf(file.path(outdir, "Scatter_plot_to_show_ORG_organizer_4_populations_clustering.pdf"), width = 6, height = 4.5)
organizer.cells$plot
dev.off()
sr.org@meta.data %>%
  rownames_to_column(var = "barcode") %>%
  mutate(CellType.organizer = case_when(barcode %in% organizer.cells$id$Organizer.1 ~ "Organizer.1",
                                        barcode %in% organizer.cells$id$Organizer.2 ~ "Organizer.2",
                                        barcode %in% organizer.cells$id$Organizer.3 ~ "Organizer.3",
                                        barcode %in% organizer.cells$id$Organizer.4 ~ "Organizer.4")) -> sr.org@meta.data
rownames(sr.org@meta.data) <- sr.org@meta.data$barcode
sr.org@meta.data$CellType.organizer[is.na(sr.org@meta.data$CellType.organizer)] <- "Others"
table(sr.org@meta.data$CellType.organizer)
sr.org@meta.data %>%
  mutate(CellType.organizer = case_when(CellType.organizer == "Others" & RNA_snn_res.1 %in% c(0,1,5,9,10) ~ "Others.C1",
                                        CellType.organizer == "Others" & RNA_snn_res.1 %in% c(3,6,11) ~ "Others.C2",
                                        CellType.organizer == "Others" & RNA_snn_res.1 %in% c(2,4,7,8,12) ~ "Others.C3",
                                        TRUE ~ CellType.organizer)) -> sr.org@meta.data
pdf(file.path(outdir, "Scatter_plot_to_show_ORG_organizer_4_populations_clustering_new.pdf"), width = 11.5, height = 4.5)
DimPlot(sr.org, reduction = "umap", group.by = c("Sample_Name", "CellType.organizer"),
        label = TRUE, repel = TRUE, pt.size = 0.5, cols = mk.col)
dev.off()
# find markers
Idents(sr.org) <- sr.org$CellType.organizer
sr.org.markers <- FindAllMarkers(sr.org, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sr.org.markers.go <- list()
for (i in as.character(unique(sr.org.markers$cluster))) {
  pd.gene <- sr.org.markers %>%
    mutate(cluster = as.character(cluster)) %>%
    filter(cluster == i) %>%
    top_n(n = 100, wt = avg_log2FC)
  prefix <- gsub("\\.", "_", i)
  sr.org.markers.go[[prefix]] <- Pipe.GO(species = "human", genelist = pd.gene$gene,
                                         basename = prefix, genetype = "SYMBOL",
                                         res.out = file.path(outdir, "GO_Marker_gene_of_ORG_organizer"))
}
# marker gene expression
pd.gene <- sr.org.markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC)
pdf(file.path(outdir, "Dot_plot_to_show_ORG_organizer_4_populations_markers_top30.pdf"), width = 35, height = 5)
cols <- colorRampPalette(c("#034896", "#52A7D5", "#FFFFFF", "#FF6262", "#B20000"))
Clustered_DotPlot(seurat_object = sr.org,
                  features = pd.gene$gene, k = 5, flip = T,
                  exp_color_min = -1, exp_color_max = 1,
                  colors_use_exp = cols(30))
dev.off()
pd.gene <- sr.org.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
pdf(file.path(outdir, "Dot_plot_to_show_ORG_organizer_4_populations_markers_top20.pdf"), width = 30, height = 5)
cols <- colorRampPalette(c("#034896", "#52A7D5", "#FFFFFF", "#FF6262", "#B20000"))
Clustered_DotPlot(seurat_object = sr.org,
                  features = pd.gene$gene, k = 5, flip = T,
                  exp_color_min = -1, exp_color_max = 1,
                  colors_use_exp = cols(30))
dev.off()
pd.gene <- sr.org.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
pdf(file.path(outdir, "Dot_plot_to_show_ORG_organizer_4_populations_markers_top10.pdf"), width = 20, height = 5)
cols <- colorRampPalette(c("#034896", "#52A7D5", "#FFFFFF", "#FF6262", "#B20000"))
Clustered_DotPlot(seurat_object = sr.org,
                  features = pd.gene$gene, k = 5, flip = T,
                  exp_color_min = -1, exp_color_max = 1,
                  colors_use_exp = cols(30))
dev.off()
# top50 ligands
pd.gene <- sr.org.markers %>%
  group_by(cluster) %>%
  filter(gene %in% unique(CellChatDB.human$interaction$ligand)) %>%
  top_n(n = 50, wt = avg_log2FC)
pdf(file.path(outdir, "Dot_plot_to_show_ORG_organizer_4_populations_markers_top50_ligands.pdf"), width = 20, height = 5)
cols <- colorRampPalette(c("#034896", "#52A7D5", "#FFFFFF", "#FF6262", "#B20000"))
Clustered_DotPlot(seurat_object = sr.org,
                  features = pd.gene$gene, k = 5, flip = T,
                  exp_color_min = -1, exp_color_max = 1,
                  colors_use_exp = cols(30))
dev.off()
# top50 LRs
pd.gene <- sr.org.markers %>%
  group_by(cluster) %>%
  filter(gene %in% lrs) %>%
  top_n(n = 50, wt = avg_log2FC)
pdf(file.path(outdir, "Dot_plot_to_show_ORG_organizer_4_populations_markers_top50_LRs.pdf"), width = 25, height = 5)
cols <- colorRampPalette(c("#034896", "#52A7D5", "#FFFFFF", "#FF6262", "#B20000"))
Clustered_DotPlot(seurat_object = sr.org,
                  features = pd.gene$gene, k = 5, flip = T,
                  exp_color_min = -1, exp_color_max = 1,
                  colors_use_exp = cols(30))
dev.off()
# correlation analysis with monkey scRNAseq data
expr.monkey <- AverageExpression(mf.embryo.organizer, assays = "RNA", slot = "data", group.by = "CellType.organizer")
expr.org <- AverageExpression(sr.org, assays = "RNA", slot = "data", group.by = "CellType.organizer")
hvg.monkey <- VariableFeatures(mf.embryo.organizer)
hvg.org <- VariableFeatures(sr.org)
homo.gene <- data.frame(species.1 = intersect(rownames(mf.embryo.organizer), rownames(sr.org)),
                        species.2 = intersect(rownames(mf.embryo.organizer), rownames(sr.org)))
homo.gene[nrow(homo.gene)+1,] <- c("TBXT", "T")
CorrComparePlot(ExpressionTableSpecies1 = expr.monkey$RNA, DEgenesSpecies1 = hvg.monkey,
                ExpressionTableSpecies2 = expr.org$RNA, DEgenesSpecies2 = hvg.org,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(outdir, "Correlation_between_monkey_organizer_and_ORG_organizer"),
                Height = 10, Width = 10)


### >>> 5. GAS D1-D5
# clustering
pdf(file.path(outdir, "Scatter_plot_to_show_GAS_clustering.pdf"), width = 11, height = 4.5)
p1 <- DimPlot(gas.merge, reduction = "umap", group.by = "Sample_Name",
              label = TRUE, repel = TRUE, pt.size = 0.5, cols = c(ORG_D2 = "#4CAF49", EPI_D0 = "#B91122",
                                                                  GAS_D1 = "#9ECAE1", GAS_D3 = "#EB8838", GAS_D5 = "#08519C"))
p2 <- DimPlot(gas.merge, reduction = "umap", group.by = "CellType",
              label = TRUE, repel = TRUE, pt.size = 0.5, cols = gas.merge.col)
p1 + p2
dev.off()
# co-expression
pdf(file.path(outdir, "Scatter_plot_to_show_GAS_gene_expression_GSC_OTX2_FOXA2_FST_CER1_DKK1_HHEX_CHRD.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = gas.merge,
                        features = c("GSC", "OTX2", "FOXA2", "FST", "CER1", "DKK1", "HHEX", "CHRD"), pt.size = 0.5)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_GAS_gene_expression_GSC_OTX2_FOXA2_FST_CER1_DKK1_HHEX.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = gas.merge,
                        features = c("GSC", "OTX2", "FOXA2", "FST", "CER1", "DKK1", "HHEX"), pt.size = 0.5)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_GAS_gene_expression_GSC_OTX2_FOXA2_FST_CER1_DKK1.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = gas.merge, features = c("GSC", "OTX2", "FOXA2", "FST", "CER1", "DKK1"), pt.size = 0.5)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_GAS_gene_expression_GSC_OTX2_FOXA2_FST_CER1.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = gas.merge, features = c("GSC", "OTX2", "FOXA2", "FST", "CER1"), pt.size = 0.5)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_GAS_gene_expression_GSC_OTX2_FOXA2_FST.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = gas.merge, features = c("GSC", "OTX2", "FOXA2", "FST"), pt.size = 0.5)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_GAS_gene_expression_GSC_OTX2_FOXA2.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = gas.merge, features = c("GSC", "OTX2", "FOXA2"), pt.size = 0.5)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_GAS_gene_expression_OTX2_FOXA2.pdf"),
    height = 3.5, width = 5)
Plot_Density_Joint_Only(seurat_object = gas.merge, features = c("OTX2", "FOXA2"), pt.size = 0.5)
dev.off()
# plot expression
sr.out <- file.path(outdir, "Marker_gene_of_GAS")
dir.create(sr.out, recursive = T)
marker.gene <- list(Organizer = c("FST", "OTX2", "GSC", "EOMES", "FOXH1", "CER1", "MIXL1", "DKK1", "NODAL",
                                  "FOXA2", "LEFTY1", "LEFTY2", "TBXT", "SP5", "NOG", "CHRD"),
                    #Ligands = intersect(unique(CellChatDB.human$interaction$ligand), rownames(gas.merge)),
                    Buchong = c("NANOG", "HHEX", "CDH1", "POU3F1", "CDX2", "CDX1", "SOX17",
                                "GATA6", "GATA4", "LUM", "NID2", "POU5F1", "SOX2", "MYC",
                                "LIN28", "GDF3", "TDGF1", "SNAIL", "SNAIL2", "TBX6",
                                "MESP1", "LHX1", "IRX3", "HAND1", "FOXF1"),
                    FMH.marker  = c("LHX5","ZNF616","GBX2","RREB1","TEAD4","EN2","HAND1",
                                    "HESX1","FOXG1","NFATC4","PAX7","DBX1","PAX5","SMARCA5",
                                    "PAX6","IRX3","EBF1","HOXA3","MEIS1","WRNIP1","MAFB"),
                    FMH = c("PAX6","OTX2","EN1","PAX8","GBX2","EGR2","PAX3","DBX2","DBX1","IRX3","OLIG2","NKX2.2"),
                    Neural = c("SOX1","SOX3","PAX6","TUBB3","ELAVL3","OLIG2","NEUROG1","NEUROD1","NIKX2.2"),
                    NMP = c("LIMD1","MARK3","CIT","MAPK14","FAT4","LATS2","STK4","WWC1","SHANK2","SOX11","TEAD1","PJA2",
                            "TEAD2","VGLL5","SAV1","MAP2K3","TIAL1","NF2","STK3","LATS1","THBS1","TEAD3","AMOT","YAP1",
                            "AMOTL1","AMOTL2","NEK8","WTIP","DLG5","TEAD4"),
                    Amnion = c("GABRP","TGFB1","VTCN1","S100A10","GADD45G","ACTC1",
                               "IGFBP3","MCM5","CLDN10","TUBG1","UCP2","BIN3","GGA2","FGFR1","GALM"),
                    EPI = c("CDH1", "POU5F1", "SOX2", "NANOG", "LIN28AP1", "MYC", "KLF4", "TDGF1", "GDF3", "GFP"),
                    Embryonic = c("TBXT", "MESP1", "PDGFRA", "LHX1", "IRX3", "GATA6", "MYL7", "MSX1",
                                  "MSX2", "FOXF1", "HAND1", "HAND2", "BMP4", "MAB21L2", "CST1", "MEF2C",
                                  "PECAM1", "CDH5", "TEK", "TFAP2A", "TFAP2C", "GATA3", "DLX5", "LHX5",
                                  "PAX6", "SOX3", "PAX3", "PAX7"),
                    Extraembryonic = c("LUM", "NID2", "ANXA1", "POSTN", "GATA2", "KRT7", "VIM", "GATA4",
                                       "DCN", "TP63", "KRT18", "NR2F2", "ISL1", "HEY1", "CDH10", "CTSV",
                                       "ARID5B", "PLAGL1", "CREB3L1", "HOXA10", "HOXA11", "HOXA9",
                                       "HOXA13", "WNT6", "GABRP", "HEY1", "BST2"),
                    anterior = c("POU5F1","NANOG","TDGF1","NES","ID1","CDH1","SOX2","POU3F1","SOX3","IRX2",
                                 "CRABP1","PBX1","BTG1","TFAP2A","SFRP1","SFRP2","OTX2","FGFR1","HESX1","LHX5"),
                    medial = c("TBXT","TBX6","MESP1","EOMES","IRX3","LHX1","FOXC1","HOXB1","HOXA1","MEIS2","HOXB2",
                               "HOXB3","HOXA3","DLL1","GATA6","GATA4","FST","MIXL1","SP5","CER1","GFP"),
                    posterior = c("CDX2","CDX4","HOXB9","HOXC4","HOXC8","HOXC6","HOXB7","HOXB8","HOXB6","HOXB5",
                                  "VIM","SHH","FOXA2","CDH2","DUSP6","HOXA10","HOXA13","HAND1","HAND2","MSX1","GFP"),
                    megakaryocyte = c("GP1BB", "ITGA2B", "NFE2"),
                    erythroid = c("GATA1", "KLF1", "GYPB", "HBE1"),
                    myeloid.progenitor = c("CD36", "CSF1R", "LYVE1", "KIT", "CSF1R", "MYB", "SPI1", "CD34", "CD45", "CD52", "NFE2"),
                    wanghongmei = c("AFP","ANXA1","ANXA8","APOE","APOE ","CDH1","CDH2","CDX1","CDX2","CER1","CHRD",
                                    "COL1A1","COL3A1","DNMT3B","DPPA3","EOMES","FGF8","FGF8 HES7","FOLR2","FOXA2",
                                    "FOXC1","FOXC2","FOXF1","GATA1","GATA4","GATA6","GP1BB","HAND1","HBE1","HBM",
                                    "HBZ","HES7","HHEX","IGFBP5","IRX3","ISL1","KDR","KRT18","KRT8","LEFTY2","LFNG",
                                    "MEF2C","MIXL1","MPO","MYL7","NANOG","NANOS3","NKX2-5","NOTO","OSR1","OTX2","PAX3",
                                    "PAX3 ","PAX6","PAX7","PAX8","PECAM1","POU3F1","POU5F1","POU5FA","PPBP","PRDM14",
                                    "RUNX1","SIX1","SOX10","SOX17","SOX2","SOX2 ","SOX9","TBXT","TBX6","TCF15","TCF21",
                                    "TFAP2A","TFAP2C","TFC15","MESP1","TTR","UCHL1","VWF","WNT3A","WNT5B","WT1"))
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(gas.merge))
}
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(gas.merge, features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 0.5, slot = "data"))
    dev.off()
  }; rm(g)
}
# correlation analysis with monkey scRNAseq data
expr.monkey <- AverageExpression(mf.embryo.organizer, assays = "RNA", slot = "data", group.by = "CellType.organizer")
expr.gas <- AverageExpression(gas.merge, assays = "RNA", slot = "data", group.by = "CellType")
hvg.monkey <- VariableFeatures(mf.embryo.organizer)
hvg.gas <- VariableFeatures(gas.merge)
homo.gene <- data.frame(species.1 = intersect(rownames(mf.embryo.organizer), rownames(gas.merge)),
                        species.2 = intersect(rownames(mf.embryo.organizer), rownames(gas.merge)))
homo.gene[nrow(homo.gene)+1,] <- c("TBXT", "T")
CorrComparePlot(ExpressionTableSpecies1 = expr.monkey$RNA, DEgenesSpecies1 = hvg.monkey,
                ExpressionTableSpecies2 = expr.gas$RNA, DEgenesSpecies2 = hvg.gas,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(outdir, "Correlation_between_monkey_gasanizer_and_GAS_cell_types"),
                Height = 10, Width = 10)
# extract cells
cell.coor <- data.frame(group = c("Organizer.1", "Organizer.2"),
                        xmin = c(-12, -13.5),
                        xmax = c(-11, -12.5),
                        ymin = c(4.5, 2.5),
                        ymax = c(5.5, 3.5))
organizer.cells <- ExtractCellByPos(object = gas.merge, object.type = "seurat",
                                    dim.name = "umap", group.coor = cell.coor, group.name = "CellType.sub",
                                    pt.size = 0.5)
pdf(file.path(outdir, "Scatter_plot_to_show_GAS_organizer_2_populations_clustering.pdf"), width = 6, height = 4.5)
organizer.cells$plot
dev.off()
gas.merge@meta.data %>%
  rownames_to_column(var = "barcode") %>%
  mutate(CellType.organizer = case_when(barcode %in% organizer.cells$id$Organizer.1 ~ "Organizer.1",
                                        barcode %in% organizer.cells$id$Organizer.2 ~ "Organizer.2")) -> gas.merge@meta.data
rownames(gas.merge@meta.data) <- gas.merge@meta.data$barcode
gas.merge@meta.data$CellType.organizer[is.na(gas.merge@meta.data$CellType.organizer)] <- "Others"
table(gas.merge@meta.data$CellType.organizer)
gas.merge@meta.data %>%
  mutate(CellType.organizer = case_when(CellType.organizer == "Others" & CellType == "EndLCs" ~ "EndLCs",
                                        CellType.organizer == "Organizer.1" ~ "Organizer.1",
                                        CellType.organizer == "Organizer.2" ~ "Organizer.2",
                                        TRUE ~ CellType)) -> gas.merge@meta.data
pdf(file.path(outdir, "Scatter_plot_to_show_GAS_organizer_2_populations_clustering_new.pdf"), width = 11.5, height = 4.5)
p1 <- DimPlot(gas.merge, reduction = "umap", group.by = "Sample_Name",
              label = TRUE, repel = TRUE, pt.size = 0.5, cols = c(ORG_D2 = "#4CAF49", EPI_D0 = "#B91122",
                                                                  GAS_D1 = "#9ECAE1", GAS_D3 = "#EB8838", GAS_D5 = "#08519C"))
gas.merge.org.col <- c(gas.merge.col, Organizer.2 = "#9e9e9e", Organizer.1 = "#212121")
p2 <- DimPlot(gas.merge, reduction = "umap", group.by = "CellType.organizer",
              label = TRUE, repel = TRUE, pt.size = 0.5, cols = gas.merge.org.col)
p1 + p2
dev.off()
# stat of cell types
pdf(file.path(outdir, "Scatter_plot_to_show_GAS_organizer_2_populations_stage_ratio.pdf"), width = 5, height = 5)
PieDonut(subset(gas.merge@meta.data, CellType.organizer %in% c("Organizer.1", "Organizer.2")),
         aes(pies = CellType.organizer, donuts = Sample_Name),
         r0 = 0.4, r1 = 0.6, r2 = 0.8, title = "Type vs Stage", showRatioThreshold = 0.001)
dev.off()
pdf(file.path(outdir, "Scatter_plot_to_show_GAS_organizer_2_populations_cell_type_ratio.pdf"), width = 5, height = 5)
PieDonut(subset(gas.merge@meta.data, CellType.organizer %in% c("Organizer.1", "Organizer.2")),
         aes(pies = CellType.organizer, donuts = CellType),
         r0 = 0.4, r1 = 0.6, r2 = 0.8, title = "Type vs Stage", showRatioThreshold = 0.001)
dev.off()
# find markers
Idents(gas.merge) <- gas.merge$CellType.organizer
gas.merge.markers <- FindAllMarkers(gas.merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gas.merge.markers.go <- list()
for (i in as.character(unique(gas.merge.markers$cluster))) {
  pd.gene <- gas.merge.markers %>%
    mutate(cluster = as.character(cluster)) %>%
    filter(cluster == i) %>%
    top_n(n = 100, wt = avg_log2FC)
  prefix <- gsub("\\.", "_", i)
  gas.merge.markers.go[[prefix]] <- Pipe.GO(species = "human", genelist = pd.gene$gene,
                                            basename = prefix, genetype = "SYMBOL",
                                            res.out = file.path(outdir, "GO_Marker_gene_of_GAS_organizer"))
}
# marker gene expression
pd.gene <- gas.merge.markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC)
pdf(file.path(outdir, "Dot_plot_to_show_GAS_organizer_2_populations_markers_top30.pdf"), width = 50, height = 5)
cols <- colorRampPalette(c("#034896", "#52A7D5", "#FFFFFF", "#FF6262", "#B20000"))
Clustered_DotPlot(seurat_object = gas.merge,
                  features = pd.gene$gene, k = 13, flip = T,
                  exp_color_min = -1, exp_color_max = 1,
                  colors_use_exp = cols(30))
dev.off()
pd.gene <- gas.merge.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
pdf(file.path(outdir, "Dot_plot_to_show_GAS_organizer_2_populations_markers_top20.pdf"), width = 45, height = 5)
cols <- colorRampPalette(c("#034896", "#52A7D5", "#FFFFFF", "#FF6262", "#B20000"))
Clustered_DotPlot(seurat_object = gas.merge,
                  features = pd.gene$gene, k = 13, flip = T,
                  exp_color_min = -1, exp_color_max = 1,
                  colors_use_exp = cols(30))
dev.off()
pd.gene <- gas.merge.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
pdf(file.path(outdir, "Dot_plot_to_show_GAS_organizer_2_populations_markers_top10.pdf"), width = 35, height = 5)
cols <- colorRampPalette(c("#034896", "#52A7D5", "#FFFFFF", "#FF6262", "#B20000"))
Clustered_DotPlot(seurat_object = gas.merge,
                  features = pd.gene$gene, k = 13, flip = T,
                  exp_color_min = -1, exp_color_max = 1,
                  colors_use_exp = cols(30))
dev.off()
# top50 ligands
pd.gene <- gas.merge.markers %>%
  group_by(cluster) %>%
  filter(gene %in% unique(CellChatDB.human$interaction$ligand)) %>%
  top_n(n = 50, wt = avg_log2FC)
pdf(file.path(outdir, "Dot_plot_to_show_GAS_organizer_2_populations_markers_top50_ligands.pdf"), width = 30, height = 5)
cols <- colorRampPalette(c("#034896", "#52A7D5", "#FFFFFF", "#FF6262", "#B20000"))
Clustered_DotPlot(seurat_object = gas.merge,
                  features = pd.gene$gene, k = 13, flip = T,
                  exp_color_min = -1, exp_color_max = 1,
                  colors_use_exp = cols(30))
dev.off()
# top50 LRs
pd.gene <- gas.merge.markers %>%
  group_by(cluster) %>%
  filter(gene %in% lrs) %>%
  top_n(n = 50, wt = avg_log2FC)
pdf(file.path(outdir, "Dot_plot_to_show_GAS_organizer_2_populations_markers_top50_LRs.pdf"), width = 50, height = 5)
cols <- colorRampPalette(c("#034896", "#52A7D5", "#FFFFFF", "#FF6262", "#B20000"))
Clustered_DotPlot(seurat_object = gas.merge,
                  features = pd.gene$gene, k = 13, flip = T,
                  exp_color_min = -1, exp_color_max = 1,
                  colors_use_exp = cols(30))
dev.off()
# correlation analysis with monkey scRNAseq data
expr.monkey <- AverageExpression(mf.embryo.organizer, assays = "RNA", slot = "data", group.by = "CellType.organizer")
expr.gas <- AverageExpression(gas.merge, assays = "RNA", slot = "data", group.by = "CellType.organizer")
hvg.monkey <- VariableFeatures(mf.embryo.organizer)
hvg.gas <- VariableFeatures(gas.merge)
homo.gene <- data.frame(species.1 = intersect(rownames(mf.embryo.organizer), rownames(gas.merge)),
                        species.2 = intersect(rownames(mf.embryo.organizer), rownames(gas.merge)))
homo.gene[nrow(homo.gene)+1,] <- c("TBXT", "T")
CorrComparePlot(ExpressionTableSpecies1 = expr.monkey$RNA, DEgenesSpecies1 = hvg.monkey,
                ExpressionTableSpecies2 = expr.gas$RNA, DEgenesSpecies2 = hvg.gas,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(outdir, "Correlation_between_monkey_organizer_and_GAS_organizer"),
                Height = 10, Width = 10)
# subset GAS-D5
gas.merge.sub <- subset(gas.merge, Sample_Name == "GAS_D5")
gas.merge.sub$GSC <- GetAssayData(gas.merge.sub, "data")["GSC", ]
gas.merge.sub$OTX2 <- GetAssayData(gas.merge.sub, "data")["OTX2", ]
gas.merge.sub$FOXA2 <- GetAssayData(gas.merge.sub, "data")["FOXA2", ]
gas.merge.sub$CellType.organizer2 <- "Negative"
gas.merge.sub$CellType.organizer2[gas.merge.sub$GSC > 0 & gas.merge.sub$OTX2 > 0 & gas.merge.sub$FOXA2] <- "Positive"
pdf(file.path(outdir, "PieDonut_plot_to_show_GASD5_organizer_3positive_cell_type_vs_organizer_ratio.pdf"), width = 5, height = 5)
PieDonut(gas.merge.sub@meta.data,
         aes(pies = CellType, donuts = CellType.organizer2),
         r0 = 0.4, r1 = 0.6, r2 = 0.8, title = "CellType vs Organizer", showRatioThreshold = 0.0001)
dev.off()
gas.merge.sub$CellType.organizer2 <- "Negative"
gas.merge.sub$CellType.organizer2[gas.merge.sub$GSC > 0] <- "Positive"
pdf(file.path(outdir, "PieDonut_plot_to_show_GASD5_organizer_GSC_positive_cell_type_vs_organizer_ratio.pdf"), width = 5, height = 5)
PieDonut(gas.merge.sub@meta.data,
         aes(pies = CellType, donuts = CellType.organizer2),
         r0 = 0.4, r1 = 0.6, r2 = 0.8, title = "CellType vs Organizer", showRatioThreshold = 0.0001)
dev.off()
pdf(file.path(outdir, "PieDonut_plot_to_show_GASD5_organizer_GSC_positive_cell_type_vs_GFP_ratio.pdf"), width = 5, height = 5)
PieDonut(gas.merge.sub@meta.data,
         aes(pies = CellType, donuts = Type),
         r0 = 0.4, r1 = 0.6, r2 = 0.8, title = "CellType vs Organizer", showRatioThreshold = 0.0001)
dev.off()
gas.merge.sub <- subset(gas.merge.sub, CellType.organizer2 == "Positive")
gas.merge.sub <- NormalizeData(gas.merge.sub) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(gas.merge.sub), verbose = FALSE)
gas.merge.sub <- RunPCA(gas.merge.sub, features = VariableFeatures(gas.merge.sub), verbose = FALSE)
ElbowPlot(gas.merge.sub)
dev.off()
dim.n <- 5
gas.merge.sub <- RunUMAP(gas.merge.sub, dims = 1:dim.n)
pdf(file.path(outdir, "Scatter_to_show_GASD5_organizer_GSC_positive_cell_clustering.pdf"), height = 4.5, width = 5.5)
DimPlot(gas.merge.sub, group.by = "CellType", cols = pd.col.fix, pt.size = 2)
dev.off()
# plot expression
sr.out <- file.path(outdir, "Marker_gene_of_GASD5_GSC_positive_cell")
dir.create(sr.out, recursive = T)
DefaultAssay(gas.merge.sub) <- "RNA"
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(gas.merge.sub))
}
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(gas.merge.sub, features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 1, slot = "data", min.cutoff = 0.1))
    dev.off()
  }; rm(g)
}


### >>> 6. Correlation analysis between monkey and ORG
# monkey expression
mf.embryo.organizer$CellType.organizer2 <- "Negative"
mf.embryo.organizer$CellType.organizer2[mf.embryo.organizer$GSC > 0 & mf.embryo.organizer$OTX2 > 0 & mf.embryo.organizer$FOXA2] <- "Positive"
org.expr.monkey <- AverageExpression(mf.embryo.organizer, group.by = "CellType.organizer2")$RNA %>% as.data.frame()
colnames(org.expr.monkey) <- "Monkey"
rownames(org.expr.monkey)[grep("^T$", rownames(mf.embryo.organizer))] <- "TBXT"
org.expr.monkey["SOX17", ] <- mean(mk.gene.expr["SOX17", colnames(mf.embryo.organizer)])
# ORG expression
sr.org$GSC <- GetAssayData(sr.org, "data")["GSC", ]
sr.org$OTX2 <- GetAssayData(sr.org, "data")["OTX2", ]
sr.org$FOXA2 <- GetAssayData(sr.org, "data")["FOXA2", ]
sr.org$CellType.organizer2 <- "Negative"
sr.org$CellType.organizer2[sr.org$GSC > 0 & sr.org$OTX2 > 0 & sr.org$FOXA2] <- "Positive"
pdf(file.path(outdir, "PieDonut_plot_to_show_ORG_organizer_type_vs_sample_ratio.pdf"), width = 5, height = 5)
PieDonut(sr.org@meta.data,
         aes(pies = Sample_Name, donuts = CellType.organizer2),
         r0 = 0.4, r1 = 0.6, r2 = 0.8, title = "Type vs Stage", showRatioThreshold = 0.001)
dev.off()
org.expr.ours <- AverageExpression(subset(sr.org, CellType.organizer2 == "Positive"), group.by = "Sample_Name")$RNA %>% as.data.frame()
# plotting
common.gene <- Reduce(intersect, list(rownames(org.expr.monkey), rownames(org.expr.ours)))
pd <- data.frame(Monkey = org.expr.monkey[common.gene, ],
                 ORG.D1 = org.expr.ours[common.gene, ]$ORG_D1,
                 ORG.D2 = org.expr.ours[common.gene, ]$ORG_D2,
                 ORG.D3 = org.expr.ours[common.gene, ]$ORG_D3)
rownames(pd) <- common.gene
label.gene <- intersect(c("GSC", "FOXA2", "OTX2", "FST", "CER1", "DKK1", "HHEX",
                          "EOMES", "TBXT", "SOX17", "NOG", "FOXH1"), rownames(pd))
pd <- pd[c(setdiff(rownames(pd), label.gene), intersect(rownames(pd), label.gene)), ]
pd %>%
  mutate(type = c(rep("Others", nrow(pd)-length(label.gene)), rep("Markers", length(label.gene)))) %>%
  select("ORG.D3", "Monkey", "type") %>%
  rownames_to_column(var = "SYMBOL") %>%
  filter(SYMBOL %in% label.gene) -> pd.label
#pdf(file.path(outdir, "Correlation_analysis_of_organizer_between_ORGD1_and_Monkey_embryo_three_positive.pdf"), height = 4.5, width = 5.5)
#pdf(file.path(outdir, "Correlation_analysis_of_organizer_between_ORGD2_and_Monkey_embryo_three_positive.pdf"), height = 4.5, width = 5.5)
#pdf(file.path(outdir, "Correlation_analysis_of_organizer_between_ORGD3_and_Monkey_embryo_three_positive.pdf"), height = 4.5, width = 5.5)
pd %>%
  mutate(type = c(rep("Others", nrow(pd)-length(label.gene)), rep("Markers", length(label.gene)))) %>%
  select("ORG.D3", "Monkey", "type") %>%
  ggplot(aes(x = log2(ORG.D3 + 0.1), y = log2(Monkey + 0.1))) +
  geom_point(aes(color = type, size = type), shape = 16) +
  scale_color_manual(values = c("Others" = "#000000", "Markers" = "#E91E78")) +
  scale_size_manual(values = c("Others" = 1, "Markers" = 3)) +
  geom_smooth(method = "lm", se = F) +
  stat_cor(method = "pearson") +
  labs(x = "ORG", y = "Monkey embryo") +
  geom_label_repel(aes(label = SYMBOL), data = pd.label,
                   box.padding = 0.2, point.padding = 0.2, segment.color = "#130f40",
                   label.size = 0.2, fill = NA, size = 3) +
  theme_bw() +
  theme(
    axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
    axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
    axis.text.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
    axis.text.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
    panel.grid = element_blank()
  )
dev.off()
# integrate datasets
sr.list.org <- list()
sr.list.org$data1 <- mf.embryo.organizer
sr.list.org$data1$CellType <- sr.list.org$data1$cell_type
sr.list.org$data1$Sample <- sr.list.org$data1$sample
sr.list.org$data1$Datasets <- "Monkey"
sr.list.org$data1@meta.data <- sr.list.org$data1@meta.data[, c("Datasets", "Sample", "CellType", "CellType.organizer")]
sr.list.org$data2 <- sr.org
sr.list.org$data2$Sample <- sr.list.org$data2$Sample_Name
sr.list.org$data2$Datasets <- "ORG"
sr.list.org$data2@meta.data <- sr.list.org$data2@meta.data[, c("Datasets", "Sample", "CellType.organizer")]
# merge
sr.list.org <- Reduce(function(x, y) merge(x, y), sr.list.org)
table(sr.list.org$Datasets)
sr.list.org <- NormalizeData(sr.list.org) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sr.list.org), verbose = FALSE)
sr.list.org <- RunPCA(sr.list.org, features = VariableFeatures(sr.list.org), verbose = FALSE)
ElbowPlot(sr.list.org)
dev.off()
dim.n <- 5
sr.list.org <- RunUMAP(sr.list.org, dims = 1:dim.n)
sr.list.org <- RunHarmony(sr.list.org, group.by.vars = c("Datasets", "Sample"))
sr.list.org <- RunUMAP(sr.list.org, reduction = "harmony", dims = 1:dim.n)
pdf(file.path(outdir, "Integrative_analysis_between_ORGD1_and_Monkey_embryo_three_positive.pdf"), height = 4.5, width = 11)
DimPlot(sr.list.org, group.by = c("Datasets", "Sample"), cols = pd.col)
dev.off()
# plot expression
sr.out <- file.path(outdir, "Marker_gene_of_ORG_and_Monkey")
dir.create(sr.out, recursive = T)
marker.gene <- list(Organizer = c("FST", "OTX2", "GSC", "EOMES", "FOXH1", "CER1", "MIXL1", "DKK1", "NODAL",
                                  "FOXA2", "LEFTY1", "LEFTY2", "T", "SP5", "NOG", "CHRD"),
                    Buchong = c("NANOG", "HHEX", "CDH1", "POU3F1", "CDX2", "CDX1", "SOX17",
                                "GATA6", "GATA4", "LUM", "NID2", "POU5F1", "SOX2", "MYC",
                                "LIN28", "GDF3", "TDGF1", "SNAIL", "SNAIL2", "TBX6",
                                "MESP1", "LHX1", "IRX3", "HAND1", "FOXF1"),
                    FMH.marker  = c("LHX5","ZNF616","GBX2","RREB1","TEAD4","EN2","HAND1",
                                    "HESX1","FOXG1","NFATC4","PAX7","DBX1","PAX5","SMARCA5",
                                    "PAX6","IRX3","EBF1","HOXA3","MEIS1","WRNIP1","MAFB"),
                    FMH = c("PAX6","OTX2","EN1","PAX8","GBX2","EGR2","PAX3","DBX2","DBX1","IRX3","OLIG2","NKX2.2"),
                    Neural = c("SOX1","SOX3","PAX6","TUBB3","ELAVL3","OLIG2","NEUROG1","NEUROD1","NIKX2.2"),
                    NMP = unique(c("LIMD1","MARK3","CIT","MAPK14","FAT4","LATS2","STK4","WWC1","SHANK2","SOX11","TEAD1","PJA2",
                                   "TEAD2","VGLL5","SAV1","MAP2K3","TIAL1","NF2","STK3","LATS1","THBS1","TEAD3","AMOT","YAP1",
                                   "AMOTL1","AMOTL2","NEK8","WTIP","DLG5","TEAD4",
                                   "NKX1-2","SP5","FGF17","WNT3A","CDX1","CDX2","FGF8","HES7","CDH2",
                                   "TBX6","FGF15","WNT5B","WNT10A","DUSP6","AXIN2","HOXA5","HOXB8","HOXA13",
                                   "HOXC13","POU5F1","SALL4","TET1","TET2","TET3","MED12","WNT8C","WNT8A",
                                   "CTNNB1","SP8","AXIN2","TCF1","LEF1","WNT3","VANGL2","WNT5A","WNT11","FGF4",
                                   "FGFR1","ALDH1A2","RARB","DELTA1","GFD8","GDF11","BMP4","POU3F1","UTF1","FOXA2",
                                   "CER1","LEFTY1","LEFTY2","MESP1","TBX3","EOMES","GATA6","MEIS1","OTX2","RSPO3",
                                   "HOXA1","IDL1","HOXA1","MIXL1","CER1","IDI1","HMGCS1","UCHL1","PLP1","KISS1","HES6",
                                   "LPAR6","CTGF","LMO2","CCNB2","EBPL","TBXT","SOX2","SALL1","BASP1","VCAN","FASN",
                                   "ID3","AURKA","MGST1","ALPL","BCAT1","SFRP2","MSX1","HOXA9","DIT4","PENK","EGR1",
                                   "MLLT3","FOSB","GFP")),
                    Amnion = c("GABRP","TGFB1","VTCN1","S100A10","GADD45G","ACTC1",
                               "IGFBP3","MCM5","CLDN10","TUBG1","UCP2","BIN3","GGA2","FGFR1","GALM"),
                    EPI = c("CDH1", "POU5F1", "SOX2", "NANOG", "LIN28AP1", "MYC", "KLF4", "TDGF1", "GDF3", "GFP"),
                    Embryonic = c("TBXT", "MESP1", "PDGFRA", "LHX1", "IRX3", "GATA6", "MYL7", "MSX1",
                                  "MSX2", "FOXF1", "HAND1", "HAND2", "BMP4", "MAB21L2", "CST1", "MEF2C",
                                  "PECAM1", "CDH5", "TEK", "TFAP2A", "TFAP2C", "GATA3", "DLX5", "LHX5",
                                  "PAX6", "SOX3", "PAX3", "PAX7"),
                    Extraembryonic = c("LUM", "NID2", "ANXA1", "POSTN", "GATA2", "KRT7", "VIM", "GATA4",
                                       "DCN", "TP63", "KRT18", "NR2F2", "ISL1", "HEY1", "CDH10", "CTSV",
                                       "ARID5B", "PLAGL1", "CREB3L1", "HOXA10", "HOXA11", "HOXA9",
                                       "HOXA13", "WNT6", "GABRP", "HEY1", "BST2"),
                    anterior = c("POU5F1","NANOG","TDGF1","NES","ID1","CDH1","SOX2","POU3F1","SOX3","IRX2",
                                 "CRABP1","PBX1","BTG1","TFAP2A","SFRP1","SFRP2","OTX2","FGFR1","HESX1","LHX5"),
                    medial = c("TBXT","TBX6","MESP1","EOMES","IRX3","LHX1","FOXC1","HOXB1","HOXA1","MEIS2","HOXB2",
                               "HOXB3","HOXA3","DLL1","GATA6","GATA4","FST","MIXL1","SP5","CER1","GFP"),
                    posterior = c("CDX2","CDX4","HOXB9","HOXC4","HOXC8","HOXC6","HOXB7","HOXB8","HOXB6","HOXB5",
                                  "VIM","SHH","FOXA2","CDH2","DUSP6","HOXA10","HOXA13","HAND1","HAND2","MSX1","GFP"),
                    megakaryocyte = c("GP1BB", "ITGA2B", "NFE2"),
                    erythroid = c("GATA1", "KLF1", "GYPB", "HBE1"),
                    myeloid.progenitor = c("CD36", "CSF1R", "LYVE1", "KIT", "CSF1R", "MYB", "SPI1", "CD34", "CD45", "CD52", "NFE2"),
                    wanghongmei = c("AFP","ANXA1","ANXA8","APOE","APOE ","CDH1","CDH2","CDX1","CDX2","CER1","CHRD",
                                    "COL1A1","COL3A1","DNMT3B","DPPA3","EOMES","FGF8","FGF8 HES7","FOLR2","FOXA2",
                                    "FOXC1","FOXC2","FOXF1","GATA1","GATA4","GATA6","GP1BB","HAND1","HBE1","HBM",
                                    "HBZ","HES7","HHEX","IGFBP5","IRX3","ISL1","KDR","KRT18","KRT8","LEFTY2","LFNG",
                                    "MEF2C","MIXL1","MPO","MYL7","NANOG","NANOS3","NKX2-5","NOTO","OSR1","OTX2","PAX3",
                                    "PAX3 ","PAX6","PAX7","PAX8","PECAM1","POU3F1","POU5F1","POU5FA","PPBP","PRDM14",
                                    "RUNX1","SIX1","SOX10","SOX17","SOX2","SOX2 ","SOX9","TBXT","TBX6","TCF15","TCF21",
                                    "TFAP2A","TFAP2C","TFC15","MESP1","TTR","UCHL1","VWF","WNT3A","WNT5B","WT1"),
                    neuron_NMP = c("CYSTM1","EM6","GRHL3","PAX6","PAX3","PAX2"),
                    meso_NMP = c("IRX3","HOXB9","CDX4","EPHA5","HES3","HBB-BH1","LMO2","ANXA5",
                                 "PMP22","GBX2","HOXB1","APELA","DNMT3B","SNAIL","SNAIL2","TTR"),
                    GUT = c("SHH","ISL1","IHH","ZG16","CDX1","URAD","TACSTD2","CDX2","TNNC1","COL9A2","MSX2","TBX3",
                            "CDH2","ADD3","EPCAM","CTHRC1","EBPL","CTNNBL1","FN1","PENK","FDX1","CXCL12","RBP1",
                            "CTSV","FOXJ1","C15H9ORF116","EGFL6","BMP7","S100A16","S100A13","TTR","AFP","GJB1",
                            "HHEX","SHISA2","CXCR4","NKX2.1","PPY","IRX1","APOA2","AKR1D1","APOB","VTN","DPYS",
                            "VEPH1","SERPINE2","DENND2C","S100A11","COL4A1","SPOCK3","AC117945.2","SIX3","CENPU",
                            "TK1","HMGN2","PCLAF","TMSB15A","SFRP5","APOM","APOC1","MEST","TEKT1","S100A1","VEGFA",
                            "LAMA1","ZIM2","LARP7","GJA1","FOS","UBE2C","COL18A1","IGDCC3","CLDN4","SLC39A2",
                            "LPAR6","BAMBI","CDKN1C","NKX2-5","HHEX","IRX2","OSR1"),
                    ExE.meso = c("PENK","APLNR","HAS2","HAPLN1","TMEM88","PMP22","BAMBI","PITX1","CDH11","KDR"),
                    Cardiac.meso = c("NKX2.5","TBX5","TNNI1","TNNT2","MYH10","MAB21I2","HOPX","TGFBI"))
DefaultAssay(sr.list.org) <- "RNA"
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(sr.list.org))
}
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(sr.list.org, features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 0.5, slot = "data", min.cutoff = 0.1))
    dev.off()
  }; rm(g)
}
# marker genes and GO
table(sr.list.org$Sample)
sr.list.org@meta.data <- sr.list.org@meta.data %>%
  mutate(CellType.go = case_when(Sample %in% c("CS8-e1", "CS8-e2", "CS9-e1", "CS9-e2") ~ "Monkey",
                                 TRUE ~ Sample))
Idents(sr.list.org) <- sr.list.org$CellType.go
organizer.markers.merge <- FindAllMarkers(sr.list.org, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
organizer.markers.merge.go <- list()
for (i in as.character(unique(organizer.markers.merge$cluster))) {
  pd.gene <- organizer.markers.merge %>%
    mutate(cluster = as.character(cluster)) %>%
    filter(cluster == i) %>%
    top_n(n = 100, wt = avg_log2FC)
  organizer.markers.merge.go[[i]] <- Pipe.GO(species = "human", genelist = pd.gene$gene,
                                             basename = i, genetype = "SYMBOL",
                                             res.out = file.path(outdir, "GO_Marker_gene_of_organizer_merged"))
}


### ===========================
### 11th part: FOXH1 KO RNA-seq ----
### ===========================


### >>> 1. Setting output directory
outdir <- file.path(getwd(), "graphs/FOXH1_KO_RNAseq")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}


### >>> 2. Load data
# gene count
gene.count <- list()
gene.count$foxh1.ko <- read.table("../bulk_rnaseq/analysis/results/featurecounts/gene/all_samples_gene_count_name_matrix.txt",
                                  header = T, sep = "\t", row.names = 1)
colnames(gene.count$foxh1.ko)[-1:-5] <- c(paste0("FOXH1.KO_rep", 1:3), paste0("FOXH1.Ctrl_rep", 1:3))
# cell meta.data
cell.meta <- list()
cell.meta$foxh1.ko <- data.frame(SampleID = colnames(gene.count$foxh1.ko)[-1:-5],
                                 SampleGroup = str_split_fixed(colnames(gene.count$foxh1.ko)[-1:-5], "_", 2)[, 1],
                                 Batch = c(rep("A", 5), "B"),
                                 row.names = colnames(gene.count$foxh1.ko)[-1:-5])


### >>> 3. Normalize data
gene.tpm <- list()
gene.tpm$foxh1.ko <- CountToTpm(count = gene.count$foxh1.ko[, -1:-5], length = gene.count$foxh1.ko$Length)


### >>> 4. PCA analysis
pdf(file.path(outdir, "PCA_plot_of_all_sample_colored_by_sample_with_top2000_high_variables_genes.pdf"), height = 5, width = 6.5)
PCAplot(data = gene.tpm$foxh1.ko, geneset = HVG.Topn(gene.tpm$foxh1.ko, 2000),
        label = cell.meta$foxh1.ko$SampleID, color.by = cell.meta$foxh1.ko$SampleGroup,
        pt.size = 3, title = "2D-plot of PCA")
dev.off()


### >>> 5. Correlation analysis
CorClust(expr = gene.tpm$foxh1.ko, cor.method = "spearman",
         geneset = HVG.Topn(gene.tpm$foxh1.ko, 2000),
         low.expr = 0,
         res.out = file.path(outdir, "Correlation"),
         prefix = "All_samples_spearman_correlation")
CorClust(expr = gene.tpm$foxh1.ko, cor.method = "pearson",
         geneset = HVG.Topn(gene.tpm$foxh1.ko, 2000),
         low.expr = 0,
         res.out = file.path(outdir, "Correlation"),
         prefix = "All_samples_pearson_correlation")


### >>> 6. Differential expression analysis
deg <- list()
# KO vs Ctrl
deg[["foxh1.ko.vs.wt"]] <- edgeR.Ge.ET.LRT(count = gene.count$foxh1.ko[, -1:-5],
                                           meta = cell.meta$foxh1.ko,
                                           g1 = "FOXH1.KO", g2 = "FOXH1.Ctrl",
                                           lfc = 1, sig = 0.05, dir = file.path(outdir, "DEGs"),
                                           nor = gene.tpm$foxh1.ko, with.rep = T,
                                           filter.data = T, lowest.cout = 5)


### >>> 7. Volcano plot
# all DEGs
pdf(file.path(outdir, "Volcano_plot_of_all_DEGs.pdf"), height = 5, width = 6)
VisDEG.volcano(deg.data = deg[["foxh1.ko.vs.wt"]]$res, geneset = NULL,
               p.col = "PValue", lfc.col = "logFC", sig = 0.05, lfc = 1,
               title = i, pt.size = 2, up.col = "#B71C1C", down.col = "#01579B")
dev.off()
# custome gene
pd.gene <- c("OTX2", "FST", "GSC", "EOMES", "CER1", "DKK1", "NODAL", "FOXA2", 
             "NOGGIN", "HHEX", "TBXT", "TBX6", "MESP1", "LHX1", "IRX3", "SOX17", 
             "LUM", "POSTN", "NID2", "MIXL1", "HAND1")
pdf(file.path(outdir, "Heatmap_of_custom_gene.pdf"), height = 8, width = 10)
VisDeseq2.heatmap(deg.data = deg[["foxh1.ko.vs.wt"]]$res,
                  meta = cell.meta$foxh1.ko, geneset = pd.gene,
                  pt.mode = "genes", pt.value = "zscore",
                  show.rownames = T, show.colnames = T,
                  up.col = "#B71C1C", down.col = "#01579B")
dev.off()


### >>> 8. GO
go <- list()
for (i in names(deg)) {
  go[[paste0(i, ".up")]] <- Pipe.GO(species = "human", genelist = deg[[i]]$up.sig$SYMBOL,
                                    basename = paste0(gsub("\\.", "_", i), "_up"), genetype = "SYMBOL",
                                    res.out = file.path(outdir, paste0("GO/", gsub("\\.", "_", i))))
  go[[paste0(i, ".down")]] <- Pipe.GO(species = "human", genelist = deg[[i]]$down.sig$SYMBOL,
                                      basename = paste0(gsub("\\.", "_", i), "_down"), genetype = "SYMBOL",
                                      res.out = file.path(outdir, paste0("GO/", gsub("\\.", "_", i))))
}
go.terms <- c("GO:0035567", "GO:0030509", "GO:0007369", "GO:0007389", "GO:0009953",
              "GO:0009952", "GO:0007389", "GO:0009798", "GO:0001708", "GO:0009948", "GO:0007492")
pdf(file.path(outdir, "Barplot_of_custom_go_terms.pdf"), height = 6, width = 6)
Pipe.GO.replot(go.data = go$foxh1.ko.vs.wt.down$GO.BP, go.terms = go.terms)
dev.off()
Pipe.GO.replot <- function(go.data, go.terms, bar.cols = "#000000") {
  library("tidyr")
  library("dplyr")
  library("ggplot2")
  
  pd <- go.data[unique(go.terms), ] %>%
    mutate(
      GeneRatio.new = sapply(GeneRatio, function(x) eval(parse(text = as.character(x)))),
      Description = as.factor(Description), `-Log10(pvalue)` = -log10(pvalue)
    ) %>%
    mutate(Description = fct_reorder(Description, `-Log10(pvalue)`))
  p <- pd %>%
    ggplot(aes(x = `-Log10(pvalue)`, y = Description)) +
    geom_bar(stat = "identity", fill = "#ffffff", linewidth = 1, color = bar.cols) +
    geom_text(aes(x = 0.05, label = Description), hjust = 0, size = 3.25, stat = "identity") +
    geom_vline(xintercept = -log10(0.05), colour = "#000000", linetype = "longdash", size = 1, alpha = 0.5) +
    xlab("-log10(pvalue)") +
    theme_bw() +
    ggtitle("Customed Terms") +
    theme(
      axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
      axis.text.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
      axis.text.y = element_blank(), axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank()
    )
  return(p)
}


### >>> 9. GSEA analysis
gsea <- list()
for (i in names(deg)) {
  gsea[[i]] <- Pipe.GSEA(deg.obj = deg[[i]]$all.sig, deg.type = "edger",
                         lfc = 1, sig = 0.05, species = "human",
                         basename = gsub("\\.", "_", i),
                         genetype = "SYMBOL", gene.col = "SYMBOL",
                         outdir = file.path(outdir, paste0("GSEA/", gsub("\\.", "_", i))))
}



### =================================
### 12th part: Deconvolution analysis ----
### =================================

### >>> 1. Setting output dir
outdir <- file.path(getwd(), "graphs/BSeq-sc")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}


### >>> 2. Prepare data
# bulk RNAseq data
bulk.expr <- gene.tpm$foxh1.ko
bulk.expr <- Reduce(cbind, list(gene.tpm$foxh1.ko,
                                gene.tpm$foxh1.ko + 0.1, gene.tpm$foxh1.ko + 0.2,
                                gene.tpm$foxh1.ko - 0.1, gene.tpm$foxh1.ko - 0.2))
bulk.expr <- bulk.expr[, c(grep("Ctrl", colnames(bulk.expr)), grep("KO", colnames(bulk.expr)))]
bulk.meta <- data.frame(SampleID = colnames(bulk.expr),
                        SampleGroup = str_split_fixed(colnames(bulk.expr), "_", 2)[,1],
                        row.names = colnames(bulk.expr))
# cell markers
library("future")
multisession(workers = 8)
markers <- FindAllMarkers(sr.pdt$gasd5.to.hs, only.pos = T, logfc.threshold = log2(0.25))
markers %>%
  group_by(cluster) %>%
  top_n(10, wt = avg_log2FC) %>%
  pull(gene, cluster) -> cell.markers
cell.markers <- split(cell.markers, names(cell.markers))
cell.markers <- lapply(cell.markers, unname)
cell.markers <- list(EPI = c("SOX2", "POU5F1", "OTX2", "CDH1"),
                     ECT = c("DLX5", "TFAP2A", "TFAP2C", "GATA3"),
                     PS = c("TBXT", "MIXL1", "SP5", "FST", "EOMES"),
                     NM = c("TBX6", "MESP1", "EOMES"),
                     EM = c("LHX1", "IRX3", "GATA6"),
                     AM = c("FOXF1", "HAND1", "HAND2"),
                     ExM = c("POSTN", "LUM", "NID2", "HAND1"),
                     End = c("SOX17", "FOXA2", "GSC", "CST1", "GATA6"),
                     HEP = c("MEF2C", "PECAM1", "CDH5", "TEK"))
test <- Pipe.BSEQsc(
  bulk.expr = bulk.expr, bulk.meta = bulk.meta,
  sc.obj = seurat.obj, cell.markers = cell.markers,
  cluster.by = "CellType", sample.by = "Sample_Name", group.by = "SampleGroup",
  cibersort = "/home/yhw/document/BSEQ-sc/CIBERSORTx/CIBERSORT_v1.04.R",
  remove.regex = "rep.\\..", res.out = outdir
)



### ==========================================
### 13th part: Comparison with other embryoids ----
### ==========================================

### >>> 1. Setting output dir
outdir <- file.path(getwd(), "graphs/embryoids")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}


### >>> 2. Prepare data
# load data
load("/home/yhw/document/public_data/GSM5621089_Integrated_24_36_48h.Rd")
sr.gse144897 <- readRDS("/home/yhw/document/public_data/GSE144897_GSE169074_human_gastruloids/R/CodeData/GSE144897_GSE169074_human_gastruloids.rds")
# rename
sr.gse185643 <- allsamples_PCA_cluster_dims20_res0.5
rm(allsamples_PCA_cluster_dims20_res0.5)
# seurat pipeline
sr.gse185643 <- NormalizeData(sr.gse185643) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sr.gse185643))
sr.gse185643$AuthorLabel <- as.character(Idents(sr.gse185643))
sr.gse144897 <- NormalizeData(sr.gse144897) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sr.gse144897))


### >>> 3. Integrating our data and GSE185643
common.gene <- Reduce(intersect, list(rownames(sr.pdt$gasd5.to.hs), rownames(hs.cs7), rownames(sr.gse185643), rownames(sr.gse144897)))
# genes filtering
sr.integrate <- list()

sr.integrate$gasd5 <- sr.pdt$gasd5.to.hs[rownames(sr.pdt$gasd5.to.hs) %in% common.gene, ]
sr.integrate$gasd5$Datasets <- "GAS.D5"

sr.integrate$gse185643.24h <- sr.gse185643[rownames(sr.gse185643) %in% common.gene, sr.gse185643$old.ident == "hESC_24h"]
sr.integrate$gse185643.24h$Datasets <- "Embryoids.24h"

sr.integrate$gse185643.36h <- sr.gse185643[rownames(sr.gse185643) %in% common.gene, sr.gse185643$old.ident == "hESC_36h"]
sr.integrate$gse185643.36h$Datasets <- "Embryoids.36h"

sr.integrate$gse185643.48h <- sr.gse185643[rownames(sr.gse185643) %in% common.gene, sr.gse185643$old.ident == "hESC_48h"]
sr.integrate$gse185643.48h$Datasets <- "Embryoids.48h"

sr.integrate$gse144897 <- sr.gse144897[rownames(sr.gse144897) %in% common.gene, ]
sr.integrate$gse144897$Datasets <- "Gastruloids"

sr.integrate$hs.cs7 <- sr.pdt.hs.mk$hs.cs7[rownames(sr.pdt.hs.mk$hs.cs7) %in% common.gene, ]
sr.integrate$hs.cs7$Datasets <- "Hs.CS7"
# annotate cells (24h)
anchors <- FindTransferAnchors(reference = sr.integrate$hs.cs7, query = sr.integrate$gse185643.24h,
                               dims = 1:25, reference.reduction = "pca")
sr.integrate$gse185643.24h <- MapQuery(anchorset = anchors, reference = sr.integrate$hs.cs7,
                                       query = sr.integrate$gse185643.24h,
                                       refdata = list(celltype = "CellType"),
                                       reference.reduction = "pca", reduction.model = "umap")
sr.integrate$gse185643.24h$CellType <- sr.integrate$gse185643.24h$predicted.celltype
# annotate cells (36h)
anchors <- FindTransferAnchors(reference = sr.integrate$hs.cs7, query = sr.integrate$gse185643.36h,
                               dims = 1:25, reference.reduction = "pca")
sr.integrate$gse185643.36h <- MapQuery(anchorset = anchors, reference = sr.integrate$hs.cs7,
                                       query = sr.integrate$gse185643.36h,
                                       refdata = list(celltype = "CellType"),
                                       reference.reduction = "pca", reduction.model = "umap")
sr.integrate$gse185643.36h$CellType <- sr.integrate$gse185643.36h$predicted.celltype
# annotate cells (48h)
anchors <- FindTransferAnchors(reference = sr.integrate$hs.cs7, query = sr.integrate$gse185643.48h,
                               dims = 1:25, reference.reduction = "pca")
sr.integrate$gse185643.48h <- MapQuery(anchorset = anchors, reference = sr.integrate$hs.cs7,
                                       query = sr.integrate$gse185643.48h,
                                       refdata = list(celltype = "CellType"),
                                       reference.reduction = "pca", reduction.model = "umap")
sr.integrate$gse185643.48h$CellType <- sr.integrate$gse185643.48h$predicted.celltype
# integrating
sr.integrate <- lapply(X = sr.integrate, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, features = rownames(x), verbose = FALSE)
  x <- RunPCA(x, features = VariableFeatures(x), verbose = FALSE)
})
ElbowPlot(sr.integrate$gasd5)
ElbowPlot(sr.integrate$gse185643.24h)
ElbowPlot(sr.integrate$gse185643.36h)
ElbowPlot(sr.integrate$gse185643.48h)
ElbowPlot(sr.integrate$gse144897)
ElbowPlot(sr.integrate$hs.cs7)
sr.integrate$gasd5@project.name <- "gasd5"
sr.integrate$gse185643.24h@project.name <- "gse185643.24h"
sr.integrate$gse185643.36h@project.name <- "gse185643.36h"
sr.integrate$gse185643.48h@project.name <- "gse185643.48h"
sr.integrate$gse144897@project.name <- "gse144897"
sr.integrate$hs.cs7@project.name <- "hs.cs7"
dim.n <- c(7, 10, 7, 7, 6, 15)
names(dim.n) <- names(sr.integrate)
sr.integrate <- lapply(X = sr.integrate, FUN = function(x) {
  dims <- unname(dim.n[x@project.name])
  x <- RunUMAP(x, dims = 1:dims, return.model = T) %>% RunTSNE(dims = 1:dims)
})
anchors <- FindIntegrationAnchors(object.list = sr.integrate, dims = 1:25)
sr.integrated <- IntegrateData(anchorset = anchors, dims = 1:25)
DefaultAssay(sr.integrated) <- "integrated"
sr.integrated <- ScaleData(sr.integrated, verbose = FALSE)
sr.integrated <- RunPCA(sr.integrated, npcs = 30, verbose = FALSE)
ElbowPlot(sr.integrated, ndims = 30)
sr.integrated <- RunUMAP(sr.integrated, reduction = "pca", dims = 1:10)
# plotting clustering
dev.off()
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASd5_and_GSE185643_GSE144897.pdf"), height = 5, width = 12)
p1 <- DimPlot(sr.integrated, reduction = "umap", group.by = "Datasets",
              pt.size = 0.25, label = TRUE, repel = TRUE, cols = mk.col)
p2 <- DimPlot(sr.integrated, reduction = "umap", group.by = "CellType",
              pt.size = 0.25, label = TRUE, repel = TRUE, cols = pd.col.fix)
p1 + p2
dev.off()
for (i in unique(sr.integrated$Datasets)) {
  sr.integrated$CellType.tmp <- sr.integrated$CellType
  sr.integrated$CellType.tmp[sr.integrated$Datasets != i] <- "Others"
  pd.col.fix.tmp <- c(pd.col.fix, c(Others = "#E7E7E7"))
  pdf(file.path(outdir, paste0("Integrated_UMAP_clustering_highlighted_", i, ".pdf")), height = 5, width = 6)
  print(DimPlot(sr.integrated, reduction = "umap", group.by = "CellType.tmp",
                pt.size = 0.25, label = TRUE, repel = TRUE, cols = pd.col.fix.tmp, order = names(pd.col.fix.tmp)))
  dev.off()
}
for (i in unique(sr.integrated$Datasets)[2:4]) {
  sr.integrated$CellType.tmp <- sr.integrated$AuthorLabel
  sr.integrated$CellType.tmp[sr.integrated$Datasets != i] <- "Others"
  pd.col.fix.tmp <- mk.col2[1:(length(unique(sr.integrated$AuthorLabel))-1)]
  names(pd.col.fix.tmp) <- sort(unique(sr.integrated$AuthorLabel))
  pd.col.fix.tmp["Others"] <- "#E7E7E7"
  pdf(file.path(outdir, paste0("Integrated_UMAP_clustering_highlighted_", i, "_raw_label.pdf")), height = 5, width = 6)
  print(DimPlot(sr.integrated, reduction = "umap", group.by = "CellType.tmp",
                pt.size = 0.25, label = TRUE, repel = TRUE, cols = pd.col.fix.tmp, order = names(pd.col.fix.tmp)))
  dev.off()
}
# plot stat
Vis.Annotation.Ratio(sr.meta = sr.integrated@meta.data, annotation = c("Datasets", "CellType"),
                     pd.title = "Integrated", pd.height = 5, pd.width = 10, pd.col = pd.col.fix,
                     res.out = file.path(outdir, "Cell_type_stat_by_datasets"))
Vis.Annotation.Ratio(sr.meta = sr.integrated@meta.data, annotation = c("CellType", "Datasets"),
                     pd.title = "Integrated", pd.height = 5, pd.width = 10,
                     res.out = file.path(outdir, "Cell_type_stat_by_celltype"))
Vis.Annotation.Relationship(sr.meta = subset(sr.integrated@meta.data, Datasets %in% c("Embryoids.24h", "Embryoids.36h", "Embryoids.48h")),
                            annotation = c("Datasets", "AuthorLabel", "CellType"),
                            pd.title = "Relationship between various annotations",
                            pd.col = NULL, split.meta = FALSE, split.by = "AuthorLabel",
                            pd.height = 15, pd.width = 15,
                            res.out = file.path(outdir, "Sankey_plot_GSE185643_sample"))
Vis.Annotation.Relationship(sr.meta = subset(sr.integrated@meta.data, Datasets %in% c("Embryoids.24h", "Embryoids.36h", "Embryoids.48h")),
                            annotation = c("Datasets", "AuthorLabel", "CellType"),
                            pd.title = "Relationship between various annotations",
                            pd.col = NULL, split.meta = TRUE, split.by = "AuthorLabel",
                            pd.height = 15, pd.width = 15,
                            res.out = file.path(outdir, "Sankey_plot_GSE185643_sample"))
# plot expression
sr.out <- file.path(outdir, "Marker_genes_four_datasets")
dir.create(sr.out, recursive = T)
DefaultAssay(sr.integrated) <- "RNA"
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(sr.integrated))
}
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(sr.integrated, features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 0.2, slot = "data", min.cutoff = 0.1))
    dev.off()
  }; rm(g)
}


### >>> 4. Correlation
# Hs CS7 VS GSE185643 (24h)
expr.1 <- AverageExpression(sr.integrate$hs.cs7, assays = "RNA", slot = "data", group.by = "CellType")
expr.2$RNA <- AverageExpression(sr.integrate$gse185643.24h, assays = "RNA", slot = "data", group.by = "CellType")
hvg.1 <- VariableFeatures(sr.integrate$hs.cs7)
hvg.2 <- VariableFeatures(sr.integrate$gse185643.24h)
homo.gene <- data.frame(species.1 = intersect(rownames(sr.integrate$hs.cs7), rownames(sr.integrate$gse185643.24h)),
                        species.2 = intersect(rownames(sr.integrate$hs.cs7), rownames(sr.integrate$gse185643.24h)))
CorrComparePlot(ExpressionTableSpecies1 = expr.1$RNA, DEgenesSpecies1 = hvg.1,
                ExpressionTableSpecies2 = expr.2$RNA, DEgenesSpecies2 = hvg.2,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(outdir, "Correlation_between_Hs_CS7_and_GSE185643.24h_cell_types"),
                pd.order = "original", Height = 15, Width = 15)
# Hs CS7 VS GSE185643 (36h)
expr.1 <- AverageExpression(sr.integrate$hs.cs7, assays = "RNA", slot = "data", group.by = "CellType")
expr.2 <- AverageExpression(sr.integrate$gse185643.36h, assays = "RNA", slot = "data", group.by = "CellType")
hvg.1 <- VariableFeatures(sr.integrate$hs.cs7)
hvg.2 <- VariableFeatures(sr.integrate$gse185643.36h)
homo.gene <- data.frame(species.1 = intersect(rownames(sr.integrate$hs.cs7), rownames(sr.integrate$gse185643.36h)),
                        species.2 = intersect(rownames(sr.integrate$hs.cs7), rownames(sr.integrate$gse185643.36h)))
CorrComparePlot(ExpressionTableSpecies1 = expr.1$RNA, DEgenesSpecies1 = hvg.1,
                ExpressionTableSpecies2 = expr.2$RNA, DEgenesSpecies2 = hvg.2,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(outdir, "Correlation_between_Hs_CS7_and_GSE185643.36h_cell_types"),
                pd.order = "original", Height = 15, Width = 15)
# Hs CS7 VS GSE185643 (48h)
expr.1 <- AverageExpression(sr.integrate$hs.cs7, assays = "RNA", slot = "data", group.by = "CellType")
expr.2 <- AverageExpression(sr.integrate$gse185643.48h, assays = "RNA", slot = "data", group.by = "CellType")
hvg.1 <- VariableFeatures(sr.integrate$hs.cs7)
hvg.2 <- VariableFeatures(sr.integrate$gse185643.48h)
homo.gene <- data.frame(species.1 = intersect(rownames(sr.integrate$hs.cs7), rownames(sr.integrate$gse185643.48h)),
                        species.2 = intersect(rownames(sr.integrate$hs.cs7), rownames(sr.integrate$gse185643.48h)))
CorrComparePlot(ExpressionTableSpecies1 = expr.1$RNA, DEgenesSpecies1 = hvg.1,
                ExpressionTableSpecies2 = expr.2$RNA, DEgenesSpecies2 = hvg.2,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(outdir, "Correlation_between_Hs_CS7_and_GSE185643.48h_cell_types"),
                pd.order = "original", Height = 15, Width = 15)
# Hs CS7 VS GSE144897
expr.1 <- AverageExpression(sr.integrate$hs.cs7, assays = "RNA", slot = "data", group.by = "CellType")
expr.2 <- AverageExpression(sr.integrate$gse144897, assays = "RNA", slot = "data", group.by = "CellType")
hvg.1 <- VariableFeatures(sr.integrate$hs.cs7)
hvg.2 <- VariableFeatures(sr.integrate$gse144897)
homo.gene <- data.frame(species.1 = intersect(rownames(sr.integrate$hs.cs7), rownames(sr.integrate$gse144897)),
                        species.2 = intersect(rownames(sr.integrate$hs.cs7), rownames(sr.integrate$gse144897)))
CorrComparePlot(ExpressionTableSpecies1 = expr.1$RNA, DEgenesSpecies1 = hvg.1,
                ExpressionTableSpecies2 = expr.2$RNA, DEgenesSpecies2 = hvg.2,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(outdir, "Correlation_between_Hs_CS7_and_GSE144897_cell_types"),
                pd.order = "original", Height = 15, Width = 15)
# Hs CS7 vs GASD5 VS
expr.1 <- AverageExpression(sr.integrate$hs.cs7, assays = "RNA", slot = "data", group.by = "CellType")
expr.2 <- AverageExpression(sr.integrate$gasd5, assays = "RNA", slot = "data", group.by = "CellType")
hvg.1 <- VariableFeatures(sr.integrate$hs.cs7)
hvg.2 <- VariableFeatures(sr.integrate$gasd5)
homo.gene <- data.frame(species.1 = intersect(rownames(sr.integrate$gasd5), rownames(sr.integrate$hs.cs7)),
                        species.2 = intersect(rownames(sr.integrate$gasd5), rownames(sr.integrate$hs.cs7)))
CorrComparePlot(ExpressionTableSpecies1 = expr.1$RNA, DEgenesSpecies1 = hvg.1,
                ExpressionTableSpecies2 = expr.2$RNA, DEgenesSpecies2 = hvg.2,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(outdir, "Correlation_between_Hs_CS7_and_GASD5_cell_types"),
                pd.order = "original", Height = 15, Width = 15)



### ==============================
### 14th part: Endoderm + Mesoderm ----
### ==============================

### >>> 1. Setting output dir
outdir <- file.path(getwd(), "graphs/Endo_Meso")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}


### >>> 2. Endoderm analysis
sr.endo <- list()
sr.endo$Mk.CS8 <- subset(sr.pdt.hs.mk$mk.cs8, CellType %in% c("DE", "Gut", "VE", "ys.Endo1", "ys.Endo2"))
sr.endo$Mk.CS8$Datasets <- "Mk.CS8"
sr.endo$Mk.CS8$CellType.mk <- sr.endo$Mk.CS8$CellType
sr.endo$Hs.CS7 <- subset(sr.pdt.hs.mk$hs.to.mk, CellType %in% c("End"))
sr.endo$Hs.CS7$Datasets <- "Hs.CS7"
sr.endo$Hs.CS7$CellType.mk <- sr.endo$Hs.CS7$predicted.celltype
sr.endo$GAS.D5 <- subset(sr.pdt.hs.mk$gas.to.mk, CellType %in% c("End"))
sr.endo$GAS.D5$Datasets <- "GAS.D5"
sr.endo$GAS.D5$CellType.mk <- sr.endo$GAS.D5$predicted.celltype
# integrating
sr.endo <- lapply(X = sr.endo, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, features = rownames(x), verbose = FALSE)
  x <- RunPCA(x, features = VariableFeatures(x), verbose = FALSE)
})
for (i in names(sr.endo)) {
  sr.endo[[i]]@project.name <- i
  print(ElbowPlot(sr.endo[[i]]))
}
dim.n <- c(5, 5, 5)
names(dim.n) <- names(sr.endo)
sr.endo <- lapply(X = sr.endo, FUN = function(x) {
  dims <- unname(dim.n[x@project.name])
  x <- RunUMAP(x, dims = 1:dims, return.model = T) %>% RunTSNE(dims = 1:dims)
})
anchors <- FindIntegrationAnchors(object.list = sr.endo, dims = 1:25)
sr.endo.merge <- IntegrateData(anchorset = anchors, dims = 1:25)
DefaultAssay(sr.endo.merge) <- "integrated"
sr.endo.merge <- ScaleData(sr.endo.merge, verbose = FALSE)
sr.endo.merge <- RunPCA(sr.endo.merge, npcs = 30, verbose = FALSE)
ElbowPlot(sr.endo.merge, ndims = 30)
sr.endo.merge <- RunUMAP(sr.endo.merge, reduction = "pca", dims = 1:10)
# plotting clustering
dev.off()
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_endoderm.pdf"), height = 5, width = 12)
p1 <- DimPlot(sr.endo.merge, reduction = "umap", group.by = "Datasets",
              pt.size = 1, label = TRUE, repel = TRUE, cols = mk.col)
p2 <- DimPlot(sr.endo.merge, reduction = "umap", group.by = "CellType",
              pt.size = 1, label = TRUE, repel = TRUE, cols = mk.col)
p1 + p2
dev.off()
sr.endo.merge$Merged.CellType <- paste0(sr.endo.merge$Datasets, "-", sr.endo.merge$CellType)
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_endoderm_CellType_with_Datasets.pdf"), height = 5, width = 6.5)
DimPlot(sr.endo.merge, reduction = "umap", group.by = "Merged.CellType",
        pt.size = 1, label = TRUE, repel = TRUE, cols = mk.col)
dev.off()
pdf(file.path(outdir, "GASD5_endoderm_UMAP_clustering_mapped_GFP_positive.pdf"), height = 5, width = 6)
DimPlot(subset(sr.endo.merge, Datasets == "GAS.D5"), reduction = "umap", group.by = c("Type"),
        label = TRUE, repel = TRUE, pt.size = 1.5, cols = mk.col)
dev.off()
# plot expression
sr.out <- file.path(outdir, "Marker_genes_of_endoderm")
dir.create(sr.out, recursive = T)
DefaultAssay(sr.endo.merge) <- "RNA"
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(sr.endo.merge))
}
table(sr.endo.merge$Datasets)
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(sr.endo.merge, features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 1, slot = "data", min.cutoff = 0.1))
    dev.off()
    # GAS.D5
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], "_in_GAS.D5.pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(subset(sr.endo.merge, Datasets == "GAS.D5"), features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 1, slot = "data", min.cutoff = 0.1))
    dev.off()
    # Hs.CS7
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], "_in_Hs.CS7.pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(subset(sr.endo.merge, Datasets == "Hs.CS7"), features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 1, slot = "data", min.cutoff = 0.1))
    dev.off()
    # Mk.CS8
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], "_in_Mk.CS8.pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(subset(sr.endo.merge, Datasets == "Mk.CS8"), features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 1, slot = "data", min.cutoff = 0.1))
    dev.off()
  }; rm(g)
}
# LUM gene
mk.gene.expr <- read.csv("/home/yhw/bioinfo/project-wb/hs_gastruloids/analysis/tables/SOX17.csv")
rownames(mk.gene.expr) <- mk.gene.expr$X
mk.gene.expr <- mk.gene.expr[,-1]
mk.gene.expr <- t(mk.gene.expr) %>% as.data.frame()
rownames(mk.gene.expr) <- gsub("\\.", "-", str_split_fixed(rownames(mk.gene.expr), "_", 2)[, 2])
mk.gene.expr <- t(mk.gene.expr)
pd <- subset(sr.endo.merge, Datasets %in% c("GAS.D5", "Hs.CS7", "Mk.CS8"))
tmp <- rbind(GetAssayData(hs.cs7, slot = "data")["LUM", intersect(colnames(pd), colnames(hs.cs7))] %>% as.data.frame(),
             GetAssayData(sr.pdt$gasd5.to.hs, slot = "data")["LUM", intersect(colnames(pd), colnames(sr.pdt$gasd5.to.hs))] %>% as.data.frame(),
             mk.gene.expr["LUM",intersect(colnames(pd), colnames(mk.gene.expr))] %>% as.data.frame())
if (all(colnames(pd) %in% rownames(tmp))) {
  tmp <- tmp[colnames(pd),]
}
pd$LUM <- tmp
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_endoderm_gene_LUM.pdf"),
    height = 4.5, width = 5)
FeaturePlot(pd, features = "LUM", order = T,
            cols = c("#ebebeb", "#7d3c98"), pt.size = 1, min.cutoff = 0.1)
dev.off()
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_endoderm_gene_LUM_in_GAS.D5.pdf"),
    height = 4.5, width = 5)
FeaturePlot(subset(pd, Datasets == "GAS.D5"), features = "LUM", order = T,
            cols = c("#ebebeb", "#7d3c98"), pt.size = 2, min.cutoff = 0.1)
dev.off()
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_endoderm_gene_LUM_in_Hs.CS7.pdf"),
    height = 4.5, width = 5)
FeaturePlot(subset(pd, Datasets == "Hs.CS7"), features = "LUM", order = T,
            cols = c("#ebebeb", "#7d3c98"), pt.size = 2, min.cutoff = 0.1)
dev.off()
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_endoderm_gene_LUM_in_Mk.CS8.pdf"),
    height = 4.5, width = 5)
FeaturePlot(subset(pd, Datasets == "Mk.CS8"), features = "LUM", order = T,
            cols = c("#ebebeb", "#7d3c98"), pt.size = 2, min.cutoff = 0.1)
dev.off()
# SOX17 gene
pd <- subset(sr.endo.merge, Datasets %in% c("GAS.D5", "Hs.CS7", "Mk.CS8"))
tmp <- rbind(GetAssayData(hs.cs7, slot = "data")["SOX17", intersect(colnames(pd), colnames(hs.cs7))] %>% as.data.frame(),
             GetAssayData(sr.pdt$gasd5.to.hs, slot = "data")["SOX17", intersect(colnames(pd), colnames(sr.pdt$gasd5.to.hs))] %>% as.data.frame(),
             mk.gene.expr["SOX17", intersect(colnames(pd), colnames(mk.gene.expr))] %>% as.data.frame())
if (all(colnames(pd) %in% rownames(tmp))) {
  tmp <- tmp[colnames(pd),]
}
pd$SOX17 <- tmp
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_endoderm_gene_SOX17.pdf"),
    height = 4.5, width = 5)
FeaturePlot(pd, features = "SOX17", order = T,
            cols = c("#ebebeb", "#7d3c98"), pt.size = 1, min.cutoff = 0.1)
dev.off()
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_endoderm_gene_SOX17_in_GAS.D5.pdf"),
    height = 4.5, width = 5)
FeaturePlot(subset(pd, Datasets == "GAS.D5"), features = "SOX17", order = T,
            cols = c("#ebebeb", "#7d3c98"), pt.size = 2, min.cutoff = 0.1)
dev.off()
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_endoderm_gene_SOX17_in_Hs.CS7.pdf"),
    height = 4.5, width = 5)
FeaturePlot(subset(pd, Datasets == "Hs.CS7"), features = "SOX17", order = T,
            cols = c("#ebebeb", "#7d3c98"), pt.size = 2, min.cutoff = 0.1)
dev.off()
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_endoderm_gene_SOX17_in_Mk.CS8.pdf"),
    height = 4.5, width = 5)
FeaturePlot(subset(pd, Datasets == "Mk.CS8"), features = "SOX17", order = T,
            cols = c("#ebebeb", "#7d3c98"), pt.size = 2, min.cutoff = 0.1)
dev.off()
# highlight cells (monkey annotation)
for (i in unique(sr.endo.merge$Datasets)) {
  sr.endo.merge$CellType.tmp <- sr.endo.merge$CellType.mk
  sr.endo.merge$CellType.tmp[sr.endo.merge$Datasets != i] <- "Others"
  pd.col.fix.tmp <- mk.col2[1:(length(unique(sr.endo.merge$CellType.tmp)))]
  names(pd.col.fix.tmp) <- sort(unique(sr.endo.merge$CellType.tmp))
  pd.col.fix.tmp["Others"] <- "#E7E7E7"
  pdf(file.path(outdir, paste0("Integrated_UMAP_endoderm_clustering_with_monkey_annotation_highlighted_", i, ".pdf")),
      height = 5, width = 6)
  print(DimPlot(sr.endo.merge, reduction = "umap", group.by = "CellType.tmp",
                pt.size = 1, label = TRUE, repel = TRUE, cols = pd.col.fix.tmp, order = names(pd.col.fix.tmp)))
  dev.off()
}
# highlight cells (human annotation)
for (i in unique(sr.endo.merge$Datasets)) {
  sr.endo.merge$CellType.tmp <- sr.endo.merge$CellType
  sr.endo.merge$CellType.tmp[sr.endo.merge$Datasets != i] <- "Others"
  pd.col.fix.tmp <- c(pd.col.fix, c(Others = "#E7E7E7"))
  pdf(file.path(outdir, paste0("Integrated_UMAP_endoderm_clustering_with_human_annotation_highlighted_", i, ".pdf")),
      height = 5, width = 6)
  print(DimPlot(sr.endo.merge, reduction = "umap", group.by = "CellType.tmp",
                pt.size = 1, label = TRUE, repel = TRUE, cols = pd.col.fix.tmp, order = names(pd.col.fix.tmp)))
  dev.off()
}
# Correlation: GAS.D5 vs Hs.CS7
expr.1 <- AverageExpression(sr.endo$GAS.D5, assays = "RNA", slot = "data", group.by = "CellType.mk")
expr.2 <- AverageExpression(sr.endo$Hs.CS7, assays = "RNA", slot = "data", group.by = "CellType.mk")
hvg.1 <- VariableFeatures(sr.endo$GAS.D5)
hvg.2 <- VariableFeatures(sr.endo$Hs.CS7)
homo.gene <- data.frame(species.1 = intersect(rownames(sr.endo$GAS.D5), rownames(sr.endo$Hs.CS7)),
                        species.2 = intersect(rownames(sr.endo$GAS.D5), rownames(sr.endo$Hs.CS7)))
CorrComparePlot(ExpressionTableSpecies1 = expr.1$RNA, DEgenesSpecies1 = hvg.1,
                ExpressionTableSpecies2 = expr.2$RNA, DEgenesSpecies2 = hvg.2,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(outdir, "Correlation_of_endoderm_between_GAS_D5_and_Hs_CS7"),
                pd.order = "original", Height = 15, Width = 15)
# Correlation: GAS.D5 vs Mk.CS8
expr.1 <- AverageExpression(sr.endo$GAS.D5, assays = "RNA", slot = "data", group.by = "CellType.mk")
expr.2 <- AverageExpression(sr.endo$Mk.CS8, assays = "RNA", slot = "data", group.by = "CellType.mk")
hvg.1 <- VariableFeatures(sr.endo$GAS.D5)
hvg.2 <- VariableFeatures(sr.endo$Mk.CS8)
homo.gene <- data.frame(species.1 = intersect(rownames(sr.endo$GAS.D5), rownames(sr.endo$Mk.CS8)),
                        species.2 = intersect(rownames(sr.endo$GAS.D5), rownames(sr.endo$Mk.CS8)))
CorrComparePlot(ExpressionTableSpecies1 = expr.1$RNA, DEgenesSpecies1 = hvg.1,
                ExpressionTableSpecies2 = expr.2$RNA, DEgenesSpecies2 = hvg.2,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(outdir, "Correlation_of_endoderm_between_GAS_D5_and_Mk_CS8"),
                pd.order = "original", Height = 15, Width = 15)
pd <- cbind(expr.1$RNA[intersect(rownames(expr.1$RNA), rownames(expr.2$RNA)),] %>% as.data.frame(),
            expr.2$RNA[intersect(rownames(expr.1$RNA), rownames(expr.2$RNA)),] %>% as.data.frame())
colnames(pd) <- c(paste0("GASD5.", colnames(expr.1$RNA)),
                  paste0("Monkey.", colnames(expr.2$RNA)))
pd <- pd[c(intersect(hvg.1, hvg.2), c("FOXA2", "CDX2", "ISL1", "HOXB9", "HOXA10", "HOXA13", "CDH2", "SHH")), ]
pd %>%
  mutate(type = c(rep("Others", nrow(pd)-8), rep("Gut", 8))) %>%
  select("GASD5.Gut", "Monkey.Gut", "type") %>%
  rownames_to_column(var = "SYMBOL") %>%
  filter(SYMBOL %in% c("FOXA2", "CDX2", "ISL1", "HOXB9", "HOXA10", "HOXA13", "CDH2", "SHH")) -> pd.label
pdf(file.path(outdir, "Correlation_analysis_of_gut_between_Hs_CS7_GASD5_and_Mk_CS8_endoderm.pdf"),
    height = 4.5, width = 5.5)
pd %>%
  mutate(type = c(rep("Others", nrow(pd)-8), rep("Gut", 8))) %>%
  select("GASD5.Gut", "Monkey.Gut", "type") %>%
  ggplot(aes(x = log2(GASD5.Gut + 0.1), y = log2(Monkey.Gut + 0.1))) +
  #ggplot(aes(x = log1p(GASD5.Gut), y = log1p(Monkey.Gut))) +
  geom_point(aes(color = type, size = type), shape = 16) +
  scale_color_manual(values = c("Others" = "#000000", "Gut" = "#E91E78")) +
  scale_size_manual(values = c("Others" = 1.5, "Gut" = 3)) +
  geom_smooth(method = "lm", se = F) +
  stat_cor(method = "pearson") +
  labs(x = "GAS.D5", y = "Monkey embryo") +
  geom_label_repel(aes(label = SYMBOL), data = pd.label,
                   box.padding = 0.2, point.padding = 0.2, segment.color = "#130f40",
                   label.size = 0.2, fill = NA, size = 3) +
  theme_bw() +
  theme(
    axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
    axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
    axis.text.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
    axis.text.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
    panel.grid = element_blank()
  )
dev.off()
# Correlation: Hs.CS7 vs Mk.CS8
expr.1 <- AverageExpression(sr.endo$Hs.CS7, assays = "RNA", slot = "data", group.by = "CellType.mk")
expr.2 <- AverageExpression(sr.endo$Mk.CS8, assays = "RNA", slot = "data", group.by = "CellType.mk")
hvg.1 <- VariableFeatures(sr.endo$Hs.CS7)
hvg.2 <- VariableFeatures(sr.endo$Mk.CS8)
homo.gene <- data.frame(species.1 = intersect(rownames(sr.endo$Hs.CS7), rownames(sr.endo$Mk.CS8)),
                        species.2 = intersect(rownames(sr.endo$Hs.CS7), rownames(sr.endo$Mk.CS8)))
CorrComparePlot(ExpressionTableSpecies1 = expr.1$RNA, DEgenesSpecies1 = hvg.1,
                ExpressionTableSpecies2 = expr.2$RNA, DEgenesSpecies2 = hvg.2,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(outdir, "Correlation_of_endoderm_between_Hs_CS7_and_Mk_CS8"),
                pd.order = "original", Height = 15, Width = 15)
# differential expression analysis
table(sr.endo.merge$Datasets, sr.endo.merge$CellType.mk)
deg.endo <- list()
deg.cells <- c("DE", "VE")
for (i in deg.cells) {
  # GAS.D5 vs Mk.CS8
  tmp.sr <- subset(sr.endo.merge, CellType.mk == i)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Datasets", comparison = c("Mk.CS8", "GAS.D5"), sample.n = NULL)
  deg.endo[[paste0(i, "_GAS.D5")]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                  sample.n = NULL, group.by = "CellType",
                                                  g1 = "1_Mk.CS8", g2 = "2_GAS.D5", lfc = 1, sig = 0.05,
                                                  res.out = file.path(outdir, paste0("DEGs/Endo/", i, "_GAS.D5_vs_Mk.CS8")))
  # Hs.CS7 vs Mk.CS8
  tmp.sr <- subset(sr.endo.merge, CellType.mk == i)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Datasets", comparison = c("Mk.CS8", "Hs.CS7"), sample.n = NULL)
  deg.endo[[paste0(i, "_Hs.CS7")]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                  sample.n = NULL, group.by = "CellType",
                                                  g1 = "1_Mk.CS8", g2 = "2_Hs.CS7", lfc = 1, sig = 0.05,
                                                  res.out = file.path(outdir, paste0("DEGs/Endo/", i, "_Hs.CS7_vs_Mk.CS8")))
}
deg.cells <- c("Gut")
for (i in deg.cells) {
  # GAS.D5 vs Mk.CS8
  tmp.sr <- subset(sr.endo.merge, CellType.mk == i)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Datasets", comparison = c("Mk.CS8", "GAS.D5"), sample.n = NULL)
  deg.endo[[paste0(i, "_GAS.D5")]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                  sample.n = NULL, group.by = "CellType",
                                                  g1 = "1_Mk.CS8", g2 = "2_GAS.D5", lfc = 1, sig = 0.05,
                                                  res.out = file.path(outdir, paste0("DEGs/Endo/", i, "_GAS.D5_vs_Mk.CS8")))
}
deg.cells <- c("ys.Endo1", "ys.Endo2")
for (i in deg.cells) {
  # Hs.CS7 vs Mk.CS8
  tmp.sr <- subset(sr.endo.merge, CellType.mk == i)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Datasets", comparison = c("Mk.CS8", "Hs.CS7"), sample.n = NULL)
  deg.endo[[paste0(i, "_Hs.CS7")]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                  sample.n = NULL, group.by = "CellType",
                                                  g1 = "1_Mk.CS8", g2 = "2_Hs.CS7", lfc = 1, sig = 0.05,
                                                  res.out = file.path(outdir, paste0("DEGs/Endo/", i, "_Hs.CS7_vs_Mk.CS8")))
}
# volcano plot
for (i in names(deg.endo)) {
  pdf(file.path(outdir, paste0("DEGs/Endo/", i, "_vs_Mk.CS8_volcano.pdf")), height = 6, width = 9)
  print(VisDEG.volcano(deg.data = deg.endo[[i]]$all, geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.05, lfc = 1, title = i,
                       pt.size = 2, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
}


### >>> 3. Mesoderm analysis
sr.meso <- list()
sr.meso$Mk.CS8 <- subset(sr.pdt.hs.mk$mk.cs8, CellType %in% c("Al", "Caud.Meso", "ExE.Meso", "LP.Meso",
                                                              "Mes", "Nas.Meso", "Pharyg.Meso",
                                                              "ys.Meso1", "ys.Meso2"))
sr.meso$Mk.CS8$Datasets <- "Mk.CS8"
sr.meso$Mk.CS8$CellType.mk <- sr.meso$Mk.CS8$CellType
sr.meso$Hs.CS7 <- subset(sr.pdt.hs.mk$hs.to.mk, CellType %in% c("NM", "EM", "AM", "ExM"))
sr.meso$Hs.CS7$Datasets <- "Hs.CS7"
sr.meso$Hs.CS7$CellType.mk <- sr.meso$Hs.CS7$predicted.celltype
sr.meso$GAS.D5 <- subset(sr.pdt.hs.mk$gas.to.mk, CellType %in% c("NM", "EM", "AM", "ExM"))
sr.meso$GAS.D5$Datasets <- "GAS.D5"
sr.meso$GAS.D5$CellType.mk <- sr.meso$GAS.D5$predicted.celltype
# integrating
sr.meso <- lapply(X = sr.meso, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, features = rownames(x), verbose = FALSE)
  x <- RunPCA(x, features = VariableFeatures(x), verbose = FALSE)
})
for (i in names(sr.meso)) {
  sr.meso[[i]]@project.name <- i
  print(ElbowPlot(sr.meso[[i]]))
}
dim.n <- c(10, 10, 10)
names(dim.n) <- names(sr.meso)
sr.meso <- lapply(X = sr.meso, FUN = function(x) {
  dims <- unname(dim.n[x@project.name])
  x <- RunUMAP(x, dims = 1:dims, return.model = T) %>% RunTSNE(dims = 1:dims)
})
anchors <- FindIntegrationAnchors(object.list = sr.meso, dims = 1:25)
sr.meso.merge <- IntegrateData(anchorset = anchors, dims = 1:25)
DefaultAssay(sr.meso.merge) <- "integrated"
sr.meso.merge <- ScaleData(sr.meso.merge, verbose = FALSE)
sr.meso.merge <- RunPCA(sr.meso.merge, npcs = 30, verbose = FALSE)
ElbowPlot(sr.meso.merge, ndims = 30)
sr.meso.merge <- RunUMAP(sr.meso.merge, reduction = "pca", dims = 1:10)
# plotting clustering
dev.off()
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_mesoderm.pdf"), height = 5, width = 12)
p1 <- DimPlot(sr.meso.merge, reduction = "umap", group.by = "Datasets",
              pt.size = 1, label = TRUE, repel = TRUE, cols = mk.col)
p2 <- DimPlot(sr.meso.merge, reduction = "umap", group.by = "CellType",
              pt.size = 1, label = TRUE, repel = TRUE, cols = mk.col)
p1 + p2
dev.off()
sr.meso.merge$Merged.CellType <- paste0(sr.meso.merge$Datasets, "-", sr.meso.merge$CellType)
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_mesoderm_CellType_with_Datasets.pdf"),
    height = 5, width = 7)
DimPlot(sr.meso.merge, reduction = "umap", group.by = "Merged.CellType",
        pt.size = 1, label = TRUE, repel = TRUE, cols = mk.col)
dev.off()
pdf(file.path(outdir, "GASD5_mesoderm_UMAP_clustering_mapped_GFP_positive.pdf"), height = 5, width = 6)
DimPlot(subset(sr.meso.merge, Datasets == "GAS.D5"), reduction = "umap", group.by = c("Type"),
        label = TRUE, repel = TRUE, pt.size = 1.5, cols = mk.col)
dev.off()
# plot expression
sr.out <- file.path(outdir, "Marker_genes_of_mesoderm")
dir.create(sr.out, recursive = T)
DefaultAssay(sr.meso.merge) <- "RNA"
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(sr.meso.merge))
}
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(sr.meso.merge, features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 1, slot = "data", min.cutoff = 0.1))
    dev.off()
    # GAS.D5
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], "_in_GAS.D5.pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(subset(sr.meso.merge, Datasets == "GAS.D5"), features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 1, slot = "data", min.cutoff = 0.1))
    dev.off()
    # Hs.CS7
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], "_in_Hs.CS7.pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(subset(sr.meso.merge, Datasets == "Hs.CS7"), features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 1, slot = "data", min.cutoff = 0.1))
    dev.off()
    # Mk.CS8
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], "_in_Mk.CS8.pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(subset(sr.meso.merge, Datasets == "Mk.CS8"), features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 1, slot = "data", min.cutoff = 0.1))
    dev.off()
  }; rm(g)
}
# LUM gene
pd <- subset(sr.meso.merge, Datasets %in% c("GAS.D5", "Hs.CS7", "Mk.CS8"))
tmp <- rbind(GetAssayData(hs.cs7, slot = "data")["LUM", intersect(colnames(pd), colnames(hs.cs7))] %>% as.data.frame(),
             GetAssayData(sr.pdt$gasd5.to.hs, slot = "data")["LUM", intersect(colnames(pd), colnames(sr.pdt$gasd5.to.hs))] %>% as.data.frame(),
             mk.gene.expr["LUM", intersect(colnames(pd), colnames(mk.gene.expr))] %>% as.data.frame())
if (all(colnames(pd) %in% rownames(tmp))) {
  tmp <- tmp[colnames(pd),]
}
pd$LUM <- tmp
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_mesoderm_gene_LUM.pdf"),
    height = 4.5, width = 5)
FeaturePlot(pd, features = "LUM", order = T,
            cols = c("#ebebeb", "#7d3c98"), pt.size = 1, min.cutoff = 0.1)
dev.off()
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_mesoderm_gene_LUM_in_GAS.D5.pdf"),
    height = 4.5, width = 5)
FeaturePlot(subset(pd, Datasets == "GAS.D5"), features = "LUM", order = T,
            cols = c("#ebebeb", "#7d3c98"), pt.size = 2, min.cutoff = 0.1)
dev.off()
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_mesoderm_gene_LUM_in_Hs.CS7.pdf"),
    height = 4.5, width = 5)
FeaturePlot(subset(pd, Datasets == "Hs.CS7"), features = "LUM", order = T,
            cols = c("#ebebeb", "#7d3c98"), pt.size = 2, min.cutoff = 0.1)
dev.off()
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_mesoderm_gene_LUM_in_Mk.CS8.pdf"),
    height = 4.5, width = 5)
FeaturePlot(subset(pd, Datasets == "Mk.CS8"), features = "LUM", order = T,
            cols = c("#ebebeb", "#7d3c98"), pt.size = 2, min.cutoff = 0.1)
dev.off()
# SOX17 gene
pd <- subset(sr.meso.merge, Datasets %in% c("GAS.D5", "Hs.CS7", "Mk.CS8"))
tmp <- rbind(GetAssayData(hs.cs7, slot = "data")["SOX17", intersect(colnames(pd), colnames(hs.cs7))] %>% as.data.frame(),
             GetAssayData(sr.pdt$gasd5.to.hs, slot = "data")["SOX17", intersect(colnames(pd), colnames(sr.pdt$gasd5.to.hs))] %>% as.data.frame(),
             mk.gene.expr["SOX17", intersect(colnames(pd), colnames(mk.gene.expr))] %>% as.data.frame())
if (all(colnames(pd) %in% rownames(tmp))) {
  tmp <- tmp[colnames(pd),]
}
pd$SOX17 <- tmp
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_mesoderm_gene_SOX17.pdf"),
    height = 4.5, width = 5)
FeaturePlot(pd, features = "SOX17", order = T,
            cols = c("#ebebeb", "#7d3c98"), pt.size = 1, min.cutoff = 0.1)
dev.off()
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_mesoderm_gene_SOX17_in_GAS.D5.pdf"),
    height = 4.5, width = 5)
FeaturePlot(subset(pd, Datasets == "GAS.D5"), features = "SOX17", order = T,
            cols = c("#ebebeb", "#7d3c98"), pt.size = 2, min.cutoff = 0.1)
dev.off()
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_mesoderm_gene_SOX17_in_Hs.CS7.pdf"),
    height = 4.5, width = 5)
FeaturePlot(subset(pd, Datasets == "Hs.CS7"), features = "SOX17", order = T,
            cols = c("#ebebeb", "#7d3c98"), pt.size = 2, min.cutoff = 0.1)
dev.off()
pdf(file.path(outdir, "Integrated_UMAP_clustering_with_Hs_CS7_GASD5_and_Mk_CS8_mesoderm_gene_SOX17_in_Mk.CS8.pdf"),
    height = 4.5, width = 5)
FeaturePlot(subset(pd, Datasets == "Mk.CS8"), features = "SOX17", order = T,
            cols = c("#ebebeb", "#7d3c98"), pt.size = 2, min.cutoff = 0.1)
dev.off()
# highlight cells (monkey annotation)
for (i in unique(sr.meso.merge$Datasets)) {
  sr.meso.merge$CellType.tmp <- sr.meso.merge$CellType.mk
  sr.meso.merge$CellType.tmp[sr.meso.merge$Datasets != i] <- "Others"
  pd.col.fix.tmp <- mk.col2[1:(length(unique(sr.meso.merge$CellType.tmp)))]
  names(pd.col.fix.tmp) <- sort(unique(sr.meso.merge$CellType.tmp))
  pd.col.fix.tmp["Others"] <- "#E7E7E7"
  pdf(file.path(outdir, paste0("Integrated_UMAP_mesoderm_clustering_with_monkey_annotation_highlighted_", i, ".pdf")),
      height = 5, width = 6)
  print(DimPlot(sr.meso.merge, reduction = "umap", group.by = "CellType.tmp",
                pt.size = 1, label = TRUE, repel = TRUE, cols = pd.col.fix.tmp, order = names(pd.col.fix.tmp)))
  dev.off()
}
# highlight cells (human annotation)
for (i in unique(sr.meso.merge$Datasets)) {
  sr.meso.merge$CellType.tmp <- sr.meso.merge$CellType
  sr.meso.merge$CellType.tmp[sr.meso.merge$Datasets != i] <- "Others"
  pd.col.fix.tmp <- c(pd.col.fix, c(Others = "#E7E7E7"))
  pdf(file.path(outdir, paste0("Integrated_UMAP_mesoderm_clustering_with_human_annotation_highlighted_", i, ".pdf")),
      height = 5, width = 6)
  print(DimPlot(sr.meso.merge, reduction = "umap", group.by = "CellType.tmp",
                pt.size = 1, label = TRUE, repel = TRUE, cols = pd.col.fix.tmp, order = names(pd.col.fix.tmp)))
  dev.off()
}
# Correlation: GAS.D5 vs Hs.CS7
expr.1 <- AverageExpression(sr.meso$GAS.D5, assays = "RNA", slot = "data", group.by = "CellType.mk")
expr.2 <- AverageExpression(sr.meso$Hs.CS7, assays = "RNA", slot = "data", group.by = "CellType.mk")
hvg.1 <- VariableFeatures(sr.meso$GAS.D5)
hvg.2 <- VariableFeatures(sr.meso$Hs.CS7)
homo.gene <- data.frame(species.1 = intersect(rownames(sr.meso$GAS.D5), rownames(sr.meso$Hs.CS7)),
                        species.2 = intersect(rownames(sr.meso$GAS.D5), rownames(sr.meso$Hs.CS7)))
CorrComparePlot(ExpressionTableSpecies1 = expr.1$RNA, DEgenesSpecies1 = hvg.1,
                ExpressionTableSpecies2 = expr.2$RNA, DEgenesSpecies2 = hvg.2,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(outdir, "Correlation_of_mesoderm_between_GAS_D5_and_Hs_CS7"),
                pd.order = "original", Height = 15, Width = 15)
# Correlation: GAS.D5 vs Mk.CS8
expr.1 <- AverageExpression(sr.meso$GAS.D5, assays = "RNA", slot = "data", group.by = "CellType.mk")
expr.2 <- AverageExpression(sr.meso$Mk.CS8, assays = "RNA", slot = "data", group.by = "CellType.mk")
hvg.1 <- VariableFeatures(sr.meso$GAS.D5)
hvg.2 <- VariableFeatures(sr.meso$Mk.CS8)
homo.gene <- data.frame(species.1 = intersect(rownames(sr.meso$GAS.D5), rownames(sr.meso$Mk.CS8)),
                        species.2 = intersect(rownames(sr.meso$GAS.D5), rownames(sr.meso$Mk.CS8)))
CorrComparePlot(ExpressionTableSpecies1 = expr.1$RNA, DEgenesSpecies1 = hvg.1,
                ExpressionTableSpecies2 = expr.2$RNA, DEgenesSpecies2 = hvg.2,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(outdir, "Correlation_of_mesoderm_between_GAS_D5_and_Mk_CS8"),
                pd.order = "original", Height = 15, Width = 15)
# Correlation: Hs.CS7 vs Mk.CS8
expr.1 <- AverageExpression(sr.meso$Hs.CS7, assays = "RNA", slot = "data", group.by = "CellType.mk")
expr.2 <- AverageExpression(sr.meso$Mk.CS8, assays = "RNA", slot = "data", group.by = "CellType.mk")
hvg.1 <- VariableFeatures(sr.meso$Hs.CS7)
hvg.2 <- VariableFeatures(sr.meso$Mk.CS8)
homo.gene <- data.frame(species.1 = intersect(rownames(sr.meso$Hs.CS7), rownames(sr.meso$Mk.CS8)),
                        species.2 = intersect(rownames(sr.meso$Hs.CS7), rownames(sr.meso$Mk.CS8)))
CorrComparePlot(ExpressionTableSpecies1 = expr.1$RNA, DEgenesSpecies1 = hvg.1,
                ExpressionTableSpecies2 = expr.2$RNA, DEgenesSpecies2 = hvg.2,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(outdir, "Correlation_of_mesoderm_between_Hs_CS7_and_Mk_CS8"),
                pd.order = "original", Height = 15, Width = 15)
# differential expression analysis
table(sr.meso.merge$Datasets, sr.meso.merge$CellType.mk)
deg.meso <- list()
deg.cells <- c("Al", "Caud.Meso", "ExE.Meso", "LP.Meso", "Nas.Meso", "Pharyg.Meso", "ys.Meso2")
for (i in deg.cells) {
  # GAS.D5 vs Mk.CS8
  tmp.sr <- subset(sr.meso.merge, CellType.mk == i)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Datasets", comparison = c("Mk.CS8", "GAS.D5"), sample.n = NULL)
  deg.meso[[paste0(i, "_GAS.D5")]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                  sample.n = NULL, group.by = "CellType",
                                                  g1 = "1_Mk.CS8", g2 = "2_GAS.D5", lfc = 1, sig = 0.05,
                                                  res.out = file.path(outdir, paste0("DEGs/Meso/", i, "_GAS.D5_vs_Mk.CS8")))
  # Hs.CS7 vs Mk.CS8
  tmp.sr <- subset(sr.meso.merge, CellType.mk == i)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Datasets", comparison = c("Mk.CS8", "Hs.CS7"), sample.n = NULL)
  deg.meso[[paste0(i, "_Hs.CS7")]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                  sample.n = NULL, group.by = "CellType",
                                                  g1 = "1_Mk.CS8", g2 = "2_Hs.CS7", lfc = 1, sig = 0.05,
                                                  res.out = file.path(outdir, paste0("DEGs/Meso/", i, "_Hs.CS7_vs_Mk.CS8")))
}
deg.cells <- c("Mes", "ys.Meso1")
for (i in deg.cells) {
  # Hs.CS7 vs Mk.CS8
  tmp.sr <- subset(sr.meso.merge, CellType.mk == i)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Datasets", comparison = c("Mk.CS8", "Hs.CS7"), sample.n = NULL)
  deg.meso[[paste0(i, "_Hs.CS7")]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                  sample.n = NULL, group.by = "CellType",
                                                  g1 = "1_Mk.CS8", g2 = "2_Hs.CS7", lfc = 1, sig = 0.05,
                                                  res.out = file.path(outdir, paste0("DEGs/Meso/", i, "_Hs.CS7_vs_Mk.CS8")))
}
# volcano plot
for (i in names(deg.meso)) {
  pdf(file.path(outdir, paste0("DEGs/Meso/", i, "_vs_Mk.CS8_volcano.pdf")), height = 6, width = 9)
  print(VisDEG.volcano(deg.data = deg.meso[[i]]$all, geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.05, lfc = 1, title = i,
                       pt.size = 2, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
}



### =======================
### 15th part: NMP analysis ----
### =======================

### >>> 1. Setting output dir
outdir <- file.path(getwd(), "graphs/NMP")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}
nmp <- list()


### >>> 2. Monkey CS8
keep.cell <- GetAssayData(mk.cs8, slot = "data")["SOX2",] > 0.1 | GetAssayData(mk.cs8, slot = "data")["TBXT",] > 0.1
nmp$mk <- mk.cs8[, keep.cell]
keep.cell <- names(sort(table(nmp$mk$cell_type))[sort(table(nmp$mk$cell_type)) >= 5])
nmp$mk <- subset(nmp$mk, cell_type %in% keep.cell)
# normalize data and find highly variables genes
nmp$mk <- NormalizeData(nmp$mk) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
# scale data
nmp$mk <- ScaleData(nmp$mk, verbose = FALSE, features = rownames(nmp$mk))
# run PCA
nmp$mk <- RunPCA(nmp$mk, npcs = 50, verbose = FALSE, features = VariableFeatures(nmp$mk))
ElbowPlot(nmp$mk, ndims = 30)
dim.n <- 10
# cluster cells to check batch effect
nmp$mk <- RunUMAP(nmp$mk, reduction = "pca", dims = 1:dim.n, return.model = T) %>%
  RunTSNE(reduction = "pca", dims = 1:dim.n)
DimPlot(nmp$mk, reduction = "umap", group.by = c("sample", "cell_type"),
        label = TRUE, repel = TRUE, pt.size = 1, cols = mk.col)
# remove batch effect
nmp$mk <- RunHarmony(nmp$mk, group.by.vars = c("sample"))
nmp$mk <- RunUMAP(nmp$mk, reduction = "harmony", dims = 1:dim.n) %>%
  RunTSNE(reduction = "harmony", dims = 1:dim.n)
# plot clustering results
pdf(file.path(outdir, "Monkey_SOX2_TBXT_cell_clustering.pdf"), height = 5, width = 12.5)
DimPlot(nmp$mk, reduction = "umap", group.by = c("sample", "cell_type"),
        label = TRUE, repel = TRUE, pt.size = 1.5, cols = mk.col)
dev.off()
pdf(file.path(outdir, "Monkey_SOX2_TBXT_cell_clustering_mapped_expression.pdf"), height = 5, width = 20)
FeaturePlot(nmp$mk, features = c("SOX2", "TBXT"), blend = T, max.cutoff = 0.5, pt.size = 2, cols = c("#E7E7E7", "#7b2c9c", "#f07f21"))
dev.off()
pdf(file.path(outdir, "Monkey_SOX2_TBXT_cell_clustering_mapped_NMP_positive.pdf"), height = 5, width = 12.5)
DimPlot(nmp$mk, reduction = "umap", group.by = c("sample", "CellType.NMP"),
        label = TRUE, repel = TRUE, pt.size = 1.5, cols = mk.col)
dev.off()
library(scCustomize)
pdf(file.path(outdir, "Monkey_SOX2_TBXT_cell_clustering_mapped_density.pdf"), height = 5, width = 7)
Plot_Density_Joint_Only(seurat_object = nmp$mk, features = c("SOX2", "TBXT"), pt.size = 2)
dev.off()
# plot NMP cell composition
nmp$mk$SOX2 <- GetAssayData(nmp$mk, slot = "data")["SOX2", ]
nmp$mk$TBXT <- GetAssayData(nmp$mk, slot = "data")["TBXT", ]
nmp$mk$CellType.NMP <- "NMP.neg"
nmp$mk$CellType.NMP[nmp$mk$SOX2 > 0 & nmp$mk$TBXT > 0] <- "NMP.pos"
Vis.Annotation.Ratio(sr.meta = nmp$mk@meta.data, annotation = c("CellType.NMP", "cell_type"),
                     pd.title = "Stat by Cell Type", pd.col = mk.col,
                     pd.height = 5, pd.width = 20,
                     res.out = file.path(outdir, "Monkey_composition"))
# average expression of positive and negative cell
nmp.expr.mkcs8 <- AverageExpression(nmp$mk, group.by = "CellType.NMP")$RNA %>% as.data.frame()
# plot expression
sr.out <- file.path(outdir, "Marker_genes_of_monkey")
dir.create(sr.out, recursive = T)
DefaultAssay(nmp$mk) <- "RNA"
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(nmp$mk))
}
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(nmp$mk, features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 1, slot = "data", min.cutoff = 0.1))
    dev.off()
  }; rm(g)
}


### >>> 3. GASD5
keep.cell <- GetAssayData(sr.pdt.mk$gasd5.to.mk, slot = "data")["SOX2",] > 0.1 | GetAssayData(sr.pdt.mk$gasd5.to.mk, slot = "data")["TBXT",] > 0.1
table(keep.cell)
nmp$gas <- sr.pdt.mk$gasd5.to.mk[, keep.cell]
keep.cell <- names(sort(table(nmp$gas$CellType))[sort(table(nmp$gas$CellType)) >= 5])
nmp$gas <- subset(nmp$gas, cell_type %in% keep.cell)
# normalize data and find highly variables genes
nmp$gas <- NormalizeData(nmp$gas) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
# scale data
nmp$gas <- ScaleData(nmp$gas, verbose = FALSE, features = rownames(nmp$gas))
# run PCA
nmp$gas <- RunPCA(nmp$gas, npcs = 50, verbose = FALSE, features = VariableFeatures(nmp$gas))
ElbowPlot(nmp$gas, ndims = 30)
dim.n <- 10
# cluster cells to check batch effect
nmp$gas <- RunUMAP(nmp$gas, reduction = "pca", dims = 1:dim.n, return.model = T) %>%
  RunTSNE(reduction = "pca", dims = 1:dim.n)
# plot clustering results
pdf(file.path(outdir, "GASD5_SOX2_TBXT_cell_clustering.pdf"), height = 5, width = 6)
DimPlot(nmp$gas, reduction = "umap", group.by = c("CellType"),
        label = TRUE, repel = TRUE, pt.size = 1, cols = mk.col)
dev.off()
pdf(file.path(outdir, "GASD5_SOX2_TBXT_cell_clustering_with_human_CS7_annotation.pdf"), height = 5, width = 6)
DimPlot(nmp$gas, reduction = "umap", group.by = c("CellType.hs"),
        label = TRUE, repel = TRUE, pt.size = 1, cols = pd.col.fix)
dev.off()
pdf(file.path(outdir, "GASD5_SOX2_TBXT_cell_clustering_mapped_expression.pdf"), height = 5, width = 20)
FeaturePlot(nmp$gas, features = c("SOX2", "TBXT"), blend = T, max.cutoff = 0.5, pt.size = 1, cols = c("#E7E7E7", "#7b2c9c", "#f07f21"))
dev.off()
pdf(file.path(outdir, "GASD5_SOX2_TBXT_cell_clustering_mapped_density.pdf"), height = 5, width = 7)
Plot_Density_Joint_Only(seurat_object = nmp$gas, features = c("SOX2", "TBXT"), pt.size = 1)
dev.off()
pdf(file.path(outdir, "GASD5_SOX2_TBXT_cell_clustering_mapped_NMP_positive.pdf"), height = 5, width = 6)
DimPlot(nmp$gas, reduction = "umap", group.by = c("CellType.NMP"),
        label = TRUE, repel = TRUE, pt.size = 1.5, cols = mk.col)
dev.off()
pdf(file.path(outdir, "GASD5_SOX2_TBXT_cell_clustering_mapped_GFP_expr.pdf"), height = 5, width = 6)
FeaturePlot(nmp$gas, features = "GFP_count", pt.size = 1.5)
dev.off()
# average expression of positive and negative cell
nmp.expr.gasd5 <- AverageExpression(nmp$gas, group.by = "CellType.NMP")$RNA %>% as.data.frame()
# plot NMP cell composition
nmp$gas$SOX2 <- GetAssayData(nmp$gas, slot = "data")["SOX2", ]
nmp$gas$TBXT <- GetAssayData(nmp$gas, slot = "data")["TBXT", ]
nmp$gas$CellType.NMP <- "NMP.neg"
nmp$gas$CellType.NMP[nmp$gas$SOX2 > 0 & nmp$gas$TBXT > 0] <- "NMP.pos"
table(nmp$gas$CellType.NMP)
nmp$gas@meta.data$CellType <- as.character(nmp$gas@meta.data$CellType)
Vis.Annotation.Ratio(sr.meta = nmp$gas@meta.data, annotation = c("CellType.NMP", "CellType"),
                     pd.title = "Stat by Cell Type", pd.col = mk.col,
                     pd.height = 5, pd.width = 20,
                     res.out = file.path(outdir, "GASD5_composition"))
nmp$gas@meta.data$CellType.hs <- as.character(nmp$gas@meta.data$CellType.hs)
Vis.Annotation.Ratio(sr.meta = nmp$gas@meta.data, annotation = c("CellType.NMP", "CellType.hs"),
                     pd.title = "Stat by Cell Type", pd.col = pd.col.fix,
                     pd.height = 5, pd.width = 20,
                     res.out = file.path(outdir, "GASD5_composition_hs"))
Vis.Annotation.Ratio(sr.meta = nmp$gas@meta.data, annotation = c("CellType.NMP", "Type"),
                     pd.title = "Stat by Cell Type", pd.col = mk.col,
                     pd.height = 5, pd.width = 20,
                     res.out = file.path(outdir, "GASD5_composition_GFP"))
# plot expression
sr.out <- file.path(outdir, "Marker_genes_of_GASD5")
dir.create(sr.out, recursive = T)
DefaultAssay(nmp$gas) <- "RNA"
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(nmp$gas))
}
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(nmp$gas, features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 1, slot = "data", min.cutoff = 0.1))
    dev.off()
  }; rm(g)
}


### >>> 4. Correlation of NMP+ cells
expr.1 <- AverageExpression(subset(nmp$mk, CellType.NMP == "NMP.pos"), assays = "RNA", slot = "data", group.by = "CellType")
expr.2 <- AverageExpression(subset(nmp$gas, CellType.NMP == "NMP.pos"), assays = "RNA", slot = "data", group.by = "CellType")
hvg.1 <- VariableFeatures(FindVariableFeatures(subset(nmp$mk, CellType.NMP == "NMP.pos")))
hvg.2 <- VariableFeatures(FindVariableFeatures(subset(nmp$gas, CellType.NMP == "NMP.pos")))
homo.gene <- data.frame(species.1 = intersect(rownames(nmp$mk), rownames(nmp$gas)),
                        species.2 = intersect(rownames(nmp$mk), rownames(nmp$gas)))
CorrComparePlot(ExpressionTableSpecies1 = expr.1$RNA, DEgenesSpecies1 = hvg.1,
                ExpressionTableSpecies2 = expr.2$RNA, DEgenesSpecies2 = hvg.2,
                HomoGeneList = homo.gene, TFlist = hs.tfs$Symbol,
                filename = file.path(outdir, "Correlation_of_NMP_between_Monkey_and_GAS_D5"),
                pd.order = "original", Height = 15, Width = 15)


### >>> 5. Integrate monkey and gas
sr.nmp <- list()
sr.nmp$Mk.CS8 <- subset(nmp$mk, CellType.NMP == "NMP.pos")
sr.nmp$Mk.CS8$Datasets <- "Mk.CS8"
sr.nmp$GAS.D5 <- subset(nmp$gas, CellType.NMP == "NMP.pos")
sr.nmp$GAS.D5$Datasets <- "GAS.D5"
# integrating
sr.nmp <- lapply(X = sr.nmp, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, features = rownames(x), verbose = FALSE)
  x <- RunPCA(x, features = VariableFeatures(x), verbose = FALSE)
})
for (i in names(sr.nmp)) {
  sr.nmp[[i]]@project.name <- i
  print(ElbowPlot(sr.nmp[[i]]))
}
dim.n <- c(10, 10)
names(dim.n) <- names(sr.nmp)
sr.nmp <- lapply(X = sr.nmp, FUN = function(x) {
  dims <- unname(dim.n[x@project.name])
  x <- RunUMAP(x, dims = 1:dims, return.model = T)
})
sr.nmp.merge <- Reduce(function(x, y) merge(x, y, add.cell.ids = c(x@project.name, y@project.name)), sr.nmp)
# normalize data and find highly variables genes
sr.nmp.merge <- NormalizeData(sr.nmp.merge) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
# scale data
sr.nmp.merge <- ScaleData(sr.nmp.merge, verbose = FALSE, features = rownames(sr.nmp.merge))
# run PCA
sr.nmp.merge <- RunPCA(sr.nmp.merge, npcs = 50, verbose = FALSE, features = VariableFeatures(sr.nmp.merge))
ElbowPlot(sr.nmp.merge, ndims = 30)
dev.off()
dim.n <- 5
# cluster cells to check batch effect
sr.nmp.merge <- RunUMAP(sr.nmp.merge, reduction = "pca", dims = 1:dim.n, return.model = T) %>%
  RunTSNE(reduction = "pca", dims = 1:dim.n)
# remove batch effect
sr.nmp.merge <- RunHarmony(sr.nmp.merge, group.by.vars = c("Datasets"))
sr.nmp.merge <- RunUMAP(sr.nmp.merge, reduction = "harmony", dims = 1:dim.n) %>%
  RunTSNE(reduction = "harmony", dims = 1:dim.n)
# plot clustering results
pdf(file.path(res.out, "Integrated_clustering_of_NMP_cells.pdf"), height = 5, width = 12)
DimPlot(sr.nmp.merge, reduction = "umap", group.by = c("Datasets", "CellType"),
        label = TRUE, repel = TRUE, pt.size = 2, cols = mk.col)
dev.off()
# plot expression
sr.out <- file.path(outdir, "Marker_genes_of_merged_NMP")
dir.create(sr.out, recursive = T)
DefaultAssay(sr.nmp.merge) <- "RNA"
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(sr.nmp.merge))
}
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(sr.nmp.merge, features = paste0("rna_", marker.gene[[type]][g]), order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 2, slot = "data", min.cutoff = 0.1))
    dev.off()
  }; rm(g)
}


### >>> 6. Differential expression analysis
# find DEGs
table(sr.nmp.merge$Datasets, sr.nmp.merge$CellType)
deg.nmp <- list()
deg.cells <- c("AM", "Caud.Meso", "ECT", "Nas.Meso")
for (i in deg.cells) {
  tmp.sr <- subset(sr.nmp.merge, CellType == i)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Datasets", comparison = c("Mk.CS8", "GAS.D5"), sample.n = NULL)
  deg.nmp[[i]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                              sample.n = NULL, group.by = "CellType",
                              g1 = "1_Mk.CS8", g2 = "2_GAS.D5", lfc = 1, sig = 0.05,
                              res.out = file.path(outdir, paste0("DEGs/GAS.D5_vs_Mk.CS8_", i)))
}
tmp.sr <- sr.nmp.merge
tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Datasets", comparison = c("Mk.CS8", "GAS.D5"), sample.n = NULL)
deg.nmp[["All"]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                sample.n = NULL, group.by = "CellType",
                                g1 = "1_Mk.CS8", g2 = "2_GAS.D5", lfc = 1, sig = 0.05,
                                res.out = file.path(outdir, paste0("DEGs/GAS.D5_vs_Mk.CS8_All_cells")))
rm(i, tmp.sr, tmp)
# volcano plot
for (i in names(deg.nmp)) {
  pdf(file.path(outdir, paste0("DEGs/", i, "_GAS.D5_vs_Mk.CS8_volcano.pdf")), height = 6, width = 9)
  print(VisDEG.volcano(deg.data = deg.nmp[[i]]$all, geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.05, lfc = 1, title = i,
                       pt.size = 2, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
}
# GO
go.nmp <- list()
for (i in names(deg.nmp)) {
  go.nmp[[paste0(i, ".up")]] <- Pipe.GO(species = "human", genelist = deg.nmp[[i]]$up$SYMBOL,
                                         basename = paste0(i, "_up"), genetype = "SYMBOL",
                                         res.out = file.path(outdir, paste0("GO/", i)))
  go.nmp[[paste0(i, ".down")]] <- Pipe.GO(species = "human", genelist = deg.nmp[[i]]$down$SYMBOL,
                                          basename = paste0(i, "_down"), genetype = "SYMBOL",
                                          res.out = file.path(outdir, paste0("GO/", i)))
}


### >>> 7. Plot heatmap
# monkey data
expr.1 <- AverageExpression(subset(nmp$mk, CellType.NMP == "NMP.pos"), assays = "RNA", slot = "data", group.by = "CellType")
expr.1 <- expr.1$RNA %>% as.data.frame()
expr.1 <- expr.1[, deg.cells]
colnames(expr.1) <- paste0("Mk.", colnames(expr.1))
# gas data
expr.2 <- AverageExpression(subset(nmp$gas, CellType.NMP == "NMP.pos"), assays = "RNA", slot = "data", group.by = "CellType")
expr.2 <- expr.2$RNA %>% as.data.frame()
expr.2 <- expr.2[, deg.cells]
colnames(expr.2) <- paste0("Gas.", colnames(expr.2))
# gene list
gene.list <- list()
gene.list$TGFb <- read.table("/home/yhw/document/ensembl/GO/GO_Gene_Biomart_export_GO0007179_TGFbeta_receptor_signalling_pathway.txt")
gene.list$TGFb <- gene.list$TGFb$V2
gene.list$FGF <- read.table("/home/yhw/document/ensembl/GO/GO_Gene_Biomart_export_GO0008543_FGF_receptor_signaling_pathway.txt")
gene.list$FGF <- gene.list$FGF$V2
gene.list$WNT <- read.table("/home/yhw/document/ensembl/GO/GO_Gene_Biomart_export_GO0016055_Wnt_signaling_pathway.txt")
gene.list$WNT <- gene.list$WNT$V2
gene.list$BMP <- read.table("/home/yhw/document/ensembl/GO/GO_Gene_Biomart_export_GO0030509_BMP_signaling_pathway.txt")
gene.list$BMP <- gene.list$BMP$V2
gene.list$Hippo <- read.table("/home/yhw/document/ensembl/GO/GO_Gene_Biomart_export_GO0035329_hippo_signaling.txt")
gene.list$Hippo <- gene.list$Hippo$V2
# plot heatmap
common.gene <- intersect(rownames(expr.1), rownames(expr.2))
pd <- cbind(expr.1[common.gene, ], expr.2[common.gene, ])
pd.meta <- data.frame(SampleID = colnames(pd),
                      SampleGroup = str_split_fixed(colnames(pd), "\\.", 2)[, 1],
                      row.names = colnames(pd))
for (i in names(gene.list)) {
  pdf(file.path(outdir, paste0("NMP_pos_cell_gene_expression_level_in_", i, "_pathway.pdf")), height = 15, width = 15)
  print(Vis.heatmap(data = pd, meta = pd.meta, geneset = gene.list[[i]], genetype = NULL,
                    pt.value = "zscore",
                    show.rownames = T, show.colnames = T,
                    row.font.size = 8, col.font.size = 8,
                    pd.width = 15, pd.height = 15,
                    up.col = "#B71C1C", down.col = "#01579B"))
  dev.off()
}
rm(pd, pd.meta, i)


### >>> 8. Monkey CS11
mf.embryo <- Read10X("/home/yhw/document/public_data/GSE193007_monkey_post-gastrulation/GSE193007_MFE-filtered_feature_bc_matrix")
meta <- read.csv("/home/yhw/document/public_data/GSE193007_monkey_post-gastrulation/MFE56636-meta.csv", row.names = 1)
mf.embryo <- mf.embryo[,str_split_fixed(rownames(meta), "_", 2)[,2]]
rownames(meta) <- str_split_fixed(rownames(meta), "_", 2)[,2]
if (all(colnames(mf.embryo) == rownames(meta))) {
  mf.embryo <- CreateSeuratObject(counts = mf.embryo, meta.data = meta)
}; rm(meta)
mf.embryo.cs11 <- subset(mf.embryo, sample %in% c("CS11-e1", "CS11-e2"))
mf.embryo.cs11 <- NormalizeData(mf.embryo) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
mf.embryo.cs11 <- ScaleData(mf.embryo.cs11, verbose = FALSE, features = rownames(mf.embryo.cs11))
# average expression of positive and negative cell
nmp.expr.mkcs11 <- AverageExpression(mf.embryo.cs11, group.by = "cell_type")$RNA %>% as.data.frame()
rownames(nmp.expr.mkcs11)[grep("^T$", rownames(nmp.expr.mkcs11))] <- "TBXT"


### >>> 9. Human
hs.nmp <- readRDS("/home/yhw/document/public_data/GSE155121_human_gastrulation_organogenesis/GSE155121_human_data_raw.rds")
hs.nmp <- subset(hs.nmp, week_stage == "W3-1")
hs.nmp <- NormalizeData(hs.nmp) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
hs.nmp <- ScaleData(hs.nmp, verbose = FALSE)
hs.nmp$SOX2 <- GetAssayData(hs.nmp, slot = "data")["SOX2", ]
hs.nmp$TBXT <- GetAssayData(hs.nmp, slot = "data")["T", ]
hs.nmp$CellType.NMP <- "NMP.neg"
hs.nmp$CellType.NMP[hs.nmp$SOX2 > 0 & hs.nmp$TBXT > 0] <- "NMP.pos"
# average expression of positive and negative cell
nmp.expr.hs3w <- AverageExpression(hs.nmp, group.by = "CellType.NMP")$RNA %>% as.data.frame()
rownames(nmp.expr.hs3w)[grep("^T$", rownames(nmp.expr.hs3w))] <- "TBXT"


### >>> 10. NMP correlation analysis
common.gene <- Reduce(intersect, list(rownames(nmp.expr.gasd5), rownames(nmp.expr.mkcs11),
                                      rownames(nmp.expr.mkcs8), rownames(nmp.expr.hs3w)))
pd <- data.frame(GASD5 = nmp.expr.gasd5[common.gene, ]$NMP.pos,
                 Hs.3W = nmp.expr.hs3w[common.gene, ]$NMP.pos,
                 MK.CS8 = nmp.expr.mkcs8[common.gene, ]$NMP.pos,
                 MK.CS11 = nmp.expr.mkcs11[common.gene, ]$NMP)
rownames(pd) <- common.gene
label.gene <- intersect(c("MIXL1", "TBXT", "SOX2", "CDX1",
                          "CDX2", "WNT3A", "FGF4", "WNT5B",
                          "CDH2", "HES7", "SP5", "WNT8A",
                          "SP8", "WNT3", "TBX6", "MESP1",
                          "GDF3", "FOXA2", "EOMES", "HOXA1"), rownames(pd))
pd <- pd[c(setdiff(rownames(pd), label.gene), intersect(rownames(pd), label.gene)), ]
pd %>%
  mutate(type = c(rep("Others", nrow(pd)-length(label.gene)), rep("NMP", length(label.gene)))) %>%
  select("GASD5", "MK.CS11", "type") %>%
  rownames_to_column(var = "SYMBOL") %>%
  filter(SYMBOL %in% label.gene) -> pd.label
#pdf(file.path(outdir, "Correlation_analysis_of_NMP_between_Hs_CS7_GASD5_and_Human_embryo_3W.pdf"), height = 4.5, width = 5.5)
#pdf(file.path(outdir, "Correlation_analysis_of_NMP_between_Hs_CS7_GASD5_and_Monkey_embryo_CS8.pdf"), height = 4.5, width = 5.5)
#pdf(file.path(outdir, "Correlation_analysis_of_NMP_between_Hs_CS7_GASD5_and_Monkey_embryo_CS11.pdf"), height = 4.5, width = 5.5)
pd %>%
  mutate(type = c(rep("Others", nrow(pd)-length(label.gene)), rep("NMP", length(label.gene)))) %>%
  select("GASD5", "MK.CS11", "type") %>%
  ggplot(aes(x = log2(GASD5 + 0.1), y = log2(MK.CS11 + 0.1))) +
  geom_point(aes(color = type, size = type), shape = 16) +
  scale_color_manual(values = c("Others" = "#000000", "NMP" = "#E91E78")) +
  scale_size_manual(values = c("Others" = 1, "NMP" = 3)) +
  geom_smooth(method = "lm", se = F) +
  stat_cor(method = "pearson") +
  labs(x = "GAS.D5", y = "Human embryo") +
  geom_label_repel(aes(label = SYMBOL), data = pd.label,
                   box.padding = 0.2, point.padding = 0.2, segment.color = "#130f40",
                   label.size = 0.2, fill = NA, size = 3) +
  theme_bw() +
  theme(
    axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
    axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
    axis.text.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
    axis.text.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
    panel.grid = element_blank()
  )
dev.off()
# plot heatmap
common.gene <- Reduce(intersect, list(rownames(nmp.expr.gasd5), rownames(nmp.expr.mkcs11),
                                      rownames(nmp.expr.mkcs8), rownames(nmp.expr.hs3w)))
pd <- data.frame(GASD5 = nmp.expr.gasd5[common.gene, ]$NMP.pos,
                 Hs.3W = nmp.expr.hs3w[common.gene, ]$NMP.pos,
                 MK.CS8 = nmp.expr.mkcs8[common.gene, ]$NMP.pos,
                 MK.CS11 = nmp.expr.mkcs11[common.gene, ]$NMP)
rownames(pd) <- common.gene
pd.meta <- data.frame(SampleID = colnames(pd),
                      SampleGroup = colnames(pd),
                      row.names = colnames(pd))
for (i in names(gene.list)) {
  pdf(file.path(outdir, paste0("All_NMP_pos_cell_gene_expression_level_in_", i, "_pathway.pdf")), height = 15, width = 15)
  print(Vis.heatmap(data = pd, meta = pd.meta, geneset = gene.list[[i]], genetype = NULL,
                    pt.value = "zscore",
                    show.rownames = T, show.colnames = T,
                    row.font.size = 8, col.font.size = 8,
                    pd.width = 15, pd.height = 15,
                    up.col = "#B71C1C", down.col = "#01579B"))
  dev.off()
}
rm(pd, pd.meta, i)



### ========================
### 16th part: NSCs analysis ----
### ========================

### >>> 1. Setting output dir
outdir <- file.path(getwd(), "graphs/NSCs")
if (! dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}
nscs <- list()



### >>> 2. GASD5: HOX+ and SOX+
nscs.gene <- c("HOXB1", "HOXB2", "HOXB5", "HOXB6", "HOXB8", "HOXB9", "SOX2")
tmp <- sr.pdt$gasd5.to.hs
for (i in nscs.gene) {
  tmp@meta.data[, i] <- GetAssayData(tmp, slot = "data")[i, ]
}
keep.cell <- tmp$SOX2 > 0 & (tmp$HOXB1 > 0 | tmp$HOXB2 > 0 | tmp$HOXB5 > 0 | tmp$HOXB6 > 0 | tmp$HOXB8 > 0 | tmp$HOXB9 > 0)
nscs$hs.gasd5 <- tmp[, keep.cell]
# normalize data and find highly variables genes
nscs$hs.gasd5 <- NormalizeData(nscs$hs.gasd5) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
# scale data
nscs$hs.gasd5 <- ScaleData(nscs$hs.gasd5, verbose = FALSE, features = rownames(nscs$hs.gasd5))
# run PCA
nscs$hs.gasd5 <- RunPCA(nscs$hs.gasd5, npcs = 50, verbose = FALSE, features = VariableFeatures(nscs$hs.gasd5))
ElbowPlot(nscs$hs.gasd5, ndims = 30)
dim.n <- 5
# cluster cells to check batch effect
nscs$hs.gasd5 <- RunUMAP(nscs$hs.gasd5, reduction = "pca", dims = 1:dim.n, return.model = T) %>%
  RunTSNE(reduction = "pca", dims = 1:dim.n)
# plot clustering results
pdf(file.path(outdir, "GASD5_SOX2_HOXBs_cell_clustering.pdf"), height = 5, width = 6)
DimPlot(nscs$hs.gasd5, reduction = "umap", group.by = c("CellType"),
        label = TRUE, repel = TRUE, pt.size = 3, cols = pd.col.fix)
dev.off()
# plot density
library(scCustomize)
for (i in nscs.gene[1:6]) {
  pdf(file.path(outdir, paste0("Scatter_plot_to_show_NSCs_gene_co-expression_SOX2_and_", i, ".pdf")),
      height = 3.5, width = 5)
  print(Plot_Density_Joint_Only(seurat_object = nscs$hs.gasd5, features = c("SOX2", i), pt.size = 3))
  dev.off()
}
# plot expression
sr.out <- file.path(outdir, "Marker_genes_of_GASD5")
dir.create(sr.out, recursive = T)
DefaultAssay(nscs$hs.gasd5) <- "RNA"
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(nscs$hs.gasd5))
}
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(nscs$hs.gasd5, features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 3, slot = "data", min.cutoff = 0.1))
    dev.off()
  }; rm(g)
}


### >>> 3. GASD5: PAX6+
nscs.gene <- c("PAX6")
tmp <- sr.pdt$gasd5.to.hs
for (i in nscs.gene) {
  tmp@meta.data[, i] <- GetAssayData(tmp, slot = "data")[i, ]
}
keep.cell <- tmp$PAX6 > 0
nscs$pax6 <- tmp[, keep.cell]
# normalize data and find highly variables genes
nscs$pax6 <- NormalizeData(nscs$pax6) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
# scale data
nscs$pax6 <- ScaleData(nscs$pax6, verbose = FALSE, features = rownames(nscs$pax6))
# run PCA
nscs$pax6 <- RunPCA(nscs$pax6, npcs = 20, verbose = FALSE, features = VariableFeatures(nscs$pax6))
ElbowPlot(nscs$pax6, ndims = 30)
dim.n <- 10
# cluster cells to check batch effect
nscs$pax6 <- RunUMAP(nscs$pax6, reduction = "pca", dims = 1:dim.n, return.model = T)
# plot clustering results
pdf(file.path(outdir, "GASD5_PAX6_cell_clustering.pdf"), height = 5, width = 6)
DimPlot(nscs$pax6, reduction = "umap", group.by = c("CellType"),
        label = TRUE, repel = TRUE, pt.size = 4, cols = pd.col.fix) +
  xlim(c(-4, 4)) + ylim(c(-3, 3))
dev.off()
# plot expression
sr.out <- file.path(outdir, "Marker_genes_of_GASD5_PAX6")
dir.create(sr.out, recursive = T)
DefaultAssay(nscs$pax6) <- "RNA"
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(nscs$pax6))
}
Idents(nscs$pax6) <- "PAX6.pos"
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    #pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
    #    height = 4.5, width = 5)
    #print(FeaturePlot(nscs$pax6, features = marker.gene[[type]][g], order = T,
    #                  cols = c("#ebebeb", "#7d3c98"), pt.size = 3, slot = "data", min.cutoff = 0.1) +
    #        xlim(c(-4, 4)) + ylim(c(-3, 3)))
    #dev.off()
    # violin plot
    pdf(file.path(tmp.out, paste0("Violin_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 4.5)
    print(VlnPlot(nscs$pax6, features = marker.gene[[type]][g], pt.size = 3))
    dev.off()
  }; rm(g)
}


### >>> 4. GASD5: PAX3+
nscs.gene <- c("PAX3")
tmp <- sr.pdt$gasd5.to.hs
for (i in nscs.gene) {
  tmp@meta.data[, i] <- GetAssayData(tmp, slot = "data")[i, ]
}
keep.cell <- tmp$PAX3 > 0
table(keep.cell)
nscs$pax3 <- tmp[, keep.cell]
# normalize data and find highly variables genes
nscs$pax3 <- NormalizeData(nscs$pax3) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
# scale data
nscs$pax3 <- ScaleData(nscs$pax3, verbose = FALSE, features = rownames(nscs$pax3))
# run PCA
nscs$pax3 <- RunPCA(nscs$pax3, npcs = 20, verbose = FALSE, features = VariableFeatures(nscs$pax3))
ElbowPlot(nscs$pax3, ndims = 30)
dim.n <- 5
# cluster cells to check batch effect
nscs$pax3 <- RunUMAP(nscs$pax3, reduction = "pca", dims = 1:dim.n, return.model = T)
# plot clustering results
pdf(file.path(outdir, "GASD5_PAX3_cell_clustering.pdf"), height = 5, width = 6)
DimPlot(nscs$pax3, reduction = "umap", group.by = c("CellType"),
        label = TRUE, repel = TRUE, pt.size = 3, cols = pd.col.fix) +
  xlim(c(-4, 4)) + ylim(c(-4, 4))
dev.off()
# plot expression
sr.out <- file.path(outdir, "Marker_genes_of_GASD5_PAX3")
dir.create(sr.out, recursive = T)
DefaultAssay(nscs$pax3) <- "RNA"
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(nscs$pax3))
}
Idents(nscs$pax3) <- "PAX3.pos"
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    #pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
    #    height = 4.5, width = 5)
    #print(FeaturePlot(nscs$pax3, features = marker.gene[[type]][g], order = T,
    #                  cols = c("#ebebeb", "#7d3c98"), pt.size = 3, slot = "data", min.cutoff = 0.1) +
    #        xlim(c(-4, 4)) + ylim(c(-3, 3)))
    #dev.off()
    # violin plot
    pdf(file.path(tmp.out, paste0("Violin_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 4.5)
    print(VlnPlot(nscs$pax3, features = marker.gene[[type]][g], pt.size = 3))
    dev.off()
  }; rm(g)
}


### >>> 4. GASD5: PAX7+
nscs.gene <- c("PAX7")
tmp <- sr.pdt$gasd5.to.hs
for (i in nscs.gene) {
  tmp@meta.data[, i] <- GetAssayData(tmp, slot = "data")[i, ]
}
keep.cell <- tmp$PAX7 > 0
table(keep.cell)
nscs$pax7 <- tmp[, keep.cell]
# normalize data and find highly variables genes
nscs$pax7 <- NormalizeData(nscs$pax7) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
# scale data
nscs$pax7 <- ScaleData(nscs$pax7, verbose = FALSE, features = rownames(nscs$pax7))
# run PCA
nscs$pax7 <- RunPCA(nscs$pax7, npcs = 10, verbose = FALSE, features = VariableFeatures(nscs$pax7))
ElbowPlot(nscs$pax7, ndims = 30)
dim.n <- 5
# cluster cells to check batch effect
nscs$pax7 <- RunUMAP(nscs$pax7, reduction = "pca", dims = 1:dim.n, return.model = T)
# plot clustering results
pdf(file.path(outdir, "GASD5_PAX7_cell_clustering.pdf"), height = 5, width = 6)
DimPlot(nscs$pax7, reduction = "umap", group.by = c("CellType"),
        label = TRUE, repel = TRUE, pt.size = 3, cols = pd.col.fix) +
  xlim(c(-20, 15)) + ylim(c(-10, 15))
dev.off()
# plot expression
sr.out <- file.path(outdir, "Marker_genes_of_GASD5_PAX7")
dir.create(sr.out, recursive = T)
DefaultAssay(nscs$pax7) <- "RNA"
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(nscs$pax7))
}
Idents(nscs$pax7) <- "PAX7.pos"
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(nscs$pax7, features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 3, slot = "data", min.cutoff = 0.1) +
            xlim(c(-20, 15)) + ylim(c(-10, 15)))
    dev.off()
    # violin plot
    pdf(file.path(tmp.out, paste0("Violin_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 4.5)
    print(VlnPlot(nscs$pax7, features = marker.gene[[type]][g], pt.size = 3))
    dev.off()
  }; rm(g)
}



### =============================
### 17th part: FOXH1 KO scRNA-seq ----
### =============================


### >>> 1. Setting output dir
res.out <- file.path(getwd(), "graphs/FOXH1_scRNA-seq")
if (! dir.exists(res.out)) {
  dir.create(res.out, recursive = T)
}


### >>> 2. Load data
# count data
ge.co <- read.csv("foxh1/output/Combined_St230806_RSEC_ReadsPerCell_for_R.csv", header = T, row.names = 1)
ge.co <- t(ge.co)
# cell meta
cell.meta <- read.csv("foxh1/output/St230806_Sample_Tag_Calls_for_R.csv", stringsAsFactors = F, row.names = 1)


### >>> 3. Create Seurat object
if (all(colnames(ge.co) == rownames(cell.meta))) {
  sr.foxh1 <- CreateSeuratObject(counts = ge.co, project = "FOXH1", meta.data = cell.meta) %>%
    PercentageFeatureSet(col.name = "Pect.mt", pattern = "^MT\\.")
}
pdf(file.path(outdir, "FOXH1_KO_quality_raw_violin_plot.pdf"), height = 4, width = 10)
VlnPlot(sr.foxh1, features = c("nFeature_RNA", "nCount_RNA", "Pect.mt"), ncol = 3, 
        group.by = "Sample_Name", pt.size = 0, cols = pd.col)
dev.off()


### >>> 4. Split data
sr.foxh1.org <- subset(sr.foxh1, Sample_Name %in% c("ORG_KO", "ORG_WT"))
sr.foxh1.gas <- subset(sr.foxh1, Sample_Name %in% c("GAS_KO", "GAS_WT"))


### >>> 5. Filter data
# GAS
sr.foxh1.gas$GFP_count <- GetAssayData(sr.foxh1.gas, slot = "count")["GFP",]
sr.foxh1.gas <- subset(sr.foxh1.gas, nFeature_RNA >= 500 & nCount_RNA >= 5000 & Pect.mt <= 25)
pdf(file.path(res.out, "FOXH1_KO_GAS_quality_after_filtering_violin_plot.pdf"), height = 4, width = 5)
VlnPlot(sr.foxh1.gas, features = c("nFeature_RNA", "nCount_RNA", "Pect.mt"), ncol = 3, 
        group.by = "Sample_Name", pt.size = 0, cols = pd.col)
dev.off()
# ORG
sr.foxh1.org$GFP_count <- GetAssayData(sr.foxh1.org, slot = "count")["GFP",]
sr.foxh1.org <- subset(sr.foxh1.org, nFeature_RNA >= 500 & nCount_RNA >= 5000 & Pect.mt <= 25)
pdf(file.path(res.out, "FOXH1_KO_ORG_quality_after_filtering_violin_plot.pdf"), height = 4, width = 5)
VlnPlot(sr.foxh1.org, features = c("nFeature_RNA", "nCount_RNA", "Pect.mt"), ncol = 3, 
        group.by = "Sample_Name", pt.size = 0, cols = pd.col)
dev.off()


### >>> 6. Process by Seurat
# GAS
sr.foxh1.gas <- NormalizeData(sr.foxh1.gas) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sr.foxh1.gas))
sr.foxh1.gas <- RunPCA(sr.foxh1.gas, features = VariableFeatures(object = sr.foxh1.gas))
ElbowPlot(sr.foxh1.gas, ndims = 30)
dev.off()
sr.foxh1.gas <- RunUMAP(sr.foxh1.gas, dims = 1:20)
pdf(file.path(res.out, "FOXH1_KO_GAS_clustering_umap_plot_before_annotation.pdf"), height = 4, width = 5)
DimPlot(sr.foxh1.gas, reduction = "umap", group.by = c("Sample_Name"), cols = pd.col)
dev.off()
sr.foxh1.gas <- RunHarmony(sr.foxh1.gas, group.by.vars = "Sample_Name")
sr.foxh1.gas <- RunUMAP(sr.foxh1.gas, reduction = "harmony", dims = 1:20)
pdf(file.path(res.out, "FOXH1_KO_GAS_clustering_umap_plot_before_annotation_after_correction.pdf"), height = 4, width = 5)
DimPlot(sr.foxh1.gas, reduction = "umap", group.by = c("Sample_Name"), cols = pd.col)
dev.off()
# ORG
sr.foxh1.org <- NormalizeData(sr.foxh1.org) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sr.foxh1.org))
sr.foxh1.org <- RunPCA(sr.foxh1.org, features = VariableFeatures(object = sr.foxh1.org))
ElbowPlot(sr.foxh1.org, ndims = 30)
dev.off()
sr.foxh1.org <- RunUMAP(sr.foxh1.org, dims = 1:10)
pdf(file.path(res.out, "FOXH1_KO_ORG_clustering_umap_plot_before_annotation.pdf"), height = 4, width = 5)
DimPlot(sr.foxh1.org, reduction = "umap", group.by = c("Sample_Name"), cols = pd.col)
dev.off()
sr.foxh1.org <- RunHarmony(sr.foxh1.org, group.by.vars = "Sample_Name")
sr.foxh1.org <- RunUMAP(sr.foxh1.org, reduction = "harmony", dims = 1:10)
pdf(file.path(res.out, "FOXH1_KO_ORG_clustering_umap_plot_before_annotation_after_correction.pdf"), height = 4, width = 5)
DimPlot(sr.foxh1.org, reduction = "umap", group.by = c("Sample_Name"), cols = pd.col)
dev.off()


### >>> 7. Annotate cells by integration
# GAS
anchors <- FindTransferAnchors(reference = hs.cs7, query = sr.foxh1.gas, 
                               dims = 1:25, reference.reduction = "pca")
sr.pdt$foxh1.to.hs <- MapQuery(anchorset = anchors, reference = hs.cs7, query = sr.foxh1.gas,
                               refdata = list(celltype = "CellType"), 
                               reference.reduction = "pca", reduction.model = "umap")
sr.pdt$foxh1.to.hs@meta.data <- sr.pdt$foxh1.to.hs@meta.data %>%
  mutate(CellType = factor(predicted.celltype, levels = names(pd.col.fix)))
# umap and tsne plot
p1 <- DimPlot(hs.cs7, reduction = "umap", group.by = "CellType", label = F, 
              pt.size = 1, cols = pd.col.fix) + ggtitle("Human embryos") + theme(aspect.ratio=4/7)
p2 <- DimPlot(sr.pdt$foxh1.to.hs, reduction = "ref.umap", group.by = "CellType", 
              label = F, pt.size = 1, cols = pd.col.fix) + ggtitle("FOXH1") + theme(aspect.ratio=4/7)
p3 <- DimPlot(sr.pdt$foxh1.to.hs, reduction = "ref.umap", group.by = "Sample_Name", 
              label = F, pt.size = 1, cols = pd.col) + ggtitle("FOXH1") + theme(aspect.ratio=4/7)
pdf(file.path(res.out, "GAS_FOXH1_and_human_CS7_embryo_visualized_separately_umap.pdf"),
    height = 3.5, width = 15)
print(p1 + p2 + p3)
dev.off()


### >>> 7. Annotate cells by hands
# GAS
sr.foxh1.gas <- FindNeighbors(sr.foxh1.gas, reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = seq(0, 3, 0.5))
if (all(colnames(sr.foxh1.gas) == colnames(sr.pdt$foxh1.to.hs))) {
  sr.foxh1.gas$predicted.CellType <- sr.pdt$foxh1.to.hs$CellType
}
pdf(file.path(res.out, "FOXH1_KO_GAS_clustering_umap_plot_before_annotation_after_correction.pdf"), height = 5, width = 12)
DimPlot(sr.foxh1.gas, reduction = "umap", group.by = c("Sample_Name", "RNA_snn_res.2", "predicted.CellType"), 
        cols = mk.col, label = T, pt.size = 1)
dev.off()
pdf(file.path(res.out, "FOXH1_KO_GAS_marker_gene_Dot_plot_before_annotation_after_correction.pdf"), height = 8, width = 15)
DotPlot(sr.foxh1.gas, features = unique(pd.ge), cols = c("#eaeaea", "#fc0330"),
        col.min = 0, col.max = 2, group.by = "RNA_snn_res.2", assay = "RNA",
        scale = T, scale.min = 0, scale.max = 100) +
  RotatedAxis() +
  scale_x_discrete(labels = unique(pd.ge)) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
dev.off()
# EPI --> cluster 0,1,2,4,5,6,7,8,9,10,13,16
FeaturePlot(sr.foxh1.gas, features = c("DLX5", "TFAP2A", "TFAP2C", "GATA3"), ncol = 4)
# PS --> cluster 11,20
FeaturePlot(sr.foxh1.gas, features = c("TBXT", "MIXL1", "SP5", "FST", "EOMES"), ncol = 4)
# ECT --> cluster 3,12
FeaturePlot(sr.foxh1.gas, features = c("DLX5", "TFAP2A", "TFAP2C", "GATA3"), ncol = 4)
# END --> cluster 17
FeaturePlot(sr.foxh1.gas, features = c("SOX17", "FOXA2", "GSC", "CST1", "GATA6"), ncol = 4)
# NM --> cluster 21
FeaturePlot(sr.foxh1.gas, features = c("TBX6", "MESP1", "EOMES", "LHX1", "IRX3", 
                                       "GATA6", "FOXF1", "HAND1", "HAND2"), ncol = 4)
# EM --> cluster 15
FeaturePlot(sr.foxh1.gas, features = c("LHX1", "IRX3", "GATA6"), ncol = 4)
# ExM --> cluster 19
FeaturePlot(sr.foxh1.gas, features = c("POSTN", "LUM", "NID2", "HAND1", "ANXA1"), ncol = 4)
# HEP --> cluster 14
FeaturePlot(sr.foxh1.gas, features = c("MEF2C", "PECAM1", "CDH5", "TEK"), ncol = 4)
# plot expression
sr.out <- file.path(res.out, "Marker_genes_of_GAS")
dir.create(sr.out, recursive = T)
DefaultAssay(sr.foxh1.gas) <- "RNA"
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(sr.foxh1.gas))
}
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(sr.foxh1.gas, features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 1, slot = "data", min.cutoff = 0.1))
    dev.off()
  }; rm(g)
}
# rename cells
sr.foxh1.gas@meta.data <- sr.foxh1.gas@meta.data %>% 
  mutate(CellType = case_when(RNA_snn_res.2 %in% c(2) ~ "neural crest progenitors",
                              RNA_snn_res.2 %in% c(14) ~ "migratory neural crest",
                              RNA_snn_res.2 %in% c(20) ~ "primitive streak",
                              RNA_snn_res.2 %in% c(11) ~ "epiblast",
                              RNA_snn_res.2 %in% c(18,7,16,4,8,0,10,13,5,9,6,1) ~ "neural stem cells",
                              RNA_snn_res.2 %in% c(12,3,22) ~ "non-neural ectoderm",
                              RNA_snn_res.2 %in% c(19) ~ "amnion ectoderm",
                              RNA_snn_res.2 %in% c(17) ~ "hypoblast/yolk-sac endoderm",
                              RNA_snn_res.2 %in% c(15) ~ "definitive endoderm",
                              RNA_snn_res.2 %in% c(21) ~ "extraembryonic mesenchymal cells"))
pdf(file.path(res.out, "FOXH1_KO_GAS_clustering_umap_plot_after_annotation_after_correction.pdf"), height = 5, width = 20)
library(Seurat)
DimPlot(sr.foxh1.gas, reduction = "umap", group.by = c("Sample_Name", "RNA_snn_res.2", "CellType"), 
        cols = mk.col, label = T, pt.size = 1)
dev.off()
Vis.Annotation.Ratio(sr.meta = sr.foxh1.gas@meta.data, 
                     annotation = c("CellType", "Sample_Name"),
                     pd.title = "Stat by Group", pd.col = NULL,
                     pd.height = 15, pd.width = 15,
                     res.out = file.path(res.out, "Composition"))
# ORG
sr.foxh1.org <- FindNeighbors(sr.foxh1.org, reduction = "harmony", dims = 1:10) %>%
  FindClusters(resolution = seq(0, 3, 0.5))
pdf(file.path(res.out, "FOXH1_KO_ORG_clustering_umap_plot_before_annotation_after_correction.pdf"), height = 5, width = 12)
DimPlot(sr.foxh1.org, reduction = "umap", group.by = c("Sample_Name", "RNA_snn_res.2"), 
        cols = mk.col, label = T, pt.size = 1)
dev.off()
pdf(file.path(res.out, "FOXH1_KO_ORG_marker_gene_Dot_plot_before_annotation_after_correction.pdf"), height = 8, width = 15)
DotPlot(sr.foxh1.org, features = unique(pd.ge), cols = c("#eaeaea", "#fc0330"),
        col.min = 0, col.max = 2, group.by = "RNA_snn_res.2", assay = "RNA",
        scale = T, scale.min = 0, scale.max = 100) +
  RotatedAxis() +
  scale_x_discrete(labels = unique(pd.ge)) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
dev.off()
# plot expression
sr.out <- file.path(res.out, "Marker_genes_of_ORG")
dir.create(sr.out, recursive = T)
DefaultAssay(sr.foxh1.org) <- "RNA"
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(sr.foxh1.org))
}
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(sr.foxh1.org, features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 1, slot = "data", min.cutoff = 0.1))
    dev.off()
  }; rm(g)
}
# integrate with GAS before
sr.foxh1.org.merge <- merge(sr.foxh1.org, y = gas.d5before.merge, add.cell.ids = c("B2", "B1"), project = "ORG")
sr.foxh1.org.merge <- NormalizeData(sr.foxh1.org.merge) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sr.foxh1.org.merge))
sr.foxh1.org.merge <- RunPCA(sr.foxh1.org.merge, features = VariableFeatures(object = sr.foxh1.org.merge))
ElbowPlot(sr.foxh1.org.merge, ndims = 30)
dev.off()
sr.foxh1.org.merge <- RunUMAP(sr.foxh1.org.merge, dims = 1:15)
pdf(file.path(res.out, "FOXH1_KO_ORG_merged_GAS_clustering_umap_plot_before_annotation.pdf"), height = 4, width = 5)
DimPlot(sr.foxh1.org.merge, reduction = "umap", group.by = c("Sample_Name"), cols = pd.col)
dev.off()
sr.foxh1.org.merge <- RunHarmony(sr.foxh1.org.merge, group.by.vars = "Sample_Name")
sr.foxh1.org.merge <- RunUMAP(sr.foxh1.org.merge, reduction = "harmony", dims = 1:15)
sr.foxh1.org.merge <- FindNeighbors(sr.foxh1.org.merge, reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = seq(0, 3, 0.5))
pdf(file.path(res.out, "FOXH1_KO_ORG_merged_GAS_clustering_umap_plot_before_annotation_after_correction.pdf"), height = 5, width = 18)
DimPlot(sr.foxh1.org.merge, reduction = "umap", group.by = c("Sample_Name", "CellType", "RNA_snn_res.3"), cols = mk.col, pt.size = 0.25)
dev.off()
sr.out <- file.path(res.out, "Marker_genes_of_ORG_merged_GAS")
dir.create(sr.out, recursive = T)
DefaultAssay(sr.foxh1.org.merge) <- "RNA"
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(sr.foxh1.org.merge))
}
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(sr.foxh1.org.merge, features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 0.25, slot = "data", min.cutoff = 'q25'))
    dev.off()
  }; rm(g)
}
# integrate with ORG before
table(sr.list$cl$Sample_Name)
sr.foxh1.org.merge2 <- merge(sr.foxh1.org, y = sr.list$cl, add.cell.ids = c("B1", "B2"), project = "ORG")
sr.foxh1.org.merge2 <- NormalizeData(sr.foxh1.org.merge2) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sr.foxh1.org.merge2))
sr.foxh1.org.merge2 <- RunPCA(sr.foxh1.org.merge2, features = VariableFeatures(object = sr.foxh1.org.merge2))
ElbowPlot(sr.foxh1.org.merge2, ndims = 30)
dev.off()
sr.foxh1.org.merge2 <- RunUMAP(sr.foxh1.org.merge2, dims = 1:10)
pdf(file.path(res.out, "FOXH1_KO_ORG_merged_ORG_clustering_umap_plot_before_annotation.pdf"), height = 4, width = 5)
DimPlot(sr.foxh1.org.merge2, reduction = "umap", group.by = c("Sample_Name"), cols = pd.col)
dev.off()
sr.foxh1.org.merge2 <- RunHarmony(sr.foxh1.org.merge2, group.by.vars = "Sample_Name")
sr.foxh1.org.merge2 <- RunUMAP(sr.foxh1.org.merge2, reduction = "harmony", dims = 1:10)
sr.foxh1.org.merge2 <- FindNeighbors(sr.foxh1.org.merge2, reduction = "harmony", dims = 1:10) %>%
  FindClusters(resolution = seq(0, 3, 0.5))
pdf(file.path(res.out, "FOXH1_KO_ORG_merged_ORG_clustering_umap_plot_before_annotation_after_correction.pdf"), height = 5, width = 13)
DimPlot(sr.foxh1.org.merge2, reduction = "umap", group.by = c("Sample_Name", "RNA_snn_res.3"), 
        cols = mk.col, pt.size = 0.25, label = T)
dev.off()
sr.out <- file.path(res.out, "Marker_genes_of_ORG_merged_ORG")
dir.create(sr.out, recursive = T)
DefaultAssay(sr.foxh1.org.merge2) <- "RNA"
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(sr.foxh1.org.merge2))
}
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeaturePlot(sr.foxh1.org.merge2, features = marker.gene[[type]][g], order = T,
                      cols = c("#ebebeb", "#7d3c98"), pt.size = 0.25, slot = "data", min.cutoff = 'q25'))
    dev.off()
  }; rm(g)
}



### ====================
### 18th part: Save data ----
### ====================

saveRDS(sr.list, "/home/yhw/bioinfo/project-wb/CodeData/sr.list.rds")
saveRDS(sr.list.org, "/home/yhw/bioinfo/project-wb/CodeData/sr.list.org.rds")
saveRDS(sr.pdt, "/home/yhw/bioinfo/project-wb/CodeData/sr.pdt.rds")
saveRDS(sr.pdt.hs.mk, "/home/yhw/bioinfo/project-wb/CodeData/sr.pdt.hs.mk.rds")
saveRDS(sr.pdt.mk, "/home/yhw/bioinfo/project-wb/CodeData/sr.pdt.mk.rds")
saveRDS(sr.merge, "/home/yhw/bioinfo/project-wb/CodeData/sr.merge.rds")
saveRDS(sr.endo.merge, "/home/yhw/bioinfo/project-wb/CodeData/sr.endo.merge.rds")
saveRDS(sr.meso.merge, "/home/yhw/bioinfo/project-wb/CodeData/sr.meso.merge.rds")
saveRDS(sr.endo, "/home/yhw/bioinfo/project-wb/CodeData/sr.endo.rds")
saveRDS(sr.meso, "/home/yhw/bioinfo/project-wb/CodeData/sr.meso.rds")
saveRDS(sr.gas, "/home/yhw/bioinfo/project-wb/CodeData/sr.gas.rds")
saveRDS(sr.org, "/home/yhw/bioinfo/project-wb/CodeData/sr.org.rds")
saveRDS(sr.foxh1, "/home/yhw/bioinfo/project-wb/CodeData/sr.foxh1.rds")
saveRDS(sr.foxh1.gas, "/home/yhw/bioinfo/project-wb/CodeData/sr.foxh1.gas.rds")
saveRDS(gas.merge, "/home/yhw/bioinfo/project-wb/CodeData/gas.merge.rds")
saveRDS(gas.d5before.merge, "/home/yhw/bioinfo/project-wb/CodeData/gas.d5before.merge.rds")
saveRDS(gas.d5before.merge.urd, "/home/yhw/bioinfo/project-wb/CodeData/gas.d5before.merge.urd.rds")
saveRDS(gene.count$foxh1.ko, "/home/yhw/bioinfo/project-wb/CodeData/foxh1.ko.RNAseq.count.rds")
save.image("/home/yhw/bioinfo/project-wb/CodeData/WB_downstream_v1.RData")
