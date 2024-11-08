

library(sigminer)
library(tidyverse)

maf_xena <- read_maf("~/October_2016_whitelist_2583.snv_mnv_indel.maf.coding.xena")

tab <- read_table(file = "~/October_2016_whitelist_2583.snv_mnv_indel.maf.coding.xena")

head(tab, 10)
table(tab$effect)


load("~/RNAsignatures/data/mergedRNAMAF.RData")
head(as.data.frame(mergedRNAmaf))

table(tab$end - tab$start)
table(mergedRNAmaf$Variant_Type)



table(mergedRNAmaf$End_Position - mergedRNAmaf$Start_Position,mergedRNAmaf$Variant_Type)

mergedRNAmaf


table(tab$effect)
table(mergedRNAmaf$Variant_Classification, mergedRNAmaf$Variant_Type)


maf_tab <- tab %>%
  mutate(Variant_Classification = effect, 
         Chromosome = paste0("chr", chr), 
         Hugo_Symbol = gene, 
         Start_Position = start,
         End_Position = end, 
         Strand = "*", 
         Tumor_Sample_Barcode = Sample, 
         Reference_Allele = reference, 
         Tumor_Seq_Allele1 = alt, 
         Tumor_Seq_Allele2 = alt) %>%
  select(Variant_Classification, Chromosome, Hugo_Symbol, 
         Start_Position, End_Position, Strand, Tumor_Sample_Barcode, 
         Reference_Allele, 
         Tumor_Seq_Allele1, 
         Tumor_Seq_Allele2)

maf_tab2 <- read_maf(maf_tab)


dd <- read_xena_variants("~/October_2016_whitelist_2583.snv_mnv_indel.maf.coding.xena", )


dt <- data.table::fread("~/October_2016_whitelist_2583.snv_mnv_indel.maf.coding.xena")
detect_name <- function(x, y, z) {
  if (x %in% z) 
    x
  else y
}
data.table::setnames(dt, 
                     old = c(detect_name("Sample_ID", "Sample", colnames(dt)), 
                                 "gene", 
                                 detect_name("chrom", "chr", colnames(dt)), 
                                 "start", 
                                 "end", 
                                 detect_name("ref", "reference", colnames(dt)), 
                                 "alt"), 
                     new = c("Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2"))
head(dt)






tally_maf <- sig_tally(df_maf, mode = "ALL")
dim(tally_maf$SBS_96)

colSums(tally_maf$ID_83)




### Read the MAF from Xena
load_maf_from_xena <- function(file) {
  dt <- data.table::fread(file)
  detect_name <- function(x, y, z) {
    if (x %in% z) 
      x
    else y
  }
  data.table::setnames(dt, 
                       old = c(detect_name("Sample_ID", "Sample", colnames(dt)), 
                               "gene", 
                               detect_name("chrom", "chr", colnames(dt)), 
                               "start", 
                               "end", 
                               detect_name("ref", "reference", colnames(dt)), 
                               "alt"), 
                       new = c("Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2"))
  dt$Variant_Type <- dplyr::case_when(nchar(dt$Reference_Allele) == 
                                        1L & nchar(dt$Tumor_Seq_Allele2) == 1L ~ "SNP", nchar(dt$Reference_Allele) < 
                                        nchar(dt$Tumor_Seq_Allele2) ~ "INS", nchar(dt$Reference_Allele) > 
                                        nchar(dt$Tumor_Seq_Allele2) ~ "DEL", nchar(dt$Reference_Allele) == 
                                        2L & nchar(dt$Tumor_Seq_Allele2) == 2L ~ "DNP", nchar(dt$Reference_Allele) == 
                                        3L & nchar(dt$Tumor_Seq_Allele2) == 3L ~ "TNP", TRUE ~ 
                                        "Unknown")
  dt$Variant_Classification <- "Unknown"
  dt$Hugo_Symbol <- "Unknown"
  df_maf <- maftools::read.maf(dt, clinicalData = NULL, removeDuplicatedVariants = TRUE, 
                               useAll = TRUE, gisticAllLesionsFile = NULL, gisticAmpGenesFile = NULL, 
                               gisticDelGenesFile = NULL, gisticScoresFile = NULL, cnLevel = "all", 
                               cnTable = NULL, isTCGA = FALSE, vc_nonSyn = "Unknown",
                               verbose = FALSE)
  return(df_maf)
}

data <- load_maf_from_xena("~/October_2016_whitelist_2583.snv_mnv_indel.maf.coding.xena")

tally_maf <- sig_tally(df_maf, mode = "ID", add_trans_bias = TRUE)
tally_maf$SBS_96

tally_maf_coding <- sig_tally(data, mode = "ID")
dim(tally_maf_coding$nmf_matrix)
colSums(tally_maf_coding$nmf_matrix)






colSums(tally_maf$SBS_96)

tumor_data <- read_tsv("~/pcawg_specimen_histology_August2016_v9_donor") %>% as.data.frame()
project_codes <- read_table("~/project_code_donor") %>% as.data.frame()




brca_id <- tumor_data$icgc_specimen_id[tumor_data$histology_abbreviation == "Breast-AdenoCA" & !is.na(tumor_data$histology_abbreviation)]

id_brca <- rownames(tally_maf$ID_83) %in% brca_id

Xindels <- t(tally_maf$ID_83[id_brca, ])


indel_nmf <- sig_auto_extract(t(Xindels))

indel_nmf$Signature.norm
show_sig_profile(indel_nmf$Signature.norm, mode = "ID", style = "cosmic", x_label_angle = 90)


  rna.maf <- read_maf(mergedRNAmaf)
  mt_tally <- sig_tally(
    rna.maf,mode = "ID",
    ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
    useSyn = TRUE,
    #add_trans_bias = TRUE
  )

  colSums(mt_tally$nmf_matrix)

  indel_nmf_RNA <- sig_auto_extract(mt_tally$nmf_matrix)
  
  indel_nmf$Signature.norm
  Sig_tr <- indel_nmf_RNA$Signature.norm %*% rWishart(1, 10, diag(4))[,,1]
  colnames(Sig_tr) <- colnames(indel_nmf_RNA$Signature.norm)
  show_sig_profile(Sig_tr, mode = "ID", style = "cosmic", x_label_angle = 90)
  
  
mergedRNAmaf




