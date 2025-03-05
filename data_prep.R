library(writexl)
library(readxl)

setwd("~/Desktop/PhD/2024 Second Term/SNP analysis")

#read the manually made metadata
metadata <- read_xlsx("snp_metadata.xlsx")

#for consistency change the names of OA and RA samples
metadata$FID[c(1:69,128:132)] <- sub("^(RA|OA)(\\d{3})$", "\\2_\\1", metadata$FID[c(1:69,128:132)])

#for single cell get the synbio
synbio <- read_xlsx("SynBio_samples_DNA_summary.xlsx")

#gender and diagnosis information
metadata_copy <- merge(metadata, synbio[,c(1,3,4)] , all.x = TRUE, by.x = "FID", by.y = "Samples_sequenced",sort =FALSE)

#give info if single cell or not
metadata_copy$single_cell <- ifelse(is.na(metadata_copy$Gender),"bulk", "single_cell")

#for ra and oa get gender and diagnosis
ra_and_oa <- read_excel("~/Desktop/PhD/SHK_Schulthess Klinik_20241009.xlsx")[,c(2,5,8)]

#(e.g., "316_RA" → "316") or SHK 120 to 120
metadata_copy$Nr. <- sub("_(RA|OA)$", "", metadata_copy$FID)
ra_and_oa$Nr. <- sub("^SHK ", "", ra_and_oa$Nr.)
#merge them
merged_df <- merge(metadata_copy[46:132,c("FID","single_cell", "Nr.")], ra_and_oa, by = "Nr.", all.x =TRUE)
merged_df$Nr. <- NULL
merged_df <- merged_df[,c(1,4,3,2)]
colnames(merged_df) <- c("FID","Diagnosis","Gender","single_cell")
metadata_copy$Nr. <- NULL
metadata_copy <- rbind(metadata_copy[1:45,],merged_df)

#get the order back, important!
metadata_final <- metadata_copy[match(metadata$FID, metadata_copy$FID), ]

#rest will be changed manually
write_xlsx(metadata_final,"metadata_final.xlsx")


# --- to find the LD related risk loci -----

il6_ld <- read.delim("rs1800795.tsv")
