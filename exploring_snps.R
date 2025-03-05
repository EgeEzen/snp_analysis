library(writexl)
library(readxl)
library(ggplot2)
library(envalysis)
library(GenomicRanges)
library(biomaRt)

setwd("~/Desktop/PhD/2024 Second Term/SNP analysis")

# RBPJ -> "chr4:26,078,378-26,129,770"
# ZFP36L1 -> "chr14:67,772,093-68,816,860"
# ZFP36L1 <- "chr14:67,772,093-67,939,102"
# 
# indicate genomic locus
string_from_ucsc <-"chr12:57,534,820-57,823,560"

# extract chromosome number and location ucsc genomic locus
string_from_ucsc <- gsub(",", "", string_from_ucsc)
# Split the string into chromosome and positions
string_from_ucsc <- strsplit(string_from_ucsc, "[:-]")[[1]]
# Extract the chromosome, start, and end positions
chrom <- string_from_ucsc[1]
start <- as.numeric(string_from_ucsc[2])
end <- as.numeric(string_from_ucsc[3])

command <- paste0("./plink --bfile GSA2022_964_025_V3_and_Ospelt1_2_filtered_hg38 --recode vcf --chr "
                  ,chrom," --from-bp ",start," --to-bp ",end, 
                  " --allow-extra-chr --out snp_outputs/test/result")
# run the plink
system(command)

command <- paste0("./plink --bfile GSA2022_964_025_V3_and_Ospelt1_2_filtered_hg38 --recode A --chr "
                  ,chrom," --from-bp ",start," --to-bp ",end, 
                  " --allow-extra-chr --out snp_outputs/test/result")

# run the plink
system(command)

# read the snp output
snp_output <- read.table(paste0("snp_outputs/test/result.raw"),header= TRUE)
# read the vcf output to get the coordinates
vcf_output <- read.table(paste0("snp_outputs/test/result.vcf"),header=FALSE)[,c(1:3)]

#to find the leading snp in the region read ra credible snps 
ra_cred_snps <- read.delim("RA_all_cred_hg38_PP.bed",header = FALSE)
colnames(ra_cred_snps) <- c("chromosome","start","end","posteriorprob")

#convert to granges and find the interested region
ra_cred_snps_gr <- GRanges(seqnames = ra_cred_snps$chromosome,
                           ranges = IRanges(start = ra_cred_snps$start, 
                                            end = ra_cred_snps$end),
                           posteriorprob = ra_cred_snps$posteriorprob)

snps_in_region <- subsetByOverlaps(ra_cred_snps_gr, GRanges(seqnames = chrom,
                                                            ranges = IRanges(start = start, end = end)))
snps_in_region <- as.data.frame(snps_in_region)

#now put the names of these snps to the snps in the region
snps_in_region <- merge(snps_in_region, vcf_output, all.x = TRUE, by.x = "start",by.y = "V2")
snps_in_region <- snps_in_region[,c(2,1,3,6,8)]
snps_in_region <- snps_in_region[!is.na(snps_in_region$V3),] #not all snps are genotyped...

#filter snp_output so that only snps of interests are there 
snp_ids_clean <- gsub("_[ATCG]", "", colnames(snp_output)[7:ncol(snp_output)])
colnames(snp_output)[7:ncol(snp_output)] <- snp_ids_clean
common_snps <- intersect(snps_in_region$V3, colnames(snp_output))
filtered_snp_output <- snp_output[, c("FID", common_snps), drop = FALSE]

#to check if they are in linkage disequlibrium
same_value_per_row <- apply(filtered_snp_output[, -1], 1, function(row) length(unique(row)) == 1)
table(same_value_per_row) 
rows_with_identical_snps <- filtered_snp_output[same_value_per_row, ]
head(rows_with_identical_snps)
rows_with_different_snps <-  filtered_snp_output[!same_value_per_row, ]
head(rows_with_different_snps)

#to see it in a manhattan plot
library(ggplot2)
library(ggrepel)

top_snps <- snps_in_region[order(snps_in_region$posteriorprob, decreasing = TRUE)[1:4], ]

# Manhattan plot for the SNPs
ggplot(snps_in_region, aes(x = start, y = posteriorprob)) +
  geom_point(color = "#2DAA9E90", size = 2) +  # SNP points
  geom_label_repel(data = top_snps, aes(label = V3), fill="white",alpha=0.7) +  # Labels for top 4 SNPs
  theme_publish() +
  labs(x = "Genomic Position (chr4)", 
       y = "Posterior Probability", 
       title = "AGAP2-AS1 Region SNPs"
      ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size=12))  # Rotate x-axis labels if needed
  
ggsave(file = "posterior_probability_agap2-as1.png",
       plot = last_plot(),
       dpi = 300, 
       bg = "transparent", 
       width = 6, height = 4
)
