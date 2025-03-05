library(writexl)
library(readxl)
library(ggplot2)
library(envalysis)
library(dplyr)
library(tidyr)

setwd("~/Desktop/PhD/2024 Second Term/SNP analysis")

# create snp file 
create_snp_file <- function(snp){

  command <- paste0("./plink --bfile GSA2022_964_025_V3_and_Ospelt1_2_filtered_hg38 --recode A --snps ",snp, 
                    " --allow-extra-chr --out snp_outputs/",snp,"/result")
  # run the plink
  system(paste0("mkdir snp_outputs/",snp)) #make the directory
  system(command)
  
  # for the allelic information 
  system(paste0("grep ", snp ," GSA2022_964_025_V3_and_Ospelt1_2_filtered_hg38.bim > snp_outputs/",snp,"/allel.txt"))
  
  # read the snp output
  snp_output <- read.table(paste0("snp_outputs/",snp,"/result.raw"),header= TRUE)
  
  # organize data
  snp_output <- snp_output[,c(1,7)]
  
  # read the snp output additional information (regarding snp bases)
  snp_output_text <- read.table(paste0("snp_outputs/",snp,"/allel.txt"),tryLogical = FALSE)
  
  # convert numbers to bases (show allele)
  snp_output$base <- ifelse(snp_output[,2] == "0",paste0( snp_output_text[,6], snp_output_text[,6]),
                            ifelse(snp_output[,2] == "1",paste0( snp_output_text[,5], snp_output_text[,6]), 
                                   paste0( snp_output_text[,5], snp_output_text[,5])))
  
  head(snp_output)
  
  #read metadata and merge with the snp outputs
  metadata <- read_xlsx("metadata_final_.xlsx")
  
  #prepare the data
  snp_output <- merge(snp_output, metadata, by = "row.names")
  snp_output <- snp_output[,c(5,3,4,6:8)]
  colnames(snp_output) <- c("ID", colnames(snp_output)[2:6])
  
  #write it
  write_xlsx(snp_output, paste0("snp_outputs/",snp,"_table.xlsx"))
  
  return(snp_output)
}

# doing eqtl within ra using all
eqtl_in_ra <- function(snp_output, use_norm_counts, gene_of_interest){
  
  if (!use_norm_counts){
    #for eQTL on RA - read the count data
    rna_seq_count_data <- read.csv("~/Desktop/PhD/2024 Second Term/RNA-seq all/count_data_all.csv",row.names = 1)
    #do CPM for eQTL
    lib_size <- colSums(rna_seq_count_data)
    crna_seq_count_cpm <- t( t(rna_seq_count_data) / lib_size * 1e6 )
    gene_of_interest_cpm <- crna_seq_count_cpm[gene_of_interest,] 
  } else{
    # !!! or read the normalized and batch effect removed count data? !!! #
    rna_seq_count_data <- read.csv("~/Desktop/PhD/2024 Second Term/RNA-seq all/vst_counts_batch_effect_removed.csv",row.names = 1)
    gene_of_interest_cpm <- rna_seq_count_data[gene_of_interest,]
    named_vector <- as.numeric(gene_of_interest_cpm)
    names(named_vector) <- colnames(gene_of_interest_cpm)
    gene_of_interest_cpm <- named_vector
  }

  #get matchings from the snp_output
  snp_output_bulk_ra <- snp_output[snp_output$single_cell == "bulk" & snp_output$Diagnosis == "RA",]
  snp_output_bulk_ra$ID <- sub("^([0-9]+)_RA", "RA_\\1", snp_output_bulk_ra$ID) #to match the cpm colnames
  matching_indices <- match(snp_output_bulk_ra$ID, names(gene_of_interest_cpm))
  
  #subset cpm values based on the matching samples
  matched_cpm <- gene_of_interest_cpm[matching_indices]
  snp_output_bulk_ra$CPM <- matched_cpm
  
  #for the plot 
  snp_plot <-snp_output_bulk_ra[!is.na(snp_output_bulk_ra$CPM),]
  
  # Convert genotype to factor
  genotype_order <- unique(snp_plot[,2])
  genotype_base <- unique(snp_plot[,3])
  genotype_list <- setNames(genotype_base, genotype_order)
  genotype_list <- genotype_list[order(as.numeric(names(genotype_list)))] #to make sure that they follow 0,1,2
  snp_plot$base <- factor(snp_plot$base, levels = genotype_list)
  snp_plot$numeric <- as.character(as.numeric(snp_plot$base)-1)
  
  #linear regression
  #lm_eqtl <- lm(CPM ~ numeric + Gender, data = snp_plot)
  #summary(lm_eqtl)
  
  # Create a violin+box+point plot
  p <- ggplot(snp_plot, aes(x = numeric, y = CPM, fill = base)) +
    geom_violin(alpha = 0.3, position = "dodge") +
    geom_boxplot(width = 0.2, position = position_dodge(0.9), alpha = 0.7,show.legend = FALSE) + 
    geom_jitter(width = 0.2, alpha = 0.7,show.legend = FALSE) + 
    labs(title = paste0("eQTL on ",gene_of_interest, " in RA\n","SNP number:",snp), x = "Genotype", y = "Gene Expression") +
    theme_publish() +
    theme(legend.position = "right")
  
  if (!use_norm_counts){
    file = paste0("figures_test_3/",gene_of_interest,"_",snp, "cpm_eQTL.png")
  }else{
    file = paste0("figures_test_3/",gene_of_interest,"_",snp, "BE_removed_nrm_eQTL.png")
  }
  ggsave(file = file,
    plot = p,
    dpi = 300, 
    bg = "transparent", 
    width = 4, height = 4
  )
  return(p)
}

# for getting the tnf vs nontnf columns 
get_tnf_cols <- function(df){
  
  ra_cols <- colnames(df)[grep("^RA_", colnames(df))]
  ra_numbers <- unique(gsub("_TNF", "", ra_cols))
  matched_cols <- unlist(lapply(ra_numbers, function(x) {
    pair <- ra_cols[grep(paste0("^", x, "(_TNF)?$"), ra_cols)]
    if (length(pair) == 2) return(pair)  # Keep only if both exist
  }))
  df_updated <- df[, matched_cols]
  return(df_updated)
}

# doing eqtl in ra tnf vs nontnf 
eqtl_in_ra_tnf_vs_non_tnf <- function(snp_output, use_norm_counts, gene_of_interest){
  
  if (!use_norm_counts){
    #for eQTL on RA - read the count data
    rna_seq_count_data <- read.csv("~/Desktop/PhD/2024 Second Term/RNA-seq all/count_data_all.csv",row.names = 1)
    #do CPM for eQTL
    lib_size <- colSums(rna_seq_count_data)
    crna_seq_count_cpm <- t( t(rna_seq_count_data) / lib_size * 1e6 )
    crna_seq_count_cpm_tnf <-get_tnf_cols(crna_seq_count_cpm)
    gene_of_interest_cpm <- crna_seq_count_cpm_tnf[gene_of_interest,] 
  } else{
    # !!! or read the normalized and batch effect removed count data? !!! #
    rna_seq_count_data <- read.csv("~/Desktop/PhD/2024 Second Term/RNA-seq all/vst_counts_batch_effect_removed.csv",row.names = 1)
    rna_seq_count_data_tnf <- get_tnf_cols(rna_seq_count_data)
    gene_of_interest_cpm <- rna_seq_count_data_tnf[gene_of_interest,]
    named_vector <- as.numeric(gene_of_interest_cpm)
    names(named_vector) <- colnames(gene_of_interest_cpm)
    gene_of_interest_cpm <- named_vector
  }
  
  #get matchings from the snp_output
  snp_output_bulk_ra <- snp_output[snp_output$single_cell == "bulk" & snp_output$Diagnosis == "RA",]
  snp_output_bulk_ra$ID <- sub("^([0-9]+)_RA", "RA_\\1", snp_output_bulk_ra$ID) #to match the cpm colnames
  matching_indices <- match(snp_output_bulk_ra$ID, names(gene_of_interest_cpm))
  
  #subset cpm values based on the matching samples
  matched_cpm <- gene_of_interest_cpm[matching_indices]
  snp_output_bulk_ra$CPM <- matched_cpm
  #tnf will be always +1 index
  matching_indices[!is.na(matching_indices)] <- matching_indices[!is.na(matching_indices)] + 1
  snp_output_bulk_ra$CPM_TNF <- gene_of_interest_cpm[matching_indices]
  
  #for the plot 
  snp_plot <-snp_output_bulk_ra[!is.na(snp_output_bulk_ra$CPM),]

  # Convert genotype to factor
  genotype_order <- unique(snp_plot[,2])
  genotype_base <- unique(snp_plot[,3])
  genotype_list <- setNames(genotype_base, genotype_order)
  genotype_list <- genotype_list[order(as.numeric(names(genotype_list)))] #to make sure that they follow 0,1,2
  snp_plot$base <- factor(snp_plot$base, levels = genotype_list)
  snp_plot$numeric <- as.character(as.numeric(snp_plot$base)-1)
  
  #linear regression
  #lm_eqtl <- lm(CPM ~ numeric + Gender, data = snp_plot)
  #summary(lm_eqtl)
  
  #for plotting
  snp_plot <- snp_plot %>%
    pivot_longer(cols = c(CPM, CPM_TNF), names_to = "Condition", values_to = "Expression") %>%
    mutate(Condition = factor(ifelse(Condition == "CPM", "Non-TNF", "TNF"), levels = c("Non-TNF", "TNF"))) 
  
  # Create plot
  p_tnf <- ggplot(snp_plot, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_violin(alpha = 0.3, position = "dodge") +  # Violin plot for distribution
    geom_boxplot(width = 0.2, position = position_dodge(0.9), alpha = 0.7, show.legend = FALSE) +  # Boxplot
    geom_jitter(width = 0.2, alpha = 0.7, show.legend = FALSE) +  # Jitter for individual points
    geom_line(aes(group = ID), color = "black", alpha = 0.5) +  # Connect same samples
    facet_wrap(~ numeric) +  # Facet by Genotype
    labs(title = paste0("eQTL on ", gene_of_interest, " in RA\n", "SNP number: ", snp),
         x = "Condition", y = "Gene Expression") +
    theme_publish() +
    theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))

  # save the file
  if (!use_norm_counts){
    file = paste0("figures_test_3/",gene_of_interest,"_",snp, "cpm_eQTL_tnf.png")
  }else{
    file = paste0("figures_test_3/",gene_of_interest,"_",snp, "BE_removed_nrm_eQTL_tnf.png")
  }
  ggsave(file = file,
         plot = p_tnf,
         dpi = 300, 
         bg = "transparent", 
         width = 5, height = 4
  )
  return(p_tnf)
}

snp_list <-c("rs3128921")

for (snp in snp_list) {

  use_norm_counts <- FALSE
  gene_of_interest <- "HLA-DRB1"
  
  snp_output <- create_snp_file(snp)
  plot <- eqtl_in_ra(snp_output, use_norm_counts, gene_of_interest)
  tnf_plot <- eqtl_in_ra_tnf_vs_non_tnf(snp_output, use_norm_counts, gene_of_interest)
  
}
