print("********************************* Linking causal region, causal variant and causal eGene *********************************")
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
init_path = getwd()
setwd(init_path)

args <- commandArgs(trailingOnly=TRUE)
p.thresh <- as.numeric(args[1])

alter_name = "/../processing/region_alter_pathogenic.csv"
ref_name = "/../processing/region_reference_pathogenic.csv"
snps_name = "/../demo_data/GWAS.txt"
eqtls_name = "/../demo_data/eqtls_200kb.txt"

greater <- read.csv(paste0(init_path, alter_name), header = T, stringsAsFactors = F)
less <- read.csv(paste0(init_path, ref_name), header = T, stringsAsFactors = F)
fdr = p.adjust(greater$p_fisher, method = 'bonferroni')
greater <- cbind(greater, fdr)
greater_sig <- greater[greater$fdr<p.thresh,]
fdr = p.adjust(less$p_fisher, method = 'bonferroni')
less <- cbind(less, fdr)
less_sig <- less[less$fdr<p.thresh,]
greater_ssig <- greater_sig[!(greater_sig$snp %in% less_sig$snp), ]
less_ssig <- less_sig[!(less_sig$snp %in% greater_sig$snp), ]
snp_pos <- read.table(paste0(init_path, snps_name), header = T, stringsAsFactors = F, sep = '\t')
colnames(snp_pos)[1] <- "snp"
colnames(snp_pos)[2] <- "chromosome"
greater_ssig <- left_join(greater_ssig, snp_pos, by = "snp")
less_ssig <- left_join(less_ssig, snp_pos, by = "snp")
less_ssig_sep <- separate_rows(less_ssig, causal_variants, sep = " ")[,1:4]
greater_ssig_sep <- separate_rows(greater_ssig, causal_variants, sep = " ")[,1:4]

eqtls <- read.table(paste0(init_path, eqtls_name), header = T, stringsAsFactors = F, sep = '\t')
colnames(eqtls)[2] <- "causal_variants"
alter_pathogenic <- inner_join(greater_ssig_sep, eqtls, by = c("snp", "causal_variants"))
alter_pathogenic <- alter_pathogenic[, -2]
alter_pathogenic <- alter_pathogenic[,c(1, 2, 4, 3)]
ref_pathogenic <- inner_join(less_ssig_sep, eqtls, by = c("snp", "causal_variants"))
ref_pathogenic <- ref_pathogenic[, -2]
ref_pathogenic <- ref_pathogenic[,c(1, 2, 4, 3)]
write.csv(alter_pathogenic, paste0(init_path, "/../alter_pathogenic_region_variant_egene.csv"), quote=F, row.names=F)
write.csv(ref_pathogenic, paste0(init_path, "/../ref_pathogenic_region_variant_egene.csv"), quote=F, row.names=F)
print("All done!")

