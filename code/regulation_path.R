print("************************************* Generating regulation paths *************************************")
suppressMessages(library(dplyr))
init_path = getwd()
setwd(paste0(init_path, "/../processing"))

alter_pathogenic = read.table("../alter_pathogenic_region_variant_egene.csv", header = T, sep = ",", stringsAsFactors = F)
ref_pathogenic = read.table("../ref_pathogenic_region_variant_egene.csv", header = T, sep = ",", stringsAsFactors = F)
genes <- read.csv("cell_type_auc_mean.csv",header = T, sep = ",", stringsAsFactors = F)

max_celltype <- apply(genes[,2:6], 1, function(t) colnames(genes)[2:6][which.max(t)])
max_celltype <- as.data.frame(cbind(heart_egenes = genes$X, max_celltype))

if (nrow(alter_pathogenic)!=0){
    alter <- inner_join(alter_pathogenic, max_celltype, by = "heart_egenes")
    alter <- cbind(alter, disease = rep("cHF", nrow(alter)))
    alter <- alter[,c(1,2,3,5,6,4)]
    alter <- cbind(alter, type = rep("alternate-allele-pathogenic", nrow(alter)))
    colnames(alter)[5] <- "max_celltype"
}else alter <- NA
                      
if (nrow(ref_pathogenic)!=0){
    ref <- inner_join(ref_pathogenic, max_celltype, by = "heart_egenes")
    ref <- cbind(ref, disease = rep("cHF", nrow(ref)))
    ref <- ref[,c(1,2,3,5,6,4)]
    ref <- cbind(ref, type = rep("reference-allele-pathogenic", nrow(ref)))
    colnames(ref)[5] <- "max_celltype"
}else ref <- NA
                      
combine <- rbind(alter, ref)
combine <- na.omit(combine)

write.table(combine, "../regulation_path.csv", sep = ',', quote = F, row.names = F)
print('All Done! You may check the regulation paths in regulation_path.csv.')