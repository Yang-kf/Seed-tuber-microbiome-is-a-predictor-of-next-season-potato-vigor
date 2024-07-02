#Citation
#Note: Field K is the same as field S (short for SPNA) in the script

library(randomForest)
require(caTools)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(gridExtra)


#load("my_workspace.RData")
#read function
RF_results_function <- readRDS("C:/Users/Y Song/OneDrive - Universiteit Utrecht/Postdoc/HZPC_project/Vitality_prediction_model/Blue_vitality_nooutlier/RF_results_function.rds")

######HEEL BAC######
{

#move to working direction
#load metadata
map = read.table("~/metadata_Merge_Heel_bac.txt",header = T, row.names = 1,sep = "\t")

setwd("~/Heel_bac") 
#input data at different phylogenetic level
#ASV
table_8 =read.table("~/Merge_Heel_bac_export_NoSruber_rar8000/Merge_Heel_bac_all_feature-table_NoSruber_rar8000.tsv",header = T, row.names = 1,sep = "\t")

# #Species
table_7 =read.table("~/Merge_Heel_bac_export_NoSruber_rar8000/Heel20182019_bac_feature_table_filtered2_NoNonsense_NoContamination_NoOutlier_NoSruber_rar8000_collapse_species.tsv",header = T, row.names = 1,sep = "\t")

# #Genus
table_6 =read.table("~/Merge_Heel_bac_export_NoSruber_rar8000/Heel20182019_bac_feature_table_filtered2_NoNonsense_NoContamination_NoOutlier_NoSruber_rar8000_collapse_genus.tsv",header = T, row.names = 1,sep = "\t")

# #Family
table_5 =read.table("~/Merge_Heel_bac_export_NoSruber_rar8000/Heel20182019_bac_feature_table_filtered2_NoNonsense_NoContamination_NoOutlier_NoSruber_rar8000_collapse_family.tsv",header = T, row.names = 1,sep = "\t")

# #Order
table_4 =read.table("~/Merge_Heel_bac_export_NoSruber_rar8000/Heel20182019_bac_feature_table_filtered2_NoNonsense_NoContamination_NoOutlier_NoSruber_rar8000_collapse_order.tsv",header = T, row.names = 1,sep = "\t")

# #Class
table_3 =read.table("~/Merge_Heel_bac_export_NoSruber_rar8000/Heel20182019_bac_feature_table_filtered2_NoNonsense_NoContamination_NoOutlier_NoSruber_rar8000_collapse_class.tsv",header = T, row.names = 1,sep = "\t")

# #Phylum
table_2 =read.table("~/Merge_Heel_bac_export_NoSruber_rar8000/Heel20182019_bac_feature_table_filtered2_NoNonsense_NoContamination_NoOutlier_NoSruber_rar8000_collapse_phylum.tsv",header = T, row.names = 1,sep = "\t")

# 
}
######HEEL FUN######
{

#move to working direction
#load metadata
map = read.table("~/metadata_Merge_Heel_fun.txt",header = T, row.names = 1,sep = "\t")
#ASV
table_8 =read.table("~/Merge_Heel_fun_export_rar4000/Merge_Heel_fun_all_feature-table_rar4000.tsv",header = T, row.names = 1,sep = "\t")
#Species
table_7 =read.table("~/Merge_Heel_fun_export_rar4000/Heel20182019_fun_feature_table_filtered2_NoNonsense_rar4000_collapse_species.tsv",header = T, row.names = 1,sep = "\t")
#Genus
table_6 =read.table("~/Merge_Heel_fun_export_rar4000/Heel20182019_fun_feature_table_filtered2_NoNonsense_rar4000_collapse_genus.tsv",header = T, row.names = 1,sep = "\t")
#Family
table_5 =read.table("~/Merge_Heel_fun_export_rar4000/Heel20182019_fun_feature_table_filtered2_NoNonsense_rar4000_collapse_family.tsv",header = T, row.names = 1,sep = "\t")
#Order
table_4 =read.table("~/Merge_Heel_fun_export_rar4000/Heel20182019_fun_feature_table_filtered2_NoNonsense_rar4000_collapse_order.tsv",header = T, row.names = 1,sep = "\t")
#Class
table_3 =read.table("~/Merge_Heel_fun_export_rar4000/Heel20182019_fun_feature_table_filtered2_NoNonsense_rar4000_collapse_class.tsv",header = T, row.names = 1,sep = "\t")
#Phylum
table_2 =read.table("~/Merge_Heel_fun_export_rar4000/Heel20182019_fun_feature_table_filtered2_NoNonsense_rar4000_collapse_phylum.tsv",header = T, row.names = 1,sep = "\t")
}
######HEEL BAC FUN Merge######
{
  
  map = read.table("~/metadata_Merge_Heel_bac.txt",header = T, row.names = 1,sep = "\t")
  #ASV
  table_8 =read.table("~/Merge_Heel_bac8000_fun4000_feature-table_ASV.txt",header = T, row.names = 1,sep = "\t")
  #Species
  table_7 =read.table("~/Merge_Heel_bac8000_fun4000_feature-table_species.txt",header = T, row.names = 1,sep = "\t")
  #Genus
  table_6 =read.table("~/Merge_Heel_bac8000_fun4000_feature-table_genus.txt",header = T, row.names = 1,sep = "\t")
  #Family
  table_5 =read.table("~/Merge_Heel_bac8000_fun4000_feature-table_family.txt",header = T, row.names = 1,sep = "\t")
  #Order
  table_4 =read.table("~/Merge_Heel_bac8000_fun4000_feature-table_order.txt",header = T, row.names = 1,sep = "\t")
  #Class
  table_3 =read.table("~/Merge_Heel_bac8000_fun4000_feature-table_class.txt",header = T, row.names = 1,sep = "\t")
  #Phylum
  table_2 =read.table("~/Merge_Heel_bac8000_fun4000_feature-table_phylum.txt",header = T, row.names = 1,sep = "\t")
}

# Define the list of input file names
phylevels <- list(table_8,table_7, table_6, table_5,table_4, table_3, table_2)
index <- list("ASV","species","genus","family","order", "class", "phylum")


corresults <- data.frame(label= character(), rho = numeric(),pvalue = numeric(), R2 = numeric())
for (i in seq_along(phylevels)) {
  # Get the current object
  otu_table <- phylevels[[i]]
  setwd("C:/Users/Y Song/OneDrive - Universiteit Utrecht/Postdoc/HZPC_project/Vitality_prediction_model/Blue_vitality_nooutlier/Heel_all") #change address
  index_value <- paste0("Heel_all_",index[[i]])#change dataset name here
  print(index_value)
  dir.create(paste0(index_value))
  setwd(index_value)
  
{
##Train with year2
train_map = map[map$Year %in% c("2019-2020"),]

# filter sample in otu_table
idx = rownames(train_map) %in% colnames(otu_table)
train_map = train_map[idx,]
train_otu = otu_table[, rownames(train_map)]  

X_train_otu <- t(train_otu)

Y_train_var <- as.factor(train_map$Variety)
Y_train_M2 <- train_map$M_52_method2
Y_test_S2 <- train_map$S_48_method2
Y_test_V2 <- train_map$V_50_method2


#use year1 as final validation (testing)
vali_map = map[map$Year %in% c("2018-2019"),]
idx = rownames(vali_map) %in% colnames(otu_table)
vali_map = vali_map[idx,]
vali_otu = otu_table[, rownames(vali_map)]

X_vali_otu <- t(vali_otu)

Y_vali_var <- as.factor(vali_map$Variety)
Y_vali_M2 <- vali_map$M_52_method2
Y_vali_S2 <- vali_map$S_48_method2
Y_vali_V2 <- vali_map$V_50_method2

#######Train#####

#make a base tree
X_train <- X_train_otu
Y_train <- Y_train_M2
set.seed(315)
rf = randomForest(X_train, Y_train, importance=TRUE, proximity=TRUE, ntree = 1000)

file_name <- paste("rf",index_value,".txt",sep = "")
sink(file_name)
rf
sink()

file_name <- paste("rf",index_value,".RData",sep = "")
save(rf, file = file_name)


set.seed(315)
result = rfcv(X_train, Y_train, cv.fold=5)
png("train_error.png")
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))
dev.off()

sink(paste0(index_value,"_train.error.txt"))
result$error.cv
sink()

# output feature importance
imp= as.data.frame(rf$importance)
imp = imp[order(imp[,1],decreasing = T),]
head(imp,n=10)
write.table(imp,file = paste0(index_value,"_importance_rf.txt"),quote = F,sep = '\t', row.names = T, col.names = T)

#######Test#######
#### same year different fields
 

###apply this function
#index_value <- "Heel_bac_species"

#M
train.p = predict( rf,X_train, type = "response")
df = data.frame(Observed = Y_train_M2, Predicted = train.p, Variety = train_map$Variety)
corresults1 <- RF_results_function (df,paste0(index_value,"_test_M2"),train_map)
corresults1
#S
#train.p = predict( rf,X_train, type = "response")
df = data.frame(Observed = Y_test_S2, Predicted = train.p)
corresults2 <-RF_results_function (df,paste0(index_value,"_test_S2"),train_map)
corresults2
#V
#train.p = predict( rf,X_train, type = "response")
df = data.frame(Observed = Y_test_V2, Predicted = train.p)
corresults3 <-RF_results_function(df,paste0(index_value,"_test_V2"),train_map)
corresults3

#with out-of-bag error
train.p = predict( rf, type = "response")
df = data.frame(Observed = Y_train_M2, Predicted = train.p)
corresults4 <-RF_results_function (df,paste0(index_value,"test_M2_oob"),train_map)
corresults4
#train.p = predict( rf, type = "response")
df = data.frame(Observed = Y_test_S2, Predicted = train.p)
corresults5 <-RF_results_function (df,paste0(index_value,"_test_S2_oob"),train_map)
corresults5
#train.p = predict( rf, type = "response")
df = data.frame(Observed = Y_test_V2, Predicted = train.p)
corresults6 <-RF_results_function(df,paste0(index_value,"_test_V2_oob"),train_map)


#######validation#####
###different year##
X_vali <- X_train_otu
Y_train <- Y_train_M2
#M
vali.p = predict( rf,X_vali_otu, type = "response")
df = data.frame(Observed = Y_vali_M2, Predicted = vali.p)
corresults7 <-RF_results_function (df,paste0(index_value,"_vali_M2"),vali_map)

#S
df = data.frame(Observed = Y_vali_S2, Predicted = vali.p)
corresults8 <-RF_results_function (df,paste0(index_value,"_vali_S2"),vali_map)

#V
df = data.frame(Observed = Y_vali_V2, Predicted = vali.p)
corresults9 <-RF_results_function (df,paste0(index_value,"_vali_V2"),vali_map)

# Initialize an empty dataframe to store the results

corresults <- rbind(corresults,corresults1,corresults2,corresults3,corresults4,corresults5,corresults6,corresults7,corresults8,corresults9)
corresults
 
save.image(file = paste0(index_value,"my_workspace.RData"))
  }

}

print(corresults)

dir.create("C:/Users/Y Song/OneDrive - Universiteit Utrecht/Postdoc/HZPC_project/Vitality_prediction_model/Blue_vitality_nooutlier/Heel_all/model_corresults") #change location
setwd("C:/Users/Y Song/OneDrive - Universiteit Utrecht/Postdoc/HZPC_project/Vitality_prediction_model/Blue_vitality_nooutlier/Heel_all/model_corresults") #change location
index_value = "Heel_all_" #change dataset name here
write.csv(corresults, file = paste0(index_value,"corresults.csv"))


{#arrange results
  # Define the vector of index values
  index_values <- c("test_M2", "test_S2", "test_V2")
  
  # Iterate over the index values and write the filtered data frames to files
  for (index_value2 in index_values) {
    # Filter the data frame
    df <- corresults %>% filter(grepl(index_value2, label) & !grepl("oob", label))
    
    # Write the filtered data frame to a CSV file
    write.csv(df, file = paste0(index_value, index_value2, "_corresults.csv"))
  }
  
  
  # Define the vector of index values
  index_values <- c("test_M2_oob", "test_S2_oob", "test_V2_oob","vali_M2","vali_S2","vali_V2")
  
  # Iterate over the index values and write the filtered data frames to files
  for (index_value2 in index_values) {
    # Filter the data frame
    df <- corresults %>% filter(grepl(index_value2, label))
    
    # Write the filtered data frame to a CSV file
    write.csv(df, file = paste0(index_value, index_value2, "_corresults.csv"))
  }
}
