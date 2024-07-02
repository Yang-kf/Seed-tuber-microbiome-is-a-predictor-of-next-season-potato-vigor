
#make a confusion matrix
#visualize it
library(dplyr)
library(caret)
library(gridExtra)
library(ggplot2)
#df is the output from RF analysis observed and predicted
setwd("C:/Users/Y Song/OneDrive - Universiteit Utrecht/Postdoc/HZPC_project/Vitality_prediction_model/Blue_vitality_nooutlier/Eye_all/Eye_all_ASV")
map = read.table("C:/Users/Y Song/OneDrive - Universiteit Utrecht/Postdoc/HZPC_project/Merge_Sample_2018_2019_rerun/metadata_with_BLUE/metadata_Merge_Eye_bac_precrop_correctSoiltype_BLUE_NoOutlier_X.txt",header = T, row.names = 1,sep = "\t")

#possible to use a loop here but I did it one by one
id = "Eye_all_ASV_train_M2_oob" #limits = c(-1.05, 1) c(-3, 2.5)
id = "Eye_all_ASV_test_M2"  #limits = c(-1.8, 1.8) c(-3, 2.5)
id = "Eye_all_ASV_test_S2"  #limits = c(-1.8, 1.8) c(-3, 2.5)
id = "Eye_all_ASV_test_V2"  #limits = c(-1.8, 1.8) c(-3, 2.8)
id = "Eye_all_ASV_vali_M2"  #limits = c(-0.7, 0.7) c(-2.5, 2.5)
id = "Eye_all_ASV_vali_S2"   #limits = c(-0.7, 0.7) c(-3, 2.5)
id = "Eye_all_ASV_vali_V2"  #limits = c(-0.7, 0.7) c(-3, 2.1)
{
filename = paste0("C:/Users/Y Song/OneDrive - Universiteit Utrecht/Postdoc/HZPC_project/Vitality_prediction_model/Blue_vitality_nooutlier/Eye_all/Eye_all_ASV/", id, ".txt")
df = read.table(file = filename,header = T, sep = "\t", row.names = 1)

# df = read.table(file = "C:/Postdoc/HZPC_project/Vitality_prediction_model/Blue_vitality_nooutlier/Eye_bac/Eye_bac_ASV/Eye_bac_ASV_test_M2.txt",header = T, sep = "\t", row.names = 1)
# df = read.table(file = "C:/Postdoc/HZPC_project/Vitality_prediction_model/Blue_vitality_nooutlier/Eye_bac/Eye_bac_ASV/Eye_bac_ASV_vali_M2.txt",header = T, sep = "\t", row.names = 1)
# df = read.table(file = "C:/Postdoc/HZPC_project/Vitality_prediction_model/Blue_vitality_nooutlier/Eye_bac/Eye_bac_ASV/Eye_bac_ASV_vali_S2.txt",header = T, sep = "\t", row.names = 1)

#define the function to classify observed and predicted into High, Low and Middle
bayesQuantile <- function(yvec) {
  n <- length(yvec)
  # labelL <- 0
  # labelM <- 1
  # labelH <- 2
  labelL <- "L"
  labelM <- "M"
  labelH <- "H"
  labels <- rep(0, n)
  yvecB1 <- quantile(yvec, probs = 1/3) # bound 1
  yvecB2 <- quantile(yvec, probs = 2/3) # bound 2
  for (ii in 1:n) {
    if (yvec[ii] <= yvecB1) {
      labels[ii] <- labelL
    } else if (yvec[ii] > yvecB1 & yvec[ii] <= yvecB2) {
      labels[ii] <- labelM
    } else if (yvec[ii] > yvecB2) {
      labels[ii] <- labelH
    }
  }
  return(labels) # return vector directly
}

#apply this function to df and generate df_class
Observed_class <- bayesQuantile(df$Observed) # extract labels vector from returned list
Predicted_class <- bayesQuantile(df$Predicted)

#generate cm for all varieties
conf_mat <- confusionMatrix(factor(Observed_class), factor(Predicted_class))
sink(paste0("confusionMatrix_statistics_",id,"_all.txt"))
conf_mat
conf_mat$byClass
sink()

#function for plotting confusion matrix
plot_confusion_matrix <- function(cm_class,var) {
  
  # calculate proportions of the confusion matrix
  cm_prop <- round(prop.table(cm_class, margin = 1), 2)
  
  # melt the confusion matrix into a long data frame
  cm_melt <- reshape2::melt(cm_prop)
  
  # plot the confusion matrix using ggplot2
  p <- ggplot(data = cm_melt,
              aes(x = Predicted_class, y = Observed_class, fill = value)) +
    geom_tile(color = "black") +
    geom_text(aes(label = ifelse(value > 0.8, "", value)), size = 9,
              color = ifelse(cm_melt$value > 0.8, "white", "black")) +  # Set font color based on condition
    #scale_fill_gradient2(low = "white", mid = "white", high = "steelblue", 
                         #limits = c(0, 1), midpoint = 0.15, space = "Lab") +
    scale_fill_gradientn(colors = c("white", "white", "steelblue","#152f52","#152f52"),
                         values = c(0, 0.15, 0.6,0.85,1),limits = c(0, 1), space = "Lab") +
    theme_minimal() +
    labs(x = "Predicted", y = "Observed", title = paste0("Variety ",var)) +
    theme(axis.text.x = element_text(size =18))+
    theme(axis.text.y = element_text(size =18))+ 
    theme(plot.title = element_text(hjust = 0.5, size = 20))+
    scale_x_discrete(limits = c("L", "M", "H"))+ 
    scale_y_discrete(limits = c("L", "M", "H"))
  
  p 
  # return the plot object
  return(p)
  
}

# Function for plotting confusion matrix with percentages
plot_confusion_matrix_percentage <- function(cm_class, var) {
  
  # Calculate row-wise percentages of the confusion matrix
  cm_prop_row <- round(prop.table(cm_class, margin = 1) * 100, 0)  # Change: Calculate row-wise percentages and multiply by 100
  
  # Melt the confusion matrix into a long data frame
  cm_melt <- reshape2::melt(cm_prop_row)
  
  # Plot the confusion matrix using ggplot2
  p <- ggplot(data = cm_melt,
              aes(x = Predicted_class, y = Observed_class, fill = value)) +
    geom_tile(color = "black") +
    geom_text(aes(label = paste0(value, "%")), size = 9,  # Change: Include "%" symbol in label
              color =  "black") +  # Set font color to black
   
    scale_fill_gradientn(colors = c("white", "white", "steelblue","#29416F","#29416F"),
                         values = c(0, 0.2, 0.6,0.9,1),limits = c(0, 100), space = "Lab") + # Gradient scale for colors
    theme_minimal() +
    labs(x = "Predicted", y = "Observed", title = paste0("Variety ", var)) +
    theme(axis.text.x = element_text(size = 18)) +
    theme(axis.text.y = element_text(size = 18)) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20)) +
    scale_x_discrete(limits = c("L", "M", "H")) + 
    scale_y_discrete(limits = c("L", "M", "H"))
  
  return(p)
  
}


#use cm_class as input for function "plot_confusion_matrix"
cm_class <- table(Observed_class,Predicted_class)
pdf(paste0("confusionMatrix_percentage_",id,"_all.pdf"),height = 4,width = 5)
p_all <- plot_confusion_matrix_percentage(cm_class,"All")
# pdf(paste0("confusionMatrix",id,"_all.pdf"),height = 4,width = 5)
# p_all <- plot_confusion_matrix(cm_class,"All Variables")
p_all
dev.off()

#add vitality data by merge df_class with metadata
df_class <- cbind(df, Observed_class, Predicted_class) # bind vectors to df
# Merge the two data frames based on their row names
merged_df <-merge(df, map, by=0, all.x=TRUE)
#save df_class_var
# Extract columns by their names and combine into a new dataframe
df_class_var <- merged_df[, c("Observed", "Predicted", "VarietyFullName", "Variety")]
row.names(df_class_var) <- merged_df$Row.names
write.table(df_class_var, file = paste0("Observed_predicted_class_",id,"_all.txt"))


# Create a vector with the variables
vars <- c("A", "B", "C", "D", "E", "F")
# Create an empty list to store the plots
plots_list <- list()
#Initialize a list to store confusion matrices for each genotype
var_conf_mats <- list()
# Loop over the variables
for(var in vars) {
  # Subset the data for the current variable
  idx <- which(df_class_var$Variety == var)
  subset_y <- df_class[idx,]
  
  # Generate the confusion matrix for the current variable
  subset_conf_mat <- table(Observed_class = subset_y$Observed_class, Predicted_class = subset_y$Predicted_class)
 
  # # Create the plot for the current variable
  # p <- plot_confusion_matrix(subset_conf_mat, var)
  # Create the plot for the current variable
  p <- plot_confusion_matrix_percentage(subset_conf_mat, var)
  
  # Add the plot to the list
  plots_list[[var]] <- p
  
  # Store the confusion matrix statistics for this genotype in the list
  subset_conf_mat_sta <- confusionMatrix(factor(subset_y$Observed_class), factor(subset_y$Predicted_class))
  var_conf_mats[[var]] <-  subset_conf_mat_sta
}


# Add the combined plot to the beginning of the list
plots_list <- c(list("All" = p_all), plots_list)

# Save all the plots in a single PDF file
# Function to extract legend from a ggplot object
get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Create a separate legend for the plots
legend <- get_legend(plots_list[[var]])

# Remove legend from all plots
plots_list_nolegend <- lapply(plots_list, function(p) {
  p + theme(legend.position = "none")
})
# Create a new PDF file
# filename <- paste0(id,"_confusion_matrix_percentage_per_variety.pdf")
# #filename <- paste0(id,"_confusion_matrix_per_variety.pdf")
# pdf(filename,width = 26, height = 4)
# 
# # Arrange the plots and the legend together
# grid.arrange(arrangeGrob(grobs = plots_list_nolegend, ncol = 6), legend, widths = c(8, 1))
# 
# # Close the PDF file
# dev.off()

# Create a new PDF file
filename <- paste0(id,"_confusion_matrix_percentage_per_variety&all.pdf")
#filename <- paste0(id,"_confusion_matrix_per_variety&all.pdf")
pdf(filename,width = 30, height = 3.5)
#including all plot
grid.arrange(arrangeGrob(grobs = plots_list_nolegend, ncol = 7), legend, widths = c(10, 1))

# Close the PDF file
dev.off()


# Print the confusion matrices for each genotype
sink(file = paste0("C:/Users/Y Song/OneDrive - Universiteit Utrecht/Postdoc/HZPC_project/Vitality_prediction_model/Blue_vitality_nooutlier/Eye_all/Eye_all_ASV/",id,"_confusion-matrix_percentage_statistics_per_variety.txt"))
for (var in vars) {
  cat(paste0("Confusion matrix for variety ", var, ":\n"))
  print(var_conf_mats[[var]])
  print(var_conf_mats[[var]]$byClass)
}
sink()

}
############################################################
#Adapted from RF_results_function
#arrange correlation plots
{
#plot corrections with updated format
p=ggplot(data = df, mapping = aes(x=Predicted,y=Observed)) +
geom_point()
p
# 预测结果评估
result = matrix(0, ncol = 7, nrow = 1)
cor = cor.test(df[,1], df[,2], method = "spearman")
df2 = df
colnames(df2) = c("x", "y")
m = lm(y ~ x, df2)
#add a new column "All" so can use facet to have the same format as the rest of varieties
df_class_var_All <- df_class_var
df_class_var_All$All <- "All"

p1 = ggplot(df_class_var_All, aes(Predicted, Observed)) +
#geom_point(size =4,aes(color=Variety)) +
geom_point(size =4,color="grey40",alpha = 0.6) +
geom_smooth(method = "lm") +
scale_color_brewer(type = "seq",palette = "Set2")+
labs(title = paste("rho = " , round(cor$estimate, digits = 3), ", P = " , signif(cor$p.value, digits = 3), ", R2 = ", round(summary(m)$r.squared, digits = 3) , sep = "")) +
theme_bw() +
theme(text = element_text(size = 17),plot.title = element_text(size = 17)) +
theme(legend.position="none") +
scale_x_continuous(limits = c(-2, 2), breaks = seq(-1, 1, by = 1)) + #labels = scales::number_format(accuracy = 1)) +
scale_y_continuous(limits = c(-3, 3), breaks = seq(-2, 2, by = 2))#labels = scales::number_format(accuracy = 1))
p1

p3 = ggplot(df_class_var_All, aes(Predicted, Observed)) +
geom_point(size =3,aes(color=Variety),alpha = 0.8) +
#geom_point(size =4,color="grey40",alpha = 0.6) +
geom_smooth(method = "lm") +
scale_color_brewer(type = "seq",palette = "Set2")+
labs(title = paste("rho = " , round(cor$estimate, digits = 3), ", P = " , signif(cor$p.value, digits = 3), ", R2 = ", round(summary(m)$r.squared, digits = 3) , sep = "")) +
theme_bw() +
theme(text = element_text(size = 17),plot.title = element_text(size = 17)) +
theme(legend.position="none")+
scale_x_continuous(limits = c(-0.7, 0.7), breaks = seq(-0.5, 0.5, by = 0.5)) + #labels = scales::number_format(accuracy = 1)) +
scale_y_continuous(limits = c(-3, 2.1), breaks = seq(-2, 2, by = 2))#labels = scales::number_format(accuracy = 1))
p3
(p2 = p3+facet_wrap(~All,ncol=6))

#make a new factor with whole variety names

VarietyFull <- c("A" = "Challenger",
                 "B" = "Colomba",
                 "C" = "Festien",
                 "D" = "Innovator",
                 "E" = "Sagitta",
                 "F" = "Seresta")
p4 = p3+facet_wrap(~Variety,ncol=6,labeller = labeller(Variety = VarietyFull))
p4
#save figures
# # Create the file name using the index value
file_name <- paste(id,"_correlation_color.pdf",sep = "")
pdf(file_name,width=18,height=2.8)
grid.arrange(p2, p4, ncol=2,widths=c(2.1, 10))
grid.arrange(p2, p4, ncol=3,widths=c(2.1,10,3.5))

dev.off()

# library(tidyverse)
# try to add col to all variety
cor <- df_class_var %>% group_by(df_class_var$Variety) %>% summarize(spearman_cor=cor(Predicted, Observed,method="spearman"),p = cor.test(Predicted, Observed, method = "spearman")$p.value,R2 = summary(lm(Observed ~ Predicted))$r.squared)
file_name <- paste(id,"_cor_variety_color.txt",sep = "")
write.table(cor,file = file_name,quote = F,sep = '\t', row.names = T, col.names = T)
}













#STOP HERE######

# 
# Start a new PDF device to save the plots
pdf("all_plots.pdf")

# Loop over the unique genotypes and generate a confusion matrix for each one
for (var in varieties) {
  # Get the subset of samples that belong to this genotype
  idx <- which(df_class_var$Variety == var)
  #idx <- which(df_class_var$Variety == "A")
  subset_y <- df_class[idx,]
  ##STOP here, need to find a step to use plot_confusion_matrix function
  #subset_y_class <- as.data.frame(cbind(subset_y$Observed_class, subset_y$Predicted_class))
  
  # Generate the confusion matrix for this genotype
  subset_conf_mat <- table(Observed_class = subset_y$Observed_class, Predicted_class = subset_y$Predicted_class)
  print(var)
  p <- plot_confusion_matrix(subset_conf_mat,var)
  pdf(p,file = "C:/Postdoc/HZPC_project/Vitality_prediction_model/Blue_vitality_nooutlier/Eye_bac/Eye_bac_ASV/Eye_bac_ASV_test_M2_confusion-matrix.pdf")
  # Store the confusion matrix for this genotype in the list
  #var_conf_mats[[var]] <-  subset_conf_mat
  
}
# Close the PDF device
dev.off()








plot_confusion_matrix(subset_conf_mat)

# p <- ggplot(data = cm_melt,
#             aes(x = Var1, y = Var2, fill = value)) +
#   geom_tile() +
#   geom_text(aes(label = value), size = 12, color = "black") +
#   scale_fill_gradient(low = "white", high = "steelblue") +
#   theme_minimal() +
#   labs(x = "Predicted", y = "Observed", title = "Confusion Matrix") +
#   theme(axis.text.x = element_text(size =15))+
#   theme(axis.text.y = element_text(size =15))
# p
# Create a confusion matrix
library(ggplot2)
cm_class <- table(Observed_class,Predicted_class)#,merged_df$Variety)
plot_confusion_matrix(cm_class,"all")
plot_confusion_matrix(subset_conf_mat)

plot_confusion_matrix(subset_y_class)
subset_conf_mat
cm_class
a = as.data.frame(subset_conf_mat) ##somehow the formart is wrong as plot function input
a
plot_confusion_matrix(a)

plot_confusion_matrix(subset_y_class)# load the caret package
install.packages("caret")
library(caret)

install.packages("DMwR")
# library("DMwR")
# # install DMwR package if not already installed
# if (!requireNamespace("DMwR", quietly = TRUE)) {
#   install.packages("DMwR")
# }
# 
# # load the DMwR package
# library(DMwR)
# create a confusion matrix

cm <- confusionMatrix(factor(Predicted_class, levels = c(0, 1, 2)), 
                      factor(Predicted_class, levels = c(0, 1, 2)))
cm
# calculate the F-score
# f_score <- fMeasure(cm)
# f_score <- F1(cm)
f_score <- cm$byClass["F1"]
print(f_score)


library(ggplot2)
library(scales)

# Create a matrix of relative frequencies
cm_matrix <- cm$table
cf_matrix_prop <- as.data.frame.matrix(cm_matrix/sum(cm_matrix))

# Convert the matrix to a data frame for ggplot2
df <- as.data.frame(cf_matrix_prop)
df$x <- rownames(df)
df <- reshape2::melt(df, id.vars = "x")

# Create the heatmap using ggplot2
ggplot(df, aes(x = variable, y = x, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("#f7fbff", "#08306b"), labels = percent) +
  theme_minimal() +
  labs(title = "Confusion Matrix", x = "Predicted", y = "True", fill = "Proportion") +
  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank())
# Visualize confusion matrix using a heatmap
heatmap(cm_df, Colv=NA, Rowv=NA, 
        col = cm$byClass$Color,
        xlab = "Predicted", ylab = "Actual")
heatmap(cm_df,
        xlab = "Predicted", ylab = "Actual")
# convert confusion matrix to data frame
cm_df <- as.data.frame.matrix(cm$table)
# add rownames and colnames to data frame
cm_df$Actual <- colnames(cm_df)
cm_df$Prediction <- rownames(cm_df)

cf_matrix_prop$Actual <- colnames(cf_matrix_prop)
cf_matrix_prop$Prediction <- rownames(cf_matrix_prop)
# convert data frame to long format
library(tidyr)
cm_df_long <- gather(cf_matrix_prop, key = "Variable", value = "Value", -Prediction, -Actual)

# plot confusion matrix using ggplot2
ggplot(cm_df_long, aes(x = Prediction, y = Actual, fill = Value)) + 
  geom_tile() +
  scale_fill_gradient(low = "red", high = "steelblue", na.value = "black") +
  labs(x = "Predicted", y = "Actual", title = "Confusion Matrix") +
  theme_minimal()
This will create a heatmap of the confusion matrix, with each cell colored according to the number of observations in that cell. You can adjust the low and high arguments of scale_fill_gradient to change the color gradient, and na.value to change the color of cells with missing values.
ggplot(cm_df_long, aes(x = Prediction, y = Actual, fill = Value))+
  scale_fill_gradient(low = "red", high = "steelblue", na.value = "black") 
cf_matrix_prop <- as.data.frame.matrix(cf_matrix_prop)  # Convert to matrix format
heatmap(cf_matrix_prop, annot = TRUE, fmt = ".2%", cmap = "Blues", cexCol = 0.8, cexRow = 0.8)
cf_matrix_prop

# Visualizing Confusion Matrix
fourfoldplot(as.table(cm_df),color=c("green","red"),main = "Confusion Matrix")

ggplot(cm_df_long, aes(x = Actual, y = Prediction, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue", na.value = "white") +
  geom_text(aes(label = paste(Value, "%")), color = "black", size = 4) +
  labs(x = "Actual Labels", y = "Predicted Labels", title = "Confusion Matrix Heatmap") +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), legend.title = element_blank())


# Load library
library(caret)

# Create example data
set.seed(123)
y_true <- sample(c("class A", "class B"), size = 100, replace = TRUE)
y_pred <- sample(c("class A", "class B"), size = 100, replace = TRUE)

# Calculate confusion matrix and F-score
conf_mat <- confusionMatrix(factor(y_pred), factor(y_true))
conf_mat
f_score <- conf_mat$byClass["F1"]
f_score
# Load the caret package
library(caret)

# Create example data with 3 classes
observed <- factor(rep(LETTERS[1:3], c(20, 30, 50)))
predicted <- factor(c(rep("A", 15), rep("B", 25), rep("C", 10), rep("B", 5), rep("C", 40), rep("A", 10)))

# Calculate the confusion matrix and F1 scores for each class
cm <- confusionMatrix(observed, predicted)
f1_scores <- cm$byClass["F1"]
