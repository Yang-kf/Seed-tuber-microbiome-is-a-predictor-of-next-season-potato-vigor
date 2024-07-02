# Data manipulation and analysis
library(dplyr) # Load the dplyr package for data manipulation
library(tibble) # Load the tibble package for data structures
library(MASS) # Load the MASS package for statistical functions
library(infotheo) # Load the infotheo package for information-theoretic analysis
library(magrittr) # Load the magrittr package for piping operator %<>% not loaded by dplyr
library(aricode) # Load the aricode package for handling mutual information calculation

# Plotting and visualization
library(ggplot2) # Load the ggplot2 package for plotting
library(patchwork) # Load the patchwork package for arranging ggplots
library(ggpubr) # Load the ggpubr package for publication-ready plots
library(ggraph) # Load the ggraph package for plotting graphs
library(ggtree) # Load the ggtree package for visualizing phylogenetic trees
library(ggtreeExtra) # Load the ggtreeExtra package for additional functionalities for ggtree
library(ggnewscale) # Load the ggnewscale package for creating secondary axes
library(ggtext) # Load the ggtext package for using formatted text in ggplot2
library(grid) # Load the grid package for grid graphics

# DNA sequence analysis
library(kmer) # Load the kmer package for DNA sequence analysis
library(ape) # Load the ape package for phylogenetic analysis

# Model building and machine learning
library(randomForest) # Load the randomForest package for random forest modeling

# Other
library(ggsci) # Load the ggsci package for scientific color palettes
library(scales) # Load the scales package for scaling and transforming data
library(extrafont) # Load the extrafont package for using external fonts

setwd("~/Postdoc/Song_et_al_2023") # Set the working directory to "~/Postdoc/Song_et_al_2023"
extrafont::font_import("/Windows/Fonts/", pattern = "Arial*", prompt = F, recursive = T) # Import Arial font files from Windows Fonts directory
extrafont::font_import("/Windows/Fonts/", pattern = "arial*", prompt = F, recursive = T) # Import arial font files from Windows Fonts directory
extrafont::loadfonts(device = "win") # Load fonts for use in Windows device
extrafont::loadfonts(device = "postscript") # Load fonts for use in postscript device

##Functions --------------------------------------------------------------------

cleanNames <- function(x) { # Define a function cleanNames taking one argument x
  x %<>% # Pipe x into the following operations
    sub("^X", "", .) %>% # Remove leading "X" from strings in x
    gsub("^S", "", .) %>% # Remove leading "S" from strings in x
    gsub("^R", "", .) %>% # Remove leading "R" from strings in x
    stringi::stri_replace(., regex = "-", replacement = ".", mode = "all") # Replace "-" with "." in strings in x using stringi package
  return(x) # Return the modified x
}

#getp ====

getp <- function(v1, v2, method = "spearman") { # Define a function getp with arguments v1, v2, and method (defaulting to "spearman")
  tryCatch({ # Try the following block of code
    ct <- cor.test(v1, v2, method = method) # Perform correlation test between v1 and v2 using specified method
    return(c(corr = unname(ct$estimate), p = ct$p.value)) # Return correlation coefficient and p-value from the correlation test
  }, error = function(e) return(c(corr = NA, p = NA))) # If an error occurs, return NA values for correlation coefficient and p-value
}

treemapFromDataFrame <- function(df) { # Define a function treemapFromDataFrame taking one argument df
  nc <- ncol(df) # Get the number of columns in df
  if (nc > 2) { # If df has more than 2 columns, do the following
    tmlist <- vector(mode = "list", length <-  nc - 1) # Create a list to store treemaps
    for (i in 1:(nc - 1)) { # Iterate through each column except the last
      tmlist[[i]] <- treemapFromDataFrame(df[i:(i + 1)]) # Recursively call treemapFromDataFrame on subsets of df
    }
    res.tm <- do.call(rbind, tmlist) # Combine the treemaps into a single treemap
    return(res.tm) # Return the resulting treemap
  } else { # If df has 2 or fewer columns, do the following
    df[[1]] %<>% as.character() # Convert the first column of df to character type
    df[[2]] %<>% as.character() # Convert the second column of df to character type
    froms <- unique(df[[1]]) # Get unique values from the first column of df
    fromslist <- sapply(froms, \(x) df[df[[1]] == x, 2], simplify = FALSE) # Create a list of 'to' values corresponding to each 'from' value
    res.df <- lapply(froms, \(x) data.frame(from = x, to = fromslist[[x]])) %>% # Create data frames for each 'from'-'to' pair
      do.call(rbind, .) %>% .[!duplicated.data.frame(.), ] # Combine data frames and remove duplicates
    return(res.df) # Return the resulting data frame
  }
}
asvLabel <- function(a, sel) { # Define a function asvLabel taking two arguments a and sel
  paste0(sel[a, "Kingdom"], " - ", sel[a, "Phylum"]) # Create labels by concatenating Kingdom and Phylum from sel dataframe
}

asvLabel2 <- function(a, slabels) { # Define a function asvLabel2 taking two arguments a and slabels
  as.character(a) %>% # Convert a to character and pipe it into the following operations
    sapply(\(x) paste0(slabels[x], # Apply paste0 function to each element x in a using slabels
                       "<br>(",
                       stringi::stri_sub(x, 1, 5), # Extract first 5 characters of x
                       "...",
                       stringi::stri_sub(x, -5, -1), # Extract last 5 characters of x
                       ")"))
}

asvLabel3 <- function(a, slabels) { # Define a function asvLabel3 taking two arguments a and slabels
  as.character(a) %>% # Convert a to character and pipe it into the following operations
    sapply(\(x) paste0(slabels[x], # Apply paste0 function to each element x in a using slabels
                       "<br>ASV ",
                       stringi::stri_sub(x, 1, 5))) # Extract first 5 characters of x and concatenate with ASV label
}

jaccard <- function(a, b) { # Define a function jaccard taking two arguments a and b
  sum((as.logical(a) & as.logical(b)))/sum((as.logical(a) | as.logical(b))) # Calculate Jaccard similarity coefficient
}

generate_label_df <- function(TUKEY, variable){ # Define a function generate_label_df taking two arguments TUKEY and variable
  Tukey.levels <- TUKEY[[variable]][,4] # Extract levels from Tukey post-hoc for the specified variable
  Tukey.labels <- data.frame(multcompView::multcompLetters(Tukey.levels)['Letters']) # Generate labels using multcompView package
  Tukey.labels$treatment=rownames(Tukey.labels) # Add 'treatment' column with row names as values
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ] # Order Tukey.labels dataframe by 'treatment'
  return(Tukey.labels) # Return the ordered Tukey.labels dataframe
}

pValBootstrapped <- function(cm, nrep = 1000, obs_class, prd_class) { # Define a function pValBootstrapped taking four arguments cm, nrep, obs_class, and prd_class
  f <- 1-((cm/rowSums(cm))[obs_class, prd_class]) # Calculate F-score
  sizes <- rep(names(rowSums(cm)), rowSums(cm)) # Create a vector of sample sizes based on row sums of cm

  fullboot <- sapply(1:nrep, \(x) { # Perform bootstrapping
    boot <- sample(sizes, sum(rowSums(cm)), replace = T) # Sample with replacement from sizes
    1 - (table(boot)/sum(table(boot)))[prd_class] # Calculate F-score for bootstrapped samples
  })

  p_value <- sum(fullboot >= f)/length(fullboot) # Calculate p-value
  return(p_value %>% round(., digits = 3)) # Return rounded p-value
}

cssNormalise <- function(X, samples_in_rows = T) { # Define a function cssNormalise taking two arguments X and samples_in_rows (defaulting to TRUE)
  X.css <- if (samples_in_rows) t(X) else X # Transpose X if samples_in_rows is TRUE
  X.css <- metagenomeSeq::newMRexperiment(X.css) # Convert X.css to MRexperiment object
  X.css <- metagenomeSeq::cumNorm(X.css, p = metagenomeSeq::cumNormStatFast(X.css)) # Perform cumulative normalization
  X.css <- data.frame(metagenomeSeq::MRcounts(X.css, norm = T, log = T)) # Convert normalized counts to data frame
  return(X.css) # Return normalized data frame
}

##Datasets ---------------------------------------------------------------------

#RF objects
rf_eye <- readRDS("RF_results/rf_eye.RDS") # Read random forest object for eye data
rf_heel <- readRDS("RF_results/rf_heel.RDS") # Read random forest object for heel data

rf_eye_df <- rf_eye$importance %>% as.data.frame() %>% .[order(.$`%IncMSE`, decreasing = T), ] # Convert importance measures to data frame and order by %IncMSE
rf_heel_df <- rf_heel$importance %>% as.data.frame() %>% .[order(.$`%IncMSE`, decreasing = T), ] # Convert importance measures to data frame and order by %IncMSE

#Taxonomy info for 16S and ITS ====
tstB <- read.delim("RF_results/taxonomy_bac.tsv", col.names = c("Feature.ID", "Taxon", "Confidence")) %>% cbind(., "Type" = "Bacteria") # Read bacterial taxonomy information
tstF <- read.delim("RF_results/taxonomy_fun.tsv", col.names = c("Feature.ID", "Taxon", "Confidence")) %>% cbind(., "Type" = "Fungi") # Read fungal taxonomy information
allTaxa <- rbind(tstB, tstF) # Combine bacterial and fungal taxonomy information
rownames(allTaxa) <- allTaxa$Feature.ID # Set row names of allTaxa to Feature.ID column values

#RF results ====
rf.b <- read.delim("RF_results/Eye_all_ASV_importance_rf.txt") %>% rownames_to_column(var = "Feature.ID") # Read ASV importance results for eye data
colnames(rf.b)[2] <- "MSE" # Rename second column to "MSE"
rownames(rf.b) <- rf.b$Feature.ID # Set row names of rf.b to Feature.ID column values
rf.b <- cbind(rf.b, allTaxa[rf.b$Feature.ID, ]) # Combine rf.b with taxonomy information
rf.b[, c("MSE", "IncNodePurity")] <- randomForest::importance(x = rf_eye, scale = T) %>% as.data.frame %>% .[rownames(rf.b), ] # Calculate importance measures for rf.b
rf.b <- rf.b[order(rf.b$MSE, decreasing = T), ] # Order rf.b by MSE in decreasing order
rf.df <- data.frame(x = 1:nrow(rf.b), y = rf.b$MSE, row.names = rownames(rf.b)) # Create data frame rf.df with x and y columns
rf.df <- rf.df[rf.df$y > 0, ] # Filter rows of rf.df where y is greater than 0
rf.df$Pur <- rf_eye$importance[rownames(rf.df), "IncNodePurity"] # Add IncNodePurity values to rf.df
rf.df$SD <- rf_eye$importanceSD[rownames(rf.df)] # Add standard deviation values to rf.df
identical(rownames(rf_eye$importance), names(rf_eye$importanceSD)) # Check if row names of rf_eye$importance are identical to names of rf_eye$importanceSD

#ASVs for 16S and ITS - Eye ====
av <- data.table::fread("RF_results/Merge_Eye_bac8000_fun4000_feature-table_ASV.txt") %>% as.data.frame # Read ASV data for eye
av_eye <- av
colnames(av)[1] <- "OTU.ID" # Rename first column to "OTU.ID"
rownames(av) <- av$OTU.ID # Set row names of av to OTU.ID column values
saveRDS(av[, -1], "av.RDS") # Save av data frame without the first column
css.av <- cssNormalise(av[, -1], samples_in_rows = F) # Normalize av data frame
saveRDS(css.av, "css.av.RDS") # Save normalized av data frame
#ASVs for 16S and ITS - Heel ====
av_heel <- data.table::fread("RF_results/Merge_Heel_bac8000_fun4000_feature-table_ASV.txt") %>% as.data.frame # Read ASV data for heel
colnames(av_heel)[1] <- "OTU.ID" # Rename first column to "OTU.ID"
rownames(av_heel) <- av_heel$OTU.ID # Set row names of av_heel to OTU.ID column values

#DNA sequences ====
dna <- ape::read.FASTA("RF_results/dna_sequences_full.fasta") # Read DNA sequences from FASTA file

#NOT RUN
#km5 <- kmer::kcount(dna, k = 5)
#avt <- t(av[, -1])
#km5 <- km5[colnames(avt), ]

#Metadata datasets ====
vitBlue <- list.files("vitality_correlation_figure/", pattern = ".csv", full.names = TRUE) %>% # List csv files in directory
  lapply(\(x) read.csv(x, header = TRUE, col.names = c("Row", "Variety", "Batch", "Blue", "Method2", "Method1")) %>%
           dplyr::mutate(field = basename(x) %>% strsplit('_') %>% unlist %>% '[['(1) %>% '['(1), # Add 'field' column based on file name
                         year = basename(x) %>% strsplit('_') %>% unlist %>% '[['(1) %>% '['(2), # Add 'year' column based on file name
                         Batchyear = paste0(Batch, "-", year)) # Add 'Batchyear' column by combining Batch and year
  ) %>% do.call(rbind, .) # Combine data frames into one

md <- read.csv("RF_results/metadata_Merge_Eye_bac_UP.csv") # Read metadata file for eye data
md$SampleID %<>% sub("^S", "", .) %>% sub("^X", "", .) %>% sub("^R", "", .) # Remove leading 'S', 'X', and 'R' from SampleID
rownames(md) <- md$SampleID # Set row names of metadata to SampleID
md$Year2 <- strsplit(md$Year, "-") %>% sapply('[[', 2) # Extract year from Year column
md %<>% dplyr::mutate(Batchyear = paste0(OriginalNo, "-", Year2)) # Add Batchyear column
md$M <- md %>% '[['("Batchyear") %>% lapply(\(x) vitBlue[vitBlue$Batchyear == x & vitBlue$field == "M", "Method2"]) %>% unlist # Add M column based on vitBlue
md$S <- md %>% '[['("Batchyear") %>% lapply(\(x) vitBlue[vitBlue$Batchyear == x & vitBlue$field == "S", "Method2"]) %>% unlist # Add S column based on vitBlue
md$V <- md %>% '[['("Batchyear") %>% lapply(\(x) vitBlue[vitBlue$Batchyear == x & vitBlue$field == "V", "Method2"]) %>% unlist # Add V column based on vitBlue
md$MBlue <- md %>% '[['("Batchyear") %>% lapply(\(x) vitBlue[vitBlue$Batchyear == x & vitBlue$field == "M", "Blue"]) %>% unlist # Add MBlue column based on vitBlue
md$SBlue <- md %>% '[['("Batchyear") %>% lapply(\(x) vitBlue[vitBlue$Batchyear == x & vitBlue$field == "S", "Blue"]) %>% unlist # Add SBlue column based on vitBlue
md$VBlue <- md %>% '[['("Batchyear") %>% lapply(\(x) vitBlue[vitBlue$Batchyear == x & vitBlue$field == "V", "Blue"]) %>% unlist # Add VBlue column based on vitBlue
md %<>%
  dplyr::group_by(Variety, Year) %>%
  dplyr::mutate(M2 = scale(MBlue) %>% as.numeric(), # Add scaled M2 column
                S2 = scale(SBlue) %>% as.numeric(), # Add scaled S2 column
                V2 = scale(VBlue) %>% as.numeric()) %>% # Add scaled V2 column
  as.data.frame() # Convert to data frame
rownames(md) <- md$SampleID # Set row names to SampleID

md2 <- read.csv("RF_results/metadata_Merge_Eye_bac_precrop_correctSoiltype_NoOutlier_X.csv") # Read metadata file for eye data
md2$SampleID %<>% sub("^S", "", .) %>% sub("^X", "", .) %>% sub("^R", "", .) # Remove leading 'S', 'X', and 'R' from SampleID
rownames(md2) <- md2$SampleID # Set row names of metadata to SampleID

#Wrong
md3 <- read.delim("RF_results/metadata_Merge_Eye_bac_precrop_correctSoiltype_NoOutlier_X.csv") # Read metadata file for eye data
md3$SampleID %<>% sub("^S", "", .) %>% sub("^X", "", .) %>% sub("^R", "", .) # Remove leading 'S', 'X', and 'R' from SampleID
rownames(md3) <- md3$SampleID # Set row names of metadata to SampleID

md_heel <- read.delim("RF_results/metadata_Merge_Heel_bac_precrop_correctSoiltype_BLUE_NoOutlier_X.txt") # Read metadata file for heel data
md_heel$SampleID %<>% sub("^S", "", .) %>% sub("^X", "", .) %>% sub("^R", "", .) # Remove leading 'S', 'X', and 'R' from SampleID
rownames(md_heel) <- md_heel$SampleID # Set row names of metadata to SampleID
md_heel$M <- md_heel$M_52_method2 # Rename column
md_heel$S <- md_heel$S_48_method2 # Rename column
md_heel$V <- md_heel$V_50_method2 # Rename column

#Prediction file ====
prd_files <- list.files(path = "RF_results/", pattern = "Observed*", full.names = TRUE) # List prediction files
prd <- lapply(prd_files, \(x) { # Read and process prediction files
  y <- read.delim(x, sep = " ") # Read file
  field <- strsplit(basename(x), "_") %>% unlist %>% '[['(8) %>% gsub("[0-9]", "", .) # Extract field information
  year <- strsplit(basename(x), "_") %>% unlist %>% '[['(7) %>% '=='("test") %>% ifelse(yes = "2019", no = "2020") # Extract year information
  y$Field <- field # Add Field column
  y$Year <- year # Add Year column
  y$Class <- y$Observed_class # Add Class column
  y$Sample <- rownames(y) # Add Sample column
  y[, c("Sample", "Class", "Field", "Year")] # Select required columns
})
prd <- do.call(rbind, prd) # Combine prediction files
prd$Sample %<>% sub("^S", "", .) %<>% sub("^X", "", .) %<>% sub("^R", "", .) # Remove leading 'S', 'X', and 'R' from Sample
#Order ASVs by MSE
o <- order(rf.b$MSE, decreasing = TRUE) # Order rf.b by MSE in descending order
rf.b %<>% .[o, ] # Reorder rf.b
av %<>% .[rf.b$Feature.ID, ] # Reorder av
#Check, should be TRUE
identical(rf.b$Feature.ID, av$OTU.ID) # Check if Feature.ID in rf.b is identical to OTU.ID in av
#Matrix of ASV relative abundances ====
m.av <- av %>% dplyr::select(!OTU.ID) %>% apply(2, \(x) x/sum(x)) # Calculate relative abundances of ASVs
#Cutoffs for top 5% and 1% contributors
q0 <- quantile(rf.b$MSE, .01) # Calculate 1st percentile of MSE values
ar0 <- which(rf.b$MSE <= q0) %>% length # Calculate length of values less than or equal to 1st percentile
q1 <- quantile(rf.df$y, .95) # Calculate 95th percentile of y values
ar1 <- which(rf.df$y >= q1) %>% length # Calculate length of values greater than or equal to 95th percentile
q2 <- quantile(rf.df$y, .99) # Calculate 99th percentile of y values
ar2 <- which(rf.df$y >= q2) %>% length # Calculate length of values greater than or equal to 99th percentile
q3 <- quantile(rf.df$y, .90) # Calculate 90th percentile of y values
ar3 <- which(rf.df$y >= q3) %>% length # Calculate length of values greater than or equal to 90th percentile

ar <- ar2 # Set ar to length of values greater than or equal to 99th percentile

#Taxonomy matrix ====
sel <- rf.b # Select rf.b for taxonomy matrix
tax <- lapply(sel$Taxon, \(x) strsplit(x, "[a-zA-z](_[0-9]){0,1}__", perl = TRUE) %>% unlist %>% lapply(strsplit, ";") %>% unlist) # Split Taxon column by "__"
taxMat <- matrix("Unclassified", nrow = nrow(sel), ncol = 7) # Create taxonomy matrix with all values "Unclassified"
colnames(taxMat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") # Set column names
for(i in seq_along(tax)) taxMat[i, 1:length(tax[[i]])] <- tax[[i]] # Fill taxonomy matrix with taxonomic information

ambiguous_names <- c("Ambiguous_taxa", "uncultured bacterium", "uncultured", "unidentified", "Unassigned") # Define ambiguous names
taxMat[taxMat %in% ambiguous_names] <- "Unclassified" # Replace ambiguous names with "Unclassified"
#add some info to unclassified taxonomies (optional)
unclInfo <- which(taxMat == "Unclassified", arr.ind = TRUE) %>% .[order(.[, "row"]), ] %>% data.frame %>% .[!duplicated(.$row), ] # Get info for unclassified taxonomies
unclInfo$col %<>% '-'(1) # Subtract 1 from col
taxMatUnlcInfo <- taxMat # Create a copy of taxMat
for (i in 1:nrow(unclInfo)) {x <- unclInfo[i, , drop = TRUE]; taxMatUnlcInfo[x[[1]], (x[[2]] + 1):ncol(taxMatUnlcInfo)] <- paste("Unclassified", taxMatUnlcInfo[x[[1]], x[[2]]])} # Fill taxMatUnlcInfo with taxonomic information
taxMatUnlcInfo <- cbind(taxMatUnlcInfo, "ASV" = sel$Feature.ID) # Add ASV column
rownames(taxMatUnlcInfo) <- taxMatUnlcInfo[, "ASV"] # Set row names to ASV

sel <- cbind(sel, taxMatUnlcInfo) # Combine sel with taxMatUnlcInfo

#Final DF with all info ====
m2 <- reshape2::melt(m.av[1:ar, ]) %>% # Reshape m.av to long format
  rownames_to_column(var = "Sample") %>% # Add Sample column
  dplyr::mutate(Sample = Var2 %>% cleanNames(), # Clean Sample names
                ASV = Var1) %>% # Add ASV column
  dplyr::select(!c(Var1, Var2)) %>% # Select required columns
  dplyr::mutate(Var = md2[Sample, "VarietyFullName"], # Add Var column
                Lot = md2[Sample, "Seedlot"], # Add Lot column
                Field = "M", # Add Field column
                Vit = md[Sample, "MBlue"], # Add Vit column
                Year = md2[Sample, "Year"] %>% strsplit("-") %>% sapply('[', 2), # Add Year column
                Class = lapply(Sample, \(x) prd[prd$Sample == x & prd$Field == "M", "Class", drop = TRUE]) %>% unlist, # Add Class column
                Type = md2[Sample, "SampleType"]) %>% # Add Type column
  rbind(., # Combine with modified versions for Field "S" and "V"
        dplyr::mutate(Field = "S", # Change Field to "S"
                      Vit = md[Sample, "SBlue"], # Add Vit column for "S"
                      Class = lapply(Sample, \(x) prd[prd$Sample == x & prd$Field == "S", "Class", drop = TRUE]) %>% unlist), # Add Class column for "S"
        dplyr::mutate(Field = "V", # Change Field to "V"
                      Vit = md[Sample, "VBlue"], # Add Vit column for "V"
                      Class = lapply(Sample, \(x) prd[prd$Sample == x & prd$Field == "V", "Class", drop = TRUE]) %>% unlist) # Add Class column for "V"
  )

m2 <- m2 %>% # Modify m2
  dplyr::group_by(Year, Field) %>% # Group by Year and Field
  dplyr::mutate(VitCentered0 = (Vit - mean(Vit, na.rm = TRUE)) / sd(Vit, na.rm = TRUE), sd0 = sd(VitCentered0, na.rm = TRUE)) %>% # Calculate VitCentered0 and sd0
  dplyr::group_by(Var, Year, Field) %>% # Group by Var, Year, and Field
  dplyr::mutate(VitCentered = (Vit - mean(Vit, na.rm = TRUE)) / sd(Vit, na.rm = TRUE), sd = sd(VitCentered, na.rm = TRUE)) %>% # Calculate VitCentered and sd
  dplyr::mutate(VitCenteredMed = (Vit - median(Vit, na.rm = TRUE)) / sd(Vit), sd = sd(VitCentered, na.rm = TRUE)) %>% # Calculate VitCenteredMed and sd
  dplyr::mutate(pa = value != 0 %>% ifelse(., "Present", "Absent") %>% factor, # Determine presence or absence and convert to factor
                qb = value <= quantile(value, .9), # Calculate quantile-based value
                qs = value >= quantile(value, .1)) %>% # Calculate quantile-based value
  dplyr::group_by(Year, Field, Var) %>% # Group by Year, Field, and Var
  dplyr::mutate(VitCentered = (Vit - mean(Vit, na.rm = TRUE)) / sd(Vit, na.rm = TRUE)) %>% # Recalculate VitCentered
  dplyr::mutate(VitRel = ifelse(abs(VitCentered) >= sd, ifelse(VitCentered >= 0, "Big", "Small"), NA) %>% factor(., levels = c("Small", "Big"))) # Calculate VitRel and convert to factor

levelnames <- c("H" = "Big", "M" = "Medium", "L" = "Small") # Define levelnames
m2$Class %<>% levelnames[.] %>% factor(., levels = c("Small", "Medium", "Big")) # Replace levels in Class column with defined levelnames and convert to factor
m2$VitRel <- m2$VitCentered # Set VitRel to VitCentered
m2$VitRel[abs(m2$VitRel) < m2$sd] <- NA # Set VitRel to NA where absolute VitCentered is less than sd
m2$VitRel[m2$VitRel < 0] <- 0 # Set VitRel to 0 where VitCentered is less than 0
m2$VitRel[m2$VitRel > 0] <- 1 # Set VitRel to 1 where VitCentered is greater than 0
m2$VitRel %<>% as.logical %>% ifelse(., "Big", "Small") %>% factor(., levels = c("Small", "Big")) # Convert VitRel to factor
mdf <- m2 %>% # Create mdf
  dplyr::filter(!is.na(VitRel)) %>% # Filter out NA values in VitRel
  dplyr::group_by(ASV, pa, VitRel, .drop = FALSE) %>% # Group by ASV, pa, and VitRel
  dplyr::count() %>% # Count occurrences
  dplyr::group_by(ASV) %>% # Group by ASV
  dplyr::mutate(prop = n/sum(n)) %>% # Calculate proportion
  dplyr::mutate(propRel = 0 + (1 - 0) * (prop - min(prop)) / (max(prop) - min(prop))) %>% # Calculate relative proportion
  dplyr::group_by(ASV, VitRel) %>% # Group by ASV and VitRel
  dplyr::mutate(nVitRel = n/sum(n)) %>% # Calculate relative proportion based on VitRel
  dplyr::group_by(ASV, pa) %>% # Group by ASV and pa
  dplyr::mutate(npaRel = n/sum(n)) # Calculate relative proportion based on pa
mdf2 <- m2 %>% # Create mdf2
  dplyr::group_by(ASV, pa, Class, .drop = TRUE) %>% # Group by ASV, pa, and Class
  dplyr::count() %>% # Count occurrences
  dplyr::group_by(ASV, Class, .drop = TRUE) %>% # Group by ASV and Class
  dplyr::mutate(nVitRel = n/sum(n)) %>% # Calculate relative proportion based on VitRel
  dplyr::filter(as.character(Class) %in% c("Big", "Small")) %>% # Filter out rows where Class is not "Big" or "Small"
  dplyr::group_by(ASV, pa, .drop = TRUE) %>% # Group by ASV and pa
  dplyr::mutate(npaRel = n/sum(n)) # Calculate relative proportion based on pa
m2_counts <- reshape2::melt(av[1:ar, ]) %>% # Reshape av to long format
  dplyr::mutate(Sample = variable %>% cleanNames(), # Clean Sample names
                ASV = OTU.ID, # Add ASV column
                value = value) %>% # Add value column
  dplyr::select(!c(variable, OTU.ID)) %>% # Select required columns
  dplyr::mutate(Var = md2[Sample, "VarietyFullName"], # Add Var column
                Lot = md2[Sample, "Seedlot"], # Add Lot column
                Field = "M", # Add Field column
                Vit = md[Sample, "MBlue"], # Add Vit column
                Year = md2[Sample, "Year"] %>% strsplit("-") %>% sapply('[', 2), # Add Year column
                Class = lapply(Sample, \(x) prd[prd$Sample == x & prd$Field == "M", "Class", drop = TRUE]) %>% unlist, # Add Class column for "M"
                Type = md2[Sample, "SampleType"]) %>% # Add Type column
  rbind(., # Combine with modified versions for Field "S" and "V"
        m2_counts %>% dplyr::mutate(
          Field = "S", # Change Field to "S"
          Vit = md[Sample, "SBlue"], # Add Vit column for "S"
          Class = lapply(Sample, \(x) prd[prd$Sample == x & prd$Field == "S", "Class", drop = TRUE]) %>% unlist), # Add Class column for "S"
        m2_counts %>% dplyr::mutate(
          Field = "V", # Change Field to "V"
          Vit = md[Sample, "VBlue"], # Add Vit column for "V"
          Class = lapply(Sample, \(x) prd[prd$Sample == x & prd$Field == "V", "Class", drop = TRUE]) %>% unlist) # Add Class column for "V"
  )

m2_counts <- m2_counts %>% # Modify m2_counts
  dplyr::group_by(Year, Field) %>% # Group by Year and Field
  dplyr::mutate(VitCentered0 = (Vit - mean(Vit, na.rm = TRUE)) / sd(Vit, na.rm = TRUE), sd0 = sd(VitCentered0, na.rm = TRUE)) %>% # Calculate VitCentered0 and sd0
  dplyr::group_by(Var, Year, Field) %>% # Group by Var, Year, and Field
  dplyr::mutate(VitCentered = (Vit - mean(Vit, na.rm = TRUE)) / sd(Vit, na.rm = TRUE), sd = sd(VitCentered, na.rm = TRUE)) %>% # Calculate VitCentered and sd
  dplyr::mutate(VitCenteredMed = (Vit - median(Vit, na.rm = TRUE)) / sd(Vit), sd = sd(VitCentered, na.rm = TRUE)) %>% # Calculate VitCenteredMed and sd
  dplyr::mutate(pa = value != 0 %>% ifelse(., "Present", "Absent") %>% factor, # Determine presence or absence and convert to factor
                qb = value <= quantile(value, .9), # Calculate quantile-based value
                qs = value >= quantile(value, .1)) %>% # Calculate quantile-based value
  dplyr::group_by(Year, Field, Var) %>% # Group by Year, Field, and Var
  dplyr::mutate(VitCentered = (Vit - mean(Vit, na.rm = TRUE)) / sd(Vit, na.rm = TRUE)) %>% # Recalculate VitCentered
  dplyr::mutate(VitRel = ifelse(abs(VitCentered) >= sd, ifelse(VitCentered >= 0, "Big", "Small"), NA) %>% factor(., levels = c("Small", "Big"))) # Calculate VitRel and convert to factor

levelnames <- c("H" = "Big", "M" = "Medium", "L" = "Small") # Define levelnames
m2_counts$Class %<>% levelnames[.] %>% factor(., levels = c("Small", "Medium", "Big")) # Replace levels in Class column with defined levelnames and convert to factor
m2_counts$VitRel <- m2_counts$VitCentered # Set VitRel to VitCentered
m2_counts$VitRel[abs(m2_counts$VitRel) < m2_counts$sd] <- NA # Set VitRel to NA where absolute VitCentered is less than sd
m2_counts$VitRel[m2_counts$VitRel < 0] <- 0 # Set VitRel to 0 where VitCentered is less than 0
m2_counts$VitRel[m2_counts$VitRel > 0] <- 1 # Set VitRel to 1 where VitCentered is greater than 0
m2_counts$VitRel %<>% as.logical %>% ifelse(., "Big", "Small") %>% factor(., levels = c("Small", "Big")) # Convert VitRel to factor

#Broader data
m3 <- reshape2::melt(m.av[1:ar3, ]) %>% # Reshape m.av to long format
  rownames_to_column(var = "Sample") %>% # Add Sample column
  dplyr::mutate(Sample = Var2 %>% cleanNames(), # Clean Sample names
                ASV = Var1) %>% # Add ASV column
  dplyr::select(!c(Var1, Var2)) %>% # Select required columns
  dplyr::mutate(Var = md2[Sample, "VarietyFullName"], # Add Var column
                Lot = md2[Sample, "Seedlot"], # Add Lot column
                Field = "M", # Add Field column
                Vit = md[Sample, "MBlue"], # Add Vit column
                Year = md2[Sample, "Year"] %>% strsplit("-") %>% sapply('[', 2), # Add Year column
                Class = lapply(Sample, \(x) prd[prd$Sample == x & prd$Field == "M", "Class", drop = TRUE]) %>% unlist, # Add Class column for "M"
                Type = md2[Sample, "SampleType"]) %>% # Add Type column
  rbind(., # Combine with modified versions for Field "S" and "V"
        m3 %>% dplyr::mutate(
          Field = "S", # Change Field to "S"
          Vit = md[Sample, "SBlue"], # Add Vit column for "S"
          Class = lapply(Sample, \(x) prd[prd$Sample == x & prd$Field == "S", "Class", drop = TRUE]) %>% unlist
        ),
        m3 %>% dplyr::mutate(
          Field = "V", # Change Field to "V"
          Vit = md[Sample, "VBlue"], # Add Vit column for "V"
          Class = lapply(Sample, \(x) prd[prd$Sample == x & prd$Field == "V", "Class", drop = TRUE]) %>% unlist
        )
  )

m3 <- m3 %>% # Modify m3
  dplyr::group_by(Year, Field) %>% # Group by Year and Field
  dplyr::mutate(VitCentered0 = (Vit - mean(Vit, na.rm = TRUE)) / sd(Vit, na.rm = TRUE), sd0 = sd(VitCentered0, na.rm = TRUE)) %>% # Calculate VitCentered0 and sd0
  dplyr::group_by(Var, Year, Field) %>% # Group by Var, Year, and Field
  dplyr::mutate(VitCentered = (Vit - mean(Vit, na.rm = TRUE)) / sd(Vit, na.rm = TRUE), sd = sd(VitCentered, na.rm = TRUE)) %>% # Calculate VitCentered and sd
  dplyr::mutate(VitCenteredMed = (Vit - median(Vit, na.rm = TRUE)) / sd(Vit), sd = sd(VitCentered, na.rm = TRUE)) %>% # Calculate VitCenteredMed and sd
  dplyr::mutate(pa = value != 0 %>% ifelse(., "Present", "Absent") %>% factor, # Determine presence or absence and convert to factor
                qb = value <= quantile(value, .9), # Calculate quantile-based value
                qs = value >= quantile(value, .1)) # Calculate quantile-based value
levelnames <- c("H" = "Big", "M" = "Medium", "L" = "Small") # Define levelnames
m3$Class %<>% levelnames[.] %>% factor(., levels = c("Small", "Medium", "Big")) # Replace levels in Class column with defined levelnames and convert to factor
m3$VitRel <- m3$VitCentered # Set VitRel to VitCentered
m3$VitRel[abs(m3$VitRel) < m3$sd] <- NA # Set VitRel to NA where absolute VitCentered is less than sd
m3$VitRel[m3$VitRel < 0] <- 0 # Set VitRel to 0 where VitCentered is less than 0
m3$VitRel[m3$VitRel > 0] <- 1 # Set VitRel to 1 where VitCentered is greater than 0
m3$VitRel %<>% as.logical %>% ifelse(., "Big", "Small") %>% factor(., levels = c("Small", "Big")) # Convert VitRel to factor

m3_counts <- reshape2::melt(av[1:1200, ]) %>% # Reshape av to long format
  dplyr::mutate(Sample = variable %>% cleanNames(), # Clean Sample names
                ASV = OTU.ID, # Add ASV column
                value = value) %>% # Add value column
  dplyr::select(!c(variable, OTU.ID)) %>% # Select required columns
  dplyr::mutate(Var = md2[Sample, "VarietyFullName"], # Add Var column
                Lot = md2[Sample, "Seedlot"], # Add Lot column
                Field = "M", # Add Field column
                Vit = md[Sample, "MBlue"], # Add Vit column
                Year = md2[Sample, "Year"] %>% strsplit("-") %>% sapply('[', 2), # Add Year column
                Class = lapply(Sample, \(x) prd[prd$Sample == x & prd$Field == "M", "Class", drop = TRUE]) %>% unlist, # Add Class column for "M"
                Type = md2[Sample, "SampleType"]) %>% # Add Type column
  rbind(., # Combine with modified versions for Field "S" and "V"
        m3_counts %>% dplyr::mutate(
          Field = "S", # Change Field to "S"
          Vit = md[Sample, "SBlue"], # Add Vit column for "S"
          Class = lapply(Sample, \(x) prd[prd$Sample == x & prd$Field == "S", "Class", drop = TRUE]) %>% unlist),
        m3_counts %>% dplyr::mutate(
          Field = "V", # Change Field to "V"
          Vit = md[Sample, "VBlue"], # Add Vit column for "V"
          Class = lapply(Sample, \(x) prd[prd$Sample == x & prd$Field == "V", "Class", drop = TRUE]) %>% unlist)
  )

m3_counts <- m3_counts %>% # Modify m3_counts
  dplyr::group_by(Year, Field) %>% # Group by Year and Field
  dplyr::mutate(VitCentered0 = (Vit - mean(Vit, na.rm = TRUE)) / sd(Vit, na.rm = TRUE), sd0 = sd(VitCentered0, na.rm = TRUE)) %>% # Calculate VitCentered0 and sd0
  dplyr::group_by(Var, Year, Field) %>% # Group by Var, Year, and Field
  dplyr::mutate(VitCentered = (Vit - mean(Vit, na.rm = TRUE)) / sd(Vit, na.rm = TRUE), sd = sd(VitCentered, na.rm = TRUE)) %>% # Calculate VitCentered and sd
  dplyr::mutate(VitCenteredMed = (Vit - median(Vit, na.rm = TRUE)) / sd(Vit), sd = sd(VitCentered, na.rm = TRUE)) %>% # Calculate VitCenteredMed and sd
  dplyr::mutate(pa = value != 0 %>% ifelse(., "Present", "Absent") %>% factor, # Determine presence or absence and convert to factor
                qb = value <= quantile(value, .9), # Calculate quantile-based value
                qs = value >= quantile(value, .1)) # Calculate quantile-based value
levelnames <- c("H" = "Big", "M" = "Medium", "L" = "Small") # Define levelnames
m3_counts$Class %<>% levelnames[.] %>% factor(., levels = c("Small", "Medium", "Big")) # Replace levels in Class column with defined levelnames and convert to factor
m3_counts$VitRel <- m3_counts$VitCentered # Set VitRel to VitCentered
m3_counts$VitRel[abs(m3_counts$VitRel) < m3_counts$sd] <- NA # Set VitRel to NA where absolute VitCentered is less than sd
m3_counts$VitRel[m3_counts$VitRel < 0] <- 0 # Set VitRel to 0 where VitCentered is less than 0
m3_counts$VitRel[m3_counts$VitRel > 0] <- 1 # Set VitRel to 1 where VitCentered is greater than 0
m3_counts$VitRel %<>% as.logical %>% ifelse(., "Big", "Small") %>% factor(., levels = c("Small", "Big")) # Convert VitRel to factor

# ASV vs taxonomy info ====
stripLabels <- taxMatUnlcInfo[, "Genus"]  # Used to be sel$Genus[1:ar2]
stripLabels[!startsWith(stripLabels, "Unclassified")] %<>% paste0("*", ., "*")
stripLabels  %<>% gsub("Unclassified ", "Unclassified<br>", .) %>%
  stringr::str_replace_all(., "\\([:print:]*\\)", "")
stripLabels[!startsWith(taxMatUnlcInfo[, "Family"], "Unclassified") & (startsWith(stripLabels, "Unclassified"))] %<>% gsub("<br>", "<br>*", .) %>% paste0("*")
stripLabels %<>% gsub(" ", "", .)
names(stripLabels) <- taxMatUnlcInfo[, "ASV"]

##### Figures -------------------------------------------------------------------
# MSE curve ====
point_col <- c(
  ggsci::pal_material(palette = "light-blue", n = 20)(20)[16],
  ggsci::pal_material(palette = "light-blue", n = 20)(20)[5],
  ggsci::pal_material(palette = "light-blue", n = 20)(20)[1],
  ggsci::pal_material(palette = "grey", n = 20)(20)[9]
)
point_str <- c(
  ggsci::pal_material(palette = "blue", n = 20)(20)[5],
  ggsci::pal_material(palette = "blue", n = 20)(20)[5],
  ggsci::pal_material(palette = "blue", n = 20)(20)[5],
  ggsci::pal_material(palette = "grey", n = 20)(20)[3]
)

rf.df$color <- ifelse(rf.df$y >= q1, yes = ifelse(rf.df$y >= q2, yes = point_col[1], no = point_col[2]), no = point_col[4])
rf.df$stroke <- ifelse(rf.df$y >= q1, yes = ifelse(rf.df$y >= q2, yes = point_str[1], no = point_str[2]), no = point_str[4])
rf.df$color2 <- ifelse(rf.df$y >= q3, yes = ifelse(rf.df$y >= q1, yes = ifelse(rf.df$y >= q2, yes = point_col[1], no = point_col[2]), no = point_col[3]), no = point_col[4])
rf.df$stroke2 <- ifelse(rf.df$y >= q3, yes = ifelse(rf.df$y >= q1, yes = ifelse(rf.df$y >= q2, yes = point_str[1], no = point_str[2]), no = point_str[3]), no = point_str[4])

rf.df$PurSizes <- cut(rf.df$Pur, breaks = c(0, 1, 2, 3, 4, 5, 6), right = TRUE, c(0, 1, 2, 3, 4, 5)) %>% factor(., levels = c(0, 1, 2, 3, 4, 5) %>% as.character)
# Create a ggplot object for the MSE curve
fig_mse <- ggplot(rf.df, aes(x = x, y = y)) +
  # Add points representing Mean Contribution to Node Purity with varying sizes
  geom_point(aes(size = rf.df$PurSizes %>% as.character() %>% as.numeric() %>% paste0(., "-", .+1) %>% as.factor()), alpha = 1) +
  # Add filled points with different sizes and colors based on the purity levels
  geom_point(size = c("0-1" = 7, "1-2" = 8, "2-3" = 9, "3-4" = 10, "4-5" = 13, "5-6" = 16)[rf.df$PurSizes],
             fill = alpha(rf.df$stroke, alpha = .6),
             color = "white",
             shape = 21, stroke = 0) +
  # Add points with different colors based on the purity levels
  geom_point(color = rf.df$color,
             size = c("0-1" = 5, "1-2" = 6, "2-3" = 7, "3-4" = 8, "4-5" = 11, "5-6" = 14)[rf.df$PurSizes]) +
  # Add text annotations for top 5% and top 1% ASVs
  geom_text(aes(x = 1, y = q1 + .35, label = paste0("Top 5%\n", ar1, " ASVs")),
            family = "Arial", size = 15, hjust = 0, lineheight = .75) +
  geom_hline(yintercept = q1, lty = "dashed", linewidth = 1, color = "grey") + # Add dashed horizontal line at q1
  geom_text(aes(x = 1, y = q2 + .35, label = paste0("Top 1%\n", ar2, " ASVs")),
            family = "Arial", size = 15, hjust = 0, lineheight = .75) +
  geom_hline(yintercept = q2, lty = "dashed", linewidth = 1, color = "grey") + # Add dashed horizontal line at q2
  scale_size_manual(values = c("0-1" = 5, "1-2" = 6, "2-3" = 7, "3-4" = 8, "4-5" = 11, "5-6" = 14)[rf.df$PurSizes %>% as.character() %>% as.numeric() %>% paste0(., "-", .+1)]) + # Scale the size of points manually
  annotation_logticks(sides = "b",
                      long = grid::unit(10, "mm"),
                      mid = grid::unit(5, "mm"),
                      short = grid::unit(3, "mm")) + # Add logarithmic scale ticks
  scale_x_log10() + # Set x-axis to log scale
  ylab("Mean Contribution to Accuracy") + # Set y-axis label
  labs(size = stringr::str_wrap("Mean Contribution to Node Purity", width = 15)) + # Set legend label
  guides(color = "none",size = ggplot2::guide_legend(title.position = "left"),
         shape = "none",
         alpha = "none") + # Hide legends
  theme(axis.title.x = element_blank(), # Remove x-axis title
        axis.title.y = element_text(), # Set y-axis title text properties
        axis.ticks.x = element_blank(), # Remove x-axis ticks
        axis.text.x = element_text(margin = ggplot2::margin(t = 20)), # Set x-axis text margin
        axis.text.y = element_text(margin = ggplot2::margin(r = 20)), # Set y-axis text margin
        text = element_text(family = "Arial", size = 60), # Set text properties
        legend.position = c(.7, .7), # Set legend position
        legend.title.align = 1, # Set legend title alignment
        legend.direction = "vertical", # Set legend direction
        legend.box = "vertical", # Set legend box style
        legend.key = element_rect(fill = "transparent"), # Set legend key style
        legend.background = element_rect(fill = "transparent"), # Set legend background style
        axis.line = element_line(color = "grey20"), # Set axis line color
        panel.grid.major = element_line(colour = "gray90"), # Set major grid line color
        panel.background = element_rect(fill = "white")) # Set panel background color

fig_mse

# Create a ggplot object for the MSE curve with additional annotations
fig_mse2 <- ggplot(rf.df, aes(x = x, y = y)) +
  geom_point(aes(size = rf.df$PurSizes %>% as.character() %>% as.numeric() %>% paste0(., "-", .+1) %>% as.factor()), alpha = 1) +
  geom_point(size = c("0-1" = 7, "1-2" = 8, "2-3" = 9, "3-4" = 10, "4-5" = 13, "5-6" = 16)[rf.df$PurSizes],
             fill = alpha(rf.df$stroke2, alpha = .6),
             color = "white",
             shape = 21, stroke = 0) +
  geom_point(color = rf.df$color2,
             size = c("0-1" = 5, "1-2" = 6, "2-3" = 7, "3-4" = 8, "4-5" = 11, "5-6" = 14)[rf.df$PurSizes]) +
  geom_text(aes(x = 1, y = q1 + .35, label = paste0("Top 5%\n", ar1, " ASVs")),
            family = "Montserrat Medium", size = 14, hjust = 0, lineheight = .75) +
  geom_hline(yintercept = q1, lty = "dashed", linewidth = 1, color = "grey") + # Add dashed horizontal line at q1
  geom_text(aes(x = 1, y = q2 + .35, label = paste0("Top 1%\n", ar2, " ASVs")),
            family = "Montserrat Medium", size = 14, hjust = 0, lineheight = .75) +
  geom_hline(yintercept = q2, lty = "dashed", linewidth = 1, color = "grey") + # Add dashed horizontal line at q2
  geom_text(aes(x = 1, y = q3 - .35, label = paste0("Top 10%\n", ar3, " ASVs")),
            family = "Montserrat Medium", size = 14, hjust = 0, lineheight = .75) +
  geom_hline(yintercept = q3, lty = "dashed", linewidth = 1, color = "grey") + # Add dashed horizontal line at q3
  scale_size_manual(values = c("0-1" = 5, "1-2" = 6, "2-3" = 7, "3-4" = 8, "4-5" = 11, "5-6" = 14)[rf.df$PurSizes %>% as.character() %>% as.numeric() %>% paste0(., "-", .+1)]) + # Scale the size of points manually
  annotation_logticks(sides = "b",
                      long = grid::unit(10, "mm"),
                      mid = grid::unit(5, "mm"),
                      short = grid::unit(3, "mm")) + # Add logarithmic scale ticks
  scale_x_log10() + # Set x-axis to log scale
  ylab("Mean Contribution to Accuracy") + # Set y-axis label
  labs(size = stringr::str_wrap("Mean Contribution to Node Purity", width = 15)) + # Set legend label
  guides(color = "none",size = ggplot2::guide_legend(title.position = "left"),
         shape = "none",
         alpha = "none") + # Hide legends
  theme(axis.title.x = element_blank(), # Remove x-axis title
        axis.title.y = element_markdown(), # Set y-axis title text properties
        axis.ticks.x = element_blank(), # Remove x-axis ticks
        axis.text.x = element_text(margin = ggplot2::margin(t = 20)), # Set x-axis text margin
        axis.text.y = element_text(margin = ggplot2::margin(r = 20)), # Set y-axis text margin
        text = element_text(family = "Montserrat SemiBold", size = 50), # Set text properties
        legend.position = c(.7, .7), # Set legend position
        legend.title.align = 1, # Set legend title alignment
        legend.direction = "vertical", # Set legend direction
        legend.box = "vertical", # Set legend box style
        legend.key = element_rect(fill = "transparent"), # Set legend key style
        legend.background = element_rect(fill = "transparent"), # Set legend background style
        axis.line = element_line(color = "grey20"), # Set axis line color
        panel.grid.major = element_line(colour = "gray90"), # Set major grid line color
        panel.background = element_rect(fill = "white")) # Set panel background color

fig_mse2

#Bubble plot ====
# Bubble plot creation:

# Define taxonomic levels and nomenclature
taxLevels <- c("Kingdom", "Phylum", "Family", "Genus", "Species")
nomenclature <- c("Firmicutes" = "Bacillota",
                  "Proteobacteria" = "Pseudomonadota",
                  "Actinobacteria" = "Actinomycetota",
                  "Bacteroidetes" = "Bacteroidota")

# Subset data to select taxonomic levels and top ASVs
ts <- sel[, c(taxLevels, "MSE")][1:ar1, ]

# Replace "Unclassified" values in Genus with NA
ts$Genus[startsWith(ts$Genus, "Unclassified")] <- NA

# Replace NA values with random strings
ts[is.na(ts)] <- paste0("random_", stringi::stri_rand_strings(sum(is.na(ts)), 20))

# Update Phylum names using nomenclature
ts$Phylum[ts$Phylum %in% names(nomenclature)] %<>% nomenclature[.]
ts$Phylum %<>% paste(ts$Kingdom, "-", .)

# Create edges and vertices for tree map
edges <-  treemapFromDataFrame(ts[, taxLevels[-1]])
vertices <- unique(c(edges[[1]], edges[[2]]))

# Prepare data for nodes
namedisp <- vertices
namedisp[!(namedisp %in% sel[1:ar, "Genus"])] <- NA
sizes <- sapply(vertices, \(x) ts[which(ts == x, arr.ind = T)[1], "MSE", drop = T] %>% sum)
phyla <- sapply(vertices, \(x) ts[which(ts == x, arr.ind = T)[1], "Phylum", drop = T])
family <- sapply(vertices, \(x) ts[which(ts == x, arr.ind = T)[1], "Family", drop = T])
family %<>% names %>% sapply(\(x) x == family[x])
kingdom <- sapply(vertices, \(x) ts[which(ts == x, arr.ind = T)[1], "Kingdom", drop = T])
linewidth <- (!is.na(namedisp))
vertices <- data.frame(name = vertices, size = sizes, phyla = phyla, family = family, nameDisp = namedisp, lineWidth = linewidth)

# Create tidy graph
g <- tidygraph::tbl_graph(vertices, edges)

# Define color order
colorOrder <- sizes[names(sizes) %in% ts$Phylum] %>% names %>% sort %>% sapply(\(x) which(x == (sizes[names(sizes) %in% ts$Phylum] %>% sort(decreasing = T) %>% names)))
pal_circles <- pal_simpsons()(16)[c(10, 8, 6, 7, 9, 1, 11, 13, 14, 16)]

# Set seed for reproducibility
set.seed(6)

# Create the bubble plot using ggraph
fig_bub <- ggraph(g, layout = "circlepack",
                  weight = size/min(size)) +
  geom_node_circle(aes(color = phyla, filter = lineWidth),
                   alpha = 0, linewidth = 2) +
  geom_node_circle(aes(color = phyla, filter = family),
                   alpha = 0, linewidth = .5, linetype = "dashed") +
  geom_node_circle(aes(fill = phyla, color = phyla, filter = !family),
                   alpha = .1, linewidth = .5) +
  geom_node_text(aes(label = nameDisp), size = 25,
                 repel = F, max.overlaps = 15,
                 min.segment.length = .25,
                 force = .01,
                 force_pull = 100,
                 position = "identity",
                 show.legend = F,
                 family = "Arial Bold",
                 fontface = "italic") +
  scale_fill_manual(values = pal_circles) +
  scale_color_manual(values = pal_circles) +
  scale_linewidth_manual() +
  coord_fixed() +
  guides(fill = guide_legend(ncol = 3, byrow = T),
         color = guide_legend(ncol = 3, byrow = T)) +
  theme(panel.background = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        legend.text = element_text(family = "Arial", size = 30),
        legend.position = "top")

# Display the bubble plot
fig_bub

#### Correlation in heatmap ====

# Filter and normalize the data
m2h <- m2 %>%
  .[(.$qb & .$qs), ] %>%
  dplyr::filter(value != 0) %>%
  dplyr::group_by(Var, ASV, Field, Year) %>%
  dplyr::mutate(valueMM = (value - min(value))/(max(value) - min(value)))

# Calculate correlations for each ASV and variable combination
m2h <- sapply(m2h$ASV %>% as.character %>% unique, \(a) {

  x <- m2h[m2h$ASV == a, "VitCentered", drop = T]
  y <- m2h[m2h$ASV == a, "valueMM", drop = T]
  ct <- getp(x, y)

  vars <- sapply(m2h$Var %>% unique, \(v) {
    a %<>% as.character
    x <- m2h[m2h$ASV == a & m2h$Var == v, "VitCentered", drop = T]
    y <- m2h[m2h$ASV == a & m2h$Var == v, "valueMM", drop = T]
    ct <- getp(x, y)
  })

  cbind("All" = ct, vars)
}, simplify = F) %>%
  lapply(data.frame) %>%
  do.call(cbind, .) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "name") %>%
  dplyr::mutate(ASV = strsplit(name, ".", fixed = T) %>% sapply('[', 1) %>% factor(., levels = rev(rownames(m.av)[1:length(unique(.))])),
                Variety = strsplit(name, ".", fixed = T) %>% sapply('[', 2),
                Label = stripLabels[ASV],
                Pformat = paste(round(p, digits = 2), cut(p,
                                                          breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                          labels = c("****", "***", "**", "*", "")),
                                sep = "") %>% ifelse(. == "NANA", yes = NA, no = .),
                Plabel = cut(p,
                             breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                             labels = c("****", "***", "**", "*", "")),
                corrformat = cut(corr, breaks = seq(-1, 1, .2)))

# Generate a phylogenetic tree for ASVs
tree_asv <- dna[rownames(m.av)[1:ar1] %>% unique %>% as.character] %>% kmer::kdistance() %>% hclust() %>% ape::as.phylo()
tree_asv <- ape::keep.tip(phy = tree_asv, tip = rownames(m.av)[1:ar])
t2 <- fortify(tree_asv) %>% subset(., isTip)
t2o <- t2$label[order(t2$y, decreasing = T)]

# Create the correlation heatmap plot
fig_cor <- ggplot(m2h[m2h$ASV %in% rev(levels(m2h$ASV))[1:ar], ] %>%
                    dplyr::mutate(facet = Variety == "All") %>%
                    dplyr::mutate(ASV = factor(as.character(ASV), levels = rev(t2o)),
                                  label = asvLabel3(ASV, stripLabels) %>% factor(levels = asvLabel3(rev(t2o), stripLabels))),
                  aes(x = Variety, y = label, fill = corr)) +
  geom_tile(aes(fill = corr)) +  # Add the heatmap tiles
  geom_text(aes(label = Plabel), size = 25, family = "Arial", vjust = 0.7, color = pal_material("grey", 20)(20)[20]) +  # Add text labels for significance
  geom_segment(x = 1.5, xend = 1.5, y = .5, yend = ar + .5, size = 1, color = pal_material("grey")(10)[8]) +  # Add a segment to separate ASVs
  scale_fill_gradient2(high = pal_material("green", n = 5)(5)[5],  # Set color gradient for correlations
                       low = pal_material("red", n = 5)(5)[5],
                       na.value = pal_material("grey")(10)[7],
                       limits = c(-1, 1),
                       breaks = seq(1, -1, -1)) +
  scale_color_gradient2(high = "black",  # Set color gradient for p-values
                        low = "transparent",
                        mid = "transparent",
                        midpoint = -log10(.1)) +
  labs(fill = "Spearman's ??") +  # Add labels for legend
  guides(fill = guide_colorbar(title.position = "top")) +  # Set legend position
  theme(axis.text.y = element_markdown(hjust = 0, halign = 0, size = 40),  # Customize axis text
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 50, family = "Arial"),
        legend.title = element_text(family = "Arial"),
        legend.text.align = .5,
        legend.title.align = .5,
        legend.key.width = unit(20, "mm"),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.box.just = "center")

# Display the correlation heatmap
fig_cor

##Mutual information in heatmap ====
# Mutual Information calculation and heatmap creation:

# Filter and normalize the data
minf <- m2 %>%
  .[(.$qb & .$qs), ] %>%
  dplyr::filter(value != 0) %>%
  dplyr::group_by(Var, ASV, Field, Year) %>%
  dplyr::mutate(valueMM = (value - min(value))/(max(value) - min(value)))

# Calculate Mutual Information for each ASV and variable combination
minf <- rbind(minf %>%
                dplyr::mutate(Var = "All") %>%
                dplyr::summarise(MI = infotheo::mutinformation(discretize(VitCentered), discretize(valueMM)),
                                 NMI = aricode::NMI(discretize(VitCentered) %>% unlist, discretize(valueMM) %>% unlist, variant = "max"),
                                 .groups = "keep"),
              minf %>%
                dplyr::group_by(ASV, Var, .drop = F) %>%
                dplyr::summarise(MI = ifelse(length(VitCentered) == 0,
                                             yes = NA,
                                             no = mutinformation(discretize(VitCentered), discretize(valueMM))),
                                 NMI = ifelse(length(VitCentered) == 0,
                                              yes = NA,
                                              no = aricode::NMI(discretize(VitCentered) %>% unlist, discretize(valueMM) %>% unlist, variant = "max")),
                                 .groups = "keep")
) %>%  dplyr::mutate(Label = asvLabel2(ASV, stripLabels))

# Create the Mutual Information heatmap plot
fig_cor_mi <- ggplot(minf[minf$ASV %in% rev(levels(minf$ASV))[1:ar], ] %>%
                       dplyr::mutate(facet = Var == "All") %>%
                       dplyr::mutate(ASV = factor(as.character(ASV), levels = rev(t2o))),
                     aes(x = Var, y = Label, fill = MI)) +
  geom_tile() +  # Add the heatmap tiles
  geom_segment(x = 1.5, xend = 1.5, y = .5, yend = ar + .5, size = 1, color = pal_material("grey")(10)[8]) +  # Add a segment to separate ASVs
  scale_y_discrete() +  # Customize y-axis labels
  scale_fill_gradient2(high = pal_material("green", n = 5)(5)[5],  # Set color gradient for Mutual Information values
                       low = pal_material("red", n = 5)(5)[5],
                       na.value = pal_material("grey")(10)[7]) +
  labs(fill = "Mutual Information") +  # Add label for legend
  guides(fill = guide_colorbar(title.position = "top")) +  # Set legend position
  theme(axis.text.y = element_markdown(hjust = 0, halign = 0, size = 40),  # Customize axis text
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 40, family = "Montserrat"),
        legend.title = element_markdown(family = "Arial"),
        legend.text.align = .5,
        legend.title.align = .5,
        legend.key.width = unit(20, "mm"),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.box.just = "center")

# Display the Mutual Information heatmap
fig_cor_mi

#Prevalence in heatmap ====
# Prevalence calculation and heatmap creation:

# Reshape the data for prevalence calculation
m2h2 <- reshape2::melt(m.av[1:ar, ]) %>%
  dplyr::filter(Field == "M") %>%
  rownames_to_column(var = "Sample") %>%
  dplyr::mutate(Sample = Var2 %>% cleanNames(),
                ASV = Var1) %>%
  dplyr::filter(Sample %in% prd$Sample) %>%
  dplyr::mutate(Var = md2[Sample, "VarietyFullName"],
                Lot = md2[Sample, "Seedlot"],
                Year = md2[Sample, "Year"] %>% strsplit("-") %>% sapply('[', 2)) %>%
  dplyr::group_by(Var, ASV, Year) %>%
  dplyr::summarise(prev = sum(value > 0)/length(value), .groups = "keep") %>%
  dplyr::group_by(Var, ASV) %>%
  dplyr::summarise(prev = mean(prev), .groups = "keep")

# Add summary statistics for all varieties
m2h2 %<>%
  rbind(., m2 %>%
          dplyr::filter(Field == "M") %>%
          dplyr::group_by(ASV, Year) %>%
          dplyr::summarise(prev = sum(value > 0)/length(value), .groups = "keep") %>%
          dplyr::group_by(ASV) %>%
          dplyr::summarise(prev = mean(prev), .groups = "keep") %>%
          dplyr::mutate(Var = "All"))

# Create the prevalence heatmap plot
fig_prv <- ggplot(m2h2[m2h2$ASV %in% rev(levels(m2h2$ASV))[1:ar], ] %>%
                    dplyr::mutate(facet = Var == "All") %>%
                    dplyr::mutate(ASV = factor(as.character(ASV), levels = rev(t2o))),
                  aes(x = Var, y = ASV, fill = prev)) +
  geom_tile(aes(fill = prev)) +  # Add the heatmap tiles
  # geom_text(aes(label = Plabel), size = 7, family = "Montserrat SemiBold", vjust = 0.7, color = pal_material("grey", 20)(20)[20]) +
  geom_segment(x = 1.5, xend = 1.5, y = .5, yend = ar + .5, size = 1, color = pal_material("grey")(10)[8]) +  # Add a segment to separate ASVs
  scale_y_discrete(labels = stripLabels) +  # Customize y-axis labels
  scale_fill_gradient2(high = pal_material("green", n = 5)(5)[5],  # Set color gradient for prevalence values
                       low = pal_material("red", n = 5)(5)[5],
                       na.value = pal_material("grey")(10)[7],
                       limits = c(0, 1), label = scales::percent,  # Set limits and format for color bar
                       breaks = seq(1, 0, -.5)) +
  labs(fill = "Prevalence (%)") +  # Add label for legend
  guides(fill = guide_colorbar(title.position = "top")) +  # Set legend position
  theme(axis.text.y = element_markdown(hjust = 0, halign = 0),  # Customize axis text
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 50, family = "Arial"),
        legend.title = element_markdown(family = "Arial"),
        legend.text.align = .5,
        legend.title.align = .5,
        legend.key.width = unit(20, "mm"),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.box.just = "center")

# Display the prevalence heatmap
fig_prv


#Abundance in heatmap ====
# Abundance calculation and heatmap creation:

# Calculate median abundance by ASV and year, then mean abundance by ASV
m2h3 <- m2 %>%
  dplyr::group_by(Var, ASV, Year) %>%
  dplyr::summarise(abn = median(value[value != 0], na.rm = TRUE), .groups = "keep") %>%
  dplyr::group_by(Var, ASV) %>%
  dplyr::summarise(abn = mean(abn, na.rm = TRUE), .groups = "keep")

# Add summary statistics for all ASVs
m2h3 %<>%
  rbind(., m2 %>%
          dplyr::group_by(ASV, Year) %>%
          dplyr::summarise(abn = median(value[value != 0], na.rm = TRUE), .groups = "keep") %>%
          dplyr::group_by(ASV) %>%
          dplyr::summarise(abn = mean(abn, na.rm = TRUE), .groups = "keep") %>%
          dplyr::mutate(Var = "All")) %>%
  dplyr::mutate(abn = ifelse(ASV %in% c("5e2715e1ff601dac60cda12b580ac600", "dbfe8b11a5d5634e4e45cc8e8547ea50"),
                             yes = abn/10, no = abn))  # Adjusting specific ASVs' abundance

# Create the abundance heatmap plot
fig_abn <- ggplot(m2h3[m2h3$ASV %in% rev(levels(m2h3$ASV))[1:ar], ] %>%
                    dplyr::mutate(facet = Var == "All") %>%
                    dplyr::mutate(ASV = factor(as.character(ASV), levels = rev(t2o))),
                  aes(x = Var, y = ASV, fill = abn)) +
  geom_tile(aes(fill = abn)) +  # Add the heatmap tiles
  geom_segment(x = 1.5, xend = 1.5, y = .5, yend = ar + .5, linewidth = 1, color = pal_material("grey")(10)[8]) +  # Add a segment to separate ASVs
  geom_segment(x = 7.6, xend = 7.6, y = .5, yend = 2.5, linewidth = 1, color = pal_material("grey")(10)[8]) +  # Add a segment to separate ASVs
  scale_y_discrete(labels = stripLabels) +  # Customize y-axis labels
  scale_fill_gradient2(high = pal_material("green", n = 5)(5)[5],  # Set color gradient for abundance values
                       low = pal_material("red", n = 5)(5)[5],
                       midpoint = 0,
                       na.value = "white",
                       label = scales::label_percent(scale = 1000, suffix = ""),  # Set format for color bar
                       breaks = seq(round(max(m2h3$abn, na.rm = TRUE), digits = 3), 0, -.001)) +
  labs(fill = "Abundance (???)") +  # Add label for legend
  guides(fill = guide_colorbar(title.position = "top")) +  # Set legend position
  theme(axis.text.y = element_markdown(hjust = 0, halign = 0),  # Customize axis text
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 50, family = "Arial"),
        legend.title = element_markdown(family = "Arial"),
        legend.text.align = .5,
        legend.title.align = .5,
        legend.key.width = unit(20, "mm"),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.box.just = "center")

# Display the abundance heatmap
fig_abn

#Tree ====
# Define a new variable 'tree_asv2' using sapply to apply a function to each element of 'tree_asv$tip.label'
tree_asv2 <- sapply(tree_asv$tip.label, \(x) {
  paste0(sel[x, "Kingdom"],    # Concatenate 'Kingdom' column value with "-"
         " - ",
         sel[x, "Phylum"])     # Concatenate 'Phylum' column value
}) %>% split(tree_asv$tip.label, .) %>%    # Split 'tree_asv$tip.label' by the values obtained from the previous step and group them
  groupOTU(tree_asv, ., group = "Phylum")  # Group the OTUs (Operational Taxonomic Units) based on 'Phylum'

# Generate color information for the tree
colorTree <- ggtree::fortify(tree_asv2) %>%    # Fortify 'tree_asv2' for plotting
  '[['("Phylum") %>%                           # Subset by "Phylum"
  as.character() %>%                           # Convert to character
  sapply(\(x) which(names(colorOrder) == x) %>%   # Apply a function to each element of 'colorOrder' matching 'x'
           ifelse(length(.) == 0, yes = NA, no = .)) %>%   # If no match found, return NA
  pal_circles[.] %>%                           # Subset 'pal_circles' based on the indices obtained from the previous step
  ifelse(is.na(.), yes = "black", no = .)      # Replace NA values with "black"

# Generate tree with colors
tree_asv_color <- ggtree::fortify(tree_asv2) %>%    # Fortify 'tree_asv2' for plotting
  ggtree(., layout = "rectangular",                 # Plot the tree with rectangular layout
         color = colorTree, linewidth = 2) +        # Customize color and line width

  # Generate tree with labels and colors
  tree_asv_color_h <- ggtree::fortify(tree_asv2) %>%  # Fortify 'tree_asv2' for plotting
  dplyr::mutate(label = ifelse(is.na(label),       # Mutate 'label' column based on conditions
                               yes = NA,           # If NA, return NA
                               no = asvLabel3(label, stripLabels))) %>%  # Otherwise, apply 'asvLabel3' function
  ggtree(., layout = "dendrogram",                  # Plot the tree with dendrogram layout
         color = colorTree, linewidth = 2)


# Density plots ====
# Filter data, group by ASV, and calculate summary statistics
m3 <- m2 %>%                                       # Subset 'm2'
  dplyr::filter((value != 0) & qb & qs) %>%        # Filter rows based on conditions
  dplyr::group_by(ASV) %>%                         # Group by 'ASV'
  dplyr::summarise(Above = (100*sum(VitCenteredMed > 0, na.rm = T)/
                              sum(VitCenteredMed != 0, na.rm = T)) %>% floor,   # Calculate 'Above' value
                   Below = 100 - Above,             # Calculate 'Below' value
                   label = asvLabel3(as.character(ASV), stripLabels)  %>%        # Generate labels
                     gsub("Verrucomicro", "Verrucomicro-<br>", .) %>%
                     gsub("Blastoca", "Blastoca-<br>", .) %>%
                     gsub("Sandaraci", "Sandaraci-<br>", .) %>%
                     gsub("Pyronema", "Pyronema-<br>", .) %>%
                     gsub("Comamona", "Comamona-<br>", .))

stripLabels3 <- asvLabel3(as.character(unique(m2$ASV)), stripLabels)  %>%
  gsub("Verrucomicro", "Verrucomicro*-<br>*", .) %>%
  gsub("Blastoca", "Blastoca*-<br>*", .) %>%
  gsub("Sandaraci", "Sandaraci*-<br>*", .) %>%
  gsub("Pyronema", "Pyronema*-<br>*", .) %>%
  gsub("Comamona", "Comamona*-<br>*", .)
# Create a ggplot object with data manipulation and aesthetics
fig_den <- ggplot(m2 %>%
                    .[(.$qb & .$qs), ] %>%
                    dplyr::filter(value != 0) %>%
                    dplyr::group_by(Var, ASV, Year) %>%
                    dplyr::mutate(valueMM = (value - min(value))/(max(value) - min(value))),
                  aes(y = VitCenteredMed, x = valueMM)) +
  # Add filled 2D density plot with specified contour variable
  geom_density_2d_filled(contour_var = "ndensity") +
  # Add robust linear regression line
  stat_smooth(method = MASS::rlm, method.args = list(psi = MASS::psi.bisquare), color = pal_material(palette = "amber", n = 10)(2)[2]) +
  # Add horizontal line at y = 0
  geom_hline(yintercept = 0, color = pal_material(palette = "amber", n = 10)(3)[3]) +
  # Set fill color using viridis color palette
  scale_fill_viridis_d(option = "B") +
  # Add text annotations for data points above a threshold
  geom_text(data = m3,
            mapping = aes(x = .8, y = 1.5, label = paste0(Above, "%"), group = label),
            color = pal_material(palette = "amber", n = 10)(2)[1],
            size = 15) +
  # Add text annotations for data points below a threshold
  geom_text(data = m3,
            mapping = aes(x = .8, y = -1.5, label = paste0(Below, "%"), group = label),
            color = pal_material(palette = "amber", n = 10)(2)[1],
            size = 15) +
  # Add correlation coefficient with significance annotations
  stat_cor(method = "spearman",
           cor.coef.name = "rho",
           mapping = aes(color = NULL,
                         label = paste(after_stat(r.label),
                                       cut(after_stat(p),
                                           breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                           labels = c("'****'", "'***'", "'**'", "'*'", "''")),
                                       sep = "*")),
           color = pal_material(palette = "amber", n = 10)(2)[1],
           p.accuracy = .01,
           r.accuracy = .01,
           label.y = 2.3,
           label.x.npc = .05,
           size = 15,
           hjust = 0) +
  # Set y-axis limits
  ylim(c(-3, 3)) +
  # Set x-axis breaks
  scale_x_continuous(n.breaks = 2) +
  # Wrap facetting by ASV with specified label stripping
  facet_wrap(~ASV, nrow = 2, labeller = labeller(ASV = stripLabels3)) +
  # Set legend title
  labs(fill = "Frequency") +
  # Set x-axis label
  xlab("ASV abundance (scaled)") +
  # Set y-axis label
  ylab("Scaled CSA") +
  # Set theme properties
  theme(text = element_text("Arial", size = 60),
        strip.text = element_markdown(size = 40),
        panel.background = element_rect(fill = "white"),
        legend.position = "none")

fig_den

# Create a second ggplot object with similar specifications but different adjustments
fig_den2 <- ggplot(m2 %>%
                     .[(.$qb & .$qs), ] %>%
                     dplyr::filter(value != 0) %>%
                     dplyr::group_by(Var, ASV, Year) %>%
                     dplyr::mutate(valueMM = (value - min(value))/(max(value) - min(value))),
                   aes(y = VitCenteredMed, x = valueMM)) +
  # Add filled 2D density plot with specified contour variable and adjusted density
  geom_density_2d_filled(contour_var = "ndensity", adjust = .8) +
  # Add robust linear regression line
  stat_smooth(method = MASS::rlm, method.args = list(psi = MASS::psi.bisquare), color = pal_material(palette = "amber", n = 10)(2)[2]) +
  # Add horizontal line at y = 0
  geom_hline(yintercept = 0, color = pal_material(palette = "amber", n = 10)(3)[3]) +
  # Set fill color using viridis color palette
  scale_fill_viridis_d(option = "B") +
  # Add text annotations for data points above a threshold
  geom_text(data = m3,
            mapping = aes(x = .8, y = 1.5, label = paste0(Above, "%"), group = label),
            color = pal_material(palette = "amber", n = 10)(2)[1],
            size = 10) +
  # Add text annotations for data points below a threshold
  geom_text(data = m3,
            mapping = aes(x = .8, y = -1.5, label = paste0(Below, "%"), group = label),
            color = pal_material(palette = "amber", n = 10)(2)[1],
            size = 10) +
  # Add correlation coefficient with significance annotations
  stat_cor(method = "spearman",
           cor.coef.name = "rho",
           mapping = aes(color = NULL,
                         label = paste(after_stat(r.label),
                                       cut(after_stat(p),
                                           breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                           labels = c("'****'", "'***'", "'**'", "'*'", "''")),
                                       sep = "*")),
           color = pal_material(palette = "amber", n = 10)(2)[1],
           p.accuracy = .01,
           r.accuracy = .01,
           label.y = 2.3,
           label.x.npc = .05,
           size = 10,
           hjust = 0) +
  # Set y-axis limits
  ylim(c(-3, 3)) +
  # Set x-axis breaks
  scale_x_continuous(n.breaks = 2) +
  # Wrap facetting by ASV with specified label stripping
  facet_wrap(~ASV, nrow = 2, labeller = labeller(ASV = stripLabels3)) +
  # Set legend title
  labs(fill = "Frequency") +
  # Set x-axis label
  xlab("ASV abundance (normalized)") +
  # Set y-axis label
  ylab("Plant vigor (scaled CSA)") +
  # Set theme properties
  theme(text = element_text("Arial", size = 60),
        strip.text = element_markdown(size = 30),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        axis.title.y = element_markdown())

fig_den2
# Define a layout structure for arranging plots
lay_f3 <- c("
AABBBB
AABBBB
CCCCCC
CCCCCC
CCCCCC
             ")

# Combine multiple ggplot objects and layout structure
figure_1_3 <- fig_mse +
  fig_bub +
  fig_cpa + #theme(plot.margin = ggplot2::margin(0, 100, 0, -800))) +
  plot_layout(design = lay_f3) +  # Arrange plots according to the layout design
  plot_annotation(tag_levels = "A")  # Set annotation tag levels

# Save the combined plot as a JPEG file
jpeg(filename = "Figure_1_full_v3.jpg", width = 3000, height = 3000, units = "px", quality = 300)
figure_1_3  # Display the combined plot
dev.off()  # Turn off JPEG device

##Partial dependence of ASVs ====
# Define unique identifiers for ASVs
cellvibrio1 <- "c205d92bea0c6489fc81746eb5be60fb"
cellvibrio2 <- "08f91cc4ac0790daa33ca1b7910934c9"
streptomyces1 <- "6c8e8f15c8b62edfc17dd5fc30ecd8cd"
streptomyces2 <- "3bc9a5ccb2e7b4122351bd8d4098515f"
acinetobacter1 <- "bcae8f7d2c42acfe3559c665f68f6320"

# Create a named vector of ASV identifiers
asvs <- c("Streptomyces" = streptomyces1,
          "Acinetobacter" = acinetobacter1,
          "Cellvibrio" = cellvibrio1,
          "Cellvibrio2" = cellvibrio2,
          "Streptomyces2" = streptomyces2)

# Copy the av dataframe
av_test <- av
av_test <- av_test[, -1]

# Create a list to store partial dependence plots
pdlist <- vector(mode = "list", length = ar2)

# Generate partial dependence plots for each ASV
for (i in 1:ar2) pdlist[[i]] <- randomForest::partialPlot(rf_eye,
                                                          av_test[, names(p_raw_eye)] %>% t,
                                                          eval(rownames(sel)[i]), plot = F)
# Store partial dependence plots in pd
pd <- pdlist
names(pd) <- rownames(sel)[1:ar2]
pd %<>% lapply(as.data.frame)

# Add ASV column to each dataframe in pd
for (n in names(pd)) pd[[n]]$ASV <- n

# Combine all dataframes in pd into a single dataframe
pd %<>% do.call(rbind, .)

# Set ASV column as a factor with specified levels
pd$ASV %<>% factor(., levels = rownames(sel)[1:ar2])

# Read saved pd dataframe from RDS file
pd <- readRDS("pdp.RDS")

# Subset pd dataframe for ASVs present in rownames(sel)
pdplot <- pd[pd$ASV %in% rownames(sel)[1:ar2], ] %>%
  # Mutate ASV column as a factor with specified levels and labels
  dplyr::mutate(ASV = factor(ASV, levels = rownames(sel)[1:ar2]),
                col = paste0(sel[ASV, "Kingdom"],
                             " - ",
                             sel[ASV, "Phylum"]))

# Define palette order for ASV phyla
pal_order <- colorOrder
names(pal_order) %<>% strsplit(" - ") %>% sapply('[[', 2)
pal_order["Proteobacteria"] <- pal_circles[6]
# Add color assignments for other phyla
# (Assign colors according to the order of appearance in the plot)
...

# Define palette for partial dependence plot
pal_pdp <- pal_circles[c(2, 6, 6, 3,
                         8, 6, 6, 6,
                         6, 3, 9, 2,
                         6, 2, 1, 3,
                         6, 2, 6, 3,
                         6, 9)]

# Create partial dependence plot
fig_pdp <- ggplot(pdplot %>%
                    dplyr::mutate(label = asvLabel3(ASV, stripLabels) %>%
                                    gsub("Verrucomicro", "Verrucomicro*-<br>*", .) %>%
                                    gsub("Blastoca", "Blastoca*-<br>*", .) %>%
                                    gsub("Sandaraci", "Sandaraci*-<br>*", .) %>%
                                    gsub("Pyronema", "Pyronema*-<br>*", .) %>%
                                    gsub("Comamona", "Comamona*-<br>*", .) %>%
                                    factor(., levels = .[(ASV %>% as.numeric) %>% match(unique(.), .)])),
                  aes(x = x, y = y, color = ASV)) +
  geom_line(linewidth = 5) +
  scale_color_manual(values = pal_pdp) +
  scale_x_log10(breaks = c(10, 100), label = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(breaks = c(-.2, -.1, 0, 1)) +
  annotation_logticks(sides = "b",
                      long = grid::unit(10, "mm"),
                      mid = grid::unit(5, "mm"),
                      short = grid::unit(3, "mm")) +
  ylab("Partial contribution to plant vigor") +
  xlab("ASV counts") +
  guides(color = "none") +
  facet_wrap(~label, nrow = 2, scales = "free_x") +
  theme(legend.text = element_markdown(),
        text = element_text(family = "Arial", size = 60),
        strip.text = element_textbox(family = "Arial",
                                     size = 30, halign = .5),
        # strip.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "grey20"),
        # axis.text.x = element_text(margin = ggplot2::margin(t = 20)),
        # axis.text.y = element_text(margin = ggplot2::margin(r = 20)),
        panel.grid.major = element_line(colour = pal_material("grey")(10)[2])
  )
#Exported like fig_den2, 35 x 15 inches in landscape as PDF
fig_pdp

# Define a layout structure for arranging plots
lay_f7 <- c("
AABBBBB
AABBBBB
CCCCCCC
             ")

# Combine multiple ggplot objects and layout structure
figure_1_7 <-
  (fig_mse) +
  (fig_bub + theme(legend.text = element_text(size = 40),
                   plot.background = element_rect(fill = "transparent"),
                   strip.clip = "off",
                   plot.margin = margin(l = -100, b = -100))) +
  (fig_pdp + theme(strip.text = element_textbox(family = "Arial",
                                                size = 35, halign = .5),
                   strip.background = element_rect(fill = "transparent"))) +
  plot_layout(design = lay_f7) +  # Arrange plots according to the layout design
  plot_annotation(tag_levels = "A")  # Set annotation tag levels

# Export figure as TIFF format
tiff("Figure_5.tif", width = 4000, height = 4000, units = "px", res = 300)
figure_1_7
dev.off()

# Save figure as TIFF format using ggsave
ggsave("Figure_5.tif", plot = figure_1_7,
       device = "tiff", width = 4000, height = 4000, units = "px")


###Figure 5 updated ----
# Figures as above, but using only the top 3 predictors for main figure 5
# Create a new ggplot for the top 3 ASVs
fig_den_top3 <- ggplot(m2 %>%
                         .[(.$qb & .$qs), ] %>%
                         dplyr::filter(value != 0) %>%
                         dplyr::group_by(Var, ASV, Year) %>%
                         dplyr::mutate(valueMM = (value - min(value))/(max(value) - min(value))) %>%
                         dplyr::filter(ASV %in% (m2$ASV[1:3] %>% as.character)),
                       aes(y = VitCenteredMed, x = valueMM)) +
  geom_density_2d_filled(contour_var = "ndensity") +
  stat_smooth(method = MASS::rlm, method.args = list(psi = MASS::psi.bisquare), color = pal_material(palette = "amber", n = 10)(2)[2]) +
  geom_hline(yintercept = 0, color = pal_material(palette = "amber", n = 10)(3)[3]) +
  scale_fill_viridis_d(option = "B") +
  geom_text(data = m3 %>% dplyr::filter(ASV %in% ((m2$ASV[1:3] %>% as.character))),
            mapping = aes(x = .8, y = 1.5, label = paste0(Above, "%"), group = label),
            color = pal_material(palette = "amber", n = 10)(2)[1],
            size = 15) +
  geom_text(data = m3 %>% dplyr::filter(ASV %in% ((m2$ASV[1:3] %>% as.character))),
            mapping = aes(x = .8, y = -1.5, label = paste0(Below, "%"), group = label),
            color = pal_material(palette = "amber", n = 10)(2)[1],
            size = 15) +
  stat_cor(method = "spearman",
           cor.coef.name = "rho",
           mapping = aes(color = NULL,
                         label = paste(after_stat(r.label),
                                       cut(after_stat(p),
                                           breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                           labels = c("'****'", "'***'", "'**'", "'*'", "''")),
                                       sep = "*")),
           color = pal_material(palette = "amber", n = 10)(2)[1],
           p.accuracy = .01,
           r.accuracy = .01,
           label.y = 2.3,
           label.x.npc = .05,
           size = 15,
           hjust = 0) +
  ylim(c(-3, 3)) +
  scale_x_continuous(n.breaks = 2) +
  facet_wrap(~ASV, nrow = 1, labeller = labeller(ASV = stripLabels3)) +
  labs(fill = "Frequency") +
  xlab("ASV abundance (normalized)") +
  ylab("Scaled CSA") +
  theme(text = element_text("Arial", size = 60),
        panel.background = element_rect(fill = "white"),
        legend.position = "none")

fig_den_top3

fig_den2_top3 <- ggplot(m2 %>%
                          .[(.$qb & .$qs), ] %>%
                          dplyr::filter(value != 0) %>%
                          dplyr::group_by(Var, ASV, Year) %>%
                          dplyr::mutate(valueMM = (value - min(value))/(max(value) - min(value))) %>%
                          dplyr::filter(ASV %in% (m2_counts$ASV[1:3] %>% as.character)),
                        aes(y = VitCenteredMed, x = valueMM)) +
  geom_density_2d_filled(contour_var = "ndensity", adjust = .8) +
  #geom_point(alpha = .7, size = 4, color = "white") +
  stat_smooth(method = MASS::rlm, method.args = list(psi = MASS::psi.bisquare),
              color = pal_material(palette = "amber", n = 10)(2)[2], linewidth = 2) +
  geom_hline(yintercept = 0, color = pal_material(palette = "amber", n = 10)(3)[3]) +
  scale_fill_viridis_d(option = "B", labels = scales::percent(seq(0.1, 1, .1))) +
  stat_cor(method = "spearman",
           cor.coef.name = "rho",
           mapping = aes(color = NULL,
                         label = paste(after_stat(r.label),
                                       cut(after_stat(p),
                                           breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                           labels = c("'****'", "'***'", "'**'", "'*'", "''")),
                                       sep = "*")),
           color = pal_material(palette = "amber", n = 10)(2)[1],
           p.accuracy = .01,
           r.accuracy = .01,
           label.y = 2.3,
           label.x.npc = .05,
           size = 17,
           hjust = 0) +
  ylim(c(-3, 3)) +
  scale_x_continuous(n.breaks = 2) +
  facet_wrap(~ASV, nrow = 1, labeller = labeller(ASV = stripLabels3)) +
  guides(fill = guide_legend(title = "Point density (%)")) +
  xlab("ASV abundance (normalized)") +
  ylab("Plant vigor<br>(scaled CSA)") +
  theme(text = element_text("Arial", size = 60),
        panel.background = element_rect(fill = "white"),
        legend.position = "bottom",
        axis.title.y = element_markdown())

fig_den2_top3


fig_den2_top3_points <- ggplot(m2 %>%
                                 .[(.$qb & .$qs), ] %>%
                                 dplyr::filter(value != 0) %>%
                                 dplyr::group_by(Var, ASV, Year) %>%
                                 dplyr::mutate(valueMM = (value - min(value))/(max(value) - min(value))) %>%
                                 dplyr::filter(ASV %in% (m2$ASV[1:3] %>% as.character)),
                               aes(y = VitCenteredMed, x = valueMM)) +
  #geom_density_2d_filled(contour_var = "ndensity") +
  geom_point(alpha = .5, size = 4, color = "black") +
  #geom_smooth(span = 1, linewidth = 2) +
  stat_smooth(method = MASS::rlm, method.args = list(psi = MASS::psi.bisquare), color = pal_material(palette = "amber", n = 10)(10)[10], linewidth = 2) +
  geom_hline(yintercept = 0, color = pal_material(palette = "blue-grey", n = 10)(3)[3]) +
  scale_fill_viridis_d(option = "B") +
  # stat_cor(method = "spearman",
  #          cor.coef.name = "rho",
  #          mapping = aes(color = NULL,
  #                        label = paste(after_stat(r.label),
  #                                      cut(after_stat(p),
  #                                          breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
  #                                          labels = c("'****'", "'***'", "'**'", "'*'", "''")),
  #                                      sep = "*")),
  #          color = pal_material(palette = "grey", n = 10)(10)[10],
  #          p.accuracy = .01,
  #          r.accuracy = .01,
#          label.y = 2.3,
#          label.x.npc = .05,
#          size = 17,
#          hjust = 0) +
ylim(c(-3, 3)) +
  scale_x_continuous(n.breaks = 2) +
  facet_wrap(~ASV, nrow = 1, labeller = labeller(ASV = stripLabels3)) +
  labs(fill = "Frequency") +
  xlab("ASV abundance (normalized)") +
  ylab("Plant vigor<br>(scaled CSA)") +
  theme(text = element_text("Arial", size = 60),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        axis.title.y = element_markdown())

fig_den2_top3_points

lab.fun <- function(x) if (length(x) == 2) 100 * x else ""

fig_den2_top3 <- ggplot(m2 %>%
                          .[(.$qb & .$qs), ] %>%
                          dplyr::filter(value != 0) %>%
                          dplyr::group_by(Var, ASV, Year) %>%
                          dplyr::mutate(valueMM = (value - min(value))/(max(value) - min(value))) %>%
                          dplyr::filter(ASV %in% (m2_counts$ASV[1:3] %>% as.character)),
                        aes(y = VitCenteredMed, x = valueMM)) +
  geom_density_2d_filled(contour_var = "ndensity", adjust = .8) +
  stat_smooth(method = MASS::rlm, method.args = list(psi = MASS::psi.bisquare),
              color = pal_material(palette = "amber", n = 10)(2)[2], linewidth = 2) +
  geom_hline(yintercept = 0, color = pal_material(palette = "amber", n = 10)(3)[3]) +
  scale_fill_viridis_d(option = "B", labels = lab.fun,
                       name = "Point density (%)",
                       guide = guide_colorsteps(show.limits = T)) +
  # stat_cor(method = "spearman",
  #          cor.coef.name = "rho",
  #          mapping = aes(color = NULL,
  #                        label = paste(after_stat(r.label),
  #                                      cut(after_stat(p),
  #                                          breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
  #                                          labels = c("'****'", "'***'", "'**'", "'*'", "''")),
  #                                      sep = "*")),
  #          color = pal_material(palette = "amber", n = 10)(2)[1],
  #          p.accuracy = .01,
  #          r.accuracy = .01,
#          label.y = 2.3,
#          label.x.npc = .05,
#          size = 17,
#          hjust = 0) +
ylim(c(-3, 3)) +
  scale_x_continuous(n.breaks = 2) +
  facet_wrap(~ASV, nrow = 1, labeller = labeller(ASV = stripLabels3)) +
  xlab("ASV abundance (normalized)") +
  ylab("Plant vigor<br>(scaled CSA)") +
  theme(text = element_text("Arial", size = 60),
        panel.background = element_rect(fill = "white"),
        legend.position = "bottom",
        axis.title.y = element_markdown(), legend.key.size = unit(5, "line"), legend.spacing = unit(3, "line"))
#For the legend, export as
fig_den2_top3

fig_pdp_top3 <- ggplot(pdplot %>%
                         dplyr::mutate(label = asvLabel3(ASV, stripLabels) %>%
                                         gsub("Verrucomicro", "Verrucomicro*-<br>*", .) %>%
                                         gsub("Blastoca", "Blastoca*-<br>*", .) %>%
                                         gsub("Sandaraci", "Sandaraci*-<br>*", .) %>%
                                         gsub("Pyronema", "Pyronema*-<br>*", .) %>%
                                         gsub("Comamona", "Comamona*-<br>*", .) %>%
                                         factor(., levels = .[(ASV %>% as.numeric) %>% match(unique(.), .)])) %>%
                         dplyr::filter(ASV %in% (m2$ASV[1:3] %>% as.character)),
                       aes(x = x, y = y, color = ASV)) +
  geom_line(linewidth = 5) +
  scale_color_manual(values = pal_pdp) +
  scale_x_log10(breaks = c(10, 100), label = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(breaks = c(-.2, -.1, 0, 1)) +
  annotation_logticks(sides = "b",
                      long = grid::unit(10, "mm"),
                      mid = grid::unit(5, "mm"),
                      short = grid::unit(3, "mm")) +
  ylab("Partial contribution<br>to plant vigor") +
  xlab("ASV counts") +
  guides(color = "none") +
  facet_wrap(~ASV, nrow = 1, labeller = labeller(ASV = stripLabels3)) +
  theme(legend.text = element_markdown(),
        text = element_text(family = "Arial", size = 60),
        strip.text = element_textbox(family = "Arial",
                                     size = 55, halign = .5),
        # strip.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "grey20"),
        # axis.text.x = element_text(margin = ggplot2::margin(t = 20)),
        # axis.text.y = element_text(margin = ggplot2::margin(r = 20)),
        axis.title.y = element_markdown(),
        panel.grid.major = element_line(colour = pal_material("grey")(10)[2])
  )

fig_pdp_top3



fig_cor_top3 <- ggplot(m2h[m2h$ASV %in% rev(levels(m2h$ASV))[1:ar], ] %>%
                         dplyr::mutate(facet = Variety == "All") %>%
                         dplyr::mutate(ASV = factor(as.character(ASV), levels = rev(t2o)),
                                       label = asvLabel3(ASV, stripLabels) %>% factor(levels = asvLabel3(t2o, stripLabels))) %>%
                         dplyr::filter(ASV %in% (m2$ASV[1:3] %>% as.character)),
                       aes(x = Variety, y = label, fill = corr)) +
  geom_tile(aes(fill = corr)) +
  geom_text(aes(label = Plabel), size = 25, family = "Arial", vjust = 0.7, color = pal_material("grey", 20)(20)[20]) +
  geom_segment(x = 1.5, xend = 1.5, y = .5, yend = 3.5, size = 2, color = pal_material("grey")(10)[9]) +
  # scale_y_discrete(labels = stripLabels) +
  scale_fill_gradient2(high = pal_material("green", n = 5)(5)[5],
                       low = pal_material("red", n = 5)(5)[5],
                       na.value = pal_material("grey")(10)[7],
                       limits = c(-1, 1),
                       breaks = seq(1, -1, -1)) +
  scale_color_gradient2(high = "black",
                        low = "transparent",
                        mid = "transparent",
                        midpoint = -log10(.1)) +
  labs(title = expression(paste("Spearman's ", italic("??")))) +
  guides(fill = guide_colorbar(title.position = "top")) +
  theme(axis.text.y = element_markdown(hjust = 0, halign = 0, size = 40),
        plot.title = element_text(hjust = .5),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 50, family = "Arial"),
        legend.title = element_blank(),
        legend.text.align = .5,
        legend.title.align = .5,
        legend.key.width = unit(20, "mm"),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.box.just = "center")
fig_cor_top3


fig_prv_top3 <- ggplot(m2h2[m2h2$ASV %in% rev(levels(m2h2$ASV))[1:ar], ] %>%
                         dplyr::mutate(facet = Var == "All") %>%
                         dplyr::mutate(ASV = factor(as.character(ASV), levels = t2o)) %>%
                         dplyr::filter(ASV %in% (m2$ASV[1:3] %>% as.character)),
                       aes(x = Var, y = ASV, fill = prev)) +
  geom_tile(aes(fill = prev)) +
  # geom_text(aes(label = Plabel), size = 7, family = "Montserrat SemiBold", vjust = 0.7, color = pal_material("grey", 20)(20)[20]) +
  geom_segment(x = 1.5, xend = 1.5, y = .5, yend = 3.5, size = 2, color = pal_material("grey")(10)[9]) +
  scale_y_discrete(labels = stripLabels) +
  scale_fill_gradient2(high = pal_material("green", n = 5)(5)[5],
                       low = pal_material("red", n = 5)(5)[5],
                       na.value = pal_material("grey")(10)[7],
                       limits = c(0, 1), label = scales::percent,
                       breaks = seq(1, 0, -.5)) +
  # scale_color_gradient2(high = "black",
  #                       low = "transparent",
  #                       mid = "transparent",
  #                       midpoint = -log10(.1)) +
  labs(title = "Prevalence (%)") +
  guides(fill = guide_colorbar(title.position = "top")) +
  theme(axis.text.y = element_markdown(hjust = 0, halign = 0),
        plot.title = element_text(hjust = .5),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 50, family = "Arial"),
        legend.title = element_blank(),
        legend.text.align = .5,
        legend.title.align = .5,
        legend.key.width = unit(20, "mm"),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.box.just = "center")
fig_prv_top3

m2h3_top3 <- m2h3[m2h3$ASV %in% rev(levels(m2h3$ASV))[1:ar], ] %>%
  dplyr::mutate(facet = Var == "All") %>%
  dplyr::mutate(ASV = factor(as.character(ASV), levels = t2o)) %>%
  dplyr::filter(ASV %in% (m2$ASV[1:3] %>% as.character))

fig_abn_top3 <- ggplot(m2h3_top3,
                       aes(x = Var, y = ASV, fill = abn)) +
  geom_tile(aes(fill = abn)) +
  geom_segment(x = 1.5, xend = 1.5, y = .5, yend = 3.5, linewidth = 2, color = pal_material("grey")(10)[9]) +
  scale_y_discrete(labels = stripLabels) +
  scale_fill_gradient2(high = pal_material("green", n = 5)(5)[5],
                       low = pal_material("red", n = 5)(5)[5],
                       midpoint = 0,
                       na.value = "white",
                       label = scales::label_percent(scale = 1000, suffix = ""),
                       breaks = seq(round(max(m2h3_top3$abn, na.rm = T), digits = 3), 0, -.001)) +
  labs(title = "Abundance (???)") +
  guides(fill = guide_colorbar(title.position = "top")) +
  theme(axis.text.y = element_markdown(hjust = 0, halign = 0),
        plot.title = element_text(hjust = .5),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 50, family = "Arial"),
        legend.title = element_blank(),
        legend.text.align = .5,
        legend.title.align = .5,
        legend.key.width = unit(20, "mm"),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.box.just = "center")
fig_abn_top3

fig_cpa_top3 <- (fig_cor_top3 + theme(axis.text.x = element_text(size = 55), text = element_text(size = 55), legend.text = element_text(size = 50))) +
  (fig_prv_top3 + theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 55), legend.text = element_text(size = 50))) +
  (fig_abn_top3 + theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 55), legend.text = element_text(size = 50))) +
  theme(plot.margin = ggplot2::margin(0, 0, 0, -100000)) +
  plot_layout(nrow = 1)
fig_cpa_top3

lay_figure_5_top3<- c("
AAAABBBBBBBBBB
AAAABBBBBBBBBB
AAAABBBBBBBBBB
AAAABBBBBBBBBB
AAAABBBBBBBBBB
AAAABBBBBBBBBB
AAAABBBBBBBBBB
AAAABBBBBBBBBB
AAAABBBBBBBBBB
AAAABBBBBBBBBB
CCCCCCCDDDDDDD
CCCCCCCDDDDDDD
EEEEEEEEEEEEEE
EEEEEEEEEEEEEE
             ")

figure_5_top3 <-
  (fig_mse) +
  (fig_bub + theme(legend.text = element_text(size = 45),
                   legend.position = "bottom",
                   plot.background = element_rect(fill = "transparent"),
                   strip.clip = "off",
                   plot.margin = ggplot2::margin(l = -100, b = -300))) +
  (fig_den2_top3 + theme(strip.text = element_textbox(family = "Arial",
                                                      size = 40, halign = .5),
                         legend.text = element_text(size = 40),
                         plot.background = element_rect(fill = "transparent"),
                         strip.background = element_rect("grey"))) +
  (fig_pdp_top3 + theme(strip.text = element_textbox(family = "Arial",
                                                     size = 40, halign = .5),
                        strip.background = element_rect(fill = "grey"))) +
  (fig_cpa_top3) +
  plot_layout(design = lay_figure_5_top3) +
  plot_annotation(tag_levels = "A")

figure_5_top3


### SUPPLEMENTARY FIGURE X ----

#For Eye ----
# Assign the av_eye dataset to the variable av
av <- av_eye

# Clean the column names of av using the cleanNames function
colnames(av) %<>% cleanNames()

# Remove the first column from av
av <- av[, -1]

# Filter av for columns corresponding to the year 2019, select non-zero rows, compute row means, and sort the result
pa.av.19 <- (av[, colnames(av) %in% rownames(md)[md$Year2 == "2019"]] != 0) %>% (\(x) x[rowSums(x) != 0, ])(.) %>% rowMeans %>% sort()

# Filter av for columns corresponding to the year 2020, select non-zero rows, compute row means, and sort the result
pa.av.20 <- (av[, colnames(av) %in% rownames(md)[md$Year2 == "2020"]] != 0) %>% (\(x) x[rowSums(x) != 0, ])(.) %>% rowMeans %>% sort()

# Compute the row means of non-zero elements in av and sort the result
pa.av <-  rowMeans(av != 0) %>% sort()

# Create a data frame pa.df with prevalence, ASV names, year, type, importance, and purity
pa.df <- data.frame(
  Prevalence = c(pa.av.19, pa.av.20),
  ASV = c(names(pa.av.19), names(pa.av.20)),
  Year = c(rep("2019", length(pa.av.19)), rep("2020", length(pa.av.20)))
) %>%
  mutate(Type = allTaxa[ASV, "Type"],
         Year = as.factor(Year),
         Importance = rf.b[ASV, "MSE"],
         Purity = rf.b[as.character(ASV), "IncNodePurity"])

# Create a ggplot object prev.sample.1 to visualize the data
prev.sample.1 <- ggplot(pa.df %>%
                          group_by(Type, Year, Prevalence, .drop = F) %>%
                          count(name = "Count") %>%
                          mutate(Year = as.character(Year)) %>%
                          group_by(Type, Year, .drop = F) %>%
                          mutate(NSamples = ifelse(Year == "2019", yes = sum(md$Year2 == "2019"), no = sum(md$Year2 == "2020")),
                                 Year = ifelse(Year == "2019", yes = "Year 1", no = "Year 2"),
                                 PrevalencePC = cut(Prevalence*100, breaks = seq(0, 100, 1),
                                                    labels = paste0(as.character(c(0, 2:100)), "%"))),
                        aes(x = PrevalencePC, y = Count + 1)) +
  facet_wrap(~Type*Year) +
  geom_col(position = "dodge") +
  scale_y_log10(breaks = c(10, 100, 1000), label = scales::trans_format("log10", scales::math_format(10^.x))) +
  ylab("Number of ASVs") +
  xlab("Prevalence in samples (%)") +
  scale_x_discrete(breaks = c("0%", "25%", "50%", "75%", "100%")) +
  annotation_logticks(sides = "l",
                      long = grid::unit(10, "mm"),
                      mid = grid::unit(5, "mm"),
                      short = grid::unit(3, "mm")) +
  theme(legend.text = element_markdown(),
        text = element_text(family = "Arial", size = 30),
        strip.text = element_textbox(family = "Arial",
                                     size = 30, halign = .5),
        panel.background = element_rect(fill = "white"),
        panel.spacing.x = unit(3, "lines"),
        axis.line = element_line(color = "grey20"),
        axis.title.y = element_markdown(),
        panel.grid.major = element_line(colour = pal_material("grey")(10)[2])
  )

# Print the ggplot object prev.sample.1
prev.sample.1


#For Heel ----
# Assign the value of 'av_heel' to the variable 'av'
av <- av_heel

# Clean the column names of 'av' and assign them back to 'av'
colnames(av) %<>% cleanNames()

# Remove the first column from 'av'
av <- av[, -1]

# Reorder 'av' to match the row names of 'rf_heel_df'
av <- av[rownames(rf_heel_df), ]

# Filter 'av' for columns corresponding to '2019' in 'md_heel', remove rows with all zeros, calculate row means, and sort
pa.av.19 <- (av[, colnames(av) %in% rownames(md_heel)[md_heel$Year2 == "2019"]] != 0) %>% (\(x) x[rowSums(x) != 0, ])(.) %>% rowMeans %>% sort()

# Filter 'av' for columns corresponding to '2020' in 'md_heel', remove rows with all zeros, calculate row means, and sort
pa.av.20 <- (av[, colnames(av) %in% rownames(md_heel)[md_heel$Year2 == "2020"]] != 0) %>% (\(x) x[rowSums(x) != 0, ])(.) %>% rowMeans %>% sort()

# Calculate row means of 'av' where 'av' is not zero, and sort
pa.av <-  rowMeans(av != 0) %>% sort()

# Read a tab-delimited file and store the data in 'heel_taxa', converting the "OTU.ID" column to row names
heel_taxa <- read.delim("RF_results/Merge_Heel_bac_all_feature-table_withSruber-with-taxonomy.tsv") %>%
  tibble::column_to_rownames(var = "OTU.ID")

# Keep only the "Taxon" column in 'heel_taxa', without dropping the dimensions
heel_taxa <- heel_taxa[, "Taxon", drop = F]

# Remove all occurrences of "; " from the "Taxon" column in 'heel_taxa'
heel_taxa$Taxon %<>% gsub("; ", "", .)

# Create a new column "Type" by splitting the "Taxon" column on ";" and extracting the first part,
# then removing the "D_[0-9]__" prefix from the extracted string
heel_taxa$Type <- strsplit(heel_taxa$Taxon, ";") %>% sapply('[', 1) %>% gsub("D_[0-9]__", "", ., perl = T)

# Create a data frame 'pa.df' with columns "Prevalence", "ASV", and "Year" using the specified vectors
pa.df <- data.frame(
  Prevalence = c(pa.av.19, pa.av.20),
  ASV = c(names(pa.av.19), names(pa.av.20)),
  Year = c(rep("2019", length(pa.av.19)), rep("2020", length(pa.av.20)))
) %>%
  # Add columns "Type", "Year", "Importance", and "Purity" to 'pa.df' by merging with 'heel_taxa' and 'rf_heel_df'
  mutate(Type = heel_taxa[ASV, "Type"],
         Year = as.factor(Year),
         Importance = rf_heel_df[as.character(ASV), "MSE"],
         Purity = rf_heel_df[as.character(ASV), "IncNodePurity"])

# Replace NA values in the "Type" column of 'pa.df' with corresponding values from 'allTaxa'
pa.df$Type[is.na(pa.df$Type)] <- allTaxa[pa.df$ASV[is.na(pa.df$Type)], "Type"]

# Create a ggplot object 'prev.sample.heel.1' to visualize the data in 'pa.df'
prev.sample.heel.1 <- ggplot(pa.df %>%
                               # Filter rows with non-NA "Type" values and group by "Type", "Year", "Prevalence"
                               filter(!is.na(Type)) %>%
                               group_by(Type, Year, Prevalence, .drop = F) %>%
                               # Count the number of occurrences in each group and rename the count column to "Count"
                               count(name = "Count") %>%
                               # Convert the "Year" column to character
                               mutate(Year = as.character(Year)) %>%
                               # Group by "Type" and "Year" again
                               group_by(Type, Year, .drop = F) %>%
                               # Add a new column "NSamples" based on the value of "Year" and a transformed "Year" column
                               mutate(NSamples = ifelse(Year == "2019", yes = sum(md$Year2 == "2019"), no = sum(md$Year2 == "2020")),
                                      Year = ifelse(Year == "2019", yes = "Year 1", no = "Year 2"),
                                      # Create a "PrevalencePC" column by categorizing "Prevalence" percentages into bins
                                      PrevalencePC = cut(Prevalence*100, breaks = seq(0, 100, 1),
                                                         labels = paste0(as.character(c(0, 2:100)), "%"))),
                             # Set the aesthetics for the plot
                             aes(x = PrevalencePC, y = Count + 1)) +
  # Create separate panels for each combination of "Type" and "Year"
  facet_wrap(~Type*Year) +
  # Create bar plots with log-scaled y-axis
  geom_col(position = "dodge") +
  # Set the y-axis to a log10 scale with specific breaks and labels
  scale_y_log10(breaks = c(10, 100, 1000), label = scales::trans_format("log10", scales::math_format(10^.x))) +
  # Customize the x-axis with specific labels
  scale_x_discrete(breaks = c("0%", "25%", "50%", "75%", "100%")) +
  # Set the y-axis label
  ylab("Number of ASVs") +
  # Set the x-axis label
  xlab("Prevalence in samples (%)") +
  # Add log scale ticks to the left side of the plot
  annotation_logticks(sides = "l",
                      long = grid::unit(10, "mm"),
                      mid = grid::unit(5, "mm"),
                      short = grid::unit(3, "mm")) +
  # Customize the plot theme
  theme(legend.text = element_markdown(),
        text = element_text(family = "Arial", size = 30),
        strip.text = element_textbox(family = "Arial",
                                     size = 30, halign = .5),
        panel.spacing.x = unit(3, "lines"),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "grey20"),
        axis.title.y = element_markdown(),
        panel.grid.major = element_line(colour = pal_material("grey")(10)[2])
  )

# Display the plot
prev.sample.heel.1

#Final figure
lay.sup.fig.x.1 <- c(
  "
  A
  B
  "
)
sup.fig.x.1 <- prev.sample.1 + prev.sample.heel.1 +
  plot_layout(design = lay.sup.fig.x.1) +
  plot_annotation(tag_levels = "a") &
  theme(plot.margin = margin(0, 20, 0, 0))

sup.fig.x.1
ggsave("Figure_supplementary_X1.jpg", plot = sup.fig.x.1, width = 6000, height = 6000, dpi = 300, units = "px")


### SUPPLEMENTARY FIGUE X2 ----
# Calculate the 10th percentile of the non-zero MSE values from rf.b
q_1 <- quantile(rf.b$MSE[rf.b$MSE > 0], .1)
# Determine the number of MSE values less than or equal to the 10th percentile
ar_1 <- which(rf.b$MSE[rf.b$MSE > 0] <= q_1) %>% length

# Calculate the 99th percentile of the non-zero MSE values from rf.b
q_50 <- quantile(rf.b$MSE[rf.b$MSE > 0], .99)
# Determine the number of MSE values greater than the 99th percentile
ar_50 <- which(rf.b$MSE[rf.b$MSE > 0] > q_50) %>% length

# Subset the av matrix for columns corresponding to rows in md where Year2 is 2019, then log-transform the data
log.rav19 <- av[, colnames(av) %in% rownames(md)[md$Year2 == "2019"]] %>% apply(2, \(x) log(x + 1)/log(sum(x + 1)))
# Subset the av matrix for columns corresponding to rows in md where Year2 is 2019, then normalize the data
rav19 <- av[, colnames(av) %in% rownames(md)[md$Year2 == "2019"]] %>% apply(2, \(x) x/sum(x))
# Calculate the median of the non-zero log-transformed values in log.rav19 for each row
m.log.rav19 <- log.rav19 %>% apply(1, \(x) median(x[x != 0]))
# Calculate the 75th percentile of the non-zero log-transformed values in log.rav19 for each row
q75.log.rav19 <- log.rav19 %>% apply(1, \(x) quantile(x[x != 0], probs = .75))
# Calculate the 25th percentile of the non-zero log-transformed values in log.rav19 for each row
q25.log.rav19 <- log.rav19 %>% apply(1, \(x) quantile(x[x != 0], probs = .25))
# Calculate the standard deviation of the non-zero log-transformed values in log.rav19 for each row
sd.log.rav19 <- log.rav19 %>% apply(1, \(x) sd(x[x != 0]))
# Calculate the median of the non-zero normalized values in rav19 for each row
m.rav19 <- rav19 %>% apply(1, \(x) median(x[x != 0]))
# Calculate the 75th percentile of the non-zero normalized values in rav19 for each row
q75.rav19 <- rav19 %>% apply(1, \(x) quantile(x[x != 0], probs = .75))
# Calculate the 25th percentile of the non-zero normalized values in rav19 for each row
q25.rav19 <- rav19 %>% apply(1, \(x) quantile(x[x != 0], probs = .25))
# Calculate the standard deviation of the non-zero normalized values in rav19 for each row
sd.rav19 <- rav19 %>% apply(1, \(x) sd(x[x != 0]))

# Create a data frame with calculated statistics for 2019 and 2020
rav.df <- data.frame(
  AbundanceMedian = c(m.rav19, m.rav20),         # Median abundance
  AbundanceQ75 = c(q75.rav19, q75.rav20),        # 75th percentile abundance
  AbundanceQ25 = c(q25.rav19, q25.rav20),        # 25th percentile abundance
  AbundanceSD = c(sd.rav19, sd.rav20),           # Standard deviation of abundance
  Prevalence = c(pa.av.19[names(m.rav19)], pa.av.20[names(m.rav20)]), # Prevalence values
  ASV = c(names(m.rav19), names(m.rav20)),       # ASV names
  Year = c(rep("2019", length(m.rav19)), rep("2020", length(m.rav20))) # Year labels
) %>%
  # Add columns for Type, Year, Importance, ImportanceQ, and Purity to the data frame
  mutate(Type = factor(allTaxa[ASV, "Type"]),
         Year = as.factor(Year),
         Importance = rf.b[as.character(ASV), "MSE"],
         ImportanceQ = cut(rf.b[as.character(ASV), "MSE"], breaks = c(-Inf, 3, 6, 9, Inf), include.lowest = T),
         Purity = rf.b[as.character(ASV), "IncNodePurity"])

# Merge two reshaped data frames from rav.df by ASV and Type columns
diff.df <- merge.data.frame(rav.df %>%
                              # Reshape data from long to wide format, with ASV and Type as identifiers and AbundanceMedian as values
                              reshape2::dcast(ASV + Type ~ Year, value.var = "AbundanceMedian") %>%
                              # Rename columns "2020" and "2019" to Abn20 and Abn19
                              rename(Abn20 = "2020", Abn19 = "2019"),
                            # Reshape data from long to wide format, with ASV and Type as identifiers and Prevalence as values
                            rav.df %>% reshape2::dcast(ASV + Type ~ Year, value.var = "Prevalence") %>%
                              # Rename columns "2020" and "2019" to Prev20 and Prev19
                              rename(Prev20 = "2020", Prev19 = "2019"),
                            # Merge the two reshaped data frames by ASV and Type columns
                            by = c("ASV", "Type")) %>%
  # Add new columns: AbnDiff, PrevDiff, Importance, and Rank
  mutate(AbnDiff = abs(Abn19 - Abn20),           # Calculate absolute difference in abundance between 2019 and 2020
         PrevDiff = abs(Prev19 - Prev20),        # Calculate absolute difference in prevalence between 2019 and 2020
         Importance = rf.b[ASV, "MSE"],          # Get Importance value from rf.b data frame
         Rank = ifelse(is.na(stripLabels_sup[ASV]), yes = "Rest of ASVs", no = "Top 10% important")) # Assign Rank based on stripLabels_sup

# Calculate standard deviation of prevalence for each row in diff.df
diff.df$SDPrev <- sapply(1:nrow(diff.df), \(i) sd(c(diff.df$Abn19[i], diff.df$Abn20[i])))

# Calculate coefficient of variation of prevalence for each row in diff.df
diff.df$CVPrev <- sapply(1:nrow(diff.df), \(i) diff.df$SDPrev[i] / mean(c(diff.df$Abn19[i], diff.df$Abn20[i])))

# Create a scatter plot of CVPrev vs. Importance
ggplot(diff.df, aes(x = CVPrev, y = Importance)) +
  geom_point()

# Merge two reshaped data frames from rav.df by ASV and Type columns and filter out rows with NA values in Abn19 or Abn20
diff.abn.df <- merge.data.frame(rav.df %>%
                                  # Reshape data from long to wide format, with ASV and Type as identifiers and AbundanceMedian as values
                                  reshape2::dcast(ASV + Type ~ Year, value.var = "AbundanceMedian") %>%
                                  # Rename columns "2020" and "2019" to Abn20 and Abn19
                                  rename(Abn20 = "2020", Abn19 = "2019"),
                                # Reshape data from long to wide format, with ASV and Type as identifiers and Prevalence as values
                                rav.df %>% reshape2::dcast(ASV + Type ~ Year, value.var = "Prevalence") %>%
                                  # Rename columns "2020" and "2019" to Prev20 and Prev19
                                  rename(Prev20 = "2020", Prev19 = "2019"),
                                # Merge the two reshaped data frames by ASV and Type columns
                                by = c("ASV", "Type")) %>%
  # Add new columns: AbnDiff, PrevDiff, and Importance
  mutate(AbnDiff = abs(Abn19 - Abn20),          # Calculate absolute difference in abundance between 2019 and 2020
         PrevDiff = abs(Prev19 - Prev20),       # Calculate absolute difference in prevalence between 2019 and 2020
         Importance = rf.b[ASV, "MSE"]) %>%     # Get Importance value from rf.b data frame
  filter(!is.na(Abn19) | !is.na(Abn20))         # Filter out rows where both Abn19 and Abn20 are NA

# Create a ggplot for important prevalence samples
imp.prev.sample.1 <- ggplot(pa.df %>%
                              # Modify Year column values for readability
                              mutate(Year = ifelse(Year == "2019", yes = "Year 1", no = "Year 2")),
                            aes(x = Prevalence, y = Importance)) +
  # Add points for data with Importance less than or equal to q_50, colored grey
  geom_point(data = pa.df %>%
               mutate(Year = ifelse(Year == "2019", yes = "Year 1", no = "Year 2")) %>%
               filter(Importance <= q_50),
             aes(size = Purity),
             alpha = .2, color = "grey") +
  # Add points for data with Importance greater than q_50, colored red
  geom_point(data = pa.df %>%
               mutate(Year = ifelse(Year == "2019", yes = "Year 1", no = "Year 2")) %>%
               filter(Importance > q_50),
             aes(size = Purity), alpha = 1, color = "red", show.legend = F) +
  # Label y-axis
  ylab("Mean Contribution to Accuracy (MSE)") +
  # Label x-axis
  xlab("Prevalence (%)") +
  # Scale x-axis to percentage format
  scale_x_continuous(label = scales::percent) +
  # Facet the plot by Type and Year
  facet_wrap(~Type*Year) +
  # Add labels to legend
  labs(size = stringr::str_wrap("Mean Contribution<br>to Node Purity", width = 15),
       alpha = stringr::str_wrap("Mean Contribution<br>to Node Purity", width = 15)) +
  # Customize theme
  theme(legend.title = element_markdown(),          # Set legend title to support markdown
        legend.position = "right",                  # Position legend to the right
        text = element_text(family = "Arial", size = 40),      # Set text family and size
        strip.text = element_textbox(family = "Arial",         # Customize facet strip text
                                     size = 40, halign = .5),
        panel.background = element_rect(fill = "white"),       # Set panel background color
        panel.spacing.x = unit(3, "lines"),                    # Set spacing between panels
        axis.line = element_line(color = "grey20"),            # Customize axis lines
        axis.title.y = element_markdown(),                     # Set y-axis title to support markdown
        panel.grid.major = element_line(colour = pal_material("grey")(10)[2])  # Set major grid lines
  )


# Create a ggplot object 'imp.abn.sample.2' using the 'rav.df' data frame
imp.abn.sample.2 <- ggplot(rav.df %>%
                             # Modify 'Year' column values: "2019" to "Year 1", others to "Year 2"
                             mutate(Year = ifelse(Year == "2019", yes = "Year 1", no = "Year 2")), # / NSamples),
                           # Set aesthetic mappings for the plot: x-axis as 'AbundanceMedian' and y-axis as 'Importance'
                           aes(x = AbundanceMedian, y = Importance)) +
  # Create facets based on 'Type' and 'Year'
  facet_wrap(~Type*Year) +
  # Add points to the plot for a subset of data where 'Importance' <= 'q_50'
  geom_point(data = rav.df %>%
               # Modify 'Year' column values: "2019" to "Year 1", others to "Year 2"
               mutate(Year = ifelse(Year == "2019", yes = "Year 1", no = "Year 2")) %>%
               # Filter data to include only rows where 'Importance' <= 'q_50'
               filter(Importance <= q_50),
             # Set aesthetic mappings for size and transparency
             aes(size = Purity),
             alpha = .2, color = "grey") +
  # Add points to the plot for a subset of data where 'Importance' > 'q_50'
  geom_point(data = rav.df %>%
               # Modify 'Year' column values: "2019" to "Year 1", others to "Year 2"
               mutate(Year = ifelse(Year == "2019", yes = "Year 1", no = "Year 2")) %>%
               # Filter data to include only rows where 'Importance' > 'q_50'
               filter(Importance > q_50),
             # Set aesthetic mappings for size and transparency, and remove legend for this geom
             aes(size = Purity), alpha = 1, color = "red", show.legend = F) +
  # Set y-axis label
  ylab("Mean Contribution to Accuracy (MSE)") +
  # Set x-axis label with HTML formatting
  xlab("Median ASV abundance<br>(non-zero values)") +
  # Set labels for size and alpha in the legend with wrapped text
  labs(size = stringr::str_wrap("Mean Contribution<br>to Node Purity", width = 15),
       alpha = stringr::str_wrap("Mean Contribution<br>to Node Purity", width = 15)) +
  # Scale x-axis to log10 and set breaks and labels
  scale_x_log10(breaks = 10^(-c(5:1)), label = scales::percent) +
  # Add log scale ticks to the bottom of the plot
  annotation_logticks(sides = "b",
                      long = grid::unit(10, "mm"),
                      mid = grid::unit(5, "mm"),
                      short = grid::unit(3, "mm")) +
  # Customize the theme of the plot
  theme(legend.title = element_markdown(),
        text = element_text(family = "Arial", size = 40),
        strip.text = element_textbox(family = "Arial",
                                     size = 40, halign = .5),
        panel.background = element_rect(fill = "white"),
        panel.spacing.x = unit(3, "lines"),
        axis.line = element_line(color = "grey20"),
        axis.title.x = element_markdown(),
        panel.grid.major = element_line(colour = pal_material("grey")(10)[2])
  )
# Combine 'imp.prev.sample.1' and 'imp.abn.sample.2' plots, adjust theme and layout, and store in 'imp.prev.abn'
imp.prev.abn <- (imp.prev.sample.1 + theme(legend.position = "none")) + imp.abn.sample.2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
# Display the combined plot
imp.prev.abn

# Create a ggplot object 'prev.diff.imp.sample.1' using the 'diff.df' data frame
prev.diff.imp.sample.1 <- ggplot(diff.df,
                                 # Set aesthetic mappings for the plot: x-axis as 'Prev20' and y-axis as 'Prev19'
                                 aes(x = Prev20, y = Prev19)) +
  # Create facets based on 'Type'
  facet_wrap(~Type) +
  # Add points to the plot for a subset of data where 'Importance' <= 'q_50'
  geom_point(data = diff.df[diff.df$Importance <= q_50, ], size = 4, alpha = .2, color = "grey") +
  # Add a reference line with a slope of 1 and intercept of 0
  geom_abline(slope = 1, intercept = 0, colour = "black") +
  # Add points to the plot for a subset of data where 'Importance' > 'q_50'
  geom_point(data = diff.df[diff.df$Importance > q_50, ], size = 4, alpha = 1, color = "red") +
  # Scale x-axis with percentage labels and set limits
  scale_x_continuous(label = scales::percent, limits = c(0, 1)) +
  # Scale y-axis with percentage labels and set limits
  scale_y_continuous(label = scales::percent, limits = c(0, 1)) +
  # Set x-axis label
  xlab("Prevalence in year 2 (%)") +
  # Set y-axis label
  ylab("Prevalence in year 1 (%)") +
  # Add a statistical correlation layer with custom label formatting
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")), size = 7) +
  # Remove legend for alpha aesthetic
  guides(alpha = "none") +
  # Set aspect ratio to ensure fixed coordinate ratio
  coord_fixed() +
  # Customize the theme of the plot
  theme(legend.title = element_markdown(),
        legend.position = "right",
        text = element_text(family = "Arial", size = 40),
        strip.text = element_textbox(family = "Arial",
                                     size = 40, halign = .5),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "grey20"),
        panel.spacing.x = unit(3, "lines"),
        axis.title.y = element_markdown(),
        panel.grid.major = element_line(colour = pal_material("grey")(10)[2])
  )

# Create a ggplot object 'abn.diff.imp.sample.1' using the 'diff.abn.df' data frame
abn.diff.imp.sample.1 <- ggplot(diff.abn.df,
                                # Set aesthetic mappings for the plot: x-axis as 'Abn20' and y-axis as 'Abn19'
                                aes(x = Abn20, y = Abn19)) +
  # Create facets based on 'Type'
  facet_wrap(~Type) +
  # Add points to the plot for a subset of data where 'Importance' <= 'q_50'
  geom_point(data = diff.abn.df[diff.abn.df$Importance <= q_50, ], size = 4, color = "grey") +
  # Add a reference line with a slope of 1 and intercept of 0
  geom_abline(slope = 1, intercept = 0, colour = "black") +
  # Add points to the plot for a subset of data where 'Importance' > 'q_50'
  geom_point(data = diff.abn.df[diff.abn.df$Importance > q_50, ], size = 4, color = "red") +
  # Set x-axis label with HTML formatting
  xlab("Median non-zero abundance<br>in year 2 (%)") +
  # Set y-axis label with HTML formatting
  ylab("Median non-zero abundance<br>in year 1 (%)") +
  # Add a statistical correlation layer with custom label formatting
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")), size = 7) +
  # Scale x-axis to log10 and set breaks and labels
  scale_x_log10(breaks = 10^(-c(5:1)), label = scales::percent) +
  # Scale y-axis to log10 and set breaks and labels
  scale_y_log10(breaks = 10^(-c(5:1)), label = scales::percent) +
  # Add log scale ticks to the bottom and left sides of the plot
  annotation_logticks(sides = "bl",
                      long = grid::unit(10, "mm"),
                      mid = grid::unit(5, "mm"),
                      short = grid::unit(3, "mm")) +
  # Remove legend for alpha aesthetic
  guides(alpha = "none") +
  # Set aspect ratio to ensure fixed coordinate ratio
  coord_fixed() +
  # Customize the theme of the plot
  theme(legend.title = element_markdown(),
        text = element_text(family = "Arial", size = 40),
        strip.text = element_textbox(family = "Arial",
                                     size = 40, halign = .5),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "grey20"),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        panel.spacing.x = unit(3, "lines"),
        panel.grid.major = element_line(colour = pal_material("grey")(10)[2])
  )

# Creating a character vector with layout design for the plot
lay.sup.fig.x.2 <- c(
  "
  AB
  AB
  CD
  "
)

# Creating a plot object by combining several plots with a specified layout and annotations
sup.fig.x.2 <- (imp.prev.abn) + prev.diff.imp.sample.1 + abn.diff.imp.sample.1 +
  plot_layout(design = lay.sup.fig.x.2, guides = "keep") +
  plot_annotation(tag_levels = "a")

sup.fig.x.2
# Saving the created plot as a JPEG file with specified dimensions and resolution
ggsave("Figure_supplementary_X2.jpg", plot = sup.fig.x.2, width = 9000, height = 7000, dpi = 300, units = "px")


