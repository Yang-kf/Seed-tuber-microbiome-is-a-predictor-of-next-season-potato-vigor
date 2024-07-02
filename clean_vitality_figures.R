# Load the ggplot2 library
library(ggplot2)
# Load the magrittr library
library(magrittr)
# Load the extrafont library
library(extrafont)
# Load the ggsci library
library(ggsci)
# Load the ggtext library
library(ggtext)
# Load the dplyr library
library(dplyr)
# Load the gridExtra library
library(gridExtra)
# Load the emmeans library
library(emmeans)
# Load the GGally library
library(GGally)
# Load the lme4 library
library(lme4)
# Load the patchwork library
library(patchwork)
# Load the performance library
library(performance)

# Set the working directory to "~/Postdoc/Song_et_al_2023"
setwd("~/Postdoc/Song_et_al_2023")
# Assign the path to Ghostscript to the variable gscript
gscript <- "/Program Files/gs/gs10.02.0/bin/gswin64c.exe"
# Set the environment variable R_GSCMD to the path of Ghostscript
# Sys.setenv(R_GSCMD = gscript)

# Create a vector pal_effects containing colors from specified palettes
pal_effects <- c(
  pal_material("green")(10)[5],
  pal_material("grey")(10)[5],
  pal_material("red")(10)[5])

# Importing fonts from the Windows Fonts directory that match the pattern "arial*"
extrafont::font_import("/Windows/Fonts/", pattern = "arial*", prompt = F)

# Loading fonts for the postscript device
extrafont::loadfonts(device = "postscript")

# Loading fonts for the Windows device
extrafont::loadfonts(device = "win")

# Getting a list of CSV files in the "vitality_correlation_figure/" directory and reading them into a list of data frames
vitBlue <- list.files("vitality_correlation_figure/", pattern = ".csv", full.names = T) %>%
  lapply(\(x) read.csv(x, header = T, col.names = c("Row", "Variety", "Batch", "Blue", "Method2", "Method1")) %>%
           dplyr::mutate(field = basename(x) %>% strsplit('_') %>% '[['(1) %>% '['(1),
                         year = basename(x) %>% strsplit('_') %>% '[['(1) %>% '['(2))
  ) %>% do.call(rbind, .)

# Getting a list of files in the "raw_vitality" directory with names matching the pattern "Timeseries*" and reading them into a list of data frames
vit_files <- list.files(path = "raw_vitality", pattern = "Timeseries*", full.names = T)
vit <- lapply(vit_files, \(x) {
  y <- read.csv(x)
  field <- strsplit(basename(x), "_") %>% unlist %>% '['(3)
  year <- strsplit(basename(x), "_") %>% unlist %>% '['(4) %>% sub(".csv", "", .)
  timepoint <- list.files("vitality_correlation_figure/",
                          pattern = paste0(field, "_", year), full.names = T) %>%
    read.csv %>% colnames %>% '['(4)
  vals <- y[, timepoint]
  y$Field <- field
  y$Year <- year
  y <- y[!startsWith(colnames(y), "X")]
  y$Value <- vals
  y
})

# Combining the list of data frames into a single data frame
vit <- do.call(rbind, vit)

# Reading metadata from a CSV file into a data frame
md2 <- read.csv("RF_results/metadata_Merge_Eye_bac_precrop_correctSoiltype_NoOutlier_X.csv")
md2$SampleID %<>% sub("^S", "", .) %>% sub("^X", "", .) %>% sub("^R", "", .)
rownames(md2) <- md2$SampleID

# Modifying values in the 'CorrectBatch' column of the 'vit' data frame
vit$CorrectBatch <- vit$Batch
vit$CorrectBatch[vit$Year == "2019" & vit$Batch > 120] %<>% '-'(80)

# Converting the 'Batch' and 'CorrectBatch' columns to factors
vit$Batch %<>% as.factor
vit$CorrectBatch %<>% as.factor

# Modifying values in the 'CorrectBatch' column of the 'vitBlue' data frame
vitBlue$CorrectBatch <- as.character(vitBlue$Batch)
vitBlue$CorrectBatch[vitBlue$year == "2019"] <- sapply(as.character(vitBlue$Batch)[vitBlue$year == "2019"], \(y) {
  ifelse(y %in% names(x), yes = x[y], no = y)
})

# Extracting 'Blue' values from 'vitBlue' based on matching conditions and assigning them to 'blue'
blue <- sapply(seq(nrow(vit)), \(x) vitBlue[vitBlue$year == vit$Year[x] &
                                              vitBlue$field == vit$Field[x] &
                                              vitBlue$Batch == vit$Batch[x], "Blue"])

# Replacing empty 'blue' values with NA
blue[sapply(blue, length) == 0] <- NA

# Unlisting 'blue'
blue %<>% unlist

# Adding a 'Blue' column to the 'vit' data frame
vit$Blue <- blue

# Grouping and summarizing data in the 'vit' data frame
vit %<>% dplyr::group_by(Variety, Year, Field, Plot.number, CorrectBatch) %>%
  dplyr::summarise(ValueMean = mean(Value, na.rm = T), Blue = unique(Blue)) %>%
  dplyr::group_by(Variety, Year, Field, CorrectBatch) %>%
  dplyr::mutate(ValueCorrected = (Blue + scale(ValueMean, scale = F)) %>% as.numeric) %>%
  dplyr::rename(Batch = CorrectBatch)

##Vitality per seed lot ====
# Load ggplot2 package and subset data for the year 2019
fig_bes <- ggplot(vit[vit$Year == "2019", ] %>%
                    # Group data by Variety, Field, and Year
                    dplyr::group_by(Variety, Field, Year) %>%
                    # Scale the 'ValueCorrected' column without centering
                    dplyr::mutate(vs = scale(ValueCorrected, scale = F)) %>%
                    # Group by Variety and Batch, and calculate mean for each batch
                    dplyr::group_by(Variety, Batch) %>%
                    dplyr::mutate(meanBatch = mean(vs)) %>%
                    # Group only by Variety
                    dplyr::group_by(Variety) %>%
                    # Arrange batches by meanBatch
                    dplyr::arrange(meanBatch, .by_group = T) %>%
                    # Normalize meanBatch values
                    dplyr::mutate(meanBatchNorm = (meanBatch - min(meanBatch))/(max(meanBatch) - min(meanBatch))) %>%
                    # Convert Batch to a factor with reversed levels
                    dplyr::mutate(Batch = factor(Batch, levels = rev(unique(Batch)))),
                  aes(x = Batch, y = vs)) +
  # Create violin plots with meanBatchNorm as fill
  geom_violin(aes(fill = meanBatchNorm)) +
  # Add jittered points colored by Field
  geom_jitter(aes(color = Field), height = 0, width = .2, size = 1) +
  # Add dashed horizontal line at y = 0
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 1, color = pal_material("grey")(10)[7]) +
  ylab("Vitality (spatially corrected)") +
  xlab("Seedlot") +
  # Set fill gradient scale
  scale_fill_gradient2(low = pal_effects[3],
                       mid = pal_effects[2],
                       high = pal_effects[1],
                       midpoint = 0.5) +
  # Set legend guides
  guides(color = guide_legend(override.aes = list(size = 7)),
         fill = "none") +
  # Facet by Variety with free scales
  facet_wrap(~Variety, scales = "free") +
  # Set theme elements
  theme(strip.text = element_markdown(size = 10),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_blank(),
        text = element_text("Montserrat Medium", size = 25),
        axis.text.x = element_text("Montserrat Medium", size = 10, angle = 30, vjust = 1, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.spacing.x = unit(3, "lines"))

##Vitality per seed lot simplified ====
# Load ggplot2 package and subset data for the year 2020
fig_vps <- ggplot(vit[vit$Year == "2020", ] %>%
                    # Mutate Variety to sentence case and create Batchyear column
                    dplyr::mutate(Variety = stringr::str_to_sentence(Variety),
                                  Batchyear = paste0(Batch, "-", Year)) %>%
                    # Filter out NA values in ValueCorrected
                    dplyr::filter(!is.na(ValueCorrected)) %>%
                    # Group data by Variety, Field, and Year
                    dplyr::group_by(Variety, Field, Year) %>%
                    # Scale the 'ValueCorrected' column with centering
                    dplyr::mutate(vs = scale(ValueCorrected, scale = T) %>% as.numeric) %>%
                    # Group by Variety and Batch, and calculate mean for each batch
                    dplyr::group_by(Variety, Batch) %>%
                    dplyr::mutate(meanBatch = mean(vs, na.rm = T)) %>%
                    # Group only by Variety
                    dplyr::group_by(Variety) %>%
                    # Arrange batches by meanBatch
                    dplyr::arrange(meanBatch, .by_group = T) %>%
                    # Normalize meanBatch values
                    dplyr::mutate(meanBatchNorm = (meanBatch - min(meanBatch, na.rm = T))/(max(meanBatch, na.rm = T) - min(meanBatch, na.rm = T))) %>%
                    # Mutate Batch as a factor with reversed levels
                    dplyr::mutate(Batch = factor(Batch, levels = rev(unique(Batch)))) %>%
                    # Group by Year, Variety, Batch, and Field, and summarise vs statistics
                    dplyr::group_by(Year, Variety, Batch, Field) %>%
                    dplyr::summarise(vsmn = mean(vs),
                                     vshi = max(vs),
                                     vslo = min(vs),
                                     meanBatchNorm = meanBatchNorm[1])  %>%
                    # Mutate Field to a factor with specific levels
                    dplyr::mutate(Field = paste0("Field ", Field %>% gsub("S", "K", .)) %>% factor(levels = c("Field M", "Field K", "Field V"))),
                  aes(x = Batch, y = vsmn, shape = Field, color = meanBatchNorm)) +
  # Add stripped columns
  ggstats::geom_stripped_cols(aes(color = NULL), odd = pal_material("grey", 30)(30)[10], even = pal_material("grey", 30)(30)[5]) +
  # Add error bars
  geom_errorbar(aes(x = Batch, ymin = vslo, ymax = vshi), size = 1, position = position_dodge(width = 1), width = 0) +
  # Add points
  geom_point(aes(y = vsmn), size = 4, position = position_dodge(width = 1)) +
  # Add dashed horizontal line at y = 0
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 1, color = pal_material("grey")(10)[7]) +
  ylab("Plant vigor (scaled CSA)") +
  xlab("Seedlot") +
  # Set color gradient scale
  scale_color_gradient2(low = pal_effects[3],
                        mid = pal_effects[2],
                        high = pal_effects[1],
                        midpoint = 0.5) +
  # Facet by Variety with free scales
  facet_wrap(~Variety, scales = "free") +
  # Set legend guides
  guides(shape = guide_legend(override.aes = list(size = 7)),
         color = "none") +
  # Set theme elements
  theme(strip.text = element_markdown(family = "Arial", size = 50),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_blank(),
        text = element_text(family = "Arial", size = 45),
        axis.title.y = element_markdown(family = "Arial", size = 45),
        axis.text.x = element_text(family = "Arial", size = 25, angle = 90, vjust = .5, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.spacing.x = unit(3, "lines"))

fig_vps

# Saved by exporting to PDF in landscape 35 x 12.9 inch using cairo device
jpeg(filename = "Figure_S5_scaled_2019.pdf", width = 2500, height = 1500, units = "px", quality = 100)
fig_vps
dev.off()


##Vitality per seed lot simplified -- all-years ====
# Define a logCenter function
logCenter <- \(x) (sign(x) * log10(abs(x))) %>% ifelse(is.nan(.), 0, .)
# Load ggplot2 package and process all years data
ggplot(vit %>%
         # Mutate Variety to sentence case and create Batchyear column
         dplyr::mutate(Variety = stringr::str_to_sentence(Variety),
                       Batchyear = paste0(Batch, "-", Year),
                       Field = gsub("S", "K", Field)) %>%
         # Group data by Variety, Field, and Year
         dplyr::group_by(Variety, Field, Year) %>%
         # Scale the 'ValueCorrected' column without centering
         dplyr::mutate(vs = scale(ValueCorrected, scale = F) %>% as.numeric) %>%
         # Group by Variety and Batchyear, and calculate mean for each batch
         dplyr::group_by(Variety, Batchyear) %>%
         dplyr::mutate(meanBatch = mean(vs, na.rm = T)) %>%
         # Group only by Variety
         dplyr::group_by(Variety) %>%
         # Arrange batches by meanBatch
         dplyr::arrange(meanBatch, .by_group = T) %>%
         # Normalize meanBatch values
         dplyr::mutate(meanBatchNorm = (meanBatch - min(meanBatch, na.rm = T))/(max(meanBatch, na.rm = T) - min(meanBatch, na.rm = T))) %>%
         # Mutate Batch as a factor with reversed levels
         dplyr::mutate(Batch = factor(Batchyear, levels = rev(unique(Batchyear)))) %>%
         # Group by Year, Variety, Batch, and Field, and summarise vs statistics
         dplyr::group_by(Year, Variety, Batch, Field) %>%
         dplyr::summarise(vsmn = mean(vs),
                          vshi = max(vs),
                          vslo = min(vs),
                          meanBatchNorm = meanBatchNorm)  %>%
         # Mutate Field to a factor with specific levels
         dplyr::mutate(Field = paste0("Field ", Field) %>%
                         factor(., levels = c("Field M", "Field K", "Field V"))),
       aes(x = Batch, y = vsmn, shape = Field, color = meanBatchNorm)) +
  # Add stripped columns
  ggstats::geom_stripped_cols(aes(color = NULL), odd = pal_material("grey", 30)(30)[10], even = pal_material("grey", 30)(30)[5]) +
  # Add error bars
  geom_errorbar(aes(x = Batch, ymin = vslo, ymax = vshi), size = 1, position = position_dodge(width = 1), width = 0) +
  # Add points
  geom_point(aes(y = vsmn), size = 3, position = position_dodge(width = 1)) +
  # Add dashed horizontal line at y = 0
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 1, color = pal_material("grey")(10)[7]) +
  ylab("Vitality (centered)") +
  xlab("Seedlot") +
  # Set color gradient scale
  scale_color_gradient2(low = pal_effects[3],
                        mid = pal_effects[2],
                        high = pal_effects[1],
                        midpoint = 0.5) +
  # Facet by Variety with free scales
  facet_wrap(~Variety, scales = "free") +
  # Set legend guides
  guides(shape = guide_legend(override.aes = list(size = 7)),
         color = "none") +
  # Set theme elements
  theme(strip.text = element_markdown(size = 20),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_blank(),
        text = element_text("Montserrat Medium", size = 25),
        axis.text.x = element_text("Montserrat Medium", size = 10, angle = 90, vjust = .5, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.spacing.x = unit(3, "lines"))

##Same but only for 1 variety ====
# Load ggplot2 package and subset data for the year 2020 and specific Variety
fig_bes <- ggplot(vit[vit$Year == "2020", ] %>%
                    # Mutate Variety to sentence case
                    dplyr::mutate(Variety = stringr::str_to_sentence(Variety)) %>%
                    # Filter data for a specific Variety
                    dplyr::filter(Variety == "Festien") %>%
                    # Group data by Variety, Field, and Year
                    dplyr::group_by(Variety, Field, Year) %>%
                    # Scale the 'ValueCorrected' column with centering
                    dplyr::mutate(vs = scale(ValueCorrected, scale = T) %>% as.numeric) %>%
                    # Group by Variety and Batch, and calculate mean for each batch
                    dplyr::group_by(Variety, Batch) %>%
                    dplyr::mutate(meanBatch = mean(vs, na.rm = T)) %>%
                    # Group only by Variety
                    dplyr::group_by(Variety) %>%
                    # Arrange batches by meanBatch
                    dplyr::arrange(meanBatch, .by_group = T) %>%
                    # Normalize meanBatch values
                    dplyr::mutate(meanBatchNorm = (meanBatch - min(meanBatch, na.rm = T))/(max(meanBatch, na.rm = T) - min(meanBatch, na.rm = T))) %>%
                    # Mutate Batch as a factor with reversed levels
                    dplyr::mutate(Batch = factor(Batch, levels = rev(unique(Batch)))) %>%
                    # Group by Year, Variety, Batch, and Field, and summarise vs statistics
                    dplyr::group_by(Year, Variety, Batch, Field) %>%
                    dplyr::summarise(vsmn = mean(vs),
                                     vshi = max(vs),
                                     vslo = min(vs),
                                     meanBatchNorm = meanBatchNorm) %>%
                    # Mutate Field to a factor with specific levels
                    dplyr::mutate(Field = paste0("Field ", Field %>% gsub("S", "K", .)) %>% factor(levels = c("Field M", "Field K", "Field V"))),
                  aes(x = Batch, y = vsmn, shape = Field, color = meanBatchNorm)) +
  # Add stripped columns
  ggstats::geom_stripped_cols(aes(color = NULL), odd = pal_material("grey", 30)(30)[10], even = pal_material("grey", 30)(30)[5]) +
  # Add error bars
  geom_errorbar(aes(x = Batch, ymin = vslo, ymax = vshi), size = 1.5, width = 0, position = position_dodge(width = .8)) +
  # Add points
  geom_point(aes(y = vsmn), size = 7, position = position_dodge(width = .8)) +
  # Add dashed horizontal line at y = 0
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 1, color = pal_material("grey")(10)[7]) +
  ylab("Scaled CSA per field") +
  xlab("Seedlot") +
  # Set color gradient scale
  scale_color_gradient2(low = pal_effects[3],
                        mid = pal_effects[2],
                        high = pal_effects[1],
                        midpoint = 0.5) +
  # Facet by Variety with free scales
  facet_wrap(~Variety, scales = "free") +
  # Set legend guides
  guides(shape = guide_legend(override.aes = list(size = 7)),
         color = "none") +
  # Set theme elements
  theme(strip.text = element_markdown(size = 50),
        legend.position = c(.8, .8),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        strip.background = element_rect(fill = "grey"),
        text = element_text(family = "Arial", size = 45),
        axis.title.y = element_markdown(family = "Arial", size = 30),
        axis.text.x = element_text(family = "Arial", size = 35, angle = 90, vjust = .5),
        panel.background = element_rect(fill = "white"),
        panel.spacing.x = unit(3, "lines"))

jpeg(filename = "Figure_1_panelE_scaled_2020.jpg", width = 1200, height = 800, units = "px", quality = 100)
fig_bes
dev.off()

vlist <- vit[vit$Year == "2019", c("Variety", "Batch", "Field", "Blue")] %>% # Subset rows where Year is 2019 and select specified columns
  dplyr::group_by(Field) %>% # Group data by the Field column
  .[!duplicated(.),] %>% # Remove duplicated rows
  dplyr::mutate(Blue = scale(Blue)) %>% # Scale the Blue column
  dplyr::ungroup() %>% # Remove grouping
  dplyr::arrange(Variety, Batch) %>% # Arrange rows by Variety and Batch
  split(., .$Field) %>% # Split data into a list by Field
  lapply(dplyr::select, !c(Field)) # Remove the Field column from each element in the list

for (x in names(vlist)) {colnames(vlist[[x]])[3] <- x} # Loop through list elements and rename the third column to its corresponding Field name
vlist <- cbind(vlist[[1]], vlist[[2]][, c("S")], vlist[[3]][, c("V")]) # Combine the first two list elements and select "S" column from the second element and "V" column from the third element
vlist[, c("M", "S", "V")] %<>% sapply(as.numeric) # Convert columns "M", "S", "V" to numeric

# Commented out line purrr::reduce(., dplyr::left_join, by = c("Variety", "Batch", "Ridge"))

vlist$Variety %<>% stringr::str_to_sentence() # Convert Variety column strings to sentence case
colnames(vlist)[colnames(vlist) == "S"] <- "K" # Rename column "S" to "K"

low <- function(data, mapping, ...) { # Define a function low
  p <- ggplot(data = data, mapping = mapping) + # Create a ggplot object
    geom_point(mapping = mapping, size = 4, show.legend = F, alpha = .6) + # Add points
    geom_smooth(method = "lm", alpha = .25, linewidth = 2, show.legend = F) + # Add smoothed line
    scale_x_continuous(limits = c(-2.4, 2.4), breaks = c(-2, 0, 2)) + # Set x-axis limits and breaks
    scale_y_continuous(limits = c(-2.1, 3.8), breaks = c(-2, 0, 2)) + # Set y-axis limits and breaks
    theme(axis.line = element_line(color = pal_material("grey", 20)(20)[18]), # Set axis line color
          axis.title = element_markdown()) # Set axis title format
  p # Return the plot
}

upper <- function(data) { # Define a function upper
  cor(vlist$V, lmer(V~S+(1|Variety), data = vlist) %>% fitted) # Calculate correlation between V and fitted values of a linear mixed-effects model
}

dia <- function(data, mapping, ...) { # Define a function dia
  ggplot(data = data, mapping = mapping) + # Create a ggplot object
    geom_density(mapping = aes(fill = Variety, y = after_stat(count)), # Add density plot
                 alpha = .2, lwd = 1, bins = 20) + # Set transparency, line width, and number of bins
    scale_x_continuous(limits = c(-2.4, 2.4), breaks = c(-2, 0, 2)) + # Set x-axis limits and breaks
    theme(axis.text.y = element_text(color = "black"), # Set y-axis text color
          axis.line = element_line(color = pal_material("grey", 20)(20)[18])) # Set axis line color
}

vitgp <- ggpairs(vlist, aes(color = Variety), columns = c(3, 4, 5), # Create a pairs plot
                 columnLabels = c("Field M", "Field K", "Field V"),
                 diag = list(continuous = dia),
                 lower = list(continuous = low),
                 upper = list(continuous = wrap(ggally_cor,
                                                title_args = list(size = 11, family = "Arial"),
                                                group_args = list(size = 11, family = "Arial"),
                                                method = "pearson",
                                                title = "All",
                                                digits = 2)))
vitgp[1, 1]

# Define linear mixed-effects models for K vs. M, V vs. M, and V vs. K
lmKM <- lmer(K~M+(1|Variety), data = vlist) %>% performance::r2_nakagawa()
lmVM <- lmer(V~M+(1|Variety), data = vlist) %>% performance::r2_nakagawa()
lmVK <- lmer(V~K+(1|Variety), data = vlist) %>% performance::r2_nakagawa()

# Calculate conditional R-squared values for K vs. M, V vs. M, and V vs. K
KvsM <- paste0("*R<sup>2</sup>* = ", round(lmKM$R2_conditional, 2))
VvsM <- paste0("*R<sup>2</sup>* = ", round(lmVM$R2_conditional, 2))
VvsK <- paste0("*R<sup>2</sup>* = ", round(lmVK$R2_conditional, 2))

# Update the plot to include label and axis information for K vs. M
vitgp[2, 1] <- vitgp[2, 1] +
  ylab("Field K") +
  xlab("Field M") +
  geom_richtext(x = -.5, y = 3, label = KvsM, color = "black",
                fill = NA,
                size = 8,
                family = "Arial",
                label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt"))
# Update the plot to include label and axis information for V vs. M
vitgp[3, 1] <- vitgp[3, 1] +
  ylab("Field V") +
  xlab("Field M") +
  geom_richtext(x = -.5, y = 3, label = VvsM, color = "black",
                fill = NA,
                size = 8,
                family = "Arial",
                label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt"))
# Update the plot to include label and axis information for V vs. K
vitgp[3, 2] <- vitgp[3, 2] +
  ylab("Field V") +
  xlab("Field K") +
  geom_richtext(x = -.5, y = 3, label = VvsK, color = "black",
                fill = NA,
                size = 8,
                family = "Arial",
                label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt"))
#geom_smooth(aes(y = pred, color = NULL), color = "black", se = F, method = "lm", lty = "dashed")
# Update the plot theme and save it as a JPG file
vitgp <- vitgp +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme(text = element_text(family = "Arial", size = 25),
        strip.text = element_markdown(size = 20),
        strip.background = element_rect(fill = "grey"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.spacing = unit(2, "lines"))
# Save the updated plot as a JPG file
ggsave(vitgp, device = "jpg", filename = "correlogram_2020.jpg", dpi = 320, width = 5000, height = 5000, units = "px")

# Reshape the data for plotting
vlist2 <- vlist %>%
  reshape2::melt(id.var = c("Variety", "Batch"),
                 variable.name = "Field",
                 value.name = "CSA") %>%
  dplyr::mutate(Field = paste("Field", Field) %>% factor(., levels = c("Field M", "Field K", "Field V")))
# Summarize data for plotting paired plots
vlist3 <- vlist2 %>%
  dplyr::group_by(Variety, Field) %>%
  dplyr::filter(!is.na(CSA)) %>%
  dplyr::mutate(q25 = quantile(CSA, .25, na.rm = T),
                q75 = quantile(CSA, .75, na.rm = T)) %>%
  dplyr::summarise(CSA = median(CSA, na.rm = T),
                   q25 = q25[1],
                   q75 = q75[1])
# Create paired plots
paired1 <- ggplot(vlist3,
                  aes(x = Field, y = CSA, color = Variety)) +
  PupillometryR::geom_flat_violin(data = vlist2, aes(fill = Variety), position = position_nudge(.2), alpha = .2, show.legend = F) +
  geom_errorbar(aes(x = Field, ymin = q25, ymax = q75),
                width = .1,
                linewidth = 1, position = position_nudge(-.2), show.legend = F) +
  geom_point(aes(y = CSA), size = 5, position = position_nudge(-.2)) +
  geom_point(data = vlist2, alpha = .5, size = 2, position = position_jitter(height = 0, width = .05), show.legend = F) +
  geom_line(data = vlist2 %>%
              dplyr::group_by(Variety, Field) %>%
              dplyr::summarise(CSA = median(CSA, na.rm = T)),
            aes(group = Variety), linewidth = 1, position = position_nudge(-.2), linetype = "dashed", show.legend = F) +
  ylab("Plant vigor (scaled CSA)") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.y = element_text(color = "black"),
        axis.title.y = element_markdown(),
        axis.title.x = element_blank(),
        axis.line = element_line(color = pal_material("grey", 20)(20)[18]),
        text = element_text(family = "Arial", size = 25),
        strip.text = element_markdown(size = 20),
        strip.background = element_rect(fill = "grey"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank())


lay_summary <- c(
  "
  AAB
  AAC
  AAD
  "
)

fig_correlogram_updated <- paired1 + vitgp[2, 1] + vitgp[3, 1] + vitgp[3, 2] + plot_layout(design = lay_summary, guides = "collect")
fig_correlogram_updated
ggsave(fig_correlogram_updated, device = "jpg", filename = "Figure_correlogram_2020_updated.jpg", dpi = 320, width = 5000, height = 5000, units = "px")

# Assigning a character vector with three strings to the variable lay_summary
lay_summary <- c(
  "
  AAB
  AAC
  AAD
  "
)

# Creating a new plot by combining elements from paired1 and vitgp matrices, and adjusting the layout using lay_summary
fig_correlogram_updated <- paired1 + vitgp[2, 1] + vitgp[3, 1] + vitgp[3, 2] + plot_layout(design = lay_summary, guides = "collect")

# Displaying the updated correlogram plot
fig_correlogram_updated

# Saving the updated correlogram plot as a JPG file with specified parameters
ggsave(fig_correlogram_updated, device = "jpg", filename = "Figure_correlogram_2020_updated.jpg", dpi = 320, width = 5000, height = 5000, units = "px")
