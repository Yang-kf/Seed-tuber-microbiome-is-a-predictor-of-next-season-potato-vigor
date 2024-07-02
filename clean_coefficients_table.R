library(pheatmap)  # Load the pheatmap package for creating heatmaps.
library(readxl)  # Load the readxl package for reading Excel files.
library(ggplot2)  # Load the ggplot2 package for creating plots.
library(magrittr)  # Load the magrittr package for the pipe operator.
library(extrafont)  # Load the extrafont package for handling fonts.
library(ggsci)  # Load the ggsci package for additional color palettes.
library(ggtext)  # Load the ggtext package for using text formatting in ggplot2.
library(tibble)  # Load the tibble package for working with data frames.
library(dplyr)  # Load the dplyr package for data manipulation.
library(stringi)  # Load the stringi package for string manipulation.

setwd("~/Postdoc/Song_et_al_2023/coef_determination_figure/")  # Set the working directory.

extrafont::font_import("/Windows/Fonts/", pattern = "Montserrat.*", prompt = F)  # Import Montserrat font.
extrafont::font_import("/Windows/Fonts/", pattern = "Merriweather.*", prompt = F)  # Import Merriweather font.
extrafont::font_import("/Windows/Fonts/", pattern = "Arial.*", prompt = F)  # Import Arial font.
extrafont::loadfonts(device = "win")  # Load the fonts for use in plotting.

# COEFFICIENT OF DETERMINATION R2 ----
cdet <- lapply(1:3, \(x) readxl::read_xlsx("coef_determination_heel.xlsx", sheet = x, col_names = TRUE)) %>%  # Read Excel sheets into a list.
  do.call(rbind, .) %>%  # Combine the list of data frames into a single data frame.
  dplyr::mutate(dplyr::across(2:ncol(.), as.numeric))  # Convert columns 2 through ncol(.) to numeric.

colnames(cdet)[1] <- "Level"  # Rename the first column to "Level".
cdet <- cdet[!stringi::stri_detect(cdet$Level, regex = "Heel_.*_HFE"), ]  # Remove rows where Level matches the specified regex.

cdet$compartment <- cdet$Level %>% strsplit(., "_") %>% sapply('[', 1)  # Extract the compartment from Level column.
cdet$taxon <- cdet$Level %>% strsplit(., "_") %>% sapply('[', 2)  # Extract the taxon from Level column.
firstup <- \(x) {substr(x, 1, 1) %<>% toupper; x}  # Define a function to capitalize the first letter.
cdet$taxLevel <- cdet$Level %>% strsplit(., "_") %>% sapply('[', 3) %>% firstup  # Extract and capitalize the taxonomic level.



# Reshaping the 'cdet' dataframe using the melt function from the reshape2 package, specifying values from columns 2 to 8.
rs.cdet <- reshape2::melt(cdet, value = c(2:8))

# Converting 'taxLevel' column in rs.cdet to a factor with levels specified from the first 7 elements of 'taxLevel' in 'cdet'.
rs.cdet$taxLevel %<>% factor(., levels = cdet$taxLevel[1:7])

# Creating a new column 'perc' in rs.cdet, representing values of 'value' column formatted as percentages.
rs.cdet$perc <- scales::percent(rs.cdet$value, accuracy = 1)

# Splitting the 'variable' column by "_" and extracting the first element, then converting to a factor with levels "Train", "Test", and "Vali".
rs.cdet$set <- strsplit(rs.cdet$variable %>% as.character, "_") %>% sapply(., '[', 1) %>% factor(., levels = c("Train", "Test", "Vali"))

# Splitting the 'variable' column by "_" and extracting the second element to create the 'field' column.
rs.cdet$field <- strsplit(rs.cdet$variable %>% as.character, "_") %>% sapply(., '[', 2)

# Converting 'field' column to a factor with levels "M", "S", and "V".
rs.cdet$field %<>% factor(., levels = c("M", "S", "V"))

# Grouping rs.cdet by 'set' and 'field', then calculating normalized values using min-max scaling and storing them in 'normVal'.
rs.cdet %<>% group_by(set, field) %>% mutate(normVal = (value - min(value)) / (max(value) - min(value)))

# Creating a logical vector 'non.sig' representing nonsignificant cases based on specified conditions.
non.sig <- with(rs.cdet, taxon == "fun" & (
  (variable == "Train_M" & taxLevel == "Phylum")
  | (variable == "Vali_M" & taxLevel == "Phylum")
  | (variable == "Vali_S" & taxLevel %in% c("Family", "Class", "Phylum"))
  | (variable == "Vali_V" & taxLevel %in% c("Species", "Genus", "Family", "Order", "Class", "Phylum"))
))

# Creating a new column 'sig' in rs.cdet, assigning "Nonsignificant" or "Significant" based on the 'non.sig' logical vector.
rs.cdet$sig <- ifelse(non.sig, yes = "Nonsignificant", no = "Significant")

# Creating a vector 'f_order' specifying order for the 'field' variable in the plots.
f_order <- c("Field V", "Field K", "Field M")

# Creating a vector 't_order' specifying order for the 'taxon' variable in the plots.
t_order <- c("bac", "fun", "all")

# Creating the first figure 'fig_det1' using ggplot with specified aesthetics and settings.
fig_det1 <- ggplot(rs.cdet %>%
                     dplyr::filter(taxLevel != "HFE") %>%
                     dplyr::mutate(field = paste0("Field ", field %>% gsub("S", "K", .)) %>% factor(., levels = f_order),
                                   taxon = factor(taxon, levels = t_order),
                                   taxLevel = factor(as.character(taxLevel), levels = rev(levels(taxLevel)))),
                   aes(y = field, x = taxLevel, fill = normVal)) +
  geom_tile(color = "grey20") +
  geom_text(aes(label = round(value, 2), color = sig),
            family = "Arial",
            size = 8) +
  facet_grid(vars(set), vars(taxon), #remove_labels = "all",
             scale = "free_y", space = "free_y", switch = "y", #margins = c("set", "taxon"),
             labeller = labeller(taxon = c("all" = "Bacteria and fungi", "bac" = "Bacteria", "fun" = "Fungi"),
                                 set = c("Train" = "OOB", "Test" = "Within-year<br>testing set", "Vali" = "Across-year<br>testing set"))) +
  scale_fill_gradient2(high = ggsci::pal_material(palette = "teal", n = 40, reverse = T)(15)[15],
                       # low = ggsci::pal_material(palette = "blue", n = 40, reverse = T)(15)[15],
                       low = "white",
                       limits = c(0, 1),
                       midpoint = .5,
                       breaks = c(seq(0, 1, .1)),
                       # labels = scales::percent,
                       name = "Coefficient of<br>determination (*R<sup>2</sup>*)<br>Scaled per set<br>and field") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_color_manual(values = c(pal_material("grey")(10)[c(5, 10)])) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(family = "Arial", size = 30),
        strip.text.x.top = element_text(),
        strip.text.y.left = element_markdown(),
        legend.position = "none",
        strip.placement = "outside",
        panel.background = element_rect(fill = "white"),
        panel.spacing.x = unit(2.5, "lines"),
        panel.spacing.y = unit(1, "lines"))

fig_det1
ggsave("Figure_determination_table.tiff", plot = fig_det1, width = 6600, height = 3000, dpi = 300, units = "px")
ggsave("Figure_determination_table.jpg", plot = fig_det1, width = 6600, height = 3000, dpi = 300, units = "px")

# Create a ggplot object with data from rs.cdet, mutate field, taxon, and taxLevel variables, and convert them to factors with specified levels
fig_det2 <- ggplot(rs.cdet %>%
                     dplyr::mutate(field = paste0("Field ", field %>% gsub("S", "K", .)) %>% factor(., levels = f_order),
                                   taxon = factor(taxon, levels = t_order),
                                   taxLevel = factor(as.character(taxLevel), levels = rev(levels(taxLevel)))),
                   aes(y = field, x = taxLevel, fill = normVal)) +
  geom_tile(color = "grey20") +  # Add tiled rectangles with specified fill color
  geom_text(aes(label = perc, color = sig),  # Add text annotations with specified labels and colors
            family = "Arial",  # Specify font family
            size = 8) +  # Specify text size
  facet_grid(vars(set), vars(taxon),  # Arrange facets in a grid based on set and taxon variables
             scale = "free_y", space = "free", switch = "y",  # Specify scaling and spacing options
             labeller = labeller(taxon = c("all" = "Bacteria and fungi", "bac" = "Bacteria-only", "fun" = "Fungi-only"),  # Label taxon categories
                                 set = c("Train" = "OOB", "Test" = "Within-year<br>testing set", "Vali" = "Across-year<br>testing set"))) +  # Label set categories
  scale_fill_gradient2(high = ggsci::pal_material(palette = "teal", n = 40, reverse = T)(15)[15],  # Define gradient fill color
                       low = "white",  # Define low end color
                       limits = c(0, 1),  # Specify color scale limits
                       midpoint = .5,  # Specify midpoint for color scale
                       breaks = c(seq(0, 1, .1)),  # Specify breaks for color scale
                       name = "Coefficient of<br>determination (*R<sup>2</sup>*)<br>Scaled per set<br>and field") +
  scale_x_discrete(expand = c(0, 0)) + # Modify expansion of x axis
  scale_y_discrete(expand = c(0, 0)) + #Modify expansion of y axis
  scale_color_manual(values = c(pal_material("grey")(10)[c(5, 10)])) + # Add colors for the table
  theme(axis.text.x = element_text(angle = 30, hjust = 1), # Modify theme aesthetics
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(family = "Arial", size = 30),
        strip.text.x.top = element_text(),
        strip.text.y.left = element_markdown(),
        legend.position = "none",
        strip.placement = "outside",
        panel.background = element_rect(fill = "white"))

fig_det2
# Save plot
ggsave("Figure_determination_table_v2.tiff", plot = fig_det2, width = 6600, height = 3000, dpi = 300, units = "px")
ggsave("Figure_determination_table_v2.jpg", plot = fig_det2, width = 6600, height = 3000, dpi = 300, units = "px")

