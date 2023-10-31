#################################### SETUP #####################################
# Loading libraries
library(geiger)
library(phytools)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(cowplot)
library(nlme)
library(dendextend)
library(ComplexHeatmap)

# Clearing the workspace
rm(list = ls())

# Set the current working directory to this script's directory. You will need
# to change this to the path were it's located on your computer
setwd("[path_to_repository]/supplementary/code/")

# Setting a seed for reproductibility
seed <- 666
set.seed(seed)

# Defining a "not in" operator 
`%!in%` <- Negate(`%in%`)

################################ CONFIGURATION #################################
# Phenotype/trait to be evaluated. Either "eusociality" or "parasitoidism"
# Use only one
# dataset <- "eusociality"
dataset <- "parasitoidism"

if (dataset %!in% c("eusociality", "parasitoidism")) {
  stop("Variable `dataset` must be either 'eusociality' or 'parasitoidism'.")
}

# This will be populated with terms that are going to be manually ignored, e.g.
# redundant terms
terms_to_ignore <- c()

# Assigning each species to its group's common name (aka major group)
# Those would be: ants, bees, sawflies and wasps
major_group <- c(
  "Ccin" = "Sawflies", "Oabi" = "Sawflies", "Bkin" = "Wasps", "Nvit" = "Wasps",
  "Tpre" = "Wasps", "Cglo" = "Wasps", "Mdem" = "Wasps", "Cins" = "Wasps", 
  "Fari" = "Wasps", "Dall" = "Wasps", "Vcan" = "Wasps", "Pdom" = "Wasps", 
  "Pfus" = "Wasps", "Pexc" = "Wasps", "Pcan" = "Wasps", "Vpen" = "Wasps", 
  "Vcra" = "Wasps", "Vvel" = "Wasps", "Acom" = "Wasps", "Mrot" = "Bees", 
  "Obic" = "Bees", "Amel" = "Bees", "Alab" = "Bees", "Aflo" = "Bees", 
  "Baff" = "Bees", "Bhun" = "Bees", "Bimp" = "Bees", "Bbif" = "Bees", 
  "Bvan" = "Bees", "Bvos" = "Bees", "Bter" = "Bees", "Bpyr" = "Bees", 
  "Fvar" = "Bees", "Mgen" = "Bees", "Dnov" = "Bees", "Cgig" = "Bees", 
  "Hsal" = "Ants", "Obru" = "Ants", "Dqua" = "Ants", "Cflo" = "Ants", 
  "Fexs" = "Ants", "Nful" = "Ants", "Aech" = "Ants", "Acep" = "Ants", 
  "Acol" = "Ants", "Tsep" = "Ants", "Tcor" = "Ants", "Tzet" = "Ants", 
  "Ccos" = "Ants", "Waur" = "Ants", "Mpha" = "Ants", "Sinv" = "Ants", 
  "Veme" = "Ants", "Pbar" = "Ants", "Lhum" = "Ants", "Pgra" = "Ants", 
  "Obir" = "Ants", "Ebur" = "Ants", "Aros" = "Sawflies", "Dsim" = "Sawflies", 
  "Nfab" = "Sawflies", "Nvir" = "Sawflies", "Nlec" = "Sawflies", 
  "Npin" = "Sawflies", "Ccal" = "Bees", "Emex" = "Bees", "Nmel" = "Bees"
)

# Cutoff for the FDR-corrected ANOVA p-values
q_value_cutoff = 0.05

# Loading table with IPR term counts
term_counts <- read.table("../data/IPR_term_counts.tsv")

# Loading translation table between IPR ID and description
term_descriptions <- read.delim("../data/IPR_term_description.tsv")

# Ultrametric phylogenetic tree
tree <- read.tree("../data/ort80_hymenoptera_67species_100mil_mcmctree.nwk")
# Rotate Ampulex compressa (Acom) to be more distant from other wasps
tree <- ape::rotate(tree, 88)
# Rotate the sawflies to be closer to the solitary/parasitoid wasps
tree <- ape::rotate(tree, 68)

# Reordering the major groups vector to match the oder in the tree
tmp_groups <- c()
for (species in tree$tip.label) {
  tmp_groups[species] = major_group[[species]]
}
major_group <- tmp_groups

# Reporting unique annotated IPR terms
unique_terms <- term_counts %>% ncol()
print(
  paste0(
    "A total of ",
    unique_terms,
    " unique IPR terms were found annotated across all proteomes"
  )
)

# Reporting total annotated gene-IPR pairs
total_annotated <- sum(colSums(term_counts, na.rm = TRUE))
print(
  paste0(
    "A total of ",
    total_annotated,
    " gene-IPR pairs were found annotated across all proteomes"
  )
)

# Loading the metadata and defining some variables based on the phenotype
# being evaluated
if (dataset == "eusociality") {
  metadata <- read.table("../data/species_metadata_eusociality.tsv", sep = "\t", header = TRUE) 
  
  factor <- "IS_EUSOCIAL"
  
  feature_labels <- c("Solitary", "Eusocial")
  
  dataset_label <- "Eusociality level"
  
  fig_name <- "Figure_3_IPR_eusociality.pdf"
  
  # For this analysis we remove the intermediately eusocial species from the
  # dataset. This includes removing them from the tree and major groups vector
  spp2remove = c(metadata$ABBREV[metadata$SOCIALITY_LEVEL == 2])
  major_group <- major_group[names(major_group) %!in% spp2remove]
  tree <- drop.tip(tree, tree$tip.label[match(spp2remove, tree$tip.label)])
  
  terms_to_ignore <- c("IPR036658")
} else if (dataset == "parasitoidism") {
  metadata <- read.table("../data/species_metadata_parasitoidism.tsv", sep = "\t", header = TRUE) 
  
  factor <- "PARASITOIDISM"
  
  dataset_label <- "Parasite/parasitoid"
  
  feature_labels <- c("No parasitism", "Parasite/parasitoid")
  
  fig_name <- "Figure_4_IPR_parasitoidism.pdf"
  
  terms_to_ignore <- c(
    "IPR008753", "IPR042089", "IPR018497", "IPR013899", "IPR036063"
  )
} else {
  stop("Unrecognized dataset option.")
}

#################################### ANOVA #####################################
# Coverting the phenotype classes to a factor for the analyses later
classes <- as.factor(metadata[,factor][match(tree$tip.label, metadata$ABBREV)])
names(classes) <- metadata$ABBREV[match(tree$tip.label, metadata$ABBREV)]

# Limiting the IPR counts data frame to the species under analysis
annotation_df <- term_counts[tree$tip.label, colnames(term_counts)]

# Stores the p-value for each ANOVA result for each IPR term
p_values <- vector()

# For the sake of pretty printing
i = 1

print("Computing phylo-ANOVA...")

# ANOVA com correção filogenética
for (annotation_ID in unique(colnames(annotation_df)))  {
  # Current IPR under analysis - current iteration number
  print(paste(annotation_ID, " - ", i))
  i <- i + 1
  
  current_ipr_count <- as.vector(annotation_df[,annotation_ID])
  names(current_ipr_count) <- metadata$ABBREV[
    match(rownames(annotation_df), metadata$ABBREV)
  ]
  
  # Skip IPR terms with either 0 count or no variation
  cv <- sd(current_ipr_count) / mean(current_ipr_count)
  if ((sum(current_ipr_count) <= 0) || (cv <= 0)) {
    next
  }
  
  # Converting the counts to a named list
  current_ipr_count_wide <- cbind.data.frame(current_ipr_count, classes)
  
  # Skip IPR terms that only occur in <= 3 different count classes and with
  # a count of < 3 in any of the less frequent classes
  if (length(unique(current_ipr_count)) <= 3) { 
    if (any(table(current_ipr_count_wide[, 1]) < 3)) {
      next
    } 
  }
  
  spp <- rownames(current_ipr_count_wide)
  # Covariance matrix with Brownian motion and phylogenetic correction
  corBM <- corBrownian(phy = tree, form = ~spp)
  # Fitting a genearalised least squares model with our covariance matrix
  ancova <- gls(current_ipr_count ~ classes, data = current_ipr_count_wide, correlation = corBM)
  # ANOVA of the fitted GLS
  anova_results <- anova(ancova)
  # Saving the p-value
  p_values[annotation_ID] <- anova_results$`p-value`[2]
}

# False discovery rate correction
print("Performing multiple-hypothesis-testing correction (FDR)...")
q_values  <- p.adjust(p_values, method = "BH")

# Sub-setting only the significant IPR terms based on FDR-corrected p-values
significant_IDs <- names(q_values[q_values < q_value_cutoff & !is.na(q_values)])
significant_IDs <- significant_IDs[significant_IDs %!in% terms_to_ignore]

# Stop early if no significant ids were found
if (length(significant_IDs) == 0) {
  stop("No significant IDs were found.")
}

# Sub-setting the annotation term dataframe for only significant terms
annotation_df <- annotation_df %>% select(all_of(significant_IDs))

print("Phylo-ANOVA completed successfuly!")

################################## PLOTTING ####################################
print("Generating plots...")

# Boolean table indicating which species have the phenotype
species_to_phenotype = metadata[
  metadata[, "ABBREV"] %in% tree$tip.label, 
  c("ABBREV", factor)
] 
colnames(species_to_phenotype) = c("short_name", "has_phenotype")
# Fixing the ordering to match the tree
species_to_phenotype <- species_to_phenotype[
  match(tree$tip.label, species_to_phenotype$short_name),
]

# Extracting column names for the significant terms dataframe
colnames_for_annotation_df <- term_descriptions %>% 
  filter(ipr %in% significant_IDs) %>% 
  mutate(term_and_description = paste0(ipr, " - ", description)) %>% 
  arrange(ipr) %>% 
  select(term_and_description) %>% 
  pull(1)

# Assigning the new column names to the significant terms dataframe. Here we
# also add text wrapping for long names so it looks better in the plot
colnames(annotation_df) <- stringr::str_wrap(
  colnames_for_annotation_df, width = 55
)

# Converting the tree to a dendogram object suitable for plotting the heatmaps
tree_as_dendrogram <- as.dendrogram(as.hclust.phylo(tree), edge.root = TRUE)

# Setting the colour palette for the IPR counts
if (dataset == "parasitoidism") {
  counts_palette <- c("white", rev(viridis::magma(40)))  
} else {
  # Here the eusociality palette is manually set to only 3 colors because it
  # was verified in a previous run that all significant IPR terms had either
  # zero, one or two copies
  counts_palette <- c("white", "#FECD90FF", "#FB845FFF")
}

# Colour palette for the presence/absence of the phenotype
phenotype_colours <- c("FALSE" = "black", "TRUE" = "red")

# Convert the phenotype information to a long matrix for the heatmap
phenotype_map <- t(rbind(setNames(
  as.character(species_to_phenotype$has_phenotype), 
  species_to_phenotype$short_name
)))

# Major groups colours
groups_colors <- c(
  "Sawflies" = "#F16E1E", "Wasps" = "#DC267F", 
  "Bees" = "#4FD633", "Ants" = "#785EF0"
)

# Converting the tree to a hierarchical cluster for the heatmap
final_species_tree <- as.hclust(set(tree_as_dendrogram, "branches_lwd", 2))

# Heatmap visual parameters
default_font <- "sans serif"
monospaced_font <- "mono"
legend_params <- list(
  "title_gp" = gpar(fontfamily = default_font),
  "labels_gp" = gpar(fontfamily = default_font)
)
if (dataset == "parasitoidism") {
  hm1_legend_params <- legend_params
} else {
  hm1_legend_params <- list(
    "at" = c(0, 1, 2), 
    "labels" = c("0", "1", "2"),
    "color_bar" = "discrete", 
    "border" = "black",
    "title_gp" = gpar(fontfamily = default_font),
    "labels_gp" = gpar(fontfamily = default_font)
  )
}
hm2_legend_params <- list(
  "at" = c("FALSE", "TRUE"), 
  "labels" = feature_labels,
  "title_gp" = gpar(fontfamily = default_font),
  "labels_gp" = gpar(fontfamily = default_font)
)

if (dataset == "parasitoidism") {
  width_factor <- 1.2
} else {
  width_factor <- 2.1
}

# Heatmap of the IPR terms significantly associated with the phenotype
ht1 = Heatmap(
  as.matrix(annotation_df),
  clustering_method_columns = c("ward.D2"),
  cluster_rows = final_species_tree, 
  name = "IPR Count", 
  col = counts_palette,
  column_names_rot = 60,
  column_names_gp = gpar(fontfamily = default_font),
  column_title_gp = gpar(fontfamily = default_font),
  row_dend_width = unit(5, "cm"), 
  row_dend_gp = gpar(lwd = 1.2, edge.root = TRUE), 
  row_names_side = "left", 
  row_names_gp = gpar(fontfamily = monospaced_font, fontface = "bold"),
  row_title_gp = gpar(fontfamily = default_font),
  heatmap_legend_param = hm1_legend_params,
)

# Presence/absence matrix of the phenotype
ht2 = Heatmap(
  phenotype_map, 
  cluster_rows = final_species_tree, 
  name = dataset_label, 
  col = phenotype_colours, 
  column_dend_height = unit(2, "cm"),
  column_dend_gp = gpar(lwd = 1.2), 
  column_labels = "", 
  column_names_gp = gpar(fontfamily = default_font),
  column_title_gp = gpar(fontfamily = default_font),
  height = 1, 
  heatmap_legend_param = hm2_legend_params, 
  row_dend_gp = gpar(lwd = 1.2, edge.root = TRUE), 
  row_dend_width = unit(5, "cm"),
  row_names_side = "left", 
  row_names_gp = gpar(fontfamily = monospaced_font, fontface = "bold"),
  row_title_gp = gpar(fontfamily = default_font),
  row_title = "Species"
)

# Major group matrix
ht3 = Heatmap(
  t(rbind(major_group)), 
  cluster_rows = final_species_tree, 
  name = "Major group", 
  col = groups_colors, 
  column_dend_height = unit(2, "cm"), 
  column_names_gp = gpar(fontfamily = default_font),
  column_title_gp = gpar(fontfamily = default_font),
  column_dend_gp = gpar(lwd = 1.2), 
  column_labels = "", 
  height = 1, 
  row_dend_gp = gpar(lwd = 1.2, edge.root = TRUE), 
  row_dend_width = unit(5, "cm"),
  row_names_side = "left", 
  row_names_gp = gpar(fontfamily = monospaced_font, fontface = "bold"),
  row_title_gp = gpar(fontfamily = default_font),
  heatmap_legend_param = legend_params,
  row_title = "Species"
)

# Concatenation of the three heatmaps
ht_list = ht3 + ht2 + ht1

# Saving the final figure
cairo_pdf(
  paste0("../../figures/", fig_name), 
  width = width_factor * ncol(annotation_df),
  height = 20
)
draw(
  ht_list, 
  padding = unit(c(50, 2, 2, 2), "mm"), 
  heatmap_legend_side = "right"
)
dev.off()
print("Done! Ran successfuly!")
