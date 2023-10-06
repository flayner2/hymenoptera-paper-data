#################################### SETUP #####################################
# Loading libraries
library(treeio)
library(ggtree)
library(RRphylo)
library(ape)
library(phytools)
library(Cairo)
library(ggplot2)
library(dplyr)
library(ggnewscale)

# Clearing the workspace
rm(list = ls())

# Set the current working directory to this script's directory. You will need
# to change this to the path were it's located on your computer
setwd("~/Documents/LAB/eusociality/work/masters_paper/doc/hymenoptera-paper-data/supplementary/code/")

# Setting a seed for reproductibility
seed <- 666
set.seed(seed)

print("Building the backbone tree...")

############################# BUILDING THE TREE ################################
# "Backbone" tree from Peters et al. 2017
backbone_treedata <- read.mcmctree("../data/stitched_tree/dated_tree_aa_inde_2_used_in_Fig1.tre")
backbone_tree <- backbone_treedata@phylo

# Bombus tree from Santos Junior et al. 2022
bombus_treedata <- read.beast("../data/stitched_tree/bombus_UCLN_BEAST_MODEL_last_fossils.tre")
bombus_tree <- bombus_treedata@phylo

# Mapping the species from Santos Junior et al. 2022 to the backbone tree
bombus_data <- data.frame(
  bind = c("Bombus_impatiens_voucher_SC060", "Bombus_bifarius_voucher_SC208", 
           "Bombus_vancouverensis_nearticus", "Bombus_vosnesenskii_voucher_SC112", 
           "Bombus_terrestris_voucher_SC003", "Bombus_pyrosoma_B04",
           "Bombus_huntii_voucher_SC151", "Bombus_affinis_voucher_SC167"), 
  reference = c("Bombus_rupestris", "Bombus_impatiens_voucher_SC060", 
                "Bombus_bifarius_voucher_SC208", 
                "Bombus_bifarius_voucher_SC208-Bombus_vancouverensis_nearticus", 
                "Bombus_impatiens_voucher_SC060-Bombus_vancouverensis_nearticus",
                "Bombus_terrestris_voucher_SC003-Bombus_vancouverensis_nearticus",
                "Bombus_vosnesenskii_voucher_SC112", "Bombus_terrestris_voucher_SC003"), 
  poly = rep(FALSE, times = 8)
)

# Tree including the species from Santos Junior et al. 2022
backbone_part_1 <- tree.merger(
  backbone = backbone_tree, 
  source.tree = bombus_tree, 
  data = bombus_data, 
  plot = FALSE
)

# Apidae tree from Bossert et al. 2019
apidae_tree <- read.newick("../data/stitched_tree/Phylobayes_80p.tre")

# Mapping the species from Bossert et al. 2019 to the backbone tree
apidae_data <- data.frame(
  bind = c("Apis_laboriosa", "GEN_apis_florea", 
           "GEN_ceratina_calcarata", "GEN_eufriesea_mexicana", 
           "GEN_megachile_rotundata"), 
  reference = c("Apis_mellifera", "Apis_mellifera-Apis_laboriosa", 
                "Ceratina_chalybea", "Euglossa_dilemma", 
                "Megachile_willughbiella"), 
  poly = rep(FALSE, times = 5)
)

# Tree including the species from Bossert et al. 2019
backbone_part_2 <- tree.merger(
  backbone = backbone_part_1, 
  source.tree = apidae_tree, 
  data = apidae_data, 
  plot = FALSE
)

# Hymenoptera tree from Branstetter et al. 2017
branstetter_treedata <- read.beast("../data/stitched_tree/hym-187t-f75-beast-rand50.tre")
branstetter_tree <- branstetter_treedata@phylo

# Removing the "-" character from the names of the species in the original tree
# from Branstetter et al. 2017 since it breaks the `tree.merger` function
branstetter_tree$tip.label <- sapply(
  branstetter_tree$tip.label, function(x) sub("-", "_", x)
)

# Mapping the species from Branstetter et al. 2017 to the backbone tree
branstetter_data <- data.frame(
  bind = c(
    "TENTHREDINIDAE_Athalia_rosae_Genome", "FORMICIDAE_Atta_cephalotes_Genome",
    "Atta_colombica", "FORMICIDAE_Wasmannia_auropunctata_Genome",
    "FORMICIDAE_Monomorium_pharaonis_Genome", "FORMICIDAE_Solenopsis_invicta_Genome",
    "FORMICIDAE_Vollenhovia_emeryi_Genome", "FORMICIDAE_Pogonomyrmex_barbatus_Genome",
    "CEPHIDAE_Cephus_cinctus_Genome", "HALICTIDAE_Dufourea_novaeangliae_Genome",
    "BRACONIDAE_Microplitis_demolitor_Genome",
    "BRACONIDAE_Fopius_arisanus_Genome", "FORMICIDAE_Linepithema_humile_Genome",
    "Neodiprion_fabricii", "Neodiprion_virginianus",
    "DIPRIONIDAE_Neodiprion_lecontei_Genome", "Neodiprion_pinetum",
    "TRICHOGRAMMATIDAE_Trichogramma_pretiosum_Genome"
  ), 
  reference = c(
    "Tenthredo_koehleri-Nematus_ribesii", "Acromyrmex_echinatior",
    "FORMICIDAE_Atta_cephalotes_Genome", "Acromyrmex_echinatior-FORMICIDAE_Atta_cephalotes_Genome",
    "FORMICIDAE_Wasmannia_auropunctata_Genome-FORMICIDAE_Atta_cephalotes_Genome",
    "FORMICIDAE_Monomorium_pharaonis_Genome",
    "FORMICIDAE_Solenopsis_invicta_Genome-FORMICIDAE_Atta_cephalotes_Genome",
    "FORMICIDAE_Vollenhovia_emeryi_Genome-FORMICIDAE_Atta_cephalotes_Genome",
    "Cephus_spinipes", "Dufourea_dentiventris", 
    "Cotesia_vestalis-Dacnusa_sibirica",
    "BRACONIDAE_Microplitis_demolitor_Genome-Dacnusa_sibirica",
    "Camponotus_floridanus-FORMICIDAE_Atta_cephalotes_Genome",
    "Diprion_pini", "Neodiprion_fabricii",
    "Neodiprion_virginianus", "DIPRIONIDAE_Neodiprion_lecontei_Genome",
    "Gonatocerus_morrilli-Nasonia_vitripennis"
  ), 
  poly = rep(FALSE, times = 18)
)

# Tree including the species from Branstetter et al. 2017
backbone_part_3 <- tree.merger(
  backbone = backbone_part_2, 
  source.tree = branstetter_tree, 
  data = branstetter_data, 
  plot = FALSE
)

# Vespidae tree from Silva 2021
silva_vespidae_tree <- read.newick("../data/stitched_tree/Vespidae all data 2.tre")

# Mapping the Vespidae species from Silva 2021 to the backbone tree
silva_vespidae_data <- data.frame(
  bind = c(
    "Polistes_fuscatus",
    "Polistes_exclamans",
    "Polistes_canadensis"
  ),
  reference = c(
    "Polistes_dominula",
    "Polistes_fuscatus",
    "Polistes_exclamans"
  ),
  poly = rep(FALSE, times = 3)
)

# Tree including the Vespidae species from Silva 2021
backbone_part_4 <- tree.merger(
  backbone = backbone_part_3, 
  source.tree = silva_vespidae_tree, 
  data = silva_vespidae_data, 
  plot = FALSE
)

# Halictidae tree from Silva 2021
silva_halictidae_tree <- read.newick("../data/stitched_tree/Halictidae all data.tre")

# Mapping the Halictidae from Silva 2021 to the backbone tree
silva_halictidae_data <- data.frame(
  bind = c(
    "Nomia_melanderi",
    "Megalopta_genalis"
  ),
  reference = c(
    "Nomia_diversipes",
    "Nomia_diversipes-Lasioglossum_xanthopus"
  ),
  poly = rep(FALSE, times = 2)
)

# Tree including the Halictidae species from Silva 2021
backbone_part_5 <- tree.merger(
  backbone = backbone_part_4, 
  source.tree = silva_halictidae_tree, 
  data = silva_halictidae_data, 
  plot = FALSE
)

# Formicidae tree from Romiguier et al. 2022
romiguier_treedata <- read.mcmctree("../data/stitched_tree/FigTree.tre")
romiguier_tree <- romiguier_treedata@phylo

# Mapping the species from Romiguier et al. 2022 to the backbone tree
romiguier_data <- data.frame(
  bind = c(
    "Formica_sanguinea",
    "Formica_exsecta",
    "Nylanderia_fulva",
    "Odontomachus_haematodus",
    "Odontomachus_brunneus",
    "Pseudomyrmex_pallidus",
    "Pseudomyrmex_gracilis",
    "Ooceraea_biroi",
    "Eciton_burchellii"
  ),
  reference = c(
    "Camponotus_floridanus",
    "Formica_sanguinea",
    "Formica_exsecta-Camponotus_floridanus",
    "Harpegnathos_saltator",
    "Odontomachus_haematodus",
    "FORMICIDAE_Linepithema_humile_Genome",
    "Pseudomyrmex_pallidus",
    "Pseudomyrmex_gracilis-FORMICIDAE_Atta_cephalotes_Genome",
    "Ooceraea_biroi"
  ),
  poly = rep(FALSE, times = 9)
)

# Tree including the species frm Romiguier et al. 2022
backbone_part_6 <- tree.merger(
  backbone = backbone_part_5, 
  source.tree = romiguier_tree, 
  data = romiguier_data, 
  plot = FALSE
)

# Removing the species used only as reference points
backbone_part_6 <- drop.tip(backbone_part_6, c(
  "Formica_sanguinea", 
  "Odontomachus_haematodus",
  "Pseudomyrmex_pallidus")
)

# Here we include species that were manually inserted into the backbone tree
# based on references from the literature (see the Supplementary Table) but
# that didn't have an associated tree file available for download
manual_data <- data.frame(
  bind = c(
    "Belonocnema_kinseiy",
    "Chelonus_insularis", "Cotesia_glomerata",
    "Venturia_canescens", "Colletes_gigas",
    "Frieseomelitta_varia", "Osmia_bicornis_bicornis",
    "Vespa_velutina", "Vespula_pensylvanica",
    "Diprion_similis", "Diachasma_alloeum",
    "Trachymyrmex_septentrionalis", "Trachymyrmex_cornetzi",
    "Trachymyrmex_zeteki", "Cyphomyrmex_costatus",
    "Dinoponera_quadriceps", "Drosophila_melanogaster"
  ),
  reference = c(
    "Andricus_quercuscalicis",
    "BRACONIDAE_Microplitis_demolitor_Genome-Dacnusa_sibirica", "Cotesia_vestalis",
    "Hyposoter_didymator", "Colletes_cunicularius",
    "Tetragonula_carbonaria", "Osmia_cornuta",
    "Vespa_crabro", "Vespula_germanica",
    "Diprion_pini", "BRACONIDAE_Fopius_arisanus_Genome",
    "Acromyrmex_echinatior-FORMICIDAE_Atta_cephalotes_Genome", 
    "Trachymyrmex_septentrionalis-FORMICIDAE_Atta_cephalotes_Genome",
    "Trachymyrmex_cornetzi-FORMICIDAE_Atta_cephalotes_Genome", 
    "Trachymyrmex_zeteki-FORMICIDAE_Atta_cephalotes_Genome",
    "Harpegnathos_saltator-Odontomachus_brunneus", 
    "Pseudomallada_prasinus-Xanthostigma_xanthostigma"
  ),
  poly = rep(FALSE, times = 17)
)

# Final tree
final_tree <- tree.merger(
  backbone = backbone_part_6, 
  data = manual_data, 
  plot = FALSE
)

# Renaming the species in the final tree to the 4-letter code used in this work
names_map <- c(
  # Species present in the original tree from Peters et al. 2017
  "Acromyrmex_echinatior" = "Aech",
  "Ampulex_compressa" = "Acom",
  "Apis_mellifera" = "Amel",
  "Camponotus_floridanus" = "Cflo",
  "Harpegnathos_saltator" = "Hsal",
  "Nasonia_vitripennis" = "Nvit",
  "Orussus_abietinus" = "Oabi",
  "Polistes_dominula" = "Pdom",
  "Tribolium_castaneum" = "Tcas",
  "Vespa_crabro" = "Vcra",
  # The inserted species
  "Bombus_bifarius_voucher_SC208" = "Bbif",
  "Bombus_impatiens_voucher_SC060" = "Bimp",
  "Bombus_pyrosoma_B04" = "Bpyr",
  "Bombus_terrestris_voucher_SC003" = "Bter",
  "Bombus_vosnesenskii_voucher_SC112" = "Bvos",
  "Bombus_vancouverensis_nearticus" = "Bvan",
  "Bombus_huntii_voucher_SC151" = "Bhun", 
  "Bombus_affinis_voucher_SC167" = "Baff",
  "GEN_apis_florea" = "Aflo",
  "GEN_ceratina_calcarata" = "Ccal",
  "GEN_eufriesea_mexicana" = "Emex",
  "GEN_megachile_rotundata" = "Mrot",
  "Apis_laboriosa" = "Alab",
  "TENTHREDINIDAE_Athalia_rosae_Genome" = "Aros",
  "FORMICIDAE_Atta_cephalotes_Genome" = "Acep",
  "FORMICIDAE_Wasmannia_auropunctata_Genome" = "Waur",
  "FORMICIDAE_Monomorium_pharaonis_Genome" = "Mpha",
  "FORMICIDAE_Solenopsis_invicta_Genome" = "Sinv",
  "FORMICIDAE_Vollenhovia_emeryi_Genome" = "Veme",
  "FORMICIDAE_Pogonomyrmex_barbatus_Genome" = "Pbar",
  "CEPHIDAE_Cephus_cinctus_Genome" = "Ccin",
  "HALICTIDAE_Dufourea_novaeangliae_Genome" = "Dnov",
  "BRACONIDAE_Microplitis_demolitor_Genome" = "Mdem",
  "BRACONIDAE_Fopius_arisanus_Genome" = "Fari",
  "FORMICIDAE_Linepithema_humile_Genome" = "Lhum",
  "DIPRIONIDAE_Neodiprion_lecontei_Genome" = "Nlec",
  "TRICHOGRAMMATIDAE_Trichogramma_pretiosum_Genome" = "Tpre",
  "Atta_colombica" = "Acol",
  "Neodiprion_virginianus" = "Nvir",
  "Neodiprion_fabricii" = "Nfab",
  "Neodiprion_pinetum" = "Npin",
  "Polistes_fuscatus" = "Pfus",
  "Polistes_exclamans" = "Pexc",
  "Nomia_melanderi" = "Nmel",
  "Megalopta_genalis" = "Mgen",
  "Polistes_canadensis" = "Pcan",
  "Formica_exsecta" = "Fexs",
  "Nylanderia_fulva" = "Nful",
  "Odontomachus_brunneus" = "Obru",
  "Pseudomyrmex_gracilis" = "Pgra",
  "Ooceraea_biroi" = "Obir",
  "Eciton_burchellii" = "Ebur",
  "Belonocnema_kinseiy" = "Bkin",
  "Chelonus_insularis" = "Cins",
  "Cotesia_glomerata" = "Cglo",
  "Venturia_canescens" = "Vcan",
  "Colletes_gigas" = "Cgig",
  "Frieseomelitta_varia" = "Fvar",
  "Osmia_bicornis_bicornis" = "Obic",
  "Vespa_velutina" = "Vvel",
  "Vespula_pensylvanica" = "Vpen",
  "Diprion_similis" = "Dsim",
  "Diachasma_alloeum" = "Dall",
  "Trachymyrmex_septentrionalis" = "Tsep",
  "Trachymyrmex_cornetzi" = "Tcor",
  "Trachymyrmex_zeteki" = "Tzet",
  "Cyphomyrmex_costatus" = "Ccos",
  "Dinoponera_quadriceps" = "Dqua",
  "Drosophila_melanogaster" = "Dmel"
)

# Renaming the tips and removing the ones that should not be in the final tree
for (tip_name in final_tree$tip.label) {
  if (tip_name %in% names(names_map)) {
    final_tree$tip.label[which(final_tree$tip.label == tip_name)] <- names_map[[tip_name]]
  } else {
    final_tree <- drop.tip(final_tree, tip_name)
  }
}

print("The tree was constructed successfuly")

# Saving the tree file. The saved tree is the one constructed here, not the 
# `main_tree` below, which is the tree from the MCMCTree results
print("Saving tree in NEWICK format...")
# Remove the outgroups
tree_to_save <- drop.tip(final_tree, c("Dmel", "Tcas"))
# Erasing branch length data as it isn't useful for us
tree_to_save$edge.length <- NULL
# Saving the final tree without the outgroups and branch lengths
write.tree(tree_to_save, "../data/stitched_tree.nwk")

############################# PLOTTING THE TREE ################################
print("Generating the tree plots...")

#### The backbone tree
# Renaming the tips to follow our convention
for (tip_name in backbone_tree$tip.label) {
  split_name <- strsplit(tip_name, "_")
  if (split_name[[1]][2] == "sp") {
    fixed_name <- paste0(substr(split_name[[1]][1], 1, 2), substr(split_name[[1]][2], 1, 2))
  }
  else {
    fixed_name <- paste0(substr(split_name[[1]][1], 1, 1), substr(split_name[[1]][2], 1, 3))
  }
  backbone_tree$tip.label[which(backbone_tree$tip.label == tip_name)] <- fixed_name
}

# Graphical parameters for plotting
tree_font_size <- 6
tree_font_family <- "mono"
tree_clade_font_family <- "sans serif"
tree_pdf_width <- 20
tree_pdf_height <- 46
tree_font_size_legend <- 18
tree_font_size_legend_title <- 20
tree_legend_margin <- margin(0, 0, 0, 0)
tree_legend_box_margin <- margin(r = 10)

# Classifying the tips as either "removed" or "anchor" (reference) points
backbone_removed_nodes <- data.frame(
  node = c(
    304, 119, 312, 118, 116, 115, 295, 111, 110, 141, 144, 325, 327, 109, 282,
    272, 252, 70, 69, 246, 233, 49, 53, 48, 221, 154, 216, 212, 211, 31, 200, 11,
    10, 13, 14, 15, 9, 196, 197, 7, 186, 8, 2, 331, 169, 339, 162, 161, 335, 120, 
    309, 117, 114, 319, 142, 143, 329, 230, 30, 203, 12, 6, 1, 168, 342, 344
  ),
  status = c(rep("removed", times = 49), rep("anchor", times = 17))
)

# Plotting
cairo_pdf(
  "../figures/Supplementary_Figure_4_tree_backbone_removed_anchor.pdf", 
  width = tree_pdf_width, 
  height = tree_pdf_height
)
ggtree(backbone_tree, linewidth = 2, branch.length = "none") + 
  geom_tiplab(
    size = tree_font_size, 
    family = tree_font_family,
    fontface = "bold"
  ) + 
  geom_highlight(
    data = backbone_removed_nodes, 
    mapping = aes(node = node, fill = status, extend = 1.3)
  ) +
  geom_rootedge(rootedge = 0.5, linewidth = 2) +
  scale_fill_manual(
    name = "",
    breaks = c("removed", "anchor"),
    labels = c("Removed", "Reference"),
    values = c("red", "blue")) +
  theme(legend.text = element_text(size = tree_font_size_legend),
        legend.title = element_text(size = tree_font_size_legend_title),
        legend.margin = tree_legend_margin,
        legend.box.margin = tree_legend_box_margin) +
  scale_y_continuous(limits = c(-1, NA))
dev.off()

#### Intermediate tree 
# This is the tree with the non-reference nodes removed, indicating the nodes
# used only for reference that will afterwards be removed as well
intermediate_to_remove <- c(
  126, 125, 124, 123, 122, 121, 119, 134, 133, 132, 131, 135, 118, 116, 115, 113, 
  112, 111, 110, 141, 144, 147, 146, 145, 150, 149, 148, 109, 108, 107, 106, 105, 
  103, 102, 104, 99, 98, 97, 96, 95, 94, 101, 100, 93, 92, 91, 90, 89, 88, 87, 
  86, 85, 84, 83, 82, 81, 79, 78, 80, 77, 76, 75, 74, 73, 72, 70, 69, 68, 67, 64, 
  63, 62, 66, 65, 61, 60, 59, 58, 56, 55, 54, 57, 49, 53, 48, 45, 44, 43, 42, 47, 
  46, 154, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 22, 21, 23, 11, 10, 13, 
  14, 15, 9, 17, 16, 19, 18, 20, 7, 5, 4, 8, 2, 1, 157, 156, 155, 158, 169, 165, 
  164, 163, 162, 161, 160, 159
)
intermediate_tree <- drop.tip(backbone_tree, intermediate_to_remove)

# The reference nodes that will be removed later
intermediate_removed_nodes <- data.frame(
  node = c(
    64, 63, 18, 17, 59, 12, 46, 2, 72, 77, 35, 39
  ),
  status = c(rep("removed", times = 12))
)

# Plotting
cairo_pdf(
  "../figures/Supplementary_Figure_5_tree_intermediate_anchor.pdf", 
  width = tree_pdf_width, 
  height = 16
)
ggtree(intermediate_tree, linewidth = 2, branch.length = "none") + 
  geom_tiplab(
    size = tree_font_size, 
    family = tree_font_family,
    fontface = "bold"
  ) + 
  geom_highlight(
    data = intermediate_removed_nodes, 
    mapping = aes(node = node, fill = status, extend = 0.7)
  ) +
  geom_rootedge(rootedge = 0.5, linewidth = 2) +
  scale_fill_manual(
    name = "",
    breaks = c("removed"),
    labels = c("Removed"),
    values = c("red")) +
  theme(legend.text = element_text(size = tree_font_size_legend),
        legend.title = element_text(size = tree_font_size_legend_title),
        legend.margin = tree_legend_margin,
        legend.box.margin = tree_legend_box_margin) +
  scale_y_continuous(limits = c(-1, NA))
dev.off()

#### Final tree
# First plot contains a simple representation of the tree and a highlight of
# the included species. This tree also includes the outgroup species

# The included species not present in the original backbone tree
highlighted_nodes <- data.frame(
  node = c(
    # ants
    128, 49, 50, 51, 52, 53, 129, 56, 57, 44, 41, 42, 45, 130, 131,
    # bees
    100, 25, 26, 23, 22, 94, 108,
    # Drosophila melanogaster
    69,
    # Sawflies
    1, 132,
    # Wasps
    86, 18, 16, 5, 3, 78
  ),
  group = c(
    rep("ants", times = 15),
    rep("bees", times = 7),
    "flies",
    rep("sawflies", times = 2),
    rep("wasps", times = 6)
  )
)

# Plotting
cairo_pdf(
  "../figures/Supplementary_Figure_6_tree_final_with_outgroups_included_nodes.pdf", 
  width = tree_pdf_width, 
  height = 22
)
ggtree(final_tree, linewidth = 2, branch.length = "none") + 
  geom_tiplab(
    size = tree_font_size, 
    family = tree_font_family, 
    fontface = "bold"
  ) +
  geom_highlight(
    data = highlighted_nodes, 
    mapping = aes(node = node, fill = group, extend = 1)
  ) +
  geom_rootedge(rootedge = 0.5, linewidth = 2) +
  scale_fill_manual(
    name = "Major group",
    breaks = c("ants", "bees", "wasps", "sawflies", "flies"),
    labels = c("Ants", "Bees", "Wasps", "Sawflies", "Flies"),
    values = c("#785EF0", "#4FD633", "#DC267F", "#F16E1E", "#FFD92F")
  ) +
  theme(legend.text = element_text(size = tree_font_size_legend),
        legend.title = element_text(size = tree_font_size_legend_title),
        legend.margin = tree_legend_margin,
        legend.box.margin = tree_legend_box_margin) +
  scale_y_continuous(limits = c(-1, NA))
dev.off()

#### Main tree
# This is the resulting tree from the MCMCTree run (it's ultrametric)
main_tree <- read.newick("../data/ort80_hymenoptera_67species_100mil_mcmctree.nwk")

# Metadata
eusociality_data <- read.delim(
  "../data/species_metadata_eusociality.tsv",
  sep = "\t"
) %>% select(ABBREV, SOCIALITY_LEVEL, GROUP) %>% 
  dplyr::rename(
    species = ABBREV,
    level = SOCIALITY_LEVEL,
    major_group = GROUP
) %>% filter(!species %in% c("Dmel", "Tcas"))

parasitoidism_data <- read.delim(
  "../data/species_metadata_parasitoidism.tsv",
  sep = "\t"
) %>% select(ABBREV, PARASITOIDISM, GROUP) %>% 
  dplyr::rename(
    species = ABBREV,
    parasitoidism = PARASITOIDISM,
    major_group = GROUP
  ) %>% filter(!species %in% c("Dmel", "Tcas"))

# This is necessary for `gheatmap` to work properly
rownames(eusociality_data) <- eusociality_data$species
eusociality_data <- eusociality_data %>% select(level)
eusociality_data$level <- as.factor(eusociality_data$level)

rownames(parasitoidism_data) <- parasitoidism_data$species
parasitoidism_data <- parasitoidism_data %>% select(parasitoidism)

# Plotting
# The clade names here are the most specific names that encompass all species
# within. We chose commonly identifiable clades and weren't strict to the
# current taxonomic view of the order. Therefore, we include the classic clades
# "Symphyta" and "Parasitica", none of which are monophyletic, but which are
# helpful in picturing the species within (i.e., "Symphyta" are all the sawflies 
# and "Parasitica" the parasitoid wasps within Apocrita)

base_tree_plot <- ggtree(main_tree, linewidth = 2) %<+% full_metadata + 
  geom_tiplab(
    size = tree_font_size, 
    family = tree_font_family, 
    fontface = "bold"
  ) +
  geom_cladelab(
    node = 131, 
    label = "\"Symphyta\"",
    image = "../pics/sawfly.png",
    geom = "image",
    imagecolor = "black", 
    alpha = 1,
    offset = 0.25,
    extend = 2.3,
    barsize = 2,
    offset.text = 0.25
  ) +
  geom_cladelab(
    node = 131,
    label = "\"Symphyta\"",
    color = "black",
    angle = 90,
    offset = 0.25,
    barsize = 0,
    fontsize = tree_font_size,
    family = tree_clade_font_family,
    hjust = 0.5,
    offset.text = 0.05
  ) +
  geom_cladelab(
    node = 72, 
    label = "\"Parasitica\"",
    image = "../pics/parasitoida.png",
    geom = "image",
    imagecolor = "black", 
    alpha = 1,
    offset = 0.25,
    extend = 0.3,
    barsize = 2,
    offset.text = 0.25
  ) +
  geom_cladelab(
    node = 72, 
    label = "\"Parasitica\"",
    color = "black",
    angle = 90,
    offset = 0.25,
    barsize = 0,
    fontsize = tree_font_size,
    family = tree_clade_font_family,
    hjust = 0.5,
    offset.text = 0.05
  ) +
  geom_cladelab(
    node = 81, 
    label = "Vespidae",
    image = "../pics/vespoidea.png",
    geom = "image",
    imagecolor = "black", 
    alpha = 1,
    offset = 0.25,
    extend = 0.3,
    barsize = 2,
    offset.text = 0.25
  ) +
  geom_cladelab(
    node = 81, 
    label = "Vespidae",
    color = "black",
    angle = 90,
    offset = 0.25,
    barsize = 0,
    fontsize = tree_font_size,
    family = tree_clade_font_family,
    hjust = 0.5,
    offset.text = 0.05
  ) +
  geom_cladelab(
    node = 88, 
    label = "Apoidea",
    image = "../pics/apoidea.png",
    geom = "image",
    imagecolor = "black", 
    alpha = 1,
    offset = 0.25,
    extend = 0.3,
    barsize = 2,
    offset.text = 0.25
  ) +
  geom_cladelab(
    node = 88, 
    label = "Apoidea",
    color = "black",
    angle = 90,
    offset = 0.25,
    barsize = 0,
    fontsize = tree_font_size,
    family = tree_clade_font_family,
    hjust = 0.5,
    offset.text = 0.05
  ) +
  geom_cladelab(
    node = 108,
    label = "Formicidae",
    image = "../pics/formicoidea.png",
    geom = "image",
    imagecolor = "black", 
    alpha = 1,
    offset = 0.25,
    extend = 0.3,
    barsize = 2,
    offset.text = 0.25,
    imagesize = 0.035
  ) +
  geom_cladelab(
    node = 108, 
    label = "Formicidae",
    color = "black",
    angle = 90,
    offset = 0.25,
    barsize = 0,
    fontsize = tree_font_size,
    family = tree_clade_font_family,
    hjust = 0.5,
    offset.text = 0.05
  ) +
  geom_rootedge(rootedge = 0.2, linewidth = 2) +
  geom_treescale(fontsize = tree_font_size, family = tree_clade_font_family, y = -1) +
  scale_y_continuous(limits = c(-1, NA))

tree_plot_with_eusociality <- gheatmap(
  base_tree_plot, 
  eusociality_data, 
  custom_column_labels = "Eusociality",
  width = 0.1,
  colnames_angle = -45,
  font.size = tree_font_size,
  offset = 0.7,
  colnames_position = "top",
  colnames_offset_y = 1.2,
  colnames_offset_x = -0.08
) + 
  scale_fill_manual(
    name = "Eusociality level",
    labels = c("Solitary", "Intermediate", "Eussocial"),
    values = c("#CACACA", "#1B9E77", "#4D497D")
  )

tree_plot_with_eusociality <- tree_plot_with_eusociality + new_scale_fill()

cairo_pdf(
  "../../figures/Figure_2_Main_tree.pdf", 
  width = 22, 
  height = 22
)
gheatmap(
    tree_plot_with_eusociality,
    parasitoidism_data,
    custom_column_labels = "Parasitoidism",
    width = 0.1,
    colnames_angle = -45,
    font.size = tree_font_size,
    offset = 1,
    colnames_position = "top",
    colnames_offset_y = 1.6,
    colnames_offset_x = -0.08
) + 
  scale_fill_manual(
    name = "Parasitism/parasitoidism",
    labels = c("Not parasitic", "Parasite/parasitoid"),
    values = c("#716C68", "#F77520")
  ) + theme(
    legend.text = element_text(size = tree_font_size_legend),
    legend.title = element_text(size = tree_font_size_legend_title),
    legend.margin = tree_legend_margin,
    legend.box.margin = tree_legend_box_margin
  )
dev.off()

print("Done! Finished successfuly!")