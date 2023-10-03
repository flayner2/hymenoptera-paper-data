#################################### SETUP #####################################
# Loading libraries
library(ggplot2)
library(grid)
library(dplyr)
library(ggpubr)
library(forcats)
library(RColorBrewer)
library(Cairo)

# Clearing the workspace
rm(list = ls())

# Set the current working directory to this script's directory. You will need
# to change this to the path were it's located on your computer
setwd("~/Documents/LAB/eusociality/work/masters_paper/doc/paper/supplementary/code/")

# Setting a seed for reproductibility
seed <- 666
set.seed(seed)

# Matrix containing the BUSCO completeness metrics for all species
busco_data <- read.table("../data/busco_metadata.tsv", header = TRUE, sep = "\t")

# Parameters for the figures
output <- "../figures/Supplementary_Figure_2_BUSCO_all_metrics.pdf"
output_single <- "../figures/Supplementary_Figure_3_BUSCO_single_copy.pdf"
width <- 24
height <- 15
width_single <- 10
height_single <- 18
monospaced_family <- "mono"
default_family <- "sans serif"
colours <- c("#785EF0", "#4FD633", "#F16E1E", "#DC267F")
labels <- c("Ants", "Bees", "Sawflies", "Wasps")
y_breaks <- c(0, 25, 50, 75, 100)
y_limits <- c(0, 100)
y_labels <- c("0", "25", "50", "75", "100")
font_size_legend <- 14
font_size_legend_title <- 16
font_size_y <- 14
title_bottom_margin <- 10
x_title_size <- 14
line_color <- "black"
line_size <- 1
line_type <- 2

print("Generating plots...")

############################# INDIVIDUAL PLOTS #################################
# Total completeness (single + duplicated)
complete <- busco_data %>% 
  arrange(C) %>% 
  mutate(short = factor(short, levels = short)) %>% 
  ggplot() + 
  geom_bar(aes(y = C, x = short, fill = group), stat = "identity") + 
  xlab(NULL) +
  ylab("Complete (%)") +
  scale_y_continuous(
    position = "right", 
    breaks = y_breaks, 
    limits = y_limits, 
    labels = y_labels
  ) +
  coord_flip() +
  labs(fill = NULL) +
  scale_fill_manual(
    name = "Major group",
    labels = labels,
    values = colours
  ) +
  theme_light() +
  theme(
    axis.text.y = element_text(
      family = monospaced_family, 
      colour = ifelse(arrange(busco_data, C)$S < 90, "red", "black"),
      face = "bold",
      size = font_size_y
    ),
    axis.title.x = element_text(size = x_title_size),
    axis.title.x.top = element_text(margin = margin(b = title_bottom_margin)),
    legend.text = element_text(size = font_size_legend),
    legend.title = element_text(size = font_size_legend_title),
    text = element_text(family = default_family)
  ) +
  geom_hline(
    yintercept = 90, 
    color = line_color, 
    size = line_size, 
    linetype = line_type
  )

# Complete single copy
single <- busco_data %>% 
  arrange(S) %>% 
  mutate(short = factor(short, levels = short)) %>% 
  ggplot() + 
  geom_bar(aes(y = S, x = short, fill = group), stat = "identity") + 
  xlab(NULL) +
  ylab("Single copy (%)") +
  scale_y_continuous(
    position = "right", 
    breaks = y_breaks, limits = y_limits, labels = y_labels
  ) +
  coord_flip() +
  scale_fill_manual(
    labels = labels,
    values = colours
  ) +
  theme_light() +
  theme(
    axis.text.y = element_text(
      family = monospaced_family, 
      colour = ifelse(arrange(busco_data, S)$S < 90, "red", "black"),
      face = "bold",
      size = font_size_y
    ),
    axis.title.x = element_text(size = x_title_size),
    axis.title.x.top = element_text(margin = margin(b = title_bottom_margin)),
    legend.title = element_text(size = font_size_legend_title),
    text = element_text(family = default_family)
  ) +
  geom_hline(
    yintercept = 90, 
    color = line_color, 
    size = line_size, 
    linetype = line_type
  )

# Complete duplicated
dup <- busco_data %>% 
  arrange(D) %>% 
  mutate(short = factor(short, levels = short)) %>% 
  ggplot() + 
  geom_bar(aes(y = D, x = short, fill = group), stat = "identity") + 
  xlab(NULL) +
  ylab("Duplicated (%)") +
  scale_y_continuous(
    position = "right", 
    breaks = y_breaks, 
    limits = y_limits, 
    labels = y_labels
  ) +
  coord_flip() +
  scale_fill_manual(
    labels = labels,
    values = colours
  ) +
  theme_light() +
  theme(
    axis.text.y = element_text(
      family = monospaced_family, 
      colour = ifelse(arrange(busco_data, D)$S < 90, "red", "black"),
      face = "bold",
      size = font_size_y
    ),
    axis.title.x = element_text(size = x_title_size),
    axis.title.x.top = element_text(margin = margin(b = title_bottom_margin)),
    legend.title = element_text(size = font_size_legend_title),
    text = element_text(family = default_family)
  ) +
  geom_hline(
    yintercept = 10, 
    color = line_color, 
    size = line_size, 
    linetype = line_type
  )

# Fragmented
frag <- busco_data %>% 
  arrange(F) %>% 
  mutate(short = factor(short, level = short)) %>% 
  ggplot() + 
  geom_bar(
    aes(y = F, x = short, fill = group), stat = "identity") + 
  xlab(NULL) +
  ylab("Fragmented (%)") +
  scale_y_continuous(
    position = "right", 
    breaks = y_breaks, 
    limits = y_limits, labels = y_labels
  ) +
  coord_flip() +
  scale_fill_manual(
    labels = labels,
    values = colours
  ) +
  theme_light() +
  theme(
    axis.text.y = element_text(
      family = monospaced_family, 
      colour = ifelse(arrange(busco_data, F)$S < 90, "red", "black"),
      face = "bold",
      size = font_size_y
    ),
    axis.title.x = element_text(size = x_title_size),
    axis.title.x.top = element_text(margin = margin(b = title_bottom_margin)),
    legend.title = element_text(size = font_size_legend_title),
    text = element_text(family = default_family)
  ) +
  geom_hline(
    yintercept = 10, 
    color = line_color, 
    size = line_size, 
    linetype = line_type
  )

# Missing
miss <- busco_data %>% 
  arrange(M) %>% 
  mutate(short = factor(short, level = short)) %>% 
  ggplot() + 
  geom_bar(aes(y = M, x = short, fill = group), stat = "identity") + 
  xlab(NULL) +
  ylab("Missing (%)") +
  scale_y_continuous(
    position = "right", 
    breaks = y_breaks, 
    limits = y_limits, 
    labels = y_labels 
  ) +
  coord_flip() +
  scale_fill_manual(
    labels = labels,
    values = colours
  ) +
  theme_light() +
  theme(
    axis.text.y = element_text(
      family = monospaced_family, 
      colour = ifelse(arrange(busco_data, M)$S < 90, "red", "black"),
      face = "bold",
      size = font_size_y
    ),
    axis.title.x = element_text(size = x_title_size),
    axis.title.x.top = element_text(margin = margin(b = title_bottom_margin)),
    legend.title = element_text(size = font_size_legend_title),
    text = element_text(family = default_family)
  ) +
  geom_hline(
    yintercept = 10, 
    color = line_color, 
    size = line_size, 
    linetype = line_type
  )

################################ FINAL PLOTS ###################################
# Supplementary 1, with all metrics
figure <- ggarrange(
  complete, single, dup, frag, miss, 
  nrow = 1, ncol = 5, legend = "right", common.legend = TRUE
)
cairo_pdf(output, width, height)
annotate_figure(
  figure, 
  left = text_grob("Species", rot = 90, family = default_family)
)
dev.off()

# Supplementary 2, only single-copy
cairo_pdf(output_single, width_single, height_single)
print(
  single + 
    ylab("Single copy completeness (%)") +
    guides(fill = guide_legend("Major group")) +
    theme(
      axis.text.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 20),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 20),
      text = element_text(family = default_family)
  )
)
dev.off()
