#Generating violin plot from scRNAseq analysis:
library(tidyverse)

infile <- "Hnf4-data-new-muscle&oe-switch order.csv"  # <- adjust path if needed
df <- read_csv(infile, show_col_types = FALSE)

gene_id <- "FBgn0004914"

# ---- Make a long table with columns: gene, value, Sample, Cell type ----
# Case 1: already long (has a gene/feature column + value)
if (any(names(df) %in% c("gene", "Gene", "feature", "Feature"))) {
  gene_col <- intersect(names(df), c("gene", "Gene", "feature", "Feature"))[1]
  df_long <- df %>%
    rename(gene = all_of(gene_col)) %>%
    filter(gene == gene_id)
  
  # Case 2: wide (gene is a column name)
} else if (gene_id %in% names(df)) {
  # If "value" already exists but is not the gene column, overwrite value from the gene column
  df_long <- df %>%
    mutate(value = .data[[gene_id]]) %>%
    mutate(gene = gene_id)
  
  # Fallback: pivot all non-metadata columns into gene/value
} else {
  meta_cols <- intersect(names(df), c("Cell", "Sample", "Cluster", "Cell type"))
  df_long <- df %>%
    pivot_longer(cols = -all_of(meta_cols), names_to = "gene", values_to = "value") %>%
    filter(gene == gene_id)
}

# Ensure ordering (optional)
df_long <- df_long %>%
  mutate(
    Sample = factor(Sample),
    `Cell type` = factor(`Cell type`, levels = c("Fat body", "Muscle", "Oenocyte"))
  )

# ---- Plot (facet by tissue/cell type; two violins per facet for each Sample) ----
p <- ggplot(df_long, aes(x = Sample, y = value, fill = Sample)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.9, color = "grey30") +
  geom_jitter(width = 0.15, height = 0, size = 0.5, alpha = 0.05) +
  facet_wrap(~`Cell type`, nrow = 1) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(
    x = NULL,
    y = paste0(gene_id, " expression (value)"),
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

print(p)
