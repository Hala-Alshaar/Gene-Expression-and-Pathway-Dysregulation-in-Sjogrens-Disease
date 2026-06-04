## Enrichment analysis 
################################################################
# 31. INTEGRATED FUNCTIONAL ENRICHMENT PIPELINE
# (g:Profiler + pathfindR)
################################################################

library(gprofiler2)
library(pathfindR)
library(dplyr)
library(ggplot2)

# ==========================================================
# 1. LIMMA RESULTS -> GENE SYMBOLS
# ==========================================================

deg_table <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH"
)

annot <- fData(gset)

gene_map <- annot$Gene.symbol
names(gene_map) <- annot$ID

deg_table$Gene.symbol <- gene_map[rownames(deg_table)]

deg_table <- deg_table[
  !is.na(deg_table$Gene.symbol) &
    deg_table$Gene.symbol != "",
]

# ==========================================================
# 2. SIGNIFICANT DEGs
# ==========================================================

sig_genes <- deg_table[
  deg_table$adj.P.Val < 0.05 &
    abs(deg_table$logFC) > 1,
]

gene_symbols <- unique(
  sig_genes$Gene.symbol
)

cat("Significant genes:", length(gene_symbols), "\n")

write.csv(
  sig_genes,
  "Significant_DEGs_GeneSymbols.csv",
  row.names = FALSE
)

# ==========================================================
# 3. G:PROFILER ENRICHMENT
# ==========================================================

gost_res <- gost(
  query = gene_symbols,
  organism = "hsapiens",
  correction_method = "fdr",
  significant = TRUE
)

gprof_res <- gost_res$result

write.csv(
  gprof_res,
  "gProfiler_Results.csv",
  row.names = FALSE
)

# ==========================================================
# FIGURE 10 STYLE
# ==========================================================
library(ggplot2)
library(viridis)

p_gprof <- ggplot(
  gprof_plot_data,
  aes(
    x = source,
    y = -log10(p_value),
    color = source
  )
) +
  geom_jitter(
    aes(size = intersection_size),
    width = 0.25,
    alpha = 0.8
  ) +
  scale_color_viridis_d(option = "turbo") +
  scale_size_continuous(range = c(3, 10)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Functional Enrichment Summary",
    x = "Database",
    y = expression(-log[10](p-value)),
    color = "Database",
    size = "Gene Count"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(p_gprof)

# ==========================================================
# 4. PATHFINDR INPUT
# ==========================================================

pathfindR_input <- sig_genes[, c(
  "Gene.symbol",
  "logFC",
  "adj.P.Val"
)]

colnames(pathfindR_input) <- c(
  "Gene.symbol",
  "logFC",
  "adj.P.Val"
)

pathfindR_input <- pathfindR_input[
  !duplicated(pathfindR_input$Gene.symbol),
]

write.csv(
  pathfindR_input,
  "pathfindR_Input.csv",
  row.names = FALSE
)

# ==========================================================
# 5. REACTOME ANALYSIS
# ==========================================================

reactome_res <- run_pathfindR(
  input = pathfindR_input,
  gene_sets = "Reactome",
  pin_name_path = "Biogrid",
  p_val_threshold = 0.05
)

write.csv(
  reactome_res,
  "Reactome_PathfindR.csv",
  row.names = FALSE
)

top_reactome <- reactome_res[
  order(reactome_res$lowest_p),
]

top_reactome <- head(top_reactome, 20)

p_reactome <- ggplot(
  top_reactome,
  aes(
    x = Fold_Enrichment,
    y = reorder(
      Term_Description,
      Fold_Enrichment
    )
  )
) +
  geom_point(
    aes(
      color = -log10(lowest_p),
      size = Up_regulated
    )
  ) +
  theme_bw() +
  labs(
    title = "Reactome Pathway Enrichment",
    x = "Fold Enrichment",
    y = ""
  )

print(p_reactome)

ggsave(
  "Figure11_Reactome_Pathways.png",
  p_reactome,
  width = 9,
  height = 6
)

# ==========================================================
# 6. GO BIOLOGICAL PROCESS
# ==========================================================

gobp_res <- run_pathfindR(
  input = pathfindR_input,
  gene_sets = "GO-BP",
  pin_name_path = "Biogrid",
  p_val_threshold = 0.05
)

write.csv(
  gobp_res,
  "GO_BP_PathfindR.csv",
  row.names = FALSE
)

top_bp <- gobp_res[
  order(gobp_res$lowest_p),
]

top_bp <- head(top_bp, 20)

p_bp <- ggplot(
  top_bp,
  aes(
    x = Fold_Enrichment,
    y = reorder(
      Term_Description,
      Fold_Enrichment
    )
  )
) +
  geom_point(
    aes(
      color = -log10(lowest_p),
      size = Up_regulated
    )
  ) +
  theme_bw() +
  labs(
    title = "GO Biological Process Enrichment"
  )

print(p_bp)

ggsave(
  "Figure13_GO_BP.png",
  p_bp,
  width = 9,
  height = 6
)

# ==========================================================
# 7. GO CELLULAR COMPONENT
# ==========================================================

gocc_res <- run_pathfindR(
  input = pathfindR_input,
  gene_sets = "GO-CC",
  pin_name_path = "Biogrid",
  p_val_threshold = 0.05
)

write.csv(
  gocc_res,
  "GO_CC_PathfindR.csv",
  row.names = FALSE
)

top_cc <- gocc_res[
  order(gocc_res$lowest_p),
]

top_cc <- head(top_cc, 20)

p_cc <- ggplot(
  top_cc,
  aes(
    x = Fold_Enrichment,
    y = reorder(
      Term_Description,
      Fold_Enrichment
    )
  )
) +
  geom_point(
    aes(
      color = -log10(lowest_p),
      size = Up_regulated
    )
  ) +
  theme_bw() +
  labs(
    title = "GO Cellular Component Enrichment"
  )

print(p_cc)

ggsave(
  "Figure14_GO_CC.png",
  p_cc,
  width = 9,
  height = 6
)

# ==========================================================
# 8. TERM-GENE NETWORK (FIGURE 12 STYLE)
# ==========================================================

try(
  term_gene_graph(
    enrichment_results = reactome_res,
    top_terms = 10
  )
)

# ==========================================================
# 9. SAVE EVERYTHING
# ==========================================================

save(
  sig_genes,
  gprof_res,
  reactome_res,
  gobp_res,
  gocc_res,
  file = "Enrichment_Analysis_Objects.RData"
)

cat("\n")
cat("====================================\n")
cat("ENRICHMENT PIPELINE COMPLETED\n")
cat("====================================\n")
cat("Input type: GENE SYMBOLS\n")
cat("g:Profiler completed\n")
cat("Reactome completed\n")
cat("GO-BP completed\n")
cat("GO-CC completed\n")
cat("Figures exported\n")
cat("====================================\n")