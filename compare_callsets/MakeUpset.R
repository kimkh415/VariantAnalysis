library(tidyverse)
library(ComplexUpset)
library(RColorBrewer)
library(tidyr)


callers <- c("Cue2", "PBSV", "Sniffles2")
sv_type_levels <- c("INS", "DEL", "DUP", "INV", "INS/DUP")

base_colors <- brewer.pal(n = 4, name = "Dark2") 
sv_type_colors <- c(
  "INS" = base_colors[1], 
  "DEL" = base_colors[2], 
  "INV" = base_colors[3], 
  "DUP" = base_colors[4], 
  "INS/DUP" = "#8c564b" # Distinct brown for mixed
)

sets_ordered <- rev(as.vector(outer(callers, c('INS','DEL','DUP','INV'), paste)))

plot_upset_sv <- function(df,
                          outdir,
                          nam = "SV_caller_agreement",
                          sets_ordered = sets_ordered, # Pass the expanded list, NOT just 'callers'
                          sv_type_levels = sv_type_levels,
                          sv_type_colors = sv_type_colors) {
  upset_df <- df %>%
    uncount(Count) %>%
    mutate(
      # Cue2
      `Cue2 INS`      = FALSE, # Must exist for the plot to work, even if empty
      `Cue2 DUP`      = str_detect(Caller_set, "DUP_C\\b"),
      `Cue2 DEL`      = str_detect(Caller_set, "DEL_C\\b"),
      `Cue2 INV`      = str_detect(Caller_set, "INV_C\\b"),
      # PBSV
      `PBSV INS`      = str_detect(Caller_set, "INS_P\\b"),
      `PBSV DUP`      = str_detect(Caller_set, "DUP_P\\b"),
      `PBSV DEL`      = str_detect(Caller_set, "DEL_P\\b"),
      `PBSV INV`      = str_detect(Caller_set, "INV_P\\b"),
      # Sniffles2
      `Sniffles2 INS` = str_detect(Caller_set, "INS_S\\b"),
      `Sniffles2 DUP` = str_detect(Caller_set, "DUP_S\\b"),
      `Sniffles2 DEL` = str_detect(Caller_set, "DEL_S\\b"),
      `Sniffles2 INV` = str_detect(Caller_set, "INV_S\\b")
    ) %>%
    mutate(SVTYPE = factor(SVTYPE, levels = sv_type_levels))

  # Plot
  pdf(file.path(outdir, paste0(nam, ".pdf")), width = 15, height = 7)
  p <- upset(
    upset_df,
    intersect = sets_ordered,
    min_size = 0,
    #sort_sets = FALSE,
    sort_intersections_by = c('degree', 'cardinality'),
    width_ratio = 0.2,
    base_annotations = list(
      'Intersection size' = (
        intersection_size(
          mapping = aes(fill = SVTYPE),
          counts = TRUE,
          text = list(vjust = -0.5, hjust = 0.1, angle = 45)
        ) + scale_fill_manual(values = sv_type_colors, name = "SV Type")
      )
    ),
    # Assign colors to rows
    queries = list(
      upset_query(set="Cue2 INS",      fill=sv_type_colors['INS']),
      upset_query(set="PBSV INS",      fill=sv_type_colors['INS']),
      upset_query(set="Sniffles2 INS", fill=sv_type_colors['INS']),
      upset_query(set="Cue2 DEL",      fill=sv_type_colors['DEL']),
      upset_query(set="PBSV DEL",      fill=sv_type_colors['DEL']),
      upset_query(set="Sniffles2 DEL", fill=sv_type_colors['DEL']),
      upset_query(set="Cue2 DUP",      fill=sv_type_colors['DUP']),
      upset_query(set="PBSV DUP",      fill=sv_type_colors['DUP']),
      upset_query(set="Sniffles2 DUP", fill=sv_type_colors['DUP']),
      upset_query(set="Cue2 INV",      fill=sv_type_colors['INV']),
      upset_query(set="PBSV INV",      fill=sv_type_colors['INV']),
      upset_query(set="Sniffles2 INV", fill=sv_type_colors['INV'])
    )
  ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  print(p)
  dev.off()
}


input_file="upset_data_final.tsv"
outdir="New_Fig2a/"
df <- read_tsv(input_file, show_col_types = FALSE)
plot_upset_sv(df, outdir, "Fig2a_final", sets_ordered, names(sv_type_colors), sv_type_colors)

