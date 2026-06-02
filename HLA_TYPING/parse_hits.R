library(jsonlite)
library(tidyverse)

input_dir <- "arcasHLA_results"
out_dir <- "arcasHLA_tables_plots"
dir.create(out_dir, showWarnings = FALSE)

files <- list.files(
  input_dir,
  pattern = "genotype.json$",
  recursive = TRUE,
  full.names = TRUE
)

hla_df <- map_dfr(files, function(f) {
  sample <- basename(dirname(f))
  dat <- fromJSON(f)
  
  map_dfr(names(dat), function(gene) {
    alleles <- dat[[gene]]
    
    tibble(
      Sample = sample,
      Gene = gene,
      Allele1 = ifelse(length(alleles) >= 1, alleles[1], NA),
      Allele2 = ifelse(length(alleles) >= 2, alleles[2], NA),
      Genotype = paste(
        ifelse(length(alleles) >= 1, alleles[1], NA),
        ifelse(length(alleles) >= 2, alleles[2], NA),
        sep = " / "
      )
    )
  })
})

write_csv(hla_df, file.path(out_dir, "arcasHLA_genotype_long.csv"))

hla_wide <- hla_df %>%
  select(Sample, Gene, Genotype) %>%
  pivot_wider(
    names_from = Gene,
    values_from = Genotype,
    names_prefix = "HLA_"
  )

write_csv(hla_wide, file.path(out_dir, "arcasHLA_genotype_wide_by_sample.csv"))

p_sample <- ggplot(hla_df, aes(x = Gene, y = Sample, fill = Gene)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Genotype), size = 3) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Arcas-HLA genotype calls by sample",
    x = "HLA gene",
    y = "Sample ID"
  )

ggsave(
  file.path(out_dir, "arcasHLA_genotype_by_sample.pdf"),
  p_sample,
  width = 12,
  height = 6
)

ggsave(
  file.path(out_dir, "arcasHLA_genotype_by_sample.png"),
  p_sample,
  width = 12,
  height = 6,
  dpi = 300
)

print(hla_wide)
