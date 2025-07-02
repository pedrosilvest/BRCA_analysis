


# Import ----
library(vcfR)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(ggtext)
library(RColorBrewer)
library(patchwork)
library(forcats) 
library(gtools)

# Data Opening ----

### snv
strelka_c2040=read.vcfR(
  "/Users/pedrosilvestre/Desktop/NFS_Pedro/BRCA_project/strelka/C2040/StrelkaSomaticWorkflow/results/variants/somatic.snvs.vcf.gz",
  verbose = F)

strelka_c2041=read.vcfR(
  "/Users/pedrosilvestre/Desktop/NFS_Pedro/BRCA_project/strelka/C2041/StrelkaSomaticWorkflow/results/variants/somatic.snvs.vcf.gz",
  verbose = F)

strelka_c2040_ind = read.vcfR(
  "/Users/pedrosilvestre/Desktop/NFS_Pedro/BRCA_project/strelka/C2040/StrelkaSomaticWorkflow/results/variants/somatic.indels.vcf.gz",
  verbose = F
)

### indel
strelka_c2041_ind = read.vcfR(
  "/Users/pedrosilvestre/Desktop/NFS_Pedro/BRCA_project/strelka/C2041/StrelkaSomaticWorkflow/results/variants/somatic.indels.vcf.gz",
  verbose = F
)



varscan_c2040 = read.table(
  file="/Users/pedrosilvestre/Desktop/NFS_Pedro/BRCA_project/samples/VarScan/C2040_mpileup_3004936/C2040.snp",
                   sep = "\t", header = TRUE)

varscan_c2041 = read.table(
  file="/Users/pedrosilvestre/Desktop/NFS_Pedro/BRCA_project/samples/VarScan/C2041_mpileup_3004935/C2041.snp",
                   sep = "\t", header = TRUE)

# Strelka Analysis SNV ----
# 1. Process data----

process_strelka <- function(df) {
  df_NAlt<-data.frame(N=1:nrow(df@fix),Chr=df@fix[,1], Pos=as.numeric(df@fix[,2]),
                      NAlt=lengths(strsplit(df@fix[,5],",")) )
  
  AF<-extract.gt(df, element='FREQ', as.numeric=TRUE)
  DP<-extract.gt(df, element='DP', as.numeric=TRUE)
  
  VAF<-extract.gt(df, element='VAF', as.numeric=TRUE)
  
  
  format_field <- as.list(strsplit(df@fix[, "INFO"], ";"))
  
  sgt_values <- sapply(format_field, function(x) {
    # Find the element that starts with "SGT="
    sgt_entry <- grep("^SGT=", x, value = TRUE)
    # Remove the "SGT=" prefix to keep only the value
    if (length(sgt_entry) > 0) {
      sub("^SGT=", "", sgt_entry)
    } else {
      NA
    }
  })
  
  
  df_NAlt<-data.frame(N=1:nrow(df@fix),Chr=df@fix[,1], Pos=as.numeric(df@fix[,2]), Ref=df@fix[,4], Alt=df@fix[,5],
                      SGT = sgt_values, NAlt=lengths(strsplit(df@fix[,5],",")), Filter=df@fix[,"FILTER"])
  
  # Display the extracted SGT values
  unique(sgt_values)
  
  
  filtered_sgt_values <- sgt_values[!grepl("^(.+)->\\1$", sgt_values)]
  
  
  # Display the filtered unique SGT values
  unique(filtered_sgt_values)
  
  
  print(colnames(df@gt[,1])[-1])
  
  
  DP<-extract.gt(df, element='DP', as.numeric=TRUE) #this collects the DP for each
  FDP<-extract.gt(df, element='FDP', as.numeric=TRUE) #Reads filtered out
  AU<-extract.gt(df, element='AU', as.numeric=TRUE)
  CU<-extract.gt(df, element='CU', as.numeric=TRUE)
  GU<-extract.gt(df, element='GU', as.numeric=TRUE)
  TU<-extract.gt(df, element='TU', as.numeric=TRUE)
  
  
  get_ref_depth <- function(ref, AU, CU, GU, TU) {
    switch(ref,
           "A" = AU,
           "C" = CU,
           "G" = GU,
           "T" = TU,
           NA) # Default to NA if no match
  }
  
  # Function to get the alternate depth for normal or tumor samples
  get_alt_depth <- function(alt, AU, CU, GU, TU) {
    switch(alt,
           "A" = AU,
           "C" = CU,
           "G" = GU,
           "T" = TU,
           NA) # Default to NA if no match
  }
  
  # Create the new columns
  df_NAlt$DP_Normal <- DP[, "NORMAL"] # Assuming 'Normal' is a column in DP
  df_NAlt$DP_Tumor <- DP[, "TUMOR"] # Assuming 'Tumour' is a column in DP
  
  df_NAlt$FDP_Normal <- FDP[, "NORMAL"]
  df_NAlt$FDP_Tumor <- FDP[, "TUMOR"]
  
  
  df_NAlt$DP_Ref_Normal <- mapply(get_ref_depth, df_NAlt$Ref, AU[, "NORMAL"], CU[, "NORMAL"], GU[, "NORMAL"], TU[, "NORMAL"])
  df_NAlt$DP_Ref_Tumor <- mapply(get_ref_depth, df_NAlt$Ref, AU[, "TUMOR"], CU[, "TUMOR"], GU[, "TUMOR"], TU[, "TUMOR"])
  
  df_NAlt$DP_Alt_Normal <- mapply(get_alt_depth, df_NAlt$Alt, AU[, "NORMAL"], CU[, "NORMAL"], GU[, "NORMAL"], TU[, "NORMAL"])
  df_NAlt$DP_Alt_Tumor <- mapply(get_alt_depth, df_NAlt$Alt, AU[, "TUMOR"], CU[, "TUMOR"], GU[, "TUMOR"], TU[, "TUMOR"])
  
  df_NAlt <- df_NAlt[grepl("chr", df_NAlt$Chr), ]
  
  df_NAlt <- df_NAlt %>%
    mutate(VAF_Normal = DP_Alt_Normal / (DP_Ref_Normal + DP_Alt_Normal),
           VAF_Tumor = DP_Alt_Tumor / (DP_Ref_Tumor + DP_Alt_Tumor),
           diff_VAF = abs(VAF_Normal - VAF_Tumor))
  
  return(df_NAlt)
  
}

#head(strelka_c2041@fix)

# Process sample
st_c2041_df <- process_strelka(strelka_c2041)
colnames(st_c2041_df)
st_c2040_df <- process_strelka(strelka_c2040)
colnames(st_c2040_df)


st_c2040_filtered <- st_c2040_df %>% filter(Filter == "PASS")
st_c2041_filtered <- st_c2041_df %>% filter(Filter == "PASS")

# 2. General Differences ----



# Mutation Classification using Ref and Alt
classify_mutation <- function(ref, alt) {
  case_when(
    ref %in% c("A", "G") & alt %in% c("A", "G") ~ "Transition",  # A <-> G
    ref %in% c("C", "T") & alt %in% c("C", "T") ~ "Transition",  # C <-> T
    ref %in% c("A", "C", "G", "T") & alt %in% c("A", "C", "G", "T") ~ "Transversion",  # Any other change
    TRUE ~ NA_character_  # Handle missing or invalid cases
  )
}

# Apply to C2040
st_c2040_filtered <- st_c2040_filtered %>%
  mutate(MutationType = classify_mutation(Ref, Alt))

# Apply to C2041
st_c2041_filtered <- st_c2041_filtered %>%
  mutate(MutationType = classify_mutation(Ref, Alt))

# Count Transitions and Transversions
table(st_c2040_filtered$MutationType)
table(st_c2041_filtered$MutationType)



# Compare total mutation counts
cat("C2040 total mutations:", nrow(st_c2040_filtered), "\n")
cat("C2041 total mutations:", nrow(st_c2041_filtered), "\n")

# Mutations per chromosome
mutations_per_chr <- function(df) {
  df %>% group_by(Chr) %>% summarise(Count = n()) %>% arrange(desc(Count))
}

c2040_chr <- mutations_per_chr(st_c2040_filtered)
c2041_chr <- mutations_per_chr(st_c2041_filtered)


# Combine and label datasets
mutation_data <- bind_rows(
  c2040_chr %>% mutate(Sample = "C2040"),
  c2041_chr %>% mutate(Sample = "C2041")
)

mutation_data <- mutation_data %>%
  mutate(Chr = factor(Chr, levels = mixedsort(unique(Chr))))  # Sort naturally (chr1, chr2, ..., chrX, chrY)

# Plot with sorted chromosomes
ggplot(mutation_data %>% filter(Chr != "chrY"), aes(x = Chr, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +  
  theme_classic(base_size = 16) +  
  labs(
    title = "Mutation Distribution Across Chromosomes",
    x = "Chromosome",
    y = "Mutation Count",
    fill = "Sample"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1),  
    legend.position = "top",  
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  guides(fill = guide_legend(reverse = TRUE))  


# Create unique identifiers for mutations based on Chr and Pos
st_c2040_filtered <- st_c2040_filtered %>% mutate(MutationID = paste(Chr, Pos, Ref, Alt, sep="_"))
st_c2041_filtered <- st_c2041_filtered %>% mutate(MutationID = paste(Chr, Pos, Ref, Alt, sep="_"))

# Identify shared and unique mutations
shared_mutations <- intersect(st_c2040_filtered$MutationID, st_c2041_filtered$MutationID)
unique_c2040 <- setdiff(st_c2040_filtered$MutationID, st_c2041_filtered$MutationID)
unique_c2041 <- setdiff(st_c2041_filtered$MutationID, st_c2040_filtered$MutationID)

cat("Shared Mutations:", length(shared_mutations), "\n")
cat("Unique to C2040:", length(unique_c2040), "\n")
cat("Unique to C2041:", length(unique_c2041), "\n")

# Plot VAF distributions
ggplot(bind_rows(st_c2040_filtered %>% mutate(Sample="WT"), 
                 st_c2041_filtered %>% mutate(Sample="BRCA -/+")), 
       aes(x=VAF_Tumor, fill=Sample)) +
  geom_density(alpha=0.5) +
  theme_minimal() +
  labs(title="VAF Distribution in Tumor Samples", x="VAF", y="Density")



# Function to count mutation signatures (Ref -> Alt)
count_mutation_signatures <- function(df) {
  df %>%
    filter(!is.na(Ref) & !is.na(Alt)) %>%  # Ensure valid values
    group_by(Ref, Alt) %>%
    summarise(Count = n(), .groups = "drop") %>%
    arrange(desc(Count))  # Sort by most common mutations
}

# Apply to C2040
mutation_counts_c2040 <- count_mutation_signatures(st_c2040_filtered)

# Apply to C2041
mutation_counts_c2041 <- count_mutation_signatures(st_c2041_filtered)

# View results
print(mutation_counts_c2040)
print(mutation_counts_c2041)

mutation_order <- c("C→A", "C→G", "C→T", "T→A", "T→C", "T→G")  # Standard order

# Function to count and classify mutation signatures
count_mutation_signatures <- function(df) {
  df %>%
    filter(!is.na(Ref) & !is.na(Alt)) %>%  # Ensure valid Ref/Alt bases
    mutate(mutation_type = case_when(
      Ref == "A" & Alt == "C" ~ "T→G",
      Ref == "A" & Alt == "G" ~ "T→C",
      Ref == "A" & Alt == "T" ~ "T→A",
      Ref == "G" & Alt == "A" ~ "C→T",
      Ref == "G" & Alt == "C" ~ "C→G",
      Ref == "G" & Alt == "T" ~ "C→A",
      TRUE ~ paste0(Ref, "→", Alt)  # Keep existing C/T mutations
    )) %>%
    filter(mutation_type %in% mutation_order) %>%
    group_by(mutation_type) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(mutation_type = factor(mutation_type, levels = mutation_order))  # Set fixed order
}


plot_mutation_signatures <- function(df, sample_name) {
  mutation_counts <- count_mutation_signatures(df)  # Get processed mutation counts
  
  ggplot(mutation_counts, aes(x = mutation_type, y = Count, fill = mutation_type)) +
    geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +  # Bars for mutation counts
    geom_text(aes(label = Count), hjust = -0.2, size = 2, fontface = "bold") +  # Add count labels
    scale_fill_manual(values = c("#D55E00", "#0072B2", "#009E73", "#E69F00", "#56B4E9", "#CC79A7")) +  # Nature-like colors
    coord_flip() +  # Flip for horizontal readability
    labs(
      title = paste("Mutation in", sample_name),
      x = "Mutation Type",
      y = "Mutation Count"
    ) +
    theme_classic(base_size = 16) +  # Clean, professional theme
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(face = "bold", size = 16),
      axis.title.y = element_text(face = "bold", size = 16)
    )
}

# Apply function to each sample dataset
plot_c2040 <- plot_mutation_signatures(st_c2040_filtered, "WT")
plot_c2041 <- plot_mutation_signatures(st_c2041_filtered, "BRCA -/+")

# Display the plots
print(plot_c2040)
print(plot_c2041)

plot_c2040 + plot_c2041


# 3. Check left and right nucleotides ----

# Load the required libraries
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)



get_flanking_nucleotides <- function(chr, pos) {
  left_nuc <- as.character(BSgenome.Hsapiens.UCSC.hg38[[chr]][pos - 1])
  right_nuc <- as.character(BSgenome.Hsapiens.UCSC.hg38[[chr]][pos + 1])
  return(c(left_nuc, right_nuc))
}

# Add the columns to your dataset
st_c2040_filtered <- st_c2040_filtered %>%
  rowwise() %>%
  mutate(
    nucs = list(get_flanking_nucleotides(Chr, Pos)),  # Store result as a list
    left_nucleotide = nucs[[1]],                      # Extract first element
    right_nucleotide = nucs[[2]]                      # Extract second element
  ) %>%
  select(-nucs) %>%                                   # Remove temporary column
  ungroup()



st_c2041_filtered <- st_c2041_filtered %>%
  rowwise() %>%
  mutate(
    nucs = list(get_flanking_nucleotides(Chr, Pos)),  # Store result as a list
    left_nucleotide = nucs[[1]],                      # Extract first element
    right_nucleotide = nucs[[2]]                      # Extract second element
  ) %>%
  select(-nucs) %>%                                   # Remove temporary column
  ungroup()


st_c2040_filtered <- st_c2040_filtered %>%
  mutate(trinucleotide_context = paste0(left_nucleotide, "-", right_nucleotide))

st_c2041_filtered <- st_c2041_filtered %>%
  mutate(trinucleotide_context = paste0(left_nucleotide, "-", right_nucleotide))

unique(st_c2041_filtered$trinucleotide_context)

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)  # For annotations

library(ggplot2)
library(dplyr)

plot_mutation_signatures <- function(data, dataset_name) {
  # Ensure required columns
  if (!all(c("Ref", "Alt", "trinucleotide_context") %in% names(data))) {
    stop("Data must contain columns: Ref, Alt, trinucleotide_context")
  }
  
  # Convert mutations to the six canonical types
  data <- data %>%
    mutate(mutation_type = case_when(
      Ref == "A" & Alt == "C" ~ "T>C",
      Ref == "A" & Alt == "G" ~ "T>G",
      Ref == "A" & Alt == "T" ~ "T>A",
      Ref == "G" & Alt == "A" ~ "C>T",
      Ref == "G" & Alt == "C" ~ "C>G",
      Ref == "G" & Alt == "T" ~ "C>A",
      TRUE ~ paste0(Ref, ">", Alt)  # Keep C/T mutations unchanged
    )) %>%
    filter(mutation_type %in% c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))  # Keep only valid types
  
  # Count mutations
  mutation_counts <- data %>%
    count(mutation_type, trinucleotide_context, name = "count")
  
  # Set custom COSMIC colors
  cosmic_colors <- c("C>A" = "#0099CC",  # Blue
                     "C>G" = "#000000",  # Black
                     "C>T" = "#CC0000",  # Red
                     "T>A" = "#999999",  # Gray
                     "T>C" = "#66CC66",  # Green
                     "T>G" = "#FF9999")  # Pink
  
  # Get total mutations for annotation
  total_mutations <- sum(mutation_counts$count)
  
  # Plot
  p <- ggplot(mutation_counts, aes(x = trinucleotide_context, y = count, fill = mutation_type)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = cosmic_colors) +
    facet_grid(. ~ mutation_type, scales = "free_x", space = "free_x") +  # Align bars for each mutation type
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_blank(),  # Remove x-axis labels
      axis.ticks.x = element_blank(),  # Remove x-axis ticks
      axis.text.y = element_text(size = 12),
      strip.text = element_text(face = "bold", size = 14, color = "white", margin = margin(t = 5, b = 5)),  # Mutation type as x-axis title
      strip.background = element_rect(fill = "black"),  # Background for mutation labels
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.title.x = element_blank()  # Removes x-axis label
    ) +
    labs(
      title = paste0(dataset_name, " (", total_mutations, " SNVs)"),
      y = "Mutation Count"  # Keep only the y-axis label
    )
  
  # Print plot
  print(p)
}

# Generate separate plots for each dataset
plot_mutation_signatures(st_c2040_filtered, "WT")
plot_mutation_signatures(st_c2041_filtered, "BRCA -/+")
plot_mutation_signatures(st_c2040_filtered, "WT") + plot_mutation_signatures(st_c2041_filtered, "BRCA -/+") + 
  plot_layout(ncol = 1, heights = c(2, 2))



# Combine both datasets with a "sample" identifier
combined_data <- bind_rows(
  st_c2040_filtered %>% mutate(sample = "BRCA +/+ (WT)"),
  st_c2041_filtered %>% mutate(sample = "BRCA -/+ (Het)")
)

# Create Ref → Alt mutation type column
combined_data <- combined_data %>%
  mutate(mutation_type = paste0(Ref, "→", Alt))

# should be this since DNA is double strand
combined_data <- combined_data %>%
  mutate(mutation_type = case_when(
    Ref == "A" & Alt == "C" ~ "T→G",
    Ref == "A" & Alt == "G" ~ "T→C",
    Ref == "A" & Alt == "T" ~ "T→A",
    Ref == "G" & Alt == "A" ~ "C→T",
    Ref == "G" & Alt == "C" ~ "C→G",
    Ref == "G" & Alt == "T" ~ "C→A",
    TRUE ~ paste0(Ref, "→", Alt)  # Keep C/T mutations unchanged
  )) %>%
  filter(mutation_type %in% c("C→A", "C→G", "C→T", "T→A", "T→C", "T→G"))  # Keep only valid types


# Count occurrences for each sample
mutation_counts <- combined_data %>%
  count(sample, mutation_type, trinucleotide_context, name = "count") %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = 0) %>%
  mutate(ratio = (`BRCA +/+ (WT)` + 1) / (`BRCA -/+ (Het)` + 1))  # Avoid division by zero

mutation_counts <- mutation_counts %>%
  mutate(ratio = ifelse(ratio == 0, 1e-3, ratio))  # Replace zeros with small nonzero values

# Define custom COSMIC colors
cosmic_colors <- c("C→A" = "#0099CC",  # Blue
                   "C→G" = "#000000",  # Black
                   "C→T" = "#CC0000",  # Red
                   "T→A" = "#999999",  # Gray
                   "T→C" = "#66CC66",  # Green
                   "T→G" = "#FF9999")  # Pink

# Define colors for BRCA conditions
palette_colors <- c("BRCA +/+ (WT)" = "#377eb8",  # Blue
                    "BRCA -/+ (Het)" = "#e41a1c") # Red

# Main bar plot with transparency
p1 <- ggplot(mutation_counts, aes(x = trinucleotide_context)) +
  geom_col(aes(y = `BRCA +/+ (WT)`, fill = "BRCA +/+ (WT)"), alpha = 0.5, width = 0.7) +
  geom_col(aes(y = `BRCA -/+ (Het)`, fill = "BRCA -/+ (Het)"), alpha = 0.5, width = 0.7) +
  facet_grid(. ~ mutation_type, scales = "free_x", space = "free_x") +  # Align by mutation type
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.text.y = element_text(size = 12),
    strip.text = element_text(face = "bold", size = 14, color = "white", margin = margin(t = 5, b = 5)),  # Mutation type as x-axis title
    strip.background = element_rect(fill = "black"),  # Background for mutation labels
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",  # Keep legend for WT and Het comparison
    axis.title.x = element_blank()  # Removes x-axis label
  ) +
  labs(
    title = "Mutation Signature Comparison: BRCA +/+ (WT) vs. BRCA -/+ (Het)",
    y = "Total Count"  # Keep only y-axis label
  ) +
  scale_fill_manual(values = palette_colors)

# Print the plot
print(p1)

p1


# Define COSMIC colors, ensuring they match exactly with `mutation_type` in the dataset
cosmic_colors <- c("C→A" = "#0099CC",  # Blue
                   "C→G" = "#000000",  # Black
                   "C→T" = "#CC0000",  # Red
                   "T→A" = "#999999",  # Gray
                   "T→C" = "#66CC66",  # Green
                   "T→G" = "#FF9999")  # Pink

# Ensure `mutation_type` is a factor with the correct levels
mutation_counts <- mutation_counts %>%
  mutate(mutation_type = factor(mutation_type, levels = names(cosmic_colors)))  # Match factor levels

# Updated plot
p2 <- ggplot(mutation_counts, aes(x = trinucleotide_context, y = ratio, group = mutation_type, color = mutation_type)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  # Reference line
  geom_hline(yintercept = 2, linetype = "dotted", color = "gray50") +  # Upper threshold
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "gray50") +  # Lower threshold
  facet_grid(. ~ mutation_type, scales = "free_x", space = "free_x") +  # Align by mutation type
  scale_y_log10() +  # Log scale for better visibility
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.text.y = element_text(size = 12),
    strip.text = element_text(face = "bold", size = 14, color = "white", margin = margin(t = 5, b = 5)),  # Mutation type as x-axis title
    strip.background = element_rect(fill = "black"),  # Background for mutation labels
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank()  # Removes x-axis label
  ) +
  labs(
    title = "WT / Het Mutation Fold-Change",
    y = "WT / Het Ratio (Log Scale)"  # Keep only y-axis label
  ) +
  scale_color_manual(values = cosmic_colors)  # Assign correct colors

# Print the plot
print(p2)



final_plot <- p1 + p2 + plot_layout(ncol = 1, heights = c(1, 3))  # Give p1 more space

final_plot





# Strelka Analysis INDEL ----
# 1. ----

process_strelka_ind <- function(df) {
  df_NAlt<-data.frame(N=1:nrow(df@fix),Chr=df@fix[,1], Pos=as.numeric(df@fix[,2]),
                      NAlt=lengths(strsplit(df@fix[,5],",")) )
  
  AF<-extract.gt(df, element='FREQ', as.numeric=TRUE)
  DP<-extract.gt(df, element='DP', as.numeric=TRUE)
  
  VAF<-extract.gt(df, element='VAF', as.numeric=TRUE)
  
  
  format_field <- as.list(strsplit(df@fix[, "INFO"], ";"))
  
  sgt_values <- sapply(format_field, function(x) {
    # Find the element that starts with "SGT="
    sgt_entry <- grep("^SGT=", x, value = TRUE)
    # Remove the "SGT=" prefix to keep only the value
    if (length(sgt_entry) > 0) {
      sub("^SGT=", "", sgt_entry)
    } else {
      NA
    }
  })
  
  
  df_NAlt<-data.frame(N=1:nrow(df@fix),Chr=df@fix[,1], Pos=as.numeric(df@fix[,2]), Ref=df@fix[,4], Alt=df@fix[,5],
                      SGT = sgt_values, NAlt=lengths(strsplit(df@fix[,5],",")), Filter=df@fix[,"FILTER"])
  
  # Display the extracted SGT values
  unique(sgt_values)
  
  
  filtered_sgt_values <- sgt_values[!grepl("^(.+)->\\1$", sgt_values)]
  
  df_NAlt$Type = "Indel"
  df_NAlt$Size_Indel = 0
  for (i in 1:nrow(df_NAlt)) {
    
    ref = str_length(df_NAlt$Ref[i])
    alt = str_length(df_NAlt$Alt[i])
    
    if (ref < alt) {df_NAlt$Type[i] = "Ins"}
    if (ref > alt) {df_NAlt$Type[i] = "Del"}
    df_NAlt$Size_Indel[i] = abs(ref - alt)
    
  }
  
  return(df_NAlt)
  
  
  DP<-extract.gt(df, element='DP', as.numeric=TRUE) #this collects the DP for each
  TIR<-extract.gt(df, element='TIR')#, as.numeric=TRUE) #Reads filtered out
  TAR<-extract.gt(df, element='TAR')#, as.numeric=TRUE)
  TOR<-extract.gt(df, element='TOR')#, as.numeric=TRUE)

  
  extract_and_sum <- function(mat, sample_col) {
    raw <- mat[, sample_col]
    sum_vals <- sapply(raw, function(x) sum(as.numeric(strsplit(x, ",")[[1]])))
    list(raw = raw, sum = sum_vals)
  }
  
  
  # Fields to process
  fields <- c("TAR", "TIR", "TOR")
  
  # Loop through each field and assign raw and sum values
  for (field in fields) {
    normal_data <- extract_and_sum(get(field), "NORMAL")
    tumor_data  <- extract_and_sum(get(field), "TUMOR")
    
    df_NAlt[[paste0(field, "_Normal")]]     <- normal_data$raw
    df_NAlt[[paste0(field, "_Tumor")]]      <- tumor_data$raw
    df_NAlt[[paste0(field, "_Normal_sum")]] <- normal_data$sum
    df_NAlt[[paste0(field, "_Tumor_sum")]]  <- tumor_data$sum
  }
  
}


st_c2040_ind_df = process_strelka_ind(strelka_c2040_ind)
st_c2041_ind_df = process_strelka_ind(strelka_c2041_ind)

st_c2040_filt_ind_df = st_c2040_ind_df[st_c2040_ind_df$Filter == "PASS", ]
st_c2041_filt_ind_df = st_c2041_ind_df[st_c2041_ind_df$Filter == "PASS", ]


table_2040 <- table(st_c2040_filt_ind_df$Type)
table_2041 <- table(st_c2041_filt_ind_df$Type)

indel_table <- rbind(Sample_2040 = table_2040,
                     Sample_2041 = table_2041)
print(indel_table)

chisq.test(indel_table)


st_c2040_filt_ind_df$Sample <- "2040"
st_c2041_filt_ind_df$Sample <- "2041"

combined_df <- rbind(st_c2040_filt_ind_df, st_c2041_filt_ind_df)
boxplot(Size_Indel ~ Sample, data = combined_df,
        main = "Indel Size Distribution by Sample",
        ylab = "Indel Size", col = c("lightblue", "lightgreen"))


ggplot(combined_df, aes(x = Sample, y = Size_Indel, fill = Sample)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Indel Size Distribution by Sample", y = "Indel Size")

wilcox.test(Size_Indel ~ Sample, data = combined_df)

aggregate(Size_Indel ~ Sample, data = combined_df, median)
by(combined_df$Size_Indel, combined_df$Sample, summary)


ggplot(combined_df, aes(x = Sample, y = Size_Indel, fill = Sample)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  facet_wrap(~Type) +
  theme_minimal() +
  labs(title = "Indel Size Distribution by Type and Sample", y = "Indel Size")

ggplot(combined_df, aes(x = Sample, y = Size_Indel, fill = Sample)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  facet_wrap(~Type) +
  theme_minimal() +
  labs(title = "Indel Size by Type and Sample", y = "Indel Size")


del_df <- del_df %>% filter(grepl("^chr", Chr))
ins_df <- ins_df %>% filter(grepl("^chr", Chr))


del_counts <- del_df %>%
  group_by(Chr, Sample) %>%
  summarise(Count = n(), .groups = "drop")

ins_counts <- ins_df %>%
  group_by(Chr, Sample) %>%
  summarise(Count = n(), .groups = "drop")

ggplot(del_counts, aes(x = Chr, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Number of Deletions per Chromosome", y = "Count", x = "Chromosome") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(ins_counts, aes(x = Chr, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Number of Insertions per Chromosome", y = "Count", x = "Chromosome") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


chrom_order <- paste0("chr", c(1:22, "X", "Y"))
del_counts$Chr <- factor(del_counts$Chr, levels = chrom_order)

ggplot(del_counts, aes(x = Chr, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("2041" = "#D95F02", "2040" = "#1B9E77")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Deletions per Chromosome",
    x = "Chromosome",
    y = "Deletion Count",
    fill = "Sample"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Ensure chromosome order
ins_counts$Chr <- factor(ins_counts$Chr, levels = chrom_order)

ggplot(ins_counts, aes(x = Chr, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("2041" = "#D95F02", "2040" = "#1B9E77")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Insertions per Chromosome",
    x = "Chromosome",
    y = "Insertion Count",
    fill = "Sample"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )



# 3. Vep annotations needed ----

colnames(st_c2040_filtered) # only pass filter


# Define output file
vcf_file <- "output.vcf"

# Write the header
vcf_header <- c(
  "##fileformat=VCFv4.2",
  "##source=Generated_by_R=VarScan_Strelka.R",
  "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
)

# Extract relevant columns and create VCF data
vcf_data <- st_c2040_filtered %>%
  dplyr::mutate(
    QUAL = ".",
    ID = ".",
    INFO = paste0("DP=", DP_Tumor),  # Using DP_Tumor for Depth
    FILTER = ifelse(Filter == "PASS", "PASS", "FAIL")
  ) %>%
  dplyr::select(Chr, Pos, ID, Ref, Alt, QUAL, FILTER, INFO)

# Write to file
writeLines(vcf_header, vcf_file)
write.table(vcf_data, vcf_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)

cat("VCF file saved as", vcf_file, "\n")






# to be able to add vep annotations i need to create a vcf with the right format.

# ----- aplicar vep e ver impacto / predicted effect das mutações

# In Vep_strelka.sh
# echo "Running Vep..."
# VEP="/mnt/beegfs/apptainer/images/ensembl-vep_release_108.2.sif"
# CACHE="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/cache/"
# INPUT=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/clinvar/hcm_gtex_tovep.vcf
# output_vep=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/clinvar/hcm_gtex_tovep.vep.vcf
# apptainer exec "$VEP" vep --cache --dir "$CACHE" --canonical --vcf --force_overwrite --hgvs --af_gnomadg -i "$INPUT" -o "$output_vep" 


#### VarScan ----


unique(varscan_c2040[varscan_c2040$somatic_status != "Germline",]$chrom)
unique(varscan_c2041[varscan_c2041$somatic_status != "Germline",]$chrom)

unique
colnames(varscan_c2041)

unique(varscan_c2041$tumor_gt)
# add DP normal and DP tumor

varscan_c2040 = varscan_c2040 %>% mutate(DP_normal = normal_reads1 + normal_reads2,
                         DP_tumor = tumor_reads1 + tumor_reads2)

varscan_c2041 = varscan_c2041 %>% mutate(DP_normal = normal_reads1 + normal_reads2,
                         DP_tumor = tumor_reads1 + tumor_reads2)

### edit var_freq

varscan_c2040$tumor_var_freq <- as.numeric(gsub("%", "", varscan_c2040$tumor_var_freq)) / 100
varscan_c2040$normal_var_freq <- as.numeric(gsub("%", "", varscan_c2040$normal_var_freq)) / 100

varscan_c2041$tumor_var_freq <- as.numeric(gsub("%", "", varscan_c2041$tumor_var_freq)) / 100
varscan_c2041$normal_var_freq <- as.numeric(gsub("%", "", varscan_c2041$normal_var_freq)) / 100

### filter out normal_var freq 100% tumor_var 100%

varscan_c2040 = varscan_c2040[varscan_c2040$normal_var_freq != 1 & varscan_c2040$tumor_var_freq != 1, ]

#which(duplicated(varscan_c2040$position) == TRUE)
#which(duplicated(paste(varscan_c2040$chrom, varscan_c2040$position, sep = "_")))

#varscan_c2040$position[772928]
varscan_c2041 = varscan_c2041[varscan_c2041$normal_var_freq != 1 & varscan_c2041$tumor_var_freq != 1, ]


# VarScan ----

# Mutation Classification using Ref and Alt
classify_mutation <- function(ref, alt) {
  case_when(
    ref %in% c("A", "G") & alt %in% c("A", "G") ~ "Transition",  # A <-> G
    ref %in% c("C", "T") & alt %in% c("C", "T") ~ "Transition",  # C <-> T
    ref %in% c("A", "C", "G", "T") & alt %in% c("A", "C", "G", "T") ~ "Transversion",  # Any other change
    TRUE ~ NA_character_  # Handle missing or invalid cases
  )
}

varscan_c2040_filtered <- varscan_c2040 %>% filter(somatic_status == "Somatic")
varscan_c2041_filtered <- varscan_c2041 %>% filter(somatic_status == "Somatic")

# Apply to C2040
varscan_c2040_filtered <- varscan_c2040_filtered %>%
  mutate(MutationType = classify_mutation(ref, var))

# Apply to C2041
varscan_c2041_filtered <- varscan_c2041_filtered %>%
  mutate(MutationType = classify_mutation(ref, var))

# Count Transitions and Transversions
table(varscan_c2040_filtered$MutationType)
table(varscan_c2041_filtered$MutationType)



# Compare total mutation counts
cat("C2040 total mutations:", nrow(varscan_c2040_filtered), "\n")
cat("C2041 total mutations:", nrow(varscan_c2041_filtered), "\n")



mutations_per_chr <- function(df) {
  df %>% group_by(chrom) %>% summarise(Count = n()) %>% arrange(desc(Count))
}

varscan_c2040_chr <- mutations_per_chr(varscan_c2040_filtered)
varscan_c2041_chr <- mutations_per_chr(varscan_c2041_filtered)

# Plot mutation distribution across chromosomes
ggplot(bind_rows(varscan_c2040_chr %>% mutate(Sample="C2040"), 
                 varscan_c2041_chr %>% mutate(Sample="C2041")), 
       aes(x = chrom, y = Count, fill = Sample)) +
  geom_bar(stat="identity", position="dodge") +
  theme_minimal() +
  labs(title="Mutation Distribution by Chromosome", y="Mutation Count")




# other ----



ggplot(varscan_c2041, aes(x = tumor_var_freq)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(title = "Tumor Variant Frequency Distribution", x = "Tumor Variant Frequency", y = "Count")

ggplot(varscan_c2040, aes(x = tumor_var_freq)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(title = "Tumor Variant Frequency Distribution", x = "Tumor Variant Frequency", y = "Count")

## somatic_status != Germline e DP do normal > 7

ggplot(varscan_c2041[varscan_c2041$somatic_status != "Germline" & varscan_c2041$DP_normal > 7, ], aes(x = tumor_var_freq)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(title = "Tumor Variant Frequency Distribution", x = "Tumor Variant Frequency", y = "Count")

ggplot(varscan_c2040[varscan_c2040$somatic_status != "Germline" & varscan_c2040$DP_normal > 7, ], aes(x = tumor_var_freq)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(title = "Tumor Variant Frequency Distribution", x = "Tumor Variant Frequency", y = "Count")

library(stringr)
ggplot(
  varscan_c2041[
    str_starts(varscan_c2041$chrom, "chr") & 
      varscan_c2041$somatic_status != "Germline" & 
      varscan_c2041$DP_normal > 7, 
  ],
  aes(x = tumor_var_freq)
) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(
    title = "Tumor Variant Frequency Distribution",
    x = "Tumor Variant Frequency",
    y = "Count"
  ) +
  facet_wrap(~ chrom)


ggplot(varscan_c2041, aes(x = normal_var_freq, y = tumor_var_freq)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Normal vs Tumor Variant Allele Frequency", x = "Normal VAF", y = "Tumor VAF")

ggplot(varscan_c2041[varscan_c2041$somatic_status != "Germline" & varscan_c2041$DP_normal > 7 & varscan_c2041$DP_tumor > 7, ], aes(x = normal_var_freq, y = tumor_var_freq)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Normal vs Tumor Variant Allele Frequency", x = "Normal VAF", y = "Tumor VAF")

ggplot(varscan_c2041[varscan_c2041$somatic_status != "Germline" & varscan_c2041$DP_normal > 7 & varscan_c2041$DP_tumor > 7, ], aes(x = log10(DP_normal), y = log10(DP_tumor))) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Normal vs Tumor Variant Allele Frequency", x = "Normal VAF", y = "Tumor VAF")



ggplot(varscan_c2041[varscan_c2041$somatic_status != "Germline" & varscan_c2041$DP_normal > 7 & varscan_c2041$DP_tumor > 7, ], 
       aes(x = normal_var_freq, y = tumor_var_freq, color = chrom)) +
  geom_point(alpha = 0.6) +
  labs(title = "Tumor vs Normal Variant Frequencies by Chromosome", 
       x = "Normal Variant Frequency", 
       y = "Tumor Variant Frequency") +
  theme_minimal() +
  scale_color_viridis_d() +
  theme(legend.position = "right")


ggplot(varscan_c2040[varscan_c2040$somatic_status != "Germline" & varscan_c2040$DP_normal > 7 & varscan_c2040$DP_tumor > 7, ], 
       aes(x = normal_var_freq, y = tumor_var_freq, color = chrom)) +
  geom_point(alpha = 0.6) +
  labs(title = "Tumor vs Normal Variant Frequencies by Chromosome", 
       x = "Normal Variant Frequency", 
       y = "Tumor Variant Frequency") +
  theme_minimal() +
  scale_color_viridis_d() +
  theme(legend.position = "right")


ggplot(
  varscan_c2041[
    str_starts(varscan_c2041$chrom, "chr") & 
      varscan_c2041$somatic_status != "Germline" & 
      varscan_c2041$DP_normal > 7 & 
      varscan_c2041$DP_tumor > 7, 
  ],
  aes(x = normal_var_freq, y = tumor_var_freq)
) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Normal vs Tumor Variant Allele Frequency",
    x = "Normal VAF",
    y = "Tumor VAF"
  ) +
  facet_wrap(~ chrom)


ggplot(
  varscan_c2040[
    str_starts(varscan_c2040$chrom, "chr") & 
      varscan_c2040$somatic_status != "Germline" & 
      varscan_c2040$DP_normal > 7 & 
      varscan_c2040$DP_tumor > 7, 
  ],
  aes(x = normal_var_freq, y = tumor_var_freq)
) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Normal vs Tumor Variant Allele Frequency",
    x = "Normal VAF",
    y = "Tumor VAF"
  ) +
  facet_wrap(~ chrom)


ggplot(varscan_c2041, aes(x = somatic_status)) +
  geom_bar(fill = "purple") +
  labs(title = "Distribution of Somatic Status", x = "Somatic Status", y = "Count")

ggplot(varscan_c2041, aes(x = "Tumor Reads Strand Bias", y = tumor_reads1_plus - tumor_reads1_minus)) +
  geom_boxplot() +
  labs(title = "Strand Bias in Tumor Reads", x = "Strand", y = "Read Count Difference")

ggplot(varscan_c2041, aes(x = "Normal Reads Strand Bias", y = normal_reads1_plus - normal_reads1_minus)) +
  geom_boxplot() +
  labs(title = "Strand Bias in Normal Reads", x = "Strand", y = "Read Count Difference")


varscan_c2041$strand_bias <- varscan_c2041$tumor_reads1_plus - varscan_c2041$tumor_reads1_minus


ggplot(varscan_c2041, aes(x = log10(strand_bias + 1))) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Strand Bias in Tumor Reads", 
       x = "Strand Bias (Forward - Reverse Reads) + 1", 
       y = "Count") +
  theme_minimal()
