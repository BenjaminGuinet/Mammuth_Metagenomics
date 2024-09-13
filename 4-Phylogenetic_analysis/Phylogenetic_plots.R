# Load required libraries
library(ggtree)
library(ggplot2)
library(phytools)
library(dplyr)
library(stringr)
library(ggimage)
library(ape)
library(ggrepel)


setwd("/Users/benjguin/Desktop/All_trees/")
csv_files <- list.files(pattern = "\\.csv$")

# Read the first CSV file to get column names
first_df <- read.csv(csv_files[1], header = TRUE)

# Sort the column names
sorted_columns <- colnames(first_df)

# Initialize an empty list to store dataframes
dataframes_list <- list()

# Loop through each CSV file, read it, reorder columns, remove "index" column if exists, and store it in the list
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change options if needed
  
  # Reorder columns
  data <- data[, sorted_columns, drop = FALSE]
  print(colnames(data))
  # Remove "index" column if exists
  if ("index" %in% colnames(data)) {
    data <- data[, colnames(data) != "index", drop = FALSE]
  }
  
  dataframes_list[[file]] <- data
}

# Check if all dataframes have the same column names and order
consistent_columns <- all(sapply(dataframes_list, identical, x = colnames(dataframes_list[[1]])))

if (!consistent_columns) {
  stop("Column names or order are not consistent across all dataframes.")
}

# Merge all dataframes into a single dataframe
Metadata_df <- do.call(rbind, dataframes_list) 

Mammuth_metadata_df <- read.csv("/Users/benjguin/Desktop/Mammuth_phylogenies_eye/All_trees/Working_table_merged_damage_results_with_Kraken.txt",header = TRUE, sep=";")
Mammuth_metadata_df$Bam_file_name3 <- gsub("_extracted.bam", "", Mammuth_metadata_df$Bam_file_name)
Mammuth_metadata_df$Bam_file_name3 <- gsub("sorted_", "", Mammuth_metadata_df$Bam_file_name3)

Tree_dir <- "/Users/benjguin/Desktop/All_trees/"


# List files in directory
tree_files <- list.files(path = Tree_dir, pattern = "\\.treefile$", full.names = TRUE)

# Identify indices of elements containing "trimmed"
indices_to_remove <- grep("trimmed", tree_files)

# Remove elements containing "trimmed"
tree_files <- tree_files[-indices_to_remove]


#####. 
library(ggtree)
library(cowplot)

ggt_list <- list()
nb_trees<-0
# Loop over each tree file
for (file in tree_files) {
  nb_trees <- nb_trees+1
  print(nb_trees)
  newick_file <- basename(file)
  
  tree <- read.tree(paste0(Tree_dir, newick_file))
  
  tree$tip.label <- gsub("Mammuthus-FK033_ALL_new_Bisgaard", "Mammuthus-FK033_Pasteurella_multocida_9vssfquZ7Z", tree$tip.label)
  tree$tip.label <- gsub("Mammuthus-MD024_ALL_new_Bisgaard", "Mammuthus-MD024_Pasteurella_multocida_9vssfquZ7Z", tree$tip.label)
  
  #tree <- read.tree("/Users/benjguin/Desktop/All_trees/Basfia_output_sequences.fasta.treefile")
  
  i <- 0
  for (label in tree$tip.label) {
    i <- i + 1
    label <- str_replace(str_replace(label, "_extracted_.*", ""), ".*sorted_", "")
    
    if (grepl("Mammuthus", label)) {
      print("Mammuthus")
      print(label)
      Region <- Mammuth_metadata_df[Mammuth_metadata_df$Bam_file_name3 == label, ]$Region
      Tissue <- Mammuth_metadata_df[Mammuth_metadata_df$Bam_file_name3 == label, ]$Material2
      Sample_name <- Mammuth_metadata_df[Mammuth_metadata_df$Bam_file_name3 == label, ]$Name 
      Type <- Mammuth_metadata_df[Mammuth_metadata_df$Bam_file_name3 == label, ]$Type
      Age <- Mammuth_metadata_df[Mammuth_metadata_df$Bam_file_name3 == label, ]$X14C.cal
      
      if (Type=="Mammuthus"){
        print("")
      } else {
        Sample_name<- Type
        print (Type)
      }
      
      Score <- Mammuth_metadata_df[Mammuth_metadata_df$Bam_file_name3 == label, ]$Ancient_score
      Country <- Mammuth_metadata_df[Mammuth_metadata_df$Bam_file_name3 == label, ]$Country
      
      if (grepl("Mammuthus-P17|Mammuthus-P23|Mammuthus-P30|Mammuthus-P21", label)){
        Tissue <- "Tooth"
        Region <- "Wrangel"
        Country <- "Russia"
      }
      
      New_label <- paste0(Sample_name, " | ", Region, " (", Country ,")", " | ", Tissue, " | ", Age, " [", Score, "]")
      
      if (nchar(New_label)[1] == 6 ){
        New_label <- label
      }
      
    } else {
      
      print("Bacteria")
      print(label)
      
      host <- Metadata_df[grepl(label, Metadata_df$Accession), ]$host
      isolation_source <- Metadata_df[grepl(label, Metadata_df$Accession), ]$isolation.source
      species <- Metadata_df[grepl(label, Metadata_df$Accession), ]$Organism
      accession <- label
      
      # Replace values based on conditions
      host <- ifelse(is.na(host) | nchar(host) > 20 | tolower(host) %in% c("missing", "available"), "missing", host)
      isolation_source <- ifelse(is.na(isolation_source) | nchar(isolation_source) > 20 | tolower(isolation_source) %in% c("missing", "available"), "missing", isolation_source)
      
      New_label <- paste0(species, " (", accession ,") | ", host," | ", isolation_source)
      #New_label <- paste0(species, " | ", host," | ", isolation_source)
      
      if (nchar(New_label)[1] == 6 ){
        New_label <- label
      }
    }
    tree$tip.label[i] <- New_label
  }
  
  labels_to_root <- c(
    "Achromobacter xylosoxidans (GCA_016728825) | Missing_host | missing",
    "Missing (GCA_900248205) | missing | missing",
    "Flavobacterium terrigena (GCA_900108955) | Not_applicable | missing",
    "Flavobacterium lacus (GCA_003268815) | Not_applicable | missing",
    "Moraxella lacunata (GCA_900453245) | Human_(child) | eye",
    "Bacillus cereus (GCA_002220285) | Salad | Salad",
    "Corynebacterium glutamicum SCgG2 (GCA_000404185) | Soil | soil",
    "Paenochrobactrum gallinarii (GCA_014205685) | Not_applicable | missing",
    "Actinobacillus pleuropneumoniae (GCA_900638445) | Pig | missing"
  )
  
  tree <- midpoint.root(tree)
  
  
  for (label in labels_to_root) {
    if (label %in% tree$tip.label) {
      # Find the common ancestor of the specified label
      # Root the tree on the common ancestor
      print(paste0("rooted to ", label))
      tree <- root(tree, outgroup = label,resolve.root = TRUE)
      # Plot the rooted tree
    } else {
      cat("Label", label, "not found in the tree.\n")
    }
  }
  
  
  if ("Haemophilus influenzae (GCA_000931575) | Homo_sapiens | missing" %in% tree$tip.label){
    tree <- root(tree, outgroup = "Haemophilus influenzae (GCA_000931575) | Homo_sapiens | missing",resolve.root = TRUE)
  }
  # Create a ggtree object with adjusted size
  ggt <- 
    ggtree(tree, size = 0.25) + 
    #geom_tippoint(aes(color = case_when(
    # grepl("Wrangel", label) ~ "red",
    # grepl("Russia", label) ~ "green",
    # grepl("Canada", label) ~ "blue",
    # TRUE ~ "black"
    #)), size = 0.2) +
    geom_tiplab(size = 1, as_ylab = TRUE, align = TRUE,
                linetype = "dotted", linesize = 0.1) +  # Adjust label size
    geom_treescale(y = -0.5, fontsize = 1, linesize = 0.2) +
    scale_color_manual(values = c("red", "green", "blue", "black")) +
    theme(legend.position = "none") +
    geom_text_repel(data = subset(fortify(tree), isTip == FALSE), aes(label = label, x = x, y = y), 
                    size = 0.5,
                    fontface = 'bold',
                    box.padding = 0.05,
                    # Add extra padding around each data point.
                    point.padding = 0.05)  
  # Store the ggplot object in the list
  ggt_list[[length(ggt_list) + 1]] <- ggt
  
}

# Create a grid of plots
ggt_combined <- cowplot::plot_grid(plotlist = ggt_list, ncol = 4, rel_heights = rep(2, length(ggt_list)), rel_widths = rep(1, length(ggt_list)))
# Save the combined plot to a PDF file
pdf_filename <- "/Users/benjguin/Desktop/All_trees/combined_phylogeny.pdf"
ggsave(pdf_filename, ggt_combined, device = "pdf",height = 5, width = 5)
