### Write a table for syntenic block stats ###
library(plyr)
library(dplyr)
library(tidyr)


# Read in the map data and the synteny data

map <- read.table("data/consensus_map.tab", header=TRUE)
data <- read.table("data/all.blocks.tab", header=TRUE)

# Read in a list of chromosomes and lengths for stats

chr_list <- read.table("data/chr_lengths.txt", header=TRUE)

# Get the total number of bases for each species

chr_by_spec <- group_by(chr_list, species)
total_bases <- dplyr::summarise(chr_by_spec, total_bases = sum(length))
total_bases <- dplyr::rename(total_bases, COMP_SPECIES = species)


# Get the number of single hits for each species

syn_all <- read.csv("data/synteny.out", header=TRUE)

species <- c("gacu", "onil", "trub", "tnig", "dlab", "lcal")

syn_list = list()
for (i in 1:length(species)) {
  subset <- select(syn_all, 1:4, starts_with(species[i]))
  col <- paste(species[i], '_pos', sep='')
  subset <- na.omit(subset)
  colnames(subset)[5:7] <- c("chr", "prop", "pop")
  subset$species <- species[i]
  syn_list[[i]] <- subset
}

hit_data <- ldply(syn_list, data.frame)

by_spec <- group_by(hit_data, species)

hit_sum <- dplyr::summarise(by_spec, num.loci = n())
colnames(hit_sum)[1] <- "COMP_SPECIES"

# Get some stats for each linkage group
by_lg <- group_by(map, lg)
lg_stats <- dplyr::summarise(by_lg, total.length = max(pos), num.loci = n_distinct(locus),na.rm=TRUE)
lg_stats <- mutate(lg_stats, mean.interval = total.length / num.loci)
tot_map_length <- sum(lg_stats$total.length)


by_spec <- group_by(data, COMP_SPECIES)

block_stats <- dplyr::summarise(by_spec, num_blocks = n(), loc_per_block = sprintf("%.3f", mean(NUM_LOCI)), mean_comp_size = mean(COMP_SIZE), mean_map_size = sprintf("%.3f", mean(MAP_SIZE)), tot_comp_size = sum(COMP_SIZE), tot_map_size = sum(MAP_SIZE))


block_stats <- left_join(block_stats, total_bases, by="COMP_SPECIES")
block_stats <- left_join(block_stats, hit_sum, by="COMP_SPECIES")

# Get the proportion of the comparison genome in syntenic blocks

block_stats <- mutate(block_stats, prop_comp = sprintf("%.3f", tot_comp_size / total_bases), prop_map = sprintf("%.3f", tot_map_size / tot_map_length))

# Put the block size stats in terms of megabases

block_stats <- mutate(block_stats, mean_comp_mb = sprintf("%.1f", mean_comp_size / 1000000), tot_comp_mb = sprintf("%.1f", tot_comp_size / 1000000))


# Arrange the table by similarity to red drum

block_stats <- arrange(block_stats, desc(num_blocks))

# Change the species to the full name and add a common name column

block_stats$common_name <- block_stats$COMP_SPECIES

block_stats$COMP_SPECIES <- gsub("gacu", "Gasterosteus aculatus", block_stats$COMP_SPECIES)
block_stats$COMP_SPECIES <- gsub("onil", "Orecochromis niloticus", block_stats$COMP_SPECIES)
block_stats$COMP_SPECIES <- gsub("dlab", "Dicentrarchus labrax", block_stats$COMP_SPECIES)
block_stats$COMP_SPECIES <- gsub("lcal", "Lates calcarifer", block_stats$COMP_SPECIES)
block_stats$COMP_SPECIES <- gsub("trub", "Takifugu rubripes", block_stats$COMP_SPECIES)
block_stats$COMP_SPECIES <- gsub("tnig", "Tetraodon nigroviridis", block_stats$COMP_SPECIES)

block_stats$common_name <- gsub("gacu", "Stickleback", block_stats$common_name)
block_stats$common_name <- gsub("onil", "Nile tilapia", block_stats$common_name)
block_stats$common_name <- gsub("dlab", "European seabass", block_stats$common_name)
block_stats$common_name <- gsub("lcal", "Barramundi", block_stats$common_name)
block_stats$common_name <- gsub("trub", "Fugu", block_stats$common_name)
block_stats$common_name <- gsub("tnig", "Green spotted puffer", block_stats$common_name)

# Subsample the columns

block_stats <- select(block_stats, COMP_SPECIES, common_name, num.loci, num_blocks, loc_per_block, mean_comp_mb, mean_map_size, tot_comp_mb, tot_map_size, prop_comp, prop_map)


colnames(block_stats) <- c("Species", "Common Name", "BLAST Hits", "Number of Blocks", "Mean Loci Per Block", "Mean Comparison Block Size", "Mean Map Block Size", "Total Comparison Block Size", "Total Map Block Size", "Proportion of Genome in Blocks", "Proportion of Map in Blocks")


write.table(block_stats, "tables/raw/table_2.txt", quote=FALSE, row.names=FALSE, sep="\t")
