### Write a table of statistics for each map ###

library(dplyr)

# Read in the map data

old_map <- read.table("data/rdmap_consensus_v2.tab", header=TRUE)
new_map <- read.table("data/consensus_map.tab", header=TRUE)
AF_map <- read.table("data/AF_map.tab", header=TRUE)
AM_map <- read.table("data/AM_map.tab", header=TRUE)
BF_map <- read.table("data/BF_map.tab", header=TRUE)
BM_map <- read.table("data/BM_map.tab", header=TRUE)
A_map <- read.table("data/FamA_consensus.tab", header=TRUE)
B_map <- read.table("data/FamB_consensus.tab", header=TRUE)

# Read in the haplotype data
fama_haps <- read.table("data/hap_summary_a.txt", header=TRUE, fill=TRUE)
famb_haps <- read.table("data/hap_summary_b.txt", header=TRUE, fill=TRUE)

haps_all <- full_join(fama_haps, famb_haps, by="Locus")
haps_all <- mutate(haps_all, num.snps = ifelse(is.na(N_SNPs.x) & !is.na(N_SNPs.y), N_SNPs.y, ifelse(is.na(N_SNPs.y) & !is.na(N_SNPs.x), N_SNPs.x, ifelse(N_SNPs.x > N_SNPs.y, N_SNPs.x, N_SNPs.y))))
colnames(haps_all)[1] <- "locus"

maps <- list(new_map, A_map, B_map, AF_map, AM_map, BF_map, BM_map)

# Loop through each map and calculate stats

results <- list()
for (i in 1:length(maps)) {

  map <- maps[[i]] %>% mutate(loc.type = ifelse(grepl("Soc", locus, perl=TRUE), "msat", ifelse(grepl("Contig", locus, perl=TRUE), "snp", "est")))
  map <- left_join(map, haps_all, by="locus")

  map.res <- c()
  # Get the total number of loci
  map.res[1] <- nrow(map)

  # Get the number of each marker type (msats, ests, snps)
  map.res[2] <- nrow(filter(map, loc.type == "msat"))
  map.res[3] <- nrow(filter(map, loc.type == "est"))
  map.res[4] <- nrow(filter(map, loc.type == "snp"))

  # Get some stats for each linkage group
  by_lg <- group_by(map, lg)
  lg_stats <- summarise(by_lg, total.length = max(pos), num.loci = n_distinct(locus), num.snps = sum(num.snps, na.rm=TRUE))
  lg_stats <- mutate(lg_stats, mean.interval = total.length / num.loci)


  # Get total number of SNPs
  map.res[5] <- sum(lg_stats$num.snps)

  # Get mean lg size
  map.res[6] <- sprintf("%.3f", mean(lg_stats$total.length))

  # Get mean loci per group
  map.res[7] <- sprintf("%.3f", mean(lg_stats$num.loci))

  # Get mean marker interval
  map.res[8] <- sprintf("%.3f", mean(lg_stats$mean.interval))

  # Get total map length (cM)
  map.res[9] <- sum(lg_stats$total.length)

  results[[i]] <- map.res
}

rows <- c("Total Loci", "Msat Loci", "EST Loci", "RAD Loci", "SNP loci", "Mean LG Size", "Mean Loci per Group", "Mean Marker Interval", "Total Map Length")

table <- do.call(cbind, results)
table <- as.data.frame(cbind(rows, table))
colnames(table)[1:8] <- list("Stat", "Consensus", "Family A", "Family B", "AF", "AM", "BF", "BM")

write.table(table, "tables/raw/table_1.txt", quote=FALSE, row.names=FALSE, sep="\t")
