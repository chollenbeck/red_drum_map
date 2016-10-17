library(dplyr)

# Read in the map data and the synteny data

map <- read.table("data/consensus_map.tab", header=TRUE)
data <- read.table("data/all.blocks.tab", header=TRUE)
syn <- read.table("data/all.synteny.mapped.out", header=TRUE)

# Write a karyotype file

# Get the length of each LG
by_lg <- group_by(map, lg)
by_lg <- top_n(by_lg, 1, pos)

kary <- ungroup(by_lg)
kary <- select(kary, lg, pos)
kary$color <- "black"
kary <- mutate(kary, id = paste("lg", lg, sep=''), label = paste("LG", lg, sep=''),
               chr = "chr", dash = "-", start = "0", end = pos * 500000, color = "grey")

kary_final <- kary[,c(6, 7, 4, 5, 8, 9, 3)]

write.table(kary_final, file = "data/circos/drum_karyotype.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Write a markers file for the linkage map

markers <- mutate(map, chr = paste("soc", lg, sep=''), start = as.integer(pos * 500000), end = as.integer(pos * 500000))
markers_final <- markers[,c(4, 5, 6, 1)]

write.table(markers_final, file = "data/circos/markers.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Write a tile file for each of the comparison species

# Get a list of the comparison species

species <- as.vector(unique(data$COMP_SPECIES))


for (i in 1:length(species)) {
  spec <- filter(data, COMP_SPECIES == species[i])
  spec <- mutate(spec, chr = paste("lg", MAP_LG, sep=''), start = as.integer(MAP_START * 500000), end = as.integer(MAP_STOP * 500000))
  tiles <- select(spec, chr, start, end)
  write.table(tiles, file = paste("data/circos/", species[i], ".txt", sep = ''), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Write a tile file for the synteny mapped loci

syn_loci <- mutate(syn, chr = paste("lg", LG, sep=''), start = as.integer(LEFT_POS * 500000), end = as.integer(RIGHT_POS * 500000))
loci <- select(syn_loci, chr, start, end)
write.table(loci, file = "data/circos/syn_loci.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# Try to run circos from R

system('perl "scripts/circos-0.66/bin/circos" -conf data/circos/circos.conf -outputdir figures -outputfile fig2')
