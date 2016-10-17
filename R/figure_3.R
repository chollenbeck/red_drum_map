# Plot syntenic blocks and synteny-mapped loci on LG22

library(dplyr)
library(ggplot2)
library(rentrez)
library(tidyr)
library(seqinr)

#### Data input ####

# Linkage map
map <- read.table("data/consensus_map.tab", header=TRUE)

# Syntenic blocks
data <- read.table("data/all.blocks.tab", header=TRUE)

# Locations for synteny mapped loci
syn <- read.table("data/all.synteny.mapped.out", header=TRUE)

# Annotations for synteny mapped loci
nr_vert <- read.table("data/blastx_nr.outfmt6", header=FALSE)

# Define some colors
col <- c('#50514F','#F25F5C','#AC80A0','#247BA0','#70C1B3','#41ab5d')


# Filter the annotations for only loci with evalues < 1E-10

nr_vert <- filter(nr_vert, V11 < 1E-10)
nr_vert_by_loc <- group_by(nr_vert, V1)
nr_vert_unique <- top_n(nr_vert_by_loc, 1)
nrow(nr_vert_unique)

# Parse out the ids from the ID column

id_table <- select(nr_vert_unique, V2)
id_table <- separate(id_table, V2, c("x", "gi", "y", "accession", "z"), sep="[|]")
id_table <- select(id_table, V1, gi, accession)

# Fetch the sequences from NCBI to get the protein names

res <- entrez_fetch("protein", id=id_table$gi, rettype="fasta")

temp <- tempfile()
write(res, temp)
parsed_res <- seqinr::read.fasta(file=temp, seqtype="AA")


names = list()
for (i in 1:length(parsed_res)) {
  items <- strsplit(attr(parsed_res[[i]], "Annot"), "[|]")
  prot_name <- items[[1]][[5]]
  prot_name <- gsub("^\\s+|\\s+$", "", prot_name)
  names[i] <- prot_name
}

names.vec <- unlist(names)

# Add the names to the table

id_table$name <- names.vec

# Merge the tables

colnames(id_table)[1] <- "LOCUS"
locus_table <- left_join(syn, id_table, by="LOCUS")


# Fix the number of reported digits for the interval size

locus_table$INTERVAL_SIZE <- sprintf("%.2f", locus_table$INTERVAL_SIZE)


# Take the subset of the data for LG22

lg22 <- subset(locus_table, LG == 22)

# Arrange the table for easy plotting


lg22 <- mutate(lg22, lowrange = ifelse(LEFT_POS <= RIGHT_POS, LEFT_POS, RIGHT_POS),
                      highrange = ifelse(RIGHT_POS >= LEFT_POS, RIGHT_POS, LEFT_POS))

lg22 <- arrange(lg22, lowrange)

lg22$name[is.na(lg22$name)] <- "unknown transcript"

lg22$name <- gsub("PREDICTED: ", '', lg22$name)
lg22$name <- gsub("\\[.*\\]", '', lg22$name, perl=TRUE)

# Add height column for each comparison species that determines how the blocks are positioned vertically on the map

species <- c("gacu", "onil", "trub", "tnig", "dlab", "lcal")

base <- 0.15
step <- 0.025
for (i in 1:length(species)) {
  value <- (base + ((i-1) * step))
  data$plot.height[data$COMP_SPECIES == species[i]] <- value
}

# Make a series of line segments for delineating gene intervals
base <- 0.01
step <- 0.011
for (i in 1:nrow(lg22)) {
  value <- base + ((i-1) * step)
  lg22$line.height[i] <- value
}

# Get a dataframe for loci representing mapped loci
lg22_markers <- filter(map, lg==22)


tiff(filename = "figures/fig3.tif",
     width = 4800, height = 3000, units = "px", pointsize = 12,
     compression = "zip+p",
     bg = "white", res = 400, family = "", restoreConsole = TRUE,
     type = "cairo")

ggplot(subset(data, MAP_LG == 22)) +
  geom_segment(x=0, xend=max(lg22_markers$pos), y=-0.01, yend=-0.01, size=3, lineend="round") +
  geom_segment(data=lg22_markers, aes(x=pos, xend=pos, y=-0.01, yend=0), size=2, lineend="round") +
  geom_rect(aes(xmin=MAP_START, xmax=MAP_STOP, ymin=plot.height - 0.01, ymax=plot.height + 0.01, fill=COMP_SPECIES), linetype="solid", size=1, colour="black", alpha=1.0) +
  geom_segment(data=lg22, aes(x = lowrange, xend=highrange, y = line.height, yend=line.height), size=1.5, colour="#984447") +
  geom_text(data=subset(lg22, lowrange==highrange), aes(label = '*', x= lowrange, y = line.height), size=8, colour="#984447") +
  geom_text(data=lg22, aes(label = name, x= lowrange - 1, y = line.height), hjust="right") +
  ylim(-0.01,0.35) +
  map_theme +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="italic", size=12)
        ) +
  scale_fill_manual(breaks=species,
                    values=col,
                    labels = c("G. aculeatus", "O. niloticus", "T. rubripes", "T. nigroviridis", "D. labrax", "L. calcarifer") ) +
  xlab("Linkage Map Position (cM)")


dev.off()

