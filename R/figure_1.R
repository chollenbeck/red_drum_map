# Plot a 3-panel linkage map summary figure

library(ggplot2)
library(grid)
library(reshape2)
library(dplyr)

#### Input ####

# Read in the map data
map <- read.table("data/consensus_map.tab", header=TRUE)

# Read in the locus size data
sizes <- read.table("data/size_dist.txt", header=TRUE)

# Read in the annotations
annot <- read.csv("data/map_annot.csv", header=FALSE)

##### Data manipulations ####




# Add a column for locus type
map <- map %>% mutate(loc.type = ifelse(grepl("Soc", locus, perl=TRUE), "msat", ifelse(grepl("Contig", locus, perl=TRUE), "snp", "est")))

by_lg <- group_by(map, lg)
lg_stats <- summarise(by_lg, total.length = max(pos), num.loci = n_distinct(locus))
lg_stats <- mutate(lg_stats, mean.interval = total.length / num.loci)

map$loc.type <- factor(map$loc.type, levels = c("snp", "msat", "est"))


# Get the total number of bases mapped
map_sizes <- left_join(map, sizes, by="locus")
total_bases <- sum(map_sizes$size, na.rm = TRUE)
total_kb <- sprintf("%.1f", total_bases / 1000)


# Get the total number of loci and number of loci of each type
total_loci <- nrow(map)
total_rad <- nrow(filter(map, loc.type=="snp"))
total_msat <- nrow(filter(map, loc.type=="msat"))
total_est <- nrow(filter(map, loc.type=="est"))

# Get the proportion of annotated loci

# Filter the BLAST hits by evalue
annot_filt <- filter(annot, V11 < 1e-10)

# Filter out duplicates
annot_filt <- distinct(annot_filt, V1)

colnames(annot_filt)[1] <- "locus"
annot_map <- left_join(map, annot_filt, by="locus")
annotated <- inner_join(map, annot_filt, by="locus")
perc_annotated <- sprintf("%.2f", (nrow(annotated) / total_loci) * 100)


#### Linkage map plot ####

# Set up a data structure for scale tick marks

y <- seq(0, 120, by=10)
maj_ticks <- data.frame(y)
maj_ticks$x <- -0.5

y <- y <- seq(5, 115, by=10)
min_ticks <- data.frame(y)
min_ticks$x <- -0.5



# Plot with colors
mapplot <- ggplot(map, aes(x=lg, y=pos)) +

  geom_path(aes(group=lg), size=10, lineend="round", col="black") +
  geom_path(aes(group=lg), size=8, lineend="round", col="white", alpha=0.8) +
  geom_rect(data=subset(map, loc.type=="snp"), aes(xmin=lg-0.28, xmax=lg+0.28, ymin=pos-0.75, ymax=pos+0.75), fill='#247BA0') +
  geom_rect(data=subset(map, loc.type=="msat"), aes(xmin=lg-0.28, xmax=lg+0.28, ymin=pos-0.75, ymax=pos+0.75), fill='#50514F') +
  geom_rect(data=subset(map, loc.type=="est"), aes(xmin=lg-0.28, xmax=lg+0.28, ymin=pos-0.75, ymax=pos+0.75), fill='#F25F5C') +
  annotate("text", x = 0.25, y = seq(0, 125, by=10), label = seq(0, 125, by=10), color="#394242") +
  annotate("text", x = 0.25, y = -10, label = "cM", color="#394242") +
  annotate("text", x = 1:24, y = -13, label = 1:24, color="black", fontface="bold") +
  geom_segment(data=maj_ticks, aes(x=x, y=y, xend=0, yend=y), color="#394242") +
  geom_segment(data=min_ticks, aes(x=x, y=y, xend=-0.25, yend=y), color="#394242") +
  geom_segment(x=-0.5, y=5, xend=-0.5, yend=-125.5, size=2, color="#394242") +
  #facet_grid(.~lg) +
  theme(line=element_blank(), rect=element_blank(), axis.text = element_blank(), axis.title=element_blank()) +
  scale_fill_manual(name="Marker Type",
                    breaks=c("snp", "msat", "est"),
                    labels=c("RAD Locus", "Microsatellite", "EST-SSR")) +
  ylab("Total Length (cM)") +
  xlab("Linkage Group") +
  scale_y_reverse()


### Donut plot ####

# First make a data frame

by_type <- group_by(map, loc.type)
counts <- summarize(by_type, count=n())

# Produce some variables for the labels

counts$y.breaks <- cumsum(counts$count) - counts$count/2
counts$percent <- (counts$count / sum(counts$count))
counts$pos <- c(cumsum(360*counts$percent)-(360*counts$percent/2))
counts$pos <- c(ifelse(counts$pos<=180,counts$pos,counts$pos-180))


pie <- ggplot(counts, aes(x=factor(0.2), y=count, fill=loc.type)) +
  geom_bar(stat="identity", color="black", width=0.4) +
  coord_polar(theta='y') +
  geom_text(aes(x=1.5, y=cumsum(counts$count) - counts$count/2,
                label=paste(sprintf("%.0f", counts$percent * 100), '%', sep=''),
                angle=360), colour=c("#247BA0","#50514F", "#F25F5C"), size=6, fontface="bold") +
  map_theme +
  scale_fill_manual(values = c("#247BA0","#50514F", "#F25F5C"),
                    labels = c("RAD contigs", "Microsatellites", "EST-SSRs")) +
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=16, family="sans"),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.key = element_blank()
  )


### Text plot ###

# Create a table with text annotations

text <- c('total markers', 'RAD contigs', 'anonymous microsatellites', 'EST-SSRs', '% annotated', 'kilobases mapped')
stats <- c(total_loci, total_rad, total_msat, total_est, perc_annotated, total_kb)
color <- c('black', '#247BA0', '#50514F', '#F25F5C', '#AC80A0','#70C1B3')
ystart <- 6:1

figtext <- as.data.frame(cbind(text, stats, color, ystart))
figtext$label <- paste(figtext$stats, figtext$text, sep = ' ')

tplot <- ggplot(figtext) +
  geom_text(aes(label = label, x=0, y=ystart), colour = color, hjust="left", size=8, fontface="bold", family="sans") +
  xlim(-0.05, 1) +
  map_theme +
  theme(line=element_blank(), rect=element_blank(), axis.text = element_blank(), axis.title=element_blank())


#### Combined plot ####

tiff(filename = "figures/fig1.tif",
     width = 4800, height = 3000, units = "px", pointsize = 12,
     compression = "zip+p",
     bg = "white", res = 400, family = "", restoreConsole = TRUE,
     type = "cairo")

multiplot(pie, tplot, mapplot, layout=matrix(c(1,2,3,3), nrow=2, byrow=TRUE))

dev.off()
