library(ggplot2)
library(grid)

#### map theme ####

# A custom ggplot theme for producing figures for the red drum map manuscript

map_theme <- ggplot2::theme(
  plot.title = element_text(size=16, face="bold", family="serif", vjust=1),
  panel.background = element_rect(fill = 'white'),
  panel.grid.major.y = element_line(size = 0.25, color = "lightgrey"),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.x = element_line(size = 0.25, color = "lightgrey"),
  panel.grid.minor.x = element_blank(),
  axis.line = element_line(size = 0.7, color = "black"),
  axis.ticks = element_line(size = 0.25, color = "black"),
  axis.text = element_text(color="black", size=15, vjust=0.5, family="serif"),
  axis.title.x = element_text(color="black", size=15, vjust=-0.35, face="bold"),
  axis.title.y = element_text(color="black", size=15, vjust=1, face="bold"),
  legend.key = element_rect(fill = "white"),
  panel.margin = unit(0.25, "lines"),
  strip.background = element_rect(fill="white")
)
