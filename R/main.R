# The main analysis script for the red drum linkage map synteny analysis

# Create the figures
system("Rscript R/figure_1.R")
system("Rscript R/figure_2.R")
system("Rscript R/figure_3.R")

# Create the tables
system("Rscript R/table_1.R")
system("Rscript R/table_2.R")




