library(ggplot2)
library(reshape2)

verification_collision_data <- read.csv("~/Dropbox/MSU/Dufour/MABE_folders/MABE_source/World/ChemotaxisWorld/verification_collision_data.csv")

line_plt <- ggplot(verification_collision_data, aes(x = X, y = Y)) + geom_segment(xend = 0, yend = 0) + scale_x_continuous(limits = c(-12,12)) + scale_y_continuous(limits = c(-12,12))
point_plt <- ggplot(verification_collision_data, aes(x = X, y = Y)) + geom_point() + scale_x_continuous(limits = c(-12,12)) + scale_y_continuous(limits = c(-12,12))
