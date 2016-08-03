#See https://cran.r-project.org/doc/contrib/Lemon-kickstart/kr_addat.html for help

chemotaxis_position_visualization_data_1 <- read.csv("~/Dropbox/MSU/Dufour/MABEplay/chemotaxis_position_visualization_data_1.csv", header=FALSE)
chemotaxis_position_visualization_data_2 <- read.csv("~/Dropbox/MSU/Dufour/MABEplay/chemotaxis_position_visualization_data_2.csv", header=FALSE)
chemotaxis_position_visualization_data_3 <- read.csv("~/Dropbox/MSU/Dufour/MABEplay/chemotaxis_position_visualization_data_3.csv", header=FALSE)
yax = range(c(chemotaxis_position_visualization_data_1$V2, chemotaxis_position_visualization_data_2$V2, chemotaxis_position_visualization_data_3$V2))
xax = range(c(chemotaxis_position_visualization_data_1$V1, chemotaxis_position_visualization_data_2$V1, chemotaxis_position_visualization_data_3$V1))
plot(chemotaxis_position_visualization_data_1$V1, chemotaxis_position_visualization_data_1$V2, col=1, xlab="X", ylab="Y", main="Paths for Generation 100", ylim=yax, xlim = xax)
par(new=T)
plot(chemotaxis_position_visualization_data_2$V1, chemotaxis_position_visualization_data_2$V2, axes=F, col=2, xlab="", ylab="", ylim=yax, xlim=xax)
par(new=T)
plot(chemotaxis_position_visualization_data_3$V1, chemotaxis_position_visualization_data_3$V2, axes=F, col=3, xlab="", ylab="", ylim=yax, xlim=xax)
abline(v=0)
abline(h=0)
