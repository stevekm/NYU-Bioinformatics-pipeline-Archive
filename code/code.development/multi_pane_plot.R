m <- rbind(c(1, 1), c(2, 3))
print(m)
layout(m)
layout.show(3)
par(mar = c(3, 3, 0, 0))
for (i in 1:3) plot(1, 1, type = "n")

z<-rbind(c(1,1),c(2,3))
print(z)
layout(z)
layout.show(3)
par(mar = c(3, 3, 0, 0))
# for (i in 1:3) plot(1, 1, type = "n")
# plot(0,0,typ="n")
# frame()
legend("top",legend="Title")
plot(1,1,legend("top",legend="Title"),type=n)
layout.show(3)


###
# example provided by Sean Anderson
# run this in RStudio
x <- 13
x2 <- 8

temp.mat <- rep(NA, 5)
for(i in 1:3){
  for(j in 1:4){
    row <- seq(i, i + x, 3)
    temp.mat <- rbind(temp.mat, row)
  }
}
m.upper <- temp.mat[-1, ]

f1 <- seq(31, 31 + x2, 2)
f2 <- seq(32, 32 + x2, 2)
f3 <- seq(41, 41 + x2, 2)
f4 <- seq(42, 42 + x2, 2)

m.lower <- m.upper + 15
m <- rbind(f1, m.upper, f2,  f3, m.lower, f4)

layout(m)
# layout.show(m)
layout.show(max(m)) # shows all of the cells in the matrix as panes in the plotting area

# plots are placed in each pane in the order from 1 to maximum value plotted
# matrix entires with the same numeric value get the same plot;
# # adjacent matrix entries with the same value become a single plotting pane
# unused matrix entries remain empty, can be used as white-space between plots
# plots must be placed in order!

dev.off() # this clears the plot region

###
# its probably easier to create the matrix outside of R in a CSV in Excel, then import and turn into a matrix
tmp_matrix<-read.csv(file="/ifs/home/kellys04/projects/R_multipane_plot_matrix_template.csv",header = F)
tmp_matrix<-as.matrix(tmp_matrix)
tmp_matrix

# looks like this:
# > tmp_matrix
# V1 V2 V3 V4
# [1,]  1  1  1  1
# [2,]  2  2  3  3
# [3,]  2  2  3  3
# [4,]  2  2  3  3
# [5,]  2  2  3  3
# [6,]  4  4  5  5
# [7,]  4  4  5  5
# [8,]  4  4  5  5
# [9,]  4  4  5  5

# recreate this object with this code:
# > dput(tmp_matrix)

# create a matrix
# its probably easier to create outside of R in a CSV in Excel, then import and turn into a matrix
tmp_matrix<-structure(c(1L, 2L, 2L, 2L, 2L, 4L, 4L, 4L, 4L, 1L, 2L, 2L, 2L, 
                        2L, 4L, 4L, 4L, 4L, 1L, 3L, 3L, 3L, 3L, 5L, 5L, 5L, 5L, 1L, 3L, 
                        3L, 3L, 3L, 5L, 5L, 5L, 5L), .Dim = c(9L, 4L), .Dimnames = list(
                          NULL, c("V1", "V2", "V3", "V4")))
tmp_matrix

# pass the matrix to layout() to set the plot panels
layout(tmp_matrix) 
layout.show(max(tmp_matrix)) # show all of the panes in the layout
# plots are placed in each pane in the order from 1 to maximum value plotted
# matrix entires with the same numeric value get the same plot;
# # adjacent matrix entries with the same value become a single plotting pane
# unused matrix entries remain empty, can be used as white-space between plots
# the entire multi-panel layout will expand to fill the entire plotting area
# plots will be placed in numeric order

par(mar=c(0,0,0,0)) # need to set this for some reason
plot(1,type='n',axes=FALSE,xlab="",ylab="") # call blank plot to fill the first panel

# set up the Legend in the first panel
legend("top",legend="This is the legend")

# fill in the rest of the plots
for (i in 2:5){
  par(mar=c(2,2,1,1))
  plot(1,1,type="n")
}

# dev.off() # this clears the plot region, and allows par() to be reset on the next plot call

#######

require(VennDiagram)
library(ggplot2)
library(grid)
library(gridExtra)

A<- draw.pairwise.venn(600, 200, 61, c("X", "Y"), 
                       col= rep("gray70", 2), lwd= rep (1, 2), fill= c("skyblue1", "yellowgreen"), 
                       cat.pos=0, fontfamily = rep("sans"), cat.fontfamily= rep("sans"), sacled=FALSE)

B<- draw.pairwise.venn(400, 200, 60, c("X", "Y"), col= rep("gray70", 2), 
                         lwd= rep (1, 2), fill= c("skyblue1", "yellowgreen"), cat.pos=0,
                         fontfamily = rep("sans"), cat.fontfamily= rep("sans"), scaled=FALSE)

C<- draw.pairwise.venn(700, 500, 75, c("X", "Y"), col= rep("gray70", 2), lwd= rep (1, 2), 
                         fill= c("skyblue1", "yellowgreen"), cat.pos=0, fontfamily = rep("sans"), 
                         cat.fontfamily= rep("sans"), scaled=FALSE)



pushViewport(viewport(layout=grid.layout(ncol=3, widths = unit(rep(1/3,3), "npc"))))
pushViewport(viewport(layout.pos.col=1))
draw.pairwise.venn(600, 200, 61, c("X", "Y"), col= rep("gray70", 2), lwd= rep (1, 2), 
                   fill= c("skyblue1", "yellowgreen"), cat.pos=0, fontfamily = rep("sans"), 
                   cat.fontfamily= rep("sans"), sacled=FALSE)
popViewport()
pushViewport(viewport(layout.pos.col=2))
draw.pairwise.venn(400, 200, 60, c("X", "Y"), col= rep("gray70", 2), lwd= rep (1, 2), 
                   fill= c("skyblue1", "yellowgreen"), cat.pos=0, fontfamily = rep("sans"), 
                   cat.fontfamily= rep("sans"), scaled=FALSE)
popViewport()
pushViewport(viewport(layout.pos.col=3))
draw.pairwise.venn(700, 500, 75, c("X", "Y"), col= rep("gray70", 2), lwd= rep (1, 2), 
                   fill= c("skyblue1", "yellowgreen"), cat.pos=0, fontfamily = rep("sans"), 
                   cat.fontfamily= rep("sans"), scaled=FALSE)
popViewport(0)

#  grid.arrange(gTree(children=A), gTree(children=B), gTree(children=C), ncol=3, main="Venn")




