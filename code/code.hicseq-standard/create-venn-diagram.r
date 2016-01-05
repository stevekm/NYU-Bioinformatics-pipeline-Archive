#!/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
bed1_peaks <- as.integer(as.numeric(args[1]))
bed2_peaks <- as.integer(as.numeric(args[2]))
common_peaks <- as.integer(as.numeric(args[3]))
sample1 <- as.character(args[4])
sample2 <- as.character(args[5])
out <- args[6]

#total_peaks <- bed1_peaks + bed2_peaks - common_peaks

#bed2_lower <- bed1_peaks - common_peaks + 1

# load the required 
# package
library(VennDiagram)


output1 <- paste(out,"Venn_diagram.pdf",sep="/")
pdf(sprintf("%s", output1), height=5, width=5);

mat <- matrix(c(1,2,0,2), 2)
layout(mat, c(3.5,1), c(1,3))
par(mar=c(0.5, 4.5, 0.5, 0.5))

plot(1, type="n", axes=FALSE, xlab="", ylab="")
plot_colors <- c("cornflowerblue", "darkorchid1")
legend(x="top", inset=0, legend = c(sample1,sample2), col=plot_colors, pch=19, cex=0.55, horiz=FALSE)

par(mar=c(4.5, 4.5, 0.5, 0.5))
venn.plot <- draw.pairwise.venn(
        area1 = bed1_peaks,
        area2 = bed2_peaks,
        cross.area = common_peaks,
        category = c("", ""),
        lwd = 4,
        lty = "blank",
        fill = c("cornflowerblue", "darkorchid1"),
        label.col = "black",
        cex = 2,
        cat.cex = 1.0,
        cat.pos = c(5, -70),
        cat.dist = 0.05,
        cat.just = list(c(0.5, 0.5), c(0.5, 0.5)),
        ext.pos = 30,
        ext.dist = -0.05,
        ext.length = 0.9,
        ext.line.lwd = 2,
        ext.line.lty = "dashed"
    );

# Plot Venn
grid.draw(venn.plot);
dev.off();


# Generate the tiff as well
output2 <- paste(out,"Venn_diagram.tiff",sep="/")
tiff(sprintf("%s", output2), compression = "lzw", height=5, width=5, res=300, units="in");

mat <- matrix(c(1,2,0,2), 2)
layout(mat, c(3.5,1), c(1,3))
par(mar=c(0.5, 4.5, 0.5, 0.5))

plot(1, type="n", axes=FALSE, xlab="", ylab="")
plot_colors <- c("cornflowerblue", "darkorchid1")
legend(x="top", inset=0, legend = c(sample1,sample2), col=plot_colors, pch=19, cex=0.55, horiz=FALSE)

par(mar=c(4.5, 4.5, 0.5, 0.5))
venn.plot <- draw.pairwise.venn(
        area1 = bed1_peaks,
        area2 = bed2_peaks,
        cross.area = common_peaks,
        category = c("", ""),
        lwd = 4,
	lty = "blank",
        fill = c("cornflowerblue", "darkorchid1"),
        label.col = "black",
        cex = 2,
	cat.cex = 1.0,
	cat.pos = c(5, -70),
 	cat.dist = 0.05,
        cat.just = list(c(0.5, 0.5), c(0.5, 0.5)),
        ext.pos = 30,
        ext.dist = -0.05,
        ext.length = 0.9,
        ext.line.lwd = 2,
        ext.line.lty = "dashed"
    );

# Plot Venn
grid.draw(venn.plot);
dev.off();
