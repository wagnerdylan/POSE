# File:     geo_oribtal_plot.r
# Author:   Dylan R. Wagner
# Version:  2019-04-21
# Description:
#   Plots geocentric orbit cartesian coordinates.
#
# Usage:
#   rscript geo_oribtal_plot.r <FILE_PATH> <OUT_PATH> <ID> <TILE>
#

library(data.table)
library(ggplot2)
library(ggforce)

args <- commandArgs(trailingOnly = TRUE)
filepath <- args[1]
out_path <- args[2]
type_body <- args[3]
plot_title <- args[4]

data <- fread(filepath)
data <- subset(data, type == get("type_body"))
data <- data[, orbit_num := NA]

earth <- data.frame(x0=0, y0=0, r=6378137.0)

p <- ggplot(data, aes(x, y, color=z))
p <- p + geom_path()
p <- p + geom_hline(yintercept=0) + geom_vline(xintercept=0)
p <- p + theme_light()
p <- p + geom_circle(aes(x0 = x0, y0 = y0, r=r),
                     data=earth, inherit.aes = FALSE,
                     show.legend = FALSE)
p <- p + coord_fixed()
p <- p + labs(x="X (meters)", y="Y (meters)", title=get("plot_title"), subtitle=paste("Object:", type_body), color="Z (meters)")

ggsave(out_path, height=7, width=8)