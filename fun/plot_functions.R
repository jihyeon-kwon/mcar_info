##################
# Plot Functions #
##################

library(RColorBrewer)




# colors:

cols0 <- brewer.pal(3, "Dark2")
cols <- paste(cols0, "50", sep = "")

## I HISTOGRAMS ----------------------------------------------------------------

# 1. One Histogram -------------------------------------------------------------


hist.one <- function(data = data,
                     legend.loc = "topright",
                     legend = "",
                     xlab = "Estimates",
                     ylab = "Density",
                     mar = c(4, 4, 1, 1) + .2,
                     cex = 1,
                     cex.lab = 1,
                     freq = FALSE,
                     ylim = ylim,
                     length.out = 100,
                     breaks = seq(min(data, na.rm = TRUE),
                                  max(data, na.rm = TRUE),
                                  length.out = length.out
                     ),
                     vline = "") {
  par(mar = mar, cex = cex, cex.lab = cex.lab)
  dens <- hist(data, breaks = breaks, plot = FALSE)$density
  ymax <- max(dens) * 1.3
  ylim <- c(0, ymax)

  hist(data,
    col = cols[1],
    border = FALSE,
    breaks = breaks,
    freq = freq,
    xlab = "Estimates",
    main = "",
    ylim = ylim
  )

  box()

  # true value
  lines(
    x = c(vline, vline),
    y = c(0, ylim[2])
  )
  
  # median
  lines(
    x = c(median(data), median(data)),
    y = c(0, ylim[2]),
    col = cols0[1]
  )
  
  legend(legend.loc,
    legend = legend,
    fill = cols[1],
    border = FALSE,
    bty = "n",
    cex = 4 / 5
  )
}



# mtesting

# data <- rnorm(10000, mean = 1, sd = 1)
# hist.one(data,
#          legend = "data",
#          vline = 3)





# 2. Two Histograms -------------------------------------------------------



hist.two <- function(data1 = data1,
                     data2 = data2,
                     legend.loc = "topright",
                     legend1 = "",
                     legend2 = "",
                     xlab = "Estimates",
                     ylab = "Density",
                     mar = c(4, 4, 1, 1) + .2,
                     cex = 1,
                     cex.lab = 1,
                     freq = FALSE,
                     ylim = ylim,
                     length.out = 100,
                     vline = "",
                     breaks = seq(min(data1, data2, na.rm = TRUE),
                                  max(data1, data2, na.rm = TRUE),
                                  length.out = length.out
                     )) {
  par(mar = mar, cex = cex, cex.lab = cex.lab)
  dens1 <- hist(data1, breaks = breaks, plot = FALSE)$density
  dens2 <- hist(data2, breaks = breaks, plot = FALSE)$density
  ymax <- max(dens1, dens2) * 1.3
  ylim <- c(0, ymax)

  # first histogram
  hist(data1,
    breaks = breaks,
    col = cols[1],
    border = FALSE,
    freq = freq,
    xlab = "Estimates",
    main = "",
    ylim = ylim
  )
  
  # first median
  lines(
    x = c(median(data1, na.rm = TRUE), median(data1, na.rm = TRUE)),
    y = c(0, ylim[2]),
    col = cols0[1]
  )
  
  # second histogram
  hist(data2,
    breaks = breaks,
    col = cols[2],
    border = FALSE,
    freq = freq,
    ylim = ylim,
    add = TRUE
  )
  
  # second median
  lines(
    x = c(median(data2, na.rm = TRUE), median(data2, na.rm = TRUE)),
    y = c(0, ylim[2]),
    col = cols0[2]
  )

  box()

  lines(
    x = c(vline, vline),
    y = c(0, ylim[2])
  )

  legend(legend.loc,
    legend = c(legend1, legend2),
    fill = cols[1:2],
    border = FALSE,
    bty = "n",
    cex = 4 / 5
  )
}



# testing

# data1 <- rnorm(10000, mean = 1, sd = 1)
# data2 <- rnorm(10000, mean = 3, sd = 1)
# hist.two(data1 = data1,
#              data2 = data2,
#              legend1 = "Data1",
#              legend2 = "Data2",
#              vline = 2)





# 3. Three Histograms -------------------------------------------------------



hist.three <- function(data1 = data1,
                       data2 = data2,
                       data3 = data3,
                       legend.loc = "topright",
                       legend1 = "",
                       legend2 = "",
                       legend3 = "",
                       xlab = "Estimates",
                       ylab = "Density",
                       mar = c(4, 4, 1, 1) + .2,
                       cex = 1,
                       cex.lab = 1,
                       freq = FALSE,
                       ylim = ylim,
                       length.out = 100,
                       vline = "",
                       breaks = seq(floor(min(data1, data2, data3, na.rm = TRUE)),
                                    ceiling(max(data1, data2, data3, na.rm = TRUE)),
                                    length.out = length.out
                       )) {
  par(mar = mar, cex = cex, cex.lab = cex.lab)

  dens1 <- hist(data1, breaks = breaks, plot = FALSE)$density
  dens2 <- hist(data2, breaks = breaks, plot = FALSE)$density
  dens3 <- hist(data2, breaks = breaks, plot = FALSE)$density

  ymax <- max(dens1, dens2, dens3) * 1.3
  ylim <- c(0, ymax)

  # first hist
  hist(data1,
    breaks = breaks,
    col = cols[1],
    border = FALSE,
    freq = freq,
    xlab = "Estimates",
    main = "",
    ylim = ylim
  )
  
  # first median
  lines(
    x = c(median(data1, na.rm = TRUE), median(data1, na.rm = TRUE)),
    y = c(0, ylim[2]),
    col = cols0[1]
  )
  
  # second hist
  hist(data2,
    breaks = breaks,
    col = cols[2],
    border = FALSE,
    freq = freq,
    ylim = ylim,
    add = TRUE
  )
  
  # second median
  lines(
    x = c(median(data2, na.rm = TRUE), median(data2, na.rm = TRUE)),
    y = c(0, ylim[2]),
    col = cols0[2]
  )
  
  # third hist
  hist(data3,
    breaks = breaks,
    col = cols[3],
    border = FALSE,
    freq = freq,
    ylim = ylim,
    add = TRUE
  )
  
  # third median
  lines(
    x = c(median(data3, na.rm = TRUE), median(data3, na.rm = TRUE)),
    y = c(0, ylim[2]),
    col = cols0[3]
  )
  
  box()
  lines(
    x = c(vline, vline),
    y = c(0, ylim[2])
  )
  legend(legend.loc,
    legend = c(legend1, legend2, legend3),
    fill = cols[1:3],
    border = FALSE,
    bty = "n",
    cex = 4 / 5
  )
}



## MAPS ------------------------------------------------------------------------


