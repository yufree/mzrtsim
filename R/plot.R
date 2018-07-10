#' Plot ROC curve for simulation data
#' @param sim simulation results from `mzrtsim` or `simmzrt`
#' @param pvalue p value results from `limmafit` function or related batch correction method
#' @param points numbers of points for ROC curve
#' @export
limmaroc <- function(sim,pvalue,points = 100){
        indexc <- sim$conp
        TPR <- FPR <- FNR <- c(0:points)
        for(i in 0:points){
                indexi <- which(pvalue < i/points,arr.ind = T)
                TPR[i+1] <- length(intersect(indexi,indexc))/length(indexc)
                FPR[i+1] <- (length(indexi)-length(intersect(indexi,indexc)))/(nrow(sim$data) - length(indexc))
                FNR[i+1] <- 1-length(intersect(indexi,indexc))/length(indexc)
        }
        graphics::plot(TPR~FPR,col = 'red',pch = 19,type = 'l',xlim = c(0,1),ylim =c(0,1))
        graphics::lines(FNR~FPR,col = 'blue',pch = 19,type = 'l')

        height = (TPR[-1]+TPR[-length(TPR)])/2
        width = diff(FPR)
        auc <- sum(height*width)
        return(auc)
}

#' Principal component analysis(PCA) for limma result with batch
#' @param list results from `limmafit` function or related batch correction method
#' @param center parameters for PCA
#' @param scale parameters for scale
#' @param lv factor vector for the group infomation
#' @return plot
#' @examples
#' \dontrun{
#' library(faahKO)
#' library(enviGCMS)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' list <- getmr(cdfpath, pmethod = ' ')
#' li <- svacor(list$data,list$group$class)
#' limmapca(li)
#' }
#' @seealso \code{\link{isvacor}},\code{\link{svacor}}, \code{\link{pcacor}},\code{\link{limmafit}}
#' @export
limmapca <- function(list, center = T, scale = T, lv = NULL) {
        data <- list$data
        Signal <- list$signal
        Batch <- list$batch
        error <- list$error
        datacor <- list$dataCorrected
        if (is.null(lv)) {
                pch = colnames(data)
        } else {
                pch = as.character(lv)
        }

        graphics::par(mfrow = c(2, 5), mar = c(4, 4, 2.6, 1))

        pcao <- stats::prcomp(t(data), center = center, scale = scale)
        pcaoVars = signif(((pcao$sdev)^2)/(sum((pcao$sdev)^2)),
                          3) * 100
        graphics::plot(pcao, type = "l", main = "PCA")

        pca <- stats::prcomp(t(Signal), center = TRUE, scale = TRUE)
        pcaVars = signif(((pca$sdev)^2)/(sum((pca$sdev)^2)),
                         3) * 100
        graphics::plot(pca, type = "l", main = "PCA-signal")

        pcab <- stats::prcomp(t(Batch), center = center, scale = scale)
        pcabVars = signif(((pcab$sdev)^2)/(sum((pcab$sdev)^2)),
                          3) * 100
        graphics::plot(pcab, type = "l", main = "PCA-batch")

        pcae <- stats::prcomp(t(error), center = center, scale = scale)
        pcaeVars = signif(((pcae$sdev)^2)/(sum((pcae$sdev)^2)),
                          3) * 100
        graphics::plot(pcae, type = "l", main = "PCA-error")

        pcac <- stats::prcomp(t(datacor), center = center,
                              scale = scale)
        pcacVars = signif(((pcac$sdev)^2)/(sum((pcac$sdev)^2)),
                          3) * 100
        graphics::plot(pcac, type = "l", main = "PCA-corrected")

        graphics::plot(pcao$x[, 1], pcao$x[, 2], xlab = paste("PC1:",
                                                              pcaoVars[1], "% of Variance Explained"), ylab = paste("PC2:",
                                                                                                                    pcaoVars[2], "% of Variance Explained"), pch = pch,
                       cex = 2, main = "PCA")

        graphics::plot(pca$x[, 1], pca$x[, 2], xlab = paste("PC1:",
                                                            pcaVars[1], "% of Variance Explained"), ylab = paste("PC2:",
                                                                                                                 pcaVars[2], "% of Variance Explained"), pch = pch,
                       cex = 2, main = "PCA-signal")

        graphics::plot(pcab$x[, 1], pcab$x[, 2], xlab = paste("PC1:",
                                                              pcabVars[1], "% of Variance Explained"), ylab = paste("PC2:",
                                                                                                                    pcabVars[2], "% of Variance Explained"), pch = pch,
                       cex = 2, main = "PCA-batch")

        graphics::plot(pcae$x[, 1], pcae$x[, 2], xlab = paste("PC1:",
                                                              pcaeVars[1], "% of Variance Explained"), ylab = paste("PC2:",
                                                                                                                    pcaeVars[2], "% of Variance Explained"), pch = pch,
                       cex = 2, main = "PCA-error")

        graphics::plot(pcac$x[, 1], pcac$x[, 2], xlab = paste("PC1:",
                                                              pcacVars[1], "% of Variance Explained"), ylab = paste("PC2:",
                                                                                                                    pcacVars[2], "% of Variance Explained"), pch = pch,
                       cex = 2, main = "PCA-corrected")
}

#' Filter the data with p value and q value and show them as heatmap
#' @param list results from `limmafit` function or related batch correction method
#' @param lv factor vector for the group infomation
#' @param pt threshold for p value, default is 0.05
#' @param qt threshold for q value, default is 0.05
#' @param index index for selected peaks
#' @return heatmap for the data
#' @examples
#' \dontrun{
#' sim <- mzrtsim()
#' li <- svacor(log(sim$data), as.factor(sim$con))
#' limmaplot(li,as.factor(sim$con))
#' }
#' @seealso \code{\link{isvacor}},\code{\link{svacor}}, \code{\link{pcacor}},\code{\link{limmafit}}
#' @export
limmaplot <- function(list, lv, pt = 0.05,
                      qt = 0.05, index = NULL) {
        data <- list$data[, order(lv)]
        signal <- list$signal[, order(lv)]
        batch <- list$batch[, order(lv)]
        error <- list$error[, order(lv)]
        datacor <- list$dataCorrected[, order(lv)]
        pValues <- list$"p-values"
        qValues <- list$"q-values"
        if (!is.null(index)) {
                data <- data[index, order(lv)]
                signal <- signal[index, order(lv)]
                batch <- batch[index, order(lv)]
                error <- error[index, order(lv)]
                datacor <- datacor[index, order(lv)]
                pValues <- pValues[index]
                qValues <- qValues[index]
        }
        # line position
        pos <- cumsum(as.numeric(table(lv)/sum(table(lv)))) -
                as.numeric(table(lv)/sum(table(lv)))/2
        posv <- cumsum(as.numeric(table(lv)/sum(table(lv))))[1:(nlevels(lv) -
                                                                        1)]
        icolors <- (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11,
                                                                             "RdYlBu"))))(100)
        # plot the scale
        plotchange <- function(zlim) {
                breaks <- seq(zlim[1], zlim[2], round((zlim[2] -
                                                               zlim[1])/10))
                poly <- vector(mode = "list", length(icolors))
                graphics::plot(1, 1, t = "n", xlim = c(0, 1), ylim = zlim,
                               xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i",
                               ylab = "", xlab = "", frame.plot = F)
                graphics::axis(4, at = breaks, labels = round(breaks),
                               las = 1, pos = 0.4, cex.axis = 0.8)
                p <- graphics::par("usr")
                graphics::text(p[2]+1, mean(p[3:4]), labels = "intensity",
                               xpd = NA, srt = -90)
                bks <- seq(zlim[1], zlim[2], length.out = (length(icolors) +
                                                                   1))
                for (i in seq(poly)) {
                        graphics::polygon(c(0.1, 0.1, 0.3, 0.3), c(bks[i],
                                                                   bks[i + 1], bks[i + 1], bks[i]), col = icolors[i],
                                          border = NA)
                }
        }
        # plot without batch
        plotimage1 <- function(data, signal, error, zlim) {
                graphics::image(t(data), col = icolors, xlab = "samples",
                                main = "peaks", xaxt = "n", yaxt = "n", zlim = zlim)
                graphics::axis(1, at = pos, labels = levels(lv),
                               cex.axis = 0.8)
                graphics::axis(2, at = seq(0, 1, 1/(nrow(data) -
                                                            1)), labels = rownames(data), cex.axis = 1,
                               las = 2)
                graphics::abline(v = posv)

                graphics::image(t(signal), col = icolors, xlab = "samples",
                                main = "signal", xaxt = "n", yaxt = "n",
                                zlim = zlim)
                graphics::axis(1, at = pos, labels = levels(lv),
                               cex.axis = 0.8)
                graphics::axis(2, at = seq(0, 1, 1/(nrow(signal) -
                                                            1)), labels = rownames(signal), cex.axis = 1,
                               las = 2)
                graphics::abline(v = posv)

                graphics::image(t(error), col = icolors, xlab = "samples",
                                main = "error", xaxt = "n", yaxt = "n",
                                zlim = zlim)
                graphics::axis(1, at = pos, labels = levels(lv),
                               cex.axis = 0.8)
                graphics::axis(2, at = seq(0, 1, 1/(nrow(error) -
                                                            1)), labels = rownames(error), cex.axis = 1,
                               las = 2)
                graphics::abline(v = posv)
        }
        # plot with batch
        plotimage2 <- function(data, signal, batch, error,
                               datacor, zlim) {
                graphics::image(t(data), col = icolors, xlab = "samples",
                                main = "peaks", xaxt = "n", yaxt = "n", zlim = zlim)
                graphics::axis(1, at = pos, labels = levels(lv))
                graphics::axis(2, at = seq(0, 1, 1/(nrow(data) -
                                                            1)), labels = rownames(data), las = 1)
                graphics::abline(v = posv)

                graphics::image(t(signal), col = icolors, xlab = "samples",
                                main = "signal", xaxt = "n", yaxt = "n",
                                zlim = zlim)
                graphics::axis(1, at = pos, labels = levels(lv),
                               cex.axis = 1)
                graphics::axis(2, at = seq(0, 1, 1/(nrow(signal) -
                                                            1)), labels = rownames(signal), las = 1)
                graphics::abline(v = posv)

                graphics::image(t(batch), col = icolors, xlab = "samples",
                                main = "batch", xaxt = "n", yaxt = "n",
                                zlim = zlim)
                graphics::axis(1, at = pos, labels = levels(lv))
                graphics::axis(2, at = seq(0, 1, 1/(nrow(batch) -
                                                            1)), labels = rownames(batch), las = 1)
                graphics::abline(v = posv)

                graphics::image(t(error), col = icolors, xlab = "samples",
                                main = "error", xaxt = "n", yaxt = "n",
                                zlim = zlim)
                graphics::axis(1, at = pos, labels = levels(lv))
                graphics::axis(2, at = seq(0, 1, 1/(nrow(error) -
                                                            1)), labels = rownames(error), las = 1)
                graphics::abline(v = posv)

                graphics::image(t(datacor), col = icolors, xlab = "samples",
                                main = "peaks-corrected", xaxt = "n", yaxt = "n",
                                zlim = zlim)
                graphics::axis(1, at = pos, labels = levels(lv))
                graphics::axis(2, at = seq(0, 1, 1/(nrow(datacor) -
                                                            1)), labels = rownames(datacor), las = 1)
                graphics::abline(v = posv)

        }

        # plot heatmap
        if (is.null(batch)) {
                if (sum(pValues < pt & qValues <
                        qt) != 0) {
                        message("No batch while p-values and q-values have results")
                        graphics::layout(matrix(rep(c(1, 1, 1, 1, 2, 2,2, 3,3,
                                                      3, 4,4), 12), 12, 12, byrow = TRUE))
                        data <- data[pValues < pt & qValues < qt, ]
                        signal <- signal[pValues < pt & qValues < qt,
                                         ]
                        error <- error[pValues < pt & qValues < qt,
                                       ]

                        zlim <- range(c(data, signal, error))
                        graphics::par(mar = c(3, 4, 2, 1))
                        plotimage1(data, signal, error, zlim)
                        graphics::par(mar = c(3, 1, 2, 4))
                        plotchange(zlim)

                        li <- list(data, pValues < pt & qValues < qt)
                        names(li) <- c("data", "pqvalues")
                        return(li)
                } else {
                        message("No batch while p-values and q-values have no results")
                        graphics::layout(matrix(rep(c(1, 1, 1, 1, 2, 2,2, 3,3,
                                                      3, 4,4), 12), 12, 12,  byrow = TRUE))
                        zlim <- range(c(data, signal, error))
                        graphics::par(mar = c(3, 5, 2, 1))
                        plotimage1(data, signal, error, zlim)
                        graphics::par(mar = c(3, 1, 2, 5))
                        plotchange(zlim)
                }
        } else {
                if (sum(pValues < pt & qValues <
                        qt) != 0) {
                        message("p-values and q-values have results")
                        zlim <- range(c(data, signal, batch, error,
                                        datacor))
                        graphics::layout(matrix(rep(c(1, 1, 1, 2, 2,
                                                      3, 3, 4, 4, 5, 5, 5, 6,6), 14), 14, 14, byrow = TRUE))
                        data <- data[pValues < pt & qValues < qt,
                                     ]
                        signal <- signal[pValues < pt & qValues <
                                                 qt, ]
                        batch <- batch[pValues < pt & qValues <
                                               qt, ]
                        error <- error[pValues < pt & qValues <
                                               qt, ]
                        datacor <- datacor[pValues < pt & qValues <
                                                   qt, ]
                        zlim <- range(c(data, signal, batch, error,
                                        datacor))
                        graphics::par(mar = c(3, 4, 2, 1))
                        plotimage2(data, signal, batch, error, datacor,
                                   zlim)
                        graphics::par(mar = c(3, 1, 2, 4))
                        plotchange(zlim)
                        li <- list(datacor, data, pValues < pt &
                                           qValues < qt)
                        names(li) <- c("dataCorrected", "data", "pqvalues")
                        return(li)
                } else {
                        message("p-values and q-values have no results")
                        graphics::layout(matrix(rep(c(1, 1, 1, 2, 2,
                                                      3, 3, 4, 4, 5, 5, 5, 6,6), 14), 14, 14, byrow = TRUE))
                        zlim <- range(c(signal, data, batch, error,
                                        datacor))
                        graphics::par(mar = c(3, 4, 2, 1))
                        plotimage2(data, signal, batch, error, datacor,
                                   zlim)
                        graphics::par(mar = c(3, 1, 2, 4))
                        plotchange(zlim)
                }
        }
}

#' Plot the heatmap for the data
#' @param data matrix for heatmap
#' @param lv factor vector for the group infomation
#' @param index row index for the data matrix to be ploted
#' @export
heatplot <- function(data,lv,index = NULL){
        icolors <- (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdYlBu"))))(100)
        zlim <- range(data)
        if (!is.null(index)) {
                data <- data[index, order(lv)]
        }else{
                data <- data[, order(lv)]
        }
        plotchange <- function(zlim) {
                breaks <- seq(zlim[1], zlim[2], round((zlim[2] -
                                                               zlim[1])/10))
                poly <- vector(mode = "list", length(icolors))
                graphics::plot(1, 1, t = "n", xlim = c(0, 1), ylim = zlim,
                               xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i",
                               ylab = "", xlab = "", frame.plot = F)
                graphics::axis(4, at = breaks, labels = round(breaks),
                               las = 1, pos = 0.4, cex.axis = 0.8)
                p <- graphics::par("usr")
                graphics::text(p[2] + 2, mean(p[3:4]), labels = "intensity",
                               xpd = NA, srt = -90)
                bks <- seq(zlim[1], zlim[2], length.out = (length(icolors) +
                                                                   1))
                for (i in seq(poly)) {
                        graphics::polygon(c(0.1, 0.1, 0.3, 0.3), c(bks[i],
                                                                   bks[i + 1], bks[i + 1], bks[i]), col = icolors[i],
                                          border = NA)
                }
        }
        pos <- cumsum(as.numeric(table(lv)/sum(table(lv)))) -
                as.numeric(table(lv)/sum(table(lv)))/2
        posv <- cumsum(as.numeric(table(lv)/sum(table(lv))))[1:(nlevels(lv) -
                                                                        1)]


        graphics::layout(matrix(rep(c(1, 1, 1, 2), 10), 10, 4, byrow = TRUE))
        graphics::par(mar = c(3, 6, 2, 1))
        graphics::image(t(data), col = icolors, xlab = "samples",
                        main = "peaks intensities on log scale", xaxt = "n", yaxt = "n", zlim = zlim)
        graphics::axis(1, at = pos, labels = levels(lv),
                       cex.axis = 0.8)
        graphics::axis(2, at = seq(0, 1, 1/(nrow(data) -
                                                    1)), labels = rownames(data), cex.axis = 1, las = 2)
        graphics::abline(v = posv)
        graphics::par(mar = c(3, 1, 2, 6))
        plotchange(zlim)
}

#' Relative Log Abundance (RLA) plots
#' @param data data as mzrt profile
#' @param lv factor vector for the group infomation
#' @param type 'g' means group median based, other means all samples median based.
#' @return Relative Log Abundance (RLA) plots
#' @examples
#' sim <- mzrtsim()
#' rlaplot(sim$data, as.factor(sim$con))
#' @export
rlaplot <- function(data, lv, type = "g") {
        data <- log(data)
        data[is.nan(data)] <- 0
        outmat = NULL

        if (type == "g") {
                for (lvi in levels(lv)) {
                        submat <- data[, lv == lvi]
                        median <- apply(submat, 1, median)
                        tempmat <- sweep(submat, 1, median, "-")
                        outmat <- cbind(outmat, tempmat)
                }
        } else {
                median <- apply(data, 1, median)
                outmat <- sweep(data, 1, median, "-")

        }

        outmat <- outmat[, order(lv)]
        graphics::boxplot(outmat, col = lv[order(lv)])
        graphics::abline(h = 0)
}
#' Relative Log Abundance Ridge(RLAR) plots
#' @param data data as mzrt profile
#' @param lv factor vector for the group infomation
#' @param type 'g' means group median based, other means all samples median based.
#' @return Relative Log Abundance (RLA) plots
#' @examples
#' sim <- mzrtsim()
#' ridgesplot(sim$data, as.factor(sim$con))
#' @export
ridgesplot <- function(data, lv, type = "g") {
        data <- log(data)
        data[is.nan(data)] <- 0
        outmat = NULL

        if (type == "g") {
                for (lvi in levels(lv)) {
                        submat <- data[, lv == lvi]
                        median <- apply(submat, 1, median)
                        tempmat <- sweep(submat, 1, median, "-")
                        outmat <- cbind(outmat, tempmat)
                }
        } else {
                median <- apply(data, 1, median)
                outmat <- sweep(data, 1, median, "-")

        }

        ov <- reshape2::melt(outmat)
        colnames(outmat) <- lv[order(lv)]
        ov2 <- reshape2::melt(outmat)
        ov2$group <- as.factor(ov2$Var2)
        ggplot2::ggplot(ov, ggplot2::aes(x = ov$value, y = ov$Var2,
                                         fill = ov2$group)) + ggridges::geom_density_ridges(stat = "binline",
                                                                                            bins = 100) + ggplot2::xlim(-0.5, 0.5) + ggplot2::scale_fill_discrete(name = "Group") +
                ggplot2::labs(x = "Relative Log Abundance", y = "Samples")
}



