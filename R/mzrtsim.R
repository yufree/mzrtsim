#' A list containing HMDB GC-EI-MS spectra database
#' @docType data
#' @usage data(hmdbcms)
#' @format A list with compounds from hmdb for GC-MS simulation
"hmdbcms"

#' A list containing MoNA LC-MS spectra database
#' @docType data
#' @usage data(monams1)
#' @format A list with compounds from MoNA LC-MS spectra database for GC-MS simulation
"monams1"

#' Generate simulated count data with batch effects for npeaks
#'
#' @param ncomp compound numbers to be selected, default 100
#' @param fc fold change of compounds with the same length of compounds numbers, default NULL
#' @param ncond Number of conditions to simulate
#' @param ncpeaks percentage of compounds influenced by conditions
#' @param nbatch Number of batches to simulate
#' @param nbpeaks percentage of peaks influenced by batchs
#' @param npercond Number of samples per condition to simulate
#' @param nperbatch Number of samples per batch to simulate
#' @param batchtype type of batch. 'm' means monotonic, 'b' means block, 'r' means random error, 'mb' means mixed mode of monotonic and block, default 'mb'
#' @param rtsim logical, simulate retention time or not
#' @param db compound database with MS1 data. e.g. hmdbcms or monams1
#' @param seed Random seed for reproducibility
#' @details the numbers of batch columns should be the same with the condition columns.
#' @return list with rtmz data matrix, row index of peaks influenced by conditions, row index of peaks influenced by batchs, column index of conditions, column of batchs, raw condition matrix, raw batch matrix, peak mean across the samples, peak rsd across the samples
#' @seealso \code{\link{simdata}}
#' @export
#' @examples
#' data(hmdbcms)
#' sim <- mzrtsim(db=hmdbcms)
mzrtsim <- function(ncomp = 100,
                    fc=NULL,
                    ncond = 2,
                    ncpeaks = 0.1,
                    nbatch = 3,
                    nbpeaks = 0.1,
                    npercond = 10,
                    nperbatch = c(8, 5, 7),
                    batchtype = 'mb',
                    rtsim = TRUE,
                    db = NULL,seed = 42){
    set.seed(seed)
    if(is.null(fc)){
        fc <- stats::runif(ncomp,0,100000)
    }
    batch <- rep(1:nbatch, nperbatch)
    condition <- rep(1:ncond, npercond)
    # check the col numbers
    if (length(batch) != length(condition)) {
        stop("Try to use the same numbers for both batch and condition columns")
    }
    ncol <- length(batch)
    # change the column order
    batch <- batch[sample(ncol)]
    condition <- condition[sample(ncol)]
    # generate the col name
    bc <- paste0("C", condition, "B", batch, "_", c(1:ncol))

    # generate the base peaks
    if(is.null(db)){
        stop("You need database to simulate.")
    }else{
        name <- sapply(db,function(x) x$name)
        nameli <- sample(unique(name),ncomp)
        z <- db[which(name %in% nameli)]
        zname <- sapply(z,function(x) x$name)
        idx <- sapply(unique(zname), function(x) sample(which(zname==x),1))
        z <- z[idx]
        name <- sapply(z, function(x) rep(x$name,length(x$spectra$mz)))
        mz <- sapply(z, function(x) round(x$spectra$mz,4))
        ins <- sapply(z, function(x) x$spectra$ins)
        nins <- Map("*", ins, fc)
        namelist <- unlist(sapply(z, function(x) x$name))
        namelenth <- sapply(z, function(x) length(x$spectra$mz))
        compname <- unlist(name)
        compmz <- unlist(mz)
        compins <- unlist(nins)
    }
    # get peak numbers
    npeaks <- length(compname)
    # generate the rsd for base peaks
    samplersd <- stats::rlnorm(npeaks,meanlog = 2, sdlog = 0.5)
    # get the matrix
    matrix0 <- matrix <- matrix(0, nrow = npeaks, ncol = ncol)
    colnames(matrix0) <- colnames(matrix) <- bc
    rownames(matrix0) <- rownames(matrix) <- compname

    for (i in 1:npeaks) {
        samplei <- abs(stats::rnorm(ncol, mean = compins[i],
                                    sd = compins[i] * samplersd[i]/100))
        matrix[i, ] <- samplei
    }

    # simulate retention time
    if(rtsim){
        rtraw <- stats::runif(ncomp,min=0,max=1200)
        rt <- rep(rtraw,namelenth)
    }

    # get the numbers of signal and batch peaks
    ncompeak <- ncomp * ncpeaks
    nbpeak <- npeaks * nbpeaks
    # simulation of condition
    index <- sample(1:ncomp, ncompeak)
    compcon <- namelist[index]
    matrixc <- matrix[compname %in% compcon, ]
    ncpeak <- nrow(matrixc)
    changec <- NULL
    for (i in unique(condition)) {
        colindex <- condition == i
        change <- exp(stats::rnorm(ncpeak))
        matrixc[, colindex] <- matrixc[, colindex] * change
        changec <- cbind(changec, change)
    }
    matrix[compname %in% compcon, ] <- matrixc
    # simulation of batch
    indexb <- sample(1:npeaks, nbpeak)
    indexb <- 1:npeaks %in% indexb
    matrixb <- matrix[indexb, ]
    matrixb0 <- matrix0[indexb, ]
    changeb <- changem <- changer <-  NULL
    # check batch mode
    if(!grepl('m|b|r',batchtype)){
        stop("Batch type should be 'm', 'r' and/or 'b'")
    }
    # generate random batch effect
    if(grepl('r',batchtype)){
        for (i in 1:nrow(matrixb)){
            change <- abs(stats::rnorm(ncol(matrixb)))
            matrixb[i,] <- matrixb[i,]*change
            matrixb0[i,] <- matrixb0[i,]*change
            changer <- change
        }
    }
    # generate increasing/decreasing batch effect
    if(grepl('m',batchtype)){
        for (i in 1:nrow(matrixb)){
            changet <- seq(1,ncol(matrixb),length.out = ncol(matrixb)) * exp(stats::rnorm(1))
            change <- if (sample(c(T,F),1)) changet else rev(changet)
            matrixb[i,] <- matrixb[i,]*change
            matrixb0[i,] <- matrixb0[i,]*change
            changem <- rbind(changem,change)
        }
    }
    # generate block batch effect
    if(grepl('b',batchtype)){
        for (i in 1:nbatch) {
            colindex <- batch == i
            change <- exp(stats::rnorm(nbpeak))
            matrixb[, colindex] <- matrixb[, colindex] * change
            matrixb0[, colindex] <- matrixb0[, colindex] *
                change
            changeb <- cbind(changeb, change)
        }
    }
    matrix[indexb, ] <- matrixb
    # get peaks mean across the samples
    means <- apply(matrix, 1, mean)
    # get peaks rsd across the samples
    sds <- apply(matrix, 1, stats::sd)
    rsds <- sds/means
    # get the group
    group <- cbind.data.frame(sample_name=bc,sample_group = condition)
    return(list(data = data.frame(matrix), mz = as.numeric(compmz), rt = rt,group = group, name = namelist, conp = compname %in% compcon, batchp = indexb, batch = batch, cmatrix = matrixc, changec = changec, bmatrix = matrixb0, changeb = changeb, changer = changer, changem = changem,
                matrix = matrix0, mean = means, rsd = rsds))

}
#' Save the simulated data as csv files
#' @param sim list from `mzrtsim`
#' @param name file name
#' @seealso \code{\link{mzrtsim}}
#' @export
simdata <- function(sim, name = "sim") {
    # for metaboanalyst
    data <- rbind.data.frame(sim$con, sim$data)
    filename1 <- paste0(name, "raw.csv")
    utils::write.csv(data, file = filename1)
    # for all
    filename2 <- paste0(name, "raw2.csv")
    data2 <- rbind.data.frame(Label = sim$con, sim$data)
    colnames(data2) <- c(1:ncol(data2))

    utils::write.csv(t(data2), file = filename2)
    # for biological influnced peaks
    filename3 <- paste0(name, "con.csv")
    data3 <- rbind.data.frame(sim$con, sim$batch, sim$cmatrix)
    utils::write.csv(data3, file = filename3)
    # for batch influnced peaks
    filename4 <- paste0(name, "bat.csv")
    data4 <- rbind.data.frame(sim$con, sim$batch, sim$bmatrix)
    utils::write.csv(data4, file = filename4)
    # for batch matrix
    filename5 <- paste0(name, "batchmatrix.csv")
    data5 <- rbind.data.frame(sim$con,sim$batch, sim$matrix)
    utils::write.csv(data5, file = filename5)
    # for compounds
    filename6 <- paste0(name, "comp.csv")
    data6 <- rbind.data.frame(sim$con, sim$batch, sim$compmatrix)
    utils::write.csv(data6, file = filename6)
    # for compounds folds change
    filename7 <- paste0(name, "compchange.csv")
    utils::write.csv(sim$changec, file = filename7)
    # for batch folds change
    filename8 <- paste0(name, "blockbatchange.csv")
    utils::write.csv(sim$changeb, file = filename8)
    # for monotonic change
    filename9 <- paste0(name, "monobatchange.csv")
    utils::write.csv(sim$changem, file = filename9)
}
