#' Compute pooled QC linear index according to run order
#' @param data peaks intensity list with row as peaks and column as samples
#' @param order run order of pooled QC samples
#' @param n samples numbers used for linear regression
#' @return vector for the peaks proportion with significant changes in linear regression after FDR control.
#' @export

getpqsi <- function(data, order, n=5){
        data <- data[,order(order)]
        porp <- 1:ncol(data)
        for(i in n:ncol(data)){
                p <- apply(data[,c((i-n+1):i)],1,function(x) summary(stats::lm(x~c((i-n+1):i)))$coefficients[2,4])
                # FDR control
                q <- stats::p.adjust(p,method='BH')
                porp[i] <- sum(q<0.1)/nrow(data)
        }
        return(porp[-c(1:(n-1))])
}
#' Centering batch correction
#' @param peak peaks intensity across samples
#' @param qc peaks intensity pooled QC samples
#' @param log log transformation
#' @return corrected peaks intensity
#' @export
bccenter <- function(peak, qc = NULL, log = T){
        if(log){
                peak <- log(peak+1)
                if(!is.null(qc)){
                        qc <- log(qc+1)
                        medianqc <- stats::median(qc)
                        meanpeak <- mean(peak)
                        corpeak <- peak-meanpeak+medianqc
                }else{
                        meanpeak <- mean(peak)
                        corpeak <- peak-meanpeak
                        }
        }else{
                if(!is.null(qc)){
                        medianqc <- stats::median(qc)
                        meanpeak <- mean(peak)
                        corpeak <- peak-meanpeak+medianqc
                }else{
                        meanpeak <- mean(peak)
                        corpeak <- peak-meanpeak
                }
        }
        return(corpeak)
}
#' Scaling batch correction
#' @param peak peaks intensity across samples
#' @param qc peaks intensity pooled QC samples
#' @param log log transformation
#' @return corrected peaks intensity
#' @export
bcscaling <- function(peak, qc = NULL, log = T){
        if(log){
                peak <- log(peak+1)
                if(!is.null(qc)){
                        qc <- log(qc+1)
                        meanqc <- mean(qc)
                        meanpeak <- mean(peak)
                        sdpeak <- stats::sd(peak)
                        sdqc <- stats::sd(qc)
                        corpeak <- (peak-meanpeak)/sdpeak * sdqc+meanqc
                }else{
                        meanqc <- mean(qc)
                        sdpeak <- stats::sd(peak)
                        corpeak <- (peak-meanpeak)/sdpeak
                }
        }else{
                if(!is.null(qc)){
                        meanqc <- mean(qc)
                        meanpeak <- mean(peak)
                        sdpeak <- stats::sd(peak)
                        sdqc <- stats::sd(qc)
                        corpeak <- (peak-meanpeak)/sdpeak * sdqc+meanqc
                }else{
                        meanqc <- mean(qc)
                        sdpeak <- stats::sd(peak)
                        corpeak <- (peak-meanpeak)/sdpeak
                }
        }
        return(corpeak)
}
#' Pareto scaling batch correction
#' @param peak peaks intensity across samples
#' @param qc peaks intensity pooled QC samples
#' @param log log transformation
#' @return corrected peaks intensity
#' @export
bcpareto <- function(peak, qc = NULL, log = T){
        if(log){
                peak <- log(peak+1)
                if(!is.null(qc)){
                        qc <- log(qc+1)
                        meanqc <- mean(qc)
                        meanpeak <- mean(peak)
                        sdpeak <- stats::sd(peak)
                        sdqc <- stats::sd(qc)
                        corpeak <- (peak-meanpeak)/sqrt(sdpeak) * sqrt(sdqc)+meanqc
                }else{
                        meanqc <- mean(qc)
                        sdpeak <- stats::sd(peak)
                        corpeak <- (peak-meanpeak)/sqrt(sdpeak)
                }
        }else{
                if(!is.null(qc)){
                        meanqc <- mean(qc)
                        meanpeak <- mean(peak)
                        sdpeak <- stats::sd(peak)
                        sdqc <- stats::sd(qc)
                        corpeak <- (peak-meanpeak)/sqrt(sdpeak) * sqrt(sdqc)+meanqc
                }else{
                        meanqc <- mean(qc)
                        sdpeak <- stats::sd(peak)
                        corpeak <- (peak-meanpeak)/sqrt(sdpeak)
                }
        }
        return(corpeak)
}
#' Range scaling batch correction
#' @param peak peaks intensity across samples
#' @param qc peaks intensity pooled QC samples
#' @param log log transformation
#' @return corrected peaks intensity
#' @export
bcrange <- function(peak, qc = NULL, log = T){
        if(log){
                peak <- log(peak+1)
                if(!is.null(qc)){
                        qc <- log(qc+1)
                        meanqc <- mean(qc)
                        meanpeak <- mean(peak)
                        corpeak <- (peak-meanpeak)/(max(peak) - min(peak)) * (max(qc) - min(qc))+meanqc
                }else{
                        meanqc <- mean(qc)
                        corpeak <- (peak-meanpeak)/(max(peak) - min(peak))
                }
        }else{
                if(!is.null(qc)){
                        meanqc <- mean(qc)
                        meanpeak <- mean(peak)
                        corpeak <- (peak-meanpeak)/(max(peak) - min(peak)) * (max(qc) - min(qc))+meanqc
                }else{
                        meanqc <- mean(qc)
                        corpeak <- (peak-meanpeak)/(max(peak) - min(peak))
                }
        }
        return(corpeak)
}

#' Vast scaling batch correction
#' @param peak peaks intensity across samples
#' @param qc peaks intensity pooled QC samples
#' @param log log transformation
#' @return corrected peaks intensity
#' @export
bcvast <- function(peak, qc = NULL, log = T){
        if(log){
                peak <- log(peak+1)
                if(!is.null(qc)){
                        qc <- log(qc+1)
                        meanqc <- mean(qc)
                        meanpeak <- mean(peak)
                        sdpeak <- stats::sd(peak)
                        sdqc <- stats::sd(qc)
                        corpeak <- (peak-meanpeak)/stats::sd(peak) * meanqc/sdqc
                }else{
                        meanqc <- mean(qc)
                        sdpeak <- stats::sd(peak)
                        corpeak <- (peak-meanpeak)/stats::sd(peak) * meanpeak/sdpeak
                }
        }else{
                if(!is.null(qc)){
                        meanqc <- mean(qc)
                        meanpeak <- mean(peak)
                        sdpeak <- stats::sd(peak)
                        sdqc <- stats::sd(qc)
                        corpeak <- (peak-meanpeak)/stats::sd(peak) * meanqc/sdqc
                }else{
                        meanqc <- mean(qc)
                        sdpeak <- stats::sd(peak)
                        corpeak <- (peak-meanpeak)/stats::sd(peak) * meanpeak/sdpeak
                }
        }
        return(corpeak)
}

#' Level scaling batch correction
#' @param peak peaks intensity across samples
#' @param qc peaks intensity pooled QC samples
#' @param log log transformation
#' @return corrected peaks intensity
#' @export
bclevel <- function(peak, qc = NULL, log = T){
        if(log){
                peak <- log(peak+1)
                if(!is.null(qc)){
                        qc <- log(qc+1)
                        meanqc <- mean(qc)
                        meanpeak <- mean(peak)
                        corpeak <- (peak-meanpeak)/meanpeak * meanqc+meanqc
                }else{
                        meanqc <- mean(qc)
                        corpeak <- (peak-meanpeak)/meanpeak
                }
        }else{
                if(!is.null(qc)){
                        meanqc <- mean(qc)
                        meanpeak <- mean(peak)
                        corpeak <- (peak-meanpeak)/meanpeak * meanqc+meanqc
                }else{
                        meanqc <- mean(qc)
                        corpeak <- (peak-meanpeak)/meanpeak
                }
        }
        return(corpeak)
}
