#' Centering batch correction
#' @param peak peaks intensity across samples
#' @param qc peaks intensity pooled QC samples
#' @param log log transformation
#' @return corrected peaks intensity
#' @export
bccenter <- function(peak, qc = NULL, log = TRUE){
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
bcscaling <- function(peak, qc = NULL, log = TRUE){
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
bcpareto <- function(peak, qc = NULL, log = TRUE){
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
bcrange <- function(peak, qc = NULL, log = TRUE){
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
bcvast <- function(peak, qc = NULL, log = TRUE){
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
bclevel <- function(peak, qc = NULL, log = TRUE){
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
