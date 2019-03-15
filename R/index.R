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
