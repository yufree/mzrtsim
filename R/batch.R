#' LIMMA to test batch effects by p-value and q-value
#' @param data data as mzrt profile
#' @param lv vector for the group information
#' @param batch vector for the batch information
#' @return list object with raw data, corrected data, signal part, batch part, random errors part,  p-values, q-values. If no batch index, corresponding part would miss.
#' @seealso \code{\link{isvacor}},\code{\link{svacor}}, \code{\link{pcacor}},\code{\link{famtcor}}
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' list <- getmr(cdfpath, pmethod = ' ')
#' li <- limmafit(list$data,list$group$class)
#' }
#' @export
limmafit <- function(data, lv, batch = NULL){
        mod <- stats::model.matrix(~lv)
        mod0 <- as.matrix(c(rep(1, ncol(data))))
        datacor <- signal <- error <- pValues <-  qValues <- NULL
        if(is.null(batch)){
                batch <- NULL
                # limma fit
                lmfit <- limma::lmFit(data, mod)
                signal <- lmfit$coef[, 1:nlevels(lv)] %*% t(mod[, 1:nlevels(lv)])
                error <- data - signal
                rownames(signal) <- rownames(error) <- rownames(data)
                colnames(signal) <- colnames(error) <- colnames(data)
                # find the peaks with significant differences by F test
                # with BH correction for fdr control without correction
                pValues = sva::f.pvalue(data, mod, mod0)
                qValues = stats::p.adjust(pValues, method = "BH")
        }else{
                modcor <- cbind(mod,batch)
                modcor0 <- cbind(mod0,batch)
                lmfit <- limma::lmFit(data, modcor)
                # data decomposition with batch
                batch <- lmfit$coef[, (nlevels(lv) + 1):(nlevels(lv) + NCOL(batch))] %*% t(modcor[, (nlevels(lv) + 1):(nlevels(lv) + NCOL(batch))])
                signal <- lmfit$coef[, 1:nlevels(lv)] %*% t(modcor[, 1:nlevels(lv)])
                error <- data - signal - batch
                datacor <- signal + error
                rownames(datacor) <- rownames(batch) <- rownames(signal) <- rownames(error) <- rownames(data)
                colnames(datacor) <- colnames(batch) <- colnames(signal) <- colnames(error) <- colnames(data)

                # find the peaks with significant differences by F test
                # with BH correction for fdr control
                pValues = sva::f.pvalue(data, modcor, modcor0)
                qValues = stats::p.adjust(pValues, method = "BH")
        }
        # get the results as list
        li <- list(data, datacor, signal, batch, error, pValues, qValues)
        names(li) <- c("data","dataCorrected","signal","batch", "error", "p-values", "q-values")
        return(li)
}

#' Use Surrogate Variable Analysis(SVA) to correct the unknown batch effects
#' @param data data as mzrt profile
#' @param lv factor vector for the group infomation
#' @details Use Surrogate Variable Analysis(SVA) to correct the unknown batch effects
#' @return list object with raw data, corrected data, signal part, batch part, random errors part,  p-values, q-values. If no surrogate variables were found, corresponding part would miss.
#' @seealso \code{\link{isvacor}}, \code{\link{pcacor}},\code{\link{limmafit}},\code{\link{famtcor}}
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' list <- getmr(cdfpath, pmethod = ' ')
#' li <- svacor(list$data,list$group$class)
#' }
#' @export
svacor <- function(data, lv) {
    mod <- stats::model.matrix(~lv)
    svafit <- sva::sva(data, mod)
    if (svafit$n.sv == 0) {
        message("No surrogate variable found")
        li <- limmafit(data,lv)
    } else {
        message("Data is correcting ...")
        batch <- svafit$sv
        li <- limmafit(data,lv,batch)
        message("Done!")
    }
    return(li)
}

#' Use Independent Surrogate Variable Analysis(ISVA) to correct the unknown batch effects
#' @param data data as mzrt profile
#' @param lv factor vector for the group infomation
#' @details Use Independent Surrogate Variable Analysis(ISVA) to correct the unknown batch effects
#' @return list object with raw data, corrected data, signal part, batch part, random errors part,  p-values, q-values. If no independent surrogate variable found, corresponding part would miss.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' list <- getmr(cdfpath, pmethod = ' ')
#' list <- isvacor(list$data,list$group$class)
#' }
#' @seealso \code{\link{svacor}}, \code{\link{pcacor}},\code{\link{limmafit}},\code{\link{famtcor}}
#' @export
isvacor <- function(data, lv) {
    isvafit <- isva::DoISVA(data, lv, factor.log = T)
    if (isvafit$nsv == 0) {
        message("No surrogate variable found")
        li <- limmafit(data,lv)
    } else {
        message("Data is correcting ...")
        batch <- isvafit$isv
        li <- limmafit(data,lv,batch)
        message("Done!")
    }
    return(li)
}

#' Use Principal component analysis(PCA) to correct the unknown batch effects
#' @param data data as mzrt profile
#' @param lv factor vector for the group infomation
#' @details Use Principal component analysis(PCA) to correct the unknown batch effects
#' @return list object with various components such raw data, corrected data, signal part, random errors part, batch part, p-values, q-values.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' list <- getmr(cdfpath, pmethod = ' ')
#' li <- pcacor(list$data,list$group$class)
#' }
#' @seealso \code{\link{isvacor}},\code{\link{svacor}}, \code{\link{limmafit}},\code{\link{famtcor}}
#' @export
pcacor <- function(data, lv) {
        batch <- svd(data - rowMeans(data))$v[,1]
        message("Data is correcting ...")
        li <- limmafit(data,lv,batch)
        message("Done!")
        return(li)
}

#' Use Factor Analysis for Multiple Testing(FAMT) to correct the unknown batch effects
#' @param data data as mzrt profile
#' @param lv factor vector for the group infomation
#' @details Use Factor Analysis for Multiple Testing(FAMT) to correct the unknown batch effects
#' @return list object with various components such raw data, corrected data, p-values, q-values.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' list <- getmr(cdfpath, pmethod = ' ')
#' li <- famtcor(list$data,list$group$class)
#' }
#' @seealso \code{\link{isvacor}},\code{\link{svacor}}, \code{\link{limmafit}},\code{\link{pcacor}}
#' @export
famtcor <- function(data, lv){                                           covFAMT  <- data.frame(id    = colnames(data),
                                  trmt  = as.factor(lv))

        dataFAMT <- FAMT::as.FAMTdata(expression = data, covariates = covFAMT, idcovar    = 1)
        fitFAMT  <- FAMT::modelFAMT(dataFAMT,                            x = 2, test = 2)
        data <- data
        datacor <- fitFAMT$adjdata$expression
        pValues <- fitFAMT$pval
        qValues <- stats::p.adjust(pValues, method = "BH")
        # get the results as list
        li <- list(data, datacor, pValues, qValues)
        names(li) <- c("data","dataCorrected","p-values", "q-values")
        return(li)
        }

#' Use Random main effect and Random compound-specific error variance with a mixture structure(RRmix) to correct the unknown batch effects
#' @param data data as mzrt profile
#' @param lv factor vector for the group infomation
#' @details Use Random main effect and Random compound-specific error variance with a mixture structure(RRmix) to correct the unknown batch effects
#' @return list object with various components such raw data, posterior probability.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' list <- getmr(cdfpath, pmethod = ' ')
#' li <- rrmixcor(list$data,list$group$class)
#' }
#' @seealso \code{\link{isvacor}},\code{\link{svacor}}, \code{\link{limmafit}},\code{\link{pcacor}}
#' @export
rrmixcor <- function(data, lv){
        n <- RRmix::nfactors(data,plot = F)
        re <- RRmix::runRRmix(t(data),lv,q.in = ifelse(n>10,1,n))
        posterior <- re[['b_g']]
# get the results as list
li <- list(data, posterior)
names(li) <- c("data","posterior")
return(li)
}
