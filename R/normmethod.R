# PMID: 16762068
# None
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


# # generalize log
#
# LogTrans <- function(data){
#         x <- log(data)
#         return(x)
# }

# # Power scaling
#
# PowerScaling <- function(data){
#         x <- sqrt(data)
#         return(x)
# }

# normalize to zero mean and unit variance

AutoScaling <- function(data,lv,log=T){
        if(log) data <- log(data)
        r <- rownames(data)
        c <- colnames(data)
        data2 <- t(apply(data, 1, function(x) (x - mean(x))/sd(x, na.rm=T)))
        rownames(data2) <- r
        colnames(data2) <- c
        li <- limmafit(data2,lv)
        return(li)
}

# normalize to zero mean and squared root variance

ParetoScaling <- function(data,lv,log=T){
        if(log) data <- log(data)
        r <- rownames(data)
        c <- colnames(data)
        data2 <- t(apply(data, 1, function(x) (x - mean(x))/sqrt(sd(x, na.rm=T))))
        rownames(data2) <- r
        colnames(data2) <- c
        li <- limmafit(data2,lv)
        return(li)
}


# normalize to zero mean but variance/SE

RangeScaling <- function(data,lv,log=T){
        if(log) data <- log(data)
        r <- rownames(data)
        c <- colnames(data)
        data2 <- t(apply(data, 1, function(x) (x - mean(x))/(max(x)-min(x))))
        rownames(data2) <- r
        colnames(data2) <- c
        li <- limmafit(data2,lv)
        return(li)
}

# vast scaling

VastScaling <- function(data,lv,log=T){
        if(log) data <- log(data)
        r <- rownames(data)
        c <- colnames(data)
        data2 <- t(apply(data, 1, function(x) (x - mean(x))/sd(x, na.rm=T) * mean(x)/sd(x, na.rm=T)))
        rownames(data2) <- r
        colnames(data2) <- c
        li <- limmafit(data2,lv)
        return(li)
}

# level scaling

LevelScaling <- function(data,lv,log=T){
        if(log) data <- log(data)
        r <- rownames(data)
        c <- colnames(data)
        data2 <- t(apply(data, 1, function(x) (x - mean(x))/mean(x)))
        rownames(data2) <- r
        colnames(data2) <- c
        li <- limmafit(data2,lv)
        return(li)
}

# total sum row

TotalSum <- function(data,lv,log=T){
        if(log) data <- log(data)
        r <- rownames(data)
        c <- colnames(data)
        data2 <- apply(data, 2, function(x) x/sum(x, na.rm=T))
        rownames(data2) <- r
        colnames(data2) <- c
        li <- limmafit(data2,lv)
        return(li)
}

# Median row

MedianNorm <- function(data,lv,log=T){
        if(log) data <- log(data)
        r <- rownames(data)
        c <- colnames(data)
        data2 <- apply(data, 2, function(x) x/median(x, na.rm=T))
        rownames(data2) <- r
        colnames(data2) <- c
        li <- limmafit(data2,lv)
        return(li)
}

# Mean row

MeanNorm <- function(data,lv,log=T){
        if(log) data <- log(data)
        r <- rownames(data)
        c <- colnames(data)
        data2 <- apply(data, 2, function(x) x/mean(x, na.rm=T))
        rownames(data2) <- r
        colnames(data2) <- c
        li <- limmafit(data2,lv)
        return(li)
}

# PQN

PQNorm <- function(data,lv,log=T){
        if(log) data <- log(data)
        r <- rownames(data)
        c <- colnames(data)
        ref <- apply(data[,lv == lv[1]],1,mean)
        data2 <- apply(data, 2, function(x) x/median(as.numeric(x/ref), na.rm=T))
        rownames(data2) <- r
        colnames(data2) <- c
        li <- limmafit(data2,lv)
        return(li)
}

# VSN

VSNNorm <- function(data,lv,log=T){
        if(log) data <- log(data)
        fit <- vsn::vsnMatrix(data)
        data2 <- fit@hx
        li <- limmafit(data2,lv)
        return(li)
}

# Quantile

QuanNorm <- function(data,lv,log=T){
        if(log) data <- log(data)
        data2 <- preprocessCore::normalize.quantiles(data, copy=FALSE)
        li <- limmafit(data2,lv)
        return(li)
}

# lumi rsn

LumiRobustSpline <- function(data,lv,log=T){
        if(log) data <- log(data)
        data2 <- lumi::rsn(data)
        li <- limmafit(data2,lv)
        return(li)
}

# Limma CyclicLoess
LimmaCyclicLoess <- function(data,lv,log=T){
        if(log) data <- log(data)
        data2 <- limma::normalizeCyclicLoess(data)
        li <- limmafit(data2,lv)
        return(li)
}
# AFFA CUBICSpline
LimmaCubicSpline <- function(data,lv,log=T){
        if(log) data <- log(data)
        data2 <- affy::normalize.qspline(data)
        li <- limmafit(data2,lv)
        return(li)
}
# SVA
svacor <- function(data,lv,log=T) {
        if(log) data <- log(data)
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
# iSVA
isvacor <- function(data, lv,log=T) {
        if(log) data <- log(data)
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
# PCR
pcacor <- function(data, lv,log=T) {
        if(log) data <- log(data)
        batch <- svd(data - rowMeans(data))$v[,1]
        message("Data is correcting ...")
        li <- limmafit(data,lv,batch)
        message("Done!")
        return(li)
}

# FMAT
famtcor <- function(data, lv,log=T){
        if(log) data <- log(data)
covFAMT  <- data.frame(id    = colnames(data), trmt  = as.factor(lv))
dataFAMT <- FAMT::as.FAMTdata(expression = data, covariates = covFAMT, idcovar = 1)
fitFAMT  <- FAMT::modelFAMT(dataFAMT, x = 2, test = 2)
data <- data
datacor <- fitFAMT$adjdata$expression
pValues <- fitFAMT$pval
qValues <- stats::p.adjust(pValues, method = "BH")
# get the results as list
li <- list(data, datacor, pValues, qValues)
names(li) <- c("data","dataCorrected","p-values", "q-values")
return(li)
}
# RRmix
rrmixcor <- function(data, lv,log=T){
        if(log) data <- log(data)
        n <- RRmix::nfactors(data,plot = F)
        re <- RRmix::runRRmix(t(data),lv,q.in = ifelse(n>10,1,n))
        posterior <- re[['b_g']]
        # get the results as list
        li <- list(data, posterior)
        names(li) <- c("data","posterior")
        return(li)
}


