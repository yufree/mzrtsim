#' A vector containing m/z from matrix blank
#' @docType data
#' @usage data(mzm)
#' @format A vector containing m/z from matrix blank for LC-HRMS simulation
"mzm"

#' Generate simulated mzml data and compounds list
#' @param db compound database with MS1 data. e.g. hmdbcms or monams1
#' @param name file name of mzml
#' @param n compound numbers from database, if compound is not NULL, n will be compound number, default 100
#' @param inscutoff intensity cutoff for MS1 spectra, default 0.05
#' @param mzrange m/z range for simulation, peaks out of the range will be removed, default c(100,1000)
#' @param rtrange retention time range for simulation, default c(0,600)
#' @param ppm m/z shift in ppm for input compounds, default 5
#' @param sampleppm m/z shift in ppm within one sample, default 5
#' @param mzdigit m/z digits, default 5
#' @param noisesd standard deviation for normal distribution of m/z shift, default 0.5
#' @param scanrate time for each full scan, default 0.2 second or 5 spectra per second
#' @param pwidth peak width for the compound. If it's a single value, simulated peaks' width will use this number as the lambda of Poisson distribution. If it's a numeric vector, it will be used as the peak width for each compounds.
#' @param baseline noise baseline, default 100
#' @param baselinesd standard deviation for noise, default 30
#' @param sn signal to noise ratio of each compound, default 100 for all compounds when baseline is 100
#' @param tailingfactor tailing factor for peaks, larger means larger tailing, default 1.2
#' @param compound numeric compound index in the database for targeted analysis, default NULL
#' @param rtime retention time for the compounds, if NULL, retention time will be simulated. Default NULL
#' @param tailingindex numeric index for tailing compounds, if NULL, all peaks will tailing. Default NULL
#' @param seed Random seed for reproducibility
#' @param unique if TRUE, one compound will have one spectra. Default FALSE
#' @param matrix if TRUE, m/z from experimental data will be used for background m/z simulation.Default FALSE
#' @param matrixmz custom matrix m/z vector, default NULL and predefined list will be used.
#' @return one mzML file for simulated data and one csv file the simulated compounds with retention time, m/z and name
#' @export
#' @examples
#' data(monams1)
#' simmzml(db=monams1, name = 'test')
simmzml <-
        function(db,
                 name,
                 n = 100,
                 inscutoff = 0.05,
                 mzrange = c(100,1000),
                 rtrange = c(0, 600),
                 ppm = 5,
                 sampleppm = 5,
                 mzdigit = 5,
                 noisesd = 0.5,
                 scanrate = 0.2,
                 pwidth = 10,
                 baseline = 100,
                 baselinesd = 30,
                 sn = 100,
                 tailingfactor = 1.2,
                 compound=NULL,
                 rtime=NULL,
                 tailingindex = NULL,
                 seed=42,
                 unique=FALSE,
                 matrix=FALSE,
                 matrixmz=NULL) {
                if(unique){
                        uniquecpidx <- sapply(db, function(x) x$name)
                        db <- db[!duplicated(uniquecpidx)]
                }

                if(is.null(compound)){
                        set.seed(seed)
                        sub <- db[sample(length(db), n)]
                }else{
                        sub <- db[compound]
                        n <- length(compound)
                }
                rtime0 <- seq(rtrange[1], rtrange[2], scanrate)
                if(is.null(rtime)){
                        set.seed(seed)
                        rtime <- sample(rtime0, n, replace = T)
                }
                if(length(rtime)!=n){
                        stop('Element numbers of retention time vector should have the same number of compounds.')
                }
                if(length(pwidth)==1){
                        set.seed(seed)
                        # uniform distributed peak width
                        peakrange <- stats::rpois(n, lambda = pwidth)
                }else{
                        peakrange <- pwidth
                }
                if(length(sn)==1){
                        sn <- rep(sn,100)
                }else if(length(sn)!=n){
                        set.seed(seed)
                        sn <- sample(sn,n,replace = T)
                }
                # chromotograghy simulation for the compound
                re <- c()
                for (i in c(1:n)) {
                        gaussian_peak <-
                                stats::dnorm(rtime0, mean = rtime[i], sd = peakrange[i] / 4)
                        gaussian_peak <- gaussian_peak/max(gaussian_peak)*100*sn[i]
                        # tailing simulation
                        tailing_peak <-
                                stats::dnorm(rtime0, mean = rtime[i], sd = (2*tailingfactor-1)*peakrange[i] / 4)
                        tailing_peak <- tailing_peak/max(tailing_peak)*100*sn[i]

                        if(is.null(tailingindex)){
                                # new peak with tailing
                                peak <- c(gaussian_peak[1:which.max(gaussian_peak)],tailing_peak[(which.max(gaussian_peak)+1):length(rtime0)])
                        }else{
                                # indexed peaks tailing
                                peak <- ifelse(i%in%tailingindex,c(gaussian_peak[1:which.max(gaussian_peak)],tailing_peak[(which.max(gaussian_peak)+1):length(rtime0)]),gaussian_peak)
                        }
                        re <- rbind(re,peak)
                }
                spd <-
                        S4Vectors::DataFrame(msLevel = 1L, rtime = rtime0)
                mz <- lapply(sub, function(x)
                        x$spectra$mz)
                intensity <- lapply(sub, function(x)
                        x$spectra$ins)
                if(!is.null(inscutoff)){
                        insidx <- lapply(intensity,function(x) (x/max(x))>inscutoff)
                        mz <- lapply(seq_along(mz),function(i) mz[[i]][insidx[[i]]])
                        intensity <- lapply(seq_along(intensity), function(i)
                                intensity[[i]][insidx[[i]]])
                }
                subname <- sapply(sub, function(x) x$name)
                mzv <- unlist(mz)
                intensityv <- unlist(intensity)
                rtimev <- rep(rtime,sapply(mz,length))
                namev <- rep(subname,sapply(mz,length))
                df2 <- cbind.data.frame(mz=mzv,rt=rtimev,ins=intensityv,name=namev)


                mzc <- rem <- c()

                for (i in 1:nrow(re)) {
                        if(length(mz[[i]])==1){

                                nret <- intensity[[i]]*re[i,]

                        }else{
                                nret <- matrix(intensity[[i]])%*%re[i,]
                        }
                        mzc <- c(mzc,round(mz[[i]],digits = mzdigit))
                        rem <- rbind(rem,nret)
                }
                alld <- stats::aggregate(rem,by=list(mzc),FUN=sum)
                mzpeak <- as.numeric(alld$Group.1)
                mzpeak <- mzpeak+stats::rnorm(length(mzpeak),sd=noisesd)*mzpeak*1e-6*ppm

                if(matrix){
                        mzm <- ifelse(is.null(matrixmz),mzm,matrixmz)
                        mzm <- round(mzm,digits = mzdigit)
                        mzmatrix <- mzm[!mzm%in%mzpeak]
                        insmatrix <- matrix(stats::rnorm(length(rtime0)*length(mzmatrix), mean = baseline, sd= baselinesd),nrow = length(mzmatrix),ncol = length(rtime0))
                        noisepeak <- matrix(stats::rnorm(length(rtime0)*length(mzpeak), mean = baseline, sd= baselinesd),nrow = length(mzpeak),ncol = length(rtime0))
                        inspeak <- alld[,-1]+noisepeak
                        allins <- rbind(inspeak,insmatrix)
                        # order mz
                        mz <- c(mzpeak,mzmatrix)
                        mz <- mz[order(mz)]
                        allins <- allins[order(mz),]
                }else{
                        mz <- mzpeak
                        noise <- matrix(stats::rnorm(length(rtime0)*length(mz), mean = baseline, sd= baselinesd),nrow = length(mz),ncol = length(rtime0))
                        allins <- alld[,-1]+noise
                }

                mzl <- intensityl <- list()
                for(i in 1:length(rtime0)){
                        ins <- allins[,i]
                        intensityl[[i]] <- ins[ins>0&mz>mzrange[1]&mz<mzrange[2]]
                        mzt <- mz[ins>0&mz>mzrange[1]&mz<mzrange[2]]
                        # add ppm shift for each scan
                        mzl[[i]] <- mzt+stats::rnorm(length(mzt),sd=noisesd)*mzt*1e-6*sampleppm
                }

                spd$mz <- mzl
                spd$intensity <- intensityl

                sp0 <- Spectra::Spectra(spd)
                sp0$scanIndex <- seq_along(rtime0)
                Spectra::export(sp0, Spectra::MsBackendMzR(), file = paste0(name, '.mzML'))
                utils::write.csv(df2,file = paste0(name, '.csv'),row.names = F)
        }

#' Generate figure for mzml raw data
#' @param mzml mzML file path
#' @param mzrange m/z range for simulation, peaks out of the range will be removed, default c(100,1000)
#' @param rtrange retention time range for simulation, default c(0,600)
#' @return figure for raw data and dataframe with mz, rt, and intensity
#' @export
mzmlviz <-
        function(mzml,
                 mzrange = c(0, 600),
                 rtrange = c(100, 1000)) {
                dt <- Spectra::Spectra(mzml)
                rt <- Spectra::rtime(dt)
                ins <- Spectra::intensity(dt)
                mz <- Spectra::mz(dt)
                mzv <- unlist(mz)
                rtimev <- rep(rt, times = sapply(mz, length))
                intensityv <- unlist(ins)
                norm <-
                        (intensityv - min(intensityv)) / (max(intensityv) - min(intensityv))
                # png('test.png',width = diff(xlim),height = diff(ylim))
                # par(mar = c(0, 0, 0, 0))
                # plot(rtimev,mzv,xlim=xlim,ylim=ylim,pch=19,cex=0.01,col=gray(1-norm), axes = FALSE,xaxs = "i", yaxs = "i")
                # dev.off()
                plot(
                        rtimev,
                        mzv,
                        xlim = mzrange,
                        ylim = rtrange,
                        pch = 19,
                        cex = 0.01,
                        col = grDevices::gray(1 - norm),
                        xlab = 'retention time(s)',
                        ylab = 'm/z'
                )
                return(cbind.data.frame(mz=mzv,rt=rtimev,intensity=intensityv))
        }
