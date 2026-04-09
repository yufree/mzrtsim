#' A vector containing m/z from matrix blank
#' @docType data
#' @usage data(mzm)
#' @format A vector containing m/z from matrix blank for LC-HRMS simulation
"mzm"

# Internal mzML writer using base64enc (no mzR/Spectra dependency)
# Supports mixed MS1/MS2/MSn spectra with precursor information.
#
# @param mz_list list of numeric vectors (m/z per scan)
# @param intensity_list list of numeric vectors (intensity per scan)
# @param rtime numeric vector of retention times in seconds
# @param file output file path
# @param ms_level integer vector of MS levels per scan, or scalar (default 1L)
# @param precursor_mz numeric vector of precursor m/z per scan (NA for MS1)
# @param precursor_charge integer vector of charge states (optional, NA ok)
# @param precursor_intensity numeric vector of precursor intensities (optional)
# @param collision_energy numeric vector of collision energies (optional, NA ok)
# @param isolation_window numeric scalar or vector for half-width of
#   isolation window in m/z (default 0.7)
# @param activation character activation method: "CID", "HCD", "ETD",
#   "ECD", "IRMPD", or "SID" (default "CID"). Can be a vector per scan.
# @param centroided logical, label spectra as centroid (TRUE, default) or
#   profile (FALSE). Can be a vector per scan.
write_mzml <- function(mz_list, intensity_list, rtime, file,
                       ms_level = 1L,
                       precursor_mz = NULL,
                       precursor_charge = NULL,
                       precursor_intensity = NULL,
                       collision_energy = NULL,
                       isolation_window = 0.7,
                       activation = "CID",
                       centroided = TRUE) {
    nscans <- length(rtime)

    # Recycle scalars to per-scan vectors
    ms_level <- rep_len(as.integer(ms_level), nscans)
    centroided <- rep_len(centroided, nscans)
    if (!is.null(precursor_mz))
        precursor_mz <- rep_len(precursor_mz, nscans)
    if (!is.null(precursor_charge))
        precursor_charge <- rep_len(as.integer(precursor_charge), nscans)
    if (!is.null(precursor_intensity))
        precursor_intensity <- rep_len(precursor_intensity, nscans)
    if (!is.null(collision_energy))
        collision_energy <- rep_len(collision_energy, nscans)
    isolation_window <- rep_len(isolation_window, nscans)
    activation <- rep_len(activation, nscans)

    # Activation CV lookup
    act_cv <- c(
        CID   = "MS:1000133", HCD   = "MS:1000422",
        ETD   = "MS:1000598", ECD   = "MS:1000250",
        IRMPD = "MS:1000262", SID   = "MS:1000136"
    )
    act_name <- c(
        CID   = "collision-induced dissociation",
        HCD   = "beam-type collision-induced dissociation",
        ETD   = "electron transfer dissociation",
        ECD   = "electron capture dissociation",
        IRMPD = "infrared multiphoton dissociation",
        SID   = "surface-induced dissociation"
    )

    encode_base64 <- function(x) {
        base64enc::base64encode(
            writeBin(as.double(x), raw(), size = 8, endian = "little")
        )
    }

    has_ms1 <- any(ms_level == 1L)
    has_msn <- any(ms_level >= 2L)

    con <- file(file, "w")
    on.exit(close(con))

    # fileContent CV terms
    fc_lines <- character()
    if (has_ms1)
        fc_lines <- c(fc_lines,
            '      <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>')
    if (has_msn)
        fc_lines <- c(fc_lines,
            '      <cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>')

    # Header
    writeLines(c(
        '<?xml version="1.0" encoding="utf-8"?>',
        '<mzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd" id="mzrtsim_output" version="1.1.0">',
        '  <cvList count="2">',
        '    <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" version="4.1.30" URI="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/>',
        '    <cv id="UO" fullName="Unit Ontology" version="09:04:2014" URI="https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo"/>',
        '  </cvList>',
        '  <fileDescription>',
        '    <fileContent>',
        fc_lines,
        '    </fileContent>',
        '  </fileDescription>',
        '  <softwareList count="1">',
        '    <software id="mzrtsim" version="0.99.0">',
        '      <cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value="mzrtsim"/>',
        '    </software>',
        '  </softwareList>',
        '  <instrumentConfigurationList count="1">',
        '    <instrumentConfiguration id="IC">',
        '      <cvParam cvRef="MS" accession="MS:1000031" name="instrument model" value=""/>',
        '    </instrumentConfiguration>',
        '  </instrumentConfigurationList>',
        '  <dataProcessingList count="1">',
        '    <dataProcessing id="dp">',
        '      <processingMethod order="0" softwareRef="mzrtsim">',
        '        <cvParam cvRef="MS" accession="MS:1000544" name="Conversion to mzML" value=""/>',
        '      </processingMethod>',
        '    </dataProcessing>',
        '  </dataProcessingList>',
        '  <run id="run1" defaultInstrumentConfigurationRef="IC">',
        sprintf('    <spectrumList count="%d" defaultDataProcessingRef="dp">', nscans)
    ), con)

    # Spectra
    for (i in seq_len(nscans)) {
        mz <- mz_list[[i]]
        int <- intensity_list[[i]]
        n <- length(mz)
        mz_b64 <- encode_base64(mz)
        int_b64 <- encode_base64(int)
        lvl <- ms_level[i]

        # Spectrum type CV
        if (lvl == 1L) {
            type_cv <- '        <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>'
        } else {
            type_cv <- '        <cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>'
        }
        # Centroid / profile
        if (centroided[i]) {
            mode_cv <- '        <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>'
        } else {
            mode_cv <- '        <cvParam cvRef="MS" accession="MS:1000128" name="profile spectrum" value=""/>'
        }

        spec_lines <- c(
            sprintf('      <spectrum index="%d" id="scan=%d" defaultArrayLength="%d">', i - 1L, i, n),
            sprintf('        <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="%d"/>', lvl),
            type_cv,
            mode_cv,
            '        <scanList count="1">',
            '          <cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/>',
            '          <scan>',
            sprintf('            <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="%f" unitCvRef="UO" unitAccession="UO:0000010" unitName="second"/>', rtime[i]),
            '          </scan>',
            '        </scanList>'
        )

        # Precursor list for MS2+
        if (lvl >= 2L && !is.null(precursor_mz) && !is.na(precursor_mz[i])) {
            pmz <- precursor_mz[i]
            iw <- isolation_window[i]
            act_key <- activation[i]
            act_acc <- act_cv[act_key]
            act_nm  <- act_name[act_key]

            prec_lines <- c(
                '        <precursorList count="1">',
                '          <precursor>',
                '            <isolationWindow>',
                sprintf('              <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="%f" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>', pmz),
                sprintf('              <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="%f" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>', iw),
                sprintf('              <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="%f" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>', iw),
                '            </isolationWindow>',
                '            <selectedIonList count="1">',
                '              <selectedIon>'
            )
            prec_lines <- c(prec_lines,
                sprintf('                <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="%f" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>', pmz))
            if (!is.null(precursor_charge) && !is.na(precursor_charge[i]))
                prec_lines <- c(prec_lines,
                    sprintf('                <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="%d"/>', precursor_charge[i]))
            if (!is.null(precursor_intensity) && !is.na(precursor_intensity[i]))
                prec_lines <- c(prec_lines,
                    sprintf('                <cvParam cvRef="MS" accession="MS:1000042" name="peak intensity" value="%f"/>', precursor_intensity[i]))
            prec_lines <- c(prec_lines,
                '              </selectedIon>',
                '            </selectedIonList>',
                '            <activation>',
                sprintf('              <cvParam cvRef="MS" accession="%s" name="%s" value=""/>', act_acc, act_nm)
            )
            if (!is.null(collision_energy) && !is.na(collision_energy[i]))
                prec_lines <- c(prec_lines,
                    sprintf('              <cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="%f" unitCvRef="UO" unitAccession="UO:0000266" unitName="electronvolt"/>', collision_energy[i]))
            prec_lines <- c(prec_lines,
                '            </activation>',
                '          </precursor>',
                '        </precursorList>'
            )
            spec_lines <- c(spec_lines, prec_lines)
        }

        # Binary data arrays
        spec_lines <- c(spec_lines,
            '        <binaryDataArrayList count="2">',
            sprintf('          <binaryDataArray encodedLength="%d">', nchar(mz_b64)),
            '            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value=""/>',
            '            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>',
            '            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>',
            sprintf('            <binary>%s</binary>', mz_b64),
            '          </binaryDataArray>',
            sprintf('          <binaryDataArray encodedLength="%d">', nchar(int_b64)),
            '            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value=""/>',
            '            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>',
            '            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>',
            sprintf('            <binary>%s</binary>', int_b64),
            '          </binaryDataArray>',
            '        </binaryDataArrayList>',
            '      </spectrum>'
        )

        writeLines(spec_lines, con)
    }

    # Footer
    writeLines(c(
        '    </spectrumList>',
        '  </run>',
        '</mzML>'
    ), con)
}

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
#' @param pwidth peak width for the compound. If it's a single value, simulated peaks' width will use this number as peak width for all compounds. If it's a numeric vector, it will be used as the peak widths for every compounds.
#' @param pheight peak height for the compound. If it's a single value, simulated peaks' height will use this number as peak height for all compounds. If it's a numeric vector, it will be used as the peak heights for every compounds.
#' @param baseline noise baseline, default 100. If it's a numeric vector, it will be used as the baseline for each scan.
#' @param baselinesd standard deviation for noise, default 30
#' @param rf response factor of each compound, default 100 for all compounds
#' @param tailingfactor tailing factor for peaks, larger means larger tailing, default 1.2
#' @param compound numeric compound index in the database for targeted analysis, default NULL
#' @param rtime retention time for the compounds, if NULL, retention time will be simulated. Default NULL
#' @param tailingindex numeric index for tailing compounds, if NULL, all peaks will tailing. Default NULL
#' @param seed Random seed for reproducibility
#' @param unique if TRUE, one compound will have one spectra. Default FALSE
#' @param matrix if TRUE, m/z from experimental data will be used for background m/z simulation.Default FALSE
#' @param matrixmz custom matrix m/z vector, default NULL and predefined list from serum blank will be used.
#' @param background mzML file path for background simulation, default NULL
#' @param ms_level integer vector of MS levels per scan, or scalar (default 1L)
#' @param precursor_mz numeric vector of precursor m/z per scan (NA for MS1)
#' @param precursor_charge integer vector of charge states (optional, NA ok)
#' @param precursor_intensity numeric vector of precursor intensities (optional)
#' @param collision_energy numeric vector of collision energies (optional, NA ok)
#' @param isolation_window numeric scalar or vector for half-width of
#'   isolation window in m/z (default 0.7)
#' @param activation character activation method: "CID", "HCD", "ETD",
#'   "ECD", "IRMPD", or "SID" (default "CID"). Can be a vector per scan.
#' @return one mzML file for simulated data and one csv file the simulated compounds with retention time, m/z and name
#' @export
#' @examples
#' data(monams1)
#' simmzml(db=monams1, name = file.path(tempdir(), 'test'))
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
                 pheight = 10,
                 baseline = 100,
                 baselinesd = 30,
                 rf = 100,
                 tailingfactor = 1.2,
                 compound=NULL,
                 rtime=NULL,
                 tailingindex = NULL,
                 seed=42,
                 unique=FALSE,
                 matrix=FALSE,
                 matrixmz=NULL,
                 background=NULL,
                 ms_level = 1L,
                 precursor_mz = NULL,
                 precursor_charge = NULL,
                 precursor_intensity = NULL,
                 collision_energy = NULL,
                 isolation_window = 0.7,
                 activation = "CID") {
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

                if(!is.null(background)){
                        rawdata <- RaMS::grabMSdata(background, grab_what = "MS1", verbosity = 0L)$MS1
                        rawdata$rt <- rawdata$rt * 60
                        rtime0 <- unique(rawdata$rt)
                }else{
                        rtime0 <- seq(rtrange[1], rtrange[2], scanrate)
                }

                if(length(baseline) != 1 && length(baseline) != length(rtime0)){
                        stop('Length of baseline must be 1 or equal to the number of scans (length of rtime0).')
                }

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
                        peakrange <- rep(pwidth,n)
                }else{
                        peakrange <- pwidth
                }
                if(length(pheight)==1){
                        set.seed(seed)
                        # uniform distributed peak width
                        peakheight <- rep(pheight,n)
                }else{
                        peakheight <- pheight
                }
                if(length(rf)==1){
                        rf <- rep(rf,100)
                }else if(length(rf)!=n){
                        set.seed(seed)
                        rf <- sample(rf,n,replace = T)
                }
                # chromotograghy simulation for the compound
                re <- c()
                if(n>0){
                        for (i in 1:n) {
                                gaussian_peak <-
                                        stats::dnorm(rtime0, mean = rtime[i], sd = peakrange[i] / 4)
                                gaussian_peak <- gaussian_peak/max(gaussian_peak)*100*rf[i]*peakheight[i]
                                # tailing simulation
                                tailing_peak <-
                                        stats::dnorm(rtime0, mean = rtime[i], sd = (2*tailingfactor-1)*peakrange[i] / 4)
                                tailing_peak <- tailing_peak/max(tailing_peak)*100*rf[i]*peakheight[i]

                                if(is.null(tailingindex)){
                                        # new peak with tailing
                                        peak <- c(gaussian_peak[1:which.max(gaussian_peak)],tailing_peak[(which.max(gaussian_peak)+1):length(rtime0)])
                                }else{
                                        # indexed peaks tailing
                                        peak <- ifelse(i%in%tailingindex,c(gaussian_peak[1:which.max(gaussian_peak)],tailing_peak[(which.max(gaussian_peak)+1):length(rtime0)]),gaussian_peak)
                                }
                                re <- rbind(re,peak)
                        }
                }
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

                # Calculate simulated total max intensities
                if(n>0){
                        sim_insv <- unlist(lapply(1:n, function(i) intensity[[i]] * max(re[i,])))
                } else {
                        sim_insv <- numeric(0)
                }

                df2 <- cbind.data.frame(mz=mzv,rt=rtimev,ins=intensityv,sim_ins=sim_insv,name=namev)


                mzc <- rem <- c()

                if(n>0){
                        for (i in 1:nrow(re)) {
                                if(length(mz[[i]])==1){

                                        nret <- intensity[[i]]*re[i,]

                                }else{
                                        nret <- matrix(intensity[[i]])%*%re[i,]
                                }
                                mzc <- c(mzc,round(mz[[i]],digits = mzdigit))
                                rem <- rbind(rem,nret)
                        }
                }
                if(length(mzc) > 0){
                        alld <- stats::aggregate(rem,by=list(mzc),FUN=sum)
                        mzpeak <- as.numeric(alld$Group.1)
                        mzpeak <- mzpeak+stats::rnorm(length(mzpeak),sd=noisesd)*mzpeak*1e-6*ppm
                }else{
                        mzpeak <- numeric()
                        alld <- matrix(nrow = 0, ncol = length(rtime0)+1)
                }

                if(!is.null(background)){
                        mz <- mzpeak
                        allins <- alld[,-1]
                        mzl <- intensityl <- list()
                        rt_factor <- factor(rawdata$rt, levels = rtime0)
                        raw_by_scan <- split(seq_len(nrow(rawdata)), rt_factor)
                        rawmz <- lapply(raw_by_scan, function(idx) rawdata$mz[idx])
                        rawint <- lapply(raw_by_scan, function(idx) rawdata$int[idx])

                        for(i in 1:length(rtime0)){
                                ins <- allins[,i]
                                intensityl[[i]] <- ins[ins>0&mz>mzrange[1]&mz<mzrange[2]]
                                mzt <- mz[ins>0&mz>mzrange[1]&mz<mzrange[2]]
                                # add ppm shift for each scan
                                mzsim <- mzt+stats::rnorm(length(mzt),sd=noisesd)*mzt*1e-6*sampleppm
                                intsim <- intensityl[[i]]

                                mzt <- c(rawmz[[i]], mzsim)
                                intt <- c(rawint[[i]], intsim)
                                order <- order(mzt)
                                mzl[[i]] <- mzt[order]
                                intensityl[[i]] <- intt[order]
                        }
                }else{
                        if(matrix){
                                if(is.null(matrixmz)){
                                        mzm_env <- new.env(parent = emptyenv())
                                        utils::data("mzm", envir = mzm_env)
                                        mzm_local <- mzm_env$mzm
                                }else{
                                        mzm_local <- matrixmz
                                }
                                mzm_local <- round(mzm_local,digits = mzdigit)
                                mzmatrix <- mzm_local[!mzm_local%in%mzpeak]
                                insmatrix <- matrix(stats::rnorm(length(rtime0)*length(mzmatrix), mean = rep(baseline, each=length(mzmatrix)), sd= baselinesd),nrow = length(mzmatrix),ncol = length(rtime0))

                                if(length(mzpeak) > 0){
                                        noisepeak <- matrix(stats::rnorm(length(rtime0)*length(mzpeak), mean = rep(baseline, each=length(mzpeak)), sd= baselinesd),nrow = length(mzpeak),ncol = length(rtime0))
                                        inspeak <- alld[,-1]+noisepeak
                                }else{
                                        inspeak <- matrix(nrow = 0, ncol = length(rtime0))
                                }

                                allins <- rbind(as.matrix(inspeak),insmatrix)
                                # order mz
                                mz <- c(mzpeak,mzmatrix)
                                idx <- order(mz)
                                mz <- mz[idx]
                                allins <- allins[idx,]
                        }else{
                                mz <- mzpeak
                                if(length(mz) > 0){
                                        noise <- matrix(stats::rnorm(length(rtime0)*length(mz), mean = rep(baseline, each=length(mz)), sd= baselinesd),nrow = length(mz),ncol = length(rtime0))
                                        allins <- alld[,-1]+noise
                                }else{
                                        allins <- matrix(nrow = 0, ncol = length(rtime0))
                                }
                        }

                        mzl <- intensityl <- list()
                        for(i in 1:length(rtime0)){
                                ins <- allins[,i]
                                intensityl[[i]] <- ins[ins>0&mz>mzrange[1]&mz<mzrange[2]]
                                mzt <- mz[ins>0&mz>mzrange[1]&mz<mzrange[2]]
                                # add ppm shift for each scan
                                mzl[[i]] <- mzt+stats::rnorm(length(mzt),sd=noisesd)*mzt*1e-6*sampleppm
                        }
                }

                write_mzml(mzl, intensityl, rtime0, paste0(name, '.mzML'),
                           ms_level = ms_level,
                           precursor_mz = precursor_mz,
                           precursor_charge = precursor_charge,
                           precursor_intensity = precursor_intensity,
                           collision_energy = collision_energy,
                           isolation_window = isolation_window,
                           activation = activation)
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
                 mzrange = c(100, 1000),
                 rtrange = c(0, 600)) {
                msdata <- RaMS::grabMSdata(mzml, grab_what = "MS1", verbosity = 0L)$MS1
                mzv <- msdata$mz
                rtimev <- msdata$rt * 60
                intensityv <- log10(msdata$int + 1)
                norm <-
                        (intensityv - min(intensityv)) / (max(intensityv) - min(intensityv))
                # png('test.png',width = diff(xlim),height = diff(ylim))
                # par(mar = c(0, 0, 0, 0))
                # plot(rtimev,mzv,xlim=xlim,ylim=ylim,pch=19,cex=0.01,col=gray(1-norm), axes = FALSE,xaxs = "i", yaxs = "i")
                # dev.off()
                plot(
                        rtimev,
                        mzv,
                        xlim = rtrange,
                        ylim = mzrange,
                        pch = 15,
                        cex = 0.1,
                        col = grDevices::gray(1 - norm),
                        xlab = 'retention time(s)',
                        ylab = 'm/z'
                )
                return(cbind.data.frame(mz=mzv,rt=rtimev,intensity=intensityv))
        }

#' Generate simulated mzml data for blank sample
#' @param name file name of mzml
#' @param mzrange m/z range for simulation, peaks out of the range will be removed, default c(100,1000)
#' @param rtrange retention time range for simulation, default c(0,600)
#' @param scanrate time for each full scan, default 0.2 second or 5 spectra per second
#' @param baseline noise baseline, default 100. If it's a numeric vector, it will be used as the baseline for each scan.
#' @param baselinesd standard deviation for noise, default 30
#' @param mzdigit m/z digits, default 5
#' @param matrixmz custom matrix m/z vector, default NULL and predefined list from serum blank will be used.
#' @return one mzML file for simulated blank data
#' @export
#' @examples
#' \donttest{
#' simmzml_blank(name = file.path(tempdir(), 'blank_test'))
#' }
simmzml_blank <- function(name,
                          mzrange = c(100,1000),
                          rtrange = c(0, 600),
                          scanrate = 0.2,
                          baseline = 100,
                          baselinesd = 30,
                          mzdigit = 5,
                          matrixmz = NULL) {
    
    # Create a dummy database with a single fake compound to bypass the required argument
    dummy_db <- list(list(name="dummy", spectra=list(mz=c(100), ins=c(100))))
    
    # Call simmzml with 0 compounds and matrix=TRUE
    simmzml(db = dummy_db,
            name = name,
            n = 0,
            mzrange = mzrange,
            rtrange = rtrange,
            scanrate = scanrate,
            baseline = baseline,
            baselinesd = baselinesd,
            mzdigit = mzdigit,
            matrix = TRUE,
            matrixmz = matrixmz)
}
