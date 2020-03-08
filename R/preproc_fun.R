#' A Function to preprocess Sentinel-1 data using SNAP
#'
#' This function allows you to preprocess Sentinel-1 data in batch mode, where the S1 data should be downloaded as ZIP file before using this function.
#' The S1 data will be processed through the following steps in sequence:
#' Read --> Apply Orbit File --> Remove Thermal Noise --> Calibration --> Speckle Filtering --> Terrain Correction --> Subset --> Convert to dB --> Write
#'
#' @param sourcePath Path to where the S1 data (.zip file) is located
#' @param wd Working directory of the outputs
#' @param selectedPolarisations =<string,string,string,...> The list of polarisations. Default to "VV,VH"
#' @param sourceBands =<string,string,string,...> The list of source bands. Default to "Sigma0_VV,Sigma0_VH"
#' @param pixelSpacingInMeter The pixel spacing in meters. This should be omitted if "pixelSpacingInDegree" is defined.
#' @param pixelSpacingInDegree The pixel spacing in degrees. This should be omitted if "pixelSpacingInMeter" is defined.
#' @param geoRegion The subset region in geographical coordinates, must be given as string charactor, e.g. "<lon1> <lat1>, <lon2> <lat2>, ..., <lon1> <lat1>". The first and the last coordinates should be identical to make the polygon closed. This should be omitted if "referenceBand" and "region" are defined.
#' @param referenceBand The band used to indicate pixel coordinates. This should be omitted if "geoRegion" is defined.
#' @param region The subset region in pixel coordinates. This should be omitted if "geoRegion" is defined.
#' @param output_BEAM_DIMAP If TRUE, the function will create output in BEAM-DIMAP(.dim) format. Default to TRUE
#' @param output_GeoTIFF If TRUE, the function will create output in GeoTIFF format. Default to TRUE
#'
#' @param continue_on_failure Default to "false".
#' @param orbit_type Value must be one of "Sentinel Precise (Auto Download)", "Sentinel Restituted (Auto Download)", "DORIS Precise VOR (ENVISAT) (Auto Download)", "DELFT Precise (ENVISAT, ERS1&2) (Auto Download)", "PRARE Precise (ERS1&2) (Auto Download)". Default to "Sentinel Precise (Auto Download)".
#' @param poly_degree Default to 3
#' @param reIntroduceThermalNoise Re-introduce thermal noise. Default to "false"
#' @param removeThermalNoise Remove thermal noise. Default to "true"
#' @param auxFile The auxiliary file. Value must be one of "Latest Auxiliary File", "Product Auxiliary File", "External Auxiliary File". Default to "Latest Auxuliary File"
#' @param createBetaBand Create beta0 virtual band. Default to "false"
#' @param createGammaBand Create gamma0 virtual band. Default to "false"
#' @param outputBetaBand Output beta0 band. Default to "false"
#' @param outputGammaBand Output gamma0 band. Default to "false"
#' @param outputSigmaBand Output sigma0 band. Default to "true"
#' @param outputImageInComplex Output image in complex. Default to "false"
#' @param outputImageScaleInDb Output image scale. Default to "false"
#' @param anSize The Adaptive Neighbourhood size. Valid internal is (1,200]. Default to 50
#' @param dampingFactor The damping factor(Frost filter only). Valid interval is (0,100]. Default to 2
#' @param enl The number of looks in Speckle Filtering. Valid interval is (0, *). Default to 1.0
#' @param estimateENL Default to "false"
#' @param filter Value must be one of "None", "Boxcar", "Median", "Frost", "Gamma Map", "Lee", "Refined Lee", "Lee Sigma", "IDAN". Default to "Lee Sigma"
#' @param filterSizeX The kernel x dimension. Valid interval is (1, 100]. Default to 3
#' @param filterSizeY the kernel y dimension. Valid interval is (1, 100]. Default to 3
#' @param numLooksStr Value must be one of "1", "2", "3", "4". Default to "1"
#' @param sigmaStr Value must be one of "0.5", "0.6", "0.7", "0.8", "0.9". Default to "0.9"
#' @param targetWindowSizeStr Value must be one of "3x3", "5x5". Default to "3x3".
#' @param windowSize Value must be one of "5x5", "7x7", "9x9", "11x11", "13x13", "15x15", "17x17". Default to "7x7"
#' @param mapProjection The coordinate reference system in well knewn text format. Default value is "WGS84(DD)"
#' @param saveSelectedSourceBand Default to "true"
#' @param alignToStandardGrid Force the image grid to be aligned with a specific point. Default to "false"
#' @param applyRadiometricNormalization Default to "false"
#' @param demName The digital elevation model. Default to "SRTM 3sec"
#' @param demResamplingMethod Value must be one of "NEAREST_NEIGHBOUR", "BILINEAR_INTERPOLATION", "CUBIC_CONVOLUTION", "BISINC_5_POINT_INTERPOLATION", "BISIC_11_POINT_INTERPOLATION", "BISIC_21_POINT_INTERPOLATION", "BICUBIC_INTERPOLATION", "DELAUNAY_INTERPOLATION". Default to "BILINEAR_INTERPOLATION"
#' @param externalDEMApplyEGM Default to "true"
#' @param externalDEMNoDataValue Default to 0
#' @param imgResamplingMethod Value must be one of "NEAREST_NEIGHBOUR", "BILINEAR_INTERPOLATION", "CUBIC_CONVOLUTION", "BISINC_5_POINT_INTERPOLATION", "BISIC_11_POINT_INTERPOLATION", "BISIC_21_POINT_INTERPOLATION", "BICUBIC_INTERPOLATION", "DELAUNAY_INTERPOLATION". Default to "BILINEAR_INTERPOLATION"
#' @param incidenceAngleForSigma0 Value must be one of "Use incidence angle from Ellipsoid", "Use local incidence angle from DEM", "Use projected local incidence angle from DEM". Default to "Use projected local incidence angle from DEM"
#' @param nodataValueAtSea Mask the sea with no data value (faster). Default to "true"
#' @param outputComplex Default to "false"
#' @param saveBetaNought Default to "false"
#' @param saveGammaNought Default to "false"
#' @param saveDEM Default to "false"
#' @param saveIncidenceAngleFromEllipsoid Default to "false"
#' @param saveLatLon Default to "false"
#' @param saveLocalIncidenceAngle Default to "false"
#' @param saveProjectedLocalIncidenceAngle Default to "false"
#' @param saveSigmaNought Default to "false"
#' @param standardGridOriginX x-coordinate of the standard grid's origin point. Default to 0
#' @param standardGridOriginY y-coordinate of the standard grid's origin point. Default to 0
#' @param copyMetadata Whether to copy the metadata of the source product. Default to "false"
#' @param fullSwath Forces the operator to extend the subset region to the full swath. Default to "false"
#' @param subSamplingX The pixel sub-sampling step in X (horizontal image direction). Default to 1
#' @param subSamplingY The pixel sub-sampling step in Y (vertical image direction). Default to 1
#' @param clear_Cache_after_row_write If true, the internal tile cache is cleared after a tile row has been written. Ignored if writeEntireTileRows=false. Default to "false"
#' @param delete_output_on_failure If true, all output files are deleted after a failed write operation. Default to "true"
#' @param write_entire_tile_rows If true, the write operation waits until an entire tile row is computed. Default to "true"
#'
#' @keywords rSNAP
#' @export
#' @examples
#' preprocessing(sourcePath = "<path to source product>",
#'               wd = "<path to directory for outputs>",
#'               selectedPolarisations = "VH",
#'               sourceBands = "Sigma0_VH",
#'               pixelSpacingInMeter=10,
#'               geoRegion = "100.01839447021484 26.000097274780273, 100.41411590576172 26.000097274780273, 100.41411590576172 25.508270263671875, 100.01839447021484 25.508270263671875, 100.01839447021484 26.000097274780273")
#'


preprocessing <- function(sourcePath, wd, selectedPolarisations="VV,VH", sourceBands="Sigma0_VV,Sigma0_VH",
                          pixelSpacingInMeter, pixelSpacingInDegree, geoRegion, referenceBand, region,
                          output_BEAM_DIMAP=TRUE, output_GeoTIFF=TRUE,
                          # Parameters to be specified

                          #source, wd,
                          continue_on_failure="false", orbit_type="Sentinel Precise (Auto Download)", poly_degree=3,
                          # Orbit file

                          #source, wd,
                          reIntroduceThermalNoise="false", removeThermalNoise="true",
                          # remove Thermal Noise

                          #source, wd, selectedPolarisations,
                          auxFile="Latest Auxiliary File",
                          createBetaBand="false", createGammaBand="false",
                          outputBetaBand="false", outputGammaBand="false", outputSigmaBand="true",
                          outputImageInComplex="false", outputImageScaleInDb="false",
                          # Calibration

                          #source, wd,
                          anSize=50,
                          dampingFactor=2, enl=1.0, estimateENL="false",
                          filter="Lee Sigma", filterSizeX=3,
                          filterSizeY=3, numLooksStr=1, sigmaStr=0.9,
                          targetWindowSizeStr="3x3", windowSize="7x7",
                          # Speckle Filtering

                          #source, wd, sourceBands="Sigma0_VV,Sigma0_VH", auxFile="Latest Auxiliary File",
                          mapProjection="WGS84(DD)",
                          saveSelectedSourceBand="true", alignToStandardGrid="false",
                          applyRadiometricNormalization="false",
                          demName="SRTM 3Sec",
                          demResamplingMethod="BILINEAR_INTERPOLATION", externalDEMApplyEGM="true",
                          externalDEMNoDataValue=0, imgResamplingMethod="BILINEAR_INTERPOLATION",
                          incidenceAngleForSigma0="Use projected local incidence angle from DEM",
                          nodataValueAtSea="true", outputComplex="false", saveBetaNought="false",
                          saveGammaNought="false", saveDEM="false", saveIncidenceAngleFromEllipsoid="false",
                          saveLatLon="false", saveLocalIncidenceAngle="false",
                          saveProjectedLocalIncidenceAngle="false", saveSigmaNought="false",
                          standardGridOriginX=0, standardGridOriginY=0,
                          # Terrain Correction

                          #source, wd, sourceBands,
                          copyMetadata="false", fullSwath="false", subSamplingX=1, subSamplingY=1,
                          # Subset

                          #source, wd, sourceBands,
                          # Linear to db

                          #source, wd,
                          clear_Cache_after_row_write="false",
                          delete_output_on_failure="true",
                          write_entire_tile_rows="true"
                          # Write
) {

  # Create subdirectories
  dir_name <- c("01_read", "02_appOrbit", "03_rmThmNoise", "04_calib", "05_speckFilt", "06_terCor", "07_subset", "08_toDb", "Output")
  dir_paths <- paste0(wd, "/", dir_name)
  myDir <- lapply(dir_paths, function(x){
    dir.create(x, showWarnings = FALSE)
  })
  source <- list.files(path = sourcePath, pattern = "*.zip")
  # Listing all ".zip" files in the source path


  for (i in 1:length(source)) {

    # Read
    input <- paste0("-Pfile=", sourcePath, "/", source[i])
    # Input is the i-th file in the source

    run(command = "gpt",
        args = c("Read",
                 input),
        wd = c(dir_paths[1]),
        echo_cmd = TRUE)


    # Apply-Orbit-File
    input <- paste0("-Ssource=", dir_paths[1], "/target.dim")
    cont_on_fail <- paste0("-PcontinueOnFail=", continue_on_failure)
    orbitType <- paste0("-PorbitType=", orbit_type)
    polyDegree <- paste0("-PpolyDegree=", poly_degree)

    run(command = "gpt",
        args = c("Apply-Orbit-File",
                 input,
                 cont_on_fail,
                 orbitType,
                 polyDegree),
        wd = c(dir_paths[2]),
        echo_cmd = TRUE)


    # Remove Thermal Noise
    input <- paste0("-SsourceProduct=", dir_paths[2], "/target.dim")
    reIntrThermNoi <- paste0("-PreIntroduceThermalNoise=", reIntroduceThermalNoise)
    rmThermNoi <- paste0("-PremoveThermalNoise=", removeThermalNoise)
    selPol <- paste0("-PselectedPolarisations=", selectedPolarisations)

    run(command = "gpt",
        args = c("ThermalNoiseRemoval",
                 input,
                 reIntrThermNoi,
                 rmThermNoi,
                 selPol),
        wd = c(dir_paths[3]),
        echo_cmd = TRUE)


    # Calibration
    input <- paste0("-Ssource=", dir_paths[3], "/target.dim")
    selPol <- paste0("-PselectedPolarisations=", selectedPolarisations)
    auxfile <- paste0("-PauxFile=", auxFile)
    crtBBand <- paste0("-PcreateBetaBand=", createBetaBand)
    crtGBand <- paste0("-PcreateGammaBand=", createGammaBand)
    opBBand <- paste0("-PoutputBetaBand=", outputBetaBand)
    opGBand <- paste0("-PoutputGammaBand=", outputGammaBand)
    opImgInCom <- paste0("-PoutputImageInComplex=", outputImageInComplex)
    opImgScaInDb <- paste0("-PoutputImageScaleInDb=", outputImageScaleInDb)
    opSBand <- paste0("-PoutputSigmaBand=", outputSigmaBand)

    run(command = "gpt",
        args = c("Calibration",
                 input,
                 auxfile,
                 crtBBand,
                 crtGBand,
                 opBBand,
                 opGBand,
                 opImgInCom,
                 opImgScaInDb,
                 opSBand,
                 selPol),
        wd=c(dir_paths[4]),
        echo_cmd = TRUE)


    # Speckle-Filter
    input <- paste0("-Ssource=", dir_paths[4], "/target.dim")
    adpSize <- paste0("-PanSize=", anSize)
    dampFac <- paste0("-PdampingFactor=", dampingFactor)
    Enl <- paste0("-Penl=", enl)
    estENL <- paste0("-PestimateENL=", estimateENL)
    flt <- paste0("-Pfilter=", filter)
    fltSzX <- paste0("-PfilterSizeX=", filterSizeX)
    fltSzY <- paste0("-PfilterSizeY=", filterSizeY)
    nlook <- paste0("-PnumLooksStr=", numLooksStr)
    tarWinSzStr <- paste0("-PtargetWindowSizeStr=", targetWindowSizeStr)
    winSz <- paste0("-PwindowSize=", windowSize)
    sourceband <- paste0("-PsourceBands=", sourceBands)

    run(command = "gpt",
        args = c("Speckle-Filter",
                 input,
                 adpSize,
                 dampFac,
                 Enl,
                 estENL,
                 flt,
                 fltSzX,
                 fltSzY,
                 nlook,
                 tarWinSzStr,
                 winSz,
                 sourceband),
        wd=c(dir_paths[5]),
        echo_cmd = TRUE)


    # Terrain Correction
    input <- paste0("-Ssource=", dir_paths[5], "/target.dim")
    sourceband <- paste0("-PsourceBands=", sourceBands)
    mapPrj <- paste0("-PmapProjection=", mapProjection)
    svSltBand <- paste0("-PsaveSelectedSourceBand=", saveSelectedSourceBand)
    algnStdGrd <- paste0("-PalignToStandardGrid=", alignToStandardGrid)
    aplRadNorm <- paste0("-PapplyRadiometricNormalization=", applyRadiometricNormalization)
    axF <- paste0("-PauxFile=", auxFile)
    demnm <- paste0("-PdemName=", demName)
    demResMeth <- paste0("-PdemResamplingMethod=", demResamplingMethod)
    extDEMappEGM <- paste0("-PexternalDEMApplyEGM=", externalDEMApplyEGM)
    extDEMNoDat <- paste0("-PexternalDEMNoDataValue=", externalDEMNoDataValue)
    ImgResMeth <- paste0("-PimgResamplingMethod=", imgResamplingMethod)
    incAngSigma <- paste0("-PincidenceAngleForSigma0=", incidenceAngleForSigma0)
    NoDataSea <- paste0("-PnodataValueAtSea=", nodataValueAtSea)
    opComp <- paste0("-PoutputComplex=", outputComplex)
    svBtNou <- paste0("-PsaveBetaNought=", saveBetaNought)
    svGmNou <- paste0("-PsaveGammaNought=", saveGammaNought)
    svDEM <- paste0("-PsaveDEM=", saveDEM)
    svIncAngfrEllip <- paste0("-PsaveIncidenceAngleFromEllipsoid=", saveIncidenceAngleFromEllipsoid)
    svLatLon <- paste0("-PsaveLatLon=", saveLatLon)
    svLocIncAng <- paste0("-PsaveLocalIncidenceAngle=", saveLocalIncidenceAngle)
    svLocPrjIncAng <- paste0("-PsaveProjectedLocalIncidenceAngle=", saveProjectedLocalIncidenceAngle)
    svSgNou <- paste0("-PsaveSigmaNought=", saveSigmaNought)
    stdGrdOriX <- paste0("-PstandardGridOriginX=", standardGridOriginX)
    stdGrdOriY <- paste0("-PstandardGridOriginY=", standardGridOriginY)

    if(missing(pixelSpacingInDegree)) {
      pxSpaInMet <- paste0("-PpixelSpacingInMeter=", pixelSpacingInMeter)
      run(command = "gpt",
          args = c("Terrain-Correction",
                   input,
                   sourceband,
                   pxSpaInMet,
                   mapPrj,
                   algnStdGrd,
                   aplRadNorm,
                   axF,
                   demnm,
                   demResMeth,
                   extDEMappEGM,
                   extDEMNoDat,
                   ImgResMeth,
                   incAngSigma,
                   NoDataSea,
                   opComp,
                   svBtNou,
                   svDEM,
                   svGmNou,
                   svIncAngfrEllip,
                   svLatLon,
                   svLocIncAng,
                   svLocPrjIncAng,
                   svSltBand,
                   svSgNou,
                   stdGrdOriX,
                   stdGrdOriY),
          wd=c(dir_paths[6]),
          echo_cmd = TRUE)
    }else if(missing(pixelSpacingInMeter)) {
      pxSpaInDeg <- paste0("-PpixelSpacingInDegree=", pixelSpacingInDegree)
      run(command = "gpt",
          args = c("Terrain-Correction",
                   input,
                   sourceband,
                   pxSpaInDeg,
                   mapPrj,
                   algnStdGrd,
                   aplRadNorm,
                   axF,
                   demnm,
                   demResMeth,
                   extDEMappEGM,
                   extDEMNoDat,
                   ImgResMeth,
                   incAngSigma,
                   NoDataSea,
                   opComp,
                   svBtNou,
                   svDEM,
                   svGmNou,
                   svIncAngfrEllip,
                   svLatLon,
                   svLocIncAng,
                   svLocPrjIncAng,
                   svSltBand,
                   svSgNou,
                   stdGrdOriX,
                   stdGrdOriY),
          wd=c(dir_paths[6]),
          echo_cmd = TRUE)
    }


    # Subset
    input <- paste0("-Ssource=",dir_paths[6], "/target.dim")
    geoR <- paste0("-PgeoRegion=POLYGON((",geoRegion,"))")
    sourceband <- paste0("-PsourceBands=", sourceBands)
    cpMtdata <- paste0("-PcopyMetadata=", copyMetadata)
    fSwath <- paste0("-PfullSwath=", fullSwath)
    subSampX <- paste0("-PsubSamplingX=", subSamplingX)
    subSampY <- paste0("-PsubSamplingY=", subSamplingY)

    if(missing(referenceBand)|missing(region)) {
      run(command = "gpt",
          args = c("Subset",
                   input,
                   geoR,
                   sourceband,
                   cpMtdata,
                   fSwath,
                   subSampX,
                   subSampY),
          wd=c(dir_paths[7]),
          echo_cmd = TRUE)
    }else if(missing(geoRegion)) {
      refBand <- paste0("-PreferenceBand=", referenceBand)
      rg <- paste0("-Pregion=", region)
      run(command = "gpt",
          args = c("Subset",
                   input,
                   sourceband,
                   refBand,
                   rg,
                   cpMtdata,
                   fSwath,
                   subSampX,
                   subSampY),
          wd=c(dir_paths[7]),
          echo_cmd = TRUE)

    }


    # Linear To/From dB
    input <- paste0("-Ssource=", dir_paths[7], "/target.dim")
    sourceband <- paste0("-PsourceBands=", sourceBands)

    run(command = "gpt",
        args = c("LinearToFromdB",
                 input,
                 sourceband),
        wd=c(dir_paths[8]),
        echo_cmd = TRUE)


    # Write
    input <- paste0("-Ssource=", dir_paths[8], "/target.dim")
    output_filename <- stri_extract(source[i], regex = '[^T]*')
      # extract the strings before the first occurance of the pattern 'T' in the source file name.
    opfile <- paste0("-Pfile=", output_filename)
    clrCache <- paste0("-PclearCacheAfterRowWrite=", clear_Cache_after_row_write)
    delopFail <- paste0("-PdeleteOutputOnFailure=", delete_output_on_failure)
    wrtTile <- paste0("-PwriteEntireTileRows=", write_entire_tile_rows)

    if(output_GeoTIFF == FALSE){
      run(command = "gpt",
          args = c("Write",
                   input,
                   opfile,
                   "-PformatName=BEAM-DIMAP",
                   clrCache,
                   delopFail,
                   wrtTile),
          wd = c(dir_paths[9]),
          echo_cmd = TRUE)
    }else if(output_BEAM_DIMAP == FALSE) {
      run(command = "gpt",
          args = c("Write",
                   input,
                   opfile,
                   "-PformatName=GeoTIFF",
                   clrCache,
                   delopFail,
                   wrtTile),
          wd = c(dir_paths[9]),
          echo_cmd = TRUE)
    }else{
      run(command = "gpt",
          args = c("Write",
                   input,
                   opfile,
                   "-PformatName=BEAM-DIMAP",
                   clrCache,
                   delopFail,
                   wrtTile),
          wd = c(dir_paths[9]),
          echo_cmd = TRUE)
      run(command = "gpt",
          args = c("Write",
                   input,
                   opfile,
                   "-PformatName=GeoTIFF",
                   clrCache,
                   delopFail,
                   wrtTile),
          wd = c(dir_paths[9]),
          echo_cmd = TRUE)
    }
  }

  # Remove the intermediate results
  unlink(dir_paths[1:8], recursive=TRUE)

}
