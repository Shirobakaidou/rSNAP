# rSNAP

This package enables users to batch pre-process Sentinel-1 GRD images with SNAP Toolbox in R. For the moment, it contains a single function "preprocessing", which executes a pre-processing workflow of "read input >> subset >> apply orbit file >> remove thermal noise >> radiometric calibration >> speckle filtering >> terrain correction >> write output". 
The pre-condition of using this package is to have SNAP toolbox incl. SNAP GPT installed and add it in system path (e.g. "...\snap\bin\").
Package "processx" and package "stringi" are needed for functions in this package to work. In case the two packages are not automatically installed when installing or loading this package, please install them manually.
