# Newspaper for the EpiLPS package #

### Version 0.1.5 ### (**2021-10-31**)

* Release of first version on Github. Version name: "EpiLPS is born".

### Version 1.0.1 ### (**2021-12-08**)

* Version name: "EpiLPS smoothing"
* The *epilps()* routine can now be fitted to count data containing NAs.
* Possibility to add dates to the *epilps()* routine.
* Version 1.0.1 first submitted to CRAN on (2021-12-17).

### Version 1.0.2 ### (**2021-12-29**)

* Version name: "EpiLPS smoothing".
* Patched cubicBspline.cpp.

### Version 1.0.5 ### (**2022-05-05**)

* Version name: "Schrödinger's cat jumps out of the box".
* Added correction factor on LPSMAP covariance matrix.
* episim() can now simulate from negative binomial distribution.

### Version 1.0.6 ### (**2022-08-08**)

* Version name: "Welcome to the Metaverse".
* Added midwindow option to the perfcheck routine. It reports EpiEstim 
  reproduction number (Rt) estimates at the middle of the sliding window.
  
### Version 1.0.7 ### (**2023-01-17**)

* Version name: "Thunderlight".
* Added epicurve() routine to plot the epidemic curve based on incidence data.
* Added minor changes to plot outputs.

### Version 1.0.8 ### (**2023-02-12**)

* Version name: "EpiLPS Kernels".
* Added a kernel structure for a new architecture.

### Version 1.1.0 ### (**2023-04-04**)

* Version name: "Reducing cyclomatic complexity".
* Simplification of routine inputs.
* Separation of old *epilps()* routine into two smaller routines, namely
  *estimR()* and *estimRmcmc()*.
* Added and simplified S3 methods.

### Version 1.2.0 ### (**2023-08-30**)

* Version name: "Extension to incubation estimation".
* Added priors to KerMCMC routine.
* Added routine for histogram smoothing with LPS.
* Added routines to estimate the incubation period distribution.

### Version 1.3.0 ### (**2024-03-07**)

* Version name: "Extension to nowcasting".
* Added nowcasting routine and associated S3 method.
* Added nowcastingR routine and associated S3 method.
* Added incidence and mortality datasets for Belgium.
* Corrected Rsd formula in KerRpostmap.cpp file.
* Added 0 as a lower bound on ylim for S3_Rt_plot.R.
* Updated KerMCMC.R to compute Dlogtar(zeta_prop, lambda_cur) only once in MALA.
* Updated CITATION file.

### Version 1.4.0 ### (**2024-10-24**)

* Version name: "Extension to serial interval estimation".
* Updated nowcasting.R routine output.
* Added estimSI.R and Kerserialint.cpp for serial interval estimation.
* Added dataset on influenza serial interval data.



