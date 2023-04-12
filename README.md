# Atlantic reef fish distribution modelling

## Differences in this branch

There's a new working version for including depth (bathymetry) in the model. This involved producing a finer mesh. In this branch we also created a mesh with smaller outer bound and less buffer before simplification of the study area.

For running this new version, run the files `lgcp_prepare_data_Vbath.R` (but this can be skipped, as the uploaded raster layers are already from this file) and `lgcp_model_acch_Vbath.R`.

### Results from this test

- With a finer mesh, I was able to include depth and the model worked pretty well. Estimates for both bathymetry and the other variables were within something normal.  
- Estimated range was much lower than a model with less mesh nodes, whilst the sigma was higher.  
- Predicted lambda was also lower than with less nodes (i.e. with less nodes the value is closer to the actual number of points (see `length(po.pts)`)).  
- Also on this direction, the cross-validation scores for the binomial part of the model (presence-absence) were marginally lower.  
- Removing depth does not alter those findings, thus the reduction on the scores or in the predicted lambda should be because of the higher number of points.  

## Step-by-step modelling

1. Data for each species and environmental layers are **all ready to use**, and were prepared with the files starting with `databases_*`, `env_data_*` and `spdata_*`. You can start from the next step.  
2. Using `lgcp_prepare_data.R` you can generate the **mesh** to be used by `INLA`/`inlabru`. Within this file is also the code to extend the covariates to cover all integration points (because the study area have to be extended a little bit to ensure the barrier model works properly, some points may fall out of the covariates coverage). Running this file will save the objects to be used by the next step.  
3. Each species is modelled separately, so we can chose the best model alternative for each one accordingly. After preparing the mesh AND the extended covariates, you can generate the model with the files starting with `lgcp_model_*`. So, for example, to model the species _Acanthurus chirurgus_ (acronym "acch") you can open the file `lgcp_model_acch.R`. During the modeling, other functions and R source files are loaded. 

***
**IMPORTANT NOTES:**  
To make the barrier model works properly you need to have the `sf` package version 1.0-10 which is not yet on CRAN [this behavior occurred with the most recent updates of the packages used in the modelling]. To install use:

```
library(remotes)
install_github("r-spatial/sf")
```

To speed up computations, change the `intest` argument on the header of the file to `'eb'`.

## Testing alternative meshs

The original files used a mesh without a propper outer bound. We are now testing an appropriate outer bound that uses the `inla.nonconvex.hull()` function to generate an outer bound that follows the study area's shape.

To generate the alternative mesh one should use the **`lgcp_prepare_data_Vnonconv.R`**.

## Bathymetry inclusion

We are also testing including the  **depth (bathymetry)** covariate, which is extremely important but was previously leading the models to crash. With the new mesh (using outer bound and with a smaller max.edge) the model is now converging. It's possible to test the bathymetry inclusion by using the **`lgcp_prepare_data_Vbath.R`** and **`lgcp_model_acch_Vbath.R`**.


## File structure

``` bash
├───codes # Contains all codes
├───data
│   ├───acch # Data for species Acanthurus chirurgus
│   ├───databases # Databases used to extract data
│   ├───env # Environmental layers
│   │   ├───bath_layers # Bathymetry
│   │   ├───crop_layers # Already masked and cropped layers
│   │   ├───fut_layers # Future CMIP6 layers, for prediction
│   │   └───ready_layers # Adjusted layers for inlabru - mesh
│   ├───gis # Shapefiles
│   ├───lujo # Data for species Lutjanus jocu
│   ├───mybo # Data for species Mycteroperca bonaci
│   ├───scze # Data for species Scarus zelindae
│   └───spam # Data for species Sparisoma amplum
├───figures # Final figures
├───functions # Functions used in the modeling
└───results # Files with results for each species
    ├───acch
    │   ├───effects # Covariate effects
    │   └───predictions # Predictions to each scenario
    ├───lujo
    │   ├───effects
    │   └───predictions
    ├───mybo
    │   ├───effects
    │   └───predictions
    ├───scze
    │   ├───effects
    │   └───predictions
    └───spam
        ├───effects
        └───predictions
```

## Packages used

```
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=Portuguese_Brazil.utf8  LC_CTYPE=Portuguese_Brazil.utf8   
[3] LC_MONETARY=Portuguese_Brazil.utf8 LC_NUMERIC=C                      
[5] LC_TIME=Portuguese_Brazil.utf8    

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] fields_14.1       viridis_0.6.2     viridisLite_0.4.1 spam_2.9-1       
 [5] fs_1.6.1          patchwork_1.1.2   ggplot2_3.4.1     raster_3.6-14    
 [9] terra_1.7-3       inlabru_2.7.0     INLA_23.02.04     sp_1.6-0         
[13] foreach_1.5.2     Matrix_1.5-1     

loaded via a namespace (and not attached):
 [1] xfun_0.37          tidyselect_1.2.0   sf_1.0-10          splines_4.2.2     
 [5] lattice_0.20-45    colorspace_2.1-0   vctrs_0.5.2        generics_0.1.3    
 [9] htmltools_0.5.4    utf8_1.2.3         rlang_1.0.6        e1071_1.7-13      
[13] pillar_1.8.1       glue_1.6.2         withr_2.5.0        DBI_1.1.3         
[17] lifecycle_1.0.3    plyr_1.8.8         MatrixModels_0.5-1 dotCall64_1.0-2   
[21] munsell_0.5.0      gtable_0.3.1       evaluate_0.20      codetools_0.2-18  
[25] knitr_1.42         labeling_0.4.2     fastmap_1.1.0      class_7.3-20      
[29] fansi_1.0.4        Rcpp_1.0.10        KernSmooth_2.23-20 scales_1.2.1      
[33] classInt_0.4-8     farver_2.1.1       gridExtra_2.3      digest_0.6.31     
[37] dplyr_1.1.0        grid_4.2.2         rgdal_1.6-4        cli_3.6.0         
[41] tools_4.2.2        magrittr_2.0.3     maps_3.4.1         proxy_0.4-27      
[45] tibble_3.1.8       pkgconfig_2.0.3    rmarkdown_2.20     rstudioapi_0.14   
[49] iterators_1.0.14   R6_2.5.1           units_0.8-1        compiler_4.2.2 
```

