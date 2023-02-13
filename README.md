# Atlantic reef fish distribution modelling

## Step-by-step modelling

1. Data for each species and environmental layers are all ready to use, and were prepared with the files starting with `databases_*`, `env_data_*` and `spdata_*`.  
2. Using `lgcp_prepare_data.R` you can generate the **mesh** to be used by `INLA`/`inlabru`. Within this file is also the code to extend the covariates to cover all integration points (because the study area have to be extended a little bit to ensure the barrier model works properly, some points may fall out of the covariates coverage).  
3. Each species is modelled separately, so we can chose the best model alternative for each one accordingly. After preparing the mesh AND the extended covariates, you can generate the model with the files starting with `lgcp_model_*`. So, for example, to model the species _Acanthurus chirurgus_ (acronym "acch") you can open the file `lgcp_model_acch.R`. During the modeling, other functions and R source files are loaded. 

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
