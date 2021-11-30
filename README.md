# SOCSEM
Soil Carbonyl Sulfide (OCS) Empirical Model v 8.0

Mary Whelan (mary.whelan@gmail.com)

Department of Environmental Sciences, Rutgers University

Last updated on 2021-11-30

## Summary

This empirical suite of equations is used to generate first-guess estimates 
of soil emission and consumption of atmospheric carbonyl sulfide (OCS). 

This repository presents code that was used to produce the main
results in ME Whelan, Mingjie Shi, Wu Sun, Linda Kooijmans-de Vries, 
Ulli Seibt, and Kadmiel Maseyk "The role of soil OCS fluxes in new estimates of GPP: an empirical model"
submitted for publication. 

## Dependencies

Python code used to produce figures in the manuscript are in
the `src` folder. You need the following prerequisites to run the code.

* Python 3
* numpy

## Variables needed

The approach only requires 3 variables. 

### Soil surface temperature
Soil surface temperature in C = soil_temp 

C units are used to further emphasize our assumption that OCS fluxes from soil are near 0 when the ground is frozen.

Soil temperatures as close to the surface as are available should be used. When available, surface skin temperature is appropriate. 

### Soil moisture
Soil moisture as volumetric water content (VWC) % = soilw 

VWC must be between 0 and 100. Practically, VWC is between 2 and about 50. 

VWC is not required for estimating OCS emission from wetlands. 

### Landcover
Landcover type will determine the equation used to calculate OCS exchange. 

Landcover classification should be translated into these categories 

* grassland
* temperate forest
* boreal forest
* tropical forest
* agricultural
* wetland
* tundra
* desert
* ice

Observed desert soil has near 0 OCS exchange, and ice is assumed to have near 0 OCS exchange. Values of "0" should be assigned for these cover types.

For tundra, there are (currently) no published data. As a first approximation, it may be appropriate to use the temperate forest equations for tundra. When data becomes available, this model will be revised.  

## Output

By convention, emission to the atmosphere is positive, uptake from the atmosphere is negative. 

The model assumes ambient COS at 0.5 ppb.

Output is in picomole OCS per meter squared per second.

## License

This repository is licensed under the [MIT License](./LICENSE). 
