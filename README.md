# SOCSEM
Soil Carbonyl Sulfide (OCS) Empirical Model v 8.0

This set of equations is a series of empirical fits
of soil temperature, soil moisture, and land use
to ground surface fluxes of the gas carbonyl sulfide (OCS) 
based on sparse soil incubation and field observation.

Variables needed: 

soil surface temperature in C = soil_temp 

soil moisture as VWC % = soilw 

Landcover Type [grassland, temperate forest, boreal forest, tropical forest, agricultural, wetland]. Landcover type will determine the equation used. 

This assumes ambient COS at 0.5 ppb, deserts have near 0 OCS exchange, and ice has near 0 OCS exchange. For tundra, there are (currently) no published data. When data becomes available, this model will be revised.

By convention, emission to the atmosphere is positive,
uptake from the atmosphere is negative. Output is in picomole per meter squared per second.
