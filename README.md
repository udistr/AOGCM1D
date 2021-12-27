# AOGCM1D.jl

A code for a simple 1D coupled ocean-atmosphere model.

Prognostic variables:

Atmosphere: Heat, horizontal momentum (1 direction) and moisture.
Ocean: Heat and momentum (salinity not yet implemented)
Grid:

height coordinates
Processess included:

* Vertical advaction, explicit integration
* Vetical diffustion:
  * Large and Yeager (1994) for the ocean
  * Holtslag and Boville (1993) for the atmosphere
* Bulk formulea for the air-sea interface based on Large and Yeager (2004)
* Solar radiation including transmision to the deeper ocean
* Condensation and deposition of super-saturated water vapor

## Issues/todo

* holtslag scheme noisy and boundary layer seems too shallow
* sum of ocean SW not equal 1
* improve performance
* restart
* documentation

## Directions

```
git clone https://github.com/udistr/AOGCM1D.jl
julia --project=AOGCM1D.jl -e 'using Pkg; Pkg.instantiate(); include("AOGCM1D.jl/AOGCM1D.jl")'
```

