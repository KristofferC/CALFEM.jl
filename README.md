# CALFEM.jl

[![Build Status](https://travis-ci.org/KristofferC/CALFEM.jl.svg?branch=master)](https://travis-ci.org/KristofferC/CALFEM.jl)

`CALFEM.jl` is an API port of the simple Matlab FE toolbox [CALFEM](http://www.solid.lth.se/education/courses/finita-elementmetoden-fhlf01-fhl064/software/) written in Julia. The purpose of this package is to ease the transition for people who want to try out Julia for FE-analysis. `CALFEM.jl` is built on top of [JuAFEM](https://github.com/KristofferC/JuAFEM.jl).

Not all of CALFEM is yet implemented. For a list of implemented functions, see [this issue](https://github.com/KristofferC/CALFEM.jl/issues/1). 


## Installation

```jl
Pkg.clone("https://github.com/KristofferC/ContMechTensors.jl")
Pkg.clone("https://github.com/KristofferC/JuAFEM.jl")
Pkg.clone("https://github.com/KristofferC/CALFEM.jl")
```

for the plotting to work you also need to run

```jl
Pkg.add("Winston")
```

## Examples

See the example folder for examples. These are generally written as Jupyter notebooks.

## Documentation:

See the official CALFEM manual which can be found at [here](http://www.solid.lth.se/education/courses/finita-elementmetoden-fhlf01-fhl064/software/).

There is also a briefer online documentation available [here](http://calfemjl.readthedocs.org/en/latest/). This also documents a few differences between this package and the MATLAB package.

