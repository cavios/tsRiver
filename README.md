# tsRiver

Reconstruction of river water level time series based om multi-mission satellite altimetry 

The model is build in the R-package [`TMB`](http://www.tmb-project.org)

The source code consists of two parts:

```riverSpline.R``` and ```riverSpline.cpp```

Besides at test data set ```Missouri_River_sioux.dat``` is available.

The code 

```R
source('riverSpline.R')
```

The code depends on the following R packages:

```R
TMB

```
please install these before running the code

