# RClonOr
A version of https://github.com/fmedina7/ClonOr_cpp that is callable from R.

To use the package:

1. Install gsl.

- On Windows, you can download the library, from https://www.stats.ox.ac.uk/pub/Rtools/libs.html, as a part of the file local323.zip. You need to unzip this file to some location on you computer, then set up an Environment Variable called "LIB_GSL" that points to this location (see instructions [here](https://docs.oracle.com/en/database/oracle/machine-learning/oml4r/1.5.1/oread/creating-and-modifying-environment-variables-on-windows.html)).

- On a Mac, you can install gsl using homebrew, macports, or from [source](https://www.gnu.org/software/gsl/).

- On Ubuntu, use the command

```console
sudo apt-get install libgsl-dev
```

2. Udpate R to the latest version.

3. Install the package in R using

```R
devtools::install_github("maugu/RClonOr")
```

4. Test that the package is installed by running a short MCMC chain using

```R
library(RClonOr)
clonal_origin(system.file("extdata", "true_tree_Toy.nwk", package = "RClonOr"),system.file("extdata", "simulatedData_Toy.xmfa", package = "RClonOr"),"test_Toy.xml",a="1,1,1,2,2,1,1,1,0,0,0",x=0,y=5000,z=10,D=50,T=10,R=5,m="0.25,100,1")
```

The arguments to the R function correspond to those used in the command line arguments in the C++ version of ClonalOrigin: the original version can be found [here](https://github.com/xavierdidelot/ClonalOrigin), with the version associated with the paper "Speeding up Inference of Homologous Recombination in Bacteria" being found [here](https://github.com/fjmedinaaguayo/ClonOr_cpp). This package provides an R interface to the latter version, with the intention of making the program more widely accessible.

Sometimes the code crashes, requiring you to close the R session in the process. If this happens, simply try running the same code again.
