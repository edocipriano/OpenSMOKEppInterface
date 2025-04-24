# OpenSMOKEppInterface
C and Python interface for the [OpenSMOKE++ library](https://doi.org/10.1016/j.cpc.2015.02.014)

## 📦 Requirements
The library relies on the following external dependencies:
* **OpenSMOKE++** (_not provided with this repository_)
* [**Boost**](https://www.boost.org)
* [**Eigen**](https://eigen.tuxfamily.org/index.php?title=Main_Page)
* [**Config++**](https://hyperrealm.github.io/libconfig/)

where OpenSMOKE++ and Eigen are headers-only, while the other libraries must be compiled or installed using your favorite package manager.

### 🚀 BLAS Support for performance
Additionally, for speed-up of mathematical operations, one of the following libraries must be installed:
* [**OpenBLAS**](http://www.openmathlib.org/OpenBLAS/)
* [**Intel MKL**](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)

## ⚙️ Environment Setup
Create the following environment variables which store the path to the OpenSMOKE++ library, and to this OpenSMOKEppInterface repository.
As an example:
```
export OPENSMOKE_LATEST=${HOME}/OpenSMOKE/OpenSMOKEppMinimal4Basilisk/source
```
```
export OPENSMOKE_INTERFACE=${HOME}/OpenSMOKE/OpenSMOKEppInterface
```

## 🛠️ Installation
The default tool for installing this library is [cmake](https://cmake.org), which attempts to automatically find the dependency libraries automatically, if they are installed in standard paths.
Assuming the use of OpenBLAS support, the compilation is carried out with the following commands.
```
mkdir build
cd build
cmake -DOPENSMOKE_USE_OPENBLAS=1 ..
make
```

If any error occurs during compilation, you need to make sure that cmake prompts the correct path. An example of the cmake output is:
```
...
-- Found Boost: /opt/homebrew/lib/cmake/Boost-1.85.0/BoostConfig.cmake (found version "1.85.0") found components: date_time filesystem program_options system regex timer chrono 
-- Boost_INCLUDE_DIRS    = /opt/homebrew/include
-- Boost_LIBRARIES       = Boost::date_time;Boost::filesystem;Boost::program_options;Boost::system;Boost::regex;Boost::timer;Boost::chrono
-- Found Eigen3: /opt/homebrew/include/eigen3 (Required is at least version "2.91.0") 
-- Eigen_LIBRARIES       = /opt/homebrew/include/eigen3
-- OpenSmoke++           = /Users/ecipriano/OpenSMOKE/OpenSMOKEppMinimal4Basilisk/source
-- CONFIG++_INCLUDE_DIR  = /Users/ecipriano/Local/libconfig/build/include
-- CONFIG++_LIBRARY      = /Users/ecipriano/Local/libconfig/build/lib/libconfig++.dylib
-- OpenBLAS_INCLUDE_DIRS = /opt/homebrew/Cellar/openblas/0.3.27/include
-- BLAS_LIBRARIES        = /opt/homebrew/Cellar/openblas/0.3.27/lib/libopenblas.dylib
...
```

If one of the libraries cannot be found automatically by cmake, i.e. if you installed it in non-standard paths, you can provide the path to the library in the `cmake` command. Considering for example the OpenBLAS library, the sytax is:
```
cmake -DOPENSMOKE_USE_OPENBLAS=1 -DOpenBLAS_ROOT=/opt/homebrew/Cellar/openblas/0.3.27 ..
```
