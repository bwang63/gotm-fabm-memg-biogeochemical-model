# gotm-fabm-memg-biogeochemical-model

This biogeochemical model is an updated version of the biogeochemical model described in Laurent et al. (2021) and
is coupled to General Ocean Turbulence Model ([GoTM](https://gotm.net/portfolio/software/)) through the Framework
for Aquatic Biogeochemical Models ([FABM](https://github.com/fabm-model/fabm/wiki/Developing-a-new-biogeochemical-model); Bruggeman and Bolding 2014).
The model is updated to include two different sinking schemes, the ballast scheme and the WLin scheme, which aim to
reproduce the increase of remineralization length scale with depth by simulating the protection by minerals,
i.e., CaCO3 and opal, from remineralization and allowing for an increase of sinking velocity, respectively.
Detailed descriptions of these two schemes follow Wang and Fennel (2022).

## Step 1: Download the GoTM code
Go to [GoTM website](https://gotm.net/portfolio/software/), download its source code in your directory ($GOTM_BASE)

## Step 2: Download the biogeochemical model code
Download the biogeochemical model code from this repository, rename the folder to be 'memg', and put it under the directory `$GOTM_BASE/extern/fabm/src/models` 

## Step 3: Making FABM aware of the new biogeochemical model
Open the master configuration file '$GOTM_BASE/extern/fabm/src/CMakeLists.txt' and add our institude into this file. 
```
# List of contributing institutes that should be included in compilation.
# When adding new institutes to the source tree, please add them here as well.
# You can exclude institute directories from compilation by commenting them out.
set(DEFAULT_INSTITUTES
    akvaplan     # Akvaplan-niva, Norway
    au           # University of Aarhus, Denmark
    bb           # Bolding & Bruggeman - formerly Bolding & Burchard
    csiro        # Commonwealth Scientific and Industrial Research Organisation, Australia
    ersem        # European Regional Seas Ecosystem Model
    examples     # Examples supplied with FABM itself
    gotm         # Models ported from original GOTM/BIO library
    iow          # Leibniz Institute for Baltic Sea Research, Germany
    jrc          # EC - Joint Research Centre, Ispra, Italy
    msi          # Marine Systems Institute, Tallinn University of Technology, Estonia
    memg         # Marine Environmental Modelling Group, Dalhousie University, Canada
    niva         # Norwegian Institute for Water Research, Norway
    pclake       # The PCLake model - reference implementation
    pml          # Plymouth Marine Laboratory, United Kingdom
    selma        # Simple EcoLogical Model for the Aquatic - PROGNOS
    su           # Swansea University, United Kingdom
    uhh          # University of Hamburg, Germany
)
```

## Step 4: Configure and compile it
- Create a build directory by the following commands
```
mkdir build
cd build
```

### Configure
- Open `/scripts/linux/gotm_configure.sh`, change the `<GOTM_BASE>` and `<compiler>`.  
```
# if not set use the suggested source code installation directories
GOTM_BASE=${GOTM_BASE:=~/code}

...

# default Fortran compiler is gfortran - overide by setting compuiler like:
# export compiler=ifort
compiler=${compiler:=gfortran}
```

- Execute this script by running the command
```
$GOTM_BASE/scripts/linux/gotm_configure.sh 
```

### Compile
- Open `/scripts/linux/gotm_build.sh`, change the `<compiler>` to make it consistent with `/scripts/linux/gotm_configure.sh`
```
# default Fortran compiler is gfortran - overide by setting compuiler like:
# export compiler=ifort
compiler=${compiler:=gfortran}
```

- Execute this script by running the command
```
$GOTM_BASE/scripts/linux/gotm_build.sh 
```
