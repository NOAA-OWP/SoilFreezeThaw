export cmake=/usr/local/Cellar/cmake/3.21.3_1/bin//cmake
export python=/Users/ahmadjan/opt/anaconda3/bin/python
export wkdir=$(pwd)
export builddir="cmake_build"

#ngen=https://github.com/noaa-owp/ngen
#sft=https://github.com/NOAA-OWP/SoilFreezeThaw
smp=https://github.com/NOAA-OWP/SoilMoistureProfiles

# clone repos
#git clone $ngen
#git clone $sft "${wkdir}/ngen/extern/SoilFreezeThaw/SoilFreezeThaw"
#git clone $smp "${wkdir}/extern/SoilMoistureProfiles/SoilMoistureProfiles"

#### Init submodules ===========================================
git submodule update --init --recursive


if [ -d ${builddir} ]; then 
  rm -rf ${builddir}
fi

# Buidling
for model in iso_c_fortran_bmi cfe evapotranspiration
do
   echo "MODEL: $model"
   rm -rf extern/$model/${builddir}
   if [ "$model" = "t-route" ]; then
      cd extern/$model/src/python_routing_v02
      export LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBRARY_PATH #to fix the -lgfortran not found.
      FC=gfortran ./compiler.sh
      cd ../../../..
   else
      if [ "$model" = "cfe" ]; then
         cd extern/${model}/$model
         git checkout cfe_soilfreezethaw 
         cd ../../..
      fi

      if [ "$model" = "iso_c_fortran_bmi" ]; then
	  cmake -B extern/${model}/${builddir} -S extern/${model}
	  make -j 4 -C extern/${model}/${builddir}
      fi
   fi
   
   if [ "$model" = "evapotranspiration" ]; then
       cmake -B extern/${model}/${builddir} -S extern/${model}/${model}
       make -j 4 -C extern/${model}/${builddir}
   fi
   
done

#for model in evapotranspiration SoilFreezeThaw SoilMoistureProfiles
for model in SoilFreezeThaw SoilMoistureProfiles
do
    echo "SFT: extern/$model/${builddir}"
    rm -rf extern/$model/${builddir}
    cmake -B extern/${model}/${buildDir} -S extern/${model}/${model} -DNGEN:BOOL=ON
    make -j 4 -C extern/${model}/${buildDir}
done

#### Build ngen and unit tests ============================================
#cmake -B $buildDir -S . \
#       -DNGEN_ACTIVATE_PYTHON:BOOL=ON \
#       -DNGEN_ACTIVATE_ROUTING:BOOL=ON \
#       -DBMI_C_LIB_ACTIVE:BOOL=ON \
#       -DBMI_FORTRAN_ACTIVE:BOOL=ON 
#make -j 4 -C $builddir 

