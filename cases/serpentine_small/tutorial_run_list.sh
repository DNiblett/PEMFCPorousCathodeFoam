#!/usr/bin/env bash
#mesh the case if not already
blockMesh
snappyHexMesh -overwrite
topoSet
createPatch -overwrite
setFields
checkMesh

#then run the solver for porous+free-flow combined (uses set permeability tensor and porosity field) 
pemfcPorousFlowFoam

#rename the recently saved converged files to 1
mv "$(ls -dt */ | head -n 1)" 1

# copy concentration and potential files to new timefolder
cp -r 0/C_O2 1/C_O2
cp -r 0/D_O2 1/D_O2
cp -r 0/fi 1/fi

#then change the current drawn, operating mode and run fuel cell solver
pemfcPorousCathodeFoam
