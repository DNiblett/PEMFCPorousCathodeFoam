#!/usr/bin/env bash
# Run the resolved-channel/GDL flow solve, then the PEMFC cathode solve.
# Source the OpenFOAM-v2512 environment before invoking this script.
set -euo pipefail

case_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$case_dir"

for command in blockMesh snappyHexMesh topoSet createPatch setFields checkMesh \
               foamDictionary foamListTimes pemfcPorousFlowFoam pemfcPorousCathodeFoam
do
    command -v "$command" >/dev/null || {
        echo "Required command '$command' is not available. Source OpenFOAM first." >&2
        exit 1
    }
done

# Mesh only once. Keep 0.orig as the chemistry baseline for subsequent runs.
if [[ ! -f constant/polyMesh/points ]]; then
    blockMesh
    snappyHexMesh -overwrite
    topoSet
    createPatch -overwrite
    setFields
    checkMesh
fi

if [[ ! -d 0.orig ]]; then
    cp -a 0 0.orig
fi

# Always start flow from initial fields, not from a previous cathode result.
cp system/controlDict_flow system/controlDict
foamDictionary system/controlDict -entry startFrom -set startTime
foamDictionary system/controlDict -entry startTime -set 0
pemfcPorousFlowFoam

flow_time="$(foamListTimes -latestTime | awk '/^[0-9][0-9.eE+-]*$/ {time=$1} END {print time}')"
if [[ -z "$flow_time" || "$flow_time" == "0" ]]; then
    echo "The flow solver did not write a numerical result time." >&2
    exit 1
fi

# Transfer flow fields, then restore chemistry fields at time 0. Starting
# cathode chemistry at time 0 avoids corrupting uniform/time metadata.
for field in U p phi porosity permeability
do
    [[ -f "$flow_time/$field" ]] || {
        echo "Missing flow field '$flow_time/$field'." >&2
        exit 1
    }
    cp "$flow_time/$field" "0/$field"
done
cp 0.orig/C_O2 0/C_O2
cp 0.orig/fi 0/fi

# Select the PEMFC control file explicitly. The legacy
# controlDict_electrochem targets electrolyisConcFoam3.
cp system/controlDict_tutorial_electrochem system/controlDict
foamDictionary system/controlDict -entry startFrom -set startTime
foamDictionary system/controlDict -entry startTime -set latestTime
pemfcPorousCathodeFoam
