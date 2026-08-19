#!/usr/bin/env python3
"""Generate the Schneider et al. (2010) single-channel validation case.

Run under Blender: the script copies the bundled fractal case input files,
replaces the geometry with a 7.6 mm wide GDL plus a 2 mm x 1 mm channel, and
exports a watertight fluid STL for snappyHexMesh.
"""
from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

import bpy


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / 'cases' / '_base_template'
DEFAULT = ROOT / 'cases' / 'schneider_validation'
GDL_H = 0.00028              # Toray TGP-H-090 nominal thickness
CHANNEL_W, CHANNEL_H = 0.002, 0.001
ACTIVE_W, ACTIVE_L = 0.0076, 0.0080  # 19 segments x 0.4 mm, stated 0.75 x 0.8 cm field


BLOCK_MESH = r'''FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }
convertToMeters 1;
vertices
(
    (-0.0038 -0.0040 0) (0.0038 -0.0040 0) (0.0038 0.0040 0) (-0.0038 0.0040 0)
    (-0.0038 -0.0040 0.00128) (0.0038 -0.0040 0.00128) (0.0038 0.0040 0.00128) (-0.0038 0.0040 0.00128)
);
blocks (hex (0 1 2 3 4 5 6 7) (76 80 20) simpleGrading (1 1 1));
edges ();
boundary
(
    inlet { type patch; faces ((0 1 5 4)); }
    outlet { type patch; faces ((3 7 6 2)); }
    wallWest { type wall; faces ((0 4 7 3)); }
    wallEast { type wall; faces ((1 2 6 5)); }
    catalyst { type wall; faces ((0 3 2 1)); }
    top { type wall; faces ((4 5 6 7)); }
);
'''

SNAPPY = r'''FoamFile { version 2.0; format ascii; class dictionary; object snappyHexMeshDict; }
castellatedMesh true;
snap false;
addLayers false;
geometry { fluid { type triSurfaceMesh; file "fluid.stl"; } }
castellatedMeshControls
{
    maxLocalCells 1000000;
    maxGlobalCells 2000000;
    minRefinementCells 5;
    nCellsBetweenLevels 1;
    features ();
    refinementSurfaces { fluid { level (1 1); patchInfo { type wall; } } }
    resolveFeatureAngle 30;
    refinementRegions {};
    locationInMesh (0 -0.0035 0.00078);
    allowFreeStandingZoneFaces false;
}
snapControls { nSmoothPatch 3; tolerance 1.5; nSolveIter 30; nRelaxIter 5; }
addLayersControls { relativeSizes true; layers {}; expansionRatio 1.0; finalLayerThickness 0.3; minThickness 0.1; nGrow 0; featureAngle 30; nRelaxIter 3; nSmoothSurfaceNormals 1; nSmoothNormals 3; nSmoothThickness 10; maxFaceThicknessRatio 0.5; maxThicknessToMedialRatio 0.3; minMedialAxisAngle 90; nBufferCellsNoExtrude 0; nLayerIter 50; }
meshQualityControls { maxNonOrtho 65; maxBoundarySkewness 20; maxInternalSkewness 4; maxConcave 80; minVol 1e-13; minTetQuality 1e-15; minArea -1; minTwist 0.02; minDeterminant 0.001; minFaceWeight 0.05; minVolRatio 0.01; minTriangleTwist -1; nSmoothScale 4; errorReduction 0.75; }
writeFlags (scalarLevels);
debug 0;
mergeTolerance 1e-6;
'''

TOPO_SET = r'''FoamFile { version 2.0; format ascii; class dictionary; object topoSetDict; }
actions
(
    { name channelInletFaces; type faceSet; action new; source patchToFace; patch inlet; }
    { name channelInletFaces; type faceSet; action subset; source boxToFace; box (-0.00101 -0.00401 0.00027) (0.00101 -0.00399 0.00129); }
    { name channelOutletFaces; type faceSet; action new; source patchToFace; patch outlet; }
    { name channelOutletFaces; type faceSet; action subset; source boxToFace; box (-0.00101 0.00399 0.00027) (0.00101 0.00401 0.00129); }
    { name catalystFaces; type faceSet; action new; source patchToFace; patch fluid; }
    { name catalystFaces; type faceSet; action subtract; source normalToFace; normal (1 0 0); cos 0.01; }
    { name catalystFaces; type faceSet; action subtract; source normalToFace; normal (-1 0 0); cos 0.01; }
    { name catalystFaces; type faceSet; action subtract; source normalToFace; normal (0 1 0); cos 0.01; }
    { name catalystFaces; type faceSet; action subtract; source normalToFace; normal (0 -1 0); cos 0.01; }
    { name catalystFaces; type faceSet; action subtract; source normalToFace; normal (0 0 1); cos 0.01; }
);
'''

CREATE_PATCH = r'''FoamFile { version 2.0; format ascii; class dictionary; object createPatchDict; }
pointSync false;
patches
(
    { name channelInlet; patchInfo { type patch; } constructFrom set; set channelInletFaces; }
    { name channelOutlet; patchInfo { type patch; } constructFrom set; set channelOutletFaces; }
    { name catalyst; patchInfo { type wall; } constructFrom set; set catalystFaces; }
);
'''

SET_FIELDS = r'''FoamFile { version 2.0; format ascii; class dictionary; object setFieldsDict; }
defaultFieldValues ( volScalarFieldValue porosity 0.999999 );
regions
(
    boxToCell
    {
        type box;
        box (-0.00379 -0.00399 0) (0.00379 0.00399 0.00028);
        fieldValues ( volScalarFieldValue porosity 0.78 );
    }
);
'''

OF_DATA = '''cellCrossArea 6.08e-5;
Tgas 343;
pGas 101325;
yO2 0.21;
stoichAir 20;
inletArea 2.0e-6;
'''


def make_box(name: str, centre: tuple[float, float, float], size: tuple[float, float, float]):
    bpy.ops.mesh.primitive_cube_add(location=centre)
    obj = bpy.context.object
    obj.name = name
    obj.dimensions = size
    bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
    return obj


def create_stl(case: Path) -> None:
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)
    gdl = make_box('gdl', (0, 0, GDL_H / 2), (ACTIVE_W, ACTIVE_L, GDL_H))
    channel = make_box('channel', (0, 0, GDL_H + CHANNEL_H / 2),
                       (CHANNEL_W, ACTIVE_L, CHANNEL_H))
    bpy.context.view_layer.objects.active = gdl
    gdl.select_set(True)
    modifier = gdl.modifiers.new('watertight_union', 'BOOLEAN')
    modifier.operation, modifier.solver, modifier.object = 'UNION', 'EXACT', channel
    bpy.ops.object.modifier_apply(modifier=modifier.name)
    bpy.data.objects.remove(channel, do_unlink=True)
    destination = case / 'constant/triSurface/fluid.stl'
    bpy.ops.object.select_all(action='DESELECT')
    gdl.select_set(True)
    bpy.context.view_layer.objects.active = gdl
    bpy.ops.wm.stl_export(filepath=str(destination), export_selected_objects=True,
                          ascii_format=True)


def configure(case: Path) -> None:
    (case / 'system/blockMeshDict').write_text(BLOCK_MESH)
    (case / 'system/snappyHexMeshDict').write_text(SNAPPY)
    (case / 'system/setFieldsDict').write_text(SET_FIELDS)
    (case / 'system/topoSetDict').write_text(TOPO_SET)
    (case / 'system/createPatchDict').write_text(CREATE_PATCH)
    (case / 'system/OFdata.txt').write_text(OF_DATA)
    (case / 'constant/OFdata.txt').write_text(OF_DATA)
    control = case / 'constant/controlProperties'
    control_text = control.read_text()
    settings = {
        'operatingMode': 'potentiostatic',
        'VcellSetpoint': '0.6',
        'initialCurrent': '0.0608',
        'potentiostaticCurrentTolerance': '1e-7',
    }
    for key, value in settings.items():
        pattern = rf'(?m)^({key}\s+)[^;]+;'
        if re.search(pattern, control_text):
            control_text = re.sub(pattern, rf'\g<1>{value};', control_text)
        else:
            control_text += f'\n{key} {value};'
    control.write_text(control_text + '\n')
    transport = case / 'constant/transportProperties'
    text = transport.read_text().replace('T \t\t330;', 'T \t\t343;').replace('C_0 \t\t7;', 'C_0 \t\t5.15;')
    transport.write_text(text)
    results = case / 'constant/Results'
    results.write_text(results.read_text().replace('Vcell           0.816683;', 'Vcell           0.6;').replace('I               0.0960399;', 'I               0.0608;'))
    for directory in (case / '0.template',):
        for field in ('U', 'p', 'C_O2'):
            path = directory / field
            text = path.read_text().replace('"inlet"', '"channelInlet"')
            text = text.replace('"outlet"', '"channelOutlet"')
            path.write_text(text)


def main() -> None:
    args = sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else []
    case = Path(args[0]).resolve() if args else DEFAULT
    if case.exists():
        raise FileExistsError(f'Refusing to overwrite existing case: {case}')
    shutil.copytree(TEMPLATE, case)
    create_stl(case)
    configure(case)
    print(case)


if __name__ == '__main__':
    main()
