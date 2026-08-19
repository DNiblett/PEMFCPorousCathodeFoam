#!/usr/bin/env python3
"""Run the Schneider et al. (2010) potentiostatic current-distribution sweep.

All tutorial tortuosities are constrained to be at least one.
"""
from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import matplotlib.pyplot as plt
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


ROOT = Path(__file__).resolve().parents[1]
CASE_DEFAULT = ROOT / 'cases' / 'schneider_validation'
DEFAULT_BASHRC = Path(os.environ.get('WM_PROJECT_DIR', '')) / 'etc' / 'bashrc'
BASHRC = Path(os.environ.get('FOAM_BASHRC', DEFAULT_BASHRC))
# Approximate manual digitisation of the paper's Fig. 3a air trace at 0.1 V.
# It is retained separately so its graphical origin and uncertainty are explicit.
AIR_DIGITISED = {
    0.1: (0.20, 0.11, 0.15, 0.12, 0.25, 0.42, 1.65, 2.50, 2.85, 2.50, 2.65, 2.70, 1.75, 0.55, 0.30, 0.20, 0.13, 0.10, 0.15),
    0.2: (0.13, 0.07, 0.10, 0.06, 0.18, 0.32, 1.45, 2.25, 2.50, 2.20, 2.33, 2.38, 1.62, 0.45, 0.23, 0.12, 0.08, 0.06, 0.08),
    0.3: (0.07, 0.04, 0.06, 0.04, 0.12, 0.25, 1.18, 1.95, 2.12, 1.85, 1.98, 2.04, 1.33, 0.35, 0.18, 0.09, 0.05, 0.03, 0.05),
    0.4: (0.05, 0.03, 0.05, 0.03, 0.08, 0.15, 0.95, 1.55, 1.67, 1.47, 1.57, 1.62, 1.05, 0.27, 0.13, 0.07, 0.04, 0.03, 0.03),
    0.5: (0.03, 0.02, 0.03, 0.02, 0.05, 0.09, 0.60, 1.10, 1.20, 1.05, 1.13, 1.15, 0.75, 0.16, 0.08, 0.04, 0.03, 0.02, 0.02),
    0.6: (0.02, 0.01, 0.02, 0.01, 0.03, 0.05, 0.40, 0.65, 0.72, 0.62, 0.68, 0.70, 0.48, 0.10, 0.05, 0.03, 0.02, 0.01, 0.02),
}
O2_0P6_DIGITISED = (0.95, 0.78, 0.90, 0.83, 1.00, 0.90, 0.95, 0.80, 0.85,
                     0.68, 0.80, 0.84, 1.08, 1.10, 1.15, 1.20, 1.15, 1.10, 0.95)


def shell(command: str, case: Path) -> str:
    prefix = f'source {BASHRC} && ' if BASHRC and BASHRC.is_file() else ''
    result = subprocess.run(['bash', '-lc', prefix + command], cwd=case,
                            text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}): {command}\n"
                           f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}")
    return result.stdout + result.stderr


def numeric_times(case: Path) -> list[Path]:
    return sorted((p for p in case.iterdir() if p.is_dir() and p.name != '0'
                   and re.fullmatch(r'[0-9.eE+-]+', p.name)),
                  key=lambda p: float(p.name))


def set_scalar(path: Path, key: str, value: float) -> None:
    text = path.read_text()
    text, count = re.subn(rf'(?m)^({re.escape(key)}\s+)[^;]+;', rf'\g<1>{value:.8g};', text)
    if count != 1:
        raise RuntimeError(f'Expected one {key} in {path}, found {count}')
    path.write_text(text)


def set_tortuosity(case: Path, value: float, in_plane: float, through_plane: float) -> None:
    if min(value, in_plane, through_plane) < 1.0:
        raise ValueError('Tortuosity must be >= 1')
    path = case / 'constant/transportProperties'
    text = path.read_text()
    text, count = re.subn(r'(tau\s+tau\s+\[[^]]+\]\s+)[0-9.eE+-]+;',
                          rf'\g<1>{value:.8g};', text)
    if count != 1:
        raise RuntimeError(f'Expected one tau entry in {path}, found {count}')
    for key, entry in (('tauInPlane', in_plane), ('tauThroughPlane', through_plane)):
        text, count = re.subn(rf'(?m)^({key}\s+)[^;]+;', rf'\g<1>{entry:.8g};', text)
        if count != 1:
            raise RuntimeError(f'Expected one {key} entry in {path}, found {count}')
    path.write_text(text)


def set_cathode_exchange_current(case: Path, value: float) -> None:
    path = case / 'constant/transportProperties'
    text = path.read_text()
    text, count = re.subn(r'(j_0c\s+)\-?[0-9.eE+-]+;', rf'\g<1>{value:.8g};', text)
    if count != 1:
        raise RuntimeError(f'Expected one j_0c entry in {path}, found {count}')
    path.write_text(text)


def ensure_flow(case: Path, archive: Path, flow_solver: str, flow_source: Path | None = None) -> Path:
    if flow_source is not None:
        if not flow_source.is_dir():
            raise FileNotFoundError(f'Flow source does not exist: {flow_source}')
        return flow_source.resolve()
    existing = archive / 'flow'
    if existing.is_dir():
        return existing
    if not (case / '0').is_dir():
        shutil.copytree(case / '0.template', case / '0')
    if not (case / 'constant/polyMesh').is_dir():
        for command in ('blockMesh', 'snappyHexMesh -overwrite', 'topoSet',
                        'createPatch -overwrite', 'setFields'):
            shell(command, case)
    # Electrochemistry writes time directories without p.  Flow must restart
    # from 0, rather than selecting one of those incomplete latest times.
    for time in numeric_times(case):
        shutil.rmtree(time)
    shell(flow_solver, case)
    flow_time = numeric_times(case)[-1]
    shutil.copytree(flow_time, existing)
    return existing


def patch_flux(case: Path, fields: Path, name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    case, fields = case.resolve(), fields.resolve()
    with tempfile.TemporaryDirectory(prefix='schneider_vtk_') as temporary:
        temp = Path(temporary)
        for directory in ('constant', 'system'):
            (temp / directory).symlink_to(case / directory)
        (temp / '1').symlink_to(fields)
        shell(f'foamToVTK -case {temp} -time 1 -fields \'(Flux)\' -patches \'(catalyst)\' -no-internal', case)
        file = next((temp / 'VTK').glob('*/boundary/catalyst.vtp'))
        reader = vtk.vtkXMLPolyDataReader(); reader.SetFileName(str(file)); reader.Update()
        poly = reader.GetOutput()
        centres = vtk.vtkCellCenters(); centres.SetInputData(poly); centres.Update()
        xyz = vtk_to_numpy(centres.GetOutput().GetPoints().GetData())
        flux = np.abs(vtk_to_numpy(poly.GetCellData().GetArray('Flux')))
        areas = []
        for index in range(poly.GetNumberOfCells()):
            cell = poly.GetCell(index)
            points = np.array([poly.GetPoint(cell.GetPointId(i)) for i in range(cell.GetNumberOfPoints())])
            centre = points.mean(axis=0)
            areas.append(sum(np.linalg.norm(np.cross(points[i] - centre,
                                                     points[(i + 1) % len(points)] - centre)) / 2
                             for i in range(len(points))))
    return xyz[:, 0], flux, np.asarray(areas)


def segment_profile(x: np.ndarray, flux: np.ndarray, area: np.ndarray) -> np.ndarray:
    edges = np.linspace(-0.0038, 0.0038, 20)
    values = []
    for left, right in zip(edges[:-1], edges[1:]):
        use = (x >= left) & ((x < right) if right < edges[-1] else (x <= right))
        values.append(np.average(flux[use], weights=area[use]) / 1e4)
    return np.asarray(values)


def run_point(case: Path, flow: Path, output: Path, gas: str, voltage: float,
              oxygen: float, solver: str) -> tuple[np.ndarray, float]:
    label = f'{gas}_V_{voltage:.1f}'.replace('.', 'p')
    destination = output / label
    for stale in [case / '0', *numeric_times(case)]:
        if stale.exists():
            shutil.rmtree(stale)
    shutil.copytree(flow, case / '0')
    for field in ('C_O2', 'fi', 'porosity', 'inletVelocity', 'inletOxygen',
                  'oxygenDiffusivityMultiplier', 'oxygenDiffusivityThroughPlaneMultiplier', 'kineticMultiplier',
                  'contactResistance'):
        if not (case / '0.template' / field).is_file():
            continue
        shutil.copy2(case / '0.template' / field, case / '0' / field)
    (case / '0/inletOxygen').write_text(f'inletOxygen {oxygen:.8g};\n')
    # The cathode-potential Newton solve is stiff at the low voltages reported
    # by Schneider.  Initialise close to the Butler--Volmer/membrane balance
    # instead of the generic zero field shipped with the galvanostatic example.
    fi = case / '0/fi'
    initial_fi = voltage - 0.85
    fi.write_text(fi.read_text().replace('internalField   uniform 0;',
                                         f'internalField   uniform {initial_fi:.8g};'))
    if (case / '0.template' / 'contactResistance').is_file():
        shell('setFields', case)
    set_scalar(case / 'constant/controlProperties', 'VcellSetpoint', voltage)
    set_scalar(case / 'constant/Results', 'Vcell', voltage)
    log = shell(solver, case)
    (destination / 'logs').mkdir(parents=True, exist_ok=True)
    (destination / 'logs/electrochem.log').write_text(log)
    final = numeric_times(case)[-1]
    shutil.copytree(final, destination / 'final')
    x, flux, area = patch_flux(case, destination / 'final', label)
    profile = segment_profile(x, flux, area)
    match = re.findall(r'Predicted current density\s+=\s+([0-9.eE+-]+)', log)
    if not match:
        raise RuntimeError('Could not find predicted current density in electrochemistry log')
    return profile, abs(float(match[-1])) / 1e4


def plot(output: Path, profiles: dict[tuple[str, float], np.ndarray]) -> None:
    segment = np.arange(1, 20)
    fig, axes = plt.subplots(2, 3, figsize=(12.2, 7.0), sharex=True, sharey=True,
                             constrained_layout=True)
    for axis, voltage in zip(axes.flat, sorted(AIR_DIGITISED)):
        model, digitised = profiles[('air', voltage)], AIR_DIGITISED[voltage]
        rmse = float(np.sqrt(np.mean((model - digitised) ** 2)))
        axis.plot(segment, model, 'o-', ms=3, lw=1.35, color='#1677b8', label='Model')
        axis.plot(segment, digitised, 's--', ms=2.8, lw=1.1, mfc='white', color='#222222',
                  label='Schneider et al. (digitised)')
        axis.axvspan(7.5, 12.5, color='0.90', zorder=-5, label='2 mm channel')
        axis.set(xlim=(.5, 19.5), ylim=(0, 4.55), xticks=range(1, 20, 2), title=f'{voltage:.1f} V')
        axis.text(10, 4.25, f'RMSE={rmse:.2f} A cm$^{{-2}}$', ha='center', fontsize=8)
        axis.grid(alpha=.25)
    for axis in axes[-1, :]: axis.set_xlabel('Segment (400 μm)')
    for axis in axes[:, 0]: axis.set_ylabel('Local current density (A cm$^{-2}$)')
    axes[0, 0].legend(fontsize=7, frameon=False, loc='upper left')
    fig.suptitle('Schneider et al. (2010): H$_2$/air current-distribution comparison', fontsize=13)
    fig.savefig(output / 'schneider_air_current_distributions.png', dpi=260)

    oxygen = plt.figure(figsize=(6.0, 4.2), constrained_layout=True)
    axis = oxygen.gca()
    axis.plot(segment, profiles[('oxygen', .6)], 'o-', ms=3, lw=1.5, label='Model, O$_2$ 0.6 V')
    axis.plot(segment, O2_0P6_DIGITISED, 's--', ms=3, mfc='white', color='#222222',
              label='Schneider et al. (digitised)')
    axis.axvspan(7.5, 12.5, color='0.90', zorder=-5)
    axis.set(xlim=(.5, 19.5), xticks=range(1, 20, 2), xlabel='Segment (400 μm)',
             ylabel='Local current density (A cm$^{-2}$)', title='Oxygen, 0.6 V')
    axis.grid(alpha=.25); axis.legend(fontsize=8, frameon=False)
    oxygen.savefig(output / 'schneider_oxygen_current_distribution.png', dpi=260)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, default=CASE_DEFAULT)
    parser.add_argument('--out', type=Path, default=ROOT / 'results' / 'schneider_validation')
    parser.add_argument('--tau', type=float, default=2.0,
                        help='legacy scalar tortuosity, constrained to >= 1 (default: 2.0)')
    parser.add_argument('--tau-in-plane', type=float, default=1.0,
                        help='V2 in-plane tortuosity, constrained to >= 1 (default: 1.0)')
    parser.add_argument('--tau-through-plane', type=float, default=2.0,
                        help='V2 through-plane tortuosity, constrained to >= 1 (default: 2.0)')
    parser.add_argument('--j0c', type=float, default=-1.5,
                        help='cathode exchange-current parameter in A m^-2 (default: -1.5)')
    parser.add_argument('--solver', default='pemfcPorousCathodeFoam')
    parser.add_argument('--flow-solver', default='pemfcPorousFlowFoam')
    parser.add_argument('--air-voltages', nargs='+', type=float,
                        default=(.1, .2, .3, .4, .5, .6),
                        help='air potentiostatic setpoints matching the digitised figure')
    parser.add_argument('--flow-source', type=Path,
                        help='reuse an archived converged flow-field directory')
    args = parser.parse_args()
    if not args.case.is_dir():
        raise FileNotFoundError(f'Missing case {args.case}; run scripts/make_schneider_case.py first')
    if args.out.exists():
        raise FileExistsError(f'Refusing to overwrite output: {args.out}')
    args.out.mkdir(parents=True)
    set_tortuosity(args.case, args.tau, args.tau_in_plane, args.tau_through_plane)
    set_cathode_exchange_current(args.case, args.j0c)
    flow = ensure_flow(args.case, args.out, args.flow_solver, args.flow_source)
    profiles, rows = {}, []
    for gas, oxygen, voltages in (('air', 5.15, tuple(args.air_voltages)),
                                  ('oxygen', 24.5, (.6,))):
        for voltage in voltages:
            profile, average = run_point(args.case, flow, args.out, gas, voltage, oxygen,
                                         args.solver)
            profiles[(gas, voltage)] = profile
            rows.extend({'gas': gas, 'Vcell_V': voltage, 'segment': index + 1,
                         'local_current_A_cm2': value, 'cell_average_A_cm2': average}
                        for index, value in enumerate(profile))
            print(f'{gas} {voltage:.1f} V: {average:.3f} A cm^-2', flush=True)
    with (args.out / 'simulated_profiles.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys()); writer.writeheader(); writer.writerows(rows)
    with (args.out / 'experimental_digitised_profiles.csv').open('w', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(('gas', 'Vcell_V', 'segment', 'local_current_A_cm2'))
        for voltage, values in AIR_DIGITISED.items():
            writer.writerows(('air', voltage, n, value) for n, value in enumerate(values, 1))
        writer.writerows(('oxygen', .6, n, value) for n, value in enumerate(O2_0P6_DIGITISED, 1))
    plot(args.out, profiles)


if __name__ == '__main__':
    main()
