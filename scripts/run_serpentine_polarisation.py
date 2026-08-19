#!/usr/bin/env python3
"""Mesh the tutorial serpentine field, run a galvanostatic V2 sweep, and plot it."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import re
import shutil
import subprocess

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
CASE_DEFAULT = ROOT / 'cases' / 'serpentine'
FLOW_SOLVER = 'pemfcPorousFlowFoam'
CATHODE_SOLVER = 'pemfcPorousCathodeFoam'
FLOW_CONTROL_DICT = 'controlDict_flow'
ELECTROCHEM_CONTROL_DICT = 'controlDict_tutorial_electrochem'


def run(command: list[str], case: Path, log: Path | None = None) -> str:
    result = subprocess.run(command, cwd=case, text=True, capture_output=True)
    text = result.stdout + result.stderr
    if log is not None:
        log.write_text(text)
    if result.returncode:
        raise RuntimeError(f"{' '.join(command)} failed:\n{text}")
    return text


def latest_time(case: Path) -> str:
    result = run(['foamListTimes', '-latestTime'], case)
    matches = re.findall(r'(?m)^([0-9][0-9.eE+-]*)\s*$', result)
    value = matches[-1] if matches else ''
    if not value or value == '0':
        raise RuntimeError('No numerical time was written')
    return value


def replace_entry(path: Path, key: str, value: float) -> None:
    text = path.read_text()
    text, count = re.subn(rf'(?m)^({re.escape(key)}\s+)[^;]+;', rf'\g<1>{value:.8g};', text)
    if count != 1:
        raise RuntimeError(f'Expected one {key} entry in {path}, found {count}')
    path.write_text(text)


def value(path: Path, key: str) -> float:
    match = re.search(rf'(?m)^{re.escape(key)}\s+([0-9.eE+-]+);', path.read_text())
    if not match:
        raise RuntimeError(f'Missing {key} in {path}')
    return float(match.group(1))


def set_inlet_fields(case: Path, velocity: float | None = None) -> None:
    transport = case / 'constant/transportProperties'
    oxygen = value(transport, 'inletOxygen')
    velocity = value(transport, 'inletVelocity') if velocity is None else velocity
    oxygen_field = case / '0/C_O2'
    text = oxygen_field.read_text()
    text, internal = re.subn(r'(?m)^(internalField\s+uniform\s+)[^;]+;',
                             rf'\g<1>{oxygen:.8g};', text)
    text, inlet = re.subn(r'(channelInlet[\s\S]*?value\s+uniform\s+)[^;]+;',
                          rf'\g<1>{oxygen:.8g};', text)
    if internal != 1 or inlet != 1:
        raise RuntimeError(f'Could not set C_O2 inlet values in {oxygen_field}')
    oxygen_field.write_text(text)
    velocity_field = case / '0/U'
    text = velocity_field.read_text()
    text, count = re.subn(r'(channelInlet[\s\S]*?refValue\s+uniform\s+)-?[0-9.eE+-]+;',
                          rf'\g<1>-{velocity:.8g};', text)
    if count != 1:
        raise RuntimeError(f'Could not set U inlet value in {velocity_field}')
    velocity_field.write_text(text)


def reported_voltage(log_text: str) -> float:
    """Read the converged terminal voltage reported by the cathode solver.

    The legacy Results dictionary is modified during an OpenFOAM time write;
    parsing the solver's convergence report avoids reading a stale value.
    """
    matches = re.findall(r'(?m)^Cell Voltage =\s*([0-9.eE+-]+)\s+V', log_text)
    if not matches:
        raise RuntimeError('Cathode solver did not report a cell voltage')
    return float(matches[-1])


def prepare(case: Path) -> None:
    """Create the mesh and a fresh field baseline on first use."""
    if (case / 'constant/polyMesh').is_dir() and (case / '0.orig').is_dir():
        return
    shutil.copytree(case / '0.template', case / '0', dirs_exist_ok=True)
    set_inlet_fields(case)
    for command in (['blockMesh'], ['snappyHexMesh', '-overwrite'], ['topoSet'],
                    ['createPatch', '-overwrite'], ['setFields'], ['checkMesh']):
        run(command, case, case / f'log.{command[0]}')
    shutil.copytree(case / '0', case / '0.orig', dirs_exist_ok=True)


def plot(rows: list[dict[str, float]], path: Path) -> None:
    current = [row['current_density_A_m2'] / 1e4 for row in rows]
    voltage = [row['Vcell_V'] for row in rows]
    fig, axis = plt.subplots(figsize=(6.3, 4.4), constrained_layout=True)
    axis.plot(current, voltage, 'o-', color='#1677b8', lw=1.8, ms=5)
    axis.set(xlabel='Current density (A cm$^{-2}$)', ylabel='Cell voltage (V)',
             title='1 cm$^2$ serpentine field: V2 polarisation curve')
    axis.grid(alpha=.28)
    fig.savefig(path, dpi=260)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, default=CASE_DEFAULT)
    parser.add_argument('--currents', nargs='+', type=float,
                        default=(1000, 3000, 5000, 7000, 9000, 11000),
                        help='galvanostatic values in A m^-2')
    parser.add_argument('--force-remesh', action='store_true',
                        help='require a clean cloned case before remeshing')
    args = parser.parse_args()
    case = args.case.resolve()
    if args.force_remesh and (case / 'constant/polyMesh').exists():
        raise RuntimeError('Refusing to delete a mesh; clone the case directory before --force-remesh')
    for application in (FLOW_SOLVER, CATHODE_SOLVER, 'foamListTimes'):
        if shutil.which(application) is None:
            raise RuntimeError(f'Missing {application}; source OpenFOAM and run scripts/build_solvers.sh')
    prepare(case)
    area = value(case / 'system/OFdata.txt', 'cellCrossArea')
    temperature = value(case / 'system/OFdata.txt', 'Tgas')
    pressure = value(case / 'system/OFdata.txt', 'pGas')
    oxygen_fraction = value(case / 'system/OFdata.txt', 'yO2')
    stoich = value(case / 'system/OFdata.txt', 'stoichAir')
    inlet_area = value(case / 'system/OFdata.txt', 'inletArea')
    output = case / 'tutorial_results'; output.mkdir(exist_ok=True)
    if any(output.iterdir()):
        raise FileExistsError(f'Refusing to overwrite existing results: {output}')
    rows: list[dict[str, float]] = []
    for current_density in args.currents:
        current = current_density * area
        velocity = ((stoich * current / (4 * 96485.33212 * oxygen_fraction))
                    * (8.314462618 * temperature / pressure) / inlet_area)
        shutil.copytree(case / '0.orig', case / '0', dirs_exist_ok=True)
        (case / 'constant/appliedCurrent').write_text(f'inputValue {current:.12g};\n')
        set_inlet_fields(case, velocity)
        label = f'J_{current_density:g}'
        shutil.copy2(case / 'system' / FLOW_CONTROL_DICT, case / 'system/controlDict')
        run([FLOW_SOLVER], case, output / f'{label}_flow.log')
        flow_time = latest_time(case)
        shutil.move(str(case / '0'), str(output / f'{label}_initial'))
        shutil.move(str(case / flow_time), str(case / '0'))
        for field in ('C_O2', 'fi'):
            source = case / '0.orig' / field
            if source.is_file(): shutil.copy2(source, case / '0' / field)
        shutil.copy2(case / 'system' / ELECTROCHEM_CONTROL_DICT, case / 'system/controlDict')
        electrochem_log = run([CATHODE_SOLVER], case, output / f'{label}_electrochem.log')
        voltage = reported_voltage(electrochem_log)
        electrochem_time = latest_time(case)
        shutil.move(str(case / electrochem_time), str(output / f'{label}_final'))
        rows.append({'current_density_A_m2': current_density, 'current_A': current,
                     'inlet_velocity_m_s': velocity, 'Vcell_V': voltage})
        if voltage <= 0.0: break
    with (output / 'polarisation.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)
    plot(rows, output / 'polarisation_curve.png')
    print(f'Wrote {output / "polarisation.csv"}')


if __name__ == '__main__':
    main()
