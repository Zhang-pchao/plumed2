"""Small driver parity and error-path checks for the optional CUDA library."""

import argparse
import json
import math
import re
import resource
import subprocess
import tempfile
from pathlib import Path


def read_numbers(path):
    rows = []
    for line in path.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            rows.append([float(value) for value in line.split()])
    return rows


def compare(reference, candidate):
    assert len(reference) == len(candidate)
    maximum = 0.0
    for left, right in zip(reference, candidate):
        assert len(left) == len(right)
        for expected, actual in zip(left, right):
            assert math.isfinite(actual)
            error = abs(expected - actual)
            assert error <= 2e-8 + 2e-10 * abs(expected), (expected, actual, error)
            maximum = max(maximum, error)
    return maximum


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plumed", default="plumed")
    parser.add_argument("--library", required=True)
    parser.add_argument("--work", type=Path)
    parser.add_argument("--cpu-only", action="store_true")
    arguments = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    work = arguments.work or Path(tempfile.mkdtemp(prefix="voronoi-backends-"))
    work.mkdir(parents=True, exist_ok=True)
    positions = [
        [0.2, 0.3, 0.4],
        [1.3, 0.6, 0.2],
        [2.2, 1.7, 0.8],
        [3.0, 2.4, 1.9],
        [0.4, 0.3, 0.4],
        [0.1, 0.5, 0.4],
        [0.2, 0.2, 0.6],
        [1.4, 0.7, 0.3],
        [2.3, 1.8, 0.9],
        [3.1, 2.5, 2.0],
    ]
    common = "CENTERS=1,2,3,4 ASSIGNED=5,6,7,8,9,10 REFERENCE=1.5,1.2,1.1,2.2 KAPPA=5"
    definitions = [
        "VORONOI_COORDINATION " + common + " POWER=1",
        "VORONOI_COORDINATION " + common + " POWER=2 SIGN=POSITIVE",
        "VORONOI_COORDINATION " + common + " POWER=2 SIGN=NEGATIVE SELECT=1,3",
        "VORONOI_COORDINATION " + common + " POWER=2 SELECT=2,4 COEFFICIENTS=1.3,-0.4",
        "VORONOI_COORDINATION CENTERS=1,2,3,4 ASSIGNED=5,6,7,8,9,10 REFERENCE=1 KAPPA=8 POWER=2",
        "VORONOI_DISTANCE " + common + " GROUP1=1,3 GROUP2=2,4",
        "VORONOI_DISTANCE " + common + " GROUP1=1,2,3,4",
        "VORONOI_POSITION " + common + " AXIS=X ORIGIN=0.35",
        "VORONOI_POSITION " + common + " AXIS=Y ORIGIN=-0.7 SIGN=NEGATIVE NORMALIZE",
        "VORONOI_POSITION " + common + " AXIS=Z ORIGIN=1.2 ABSOLUTE NORMALIZE",
        "VORONOI_POSITION " + common + " AXIS=X ORIGIN=-0.4 SELECT=1,4 NORMALIZE",
    ]
    modes = {
        "independent": " NO_SHARE",
        "shared": " BACKEND=CPU",
        "auto-missing": " BACKEND=AUTO CUDA_LIBRARY=missing-backend.so",
    }
    if not arguments.cpu_only:
        modes.update(
            {
                "cuda": " BACKEND=CUDA CUDA_LIBRARY=" + arguments.library,
                "auto": " BACKEND=AUTO CUDA_LIBRARY=" + arguments.library,
            }
        )
    results = {"work": str(work), "positive": [], "negative": []}

    def run(name, lines, coordinates=None, box="7,8,9"):
        directory = work / name
        directory.mkdir()
        coordinates = coordinates or positions
        xyz = "10\nsmall backend regression\n"
        xyz += "\n".join("H " + " ".join(map(str, row)) for row in coordinates) + "\n"
        (directory / "frame.xyz").write_text(xyz)
        (directory / "plumed.dat").write_text("\n".join(lines) + "\n")
        command = [
            arguments.plumed,
            "driver",
            "--ixyz",
            "frame.xyz",
            "--length-units",
            "nm",
            "--box",
            box,
        ]
        with (directory / "stdout.log").open("w") as output:
            process = subprocess.run(
                command, cwd=directory, stdout=output, stderr=subprocess.STDOUT
            )
        return directory, process.returncode

    for geometry in ["pbc", "nopbc", "nlist", "triclinic"]:
        reference = None
        for mode, backend in modes.items():
            if mode == "cuda" and geometry == "triclinic":
                continue
            extra = " NOPBC" if geometry == "nopbc" else ""
            if geometry == "nlist":
                extra = " NLIST NL_CUTOFF=1.85 NL_STRIDE=3 NL_SKIN=0.4"
            labels = ["v" + str(index) for index in range(len(definitions))]
            lines = [
                label + ": " + definition + extra + backend
                for label, definition in zip(labels, definitions)
            ]
            lines += [
                "PRINT ARG=" + ",".join(labels) + " FILE=values FMT=%24.16e",
                "DUMPDERIVATIVES ARG="
                + ",".join(labels)
                + " FILE=derivatives FMT=%24.16e",
            ]
            box = "7,0,0,0.8,8,0,0.4,0.7,9" if geometry == "triclinic" else "7,8,9"
            directory, code = run(geometry + "-" + mode, lines, box=box)
            assert code == 0, directory
            data = {
                key: read_numbers(directory / key) for key in ["values", "derivatives"]
            }
            if reference is None:
                reference = data
            errors = {key: compare(reference[key], data[key]) for key in data}
            log = (directory / "stdout.log").read_text()
            if mode == "auto" and geometry == "triclinic":
                assert "selected CPU" in log
            if mode in ["cuda", "auto"] and geometry != "triclinic":
                assert "selected CUDA" in log
            results["positive"].append(
                dict(geometry=geometry, mode=mode, errors=errors)
            )

    invalid = [
        (
            "bad-kappa",
            "VORONOI_COORDINATION " + common.replace("KAPPA=5", "KAPPA=0"),
            "KAPPA",
        ),
        (
            "bad-reference",
            "VORONOI_COORDINATION "
            + common.replace("REFERENCE=1.5,1.2,1.1,2.2", "REFERENCE=1,2"),
            "REFERENCE",
        ),
        (
            "duplicate",
            "VORONOI_COORDINATION "
            + common.replace("CENTERS=1,2,3,4", "CENTERS=1,1,3,4"),
            "duplicate",
        ),
        (
            "overlap",
            "VORONOI_COORDINATION "
            + common.replace("ASSIGNED=5,6,7,8,9,10", "ASSIGNED=1,6,7,8,9,10"),
            "disjoint",
        ),
        ("bad-power", "VORONOI_COORDINATION " + common + " POWER=3", "POWER"),
        ("bad-select", "VORONOI_COORDINATION " + common + " SELECT=5", "SELECT"),
        (
            "bad-groups",
            "VORONOI_DISTANCE " + common + " GROUP1=1,2 GROUP2=2,3",
            "disjoint",
        ),
        ("bad-axis", "VORONOI_POSITION " + common + " AXIS=W ORIGIN=0", "AXIS"),
        (
            "bad-tolerance",
            "VORONOI_POSITION " + common + " AXIS=X ORIGIN=0 TOLERANCE=0",
            "TOLERANCE",
        ),
        (
            "zero-weight",
            "VORONOI_POSITION CENTERS=1 ASSIGNED=5 REFERENCE=1 KAPPA=5 AXIS=X ORIGIN=0 NORMALIZE",
            "TOLERANCE",
        ),
        (
            "bad-backend",
            "VORONOI_COORDINATION " + common + " BACKEND=UNKNOWN",
            "BACKEND",
        ),
    ]
    for name, definition, expected in invalid:
        for mode, backend in modes.items():
            suffix = "" if name == "bad-backend" else backend
            directory, code = run(
                name + "-" + mode,
                ["v: " + definition + suffix, "PRINT ARG=v FILE=values"],
            )
            assert code != 0, directory
            assert re.search(expected, (directory / "stdout.log").read_text(), re.I), (
                directory
            )
            results["negative"].append(dict(case=name, mode=mode, returncode=code))
    for name, changed_atom, definition in [
        ("coincident-assigned", 4, "VORONOI_COORDINATION " + common),
        ("coincident-centers", 1, "VORONOI_DISTANCE " + common + " GROUP1=1 GROUP2=2"),
    ]:
        coordinates = [row.copy() for row in positions]
        coordinates[changed_atom] = positions[0].copy()
        for mode, backend in modes.items():
            directory, code = run(
                name + "-" + mode,
                ["v: " + definition + backend, "PRINT ARG=v FILE=values"],
                coordinates=coordinates,
            )
            assert code != 0, directory
            assert re.search(
                "zero|separation|forward failed|distance",
                (directory / "stdout.log").read_text(),
                re.I,
            ), directory
            results["negative"].append(dict(case=name, mode=mode, returncode=code))
    if not arguments.cpu_only:
        for name, suffix, box in [
            ("missing-library", " CUDA_LIBRARY=missing-backend.so", "7,8,9"),
            (
                "empty-nlist",
                " CUDA_LIBRARY=" + arguments.library + " NLIST NL_CUTOFF=0.02 NL_STRIDE=3",
                "7,8,9",
            ),
            (
                "strict-triclinic",
                " CUDA_LIBRARY=" + arguments.library,
                "7,0,0,0.8,8,0,0.4,0.7,9",
            ),
        ]:
            lines = [
                "v: VORONOI_COORDINATION " + common + " BACKEND=CUDA" + suffix,
                "PRINT ARG=v FILE=values",
            ]
            directory, code = run(name, lines, box=box)
            assert code != 0, directory
            if name == "empty-nlist":
                assert "no CENTER inside the candidate list" in (
                    directory / "stdout.log"
                ).read_text()
            results["negative"].append(dict(case=name, returncode=code))
    results["status"] = "PASSED"
    (work / "RESULTS.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
