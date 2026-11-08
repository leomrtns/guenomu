#!/usr/bin/env python3
"""Check mixture output, independent RF constraints, and unsupported modes."""
import csv
import importlib.util
import math
import os
from pathlib import Path
import re
import shutil
import tempfile

spec = importlib.util.spec_from_file_location("output_checks", Path(__file__).with_name("check-output.py"))
checks = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checks)


def splits(line):
    """Parse the generated Newick and independently collect unrooted splits."""
    tokens = re.findall(r"[(),;]|[^\s(),;]+", re.sub(r"\[[^]]*\]", "", line))
    position = 0
    clades = []

    def node():
        nonlocal position
        if tokens[position] != "(":
            label = tokens[position]
            position += 1
            return frozenset([label])
        position += 1
        left = node()
        assert tokens[position] == ","
        position += 1
        right = node()
        assert tokens[position] == ")"
        position += 1
        clade = left | right
        clades.append(clade)
        return clade

    leaves = node()
    assert tokens[position:] == [";"], tokens[position:]
    result = set()
    for clade in clades:
        other = leaves - clade
        if len(clade) > 1 and len(other) > 1:
            result.add(min(tuple(sorted(clade)), tuple(sorted(other))))
    return leaves, result


def configure(directory, trees=2, alpha=0.4):
    checks.fixture(directory)
    text = (directory / "run.ctrl").read_text()
    text = text.replace("10 12", "10 120").replace("param_n_samples = 4", "param_n_samples = 6")
    text = text.replace("end_of_control_file", f"""param_n_species_trees = {trees}
param_species_rf_max = 0.8
param_species_weight_alpha = {alpha}
param_reconciliation_prior = 0.2
end_of_control_file""")
    (directory / "run.ctrl").write_text(text)


def check(directory, trees):
    with (directory / "multiple.params.tsv").open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        header = reader.fieldnames
        rows = list(reader)
    expected = ["iteration", "log_target"]
    for gene in range(6):
        for tree in range(trees):
            expected += [f"weight_{gene}_{tree}", f"support_{gene}_{tree}"]
    assert header == expected, header
    assert [int(row["iteration"]) for row in rows] == [20, 40, 60, 80, 100, 120]
    assert all(math.isfinite(float(value)) for row in rows for value in row.values())
    for row in rows:
        for gene in range(6):
            for kind in ("weight", "support"):
                values = [float(row[f"{kind}_{gene}_{tree}"]) for tree in range(trees)]
                assert all(0 <= value <= 1 for value in values), values
                assert math.isclose(sum(values), 1, abs_tol=1e-10), values
            assert all(float(row[f"weight_{gene}_{tree}"]) > 0 for tree in range(trees))
    states = []
    for tree in range(trees):
        lines = (directory / f"multiple.tree{tree}.trees").read_text().splitlines()
        assert len(lines) == len(rows)
        states.append([splits(line) for line in lines])
    for sample in range(len(rows)):
        for a in range(trees):
            leaves, first = states[a][sample]
            assert len(leaves) == 5
            for b in range(a):
                other, second = states[b][sample]
                assert leaves == other
                rf = len(first ^ second)
                assert 0 < rf <= 0.8 * (2 * len(leaves) - 6)
                assert first & second, "grand trees must share an internal split"
    assert not (directory / "job0.checkpoint.bin").exists(), "experimental mode wrote legacy checkpoint"
    return header


def exercise(binary, root, mpi=None):
    prefix = [mpi, "-n", "2", binary] if mpi else [binary]
    directory = root / ("mpi" if mpi else "serial")
    configure(directory)
    output = checks.run(prefix + ["run.ctrl"], directory)
    assert "generalized-Bayes" in output and "partition functions omitted" in output, output
    header = check(directory, 2)
    if not mpi:
        helper = os.environ.get("GUENOMU_MULTIPLE_TEST")
        if helper:
            checks.run([helper, "run.ctrl"], directory)
        three = root / "three-trees"
        configure(three, 3)
        checks.run([binary, "run.ctrl"], three)
        check(three, 3)
    for command in (["-z", "1", "run.ctrl"], ["-p", "1", "run.ctrl"]):
        output = checks.run(prefix + command, directory, success=False)
        assert "unsupported" in output, output
    text = (directory / "run.ctrl").read_text()
    (directory / "impossible.ctrl").write_text(text.replace("species_rf_max = 0.8", "species_rf_max = 0.1"))
    output = checks.run(prefix + ["impossible.ctrl"], directory, success=False)
    assert "RF prior is empty" in output, output
    (directory / "invalid.ctrl").write_text(text.replace("n_species_trees = 2", "n_species_trees = 9"))
    output = checks.run(prefix + ["invalid.ctrl"], directory, success=False)
    assert "2..8 trees" in output, output
    return header


def main():
    serial = os.environ.get("GUENOMU_SERIAL")
    mpi = os.environ.get("GUENOMU_MPI")
    assert serial or mpi, "GUENOMU_SERIAL or GUENOMU_MPI must identify a built executable"
    with tempfile.TemporaryDirectory(prefix="guenomu-multiple-") as temporary:
        root = Path(temporary)
        serial_header = exercise(serial, root) if serial else None
        if mpi:
            launcher = shutil.which("mpiexec")
            assert launcher, "mpiexec is required for MPI checks"
            mpi_header = exercise(mpi, root, launcher)
            if serial_header:
                assert serial_header == mpi_header
            helper = os.environ.get("GUENOMU_MULTIPLE_TEST")
            if helper:
                checks.run([helper, "run.ctrl"], root / "mpi")
    print("PASS: multi-tree sampling, ordered weights, independent RF checks, math and explicit errors")


if __name__ == "__main__":
    main()
