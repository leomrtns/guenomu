#!/usr/bin/env python3
"""Exercise checkpoint summaries, original gene IDs, and MPI collective errors."""
import argparse
import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def run(command, directory, success=True):
    result = subprocess.run(command, cwd=directory, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=45,
                            env=dict(os.environ, OMP_NUM_THREADS="1"))
    if (result.returncode == 0) != success:
        raise AssertionError(f"Unexpected exit {result.returncode}: {command}\n{result.stdout}")
    return result.stdout


def fixture(directory, genes=6):
    directory.mkdir()
    trees = ["((A_0,B_0),(C_0,(D_0,E_0)))",
             "(((A_0,A_1),B_0),(C_0,(D_0,E_0)))",
             "((A_0,B_0),(C_0,D_0))",
             "(((A_0,A_1),(B_0,B_1)),(C_0,(D_0,E_0)))",
             "((A_0,C_0),(B_0,E_0))",
             "((A_0,C_0),(B_0,(D_0,E_0)))"]
    for i in range(genes):
        tree = trees[i]
        for species in "ABCDE":
            tree = tree.replace(f"{species}_", f"Sp{species}_")
        (directory / f"g{i}.tre").write_text(f"#NEXUS\nBegin trees;\ntree t = {tree};\nEnd;\n")
    (directory / "genes.txt").write_text("".join(f"g{i}.tre\n" for i in range(genes)))
    (directory / "species.txt").write_text("SpA\nSpB\nSpC\nSpD\nSpE\n")
    (directory / "run.ctrl").write_text("""param_execute_action = importance
param_file_with_tree_files = genes.txt
param_file_with_species_names = species.txt
param_n_generations = 10 12
param_n_samples = 4
param_n_output = 2
param_n_mc_samples = 0
param_anneal = 1 2 0.5 1
param_use_distances = 1111000
param_random_number_seed = 12345
end_of_control_file
""")


def table(path):
    with path.open() as stream:
        rows = list(csv.reader(stream, delimiter="\t"))
    header = rows[0]
    assert len(header) == len(set(header)), (path, "duplicate columns")
    assert len(rows) > 1, (path, "no samples")
    assert all(len(row) == len(header) for row in rows[1:]), (path, "ragged columns")
    values = [dict(zip(header, map(float, row))) for row in rows[1:]]
    assert all(math.isfinite(value) for row in values for value in row.values())
    return header, values


def check(directory, genes, jobs):
    header, values = table(directory / "params.txt")
    gene_names = [name for name in header if name.startswith("lnl_")]
    assert gene_names == [f"lnl_{i}" for i in range(genes)], gene_names
    assert [row["iter"] for row in values] == [0, 3, 6, 9, 12]
    seen = set()
    for rank in range(jobs):
        local_header, local_values = table(directory / f"job{rank}.params.txt")
        local_ids = [int(name[4:]) for name in local_header if name.startswith("lnl_")]
        assert local_ids == sorted(local_ids), local_ids
        assert not (seen & set(local_ids)), "gene assigned to more than one rank"
        seen.update(local_ids)
        assert len(local_values) == len(values)
        for global_row, local_row in zip(values, local_values):
            for name in local_header:
                assert global_row[name] == local_row[name], (name, global_row[name], local_row[name])
    assert seen == set(range(genes)), seen
    for row in values:
        for distance in ("dup", "los", "dco", "rfd", "lnl"):
            assert math.isclose(row[distance], sum(row[f"{distance}_{i}"] for i in range(genes)), abs_tol=1e-6)


def exercise(binary, root, mpi=None):
    directory = root / ("mpi" if mpi else "serial")
    fixture(directory)
    prefix = [mpi, "-n", "2", binary] if mpi else [binary]
    run(prefix + ["run.ctrl"], directory)
    run(prefix + ["-z", "1", "run.ctrl"], directory)
    check(directory, 6, 2 if mpi else 1)
    if mpi:
        # Analysis must reject inconsistent checkpoints without hanging in gathers.
        checkpoint = directory / "job1.checkpoint.bin"
        contents = checkpoint.read_bytes()
        checkpoint.write_bytes(contents[:len(contents) * 4 // 5])
        output = run(prefix + ["-z", "1", "run.ctrl"], directory, success=False)
        assert "different numbers of complete samples" in output, output
        invalid = root / "mpi-invalid-gene"
        fixture(invalid)
        tree = invalid / "g0.tre"
        tree.write_text(tree.read_text().replace("A_0", "UNKNOWN"))
        run(prefix + ["run.ctrl"], invalid, success=False)
        single = root / "mpi-single"
        fixture(single, 1)
        run([mpi, "-n", "1", binary, "run.ctrl"], single)
        run([mpi, "-n", "1", binary, "-z", "1", "run.ctrl"], single)
        check(single, 1, 1)
        output = run([mpi, "-n", "2", binary, "run.ctrl"], single, success=False)
        assert "no more MPI processes than gene families" in output, output
    else:
        for distances in ("1", "11110000", "111x000", "0000000"):
            run([binary, "-D", distances, "run.ctrl"], directory, success=False)
        (directory / "job0.checkpoint.bin").write_bytes(b"")
        output = run([binary, "-z", "1", "run.ctrl"], directory, success=False)
        assert "no complete samples" in output, output
    return header(directory / "params.txt")


def header(path):
    return path.read_text().splitlines()[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--serial", default=os.environ.get("GUENOMU_SERIAL"))
    parser.add_argument("--mpi", default=os.environ.get("GUENOMU_MPI"))
    parser.add_argument("--mpiexec", default=shutil.which("mpiexec"))
    args = parser.parse_args()
    assert args.serial or args.mpi, "provide --serial/--mpi or GUENOMU_SERIAL/GUENOMU_MPI"
    with tempfile.TemporaryDirectory(prefix="guenomu-output-") as temporary:
        root = Path(temporary)
        serial_header = exercise(str(Path(args.serial).resolve()), root) if args.serial else None
        if args.mpi:
            assert args.mpiexec, "mpiexec is required for MPI checks"
            mpi_header = exercise(str(Path(args.mpi).resolve()), root, args.mpiexec)
            if serial_header:
                assert mpi_header == serial_header, "serial/MPI columns differ"
    print("PASS: input-ordered parameter summaries, numeric columns, and checkpoint error handling")


if __name__ == "__main__":
    main()
