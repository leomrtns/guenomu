#!/usr/bin/env python3
"""Exercise opt-in reference-tree expansion without altering supplied ensembles."""
import importlib.util
import os
from pathlib import Path
import shutil
import tempfile


def load(name):
    spec = importlib.util.spec_from_file_location(name.replace('-', '_'), Path(__file__).with_name(name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


checks = load('check-output')
mixture = load('check-multiple-species')


def configure(directory, samples=90, tau=0.2):
    mixture.configure(directory)
    text = (directory / 'run.ctrl').read_text().replace('end_of_control_file', f'''param_branch_uncertainty_samples = {samples}
param_branch_uncertainty_tau = {tau}
end_of_control_file''')
    (directory / 'run.ctrl').write_text(text)
    # Common phylogenetic input: only internal edges have lengths. Terminal
    # lengths are unnecessary because they carry no informative tree split.
    for tree in directory.glob('g*.tre'):
        tree.write_text(tree.read_text().replace(')', '):0.2'))


def exercise(binary, root, launcher=None):
    directory = root / ('mpi' if launcher else 'serial')
    configure(directory)
    command = [launcher, '-n', '2', binary] if launcher else [binary]
    output = checks.run(command + ['run.ctrl'], directory)
    assert 'expanded to 90 independent draws' in output, output
    generated = {}
    for gene in range(6):
        path = directory / f'uncertainty.gene{gene}.trees'
        lines = path.read_text().splitlines()
        assert len(lines) == 90, (gene, len(lines))
        generated[gene] = path.read_text()
        leafsets = [mixture.splits(line)[0] for line in lines]
        assert all(leaves == leafsets[0] for leaves in leafsets)
        assert len(leafsets[0]) >= 4
        assert len({frozenset(mixture.splits(line)[1]) for line in lines}) > 1
    # Repeated runs use per-original-family seeds and reproduce the input
    # evidence even if MPI has distributed families by a different order.
    checks.run(command + ['run.ctrl'], directory)
    for gene in range(6):
        assert (directory / f'uncertainty.gene{gene}.trees').read_text() == generated[gene]
    mixture.check(directory, 2)

    ensemble = root / ('mpi-ensemble' if launcher else 'serial-ensemble')
    configure(ensemble)
    p = ensemble / 'g0.tre'
    p.write_text(p.read_text().replace('End;', 'tree second = ((SpA_0,SpC_0),(SpB_0,(SpD_0,SpE_0)));\nEnd;'))
    checks.run(command + ['run.ctrl'], ensemble)
    assert not (ensemble / 'uncertainty.gene0.trees').exists(), 'supplied uncertainty must be preserved'

    identity = root / ('mpi-identity' if launcher else 'serial-identity')
    configure(identity, tau=0)
    checks.run(command + ['run.ctrl'], identity)
    original = (identity / 'multiple.params.tsv').read_bytes()
    assert not list(identity.glob('uncertainty.gene*.trees'))
    text = (identity / 'run.ctrl').read_text().replace('branch_uncertainty_samples = 90', 'branch_uncertainty_samples = 0')
    (identity / 'run.ctrl').write_text(text)
    checks.run(command + ['run.ctrl'], identity)
    assert original == (identity / 'multiple.params.tsv').read_bytes(), 'zero tau must leave sampler evidence and streams unchanged'

    invalid = root / ('mpi-missing' if launcher else 'serial-missing')
    configure(invalid)
    p = invalid / 'g0.tre'
    p.write_text(p.read_text().replace(':0.2', ''))
    checks.run(command + ['run.ctrl'], invalid, success=False)
    return generated


def main():
    serial = os.environ.get('GUENOMU_SERIAL')
    mpi = os.environ.get('GUENOMU_MPI')
    assert serial or mpi, 'provide a built GUENOMU_SERIAL and/or GUENOMU_MPI'
    with tempfile.TemporaryDirectory(prefix='guenomu-branch-') as name:
        root = Path(name)
        serial_trees = exercise(serial, root) if serial else None
        if mpi:
            launcher = shutil.which('mpiexec')
            assert launcher, 'MPI checks need mpiexec'
            mpi_trees = exercise(mpi, root, launcher)
            if serial_trees is not None:
                assert mpi_trees == serial_trees, 'expansion depends on MPI partition'
    print('PASS: branch uncertainty, deterministic family streams, input ensembles, zero-tau identity and invalid lengths')


if __name__ == '__main__':
    main()
