#!/usr/bin/env python3
"""Test Newick/NEXUS formats, uncertainty weights and rejected input errors."""
import math
import os
from pathlib import Path
import subprocess
import tempfile

helper = os.environ.get("GUENOMU_TREE_INPUT_TEST")
assert helper, "GUENOMU_TREE_INPUT_TEST must point to the built inspection helper"


def run(directory, text, expected=None, weights=None, diagnostic=None):
    path = directory / "input.trees"
    path.write_text(text)
    result = subprocess.run([helper, str(path)], text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, timeout=15)
    if diagnostic:
        assert result.returncode != 0, (text, result.stdout)
        assert diagnostic in result.stderr, (text, result.stderr)
        return
    assert result.returncode == 0, (text, result.stderr)
    rows = result.stdout.splitlines()
    summary = list(map(int, rows[0].split()[1:]))
    assert summary == expected, (text, summary, expected)
    if weights:
        actual = [float(row.split()[1]) for row in rows if row.startswith("TREE ")]
        assert len(actual) == len(weights)
        assert all(math.isclose(a, b) for a, b in zip(actual, weights)), actual
    return result.stdout


with tempfile.TemporaryDirectory(prefix="guenomu-tree-inputs-") as temporary:
    root = Path(temporary)
    run(root, "((A,B),(C,D));\n", [1, 1, 0, 0, 0, 4], [1.])
    run(root, "((A,B),(C,D)); ((A,B),(C,D));\n((A,C),(B,D));",
        [3, 2, 0, 0, 0, 4], [2 / 3, 1 / 3])
    run(root, "[leading [nested] comment]\n(\n(A,[across\nlines]B),\n(C,D)99\n)100;",
        [1, 1, 0, 0, 0, 4])
    run(root, "#NEXUS\nBegin trees; tree one = ((A,B),(C,D)); End;",
        [1, 1, 0, 0, 0, 4])
    out = run(root, "#NEXUS\nBegin trees; Translate 1 'A taxon',2 'O''Brien',3 'C:tip';\n"
              "tree 'name=1' = [&W 1/4] (1,(2,3)); tree two = [&W 3/4] (2,(1,3)); End;",
              [2, 2, 0, 0, 1, 3], [.25, .75])
    assert "'A taxon'" in out and "'O''Brien'" in out and "'C:tip'" in out
    out = run(root, "('A[;]taxon',('B,C','D)tip'));", [1, 1, 0, 0, 1, 3])
    assert "'A[;]taxon'" in out and "'B,C'" in out and "'D)tip'" in out
    run(root, "[&W .2] ((A,B),(C,D)); [&W .8] ((A,C),(B,D));",
        [2, 2, 0, 0, 0, 4], [.2, .8])
    out = run(root, "(A:1,B:2,C:3);", [1, 1, 1, 1, 1, 3])
    assert "LENGTH 6" in out
    run(root, "((A:0,B:0):0,(C:0,D:0):0);", [1, 1, 1, 1, 1, 4])
    run(root, "((A:1,B:1):.2,(C:1,D:1));", [1, 1, 1, 0, 0, 4])
    run(root, "((A,B):.2,(C,D):.3);", [1, 1, 1, 0, 1, 4])
    run(root, "(A,((B,C):.2,(D,E):.3));", [1, 1, 1, 0, 1, 5])
    run(root, "((A,B),(C,D)); ((A:1,B:1):1,(C:1,D:1):1);",
        [2, 1, 1, 0, 0, 4], [1.])
    for bad, message in [
        ("", "no trees"), ("#NEXUS\nBegin trees; End;", "no trees"),
        ("#NEXUS\nBegin trees; tree x=; End;", "expected a Newick"),
        ("((A,B),(C,D))", "semicolon"),
        ("((A,B),(C,D);", "unbalanced"),
        ("((A,B,C),D);", "bifurcating"),
        ("(A,B,C,D);", "bifurcating"),
        ("((A,A),(C,D));", "duplicate"),
        ("((A,B),(C,D)); ((A,B),(C,E));", "disagree"),
        ("((A:-1,B),(C,D));", "finite nonnegative"),
        ("((A:nan,B),(C,D));", "finite nonnegative"),
        ("((A:inf,B),(C,D));", "finite nonnegative"),
        ("((A:1foo,B),(C,D));", "characters after branch length"),
        ("[&W -1] ((A,B),(C,D));", "weights"),
        ("[&W 0] ((A,B),(C,D));", "positive sum"),
        ("[&W .5] ((A,B),(C,D)); ((A,C),(B,D));", "all or none"),
        ("[comment ((A,B),(C,D));", "unterminated"),
        ("(('bad,B),(C,D));", "unterminated"),
        ("#NEXUS\nBegin trees; Translate 1 A,2 B,3 C; tree x=(1,(1,3)); End;", "duplicate"),
        ("#NEXUS\nBegin trees; Translate 1 A,3 C; tree x=(1,3); End;", "cover"),
    ]:
        run(root, bad, diagnostic=message)
print("PASS: Newick/NEXUS tree inputs, quoted labels, uncertainty weights and branch provenance")
