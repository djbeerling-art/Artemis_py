"""List the MATLAB functions reachable from given entry points, with line counts.

Usage: python tools/matlab_callgraph.py DIR [DIR ...] ENTRY1,ENTRY2,...

Directories are searched in order, so put the example directory (whose
defaultsettings.m etc. shadow the top-level ones) first. A function counts as
called if its name appears as a whole word outside comments, which slightly
over-reports.
"""

import os
import re
import sys


def main(argv):
    dirs, roots = argv[:-1], argv[-1].split(",")
    funcs = {}
    for d in dirs:
        for f in sorted(os.listdir(d)):
            if f.endswith(".m"):
                funcs.setdefault(f[:-2], os.path.join(d, f))
    src = {n: re.sub(r"%.*", "", open(p, errors="replace").read()) for n, p in funcs.items()}
    calls = {n: sorted(m for m in funcs if m != n and re.search(r"\b" + re.escape(m) + r"\b", s))
             for n, s in src.items()}
    seen, queue = [], list(roots)
    while queue:
        n = queue.pop(0)
        if n in seen or n not in funcs:
            continue
        seen.append(n)
        queue += calls[n]
    total = 0
    for n in seen:
        lines = len(open(funcs[n], errors="replace").read().splitlines())
        total += lines
        print(f"{lines:5d} {os.path.relpath(funcs[n], dirs[-1])}  -> {','.join(calls[n])}")
    print("total", total, "files", len(seen))


if __name__ == "__main__":
    main(sys.argv[1:])
