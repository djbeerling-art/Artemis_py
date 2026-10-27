"""Run the PHREEQC executable the way getphreeqcrun.m does.

getphreeqcrun.m calls ``phreeqc <fname>.in <fname>.out <dbfile> <fname>.screenout``
and then greps the .out file for ERROR, "Could not open" and WARNING lines.
"""

import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path


def find_phreeqc():
    """Locate the PHREEQC executable: $PHREEQC, then PATH, then $HOME/bin (the MATLAB default)."""
    candidates = [os.environ.get("PHREEQC"), shutil.which("phreeqc"),
                  str(Path.home() / "bin" / "phreeqc")]
    for c in candidates:
        if c and Path(c).is_file() and os.access(c, os.X_OK):
            return c
    raise FileNotFoundError("PHREEQC executable not found; set $PHREEQC or put phreeqc on PATH")


@dataclass
class RunResult:
    returncode: int
    out: Path
    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

    @property
    def ok(self):
        return self.returncode == 0 and not self.errors


def run(infile, dbfile, outfile=None, screenout=None, exe=None, cwd=None):
    """Run PHREEQC on ``infile`` and collect errors and warnings from the .out file."""
    infile = Path(infile)
    stem = infile.with_suffix("")
    outfile = Path(outfile) if outfile else stem.with_suffix(".out")
    screenout = Path(screenout) if screenout else stem.with_suffix(".screenout")
    exe = exe or find_phreeqc()
    proc = subprocess.run([exe, str(infile), str(outfile), str(dbfile), str(screenout)],
                          cwd=cwd or infile.parent, capture_output=True, text=True)
    text = outfile.read_text(errors="replace") if outfile.exists() else ""
    errors = [l for l in text.splitlines() if "ERROR" in l or "Could not open" in l]
    warnings = [l for l in text.splitlines() if "WARNING:" in l]
    return RunResult(proc.returncode, outfile, errors, warnings)


_BLOCK = re.compile(r"USE solution (\d+); USE equilibrium_phases")
_END = re.compile(r"(^|;)\s*END\s*($|[#;])")


def _split_day(text):
    """Split a getphreeqcrun day input into (head, preamble, mixing part, cell blocks, sel, dmp in, dmp out).

    Returns None unless the input has the daily layout: one KINETICS_MODIFY
    simulation per cell after the mixing, then a DUMP of everything."""
    L = text.splitlines(keepends=True)
    try:
        ititle = next(i for i, l in enumerate(L) if l.startswith("TITLE"))
        idump = max(i for i, l in enumerate(L) if l.startswith("DUMP"))
    except (StopIteration, ValueError):
        return None
    starts = [i for i in range(ititle, idump - 1) if _BLOCK.match(L[i]) and L[i + 1].startswith("KINETICS_MODIFY")]
    sel = re.search(r"^SELECTED_OUTPUT[^\n]*\n\s+-file (\S+)\n", text, re.M)
    dmpin = re.findall(r"^INCLUDE\$ (\S+\.dmp)\n", "".join(L[:ititle]), re.M)
    dmpout = re.match(r"\s+-file (\S+)", L[idump + 1]) if idump + 1 < len(L) else None
    if len(starts) < 2 or not sel or len(dmpin) != 1 or not dmpout:
        return None
    ipre = ititle + 1
    while ipre < starts[0] and not L[ipre].startswith(("SOLUTION", "USE ", "MIX", "EQUILIBRIUM", "SAVE")):
        ipre += 1
    blocks = [(int(_BLOCK.match(L[b]).group(1)), "".join(L[b:e])) for b, e in zip(starts, starts[1:] + [idump])]
    return ("".join(L[:ititle]), "".join(L[ititle:ipre]), "".join(L[ititle:starts[0]]), blocks,
            sel.group(1), dmpin[0], dmpout.group(1))


def run_cells_parallel(infile, dbfile, outfile=None, screenout=None, exe=None, cwd=None, workers=4):
    """Run one day's input with each cell's kinetics in its own PHREEQC process.

    After the day's mixing, each cell's KINETICS simulation uses only that
    cell's solution, so the cells can run at the same time. The mixing runs
    first and dumps its state; each cell then starts from that dump and dumps
    its own cell. The selected output (with serial simulation numbers), the
    dump and the .out file are merged into the files a serial run writes.
    Results match a serial run to rounding (the extra dump), except in cells
    with almost no water, where differences of ~1e-5 are possible. Falls back
    to :func:`run` for inputs without the daily layout (e.g. day 1).
    """
    from concurrent.futures import ThreadPoolExecutor

    infile = Path(infile)
    stem = infile.with_suffix("")
    outfile = Path(outfile) if outfile else stem.with_suffix(".out")
    screenout = Path(screenout) if screenout else stem.with_suffix(".screenout")
    parts = _split_day(infile.read_text())
    if parts is None:
        return run(infile, dbfile, outfile, screenout, exe, cwd)
    head, pre, mixing, blocks, sel, dmpin, dmpout = parts
    base = Path(cwd or infile.parent)
    selpath, dmppath = base / sel, base / dmpout
    exe = exe or find_phreeqc()
    tmp = infile.with_name(infile.name + ".split")
    tmp.mkdir(exist_ok=True)

    def go(name, body):
        f = tmp / f"{name}.in"
        f.write_text(body)
        proc = subprocess.run([exe, str(f), str(tmp / f"{name}.out"), str(dbfile), str(tmp / f"{name}.screenout")],
                              cwd=cwd or infile.parent, capture_output=True, text=True)
        return proc.returncode

    def selfile(name):
        return head.replace(f"-file {sel}\n", f"-file {tmp / name}.sel\n", 1)

    dumps = "".join(f"DUMP\n    -file {tmp}/A{k}.dmp\n    -cells {k}\nEND\n" for k, _ in blocks)
    rcA = go("A", selfile("A") + mixing + f"DUMP\n    -file {tmp / 'A.dmp'}\n    -all\nEND\n" + dumps)
    names = [f"K{k}" for k, _ in blocks]
    if rcA == 0:
        # each cell starts from the dump of its own cell only
        bodies = [selfile(n).replace(f"INCLUDE$ {dmpin}\n", f"INCLUDE$ {tmp}/A{k}.dmp\n", 1) + pre + b
                  + f"DUMP\n    -file {tmp / n}.dmp\n    -cells {k}\nEND\n" for n, (k, b) in zip(names, blocks)]
        with ThreadPoolExecutor(workers) as ex:
            rcs = list(ex.map(go, names, bodies))
    else:
        names, rcs = [], []
    # merge the .out files (for errors and warnings) and the selected output
    with open(outfile, "w") as fo:
        for n in ["A"] + names:
            fo.write((tmp / f"{n}.out").read_text(errors="replace") if (tmp / f"{n}.out").exists() else "")
    nsim = sum(1 for l in mixing.splitlines() if _END.search(l))
    with open(selpath, "w") as fo:
        if (tmp / "A.sel").exists():
            fo.write((tmp / "A.sel").read_text())
        for j, n in enumerate(names):
            f = tmp / f"{n}.sel"
            if not f.exists():
                continue
            for line in f.read_text().splitlines(keepends=True)[1:]:
                sim = line.split("\t", 1)[0]
                fo.write(f"{nsim + 1 + j:{len(sim)}d}" + line[len(sim):])
    if rcA == 0 and not any(rcs):
        with open(dmppath, "w") as fo:
            for n in ["A"] + names:
                fo.write((tmp / f"{n}.dmp").read_text())
    shutil.copy(tmp / "A.screenout", screenout) if (tmp / "A.screenout").exists() else None
    text = outfile.read_text(errors="replace")
    errors = [l for l in text.splitlines() if "ERROR" in l or "Could not open" in l]
    # every process repeats the warnings raised while reading the shared input; count those once
    first = (tmp / "A.out").read_text(errors="replace").splitlines() if (tmp / "A.out").exists() else []
    seen = {l for l in first if "WARNING:" in l}
    warnings = [l for l in first if "WARNING:" in l]
    warnings += [l for n in names if (tmp / f"{n}.out").exists()
                 for l in (tmp / f"{n}.out").read_text(errors="replace").splitlines() if "WARNING:" in l and l not in seen]
    return RunResult(max([rcA] + rcs), outfile, errors, warnings)


def relocate_input(text, replacements):
    """Rewrite absolute paths in a published PHREEQC input.

    ``replacements`` maps old directory prefixes (e.g. the authors' HPC
    scratch directory) to new ones. Prefixes are matched literally, including
    the doubled slashes MATLAB's path concatenation leaves behind.
    """
    for old, new in replacements.items():
        text = re.sub(re.escape(old), new.replace("\\", "\\\\"), text)
    return text
