"""Run the daily PHREEQC inputs in memory with IPhreeqc instead of the executable.

``phreeqc.run`` starts a new PHREEQC process each day, which reloads the
database, recompiles the RATES BASIC programs, and rebuilds the whole state
from the previous day's dump file. :class:`Engine` keeps one IPhreeqc instance
alive for the run instead: the state stays in memory, so the daily
``INCLUDE$ <run>.dmp`` and the closing ``DUMP`` are left out, and the input's
head (SELECTED_OUTPUT, USER_PUNCH and the PHASES/RATES includes) is read again
only when it changes. With ``workers > 1`` the cells' KINETICS simulations run
at the same time (as in ``phreeqc.run_cells_parallel``), each cell in its own
instance, with the cells' states passed between instances as in-memory dumps.

An Engine is called like ``phreeqc.run``. It writes the .sel file where the
executable would, and a .out file with the full output, or with just the
warnings and errors when ``printoutput`` is off.

The IPhreeqc shared library is found from $IPHREEQC_LIB, or as libIPhreeqc.so
on the library path. Build it from usgs-coupled/IPhreeqc; the authors' runs
used PHREEQC 3.7.3.
"""

import ctypes
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .phreeqc import _BLOCK, _END, RunResult

_lib = None


def load_library(path=None):
    """Load libIPhreeqc and declare the functions used here."""
    global _lib
    if _lib is not None and path is None:
        return _lib
    lib = ctypes.CDLL(path or os.environ.get("IPHREEQC_LIB") or "libIPhreeqc.so")
    sig = {
        "CreateIPhreeqc": ([], ctypes.c_int),
        "DestroyIPhreeqc": ([ctypes.c_int], ctypes.c_int),
        "LoadDatabase": ([ctypes.c_int, ctypes.c_char_p], ctypes.c_int),
        "RunString": ([ctypes.c_int, ctypes.c_char_p], ctypes.c_int),
        "GetErrorString": ([ctypes.c_int], ctypes.c_char_p),
        "GetWarningString": ([ctypes.c_int], ctypes.c_char_p),
        "GetOutputString": ([ctypes.c_int], ctypes.c_char_p),
        "GetSelectedOutputString": ([ctypes.c_int], ctypes.c_char_p),
        "GetDumpString": ([ctypes.c_int], ctypes.c_char_p),
        "GetVersionString": ([], ctypes.c_char_p),
    }
    for name in ("SetOutputStringOn", "SetSelectedOutputStringOn", "SetDumpStringOn", "SetOutputFileOn",
                 "SetSelectedOutputFileOn", "SetDumpFileOn", "SetErrorFileOn", "SetLogFileOn",
                 "SetErrorStringOn", "SetErrorOn"):
        sig[name] = ([ctypes.c_int, ctypes.c_int], ctypes.c_int)
    for name, (args, res) in sig.items():
        f = getattr(lib, name)
        f.argtypes, f.restype = args, res
    _lib = lib
    return lib


class Instance:
    """One IPhreeqc instance with string output on and file output off."""

    def __init__(self, dbfile, lib=None):
        self.lib = lib or load_library()
        self.id = self.lib.CreateIPhreeqc()
        if self.id < 0:
            raise RuntimeError("CreateIPhreeqc failed")
        # LoadDatabase resets the output switches, so it comes first
        if self.lib.LoadDatabase(self.id, str(dbfile).encode()):
            raise RuntimeError(f"IPhreeqc could not load {dbfile}: {self.errors()}")
        for f in ("SetOutputFileOn", "SetSelectedOutputFileOn", "SetDumpFileOn", "SetErrorFileOn", "SetLogFileOn"):
            getattr(self.lib, f)(self.id, 0)
        self.lib.SetErrorStringOn(self.id, 1)
        self.lib.SetSelectedOutputStringOn(self.id, 1)
        self.lib.SetDumpStringOn(self.id, 1)
        self.lib.SetOutputStringOn(self.id, 0)
        self.output = False
        self.head = None  # the SELECTED_OUTPUT and USER_PUNCH definitions last read
        self.heading = ""  # the selected output's heading line
        self.head_warnings = ""

    def run_with(self, head, text):
        """Run ``text``, reading ``head`` first only if it has changed.

        Reading SELECTED_OUTPUT and USER_PUNCH again takes about 0.04 s (the
        BASIC is recompiled); without them the selected output has no heading
        line, so :meth:`selected` adds the one last written."""
        if head == self.head:
            return self.run(text)
        n = self.run(head + text)
        self.head, self.heading = head, ""
        s = self.selected()
        self.heading = s.splitlines(keepends=True)[0] if s else ""
        self.head_warnings = self.warnings()
        return n

    def set_output(self, on):
        if on != self.output:
            self.lib.SetOutputStringOn(self.id, int(on))
            self.output = on

    def run(self, text):
        """Run ``text``; returns the number of errors."""
        return self.lib.RunString(self.id, text.encode())

    def _str(self, f):
        s = f(self.id)
        return s.decode(errors="replace") if s else ""

    def errors(self):
        return self._str(self.lib.GetErrorString)

    def warnings(self):
        return self._str(self.lib.GetWarningString)

    def output_string(self):
        return self._str(self.lib.GetOutputString)

    def selected(self):
        s = self._str(self.lib.GetSelectedOutputString)
        return s if not self.heading or s.startswith(self.heading) else self.heading + s

    def dump(self):
        return self._str(self.lib.GetDumpString)

    def close(self):
        if self.id >= 0:
            self.lib.DestroyIPhreeqc(self.id)
            self.id = -1


_INCLUDE = re.compile(r"^INCLUDE\$\s+(\S+)[ \t]*$", re.M)
_DUMPBLOCK = re.compile(r"^DUMP\n(?:\s+-[^\n]*\n)*END\n?", re.M)
_SELFILE = re.compile(r"^SELECTED_OUTPUT[^\n]*\n\s+-file (\S+)\n", re.M)
_RESTORE = re.compile(r"^INCLUDE\$\s+\S+\.dmp[ \t]*\n", re.M)
_RAW = re.compile(r"^\w+_RAW\s+(\d+)", re.M)


def _absolute_includes(text, base):
    """Make relative INCLUDE$ paths absolute: IPhreeqc reads them from the process's directory."""
    return _INCLUDE.sub(lambda m: f"INCLUDE$ {m.group(1) if os.path.isabs(m.group(1)) else os.path.join(base, m.group(1))}",
                        text)


def _split_dump(text):
    """Split a dump into {cell number: text} by its *_RAW blocks."""
    cells = {}
    pos = [(m.start(), int(m.group(1))) for m in _RAW.finditer(text)]
    for (s, k), (e, _) in zip(pos, pos[1:] + [(len(text), 0)]):
        cells.setdefault(k, []).append(text[s:e])
    return {k: "".join(v) for k, v in cells.items()}


_TOTAL_H = re.compile(r"^\s+-total_h\s+(\S+)", re.M)


def _water(dump):
    """Hydrogen in a cell's dumped solution (mol), a measure of its water; inf if there is none."""
    m = _TOTAL_H.search(dump)
    return float(m.group(1)) if m else float("inf")


def _split_body(body):
    """(preamble, mixing, [(cell, kinetics block)]) of a day's input from TITLE on; None without that layout."""
    L = body.splitlines(keepends=True)
    starts = [i for i in range(len(L) - 1) if _BLOCK.match(L[i]) and L[i + 1].startswith("KINETICS_MODIFY")]
    if len(starts) < 2:
        return None
    ipre = 1
    while ipre < starts[0] and not L[ipre].startswith(("SOLUTION", "USE ", "MIX", "EQUILIBRIUM", "SAVE")):
        ipre += 1
    blocks = [(int(_BLOCK.match(L[b]).group(1)), "".join(L[b:e])) for b, e in zip(starts, starts[1:] + [len(L)])]
    return "".join(L[:ipre]), "".join(L[:starts[0]]), blocks


class Engine:
    """A persistent IPhreeqc run for one PhreeqcRun; call it like ``phreeqc.run``.

    ``printoutput`` False leaves the .out file with the warnings and errors only.
    ``workers`` > 1 runs the cells' kinetics in parallel instances.
    """

    def __init__(self, dbfile, workers=1, printoutput=True, lib=None):
        self.dbfile = str(dbfile)
        self.workers = workers
        self.printoutput = printoutput
        self.lib = lib or load_library()
        self.main = Instance(dbfile, self.lib)
        self.pool = []
        self.cellinst = {}  # cell number: its instance
        self.includes = None  # the head's INCLUDE$ lines the instances have read
        self.rest = ""  # the rest of the head (SELECTED_OUTPUT, USER_PUNCH)
        self.dumpfile = None  # the dump the daily input includes; the state is in memory instead
        self.selfile = None
        self.executor = None
        self.cost = {}  # seconds each cell's kinetics took last time

    def close(self):
        if self.executor is not None:
            self.executor.shutdown()
        for p in [self.main] + self.pool:
            p.close()

    def write_dump(self, path):
        """Write the whole state as the executable's DUMP -all would."""
        self.main.run("DUMP\n    -all\nEND\n")
        Path(path).write_text(self.main.dump())

    def _head(self, head):
        """Split ``head`` into its definitions (SELECTED_OUTPUT, USER_PUNCH) and the INCLUDE$s still to read.

        The INCLUDE$d PHASES and RATES are read only when they change; the
        instances read the definitions only when they change (Instance.run_with)."""
        L = [l for l in head.splitlines(keepends=True)
             if not l.startswith("#") and not (self.dumpfile and l.strip() == f"INCLUDE$ {self.dumpfile}")]
        inc = "".join(l for l in L if l.startswith("INCLUDE$"))
        self.rest = "".join(l for l in L if not l.startswith("INCLUDE$"))
        if inc == self.includes:
            return ""
        self.includes = inc
        for p in self.pool:
            p.run(inc)
        return inc

    def __call__(self, infile, dbfile=None, outfile=None, screenout=None, exe=None, cwd=None):
        infile = Path(infile)
        base = str(Path(cwd or infile.parent))
        stem = infile.with_suffix("")
        outfile = Path(outfile) if outfile else stem.with_suffix(".out")
        text = infile.read_text()
        sel = _SELFILE.search(text)
        if sel:
            self.selfile = Path(base) / sel.group(1)
        # the state stays in memory: no DUMP, and no INCLUDE$ of yesterday's dump once running
        dmp = re.search(r"^DUMP\n\s+-file (\S+)", text, re.M)
        text = _DUMPBLOCK.sub("", text)
        ititle = re.search(r"^TITLE", text, re.M)
        head, body = (text[:ititle.start()], text[ititle.start():]) if ititle else (text, "")
        if self.includes is not None and self.dumpfile:
            head = head.replace(f"INCLUDE$ {self.dumpfile}\n", "")
        # a restart's dump is read once, by the main instance, ahead of the day's simulations
        restore = "".join(_RESTORE.findall(head))
        head = _RESTORE.sub("", head)
        head = self._head(_absolute_includes(head, base)) + _absolute_includes(restore, base)
        self.dumpfile = dmp.group(1) if dmp else self.dumpfile
        body = _absolute_includes(body, base)
        parts = _split_body(body) if self.workers > 1 else None
        res = self._serial(head + body) if parts is None else self._parallel(head, *parts)
        rc, sel, out, errors, warnings = res
        if self.selfile is not None:
            self.selfile.write_text(sel)
        if self.printoutput:
            outfile.write_text(out + "".join(e + "\n" for e in errors))
        else:
            outfile.write_text("".join(l + "\n" for l in warnings + errors))
        return RunResult(rc, outfile, errors, warnings)

    def _collect(self, inst, nerr):
        errors = [l for l in inst.errors().splitlines() if l.strip()] if nerr else []
        warnings = [l for l in inst.warnings().splitlines() if "WARNING:" in l]
        return errors, warnings

    def _main(self, text):
        """Run ``text`` in the main instance; warnings include those from reading the head, as each day's run does."""
        m = self.main
        m.set_output(self.printoutput)
        before = m.head
        nerr = m.run_with(self.rest, text)
        errors, warnings = self._collect(m, nerr)
        if before == m.head:
            warnings = [l for l in m.head_warnings.splitlines() if "WARNING:" in l] + warnings
        return nerr, errors, warnings

    def _serial(self, text):
        m = self.main
        nerr, errors, warnings = self._main(text)
        return (1 if nerr else 0), m.selected(), m.output_string() if self.printoutput else "", errors, warnings

    def _parallel(self, head, pre, mixing, blocks):
        # one instance per cell, so that each cell's results don't depend on which cells ran before
        # in its instance, and so on the number of workers
        for k, _ in blocks:
            if k not in self.cellinst:
                p = Instance(self.dbfile, self.lib)
                if self.includes:
                    p.run(self.includes)
                self.cellinst[k] = p
                self.pool.append(p)
        m = self.main
        cells = [k for k, _ in blocks]
        # the mixing, then a dump of the cells whose kinetics run next
        nerr, errors, warnings = self._main(head + mixing)
        sel = [m.selected()]
        out = [m.output_string() if self.printoutput else ""]
        if nerr:
            return 1, sel[0], out[0], errors, warnings
        nsim = sum(1 for l in mixing.splitlines() if _END.search(l))
        # The cell with the least water runs in the main instance, which already holds its state:
        # started from a dump, such a cell (e.g. frozen) can need several times as long to converge.
        # The choice depends only on the state, so a run gives the same results every time.
        # The other cells start, slowest first (by their last times), from dumps in the next free
        # pool instance.
        m.run("DUMP\n    -cells " + " ".join(map(str, cells)) + "\nEND\n")
        state = _split_dump(m.dump())
        inmain = min(cells, key=lambda k: (_water(state.get(k, "")), k))
        jobs = sorted(((j, k, b) for j, (k, b) in enumerate(blocks)), key=lambda job: -self.cost.get(job[1], -job[1]))

        def go(job):
            j, k, block = job
            p = m if k == inmain else self.cellinst[k]
            p.set_output(self.printoutput)
            t = time.perf_counter()
            if p is m:
                n = p.run(pre + block)
            else:
                n = p.run_with(self.rest, state.get(k, "") + pre + block + f"DUMP\n    -cells {k}\nEND\n")
                self.cost[k] = time.perf_counter() - t
            return j, (n, p.selected(), p.output_string() if self.printoutput else "",
                       p.errors() if n else "", p.warnings(), "" if p is m else p.dump())

        if self.executor is None:
            self.executor = ThreadPoolExecutor(self.workers)
        results = dict(self.executor.map(go, jobs))
        rc = 0
        dumps = []
        lines = sel[0].splitlines(keepends=True)
        cellwarn = [[l for l in results[j][4].splitlines() if "WARNING:" in l] for j in range(len(blocks))]
        # a warning every cell's run raises (from reading the head) is already in the main run's
        common = set.intersection(*(set(w) for w in cellwarn)) & set(warnings)
        for j in range(len(blocks)):
            n, s, o, e, _, d = results[j]
            rc = rc or (1 if n else 0)
            for line in s.splitlines(keepends=True)[1:]:
                sim = line.split("\t", 1)[0]
                lines.append(f"{nsim + 1 + j:{len(sim)}d}" + line[len(sim):])
            out.append(o)
            errors += [l for l in e.splitlines() if l.strip()]
            warnings += [l for l in cellwarn[j] if l not in common]
            dumps.append(d)
        if not rc:
            m.run("".join(dumps))
        return rc, "".join(lines), "".join(out), errors, warnings
