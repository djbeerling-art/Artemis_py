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


def relocate_input(text, replacements):
    """Rewrite absolute paths in a published PHREEQC input.

    ``replacements`` maps old directory prefixes (e.g. the authors' HPC
    scratch directory) to new ones. Prefixes are matched literally, including
    the doubled slashes MATLAB's path concatenation leaves behind.
    """
    for old, new in replacements.items():
        text = re.sub(re.escape(old), new.replace("\\", "\\\\"), text)
    return text
