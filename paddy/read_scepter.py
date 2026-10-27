"""Read a SCEPTER run's profile records (prof/prof_*-NNN.txt) into {field: array(nrec, nz)}, time (yr) and z (m)."""
import glob, os, re
import numpy as np

def read(rundir, kinds=("aq", "gas", "sld", "aq(tot)", "ex(tot)")):
    out = {}
    for kind in kinds:
        files = [f for f in glob.glob(os.path.join(rundir, "prof", f"prof_{kind}-*.txt")) if re.search(r"-\d{3}\.txt$", f)]
        files.sort()
        recs = []
        for f in files:
            lines = open(f).read().split("\n")
            head = lines[0].replace("amfe3sid", "amfe3 sid").split()
            rows = np.array([[float(x) for x in l.split()] for l in lines[1:] if l.strip()])
            recs.append((head, rows))
        if not recs:
            continue
        head = recs[0][0]
        data = np.stack([r for _, r in recs])  # (nrec, nz, ncol)
        out.setdefault("z", data[0, :, 0])
        out.setdefault("time", data[:, 0, -1])
        if "time" in out and data.shape[0] == len(out["time"]) - 1:  # no record 000
            data = np.concatenate([np.full((1,) + data.shape[1:], np.nan), data])
        for j, h in enumerate(head[1:-1], start=1):
            key = h if kind in ("aq", "gas", "sld") else f"{h}{kind[2:] if kind.startswith('aq') else kind}"
            out[key] = data[:, :, j]
    return out

if __name__ == "__main__":
    import sys
    d = read(sys.argv[1])
    print(len(d["time"]), d["time"][:3], d["time"][-3:])
    print(sorted(d))
