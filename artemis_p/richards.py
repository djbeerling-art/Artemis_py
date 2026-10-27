"""Daily soil water movement by the Richards equation. Python counterpart of
getrichquick.m and soilwatermovement.m (default zengdecker_2009 method), which
follow the CLM5 hydrology code (SoilWaterMovementMod.F90, SoilHydrologyMod.F90).

The MATLAB code keeps the state of the soil column in a structure ``c`` that is
handed from one day to the next; here it is a dict of 1-D numpy arrays (layers)
and scalars. Layer loops keep MATLAB's 1-based ``j`` and index with ``j - 1``.

Behaviour reproduced on purpose (see docs/PORTING_PLAN.md):

* On day 1 the ``c`` built by getrichquick has no ``zwtcm`` field, so
  soilwatermovement ignores it and starts from its own defaults.
* The extra CLM5 layers below the soil cells are only added on day 1 (later
  days see ``c.nwadded``), and the overlap fix for the first added layer is
  applied at index ``nlevsoi + 1`` instead of ``ncells + 1``.
* Initial water in the added layers uses ``clm5.dz(nwadd)`` (one value) rather
  than ``clm5.dz(wadd)``.
* When the water table is below the soil column, the inflow to the bottom soil
  layer is stored in a local ``qin`` rather than ``c.qin``, so it is left out of
  that layer's water balance in the tridiagonal system.
* WaterTable sets ``kfrz`` instead of ``k_frz`` when looking for the frost table.
"""

import math

import numpy as np

# CLM5 Technical Note Table 2.3: interface (bottom) and node depths (m)
CLM5_ZI = np.array([0.01, 0.04, 0.09, 0.16, 0.26, 0.40, 0.58, 0.8, 1.06, 1.36, 1.7, 2.08, 2.5, 2.99, 3.58,
                    4.27, 5.06, 5.95, 6.94, 8.03, 9.795, 13.328, 19.483, 28.871, 41.998])
CLM5_Z = np.array([0.02, 0.06, 0.12, 0.2, 0.32, 0.48, 0.68, 0.92, 1.2, 1.52, 1.88, 2.28, 2.72, 3.26, 3.9,
                   4.64, 5.48, 6.42, 7.46, 8.6, 10.99, 15.666, 23.301, 34.441, 49.556])
CLM5_ZTOP = np.concatenate([[0], CLM5_Z[:-1]])
CLM5_DZ = CLM5_Z - CLM5_ZTOP

DENICE, DENH2O = 917.0, 1000.0
HFUS = 3.337e5
GRAV = 9.80616
TFRZ = 273.15
CPICE, CPLIQ = 2.1172e3, 4.188e3
SMPMIN = -1e8
WATMIN = 0.01
AQUIFER_WATER_BASELINE = 5000.0
PONDMX = 0.0
SAT_LEV = 0.9
HKDEPTH = 1 / 2.5
PERCHED_BASEFLOW_SCALAR = 0.01
ORGANIC_MAX = 130.0


def tridiag(a, b, c, r):
    """Numerical Recipes tridiagonal solver as in soilwatermovement.m."""
    n = len(b)
    u = np.zeros(n)
    gam = np.zeros(n)
    bet = b[0]
    u[0] = r[0] / bet
    for j in range(1, n):
        gam[j] = c[j - 1] / bet
        bet = b[j] - a[j] * gam[j]
        u[j] = (r[j] - a[j] * u[j - 1]) / bet
    for j in range(n - 2, -1, -1):
        u[j] = u[j] - gam[j + 1] * u[j + 1]
    return u


def _grow(x, n):
    """MATLAB assignment beyond the end of an array pads it with zeros."""
    if len(x) >= n:
        return x
    return np.concatenate([x, np.zeros(n - len(x))])


class _SoilWaterMovement:
    """One call of soilwatermovement.m: the nested functions share ``self``."""

    def __init__(self, ds, c, idx, addice=1, etscale=1.0):
        parm = ds["parm"]
        self.idx = idx  # 1-based day
        self.addice, self.etscale = addice, etscale
        self.e_ice = parm.get("e_ice", 6.0)
        self.aq_sp_yield_min = parm.get("aq_sp_yield_min", 0.1)
        self.bedrockdepthm = parm.get("bedrockdepthm", 5.0)
        self.nlevsoi = 20
        if c is None or "zwtcm" not in c:  # fresh structure (day 1)
            c = {"soilflag": np.r_[np.ones(20), np.zeros(5)], "nlevsoi": 20, "nlevgrnd": 25, "ncells": 20,
                 "wa": 4000.0, "snl": 0, "frac_h2osfc": 0.0, "qflx_liqdew_to_top_layer": 0.0,
                 "qflx_soliddew_to_top_layer": 0.0, "qflx_solidevap_from_top_layer": 0.0,
                 "qflx_ev_snow": 0.0, "topo_slope": 0.0, "draindepthm": 0.0}
        self.c = c
        self.getds(ds)

    # ------------------------------------------------------------------ getds
    def getds(self, ds):
        c, idx, parm, cell, run = self.c, self.idx, ds["parm"], ds["cell"], ds["run"]
        if c["draindepthm"] <= 0:
            if "draindepthm" in parm:
                c["draindepthm"] = parm["draindepthm"]
        nwadd = 0
        c["ncells"] = nc = ds["ncells"]
        bottom = np.asarray(cell["bottomdepth"], float)
        if "nwadded" not in c:
            wadd = np.flatnonzero(CLM5_ZI * 100 > bottom[-1])
            nwadd = len(wadd)
            c["nwadded"] = nwadd
            if nwadd > 0:
                c["soilflag"] = np.r_[np.ones(nc), c["soilflag"][wadd]]
                c["nlevsoi"] = int(c["soilflag"].sum())
                c["nlevgrnd"] = len(c["soilflag"])
        self.nlevsoi = nlevsoi = c["nlevsoi"]
        addaquifer = nwadd > 0
        c.setdefault("dtime", 86400.0)
        if "zmm" not in c:
            c["z"] = np.asarray(cell["middepth"], float) / 100
            c["zi"] = bottom / 100
            c["ztop"] = np.asarray(cell["topdepth"], float) / 100
            c["dz"] = np.asarray(cell["thickness"], float) / 100
            if addaquifer:
                c["z"] = np.r_[c["z"], CLM5_Z[wadd]]
                c["ztop"] = np.r_[c["ztop"], CLM5_ZTOP[wadd]]
                c["dz"] = np.r_[c["dz"], CLM5_DZ[wadd]]
                c["zi"] = np.r_[c["zi"], CLM5_ZI[wadd]]
                # meant for the first added layer, but nlevsoi already includes the added soil layers
                c["ztop"][nlevsoi] = c["zi"][nlevsoi - 1]
                c["dz"][nlevsoi] = c["zi"][nlevsoi - 1] - c["ztop"][nlevsoi - 1]
            c["zmm"], c["zimm"], c["dzmm"] = c["z"] * 1000, c["zi"] * 1000, c["dz"] * 1000
        for k in ("watsat", "hksat", "sucsat", "bsw"):
            if k not in c:
                v = np.asarray(cell[k], float).ravel()
                c[k] = np.r_[v, np.zeros(nwadd) + v[-1]] if addaquifer else v.copy()
        if "zwt" not in c or c["zwt"] <= 0:
            if "watertabledepth" in run and idx > 1:
                c["zwt"] = run["watertabledepth"][idx - 2] / 100
            elif c.get("draindepthm", 0) > 0:
                c["zwt"] = c["draindepthm"]
            elif parm.get("draindepthm", 0) > 0:
                c["zwt"] = parm["draindepthm"]
                c["draindepthm"] = parm["draindepthm"]
            else:
                c["zwt"] = c["zi"][nlevsoi - 1] + 1
        wsat = np.flatnonzero(c["zi"] > c["zwt"])
        if "om_frac" not in c:
            c["om_frac"] = (np.asarray(cell["fSOM"], float) * np.asarray(cell["gsoil"], float)
                            / np.asarray(cell["vsoil"], float) / ORGANIC_MAX)
        c.setdefault("sand", np.asarray(cell["sandpct"], float).copy())
        c.setdefault("clay", np.asarray(cell["claypct"], float).copy())
        c["qflx_infl"] = run["infil"][idx - 1] / 86400
        temp = np.array(run["soiltemp"][:, idx - 1], float)
        wnan = np.flatnonzero(np.isnan(temp))
        if len(wnan) > 0 and wnan[0] == 0:
            print(f"(soilwatermovement) NaN for topsoil temp at idx={idx}")
        elif len(wnan) > 0:
            temp[wnan] = temp[wnan[0] - 1]
        t = c.get("t_soisno", np.zeros(nc))
        t = _grow(np.array(t, float), nc)
        t[:nc] = temp + 273.15
        c["t_soisno"] = t
        q = _grow(np.array(c.get("qflx_rootsoi_col", np.zeros(nc)), float), nc)
        q[:nc] = run["cellET"][:, idx - 1] * self.etscale / 86400
        if "qdrain" in run:
            q[:nc] = q[:nc] + run["qdrain"][:, max(idx - 2, 0)] / 86400
        c["qflx_rootsoi_col"] = q
        if addaquifer:
            c["om_frac"] = np.r_[c["om_frac"], np.zeros(nwadd)]
            c["sand"] = np.r_[c["sand"], np.zeros(nwadd) + cell["sandpct"][-1]]
            c["clay"] = np.r_[c["clay"], np.zeros(nwadd) + cell["claypct"][-1]]
            c["t_soisno"] = np.r_[c["t_soisno"], np.zeros(nwadd) + c["t_soisno"][-1]]
            c["qflx_rootsoi_col"] = np.r_[c["qflx_rootsoi_col"], np.zeros(nwadd)]
        if "h2osoi_liq" not in c:
            soiltemp = run["soiltemp"][:, idx - 1]
            water = np.asarray(cell["targH2O"], float) if idx == 1 else run["targH2O"][:, idx - 2]
            c["h2osoi_liq"] = water * (soiltemp > 0)
            c["h2osoi_ice"] = water * (soiltemp <= 0)
            if addaquifer:
                fc = np.asarray(cell["fc"], float)
                c["h2osoi_liq"] = np.r_[c["h2osoi_liq"], np.zeros(nwadd) + fc[-1] * CLM5_DZ[nwadd - 1] * 1000]
                if len(wsat) > 0:
                    c["h2osoi_liq"][wsat] = c["watsat"][wsat] * c["dz"][wsat] * 1000
                c["h2osoi_ice"] = np.r_[c["h2osoi_ice"], np.zeros(nwadd)]
            c["h2osoi_vol"] = (c["h2osoi_liq"] + c["h2osoi_ice"]) / c["dz"] / 1000
        if "sublimation" in run:
            c["h2osno"] = run["snowpack"][idx - 1]
            c["h2osno0"] = run["snowpack"][idx - 2] if idx > 1 else c["h2osno"]
            c["h2osno_no_layers"] = run["snowpack"][idx - 1] + run["snowmelt"][idx - 1]
            if c["h2osno_no_layers"] > 10:
                c["h2osno_no_layers"] = 0.0
            c["frac_h2osfc"] = 0.0
            c["frac_sno_eff"] = max(1, c["h2osno"])
            c["qflx_ev_snow"] = run["sublimation"][idx - 1]
            if run["snowpack"][idx - 1] > 0:
                c["snl"] = -1
        parm.setdefault("soilwater_movement_method", "zengdecker_2009")
        if parm["soilwater_movement_method"] != "zengdecker_2009":
            raise NotImplementedError("only the zengdecker_2009 soil water method is ported")
        parm.setdefault("lower_boundary_condition", "aquifer")
        parm.setdefault("upper_boundary_condition", "flux")

    # ---------------------------------------------------------------- helpers
    def getjwt(self):
        c = self.c
        for j in range(1, self.nlevsoi + 1):
            if c["zwt"] <= c["zi"][j - 1]:
                return j - 1
        return self.nlevsoi

    def set_soil_water_fractions(self):
        c = self.c
        c["vol_ice"] = np.minimum(c["watsat"], c["h2osoi_ice"] / (c["dz"] * DENICE))
        c["eff_porosity"] = np.maximum(0.01, c["watsat"] - c["vol_ice"])
        c["icefrac"] = np.minimum(1, c["vol_ice"] / c["watsat"])

    def heat_capacity(self):
        """cv from paramsforfreezing (the other outputs are not used)."""
        c = self.c
        om, sand, clay = c["om_frac"], c["sand"], c["clay"]
        om_csol = 2.5
        csol = ((1 - om) * (2.128 * sand + 2.385 * clay) / (sand + clay) + om_csol * om) * 1e6
        return csol * (1 - c["watsat"]) * c["dz"] + (c["h2osoi_ice"] * CPICE + c["h2osoi_liq"] * CPLIQ)

    # ------------------------------------------------------------ Phasechange
    def phasechange(self):
        c, nlevsoi = self.c, self.nlevsoi
        cv = self.heat_capacity()
        fact = c["dtime"] / cv
        # paramsforfreezing overrides the value set in getds
        c["h2osno_no_layers"] = c["h2osno"] if 0 < c.get("h2osno", 0) < 10 else 0.0
        imelt = np.zeros(nlevsoi, int)
        hm, xm, supercool, tinc = (np.zeros(nlevsoi) for _ in range(4))
        wice0 = c["h2osoi_ice"].copy()
        wliq0 = c["h2osoi_liq"].copy()
        wmass0 = wice0 + wliq0
        t = c["t_soisno"]
        for j in range(1, nlevsoi + 1):
            k = j - 1
            if c["h2osoi_ice"][k] > 0 and t[k] >= TFRZ:
                imelt[k] = 1
                tinc[k] = TFRZ - t[k]
            if t[k] < TFRZ:
                smp = HFUS * (TFRZ - t[k]) / (GRAV * t[k]) * 1000
                supercool[k] = c["watsat"][k] * (smp / c["sucsat"][k]) ** (-1 / c["bsw"][k])
                supercool[k] = supercool[k] * c["dz"][k] * 1000
            if c["h2osoi_liq"][k] > supercool[k] and t[k] <= TFRZ:
                imelt[k] = 2
                tinc[k] = TFRZ - t[k]
            if c["h2osno_no_layers"] > 0 and j == 1:
                if t[k] >= TFRZ:
                    imelt[k] = 1
                    tinc[k] = TFRZ - t[k]
        for j in range(1, nlevsoi + 1):
            k = j - 1
            if imelt[k] > 0:
                hm[k] = -tinc[k] / fact[k]
            if imelt[k] == 1 and hm[k] < 0:
                hm[k], imelt[k] = 0, 0
            if imelt[k] == 2 and hm[k] > 0:
                hm[k], imelt[k] = 0, 0
            if imelt[k] > 0 and abs(hm[k]) > 0:
                xm[k] = hm[k] * c["dtime"] / HFUS
            if j == 1 and c["h2osno_no_layers"] > 0 and xm[k] > 0:
                temp1 = c["h2osno_no_layers"]
                c["h2osno_no_layers"] = max(0, temp1 - xm[k])
                heatr = hm[k] - HFUS * (temp1 - c["h2osno_no_layers"]) / c["dtime"]
                if heatr > 0:
                    xm[k] = heatr * c["dtime"] / HFUS
                    hm[k] = heatr
                else:
                    xm[k], hm[k] = 0, 0
            if xm[k] > 0:
                c["h2osoi_ice"][k] = max(0, wice0[k] - xm[k])
            elif xm[k] < 0:
                if wmass0[k] < supercool[k]:
                    c["h2osoi_ice"][k] = 0
                else:
                    c["h2osoi_ice"][k] = min(wmass0[k] - supercool[k], wice0[k] - xm[k])
            c["h2osoi_liq"][k] = max(0, wmass0[k] - c["h2osoi_ice"][k])

    # ---------------------------------------------------------- Richards eqn
    def soilwater_zengdecker2009(self):
        c, nlevsoi, e_ice = self.c, self.nlevsoi, self.e_ice
        self.set_soil_water_fractions()
        c["vwc_liq"] = np.maximum(c["h2osoi_liq"], 1.0e-6) / (c["dz"] * DENH2O)
        jwt = self.getjwt()
        n = nlevsoi
        watsat, sucsat, bsw, hksat = c["watsat"], c["sucsat"], c["bsw"], c["hksat"]
        t = c["t_soisno"]
        c["vwc_zwt"] = watsat[n - 1]
        if t[jwt] < TFRZ:  # t_soisno(jwt+1)
            c["vwc_zwt"] = c["vwc_liq"][n - 1]
            for j in range(n, 0, -1):
                if c["zwt"] <= c["zi"][j - 1]:
                    smp1 = HFUS * (TFRZ - t[j - 1]) / (GRAV * t[j - 1]) * 1000
                    smp1 = max(smp1, sucsat[n - 1])
                    c["vwc_zwt"] = watsat[n - 1] * (smp1 / sucsat[n - 1]) ** (-1 / bsw[n - 1])
                    c["vwc_zwt"] = min(c["vwc_zwt"], 0.5 * (watsat[n - 1] + c["h2osoi_vol"][n - 1]))
                    break
        c["zimm0"] = c["zimm"] - c["dzmm"]
        zwtmm, zimm, zimm0 = c["zwtmm"], c["zimm"], c["zimm0"]
        vol_eq, zq = np.zeros(n + 1), np.zeros(n + 1)
        for j in range(1, n + 1):
            k = j - 1
            if zwtmm <= zimm0[k]:
                vol_eq[k] = watsat[k]
            elif zwtmm < zimm[k] and zwtmm > zimm0[k]:
                tempi = 1.0
                temp0 = ((sucsat[k] + zwtmm - zimm0[k]) / sucsat[k]) ** (1 - 1 / bsw[k])
                voleq1 = -sucsat[k] * watsat[k] / (1 - 1 / bsw[k]) / (zwtmm - zimm0[k]) * (tempi - temp0)
                vol_eq[k] = (voleq1 * (zwtmm - zimm0[k]) + watsat[k] * (zimm[k] - zwtmm)) / (zimm[k] - zimm0[k])
                vol_eq[k] = max(0, min(watsat[k], vol_eq[k]))
            else:
                tempi = ((sucsat[k] + zwtmm - zimm[k]) / sucsat[k]) ** (1 - 1 / bsw[k])
                temp0 = ((sucsat[k] + zwtmm - zimm0[k]) / sucsat[k]) ** (1 - 1 / bsw[k])
                vol_eq[k] = -sucsat[k] * watsat[k] / (1 - 1 / bsw[k]) / (zimm[k] - zimm0[k]) * (tempi - temp0)
                vol_eq[k] = min(watsat[k], max(0, vol_eq[k]))
            zq[k] = -sucsat[k] * (max(0.01, vol_eq[k] / watsat[k])) ** (-bsw[k])
            zq[k] = max(SMPMIN, zq[k])
        if jwt == n:
            k = n - 1
            tempi = 1.0
            temp0 = ((sucsat[k] + zwtmm - zimm[k]) / sucsat[k]) ** (1 - 1 / bsw[k])
            vol_eq[n] = -sucsat[k] * watsat[k] / (1 - 1 / bsw[k]) / (zwtmm - zimm[k]) * (tempi - temp0)
            vol_eq[n] = min(watsat[k], max(0, vol_eq[n]))
            zq[n] = -sucsat[k] * (max(0.01, vol_eq[n] / watsat[k])) ** (-bsw[k])
            zq[n] = max(SMPMIN, zq[n])

        vwc = c["vwc_liq"]
        imped, hk, dhkdw, smp, dsmpdw = (np.zeros(n) for _ in range(5))
        hk = _grow(np.array(c.get("hk", np.zeros(n)), float), n)
        for j in range(1, n + 1):
            k = j - 1
            km = min(n, j + 1) - 1
            s1 = 0.5 * (vwc[k] + vwc[km]) / (0.5 * (watsat[k] + watsat[km]))
            s1 = min(1, s1)
            s2 = hksat[k] * s1 ** (2 * bsw[k] + 2)
            imped[k] = 10 ** (-e_ice * (0.5 * (c["icefrac"][k] + c["icefrac"][km])))
            hk[k] = imped[k] * s1 * s2
            dhkdw[k] = imped[k] * (2 * bsw[k] + 3) * s2 * (1 / (watsat[k] + watsat[km]))
            s_node = min(1, max(0.01, vwc[k] / watsat[k]))
            smp[k] = -sucsat[k] * s_node ** (-bsw[k])
            smp[k] = max(SMPMIN, smp[k])
            dsmpdw[k] = -bsw[k] * smp[k] / vwc[k]
        c["hk"] = hk

        zmm, dzmm = c["zmm"], c["dzmm"]
        zmm[n] = 0.5 * (zwtmm + zmm[n - 1])
        dzmm[n] = dzmm[n - 1] if jwt < n else zwtmm - zmm[n - 1]
        dtime, et = c["dtime"], c["qflx_rootsoi_col"]
        qin, qout = c["qin"], c["qout"]
        rmx, amx, bmx, cmx = (np.zeros(n + 1) for _ in range(4))
        dqidw0, dqidw1, dqodw1, dqodw2 = (np.zeros(n + 1) for _ in range(4))
        numj, denj = np.zeros(n), np.zeros(n)
        sdamp = 0.0
        # top layer
        qin[0] = c["qflx_infl"]
        den = zmm[1] - zmm[0]
        dzq = zq[1] - zq[0]
        num = (smp[1] - smp[0]) - dzq
        numj[0], denj[0] = num, den
        qout[0] = -hk[0] * num / den
        dqodw1[0] = -(-hk[0] * dsmpdw[0] + num * dhkdw[0]) / den
        dqodw2[0] = -(hk[0] * dsmpdw[1] + num * dhkdw[0]) / den
        rmx[0] = qin[0] - qout[0] - et[0]
        amx[0] = 0
        bmx[0] = dzmm[0] * (sdamp + 1 / dtime) + dqodw1[0]
        cmx[0] = dqodw2[0]
        for j in range(2, n):
            k = j - 1
            den = zmm[k] - zmm[k - 1]
            dzq = zq[k] - zq[k - 1]
            num = (smp[k] - smp[k - 1]) - dzq
            numj[k], denj[k] = num, den
            qin[k] = -hk[k - 1] * num / den
            dqidw0[k] = -(-hk[k - 1] * dsmpdw[k - 1] + num * dhkdw[k - 1]) / den
            dqidw1[k] = -(hk[k - 1] * dsmpdw[k] + num * dhkdw[k - 1]) / den
            den = zmm[k + 1] - zmm[k]
            dzq = zq[k + 1] - zq[k]
            num = (smp[k + 1] - smp[k]) - dzq
            qout[k] = -hk[k] * num / den
            dqodw1[k] = -(-hk[k] * dsmpdw[k] + num * dhkdw[k]) / den
            dqodw2[k] = -(hk[k] * dsmpdw[k + 1] + num * dhkdw[k]) / den
            rmx[k] = qin[k] - qout[k] - et[k]
            amx[k] = -dqidw0[k]
            bmx[k] = dzmm[k] / dtime - dqidw1[k] + dqodw1[k]
            cmx[k] = dqodw2[k]
        k = n - 1  # bottom layer
        if n > jwt:  # water table in soil column
            den = zmm[k] - zmm[k - 1]
            dzq = zq[k] - zq[k - 1]
            num = (smp[k] - smp[k - 1]) - dzq
            numj[k], denj[k] = num, den
            qin[k] = -hk[k - 1] * num / den
            dqidw0[k] = -(-hk[k - 1] * dsmpdw[k - 1] + num * dhkdw[k - 1]) / den
            dqidw1[k] = -(hk[k - 1] * dsmpdw[k] + num * dhkdw[k - 1]) / den
            qout[k] = 0
            dqodw1[k] = 0
            rmx[k] = qin[k] - qout[k] - et[k]
            amx[k] = -dqidw0[k]
            bmx[k] = dzmm[k] / dtime - dqidw1[k] + dqodw1[k]
            cmx[k] = 0
            rmx[n], amx[n], bmx[n], cmx[n] = 0, 0, dzmm[n] / dtime, 0
        else:  # water table below bottom soil layer
            s_node = min(1, max(0.01, 0.5 * ((c["vwc_zwt"] + vwc[k]) / watsat[k])))
            smp1 = max(SMPMIN, -sucsat[k] * s_node ** (-bsw[k]))
            dsmpdw1 = -bsw[k] * smp1 / (s_node * watsat[k])
            den = zmm[k] - zmm[k - 1]
            dzq = zq[k] - zq[k - 1]
            num = (smp[k] - smp[k - 1]) - dzq
            numj[k], denj[k] = num, den
            # MATLAB assigns this to a local qin, so c.qin(nlevsoi) stays 0 below
            dqidw0[k] = -(-hk[k - 1] * dsmpdw[k - 1] + num * dhkdw[k - 1]) / den
            dqidw1[k] = -(hk[k - 1] * dsmpdw[k] + num * dhkdw[k - 1]) / den
            den = zmm[k + 1] - zmm[k]
            dzq = zq[k + 1] - zq[k]
            num = (smp1 - smp[k]) - dzq
            numj[k], denj[k] = num, den
            qout[k] = -hk[k] * num / den
            dqodw1[k] = -(-hk[k] * dsmpdw[k] + num * dhkdw[k]) / den
            dqodw2[k] = -(hk[k] * dsmpdw1 + num * dhkdw[k]) / den
            rmx[k] = qin[k] - qout[k] - et[k]
            amx[k] = -dqidw0[k]
            bmx[k] = dzmm[k] / dtime - dqidw1[k] + dqodw1[k]
            cmx[k] = dqodw2[k]
            qin[n] = qout[k]
            dqidw0[n] = -(-hk[k] * dsmpdw[k] + num * dhkdw[k]) / den
            dqidw1[n] = -(hk[k] * dsmpdw1 + num * dhkdw[k]) / den
            qout[n] = 0
            dqodw1[n] = 0
            rmx[n] = qin[n] - qout[n]
            amx[n] = -dqidw0[n]
            bmx[n] = dzmm[n] / dtime - dqidw1[n] + dqodw1[n]
            cmx[n] = 0
        dwat2 = tridiag(amx, bmx, cmx, rmx)
        for k in range(n):
            c["h2osoi_liq"][k] = c["h2osoi_liq"][k] + dwat2[k] * dzmm[k]
        if jwt < n:
            wh_zwt = 0.0
            s_node = max(0.01, c["h2osoi_vol"][jwt] / watsat[jwt])
            s1 = min(1, s_node)
            ka = imped[jwt] * hksat[jwt] * s1 ** (2 * bsw[jwt] + 3)
            smp1 = max(SMPMIN, smp[max(1, jwt) - 1])
            wh = smp1 - zq[max(1, jwt) - 1]
            if jwt == 0:
                c["qcharge"] = -ka * (wh_zwt - wh) / ((c["zwt"] + 0.001) * 1000)
            else:
                c["qcharge"] = -ka * (wh_zwt - wh) / ((c["zwt"] - c["z"][jwt - 1]) * 1000 * 2.0)
            c["qcharge"] = min(10 / dtime, max(-10 / dtime, c["qcharge"]))
        else:
            c["qcharge"] = dwat2[n] * dzmm[n] / dtime
        c["qflx_deficit"] = 0.0
        for k in range(n):
            if c["h2osoi_liq"][k] < 0:
                c["qflx_deficit"] -= c["h2osoi_liq"][k]
        c["smp"], c["imped"], c["num"], c["den"] = smp, imped, numj, denj
        c["dwat2"] = dwat2
        c["h2osoi_liq0"] = c["h2osoi_liq"].copy()

    # -------------------------------------------------------------- WaterTable
    def _specific_yield(self, k):
        c = self.c
        s_y = c["watsat"][k] * (1 - (1 + 1000 * c["zwt"] / c["sucsat"][k]) ** (-1 / c["bsw"][k]))
        return max(s_y, self.aq_sp_yield_min)

    def _perched_level(self, k_frz):
        """Locate k_perch from the frost table up, updating h2osoi_vol on the way."""
        c = self.c
        k_perch = 1
        for k in range(k_frz, 0, -1):
            c["h2osoi_vol"][k - 1] = (c["h2osoi_liq"][k - 1] / (c["dz"][k - 1] * DENH2O)
                                      + c["h2osoi_ice"][k - 1] / (c["dz"][k - 1] * DENICE))
            if c["h2osoi_vol"][k - 1] / c["watsat"][k - 1] <= SAT_LEV:
                k_perch = k
                break
        if c["t_soisno"][k_frz - 1] > TFRZ:
            k_perch = k_frz
        return k_perch

    def _zwt_perched(self, k_perch):
        c = self.c
        sat = [(c["h2osoi_liq"][k] / (c["dz"][k] * DENH2O) + c["h2osoi_ice"][k] / (c["dz"][k] * DENICE))
               / c["watsat"][k] for k in (k_perch - 1, k_perch)]
        m = (c["z"][k_perch] - c["z"][k_perch - 1]) / (sat[1] - sat[0])
        b = c["z"][k_perch] - m * sat[1]
        return max(0, m * SAT_LEV + b)

    def water_table(self):
        c, nlevsoi = self.c, self.nlevsoi
        c["dzmm"][:nlevsoi] = c["dz"][:nlevsoi] * 1000
        jwt = self.getjwt()
        rous = max(self._specific_yield(nlevsoi - 1), self.aq_sp_yield_min)
        if jwt == nlevsoi:
            c["wa"] = c["wa"] + c["qcharge"] * c["dtime"]
            c["zwt"] = c["zwt"] - (c["qcharge"] * c["dtime"]) / 1000 / rous
        else:
            qcharge_tot = c["qcharge"] * c["dtime"]
            if qcharge_tot > 0:
                for j in range(jwt + 1, 0, -1):
                    s_y = self._specific_yield(j - 1)
                    ziup = 0 if j == 1 else c["zi"][j - 2]
                    qcharge_layer = max(min(qcharge_tot, s_y * (c["zwt"] - ziup) * 1000), 0)
                    if s_y > 0:
                        c["zwt"] = c["zwt"] - qcharge_layer / s_y / 1000
                    qcharge_tot -= qcharge_layer
                    if qcharge_tot <= 0:
                        break
            else:
                for j in range(jwt + 1, nlevsoi + 1):
                    s_y = self._specific_yield(j - 1)
                    qcharge_layer = min(0, max(qcharge_tot, -(s_y * (c["zi"][j - 1] - c["zwt"]) * 1000)))
                    qcharge_tot -= qcharge_layer
                    if qcharge_tot >= 0:
                        c["zwt"] = c["zwt"] - qcharge_layer / s_y / 1000
                        break
                    c["zwt"] = c["zi"][j - 1]
                if qcharge_tot > 0:
                    c["zwt"] = c["zwt"] - qcharge_tot / 1000 / rous
        t = c["t_soisno"]
        k_frz = nlevsoi if t[0] > TFRZ else 1
        # the MATLAB loop sets kfrz, not k_frz, so the frost table search has no effect
        frost_table = c["z"][k_frz - 1]
        if not (c["zwt"] < frost_table and t[k_frz - 1] <= TFRZ):
            k_perch = self._perched_level(k_frz)
            if k_frz > k_perch:
                self._zwt_perched(k_perch)
            if c["snl"] + 1 >= 1:  # dew and sublimation terms are all zero here
                c["h2osoi_liq"][0] += (1 - c["frac_h2osfc"]) * c["qflx_liqdew_to_top_layer"] * c["dtime"]
                c["h2osoi_ice"][0] += (1 - c["frac_h2osfc"]) * c["qflx_soliddew_to_top_layer"] * c["dtime"]
                if c["qflx_solidevap_from_top_layer"] * c["dtime"] > c["h2osoi_ice"][0]:
                    c["qflx_solidevap_from_top_layer"] = c["h2osoi_ice"][0] / c["dtime"]
                    c["h2osoi_ice"][0] = 0
                else:
                    c["h2osoi_ice"][0] -= (1 - c["frac_h2osfc"]) * c["qflx_solidevap_from_top_layer"] * c["dtime"]

    # ---------------------------------------------------------------- Drainage
    def drainage(self, jdrain):
        c, nlevsoi, e_ice = self.c, self.nlevsoi, self.e_ice
        for k in range(nlevsoi):
            vol_ice = min(c["watsat"][k], c["h2osoi_ice"][k] / (c["dz"][k] * DENICE))
            c["icefrac"][k] = min(1, vol_ice / c["watsat"][k])
        c["rsub_top_layer"] = c["watsat"] * 0
        dtime = c["dtime"]
        jwt = self.getjwt()
        q_perch_max = PERCHED_BASEFLOW_SCALAR * math.sin(c["topo_slope"] * (math.pi / 180))
        t = c["t_soisno"]
        k_frz = nlevsoi if t[0] > TFRZ else 1
        for k in range(2, nlevsoi + 1):
            if t[k - 2] > TFRZ and t[k - 1] <= TFRZ:
                k_frz = k
                break
        frost_table = c["z"][k_frz - 1]
        rsub_top = 0.0

        def imped_k(k):
            return 10 ** (-e_ice * (0.5 * (c["icefrac"][k - 1] + c["icefrac"][min(nlevsoi, k + 1) - 1])))

        def remove_perched(krange, rsub_top_tot, zwt_key):
            zwt = c["zwt"] if zwt_key else self._zwt
            for k in krange:
                layer = min(max(rsub_top_tot, -(c["h2osoi_liq"][k - 1] - WATMIN)), 0)
                c["rsub_top_layer"][k - 1] = layer
                rsub_top_tot -= layer
                c["h2osoi_liq"][k - 1] += layer
                if rsub_top_tot >= 0:
                    zwt = zwt - layer / c["eff_porosity"][k - 1] / 1000
                    break
                zwt = c["zi"][k - 1]
            if zwt_key:
                c["zwt"] = zwt
            return rsub_top_tot

        if c["zwt"] < frost_table and t[k_frz - 1] <= TFRZ:
            wtsub = q_perch = 0.0
            for k in range(jwt + 1, k_frz + 1):
                q_perch += imped_k(k) * c["hksat"][k - 1] * c["dzmm"][k - 1]
                wtsub += c["dzmm"][k - 1]
            if wtsub > 0:
                q_perch = q_perch / wtsub
            drain_perched = q_perch_max * q_perch * (frost_table - c["zwt"])
            remove_perched(range(jwt + 1, k_frz + 1), -drain_perched * dtime, True)
            jwt = self.getjwt()
        else:
            k_perch = self._perched_level(k_frz)
            if k_frz > k_perch:
                self._zwt = self._zwt_perched(k_perch)
                wtsub = q_perch = 0.0
                for k in range(k_perch, k_frz + 1):
                    q_perch += imped_k(k) * c["hksat"][k - 1] * c["dzmm"][k - 1]
                    wtsub += c["dzmm"][k - 1]
                if wtsub > 0:
                    q_perch = q_perch / wtsub
                drain_perched = q_perch_max * q_perch * (frost_table - self._zwt)
                remove_perched(range(k_perch + 1, k_frz + 1), -drain_perched * dtime, False)
            # topographic runoff
            fff = 1 / HKDEPTH
            dzsum = icefracsum = 0.0
            for j in range(max(jwt, 1), nlevsoi + 1):
                dzsum += c["dzmm"][j - 1]
                icefracsum += c["icefrac"][j - 1] * c["dzmm"][j - 1]
            imped = 10 ** (-e_ice * (icefracsum / dzsum))
            rsub_top_max = 10 * math.sin((math.pi / 180) * c["topo_slope"])
            rsub_top = imped * rsub_top_max * math.exp(-fff * c["zwt"])
            rous = max(self._specific_yield(nlevsoi - 1), self.aq_sp_yield_min)
            if jwt == nlevsoi:
                c["wa"] = c["wa"] - rsub_top * dtime
                c["zwt"] = c["zwt"] + (rsub_top * dtime) / 1000 / rous
                c["h2osoi_liq"][nlevsoi - 1] += max(0, c["wa"] - AQUIFER_WATER_BASELINE)
                c["wa"] = min(c["wa"], AQUIFER_WATER_BASELINE)
            else:
                rsub_top_tot = -rsub_top * dtime
                if rsub_top_tot > 0:
                    print("(soilwatermovement/Drainage) RSUB_TOP IS POSITIVE in Drainage!")
                else:
                    for j in range(jwt + 1, nlevsoi + 1):
                        s_y = self._specific_yield(j - 1)
                        layer = min(max(rsub_top_tot, -(s_y * (c["zi"][j - 1] - c["zwt"]) * 1e3)), 0)
                        c["rsub_top_layer"][j - 1] = layer
                        c["h2osoi_liq"][j - 1] += layer
                        rsub_top_tot -= layer
                        if rsub_top_tot >= 0:
                            c["zwt"] = c["zwt"] - layer / s_y / 1000
                            break
                        c["zwt"] = c["zi"][j - 1]
                    c["zwt"] = c["zwt"] - rsub_top_tot / 1000 / rous
                    c["wa"] = c["wa"] + rsub_top_tot
                jwt = self.getjwt()
            c["zwt"] = min(80, max(0, c["zwt"]))
        # excess water above saturation bubbles up to the layer above
        c["excess"] = 0.0
        xsi_arr = _grow(np.array(c.get("xsi", np.zeros(nlevsoi)), float), nlevsoi)
        liq, eff, dzmm = c["h2osoi_liq"], c["eff_porosity"], c["dzmm"]
        for j in range(nlevsoi, 1, -1):
            k = j - 1
            xsi = max(liq[k] - eff[k] * dzmm[k], 0)
            liq[k] = min(eff[k] * dzmm[k], liq[k])
            xsi_arr[k] = xsi
            if j == jdrain:
                c["excess"] = xsi
                continue
            liq[k - 1] += xsi
        c["xsi"] = xsi_arr
        xs1 = max(max(liq[0] - WATMIN, 0) - max(0, (PONDMX + c["watsat"][0] * dzmm[0] - c["h2osoi_ice"][0] - WATMIN)), 0)
        liq[0] -= xs1
        qflx_rsub_sat = xs1 / dtime
        c["xs1_liq"] = xs1
        xs1 = max(max(c["h2osoi_ice"][0], 0) - max(0, (PONDMX + c["watsat"][0] * dzmm[0] - liq[0])), 0)
        c["h2osoi_ice"][0] = min(max(0, PONDMX + c["watsat"][0] * dzmm[0] - liq[0]), c["h2osoi_ice"][0])
        c["xs1_ice"] = xs1
        for j in range(1, nlevsoi):
            k = j - 1
            if liq[k] < WATMIN:
                xs = WATMIN - liq[k]
                if j == jwt:
                    c["zwt"] = c["zwt"] + xs / eff[k] / 1000
            else:
                xs = 0.0
            liq[k] += xs
            liq[k + 1] -= xs
        k = nlevsoi - 1
        if liq[k] < WATMIN:
            xs = WATMIN - liq[k]
            for i in range(nlevsoi - 1, 0, -1):
                available = max(liq[i - 1] - WATMIN - xs, 0)
                if available >= xs:
                    liq[k] += xs
                    liq[i - 1] -= xs
                    xs = 0.0
                    break
                liq[k] += available
                liq[i - 1] -= available
                xs -= available
        else:
            xs = 0.0
        liq[k] += xs
        rsub_top = rsub_top - xs / dtime
        c["qflx_drain"] = qflx_rsub_sat + rsub_top

    # ------------------------------------------------------------------- main
    def run(self):
        c, nlevsoi = self.c, self.nlevsoi
        jdrain = 0
        if c["draindepthm"] > 0:
            wdrain = np.flatnonzero(c["zi"] >= c["draindepthm"])
            if len(wdrain) > 0:
                jdrain = wdrain[0] + 1
        c["jdrain"] = jdrain
        c["excess"] = 0.0
        self.set_soil_water_fractions()
        c["vwc_liq"] = np.maximum(c["h2osoi_liq"], 1.0e-6) / (c["dz"] * DENH2O)
        c["zwtmm"] = c["zwt"] * 1000
        c["qin"], c["qout"] = c["z"] * 0, c["z"] * 0
        if self.addice:
            self.phasechange()
        self.soilwater_zengdecker2009()
        self.water_table()
        self.drainage(jdrain)
        nc = c["ncells"]
        c["h2oday"] = c["h2osoi_liq"][:nc].copy()
        c["iceday"] = c["h2osoi_ice"][:nc].copy()
        c["h2ototday"] = c["h2osoi_liq"][:nc] + c["h2osoi_ice"][:nc]
        c["qinday"] = c["qin"][:nc] * c["dtime"]
        c["qoutday"] = c["qout"][:nc] * c["dtime"]
        c["qflx_rootsoi_day"] = c["qflx_rootsoi_col"][:nc] * c["dtime"]
        c["zwtcm"] = c["zwt"] * 100
        c["aqh2omm"] = c["wa"]
        return c


def soil_water_movement(ds, c, idx, addice=1, etscale=1.0):
    """soilwatermovement.m for day ``idx`` (1-based); returns the updated state."""
    return _SoilWaterMovement(ds, c, idx, addice, etscale).run()


_DAY_FIELDS = {  # ds.run field: c field, saved for every day by getrichquick
    "qinday": "qinday", "qoutday": "qoutday", "qflx_rootsoi_day": "qflx_rootsoi_day", "ice": "iceday",
    "H2OSOI": "h2ototday", "h2osoi_liq0": "h2osoi_liq0", "h2osoi_liq": "h2osoi_liq",
    "h2osoi_ice": "h2osoi_ice", "eff_porosity": "eff_porosity", "xsi": "xsi", "hk": "hk", "smp": "smp",
    "imped": "imped", "rsub_top_layer": "rsub_top_layer", "vwc_liq": "vwc_liq", "den": "den", "num": "num",
    "dwat2": "dwat2",
}


def get_rich_quick(ds, addice=1):
    """getrichquick.m: run soilwatermovement day by day and rebuild the cell fluxes.

    The tile-drain code is disabled in the MATLAB (``if (0 & ...)``), so
    ``ds.run.qdrain`` is all zeros and ``ds.parm.draindepthm`` is set to -1.
    """
    run, cell, parm = ds["run"], ds["cell"], ds["parm"]
    n, nc = ds["nrun"], ds["ncells"]
    run["qdrain"] = np.zeros((nc, n))
    parm["draindepthm"] = -1
    vsoil = np.asarray(cell["vsoil"], float)
    out = {}
    for k in ("aquifer", "qflx_drain_day", "qcharge", "watertabledepth", "excess"):
        run[k] = np.zeros(n)
    run.setdefault("Sws", np.zeros((nc, n)))
    run.setdefault("Swstot", np.zeros((nc, n)))
    c = None  # getrichquick's own c is not recognised by soilwatermovement
    for i in range(1, n + 1):
        c = soil_water_movement(ds, c, i, addice=addice)
        for k, ck in _DAY_FIELDS.items():
            v = np.asarray(c[ck], float)
            if k not in out:
                out[k] = np.zeros((len(v), n))
            elif len(v) > out[k].shape[0]:
                out[k] = np.vstack([out[k], np.zeros((len(v) - out[k].shape[0], n))])
            out[k][: len(v), i - 1] = v
        run["aquifer"][i - 1] = c["aqh2omm"]
        run["qflx_drain_day"][i - 1] = c["qflx_drain"] * 86400
        run["qcharge"][i - 1] = c["qcharge"]
        run["watertabledepth"][i - 1] = c["zwtcm"]
        run["excess"][i - 1] = c["excess"]
        if i == n:
            continue
        run["targH2O"][:, i] = c["h2oday"]
        run["Sws"][:, i] = c["h2oday"] / vsoil
        run["Swstot"][:, i] = c["h2ototday"] / vsoil
    run.update(out)
    run["airpor"] = np.asarray(cell["por"], float)[:, None] - run["Sws"]
    infil = np.asarray(run["infil"], float)
    run["inperc"] = np.zeros((nc, n))
    run["outperc"] = np.zeros((nc, n))
    run["inperc"][0] = infil - run["qflx_drain_day"]
    liq, ice = run["h2osoi_liq"], run["h2osoi_ice"]
    dliq = liq - np.roll(liq, 1, axis=1)
    dliq[:, 0] = 0
    dice = ice - np.roll(ice, 1, axis=1)
    dice[:, 0] = 0
    cellET = run["cellET"]
    bal = (run["inperc"][0] - dice[0] - cellET[0]) - dliq[0]
    run["outperc"][0] = bal
    for j in range(1, nc):
        prev = bal
        bal = (prev - dice[j] - cellET[j]) - dliq[j]
        run["outperc"][j] = bal
        run["inperc"][j] = prev
    return ds
