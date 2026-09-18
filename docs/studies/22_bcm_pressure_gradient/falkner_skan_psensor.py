#!/usr/bin/env python
"""Pressure-gradient sensor calibrated on Falkner-Skan similarity: lambda_p = K (d^2/nu) dU_e/ds evaluated at
the Re_v peak d = eta* delta_s, with K = (thetahat0/eta0*)^2 from Blasius. Exact lambda_theta = thetahat^2 delta_s^2 U'/nu.
Ratio r(beta) = lambda_theta/lambda_p = (thetahat(beta)/eta*(beta))^2 / K. No offset needed (U' = 0 for Blasius)."""
import json, os, numpy as np
from falkner_skan_gain import falkner_skan, evaluate
H = os.path.dirname(os.path.abspath(__file__))
rows = []; sol = None
for b in np.linspace(0, 1.0, 21): sol = falkner_skan(b, sol); rows.append(evaluate(b, sol))
sol = None
for b in np.linspace(0, -0.1988, 20)[1:]: sol = falkner_skan(b, sol); rows.append(evaluate(b, sol))
rows.sort(key=lambda r: r["beta"]); b0 = [r for r in rows if r["beta"] == 0][0]
K0 = (b0["thetahat"] / b0["eta_star"])**2
print(f"Blasius: thetahat = {b0['thetahat']:.5f}, eta* = {b0['eta_star']:.4f}  ->  K0 = (thetahat/eta*)^2 = {K0:.5f}")
print(f"{'beta':>8} {'m':>8} {'lam_theta':>10} {'lam_p(K0)':>10} {'ratio':>7}")
out = []
for r in rows:
    m = r["m"]; lt = r["lam_theta"]
    lam_p = K0 * r["eta_star"]**2 * 2 * m / (m + 1)   # (eta* delta_s)^2 U'/nu with delta_s^2 U'/nu = 2m/(m+1)
    out.append(dict(beta=r["beta"], lam_theta=lt, lam_p=lam_p, ratio=(lt / lam_p if abs(lam_p) > 1e-12 else 1.0)))
    if r["beta"] < 0 or abs(r["beta"] * 10 - round(r["beta"] * 10)) < 1e-6:
        print(f"{r['beta']:8.4f} {m:8.4f} {lt:10.4f} {lam_p:10.4f} {out[-1]['ratio']:7.3f}")
json.dump(dict(K0=K0, rows=out), open(os.path.join(H, "falkner_skan_psensor.json"), "w"), indent=1)
