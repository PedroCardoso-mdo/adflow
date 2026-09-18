#!/usr/bin/env python
"""Exact Falkner-Skan map lambda_thetaL(eta*) -> lambda_theta, represented analytically for the code.

From falkner_skan_gain.json (similarity solutions): for each beta, u = lamL_raw(eta*) + off (off =
Blasius-neutral offset) and v = lambda_theta(beta). The map is one-to-one and nonlinear:
  adverse   (u < 0): saturating  ->  v = -A tanh(-u/B)
  favourable(u > 0): convex      ->  v =  C (exp(u/D) - 1)
The four constants are least-squares fits to the SIMILARITY SOLUTION (no CFD or target data), the two
branches blended with the same smooth switch used in the code (0.5(1+tanh(u/eps))).
Prints the constants, the residual of the representation, and writes falkner_skan_map.json."""
import json, os, numpy as np
from scipy.optimize import least_squares
H = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(H, "falkner_skan_gain.json"))); off = d["off_calibrated"]
u = np.array([r["lamL_raw"] + off for r in d["rows"]]); v = np.array([r["lam_theta"] for r in d["rows"]])
adv = u < -1e-9; fav = u > 1e-9
fa = least_squares(lambda p: -p[0] * np.tanh(-u[adv] / p[1]) - v[adv], [0.07, 0.03])
ff = least_squares(lambda p: p[0] * (np.exp(u[fav] / p[1]) - 1) - v[fav], [0.13, 0.05])
A, B = fa.x; C, D = ff.x
def gmap(x, eps=0.003):
    va = -A * np.tanh(-x / B); vf = C * (np.exp(x / D) - 1); w = 0.5 * (1 + np.tanh(x / eps))
    return (1 - w) * va + w * vf
print(f"adverse:   lambda_theta = -A tanh(-u/B),   A = {A:.5f}  B = {B:.5f}   (slope at 0: {A/B:.3f})")
print(f"favourable: lambda_theta = C (exp(u/D) - 1), C = {C:.5f}  D = {D:.5f}   (slope at 0: {C/D:.3f})")
res = gmap(u) - v
print(f"representation error: max |dv| = {np.abs(res).max():.4f}, rms = {np.sqrt(np.mean(res**2)):.4f} (lambda units); linear gain 2: max |dv| = {np.abs(2*u - v).max():.4f}")
print(f"{'beta':>8} {'u=lamL+off':>11} {'lam_theta':>10} {'map':>8} {'2u':>8}")
for r, uu, vv in zip(d["rows"], u, v):
    if r["beta"] < 0 or abs(r["beta"] * 10 - round(r["beta"] * 10)) < 1e-6:
        print(f"{r['beta']:8.4f} {uu:11.4f} {vv:10.4f} {gmap(uu):8.4f} {2*uu:8.4f}")
json.dump(dict(off=off, A=float(A), B=float(B), C=float(C), D=float(D), eps=0.003,
               note="lambda_theta = blend(-A tanh(-u/B) [u<0], C (exp(u/D)-1) [u>0]), u = lamL_raw + off; fit to the Falkner-Skan similarity map only"),
          open(os.path.join(H, "falkner_skan_map.json"), "w"), indent=1)
