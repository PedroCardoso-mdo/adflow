#!/usr/bin/env python
"""Effective-threshold calibration of the SA-BCM pressure-gradient sensor.

The BCM fires in the first cell where Re_v/(Re_thetac F(lambda)) > 1. With the sensor, F varies
across the boundary layer (lambda_thetaL(eta) = c (2 eta^2/(m+1)) [m f' - (1-m)/2 eta f''] + off,
before gain), so what matters is not F at the Re_v peak eta* but

    G_eff(beta) = max_eta [ eta^2 f''(eta) / F(gain*(lamL_raw(eta) + off)) ] / max_eta [ eta^2 f''(eta) ]

i.e. the effective threshold factor is 1/G_eff. The reference model applies F(lambda_theta) with the
exact (edge) lambda_theta of the profile. Calibrate (off, gain) so that 1/G_eff(beta) ~ F(lambda_theta(beta))
over the LM2009 window, with the Blasius case neutral (G_eff(0) = 1) as a hard constraint.
Tu enters F: calibrate at Tu = 0.1, 0.5 % and report both.
"""
import json, os, sys
import numpy as np
from scipy.optimize import minimize_scalar, minimize
from falkner_skan_gain import falkner_skan, C_MENTER, ETA_MAX

def F_langtry(lam, Tu):
    lam = np.clip(lam, -0.1, 0.1)
    F1 = 1 + 0.275 * (1 - np.exp(-35 * lam)) * np.exp(-Tu / 0.5)
    F3 = 1 - (-12.986 * lam - 123.66 * lam**2 - 405.689 * lam**3) * np.exp(-(Tu / 1.5)**1.5)
    return np.minimum(np.maximum(F1, 1.0), F3)

def profile(beta, sol):
    m = beta / (2 - beta); eta = np.linspace(0, ETA_MAX, 8001); f, fp, fpp = sol.sol(eta)
    rev = eta**2 * np.abs(fpp)
    lamL_raw = C_MENTER * (2 * eta**2 / (m + 1)) * (m * fp - 0.5 * (1 - m) * eta * fpp)
    thetahat = np.trapz(fp * (1 - fp), eta); lam_theta = 2 * m * thetahat**2 / (m + 1)
    return eta, rev, lamL_raw, lam_theta

def G_eff(rev, lamL_raw, off, gain, Tu):
    return np.max(rev / F_langtry(gain * (lamL_raw + off), Tu)) / np.max(rev)

def main():
    betas = np.concatenate([np.linspace(-0.1988, 0, 20), np.linspace(0.05, 1.0, 20)])
    profs = {}; sol = None
    for b in np.linspace(0, 1.0, 21): sol = falkner_skan(b, sol); profs[round(b, 4)] = profile(b, sol)
    sol = None
    for b in np.linspace(0, -0.1988, 20): sol = falkner_skan(b, sol); profs[round(b, 4)] = profile(b, sol)
    out = {}
    for Tu in (0.1, 0.5):
        eta0, rev0, lamL0, _ = profs[0.0]
        # (1) Blasius-neutral offset for a given gain: G_eff(0) = 1  (F >= 1 everywhere along the profile)
        def off_neutral(gain):
            # the smallest offset such that F(gain*(lamL+off)) >= 1 where rev/max(rev) is within the top; require G=1
            lo, hi = 0.0, 0.2
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if G_eff(rev0, lamL0, mid, gain, Tu) > 1.0 + 1e-9: lo = mid
                else: hi = mid
            return hi
        # (2) fit gain so that 1/G_eff(beta) matches F(lambda_theta) for adverse betas (the regime of interest)
        def cost(gain):
            off = off_neutral(gain); c = 0.0
            for b, (eta, rev, lamL, lt) in profs.items():
                if b >= 0 or lt < -0.1: continue
                c += (1.0 / G_eff(rev, lamL, off, gain, Tu) - F_langtry(lt, Tu))**2
            return c
        res = minimize_scalar(cost, bounds=(0.2, 6.0), method="bounded")
        gain = float(res.x); off = float(off_neutral(gain))
        rows = []
        for b in sorted(profs):
            eta, rev, lamL, lt = profs[b]
            rows.append(dict(beta=b, lam_theta=lt, F_ref=float(F_langtry(lt, Tu)), F_eff=float(1.0 / G_eff(rev, lamL, off, gain, Tu)),
                             F_eff_menter=float(1.0 / G_eff(rev, lamL, 0.0128, 1.0, Tu)),
                             F_eff_first=float(1.0 / G_eff(rev, lamL, 0.01623, 2.0, Tu))))
        out[f"Tu{Tu}"] = dict(gain=gain, off=off, rows=rows)
        print(f"\n=== Tu = {Tu} %  -> calibrated gain = {gain:.3f}, offset = {off:.5f}")
        print(f"{'beta':>8} {'lam_th':>8} {'F_ref':>7} {'F_eff':>7} | {'Menter(0.0128,1)':>16} {'first(0.01623,2)':>16}")
        for r in rows:
            if abs(r['beta'] * 10 - round(r['beta'] * 10)) < 1e-6 or r['beta'] < 0:
                print(f"{r['beta']:8.4f} {r['lam_theta']:8.4f} {r['F_ref']:7.3f} {r['F_eff']:7.3f} | {r['F_eff_menter']:16.3f} {r['F_eff_first']:16.3f}")
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "falkner_skan_effective.json"), "w"), indent=1)

if __name__ == "__main__":
    main()
