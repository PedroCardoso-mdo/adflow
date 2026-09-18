#!/usr/bin/env python
"""Falkner-Skan calibration of the SA-BCM pressure-gradient sensor.

For each wedge parameter beta (U_e = C x^m, m = beta/(2-beta)) the similarity
profile f(eta) is solved and evaluated at eta* = argmax(eta^2 f''), i.e. the
wall-normal location of the Re_v = rho d^2 |du/dy| / mu maximum where the BCM
Term1 fires. There:

  lambda_theta (exact, LM2009 definition)  = theta^2/nu dU_e/dx = 2 m thetahat^2/(m+1)
  lambda_thetaL (Menter 2015, offset = 0)  = -c (d^2/nu) dV/dy   with dV/dy = -du/dx
      (d^2/nu) du/dx = (2 eta^2/(m+1)) [ m f' - (1-m)/2 eta f'' ]
  =>  lambda_thetaL_raw = c (2 eta^2/(m+1)) [ m f' - (1-m)/2 eta f'' ],  c = 7.57e-3

Outputs: the offset that zeroes lambda_thetaL for Blasius at eta*, the linear
gain lambda_theta ~ gain * (lambda_thetaL_raw + off) fitted over the LM2009
validity window |lambda_theta| <= 0.1, a table and a figure.
"""
import json, sys, os
import numpy as np
from scipy.integrate import solve_bvp, trapezoid

C_MENTER = 7.57e-3
OFF_MENTER = 0.0128
ETA_MAX = 14.0

def falkner_skan(beta, guess=None):
    def ode(eta, y):
        f, fp, fpp = y
        return np.vstack([fp, fpp, -f * fpp - beta * (1.0 - fp**2)])
    def bc(ya, yb):
        return np.array([ya[0], ya[1], yb[1] - 1.0])
    eta = np.linspace(0.0, ETA_MAX, 1200)
    if guess is None:
        y0 = np.vstack([eta - 1.7 * (1 - np.exp(-eta / 1.7)), 1 - np.exp(-eta / 1.7), np.exp(-eta / 1.7) / 1.7])
    else:
        y0 = guess.sol(eta)
    sol = solve_bvp(ode, bc, eta, y0, tol=1e-8, max_nodes=200000)
    assert sol.success, f"FS failed at beta={beta}: {sol.message}"
    return sol

def evaluate(beta, sol):
    m = beta / (2.0 - beta)
    eta = np.linspace(0.0, ETA_MAX, 20001)
    f, fp, fpp = sol.sol(eta)
    thetahat = trapezoid(fp * (1.0 - fp), eta)
    lam_theta = 2.0 * m * thetahat**2 / (m + 1.0)
    rev = eta**2 * np.abs(fpp)                      # ~ Re_v (BL approx: omega = du/dy)
    i = int(np.argmax(rev)); eta_s = eta[i]
    lamL_raw = C_MENTER * (2.0 * eta_s**2 / (m + 1.0)) * (m * fp[i] - 0.5 * (1.0 - m) * eta_s * fpp[i])
    return dict(beta=beta, m=m, thetahat=thetahat, lam_theta=lam_theta, eta_star=eta_s,
                fpp0=float(fpp[0]), lamL_raw=lamL_raw)

def main(out_dir):
    betas_pos = np.concatenate([np.linspace(0.0, 0.5, 26), np.linspace(0.55, 1.0, 10)])
    betas_neg = np.linspace(0.0, -0.1988, 26)[1:]
    rows = []
    sol = None
    for b in betas_pos:
        sol = falkner_skan(b, sol); rows.append(evaluate(b, sol))
    sol = None
    for b in betas_neg:
        sol = falkner_skan(b, sol); rows.append(evaluate(b, sol))
    rows.sort(key=lambda r: r["beta"])
    blas = [r for r in rows if r["beta"] == 0.0][0]
    off = -blas["lamL_raw"]
    # fit gain on the LM2009 window |lambda_theta| <= 0.1 with the corrected offset
    lt = np.array([r["lam_theta"] for r in rows]); ll = np.array([r["lamL_raw"] + off for r in rows])
    w = np.abs(lt) <= 0.1
    gain = float(np.sum(lt[w] * ll[w]) / np.sum(ll[w]**2))
    gain_neg = float(np.sum(lt[w & (lt < 0)] * ll[w & (lt < 0)]) / np.sum(ll[w & (lt < 0)]**2))
    gain_pos = float(np.sum(lt[w & (lt > 0)] * ll[w & (lt > 0)]) / np.sum(ll[w & (lt > 0)]**2))
    # what Menter's own offset gives for Blasius
    res = dict(c=C_MENTER, off_menter=OFF_MENTER, blasius=dict(eta_star=blas["eta_star"], fpp0=blas["fpp0"],
               lamL_raw=blas["lamL_raw"], lamL_with_menter_offset=blas["lamL_raw"] + OFF_MENTER),
               off_calibrated=off, gain=gain, gain_adverse=gain_neg, gain_favorable=gain_pos,
               fit_window="|lambda_theta| <= 0.1", rows=rows)
    os.makedirs(out_dir, exist_ok=True)
    json.dump(res, open(os.path.join(out_dir, "falkner_skan_gain.json"), "w"), indent=1)
    print(f"Blasius: eta* = {blas['eta_star']:.3f}  f''(0) = {blas['fpp0']:.5f}  lamL_raw = {blas['lamL_raw']:+.5f}"
          f"  -> with Menter offset 0.0128: {blas['lamL_raw']+OFF_MENTER:+.5f}")
    print(f"calibrated offset (lamL=0 at Blasius eta*): {off:.5f}")
    print(f"gain (|lam|<=0.1): {gain:.4f}   adverse-only: {gain_neg:.4f}   favorable-only: {gain_pos:.4f}")
    print(f"{'beta':>8} {'m':>8} {'eta*':>6} {'lam_theta':>10} {'lamL_raw':>10} {'lamL+off':>10} {'gain*(..)':>10}")
    for r in rows:
        if abs(r["beta"] * 20 - round(r["beta"] * 20)) < 1e-9 or r["beta"] < 0:
            print(f"{r['beta']:8.4f} {r['m']:8.4f} {r['eta_star']:6.2f} {r['lam_theta']:10.4f} {r['lamL_raw']:10.4f} "
                  f"{r['lamL_raw']+off:10.4f} {gain*(r['lamL_raw']+off):10.4f}")
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
        ax[0].plot([r["beta"] for r in rows], lt, "k-", label=r"$\lambda_\theta$ (exact)")
        ax[0].plot([r["beta"] for r in rows], gain * ll, "r--", label=fr"gain$\cdot(\lambda_{{\theta L}}+$off$)$, gain={gain:.2f}")
        ax[0].plot([r["beta"] for r in rows], [r["lamL_raw"] + OFF_MENTER for r in rows], "b:", label=r"$\lambda_{\theta L}$ Menter (off 0.0128)")
        ax[0].axhline(0.1, color="0.7", lw=0.8); ax[0].axhline(-0.1, color="0.7", lw=0.8)
        ax[0].set_xlabel(r"$\beta$"); ax[0].set_ylabel(r"$\lambda$"); ax[0].legend(fontsize=7)
        ax[1].plot(ll, lt, "ko", ms=3); xx = np.linspace(ll.min(), ll.max(), 50); ax[1].plot(xx, gain * xx, "r--")
        ax[1].set_xlabel(r"$\lambda_{\theta L}$ + off (at $\eta^*$)"); ax[1].set_ylabel(r"$\lambda_\theta$")
        fig.tight_layout(); fig.savefig(os.path.join(out_dir, "falkner_skan_gain.png"), dpi=150)
    except Exception as e:
        print("plot skipped:", e)

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__)))
