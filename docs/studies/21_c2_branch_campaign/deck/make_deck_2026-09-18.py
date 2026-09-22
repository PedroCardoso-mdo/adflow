# Status deck: SA-BCM 2D optimisation reruns for the article + the branch problem (2026-09-14)
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
FIG = "../figures/cf_tu_c2c3c4.png"; FIG_C1 = "../figures/c1_meshes_cp_profile.png"; FIG_MESH = "../figures/mesh_effect_cp_cf_tu05.png"
prs = Presentation(); prs.slide_width = Inches(13.33); prs.slide_height = Inches(7.5)
BL = prs.slide_layouts[6]
def title(s, t, sub=None):
    tb = s.shapes.add_textbox(Inches(0.4), Inches(0.25), Inches(12.5), Inches(0.9)); p = tb.text_frame.paragraphs[0]; p.text = t; p.font.size = Pt(28); p.font.bold = True
    if sub:
        p2 = tb.text_frame.add_paragraph(); p2.text = sub; p2.font.size = Pt(14); p2.font.color.rgb = RGBColor(90, 90, 90)
def bullets(s, items, x=0.5, y=1.3, w=12.3, h=5.8, size=16):
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h)); tf = tb.text_frame; tf.word_wrap = True
    for i, it in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        lvl = 0
        while it.startswith("  "): it = it[2:]; lvl += 1
        p.text = ("• " if lvl == 0 else "– ") + it; p.level = lvl; p.font.size = Pt(size - 2 * lvl); p.space_after = Pt(4)
def table(s, rows, x=0.5, y=1.3, w=12.3, colw=None, size=12):
    n, m = len(rows), len(rows[0])
    t = s.shapes.add_table(n, m, Inches(x), Inches(y), Inches(w), Inches(0.35 * n)).table
    for i, r in enumerate(rows):
        for j, v in enumerate(r):
            c = t.cell(i, j); c.text = str(v)
            for p in c.text_frame.paragraphs: p.font.size = Pt(size); p.font.bold = (i == 0)
    if colw:
        for j, cw in enumerate(colw): t.columns[j].width = Inches(cw)
    return t
def note(s, txt, y=6.7):
    tb = s.shapes.add_textbox(Inches(0.5), Inches(y), Inches(12.3), Inches(0.6)); p = tb.text_frame.paragraphs[0]; p.text = txt; p.font.size = Pt(11); p.font.italic = True; tb.text_frame.word_wrap = True

# 1 title
s = prs.slides.add_slide(BL)
title(s, "SA-BCM 2D optimisation — reruns for the article and the branch problem", "Status deck, 2026-09-14, updated 2026-09-18 (dense mesh + cross-evaluation of the optima + pressure-gradient sensor, slides 22–27) · Deucalion SA_BCM_ARTICLE/2d_opt · repo docs/studies/21_c2_branch_campaign")
bullets(s, ["Question: C2 (NACA0012 → SA-BCM, cl = 0.3, Tu = 0.5 %) gave −18 % with the buggy code, −4…−6 % with the corrected one; C3 gives −24 %.",
            "Constraint set by the user: keep the C2 definition (NACA start, SA-BCM only, article constraints, fs ≤ 0.15); only the optimizer/solver strategy may change.",
            "All SA-BCM numbers here: cold start (fresh process, 16 ranks), trimmed to cl = 0.3, ev = 1e-7, fs = 0.08. SA numbers: run2d, 32 ranks.",
            "Contents: diagnosis (why the optimizer stalls) → every strategy tried → results per Tu → the branch mechanism → proposals."], y=1.6)

# 2 baseline not unique
s = prs.slides.add_slide(BL); title(s, "1. The baseline itself is not unique", "NACA0012, cl = 0.3, Tu = 0.5 %, same code, same mesh")
table(s, [["Solver path", "cd", "note"],
          ["cold, 16 or 32 ranks", "0.006149", "bit-for-bit reproducible at fixed rank count"],
          ["cold, 48 ranks", "0.006202", "+0.9 %"], ["cold, 64 ranks (opt_c2)", "0.006456", "+5 %"],
          ["trim recipe (ANK CFL 3e5)", "0.008061", "+31 %"], ["ev = 1e-10", "0.005197", "−15 % (different farfield BC)"],
          ["article C2 geometry, today's code", "0.0064", "worse than the NACA — the old −18 % came from the bug"]], colw=[3.5, 1.5, 7.3])
note(s, "The partition/pseudo-transient path selects the converged state. Deterministic for a fixed rank count.")

# 3 line test
s = prs.slides.add_slide(BL); title(s, "2. Line test NACA → C3 in FFD space: no barrier, but discrete jumps", "11 points + fine scan t = 0.80–1.00; adjoint directional derivative at every point")
table(s, [["t", "cd (fs 0.08)", "cd (fs 0.15)", "adjoint dcd/dt", "x_tr upper"],
          ["0.0 (NACA)", "0.006309", "0.006309", "−2.1e-4", "0.317"], ["0.5", "0.005765", "0.005766", "−1.5e-4", "0.352"],
          ["0.80–0.87", "0.005638 → 0.005634 (ramp −6.7e-5/unit)", "same", "−7.3e-5 … −6.1e-5 ✔", "0.369"],
          ["0.88–0.96", "0.005733 / 0.005631 alternating", "same", "−4e-5", "0.352 / 0.369"],
          ["0.97–1.00 (C3)", "0.004703 → 0.004695", "0.004977", "−2.3e-4", "0.550"]], colw=[1.6, 3.6, 1.8, 2.3, 1.8])
bullets(s, ["Within a plateau the adjoint is exact (measured slope = adjoint). The drop is carried by branch switches: ±1.8 % (transition moving one cell) and −16 % (x_tr 0.37 → 0.55).",
            "fs = 0.15 gives identical plateaus and jumps (and +5.8 % on the good branch) → smoothing is not the lever."], y=4.0, size=14)

# 4 hysteresis + Tu ramp
s = prs.slides.add_slide(BL); title(s, "3. Hysteresis chains and Tu-ramp probe: several steady states per geometry")
table(s, [["t", "cold", "warm from NACA (fwd)", "warm from C3 (bwd)", "Tu₀ = 0.1 % → 0.5 %", "Tu₀ = 3 % → 0.5 %"],
          ["0.0", "0.006149", "0.006149", "0.006255", "0.006248", "0.009386"], ["0.5", "0.005765", "0.006069", "0.005963", "0.005863", "0.009387"],
          ["0.8", "0.005639", "0.006047", "0.005639", "0.005738", "0.009446"], ["0.96", "0.005630", "0.006045", "0.005244", "0.005244", "0.009459"],
          ["1.0 (C3)", "0.004695", "0.006045", "0.004643", "0.004643", "0.009464"]], colw=[1.2, 1.6, 2.4, 2.4, 2.4, 2.3])
bullets(s, ["Warm from the NACA carries the tripped branch all the way to the C3 shape (0.00605 instead of 0.00464) — a classic warm-started optimization can never see the C3 gain.",
            "Starting turbulent (Tu₀ = 3 %) → fully turbulent state (0.0094 ≈ SA) is stable at Tu = 0.5 % for every geometry: a spurious self-sustained mode.",
            "The good branch (x_tr 0.55) is a narrow basin around C3 (collapses within Δt = 0.04)."], y=3.9, size=14)

# 5 mechanism
s = prs.slides.add_slide(BL); title(s, "4. The three stable modes and why the gradient cannot see the gain")
bullets(s, ["Mode 1 — turbulent from the leading edge (cd 0.0094): the ν_t/ν term keeps γ = 1, which keeps ν_t high. Self-sustained; only reached by starting turbulent. Spurious.",
            "Mode 2 — transition mid-chord (0.32–0.37 c) in pairs of neighbouring cells: γ switches where Re_v crosses the Tu threshold; on a flat pressure region Re_v grows slowly, so the crossing is ill-conditioned; transition at cell i changes the downstream boundary layer and the pressure, which moves the crossing — two adjacent cells are both self-consistent (±1.8 % cd). Discrete because of the mesh; fsmooth does not remove it.",
            "Mode 3 — long laminar run (x_tr 0.55 c): only for shapes with a pressure plateau (C3 family). Narrow basin.",
            "SLSQP: the adjoint measures the slope inside a step (correct, tiny); the line search samples 1–5 points on one line with shrinking steps; a +1.8 % landing at every trial ⇒ 'converged'. The accepted point can also switch branch (objective drifts up).",
            "A denser mesh makes the steps smaller (more, lower equilibria); the ill-conditioning of the crossing on a flat pressure plateau remains."], size=15)

# 6 strategies
s = prs.slides.add_slide(BL); title(s, "5. Optimizer strategies tried (NACA0012 start, Tu = 0.5 %, fs 0.08)", "best cd cold-verified; Δ vs baseline 0.006149; C4 = −11.4 %, C3 = −24.5 %")
table(s, [["Strategy", "Idea", "cd", "Δ", "Verdict"],
          ["SLSQP as in the article (scale 10)", "—", "0.00617 / 0.00583 (fs 0.04)", "−4.4 / −5.2 %", "stalls on branch noise"],
          ["B1/B1′ warm start + guard", "continuity", "0.00756 / 0.00631", "+17 / +2.6 %", "hysteresis poisons the base point"],
          ["B cold, DV scale 2", "smaller steps", "0.00571 (not reproducible)", "−7 % then drifts", "—"],
          ["B4 / B5 cold + warm-from-best; laminar-first", "branch selection", "0.00591 / 0.00590", "−4 %", "false minima"],
          ["B6 branch envelope (cold | lam | best) + SLSQP restarts", "lower envelope, clean BFGS", "0.005571", "−9.4 %", "first method that keeps its gains"],
          ["B7/B9 Tu continuation from below (0.1 → 0.5)", "laminar shapes first", "0.0059–0.0061", "−2…−5 %", "low-Tu shapes lose laminarity at 0.5 %"],
          ["B10 Tu continuation from above (1.0/0.75 → 0.5)", "plateau shapes first", "0.0062–0.0069", "worse than NACA", "BCM stalls at any Tu"],
          ["B11 SLSQP restarted from best ±0.01 probe point", "probe + descend", "0.005484", "−10.8 %", "adjoint result, V/V0 0.994"],
          ["Cold pattern search (derivative-free), 7 rounds", "compare cold values, no gradient", "0.005237 / 0.005052 / 0.005035", "−14.8 / −17.8 / −18.1 %", "V/V0 0.9997 / 0.988 / 0.973 — NOT adjoint"]],
      colw=[3.6, 2.2, 2.6, 1.8, 2.1], size=11)

# 7 results per Tu (drag decomposition)
s = prs.slides.add_slide(BL); title(s, "6. Article-style tables per free-stream Tu (cold, cl = 0.3)", "ΔCd vs the SA-BCM baseline at the same Tu; C1 (SA) = −3.0 % vs the SA baseline at any Tu")
table(s, [["Tu", "Baseline SA-BCM cd", "C2 cd (Δ)", "C3 cd (Δ)", "C4 cd (Δ)", "x_tr C2 up / lo"],
          ["0.5 % (article)", "6.15e-3", "5.48e-3 (−10.8 %) adjoint  ·  5.05e-3 (−17.8 %) derivative-free", "4.64e-3 (−24.5 %)", "5.45e-3 (−11.4 %)", "0.352 / 0.658  ·  0.420 / 0.713"],
          ["0.25 %", "5.12e-3", "4.04e-3 (−21.1 %) adjoint", "4.10e-3 (−19.9 %)", "4.59e-3 (−10.4 %)", "0.485 / 0.971"],
          ["0.1 %", "4.30e-3", "3.49e-3 (−18.9 %) adjoint", "3.84e-3 (−10.7 %)", "4.25e-3 (−1.3 %)", "0.634 / 0.982"]], colw=[1.4, 1.9, 4.2, 1.6, 1.6, 1.6], size=11)
bullets(s, ["Tu = 0.25 %: the laminar basin is wide enough for the gradient; C2 reaches the C3 level with a clear aft movement of transition (upper +0.08 c, lower +0.16 c; small laminar separation bubble near the TE on the lower surface at x/c ≈ 0.94).",
            "Tu = 0.1 %: the baseline is already very laminar (x_tr 0.52/0.92); C4 ≈ baseline, the 'SA misleads' message weakens.",
            "Tu = 0.5 %: with adjoint only, C2 (−10.8 %) is slightly below C4 (−11.4 %); ≥ −12 % needs the derivative-free stage. Full tables (Cd/Cdv/Cdp, %Cdv/%Cdp, x_tr): tables_by_tu.md."], y=3.4, size=13)

# 8 figure
s = prs.slides.add_slide(BL); title(s, "7. Skin friction, upper / lower surface, Tu = 0.1 / 0.25 / 0.5 %", "Baseline, C2 (adjoint result at 0.5 %), C3, C4 — cold, cl = 0.3")
s.shapes.add_picture(FIG, Inches(3.4), Inches(1.15), height=Inches(6.2))

# 9 what is/isn't adjoint
s = prs.slides.add_slide(BL); title(s, "8. What is reportable as C2 at Tu = 0.5 %")
table(s, [["Option", "C2", "Method", "How to describe it"],
          ["1", "−10.8 % (x_tr 0.35 / 0.66, V/V0 0.994)", "SLSQP + adjoint, branch-envelope evaluation + restarts", "gradient-based with a branch-consistent evaluation — clean"],
          ["2", "−4…−5 %", "SLSQP + adjoint exactly as in the article", "the original setup; shows the fragility"],
          ["3", "≈ −17 % (B12 running: SLSQP + adjoint restarted from the pattern-search points)", "hybrid", "must state an explicit derivative-free exploration stage"],
          ["4", "−21 % at Tu = 0.25 %", "SLSQP + adjoint (same evaluation)", "100 % adjoint; changes the article's Tu"]], colw=[0.9, 3.4, 3.8, 4.2], size=12)
bullets(s, ["Every optimum must be re-evaluated cold in a fresh process before it is quoted: values living on warm branches are not reproducible (this cost the HOP/B11 'gains').",
            "The pattern search (108 cold solves per round) is cheap only because the case is 2D (16 k cells) — in 3D the adjoint is the only option, and the adjoint itself is verified exact within each branch."], y=4.3, size=13)

# 10 proposals
s = prs.slides.add_slide(BL); title(s, "9. Proposals to remove the branch problem (ranked)")
bullets(s, ["1. Chordwise mesh refinement in the transition zone (0.25–0.65 c): the one-cell step ±1.8 % → ±0.5 % or less; SLSQP stops mistaking a step for a minimum. Most defensible in a paper ('grid-induced multiple equilibria'). Cost: new mesh + redo the 2D cases (~half a day). Untested.",
            "2. Branch-averaged objective: solve each design twice (cold and laminar-first), use the mean of cd and of the two adjoints → the staircase becomes a ramp consistent with the gradient. Still pure adjoint; ~1 h to implement in the envelope script; 2×2 solves per evaluation. Untested.",
            "3. Adjoint direction + sampled line search: keep the (good) gradient direction, replace SLSQP's line search by 6–8 cold evaluations along it. ~8 solves/iteration; custom loop, half a day.",
            "4. Report at Tu = 0.25 % (−21 %): zero work, changes the article's condition.",
            "5. Fixed cold protocol (done): reproducibility, not the optimizer.",
            "Outside the user's rules: soften the ν_t/ν feedback term of BCM (the source of the hysteresis) — a physics decision."], size=15)
note(s, "Recommendation: 1 for the article; 2 to test now with the existing scripts. Nothing launched without approval.")


# 11 dense mesh: C1 (SA) is mesh independent
s = prs.slides.add_slide(BL); title(s, "10. Dense mesh L0 (193k): SA optimisation is robust", "C1 = NACA0012 → SA, cl = 0.3, SLSQP + adjoint, article script; 16k article mesh vs L0 (one level below the NLF article mesh)")
s.shapes.add_picture(FIG_C1, Inches(0.4), Inches(1.3), height=Inches(5.9))
table(s, [["", "16k (article)", "L0 (193k)"],
          ["baseline SA cd", "9.77e-3", "9.42e-3"], ["C1 cd", "9.48e-3", "9.22e-3"], ["ΔCd", "−3.0 %", "−2.1 %"],
          ["α", "1.175°", "1.024°"], ["t_max (x/c)", "0.129 (0.38)", "0.129 (0.38)"], ["Cp,min (x/c)", "−0.60 (0.32)", "−0.61 (0.27)"],
          ["SLSQP exit", "0, 58 evals", "0, 119 evals"]], x=6.6, y=1.4, w=6.2, colw=[2.0, 2.1, 2.1], size=12)
bullets(s, ["Same optimum on both meshes: identical thickness distribution, camber and lower-surface flattening; only the Cp details differ (suction peak 5 % c further forward, softer TE recovery on L0).",
            "The coarse mesh over-predicts the baseline SA drag by 3.7 % and the gain by 1 point; the design itself is mesh independent → for fully turbulent flow the 16k mesh was enough."], x=6.6, y=4.6, w=6.3, h=2.2, size=12)

# 12 SA-BCM cases on the two meshes
s = prs.slides.add_slide(BL); title(s, "11. SA-BCM cases on the two meshes (Tu = 0.5 %, cl = 0.3)", "fs 0.08, ev 1e-7; drag counts; Δ vs the SA-BCM baseline on the same mesh; L0 C2/C3 still running (job 1926969, 32 ranks each)")
table(s, [["Case", "16k cd", "16k Δ", "16k x_tr up/lo", "L0 cd", "L0 Δ", "L0 x_tr up/lo", "L0 status"],
          ["Baseline SA", "97.7", "—", "—", "94.21", "—", "—", "done (NK recipe)"],
          ["Baseline SA-BCM", "61.49", "—", "0.317 / 0.624", "68.20", "—", "0.216 / 0.486", "done (cold trim)"],
          ["C1 (SA opt)", "94.8", "−3.0 % vs SA", "—", "92.19", "−2.1 % vs SA", "—", "done, SLSQP exit 0"],
          ["C2 (SA-BCM opt from NACA)", "54.84", "−10.8 %", "0.352 / 0.658", "58.80", "−13.8 %", "pending", "iter 18, plateau since iter 12; 13 failed trial solves"],
          ["C3 (C1 re-optimised with SA-BCM)", "46.43", "−24.5 %", "0.550 / 0.658", "51.92 (cl 0.313, not trimmed)", "≈ −24 %", "pending", "iter 3"],
          ["C4 (C1 geometry, SA-BCM)", "54.50", "−11.4 %", "0.469 / 0.590", "57.14", "−16.2 %", "0.419 / 0.523", "done (cold trim)"]],
      colw=[2.6, 1.0, 1.3, 1.5, 1.7, 1.2, 1.4, 1.6], size=11)
bullets(s, ["The SA baseline moves −3.7 % from 16k to L0; the SA-BCM baseline moves +11 % (transition on the NACA0012 moves forward: 0.32 → 0.22 c upper, 0.62 → 0.49 c lower). The transition location was not mesh-converged on the 16k mesh.",
            "On L0 the ordering C2 < C4 appears with the plain adjoint run (−13.8 % vs −16.2 % after 18 iterations, still going) — no derivative-free stage. C3 keeps ≈ −24 % on both meshes.",
            "L0 SA-BCM solves need L2Convergence 1e-8 (ANK stalls at ~1e-8); every second SLSQP trial on L0 lands on a failed/far solve (cl 0.1–0.5) and is rejected by the line search — same staircase mechanism as on 16k, smaller steps."], y=4.5, size=12)

# 13 mesh effect on Cp / Cf
s = prs.slides.add_slide(BL); title(s, "12. SA-BCM is mesh dependent: Cp, Cf 16k vs L0", "Tu = 0.5 %, cold trimmed solves at cl = 0.3; 16k dashed, L0 solid; C2/C3 on L0 added when the optimisations finish")
s.shapes.add_picture(FIG_MESH, Inches(0.3), Inches(1.25), height=Inches(4.9))
bullets(s, ["Baseline NACA0012: the L0 mesh trips the upper surface at 0.22 c instead of 0.32 c and the lower at 0.49 instead of 0.62; Cf on the laminar run is the same, the turbulent part is longer → +11 % cd.",
            "C4 (C1 geometry): transition on L0 at 0.42 / 0.52 c vs 0.47 / 0.59 c on 16k — the coarse mesh sustained a longer laminar run on both surfaces (x_tr = first cell with γ > 0.5 ten cells above the wall; on L0 that height is smaller).",
            "Per Tu (slide 7): the same figure on the article mesh at Tu = 0.1 / 0.25 / 0.5 %; the L0 baselines at 0.25 and 0.1 % are queued in the same verification job."], x=8.3, y=1.3, w=4.8, h=5.8, size=12)

# 13 GR cases on the two meshes -- same layout as slide 11 (SA-BCM), 2026-09-15
s = prs.slides.add_slide(BL); title(s, "13. SA-γ-R̅eθt cases on two meshes (Tu 0.5 %, cl 0.3)",
                                    "counts; Δ vs GR baseline on the same mesh; L1 = original GR setup (Re 6e6, ±0.03), L0 = article condition (Re 2.92e6, ±0.05); all finished 2026-09-16")
table(s, [["Case", "L1 cd", "L1 Δ", "L1 x_tr up/lo", "L0 cd", "L0 Δ", "L0 x_tr up/lo", "status"],
          ["Baseline SA", "—", "—", "—", "94.21", "—", "—", "L0: article campaign (NK recipe)"],
          ["Baseline SA-γ-R̅eθt", "73.41", "—", "0.098 / 0.363", "69.36", "—", "0.192 / 0.584", "done (cold trim, α 2.594° / 2.700°)"],
          ["C1 (SA opt)", "—", "—", "—", "92.19", "−2.1 % vs SA", "—", "L0: article campaign, SLSQP exit 0"],
          ["C2 (GR opt from NACA)", "42.91", "−41.5 %", "0.429 / 0.766", "34.94", "−49.6 %", "0.760 / 0.889", "L1: SLSQP exit 0, 82 majors, 14.8 h; L0: stopped flat at 126 solves, trimmed (job 1929206)"],
          ["C3 (C1 re-optimised with GR)", "—", "—", "—", "34.02", "−49.0 %", "0.774 / 0.900", "L0: job 1929321, flat after 122 majors, trimmed (1930514); article FFD"],
          ["C4 (C1 geometry, GR)", "—", "—", "—", "53.12", "−23.4 %", "—", "L0: cold trim, α 0.902° (job 1929321)"]],
      colw=[2.6, 0.8, 0.9, 1.4, 0.8, 1.1, 1.4, 3.3], size=11)
bullets(s, ["The GR baseline moves −5.5 % from L1 to L0 (73.4 → 69.4 counts) — but the two runs are NOT the same Re (6e6 vs 2.92e6): the lower Re on L0 delays transition (0.10 → 0.19 c upper, 0.36 → 0.58 c lower) and lowers cdv; a like-for-like mesh check needs the same condition on both (prepared: 2d_tu05_meshcmp, 16k / L1 / L0 at Re 2.92e6, not submitted).",
            "C2 on both: same mechanism — transition pushed to 0.43 / 0.77 c (L1, α 1.996°) and 0.76 / 0.89 c (L0, α 1.118°), gain mostly cdv (−42 % and −53 %); the L0 optimum is flatter and aft-loaded (cp peak at 0.70 c), the L1 one keeps a strong peak at 0.39 c and is thicker (14.5 % vs 12.2 %). L0 values are trimmed re-solves (cl = 0.3 ± 1e-5); C2 ≈ C3 on L0 (Δy ≤ 0.005 c): the GR optimum does not depend on the start point or the FFD.",
            "GR solves on both meshes need the CSANK-only recipe (nkswitchtol 1e-8, CFL 1e8, L2 1e-8): NK@1e-6 stalled at ~1e-6 abs and was flagged fail in every solve — the same stall the original Tu 0.0937 % runs had silently (96/99 fails in job 1803832). x_tr = steepest c_f rise per surface."],
        y=4.4, size=12)
note(s, "Data: $R/08_optimization/2d_tu05_L1 and 2d_tu05_L0 (trim_gammaretheta/trim.json, logs/opt_c2_*.log, surf/); x_tr from make_mesh_effect_gr.py.", y=7.05)

# 14 GR mesh effect figure -- same layout as slide 12 (SA-BCM)
FIG_GRMESH = "../figures/mesh_effect_gr_cp_cf_tu05.png"
s = prs.slides.add_slide(BL); title(s, "14. SA-γ-R̅eθt on the two meshes: Cp, Cf L1 vs L0",
                                    "Tu = 0.5 %, cl = 0.3; L1 dashed (Re 6e6), L0 solid (Re 2.92e6); baseline = cold trims, C2 = final best point; ▼ = x_tr")
s.shapes.add_picture(FIG_GRMESH, Inches(0.3), Inches(1.25), height=Inches(5.9))
bullets(s, ["Baseline NACA0012: identical cp on both meshes; the L0/Re 2.92e6 run keeps laminar flow to 0.19 c upper and 0.58 c lower vs 0.10 / 0.36 c on L1/Re 6e6 — a Re effect, not a mesh effect (transition Reθ,crit reached later at lower Re); the laminar cf level is the same.",
            "C2 optima: both optimizers push transition aft and flatten the lower surface; on L0 the upper-surface laminar run reaches 0.76 c with cf ≈ 0.001 up to there (lower 0.89 c), on L1 it ends at 0.43 c with a separation-bubble signature (cf dip to 0 then peak) — 42.9 vs 34.9 counts.",
            "Like-for-like mesh dependence (same Re on 16k / L1 / L0) is prepared in 2d_tu05_meshcmp and not yet run; GR C3/C4 on L0 prepared (article FFD + C1 shape)."],
        x=5.4, y=1.3, w=7.7, h=5.8, size=12)
note(s, "Data: same as slide 13; figure 08_optimization/make_mesh_effect_gr.py (surf CGNS of the trims and of the best C2 evaluations).", y=7.05)

# 15 GR C2 convergence histories (L1 and L0), 2026-09-15
FIG_HIST = "../figures/c2_history_L1_L0_tu05.png"
s = prs.slides.add_slide(BL); title(s, "15. SA-γ-R̅eθt C2 convergence: L1 vs L0",
                                    "Δcd vs the trimmed GR baseline per objective evaluation; red = best point so far; × = line-search trials off the cl constraint")
s.shapes.add_picture(FIG_HIST, Inches(0.3), Inches(1.3), width=Inches(12.7))
bullets(s, ["L1 (Re 6e6, ±0.03): −25.7 % at 20 evals, −39.5 % at 50, −41.5 % at 75, then flat for 140 more evaluations (SLSQP exit 0 at 82 majors; last 10 majors within ±2e-5 counts). L0 (Re 2.92e6, ±0.05): −23.5 % at 20, −29.9 % at 50, −36.9 % at 75, −46.6 % at 100, −49.6 % at 125 with a step at 90–100; stopped by the user at 126 solves (still creeping, ~0.1 count per 10 evals).",
            "Neither run has a plateau/branch problem: the objective is smooth in the DVs, failed solves are only the far line-search trials (cl 0.1–0.5), rejected by SLSQP."],
        y=5.35, h=1.8, size=12)
note(s, "Data: $R/08_optimization/2d_tu05_L1/logs/opt_c2_1927229.log, 2d_tu05_L0/logs/opt_c2_1927080.log; figure 08_optimization/plot_c2_history.py.", y=7.1)

# 16-20 cross-evaluation of the optima on L0 (2026-09-16) -- replaces the "current state" slide
import json
CE = "/home/mdo/Desktop/Run/MDO_PhD/Transition/gama_rethetha/08_optimization/2d_tu05_L0/cross_eval"
L0D = "/home/mdo/Desktop/Run/MDO_PhD/Transition/gama_rethetha/08_optimization/2d_tu05_L0"
def jr(path):
    d = json.load(open(path)); return d.get("trimmedAlpha", d.get("alpha")), d["cd"], d["cdp"], d["cdv"]
X = {"NACA0012 (baseline)":         (jr(f"{L0D}/trim_gammaretheta/trim.json"),        jr(f"{CE}/bcm_on_naca_base32_result.json")),
     "C1 = SA optimum (→ C4)":      (jr(f"{L0D}/trim_gr_on_c1/trim.json"),            jr(f"{L0D}/c1_artigo/bcm_c4_result.json")),
     "C2 GR (from NACA)":           (jr(f"{CE}/trim_gr_on_grc2/trim.json"),           jr(f"{CE}/bcm_on_grc2/out/result.json")),
     "C3 GR (from C1)":             (jr(f"{CE}/trim_gr_on_grc3/trim.json"),           jr(f"{CE}/bcm_on_grc3/out/result.json")),
     "C2 BCM (eval 57)":            (jr(f"{CE}/trim_gr_on_bcmc2/trim.json"),          jr(f"{CE}/bcm_on_bcmc2_e57_result.json")),
     "C2 BCM (restart, eval 50)":   (jr(f"{CE}/trim_gr_on_bcm_c2r1_e50/trim.json"),   jr(f"{CE}/bcm_on_bcm_c2r1_e50/out/result.json")),
     "C3 BCM (eval 40, cl 0.312 in opt)": (jr(f"{CE}/trim_gr_on_bcmc3/trim.json"),    jr(f"{CE}/bcm_on_bcmc3/out/result.json")),
     "C3 BCM (eval 118)":           (jr(f"{CE}/trim_gr_on_bcm_c3_e118/trim.json"),    jr(f"{CE}/bcm_on_bcm_c3_e118/out/result.json"))}
BG, BB = X["NACA0012 (baseline)"]
cnt = lambda v: f"{v*1e4:.2f}"; pct = lambda a, b: f"{100*(a/b-1):+.1f} %"
rows = [["Shape (row) ↓ · analysed with (column) →", "α", "cd  SA-GR", "cdp / cdv", "Δ vs 69.36", "α", "cd  SA-BCM", "cdp / cdv", "Δ vs 68.20"]]
for k, (g, b) in X.items():
    base = k.startswith("NACA")
    rows.append([k, f"{g[0]:.3f}", cnt(g[1]), f"{cnt(g[2])} / {cnt(g[3])}", "—" if base else pct(g[1], BG[1]),
                    f"{b[0]:.3f}", cnt(b[1]), f"{cnt(b[2])} / {cnt(b[3])}", "—" if base else pct(b[1], BB[1])])
s = prs.slides.add_slide(BL); title(s, "16. Cross-evaluation of the optima (L0, Tu 0.5 %)",
                                    "row = geometry, column = analysing model · cl 0.3 · counts · Δ vs NACA baseline of the SAME model · blue = own model")
t = table(s, rows, x=0.4, y=1.35, w=12.5, colw=[3.3, 0.7, 1.0, 1.4, 1.0, 0.7, 1.0, 1.4, 1.0], size=11)
for i, k in enumerate(X, 1):
    own = 2 if " GR" in k else (6 if "BCM" in k else None)
    if own: c = t.cell(i, own); c.fill.solid(); c.fill.fore_color.rgb = RGBColor(0xDD, 0xEB, 0xF7)
bullets(s, ["Method: optimum shape (FFD vector) → warped L0 mesh → trim to cl = 0.3 with SA-GR (CSANK-only, L2 1e-8) and with SA-BCM (cold, ev 1e-7, fs 0.08, L2 1e-8, 16 ranks = baseline branch). Same mesh, same Tu, same stopping criterion.",
            "Baselines agree (NACA −1.7 %, C1 +7 %): the primal is comparable. The optima do NOT: GR optima keep most of their gain in BCM (−25 % C2, −42 % C3 — the GR C3 is the best airfoil BCM has seen, 39.2 < 46.4); BCM optima are WORSE than the NACA in GR (+5.5 %, +14.6 %, +2.6 %) — only the intermediate C3 eval 40 survives (−28 %).",
            "The later the SLSQP-BCM point (e57 → restart e50, e40 → e118), the worse it is in GR: the last increments are gains of the transition criterion, not of the shape."],
        y=4.75, h=2.3, size=12)
note(s, "Data: jobs 1929206, 1929323, 1930476, 1930514 (trims), 1927080/1929321 (GR opts), 1926969/1928596 (BCM opts + cold verif); $R/08_optimization/2d_tu05_L0/cross_eval/{trim_gr_on_*/trim.json, bcm_on_*/out/result.json, PURPOSE.md}.", y=7.05)

for case, txt in (("C2", "GR shape (left): same cp in both models; SA-BCM trips the lower surface at 0.34 c (cp dip, cf → 0) where GR stays laminar to 0.89 c → 51 vs 35 counts. BCM shape (right): full nose with suction peak; GR transitions in the recovery at 0.16 c (cf 0.007), BCM only at 0.33 c → 80 vs 54 counts."),
                  ("C3", "GR shape (left): no lower-surface dip → SA-BCM stays laminar to 0.70 / 0.78 c → 39 counts (−42 %). BCM shape (right): suction peak at 0.02 c followed by an adverse gradient to 0.2 c — tolerated by BCM (laminar to 0.48 c), not by GR (0.17 c) → 71 vs 46 counts.")):
    s = prs.slides.add_slide(BL); title(s, f"{'17' if case == 'C2' else '18'}. {case} cross-evaluated: GR optimum | BCM optimum",
                                        "red = analysed with SA-GR, blue = with SA-BCM; dotted = baseline cf; ▼ = x_tr (steepest cf rise, same criterion for both models)")
    s.shapes.add_picture(f"../figures/cross_{case}.png", Inches(2.05), Inches(1.2), height=Inches(5.45))
    note(s, txt, y=6.7)

s = prs.slides.add_slide(BL); title(s, "19. Robustness: GR C2 ≈ C3 · BCM C3 e40 → e118",
                                    "left: GR C2 (from NACA) vs C3 (from C1) · right: BCM C3 eval 40 vs 118 · both models each · Δy = wall difference")
s.shapes.add_picture("../figures/gr_c2_vs_c3.png", Inches(0.5), Inches(1.25), height=Inches(5.45))
s.shapes.add_picture("../figures/c3_e40_vs_e118.png", Inches(6.9), Inches(1.25), height=Inches(5.45))
note(s, "GR: C2 and C3 differ by ≤ 0.005 c and give the same cp, cf, x_tr (0.76/0.89 vs 0.77/0.90) — a well-defined optimum, independent of start and FFD; the small difference (no lower-surface dip) is worth 12 counts in BCM (51 → 39). "
        "BCM: eval 40 → 118 gains 2 counts in its own model by adding a round nose with a suction peak; in GR that costs +22 counts (49.6 → 71.2).", y=6.75)

s = prs.slides.add_slide(BL); title(s, "20. Conclusion: the gain does not transfer in %",
                                    "L0, Tu 0.5 %, cl 0.3 — what the 2×2 matrix says")
bullets(s, ["Primal agreement on ordinary geometries (NACA −1.7 %, C1 +7 %) → comparing the models is legitimate; what diverges is what each optimizer exploits.",
            "Each optimizer pushes the laminar boundary layer to the separation threshold (cf → 0) that ITS transition criterion tolerates. SA-GR (Reθt transport, sensitive to the λθ history) rejects suction peaks followed by an early adverse gradient; SA-BCM (algebraic, local) rejects mid-chord cp dips on the lower surface.",
            "SA-GR has one well-defined optimum (C2 ≈ C3, two starts, two FFDs). SA-BCM has several (branches, rank dependence, restart dependence) and its late increments are criterion gains, not shape gains.",
            "GR optimum seen by BCM: −49 % → −42 % (keeps most). BCM optimum seen by GR: −32 % → +2.6 % (loses everything). The best airfoil in SA-BCM is the SA-GR optimum (39.2 vs BCM's own 46.4).",
            "For the article: report each model's gain only next to this cross matrix; the GR optimum is the robust candidate. Open: BCM C3 still running (eval 118 = last feasible); like-for-like mesh check (2d_tu05_meshcmp) not run."],
        y=1.4, size=15)
note(s, "Details, incidents and replication: $R/08_optimization/2d_tu05_L0/cross_eval/PURPOSE.md (2026-09-16).", y=7.05)

# 16 GR complete table at Tu 0.0937 % (the original ../2d campaign, 2026-08-04)
FIG_GR3 = "../figures/gr_2d_cases_tu00937.png"
s = prs.slides.add_slide(BL); title(s, "21. γ-R̅eθt full table, Tu = 0.0937 % (L1, Re 6e6)",
                                    "Original GR campaign (2026-08-04), cl = 0.3, trimmed re-solves (except C3*); counts; Δ vs same-model baseline, C3 vs its start C4")
table(s, [["Case", "Geometry", "Model", "α", "cd", "cdp", "cdv", "Δcd %", "Δcdp %", "Δcdv %"],
          ["—", "NACA 0012", "SA", "2.579", "88.75", "20.39", "68.36", "—", "—", "—"],
          ["C1", "NACA opt.", "SA", "0.710", "85.66", "16.93", "68.73", "−3.5", "−17.0", "+0.5"],
          ["—", "NACA 0012", "GR", "2.616", "68.75", "16.02", "52.73", "—", "—", "—"],
          ["C2", "NACA opt.", "GR", "2.170", "36.10", "8.88", "27.23", "−47.5", "−44.6", "−48.4"],
          ["C4", "C1 analysed", "GR", "0.611", "51.19", "8.85", "42.34", "−25.5", "−44.7", "−19.7"],
          ["C3", "C1 re-opt.*", "GR", "—", "36.53", "8.66", "27.88", "−28.6", "−2.2", "−34.2"],
          ["", "C3 vs base", "", "", "", "", "", "−46.9", "−46.0", "−47.1"]],
      x=0.4, y=1.35, w=7.6, colw=[0.55, 1.3, 0.65, 0.7, 0.75, 0.7, 0.7, 0.75, 0.75, 0.75], size=10)
s.shapes.add_picture(FIG_GR3, Inches(8.15), Inches(1.3), width=Inches(5.0))
bullets(s, ["Tu = 0.0937 % was calibrated at α = 0 on this mesh (set A of the α-sweep) — five times lower than the 0.5 % of slides 11–14; same FFD, bounds ±0.03, ladder with NK@1e-6 (which, as found on 2026-09-14, stalled at ~1e-6 abs in 96/99 solves of the C3 job — SLSQP ignored the flag).",
            "C4 shows the SA optimum is already a good transition shape for the wrong reason: −25.5 % almost entirely cdp; the γ-R̅eθt optimizer then only has cdv left (−34 %). Starting from NACA (C2) ends 1.2 % lower than from C1 (C3).",
            "* C3: best feasible point of job 1803832's log (cl 0.30004); its trim (job 1804386) was cancelled, α not recorded."],
        x=0.4, y=4.4, w=7.6, h=2.6, size=10)
note(s, "Data: $R/08_optimization/2d/{trim_*/trim.json, RESULTS.md, REPORT.md}; figure 2d/profiles/cp_cf_cases.png (make_cp_cf.py).", y=7.1)


# ---------------- 22-27: SA-BCM with a pressure-gradient sensor (2026-09-18) ----------------
import os, glob
PGD = "/home/mdo/Desktop/Run/MDO_PhD/Transition/gama_rethetha/08_optimization/2d_tu05_L0/cross_eval_pg"
s = prs.slides.add_slide(BL); title(s, "22. SA-BCM + pressure-gradient sensor (SABCM_PG)",
                                    "why: the BCM threshold Re_θc(Tu) ignores the gradient GR penalises (slides 16–20) · branch sa-bcm-pg, machV4 2026-09-18")
bullets(s, ["Mechanism found on slides 17–18: same cp in both models, only the cf jump moves. The BCM optimizer parks laminar flow in adverse gradients (nose suction peak + recovery) that GR trips through Re_θt(Tu)·F(λ_θ).",
            "Model (handoff option E, kept algebraic and local): λ_θL = −7.57e-3 (d²/ν) nᵀ∇u n + off (Menter 2015), λ = smooth clip(gain·λ_θL, ±0.1), F = Langtry–Menter F(λ) in the Piotrowski–Zingg smooth form (Eqs. 54–57) → Re_θc ← Re_θc^BCM · F. Everything else (KS max, tanh/fsmooth 0.08) untouched; PG off = bit-identical.",
            "n = exact unit wall normal (foot point → cell centre) from the approximate wall-distance update, stored per cell (new block array nWall) and differentiated (new independent in both Tapenade heads); helpers smoothMinMax/bcmFlambda ported from the GR branch.",
            "Options: SABCM_PG (False), _coef 7.57e-3, _off 0.0128 (Menter), _gain 1, _lamMax 0.1, _p 300. Runs use the calibrated off = 0.01623 and gain = 2 (next slide). Volume CGNS: BCM_lambda, BCM_Flambda.",
            "Not adopted from the handoff: the Ψ/tanh(k, δ) intermittency (the sa-bcm smooth form already removes those kinks; one effect at a time) and Menter's own F_PG (too weak at low Tu)."],
        y=1.45, size=14)
note(s, "Code: adflow_sabcm sa.F90 (saSource), blockette.F90 (mirror), wallDistance.F90 (nWall), turbUtils.F90 (helpers), Makefile_tapenade heads; docs/studies/22_bcm_pressure_gradient/PURPOSE.md.", y=7.05)

s = prs.slides.add_slide(BL); title(s, "23. Calibration on Falkner–Skan profiles (no CFD)",
                                    "for each β: exact λ_θ = θ²U'/ν vs the local sensor λ_θL evaluated at η* = argmax(η² f'') — where Term1 fires")
s.shapes.add_picture("../figures/pg_falkner_skan_gain.png", Inches(0.4), Inches(1.3), width=Inches(7.6))
bullets(s, ["Blasius at η*: λ_θL = −0.0162. Menter's offset 0.0128 leaves −0.0034 → F ≈ 0.96 on a flat plate (not neutral). Calibrated offset 0.01623 → neutral (the effective max over the layer changes it by < 1 %).",
            "Linear map λ_θ ≈ gain·(λ_θL + off): LS over |λ_θ| ≤ 0.1 gives gain 2.0 (adverse-only 1.4, favourable 3.0) — exact for |λ| ≲ 0.03 (suction-peak regime), over-predicts near separation (clipped at −0.1 anyway).",
            "Gains 1, 2, 3 swept on the cross set (slide 25); 2 kept."],
        x=8.2, y=1.3, w=4.8, h=5.5, size=12)
note(s, "docs/studies/22_bcm_pressure_gradient/falkner_skan_gain.py, falkner_skan_effective.py.", y=7.05)

def _cs_line(path):
    try: return [l for l in open(path).read().splitlines() if l.startswith("MAX")][0]
    except Exception: return "pending"
s = prs.slides.add_slide(BL); title(s, "24. Verification of the differentiated model",
                                    "reg suite (tutorial wing, 2 ranks) + NACA0012 L2 adjoint vs complex step · three sensors tried, one kept")
table(s, [["check", "PG off (smooth, hard)", "sensor 1 (Menter λθL)", "sensor 2 (wall dp/ds, kept)"],
          ["fwd/rev/fast AD vs references", "pass (AD unchanged)", "pass", "pass"],
          ["forward AD vs complex step", "pass", "pass", "pass"],
          ["fwd ↔ rev dot products, _b vs _fast_b", "pass", "pass (after in-place fix)", "pass"],
          ["adjoint vs CS, NACA L2, max |Δ| cl / cd", "0.008 % / 0.26 %", "187 % / 181 % — primal not smooth (FD: sign flips with step)", "0.000 % / 0.004 % (α), 0.001 % (DV12), 0.4 % (DV8, cd ≈ 1e-4) at the same state"],
          ["independent code audit (3 lenses + refutation)", "—", "—", "signs, stencil, Bernoulli, non-dim, mirror: correct; 3 fixes applied"]],
      x=0.4, y=1.4, w=12.5, colw=[3.4, 2.0, 3.4, 3.7], size=11)
bullets(s, ["Bugs found and fixed on the way: in-place Re_θc·F broke the fast-reverse state Jacobian (no push/pop stack); Tapenade head had the reference-state scalars as input-only (seeds zeroed at entry); F(0) = 0.9986 at p = 300 (smoothMax/Min bias) → exact tanh blend; ŝ = local velocity flips sign in bubbles → wall tangent oriented by the free stream.",
            "The 10–35 % adjoint-vs-CS 'gap' of sensor 2 was the two sides sitting on different steady states of the multi-state BCM primal (ANK/NK vs DADI path); formed at the same state they agree to 0.004 %.",
            "Sensor 3 (θ = Re_θc ν/U_e, the reference's own definition) is not smooth near stagnation (λ ∝ 1/U_e²): adjoint vs CS 46–105 %, trips both surfaces at 0.1 c — dropped."],
        y=4.15, h=2.8, size=12)
note(s, "Jobs: 1935762/1936305/1936446 (reg suite), 1935763/1935808 (sensor 1 CS + FD), 1936308/1936311/1936444 (loop it0–it2); adflow_sabcm docs/studies/22_bcm_pressure_gradient/PURPOSE.md §3e, §5.1.", y=7.05)

def _pg(g, shp):
    f = f"{PGD}/g{g}/bcmpg_on_{shp}/out/result.json"
    if not os.path.isfile(f): return None
    d = json.load(open(f)); return d["cd"], d.get("fail")
GRV = {"naca": 69.36, "gr_c2_best": 34.94, "gr_c3_best": 34.02, "bcm_c2r1_e50": 79.50, "bcm_c3_e118": 71.18}
BCV = {"naca": 68.20, "gr_c2_best": 51.38, "gr_c3_best": 39.22, "bcm_c2r1_e50": 53.55, "bcm_c3_e118": 46.39}
LBL = {"naca": "NACA 0012", "gr_c2_best": "C2 GR", "gr_c3_best": "C3 GR", "bcm_c2r1_e50": "C2 BCM (restart e50)", "bcm_c3_e118": "C3 BCM (e118)"}
GS = ["1.0", "2.0"]
LOOP = "/home/mdo/Desktop/Run/MDO_PhD/Transition/gama_rethetha/08_optimization/2d_tu05_L0/pg_loop/it2_s2_sdir2"
def _s2(shp):
    f = f"{LOOP}/trim_{shp}/out/result.json"
    return json.load(open(f))["cd"] * 1e4 if os.path.isfile(f) else None
rows = [["Shape", "SA-GR (target)", "SA-BCM", "Δ vs GR"] + sum([[f"sensor 1 g={g}", "Δ vs GR"] for g in GS], []) + ["sensor 2 (kept)", "Δ vs GR"]]
for shp in LBL:
    r = [LBL[shp], f"{GRV[shp]:.2f}", f"{BCV[shp]:.2f}", pct(BCV[shp], GRV[shp])]
    for g in GS:
        v = _pg(g, shp); r += ([f"{v[0]*1e4:.2f}" + ("*" if v[1] else ""), pct(v[0]*1e4, GRV[shp])] if v else ["—", "—"])
    v2 = _s2(shp); r += ([f"{v2:.2f}", pct(v2, GRV[shp])] if v2 else ["—", "—"])
    rows.append(r)
s = prs.slides.add_slide(BL); title(s, "25. Cross-evaluation with the sensor on (L0, Tu 0.5 %)",
                                    "BCM-PG vs GR, cl 0.3 · trims as slide 16 · sensor 1 (Menter, gain 1/2) and sensor 2 (wall dp/ds, K Blasius) · counts")
table(s, rows, x=0.4, y=1.35, w=12.5, colw=[2.4, 1.3, 1.1, 1.0, 1.3, 1.0, 1.3, 1.0, 1.4, 1.0], size=11)
bullets(s, ["Sensor 1 g=1 gets closest to GR's numbers (mean |Δcd| 26 % → 13 %) — by under-reading adverse gradients ~2× (Falkner–Skan map) — and is not differentiable; sensor 2 reads the LM edge λ_θ exactly (checked against the GR's own transported R̅eθt outside the layer: 950/560/966/667 vs 880·F = 970/537/970/560) and over-trips.",
            "What every variant achieves: the ranking becomes the GR one (GR C3 < GR C2 < NACA < BCM optima) — an optimiser driven by BCM-PG no longer exploits the nose suction-peak loophole; the BCM optimum that GR rejects (+14.6 %) is rejected by BCM-PG (+5…+7 %).",
            "What none can: GR C2 (+50 %) and the absolute levels — GR's onset fires at Re_S ≈ 2.46·R̅eθt (transported, lagged near the wall) vs the BCM's Re_v ≈ 2.19·Re_θc·F: with F = 0.61 in the short adverse zone at x = 0.06–0.10 the BCM-PG trips at 0.12 c, GR at 0.42 c. That gap is the base onset criterion, not the gradient sensitivity."],
        y=4.05, h=2.9, size=12)
note(s, "Data: jobs 1934907/08 (sensor 1), 1936444 (sensor 2, pg_loop/it2_s2_sdir2); $R/08_optimization/2d_tu05_L0/{cross_eval_pg,pg_loop}; study 22 PURPOSE.md §3d.", y=7.05)

s = prs.slides.add_slide(BL); title(s, "26. cf along the chord with the sensor on",
                                    "GR red · BCM blue dotted · BCM-PG green per gain · upper thick, lower thin · C3 GR and the BCM optimum move toward GR")
s.shapes.add_picture("../figures/pg_cross_cf.png", Inches(0.4), Inches(1.9), width=Inches(12.5))
note(s, "cross_eval_pg/make_cross_pg.py → cross_pg_cf.png (surfaces from the trims).", y=7.05)

s = prs.slides.add_slide(BL); title(s, "27. Direct predictions with the sensor, and the decision",
                                    "NLF0416 L1 (Tu 0.15 %) / S809 L1 (Tu 0.07 %), 3 α each, 11 ranks, rNK · counts · x_tr upper/lower")
table(s, [["case", "SA-BCM (no sensor)", "SA-BCM + sensor 2", "SA-γ-R̅eθt", "experiment x_tr,up"],
          ["NLF0416 α2", "58.2 · 0.38/0.60", "65.1 · 0.31/0.60", "57.0 · 0.38/0.60", "—"],
          ["NLF0416 α4", "65.5 · 0.33/0.61", "79.4 · 0.21/0.60", "64.4 · 0.34/0.61", "≈ 0.34"],
          ["NLF0416 α6", "85.6 · 0.21/0.61", "102.3 · 0.11/0.61 (not conv.)", "79.3 · 0.25/0.61", "—"],
          ["S809 α4", "47.6 (not conv.)", "48.7 (not conv.)", "67.5", "—"],
          ["S809 α6", "75.0 · 0.49/0.49", "74.6 · 0.50/0.49", "69.8 · 0.50/0.50", "—"],
          ["S809 α8", "149.1 · 0.96/0.50", "149.9 · 0.96/0.49", "131.8", "—"]],
      x=0.4, y=1.4, w=12.5, colw=[1.8, 2.6, 3.0, 2.6, 2.5], size=11)
bullets(s, ["S809 (transition by laminar separation): the sensor is neutral — as it should be. NLF0416 (attached transition in a mild adverse gradient): the local F(λ) trips the upper surface too early (α4: 0.21 vs 0.34 GR/experiment, cd +21 %); plain SA-BCM already matched experiment there.",
            "Same mechanism as slide 25: the local LM correlation applies F ≈ 0.8 where the reference's lagged near-wall R̅eθt and higher onset level do not. Sensor 1 (Menter, gain 2) behaved the same on this polar (0.19 at α4) and limit-cycled in ANK at NLF α6 / S809 α4–6 (no recipe helped: p, gain, CFL cap, early NK, ν̃∞ 1e-8/1e-10, restart from the converged PG-off state).",
            "Verdict on the direct-prediction criterion: the local sensor is NOT a net improvement — it fixes the ranking of optimised shapes at Tu 0.5 % but degrades the NLF polar. It stays in the code as an option (SABCM_PG, sensor 2 default) with a verified adjoint; closing the gap to GR for real would mean adopting GR's onset/transport, i.e. using GR.",
            "Decision (2026-09-19): line closed — \"worse, and no significant improvement\". The C2/C3 BCM-PG optimisation was cancelled at SLSQP major 1. Recommendation: optimise with SA-γ-R̅eθt (the robust model, slide 20), keep SA-BCM for cheap primal evaluation, leave SABCM_PG off. Kept from the work: exact wall-normal array + adjoint plumbing, two Tapenade traps documented, the one-state rule for adjoint-vs-CS checks on this multi-state primal."],
        y=4.25, h=2.75, size=11)
note(s, "Data: 15_sabcm_polars/{results/*_ref11,*_pg,*_pgs2, plot_polars_pg.py}, jobs 1934910/11, 1935125, 1935764, 1935776, 1936549; study 22 PURPOSE.md §5.3, §5.3b.", y=7.05)

# ---------------- 28-31 multipoint C2, GR vs BCM (added 2026-09-22; data 08_optimization/2d_tu05_L0_multipoint)
s = prs.slides.add_slide(BL); title(s, "28. Multipoint C2: does robustness rescue SA-BCM?",
                                    "L0, Re 2.9e6; one design, Σ cd/n over the points, one α per point; best feasible point of each history (36 h), counts")
table(s, [["case", "points (Tu %, cl)", "baseline cd", "SA-γ-R̅eθt optimum", "Δ GR", "SA-BCM optimum", "Δ BCM"],
          ["single-point C2", "(0.5, 0.3)", "69.4 / 68.2", "34.9", "−49.6 %", "48.7", "−28.6 %"],
          ["multi-Tu C2", "(0.5, 0.3) + (0.1, 0.3)", "69.4 · 62.7 / 68.2 · 51.4", "34.9 · 34.6", "−49.7 · −44.8 %", "60.2 · 50.9", "−11.7 · −1.0 %"],
          ["multi-cl C2", "(0.5, 0.2) + (0.5, 0.3) + (0.5, 0.4)", "65.4 · 69.4 · 73.9 / 66.2 · 68.2 · 71.0", "36.4 · 36.7 · 38.4", "−44 · −47 · −48 %", "65.6 · 49.0 · 49.8", "−1 · −28 · −30 %"]],
      x=0.4, y=1.4, w=12.5, colw=[1.6, 2.4, 2.9, 1.7, 1.5, 1.5, 0.9], size=10.5)
bullets(s, ["Story: the SA-BCM weakness is its response to pressure gradients (slides 4, 20, 25) — the optimiser can only exploit the nose suction peak. Putting several points in the objective does not cure that: multi-Tu BCM gains 12 % at Tu 0.5 and nothing at Tu 0.1; multi-cl BCM gives up the cl 0.2 point entirely (−1 %) to keep −28/−30 % at 0.3/0.4 — i.e. its single-point gain, no bucket.",
            "SA-γ-R̅eθt: the multi-Tu optimum IS the single-point optimum (34.9 at Tu 0.5, and 34.6 at Tu 0.1 for free); the multi-cl optimum forms a real drag bucket, 36.4–38.4 counts over cl 0.2–0.4, at a cost of 1.8 counts (5 %) at the design point.",
            "Why: GR moves transition with the whole pressure distribution (x_tr 0.55–0.63 / 0.9 on every point), so one aft-loaded shape works everywhere; BCM keeps x_tr ≈ 0.3 on the upper surface in every case and trades the points against each other. Same three-family picture as slide 20: GR optima are one family, BCM optima another.",
            "Baselines per point are cold trims of the NACA 0012 with each model (base 0.1 %: GR −10 %, BCM −25 % vs Tu 0.5). All 5 runs hit the 36 h walltime (SLSQP still running); multi-Tu GR was still creeping (−0.4 counts / 13 evals), the rest flat for > 20 evaluations. Second Tu = 0.1 %: lower end of both models' credible range (Mack correlation valid 0.1–2 %; LM Reθt calibration from ≈ 0.1 %), free-flight vs tunnel."],
        y=3.55, h=3.5, size=11)
note(s, "Data: $R/08_optimization/2d_tu05_L0_multipoint/{out_*/opt.hst, logs/opt_*_19308{55,56,57}.log, trim_base_*/trim.json, best_evals.txt}, jobs 1930855/56/57 (3 packed nodes, 32 ranks per point); PURPOSE.md there.", y=7.05)

s = prs.slides.add_slide(BL); title(s, "29. Multi-Tu C2: baselines vs GR and BCM optima",
                                    "Tu 0.5 % (left), 0.1 % (right), cl 0.3 · dotted = baselines · red GR (eval 111), blue BCM (eval 140) · ▼ x_tr")
s.shapes.add_picture("../figures/mp_overlay_models_tu.png", Inches(1.6), Inches(1.3), height=Inches(5.7))
note(s, "08_optimization/2d_tu05_L0_multipoint/plot_mp_overlay_models.py (surfaces of the best evaluations, surf/).", y=7.05)

s = prs.slides.add_slide(BL); title(s, "30. Multi-cl C2: baselines vs GR and BCM optima",
                                    "Tu 0.5 %, cl 0.2 / 0.3 / 0.4 · dotted = baselines · red GR (eval 124), blue BCM (eval 144) · ▼ x_tr")
s.shapes.add_picture("../figures/mp_overlay_models_cl.png", Inches(0.3), Inches(1.3), height=Inches(5.7))
note(s, "GR keeps x_tr ≈ 0.63 / 0.91 at all three cl; BCM ≈ 0.34 upper and, at cl 0.2, 0.29 on the lower surface → that point stays at its baseline. Same script as slide 29.", y=7.05)

s = prs.slides.add_slide(BL); title(s, "31. The optimum airfoils: GR family vs BCM family",
                                    "top row GR (single / multi-Tu / multi-cl), bottom row BCM · dotted NACA 0012 · thickness and camber under each")
s.shapes.add_picture("../figures/mp_foils_all.png", Inches(2.9), Inches(1.3), height=Inches(4.2))
bullets(s, ["GR: t_max 11.6–12.4 % at 0.54–0.58 c, camber 2.5–2.8 % at ≈ 0.5 c, thin nose, positive camber to 0.9 c — the same airfoil whatever the objective. BCM: t_max 13.4–14.3 % at 0.25–0.28 c, camber 0.7–1.7 % at 0.27 c, full nose, S-shaped camber (negative aft) — also the same whatever the objective.",
            "Conclusion for the article: SA-BCM is not a usable optimisation model; its pressure-gradient weakness is a property of the model, and multipoint formulations do not remove it. Optimise with SA-γ-R̅eθt; SA-BCM stays a cheap analysis model."],
        y=5.55, h=1.6, size=11)
note(s, "08_optimization/2d_tu05_L0_multipoint/plot_final_foils.py → foil_*.png; montage figures/mp_foils_all.png.", y=7.05)

# ---------------- 32-38 SA-G-s (one-equation smoothed SA-gamma) and the main conclusion (added 2026-09-22; data 22_sa_g_model)
def title2(s, t, sub=None):
    """wrapped title (the plain title() box does not wrap: long titles overflow in LibreOffice)"""
    tb = s.shapes.add_textbox(Inches(0.4), Inches(0.2), Inches(12.5), Inches(1.0)); tf = tb.text_frame; tf.word_wrap = True
    p = tf.paragraphs[0]; p.text = t; p.font.size = Pt(20); p.font.bold = True
    if sub:
        p2 = tf.add_paragraph(); p2.text = sub; p2.font.size = Pt(11); p2.font.color.rgb = RGBColor(90, 90, 90)
def pic_fit(s, path, x, y, w, h):
    from PIL import Image
    iw, ih = Image.open(path).size; r = min(w / iw, h / ih); pw, ph = iw * r, ih * r
    s.shapes.add_picture(path, Inches(x + (w - pw) / 2), Inches(y + (h - ph) / 2), width=Inches(pw))
KPI = RGBColor(0, 90, 160)
def kpi(s, x, y, w, big, small):
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(1.1)); tf = tb.text_frame; tf.word_wrap = True
    p = tf.paragraphs[0]; p.text = big; p.font.size = Pt(30); p.font.bold = True; p.font.color.rgb = KPI
    p2 = tf.add_paragraph(); p2.text = small; p2.font.size = Pt(11); p2.font.color.rgb = RGBColor(70, 70, 70)

# 32 the model
s = prs.slides.add_slide(BL); title2(s, "32. A third model: SA-G-s — smoothed γ equation, local onset, no R̅eθt",
    "transition-models 281a6d0a, option transitionLocalReTheta · machV3 · 22_sa_g_model/ · 2D opts at 32 ranks")
bullets(s, ["Why: no local sensor on the algebraic BCM closes the BCM-vs-GR gap (slides 22–27), and a transported threshold inside BCM fixes the BCM optima but breaks the GR optima and the NLF polar (transport in BCM vetoed). Question: does a one-equation SA-γ exist, and does it work?",
            "Literature (9-agent sweep, 22_sa_g_model/literature_memo.md): Menter-2015 one-equation γ + SA exists — Nichols 2019 (Kestrel), Lee & Baeder 2021, Jung 2022, Liu 2020, D'Alessandro 2025, STAR-CCM+ 'SA Gamma'. None smoothed / adjoint-ready; Jung: weaker than γ-R̅eθt at Re ≳ 3e6. P&Z's 'sLM2015' is the TWO-equation model.",
            "Implemented: the GR-s γ equation unchanged; onset (Re_θc, Flength) from the local LM Re_θt(Tu∞, λ_θL), Menter's wall-distance λ_θL = 7.57e-3 d²/ν dU/ds + 0.0128 (smooth clip ±0.1); R̅eθt kept as a diagnostic (option to make it inert). Names: G / G-s as GR / GR-s. Tapenade regenerated."],
        x=0.5, y=1.45, w=7.6, h=5.4, size=12)
kpi(s, 8.5, 1.5, 4.4, "0.008 %", "adjoint vs complex step (cd), GR reference 0.012 % — gate passed at the GR level")
kpi(s, 8.5, 2.8, 4.4, "251 its", "NLF0416 a4 with the opt ladder (GR 604); converges with no coupled ANK at all")
kpi(s, 8.5, 4.1, 4.4, "2.5 %", "NACA α 0–6 mean |Δcd| vs GR, x_tr ≈ GR; NLF polar 57.9 / 64.5 / 78.6 vs GR 57.0 / 64.4 / 79.3")
kpi(s, 8.5, 5.4, 4.4, "≈ GR", "cost per evaluation (same 3-eq block, adjoint 2× slower) — not cheaper as implemented")
note(s, "Data: 22_sa_g_model/{01_convergence_strategy (jobs 1937426/64/1939932), 02_derivative_check (1937469/70/78), 03_naca_sweep (1937402), 05_polars (1937404)}; docs CORE_BASE/CORE_01 option table.", y=7.0)

# 33 cross on given shapes
s = prs.slides.add_slide(BL); title2(s, "33. Cross-evaluation: on given shapes SA-G-s reproduces GR",
    "NACA0012 L0, Tu 0.5 %, cl 0.3, cold trims · red GR · blue BCM · green G-s · counts · job 1937403")
pic_fit(s, "../figures/gs_cross.png", 0.3, 1.3, 12.7, 5.55)
note(s, "GR optima: G-s 34.5 / 33.8 (GR 34.9 / 34.0; BCM 51.4 / 39.2). BCM optima: G-s 68.9 / 62.2 (GR 79.5 / 71.2; BCM 53.5 / 46.4) → the BCM loophole reads +3 % vs NACA in G-s (GR +15 %, BCM −21 %). Mean |Δcd| vs GR over the 5 shapes: BCM 26 %, G-s 6.3 %. 22_sa_g_model/04_cross_eval/plot_cross_gs.py.", y=6.9)

# 34 C3
s = prs.slides.add_slide(BL); title2(s, "34. G-s C3 vs GR C3: same airfoil, opposite verdict",
    "eval 122, 65 majors, 35 h · red GR opt in GR · green G-s opt in G-s · grey = G-s opt solved with GR (job 1939938)")
pic_fit(s, "../figures/gs_c3_vs_gr.png", 0.3, 1.35, 7.9, 5.5)
kpi(s, 8.6, 1.4, 4.4, "32.6 → 64.4", "counts: in G-s → in GR (GR C3 optimum: 34.0)")
kpi(s, 8.6, 2.7, 4.4, "−51 % → −7 %", "gain vs the NACA claimed by G-s → seen by GR")
kpi(s, 8.6, 4.0, 4.4, "0.78 → 0.43", "upper-surface x_tr: G-s keeps it laminar, GR trips right after the nose plateau")
bullets(s, ["The two optima differ by ≤ 1 % of chord; in their own models the cf curves coincide (x_tr 0.78/0.91 vs 0.77/0.90). The difference is entirely in what the other model says."], x=8.6, y=5.3, w=4.4, h=1.5, size=11)
note(s, "22_sa_g_model/06_optimization/{output_c3 (1937508), cross_gs_optima (1939938)}; surface files from the opt log (cruise_130) and the GR trim.", y=7.0)

# 35 C2
s = prs.slides.add_slide(BL); title2(s, "35. SA-G-s C2 vs GR C2 optimum: same verdict",
    "restart (no-CANK), best so far eval 39 (37.3) · red GR opt · green G-s in G-s · grey = G-s design solved with GR (job 1941002)")
pic_fit(s, "../figures/gs_c2_vs_gr.png", 0.3, 1.35, 7.9, 5.5)
kpi(s, 8.6, 1.4, 4.4, "37.3 → 66.0", "counts: in G-s → in GR (GR C2 optimum: 34.9)")
kpi(s, 8.6, 2.7, 4.4, "−44 % → −5 %", "gain vs the NACA claimed by G-s → seen by GR")
kpi(s, 8.6, 4.0, 4.4, "0.65 → 0.39", "upper-surface x_tr in G-s → in GR")
bullets(s, ["The first C2 run stalled at eval 50 (CANK plateau at α 3.9°, 12 adjoint failures; GR had 0); the restart with SANK → NK (no coupled ANK, R̅eθt inert) descends 50 → 37 counts in 40 evaluations. Whatever it reaches, GR will judge it like C3."], x=8.6, y=5.2, w=4.4, h=1.7, size=11)
note(s, "22_sa_g_model/06_optimization/{output_c2 (1937508), output_c2_r1 (1939986, 32 ranks), resolve_now/gr_on_gs_c2_r1_e39 (1941002)}; 01_convergence_strategy pack 3 (1939932): the stalled design converges with every ladder at fixed α — the stall was the optimiser's trials, not the mesh.", y=6.95)

# 36 why
s = prs.slides.add_slide(BL); title2(s, "36. Why: optimisers exploit any model without gradient history",
    "upper-surface cp / dcp/dx (> 0 adverse) / y(G-s) − y(GR) · left C3, right C2")
pic_fit(s, "../figures/gs_c3_why.png", 0.2, 1.35, 4.6, 5.5); pic_fit(s, "../figures/gs_c2_why.png", 4.7, 1.35, 4.6, 5.5)
bullets(s, ["GR optimum: dcp/dx < 0 continuously from 0.05 to 0.65 c — one monotone favourable ramp; transition 0.77.",
            "G-s optimum: 1 % more thickness at 0.15 c → suction peak (cp −0.55 / −0.72) followed by a MILD ADVERSE plateau to 0.35 c.",
            "A local onset (λ_θL) only sees the instantaneous gradient: a mild plateau never trips it. GR's transported R̅eθt is pushed down along those 20 % of chord (history) and trips at 0.43.",
            "BCM: same loophole with an instant switch and a stronger nose peak (Cp_min −1.45, slides 20/25). Both = 'peak, then mild adverse run'."],
        x=9.45, y=1.3, w=3.7, h=5.6, size=10.5)
note(s, "22_sa_g_model/06_optimization/{gs_c3_why.png, gs_c2_why.png}; cp from the opt surface files (cruise_130 / cruise_072) and the GR cross trims.", y=6.95)

# 37 main conclusion
s = prs.slides.add_slide(BL); title2(s, "37. Main conclusion: only GR optima can be trusted",
    "three transition models, one problem: NACA0012 → min cd at cl 0.3, Tu 0.5 %, Re 2.9e6, L0 mesh · what each is good for")
table(s, [["", "SA-BCM (algebraic)", "SA-G-s (1 eq., local onset)", "SA-γ-R̅eθt = GR (2 eq.)"],
          ["analysis of given shapes", "OK (NACA 2.4 %, polars = exp.)", "OK (2.5 %, cross 6 %)", "reference"],
          ["own optimum, in own model", "53.5 (−21 %)", "32.6 (−51 %)", "34.9 (−50 %)"],
          ["that optimum evaluated with GR", "79.5 → WORSE than NACA (+15 %)", "64.4 → −7 %", "34.9 → −50 % (survives cold re-trim, restarts, multipoint)"],
          ["cost per evaluation", "≈ 1/3", "≈ 1", "1"],
          ["usable for optimisation", "no", "no", "yes"]],
      x=0.4, y=1.35, w=12.5, colw=[2.6, 3.0, 3.0, 3.9], size=11)
bullets(s, ["The optimiser is a weakness detector. BCM (switch at Re_v ≥ 2.19 Re_θc(Tu)) and G-s (onset from the local λ_θL) both judge transition from the LOCAL state and carry no upstream pressure-gradient history — so the optimiser learns to build a suction peak followed by a mild adverse run: lower cdp, still 'laminar' for the model, tripped by a real boundary layer and by GR's transported R̅eθt.",
            "Every fix on the local side failed: local sensors on BCM (F=1 over-trips the NLF, F=2 does not close the peak), integral Thwaites criteria, a transported threshold inside BCM (breaks the GR optima), multipoint objectives (slides 28–31), Cp constraints (running). The history is the physics, and it is exactly the term the optimiser cannot game — it costs one more equation and the CANK phase.",
            "Recommendation: optimise with SA-γ-R̅eθt; use SA-BCM (and G-s) for cheap polars, screening and cross-checks; never trust an optimum from a local-onset model without re-evaluating it with GR."],
        x=0.5, y=3.85, w=12.3, h=3.1, size=11.5)
note(s, "Evidence: slides 20, 25–27 (BCM + sensors), 28–31 (multipoint), 33–36 (G-s); 22_sa_g_model/README.md; adflow_sabcm docs/studies/22 §6–7.", y=7.05)

prs.save("SA-BCM_2D_reruns_and_branches_2026-09-18.pptx"); print("ok")
