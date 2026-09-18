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
                                    "reg suite tests/reg_tests (tutorial wing, 2 ranks) + NACA0012 L2 adjoint vs complex step (20 shape DVs + α)")
table(s, [["check", "PG off (smooth, hard)", "PG on (new variant bcm_pg)"],
          ["fwd/rev/fast AD vs 2026-08-25 references", "pass (AD unchanged)", "n/a (refs trained now)"],
          ["forward AD vs complex step (5 tests)", "pass", "pass"],
          ["fwd ↔ rev dot products", "pass", "pass"],
          ["rev (_b) vs fast-rev (_fast_b)", "pass", "FAILED in build 1 → fixed (in-place Re_θc·F, no stack in _fast_b) → pass"],
          ["block == blockette residual", "pass", "(not parameterised)"],
          ["AD vs FD wDot (h 1e-8, rtol 8e-4)", "fails — also on the pre-change baseline", "fails (same test)"],
          ["NACA L2 adjoint vs CS, max |Δ| cl / cd", "0.008 % / 0.26 %", "build 1: 187 % / 181 % → after fix: " + _cs_line(f"{PGD}/cs_pg_fixed/adjoint_vs_cs_table_pg.txt").replace("MAX PG: ", "")]],
      x=0.4, y=1.4, w=12.5, colw=[4.0, 3.5, 5.0], size=11)
bullets(s, ["The fast-reverse bug is the classic one the GR code documents for its clamps: Tapenade guards `x = x·f` with pushreal8/popreal8 in _b, the fast-reverse autoEdit strips the real stack. Forward AD and CS were right all along — only the adjoint's state Jacobian was wrong.",
            "Two more Tapenade traps hit and recorded in the source: no local may end in _c (autoEdit strips _cd/_cb), module variables must be copied to locals before differentiated calls."],
        y=4.6, h=2.2, size=12)
note(s, "Jobs: 1934943/1935117 (reg suite, before fix), 1935762 (after), 1934905 → 1935763 (NACA CS); cross_eval_pg/cs_pg*/adjoint_vs_cs_table_pg.txt.", y=7.05)

def _pg(g, shp):
    f = f"{PGD}/g{g}/bcmpg_on_{shp}/out/result.json"
    if not os.path.isfile(f): return None
    d = json.load(open(f)); return d["cd"], d.get("fail")
GRV = {"naca": 69.36, "gr_c2_best": 34.94, "gr_c3_best": 34.02, "bcm_c2r1_e50": 79.50, "bcm_c3_e118": 71.18}
BCV = {"naca": 68.20, "gr_c2_best": 51.38, "gr_c3_best": 39.22, "bcm_c2r1_e50": 53.55, "bcm_c3_e118": 46.39}
LBL = {"naca": "NACA 0012", "gr_c2_best": "C2 GR", "gr_c3_best": "C3 GR", "bcm_c2r1_e50": "C2 BCM (restart e50)", "bcm_c3_e118": "C3 BCM (e118)"}
GS = ["1.0", "2.0", "3.0"]
rows = [["Shape", "SA-GR (target)", "SA-BCM", "Δ vs GR"] + sum([[f"BCM-PG g={g}", "Δ vs GR"] for g in GS], [])]
for shp in LBL:
    r = [LBL[shp], f"{GRV[shp]:.2f}", f"{BCV[shp]:.2f}", pct(BCV[shp], GRV[shp])]
    for g in GS:
        v = _pg(g, shp); r += ([f"{v[0]*1e4:.2f}" + ("*" if v[1] else ""), pct(v[0]*1e4, GRV[shp])] if v else ["—", "—"])
    rows.append(r)
s = prs.slides.add_slide(BL); title(s, "25. Cross-evaluation with the sensor on (L0, Tu 0.5 %)",
                                    "BCM-PG vs GR, cl 0.3 · trims as slide 16 · gain 1 / 2 / 3, off 0.01623 · counts · * = solve stopped at ~1e-6")
table(s, rows, x=0.4, y=1.35, w=12.5, colw=[2.6, 1.3, 1.1, 1.0, 1.3, 1.0, 1.3, 1.0, 1.3, 1.0], size=11)
bullets(s, ["Distance to GR over the four shapes: mean |Δcd/cd_GR| 24 % (BCM) → 16 % (BCM-PG, g 1–2).",
            "The case that motivated the work is fixed: the BCM optimum that GR rejects (+14.6 % vs NACA) is now rejected by BCM-PG too (g=2: 83.6 vs its NACA 73.3 = +14 %); the GR C3 ends within 4 % of GR (32.5 vs 34.0). Gain 2 reproduces the GR ranking best.",
            "Not fixed: GR C2 (+50 %): every BCM trips the LOWER surface at 0.3 c where GR stays laminar to 0.9 c (next slide). F(λ) can only lower the threshold (adverse) or raise it by ≤ 10 % at Tu 0.5 % — this gap is the base threshold itself (BCM Re_θc(0.5 %) = 725 vs Langtry Re_θt = 880), not the gradient.",
            "NACA moves +6 % with g=2 (upper x_tr 0.16 vs GR 0.19 — closer than BCM's 0.22; lower 0.45 vs 0.58)."],
        y=4.05, h=2.9, size=12)
note(s, "Data: jobs 1934907 (g2), 1934908 (g1), 1934909 (g3), 1935767 (e118 sweep); $R/08_optimization/2d_tu05_L0/cross_eval_pg/g*/bcmpg_on_*/out/result.json; PURPOSE.md there.", y=7.05)

s = prs.slides.add_slide(BL); title(s, "26. cf along the chord with the sensor on",
                                    "GR red · BCM blue dotted · BCM-PG green per gain · upper thick, lower thin · C3 GR and the BCM optimum move toward GR")
s.shapes.add_picture("../figures/pg_cross_cf.png", Inches(0.4), Inches(1.9), width=Inches(12.5))
note(s, "cross_eval_pg/make_cross_pg.py → cross_pg_cf.png (surfaces from the trims).", y=7.05)

s = prs.slides.add_slide(BL); title(s, "27. Polar at low Tu and the flat plate: what is still open",
                                    "NLF0416 L1 (Tu 0.15 %) / S809 L1 (Tu 0.07 %), 3 α each, 11 ranks, rNK · TMR plate L1/L2, Tu 0.10 %")
s.shapes.add_picture("../figures/pg_profile_nlf_a4.png", Inches(0.3), Inches(1.6), width=Inches(5.7))
bullets(s, ["Polar: with PG on, NLF α6 and S809 α4/α6 limit-cycle in ANK (resrho 0.2–0.7 for 5000+ its, NK never engages) while the PG-off twins converge; none of {p 100, gain 1, CFL 1e4, NK@1e-3} helps. Continuation from the converged PG-off state: job 1935764.",
            "Where it converges, PG trips NLF α4 earlier than GR/experiment (x_tr,up 0.19 vs 0.33/0.34; 81 vs 65 counts): at the Re_v peak (d/c ≈ 5e-4) the local sensor reads λ ≈ −0.005 (figure: profiles at x/c 0.05–0.20, λ = offset at the wall, negative through the layer, positive outside), a mild adverse gradient the GR's transported edge-λ_θ does not see. S809 α8: unchanged.",
            "Flat plate: the plain SA-BCM is bistable there at Tu 0.1 % (rNK/CANK → fully turbulent, cd 0.00317 = SA; ANK-only → laminar, no convergence), so it cannot test neutrality; the Blasius neutrality rests on the Falkner–Skan analysis and on λ = offset at the wall (Snn → 0).",
            "Next: C2/C3 with BCM-PG (gain 2) on L0 (dense_L0v2_pg) once the HPC adjoint check is in; low-Tu convergence needs its own study (continuation in gain, or the R̅eθt-equation plan B)."],
        x=6.2, y=1.35, w=6.8, h=5.5, size=11)
note(s, "Data: 15_sabcm_polars/{results/*_pg,*_ref11, plot_polars_pg.py, analyze_pg_nearwall.py}, jobs 1934910/11, 1935125, 1935764; 10_tmr_flatplate/10_bcm_pg jobs 1934951, 1935139.", y=7.05)

prs.save("SA-BCM_2D_reruns_and_branches_2026-09-18.pptx"); print("ok")
