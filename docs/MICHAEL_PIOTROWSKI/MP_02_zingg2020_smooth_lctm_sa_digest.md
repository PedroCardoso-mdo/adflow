# MP_02 digest — Piotrowski & Zingg (2020), SA-sLM2015 (the canonical paper)

Companion to [`MP_02_zingg2020_smooth_lctm_sa_full.md`](MP_02_zingg2020_smooth_lctm_sa_full.md).

## Status in this repo — read this first

This is **the paper the whole branch implements**. It already exists in the KB
as the physics source of truth:

> [`../SA_GAMMA_RETHETHA_BASE/SAGR_01_paper_piotrowski_zingg_2020.md`](../SA_GAMMA_RETHETHA_BASE/SAGR_01_paper_piotrowski_zingg_2020.md)

That curated copy has **hand-transcribed LaTeX equations**; the `_full.md` in
this folder is a raw `pdftotext` extraction of the same PDF, where superscripts
and fractions are flattened. **For any equation, constant, or algorithm, use
`SAGR_01`, not the file next to this one.** Per CLAUDE.md rule 10, SAGR_01 wins
over code. The transcript here exists only so this folder is a complete record
of the PDF set.

## What this paper contributes that the others do not

- The full SA-sLM2015 model: smooth `F_onset`, smooth `F_length`, the vorticity
  limiting procedure, the SA coupling (Eq. 41, γ multiplies SA production).
- **Algorithm 2** (per-node exponential back-off damping on γ and Re̅θt) —
  ported to our NK path as `applyNKAlgorithm2Damping`.
- **Eq. 58** — the three-part scaling of the linear system (`S_a` autoscale,
  `S_r` row/equation, `S_c` column/variable).
- **Eq. 59** — the source-term Δt restriction.
- Subsonic validation set (NLF0416, S809, NLF2-0415, Sickle wing).

## Where the *other* files in this folder correct or extend it

This is the reason to keep this folder rather than relying on SAGR_01 alone —
four of the six siblings supersede parts of this paper:

| Topic | 2020 says | Superseded by |
|---|---|---|
| **`S_r` row scaling exponents** | given symbolically in terms of the SBP metric Jacobian `J` | **[`MP_06`](MP_06_piotrowski_phd_thesis_digest.md) Eq. 3.7 writes the matrix out explicitly** — `J^{2/3}` on the 5 mean-flow rows, `(scale·J^{1/3})^{-1}` on the 3 turb/transition rows. This resolves the `volRef`-vs-`J` ambiguity that made our `transitionRowVolScale` stall. |
| **Crossflow via `D_scf`** (sink in the Re̅θt equation) | the implemented form | **[`MP_05`](MP_05_zingg2023_compressibility_corrections_digest.md) replaces it** with `F_onset,scf` in the γ equation, plus a staged activation rule. |
| **Transonic applicability** | subsonic calibration only | [`MP_05`](MP_05_zingg2023_compressibility_corrections_digest.md) adds the `ψ` / `ψ_scf` compressibility corrections. |
| **Which linearization/coupling to use, and when** | fully coupled | **[`MP_03`](MP_03_zingg2022_numerical_behaviour_digest.md)** measures all four coupling strategies × source-Δt on/off, on three geometries — the definitive answer. |
| Solution update damping, physicality check, line search | described | [`MP_06`](MP_06_piotrowski_phd_thesis_digest.md) gives the **pseudocode** (Algorithms 2–4). |

## Practical note

If you are chasing a convergence problem, the 2020 paper is *not* the most
useful file in this folder — MP_03 (numerical behaviour) and MP_06 §3.1
(solver chapter) are. Come back here for the physics.
