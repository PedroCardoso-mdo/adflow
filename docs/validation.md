# Validation

What was computed, against which reference, and how far each result can be trusted.
Every number is quoted from `showcase/_manifest_2d.md` or `showcase/_manifest_3d.md`, which name the
file each one comes from. Solver settings: [recipes.md](recipes.md). Runnable versions of the cases:
`examples/transition/`.

**Provenance caveat that applies to every number on this page.** The runs were made over four months with
development builds of this code (seven different revisions, most with uncommitted changes at install
time). None of the converged solutions has been recomputed with the cleaned branch. What has been checked
on the cleaned branch is listed in [verification.md](verification.md).

Status labels: **validated** = compared with a reference and converged; **converged N orders** = residual
reduced N orders, no reference agreement claimed; **provisional** = will be redone; **attempt** = unfinished.

## 1. Overview

| Case                                            | Model(s)               | Reference                                   | Status                                                                                            |
| ----------------------------------------------- | ---------------------- | ------------------------------------------- | ------------------------------------------------------------------------------------------------- |
| NLF(1)-0416 and S809 polars, matrix dissipation | GR                     | Piotrowski and Zingg 2020 (P&Z), experiment | validated; 21 points, 8 orders                                                                    |
| Same polars, scalar dissipation, 193 224 cells  | GR, SA-BCM, SA-γ, SA   | P&Z, experiment                             | GR and SA-BCM validated (23/23 points, 8 orders); SA-γ provisional (16–17 of 23 converged)        |
| NACA 0012 mesh study, 6 levels                  | GR, SA-γ, SA-BCM       | Richardson extrapolation                    | converged 8 orders on all 18 runs; observed order irregular                                       |
| Free-stream ν̃ sweep                             | GR, SA-γ, SA-BCM       | —                                           | converged except where noted                                                                      |
| NLF(2)-0415 infinite swept wing, crossflow      | GR + crossflow, matrix | P&Z, experiment                             | validated; 7 Reynolds numbers, 8 orders                                                           |
| Same wing, scalar dissipation                   | GR + crossflow, scalar | P&Z                                         | negative result: converges, no crossflow front                                                    |
| Sickle wing, grid of the paper's authors        | GR + crossflow, matrix | P&Z                                         | converged 6.25 orders (target 8 not reached); forces near the paper; front comparison qualitative |
| AR5 wing with winglet                           | GR                     | none (mesh study)                           | converged 8 orders on 4 of 5 levels; not mesh-converged                                           |
| Flat plate, fully turbulent                     | SA                     | CFL3D / FUN3D (NASA TMR)                    | validated; +0.08 % in drag on the finest grid                                                     |
| Flat plate, free transition                     | GR                     | reference band from other codes             | converged 8 orders on the three grids; at the upper edge of the band                                            |
| 2D shape optimisation, single-point             | GR                     | —                                           | finished                                                                                          |
| 2D shape optimisation, single-point             | SA-γ, SA-BCM           | re-analysis with GR                         | settings not uniform (Section 8.2)                                                                |
| 2D shape optimisation, multipoint               | GR                     | —                                           | one finished; matrix series open                                                                  |
| 3D wing optimisation                            | GR (SA control)        | —                                           | attempt, unfinished; SA control converged                                                         |
| SA-BCM article cases (2D and 3D optimisation)   | SA-BCM                 | —                                           | numbers reproducible; central claims withdrawn (Section 9)                                        |

## 2. 2D airfoils against the paper (GR)

Conditions: NLF(1)-0416 M 0.10, Re 4.0e6, Tu 0.15 %; S809 M 0.10, Re 2.0e6, Tu 0.07 %.

### 2.1 Matrix dissipation (campaign `12_paper_conditions`, series 12-A)

`vis4` 0.04, `epsShear` 0, `epsAcoustic` 0.25, ν̃_∞ 5e-7; NLF(1)-0416 on 193 224 cells, S809 on 98 736 cells.
Jobs 1827893 and 1828015. All 21 points reduced the residual 8 orders (S809 α 7°: 1.03e-8).

| Airfoil     | α [deg] |     cl | cd [counts] |   cdp |   cdv | residual / initial |
| ----------- | ------: | -----: | ----------: | ----: | ----: | -----------------: |
| NLF(1)-0416 |      −4 | 0.0108 |       66.17 | 15.77 | 50.40 |            6.0e-09 |
| NLF(1)-0416 |       0 | 0.4957 |       52.59 | 11.16 | 41.43 |            7.9e-09 |
| NLF(1)-0416 |       2 | 0.7314 |       56.98 | 14.40 | 42.58 |            8.9e-09 |
| NLF(1)-0416 |       4 | 0.9610 |       65.32 | 20.45 | 44.86 |            6.9e-09 |
| NLF(1)-0416 |       6 | 1.1810 |       79.57 | 31.11 | 48.47 |            9.9e-09 |
| S809        |       0 | 0.1700 |       63.08 | 19.53 | 43.54 |            6.9e-09 |
| S809        |       2 | 0.4134 |       63.41 | 21.10 | 42.31 |            9.1e-09 |
| S809        |       4 | 0.6566 |       63.71 | 22.72 | 40.98 |            6.3e-09 |
| S809        |       5 | 0.7751 |       65.91 | 24.00 | 41.91 |            6.0e-09 |
| S809        |       6 | 0.8768 |       82.32 | 34.18 | 48.14 |            8.8e-09 |
| S809        |       8 | 1.0711 |      118.52 | 61.93 | 56.59 |            9.2e-09 |

Comparison with the paper:

| Quantity                                  | This code                       | P&Z                                           |
| ----------------------------------------- | ------------------------------- | --------------------------------------------- |
| NLF(1)-0416, α 4°, cd                     | 65.32 (cl 0.961)                | 65.40 (cl 0.964); experiment 64.36 (cl 0.893) |
| S809 drag bucket, α 0–5°                  | 63.1–65.9                       | 63–66                                         |
| S809 α 6° (past the corner of the bucket) | 82.3                            | 85                                            |
| Mean absolute difference from the paper   | NLF 1.1 counts (max 1.7), S809 1.2 counts (max 2.9) | digitised cd(α), 10 and 9 points |

The paper values were digitised from its figures; the table per α is
`showcase/validation_2d/01_gr_matrix_vs_paper/delta_vs_paper.md`. One digitisation of the S809 cd(α) curve had
its α axis one degree too high (the paper's points are at α 1…9°); it was corrected on 2026-11-13 and the
figure of the paper used as proof is stored beside the reference files.

### 2.2 Same polars on 193 224 cells, scalar and matrix (campaign `24_minmodel_2d`, 24-B and 24-C)

ν̃_∞ 2.79e-7. GR: 23 of 23 points reduced 8 orders with both schemes (S809 α 7° matrix from a separate
6-rank run).

| Case      | Scalar: cd | Matrix: cl | Matrix: cd |   cdp |   cdv | Matrix, series 12-A |
| --------- | ---------: | ---------: | ---------: | ----: | ----: | ------------------: |
| NLF α 0°  |      52.85 |     0.4958 |      52.52 | 11.17 | 41.35 |               52.59 |
| NLF α 4°  |      65.58 |     0.9611 |      65.22 | 20.43 | 44.79 |               65.32 |
| S809 α 0° |      63.43 |     0.1691 |      63.65 | 20.36 | 43.29 |               63.08 |
| S809 α 4° |      66.42 |     0.6557 |      63.86 | 22.96 | 40.90 |               63.71 |
| S809 α 5° |      68.13 |     0.7752 |      65.16 | 23.95 | 41.21 |               65.91 |
| S809 α 6° |      78.87 |     0.8771 |      81.34 | 33.58 | 47.76 |               82.32 |
| S809 α 7° |      97.10 |     0.9735 |     100.89 | 47.16 | 53.73 |               98.98 |
| S809 α 9° |     153.72 |     1.1629 |     138.23 | 78.45 | 59.78 |              137.72 |

Scalar dissipation against the paper on S809 at equal cl: +2 to +3 counts in the bucket, −8 counts at the
corner (cl 0.88–0.98), +12 to +15 counts above cl 1.06. The matrix results are the ones to compare with
the paper.

Transition location on NLF(1)-0416, distance to the measured band (fraction of chord, upper / lower
surface), scalar: GR 0.041 / 0.006, SA-γ 0.017 / 0.016, SA-BCM 0.024 / 0.013.

### 2.3 Superseded series (do not quote)

An earlier scalar series at ν̃_∞ 1.5e-3 and a series at the ADflow default 0.009 are not converged or are
8–28 counts high on the polar arms. The S809 scalar series at 2.79e-7 on 98 736 cells has inconsistent
convergence records.

## 3. Three models on one set-up (campaign `24_minmodel_2d`)

Scalar dissipation (`vis4` 0.01), ν̃_∞ per model (GR 2.79e-7, SA-γ 1e-7 or 2.79e-7, SA-BCM 1e-9).

### 3.1 NACA 0012, α 2.7°, Re 2.92e6, Tu 0.5 %, mesh study (24-A)

All 18 runs reduced 8 orders (jobs 1962375–79).

| Model           |   Cells |      cl | cd [counts] |   cdp |   cdv |
| --------------- | ------: | ------: | ----------: | ----: | ----: |
| GR              | 753 989 | 0.30169 |       63.45 | 12.93 | 50.52 |
| GR              | 381 240 | 0.30157 |       63.36 | 13.45 | 49.91 |
| GR              | 193 224 | 0.30113 |       62.19 | 14.24 | 47.95 |
| GR              |  98 736 | 0.29987 |       60.20 | 15.65 | 44.55 |
| SA-γ            | 753 989 | 0.30278 |       69.69 | 13.66 | 56.02 |
| SA-γ            | 193 224 | 0.30239 |       67.13 | 14.71 | 52.42 |
| SA-BCM (smooth) | 753 989 | 0.30445 |       67.69 | 13.24 | 54.45 |
| SA-BCM (smooth) | 193 224 | 0.30417 |       66.01 | 14.39 | 51.62 |

Richardson extrapolation (order 2 assumed from the two finest levels, safety factor 3): error of the
193 224-cell mesh in cd is −2.1 % (GR), −4.7 % (SA-γ), −3.0 % (SA-BCM). The observed order on cd is
irregular (2.6–7.4, because cdp and cdv compensate); a least-squares re-analysis gives an error of 3–6 %
and an uncertainty of 5–8 % on that mesh. Read the level as "a few per cent", not as a certified value.

### 3.2 Validation polars, all models (24-B), cd in counts

| Case      |     GR |                SA-γ | SA-BCM (smooth) | SA (fully turbulent) |
| --------- | -----: | ------------------: | --------------: | -------------------: |
| NLF α 0°  |  52.85 |               57.16 |           53.79 |                98.84 |
| NLF α 4°  |  65.58 |               69.28 |           64.87 |               112.76 |
| NLF α 6°  |  80.16 |  84.09 (5.6 orders) |           78.78 |               125.82 |
| S809 α 0° |  63.43 |               65.93 |           62.96 |               116.06 |
| S809 α 4° |  66.42 |               69.54 |           66.62 |               126.25 |
| S809 α 6° |  78.87 |               74.79 |           72.02 |               140.05 |
| S809 α 7° |  97.10 |               78.82 |           76.59 |               150.07 |
| S809 α 9° | 153.72 | 169.11 (3.2 orders) |          159.12 |               179.25 |

Convergence: GR 23/23, SA 23/23, SA-BCM 23/23 (one point at 7.5 orders), SA-γ 16–17 of 23 (the others at
3.2–6.0 orders). With matrix dissipation (24-C): GR 23/23, SA-BCM original 23/23, SA-γ 22/23 (two points
needed `sgammaArtVisc` 10; S809 α 8° stays at 1.3e-5).
At S809 α 6–7° only GR moves transition forward on the upper surface; the two local models stay 20 counts
lower.

Model forms of these polars. SA-γ: without artificial viscosity (`sgammaArtVisc` 0) except the two matrix
points noted above; this is the accepted form for the polars. SA-BCM: smooth form in scalar and original
form in matrix in the main set; the original form in scalar also exists on the same meshes (job 2168686,
ν̃_∞ ratio 1e-9: 22 of 23 points at 8 orders, NLF α 4° 64.73 counts, S809 α 6° 71.76; within 0.4 count of
the smooth form where both converged), in `showcase/validation_2d/03_four_model_polars/sabcm_original_scalar/`.
The later optimisations used the original form.

### 3.3 Free-stream ν̃ (24-D), NLF(1)-0416 α 4°, cd in counts

| `eddyVisInfRatio` |    GR | Note                                             |
| ----------------- | ----: | ------------------------------------------------ |
| 9e-3 (default)    | 75.51 | laminar region largely lost                      |
| 1e-5              | 66.49 |                                                  |
| 2.79e-7           | 65.58 | value used                                       |
| 1e-7              | 65.33 | S809 α 2° falls into a limit cycle at this value |

SA-BCM goes from 112.73 (9e-3) to 64.87 (1e-9); SA-γ from 72.08 (9e-3, not converged) to 68.16.
An optimised laminar shape is far more sensitive than the baseline: GR optimum 31.25 (2.79e-7) → 38.20
(1e-5) → 48.69 counts (9e-3), against 62.19 → 63.48 → 69.29 for the NACA 0012.

## 4. NLF(2)-0415 infinite swept wing: crossflow (campaign `13_nlf0415_swept`)

P&Z Section V.C. GR with `transitionCrossflow=True`, matrix dissipation (`vis4` 0.04, `epsShear` 0,
`epsAcoustic` 0.25), M 0.10, α −4°, sweep 45° (through the side-slip angle, periodic span), Tu 0.20 %,
roughness 3.3e-6, ν̃_∞ 2.79e-7, Reynolds number on the chord normal to the leading edge, 64 800 cells
(541 × 121 nodes; paper grid 81.6k nodes). Job 1826778, 9 ranks per case.

| Re / 1e6 | x_tr/c, this code | x_tr/c, paper | x_tr/c, experiment | vs paper |       cl | cd [counts] |   cdp |   cdv | Iterations | residual / initial |
| -------: | ----------------: | ------------: | -----------------: | -------: | -------: | ----------: | ----: | ----: | ---------: | -----------------: |
|    1.796 |             0.680 |         0.710 |              0.780 |   −4.2 % | −0.09793 |       71.08 | 17.63 | 53.45 |        591 |            8.3e-10 |
|    2.000 |             0.604 |         0.623 |              0.731 |   −3.1 % | −0.10154 |       72.85 | 17.76 | 55.09 |        615 |            6.7e-10 |
|    2.204 |             0.531 |         0.541 |              0.578 |   −1.9 % | −0.10399 |       73.98 | 17.74 | 56.24 |        634 |            2.1e-10 |
|    2.370 |             0.484 |         0.486 |              0.504 |   −0.4 % | −0.10552 |       74.69 | 17.66 | 57.03 |        627 |            3.7e-10 |
|    2.498 |             0.447 |         0.444 |              0.446 |   +0.7 % | −0.10658 |       75.24 | 17.57 | 57.67 |        633 |            1.6e-09 |
|    3.000 |             0.339 |         0.349 |              0.327 |   −2.9 % | −0.10888 |       76.06 | 16.96 | 59.11 |        495 |            2.1e-10 |
|    3.498 |             0.266 |         0.266 |              0.297 |   +0.1 % | −0.10984 |       75.76 | 15.91 | 59.84 |        582 |            2.1e-09 |

Mean absolute difference in x_tr: 1.9 % against the paper, 8.1 % against the experiment. x_tr is the
steepest rise of the smoothed upper-surface skin friction in 0.10 < x/c < 0.85. No force data is
available for comparison. The paper and experiment values were digitised.

| Check                            | Result                                                                                                                                                                                                                           |
| -------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Scalar dissipation (`vis4` 0.01) | Converges 8 orders in 7–14 min but x_tr ≈ 0.83–0.84 at every Reynolds number (bubble-triggered). Helicity parameter mid boundary layer 0.063–0.083 against the threshold 0.1066. Those runs used the wrong Reynolds length (√2). |
| Mesh refinement                  | 64 800 / 129 200 / 259 200 cells all converged; cd changes by 0.5–0.6 % (70.15 / 70.23 / 70.58 counts and 75.08 / 75.31 / 75.48). x_tr per level was not tabulated, and this check ran at the wrong Reynolds length.             |
| Reynolds length                  | Until the final sweep the Reynolds number was based on the streamwise chord (factor √2 too low); the resulting error grew with Re and looked like weak crossflow.                                                                |

## 5. Sickle wing: crossflow in 3D (campaign `07_sickle_wing`)

P&Z Section V.D. M 0.16, α −2.6°, Re 2.75e6, Tu 0.17 %, crossflow on, matrix dissipation, ν̃_∞ 5e-7,
medium grid of the paper's authors (4 404 960 cells). Run v14: 64 ranks, 6426 iterations, 17.7 h.

| Quantity    | This code                                                                                     | Paper (read from its figure) |
| ----------- | --------------------------------------------------------------------------------------------- | ---------------------------- |
| CL          | −0.0664                                                                                       | −0.0640 … −0.0660            |
| CD          | 64.72 counts (cdp 16.89, cdv 47.83)                                                           | 60 … 63 counts               |
| Convergence | total residual reduced 6.25 orders; density residual 3.4e-8 relative; target 1e-8 not reached | —                            |

Status: **converged 6.25 orders, comparison qualitative.** The transition front follows the experimental
line and moves forward at the kinks; with crossflow off on the same grid it misses those by about half a
chord, and CD is 12 counts lower. A station-by-station comparison of x_tr with the paper was not done.
The residual stayed flat for about 15 h while the crossflow front settled, then dropped 20× in 80 iterations.

Not validation: in-house meshes of the same wing give converged solutions (five levels, 1.3–22.5 M cells,
6.6–9.5 orders) but CL of the opposite sign (+0.054) and CD 47 % above the paper. The cause was not
closed; the geometry of the in-house meshes is the suspect. The authors' grid was received privately and
is not redistributed with this repository.

## 6. AR5 wing with winglet (campaign `11_ar5_corrected_foil`)

GR, crossflow off, M 0.0816, Re 1.127e6, α 1.8°, Tu 0.25 %, scalar dissipation at the ADflow defaults
(`vis4` 0.0156) and ν̃_∞ at the ADflow default 0.009. There is no reference solution; this is a solver and
mesh study.

|     Cells |       Orders reduced |      cl | cd [counts] |    cdp |   cdv | Onset x/c at 5 % span |
| --------: | -------------------: | ------: | ----------: | -----: | ----: | --------------------: |
|   459 452 |                 8.00 | 0.90955 |      459.28 | 423.43 | 35.84 |                  0.23 |
|   904 134 |                 8.07 | 0.93011 |      416.80 | 384.44 | 32.36 |                  0.28 |
| 1 783 442 |                 8.04 | 0.94849 |      388.56 | 356.16 | 32.40 |                  0.28 |
| 3 667 320 |                 8.05 | 0.95822 |      359.07 | 326.26 | 32.81 |                  0.33 |
| 7 415 408 | 7.61 (not converged) | 0.96054 |      339.05 | 304.29 | 34.76 |                  0.33 |

Status: **converged 8 orders on four levels; not mesh-converged.** cd still moves 5.6 % between the two
finest levels (all in cdp), so Richardson extrapolation is not usable. A separate chordwise ladder at
fixed spanwise resolution (254 to 752 points around the section, 8.6–27.2 M cells) shows that the
transition front converges only from 397 points up, at x/c = 0.62 ± 0.015; the family above (181 points)
gives 0.33. Over that ladder cd moves 4 % while the onset moves 63 %. Convergence depth of the chordwise
ladder runs was not tabulated (the 397-point mesh is quoted at 7.04 orders). The free-stream ν̃ of this
study predates the sensitivity study of Section 3.3.

## 7. Flat plate (campaign `10_tmr_flatplate`)

| Case                                                                                            | Result                                                                                                                                                                                                   | Status                                  |
| ----------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------- |
| SA, fully turbulent, NASA TMR conditions (M 0.2, Re 5e6 per unit length), total-pressure inflow | cd against CFL3D: +0.25 %, +0.16 %, +0.08 % on 137×97, 273×193, 545×385 points; upwind on the finest grid +0.02 %                                                                                        | validated                               |
| GR, free transition, Re 2.5e6, Tu 0.10 %, matrix dissipation, ν̃_∞ 5e-7, total-pressure inflow   | cd 20.05 / 19.67 / 19.51 counts on the three grids; reference band 18.4–19.5 counts with x_tr 1.0–1.2. Total residual reduced 8 orders on all three (ratios 1.8e-9 / 5.4e-9 / 9.7e-9 in the logs); job notes report x_tr 1.01–1.03. | converged 8 orders; upper edge of the reference band |
| GR, same conditions, scalar dissipation                                                         | coarse grid 18.09 counts, x_tr 1.055; medium grid fell to a laminar solution; fine grid stalled                                                                                                          | partial                                 |

An earlier +0.7 % error of the SA case was traced to a far-field inflow instead of total pressure.
A separate note reports x_tr ≈ 0.83 for a matrix setting on the coarse grid; it probably refers to a
far-field round and was not resolved.

## 8. Optimisation

All 2D optimisations: NACA 0012 start, free-form deformation box with 20 shape variables plus α, lift
constraint, volume ≥ initial, thickness bounds, SLSQP, 193 224-cell mesh unless stated.

### 8.1 Single-point, GR (finished)

| Design point                                                  | Scheme | Start cd [counts] |                   Optimum cd |  Change | Source                                                  |
| ------------------------------------------------------------- | ------ | ----------------: | ---------------------------: | ------: | ------------------------------------------------------- |
| Tu 0.3 %, cl 0.5                                              | scalar |              65.1 |                         31.6 | −51.5 % | `24_minmodel_2d`, 24-F                                  |
| Tu 0.3 %, cl 0.5                                              | matrix |              66.0 |                         30.7 | −53.6 % | 24-F                                                    |
| Tu 0.5 %, cl 0.3                                              | scalar |             62.19 |  31.24 (cdp 8.41, cdv 22.83) | −49.8 % | 24-E                                                    |
| Tu 1 %, cl 0.5                                                | scalar |              81.5 |                         44.8 | −45.0 % | 24-F                                                    |
| Tu 0.3 %, cl 0.6                                              | scalar |              70.4 |                         32.3 | −54.2 % | 24-F                                                    |
| Tu 0.3 %, cl 0.7                                              | scalar |              77.1 |                         33.3 | −56.8 % | 24-F                                                    |
| Tu 0.5 %, cl 0.3, earlier set-up (default ν̃_∞, `vis4` 0.0156) | scalar |             69.36 | 34.93 (cdp 10.34, cdv 24.59) | −49.6 % | `08_optimization`, job 1927080                          |
| Re 6e6, Tu 0.094 %, cl 0.3, 98 736 cells, earlier set-up      | scalar |             68.75 |                        36.10 | −47.5 % | `08_optimization/2d`; SA control 88.75 → 85.66 (−3.5 %) |

Caveats: all stopped on a stall rule (less than 0.1 count over 10 major iterations) or on the time limit,
not on the optimiser's criterion. The optimum from a second start (the SA optimum) reaches the same drag
(31.25 counts). The scalar optimum at Tu 0.3 % does not hold on the next finer mesh (31.6 → 40.7 counts on
381 240 cells) or in matrix dissipation (41.7); the matrix optimum does (30.7 on both meshes, the finer
one at five orders). The earlier-set-up optima were never re-trimmed.

### 8.2 Single-point, SA-γ and SA-BCM

Own model, Tu 0.3 %, cl 0.5; cd in counts and change against the NACA 0012 in the same model and scheme.

| Model              | Start (scalar / matrix) | Optimum, scalar | Optimum, matrix |
| ------------------ | ----------------------- | --------------- | --------------- |
| GR (for reference) | 65.1 / 66.0             | 31.6 (−51.5 %)  | 30.7 (−53.6 %)  |
| SA-γ               | 68.4 / 72.0             | 32.3 (−52.7 %)  | 46.5 (−35.4 %)  |
| SA-BCM (original)  | 65.5 / 68.3             | 31.3 (−52.2 %)  | 48.8 (−28.5 %)  |

**Settings are not uniform**: the SA-γ scalar optima ran with the tanh onset and no artificial viscosity and the matrix ones with coefficient 5; the SA-BCM matrix optimisation was still running when the numbers were taken; two sets
of scalar optima exist.

### 8.3 Do the optima of the simpler models hold in GR?

Optimum of each model re-analysed with GR at the same lift (old design point Tu 0.5 %, cl 0.3, scalar):

| Optimum from       | Gain in its own model |      Gain in GR |
| ------------------ | --------------------: | --------------: |
| GR                 |               −49.8 % |          (same) |
| SA-γ               |               −50.5 % | −37.5 % (38.85) |
| SA-γ, second start |               −46.0 % | −17.5 % (51.31) |
| SA-BCM (smooth)    |               −27.8 % | −26.0 % (46.05) |

At Tu 0.3 %, cl 0.5, matrix dissipation, judged by GR (NACA 66.0 counts): GR optimum 30.7 (−53.6 %);
SA-γ optimum 78.3 (+18.6 %); SA-BCM optimum 59.2 (−10.3 %, provisional). The comparison is circular in the
sense that GR is the judge; no independent method (for example e^N) has been applied.

GR optimum away from its design point (Tu 0.5 %, cl 0.3), gain in GR: −17.3 % at Tu 0.8 %, −6.9 % at Tu 1 %,
+1.4 % at Tu 1.5 %; −33.6 % at cl 0.4, −18.3 % at cl 0.5.

### 8.4 Multipoint

| Run                                                            | Model         | Result                                                                            | Status                |
| -------------------------------------------------------------- | ------------- | --------------------------------------------------------------------------------- | --------------------- |
| (Tu %, cl) = (0.5, 0.3), (1.0, 0.3), (0.5, 0.4), equal weights | GR, scalar    | −41.5 % / −40.7 % / −43.4 % per point; lift constraints met to 2e-5 (job 1980228) | finished (stall rule) |
| same                                                           | SA-γ, SA-BCM  | first runs lost in a node failure                                                 | provisional           |
| (0.3, 0.5), (1.0, 0.5), (0.3, 0.7)                             | GR            | −38.9 % / −39.9 % / −45.5 %                                                       | kept                  |
| same                                                           | SA-BCM / SA-γ | −28.5 / −15.5 / −37.8 %; SA-γ −40.5 / −27.3 / missing                             | SA-γ to redo          |
| (0.3, 0.5), (0.6, 0.5), (0.3, 0.7), matrix                     | GR            | no result written up                                                              | open                  |

A multipoint with a clean-air point (Tu 0.07 %) failed: the GR solve at that point did not converge in 40
of 69 evaluations.

### 8.5 3D wing optimisation with GR (attempt, unfinished)

Campaign `08_optimization/3d_c376`. AR5 wing, 13.68 M cells (397 points around the section), 182 shape
variables plus α, minimise cd at cl = 0.8932 with volume and thickness constraints, GR without crossflow,
scalar `vis4` 0.01, ν̃_∞ 5e-7, 4 nodes × 32 ranks, each solve to 7 orders, adjoint to 1e-3.

| Run                    | Tu     |             Start cd [counts] |                       Best cd |                         Change | cl at best (target 0.8932) | Major iterations | State                                  |
| ---------------------- | ------ | ----------------------------: | ----------------------------: | -----------------------------: | -------------------------: | ---------------: | -------------------------------------- |
| GR                     | 0.4 %  | 287.0 (cdp 255.28, cdv 31.73) | 269.0 (cdp 239.29, cdv 29.74) |                         −6.3 % |                    0.88728 |                6 | lift constraint not met                |
| GR                     | 0.5 %  |                         287.1 |                         269.5 |                         −6.1 % |                    0.88834 |               12 | still moving; last evaluation diverged |
| GR                     | 0.25 % |                         286.8 |                         275.3 |                         −4.0 % |                    0.88988 |                4 | incomplete                             |
| GR from the SA optimum | 0.5 %  |                         303.6 |                         276.3 | −3.8 % against the GR baseline |                    0.89865 |                5 | incomplete                             |
| SA control             | —      |                         362.8 |                         345.2 |                         −4.9 % |                    0.89320 |               22 | converged                              |

Status: **attempt.** The lift constraint is open by 0.5–0.7 %, so part of the drag reduction is less lift.
The work stopped when the compute allocation ran out. No optimum was re-analysed from a cold start and no
adjoint-against-complex-step check was made on this mesh. Cost: about 5 h per flow solve plus about 2 h
per adjoint. An earlier 3D run on a placeholder geometry (−3.9 %, stopped by hand) is superseded.

## 9. SA-BCM article cases and the retraction

A package for an article on a differentiable reformulation of SA-BCM exists (`showcase/bcm_article/`),
computed with an earlier code version and the old option names. The numbers reproduce:

| Case (cold, trimmed re-solves)                     |     cd [counts] |           cdv |             cdp | Change                                |
| -------------------------------------------------- | --------------: | ------------: | --------------: | ------------------------------------- |
| 2D NACA 0012, SA baseline                          |           94.21 |         75.50 |           18.71 | —                                     |
| 2D, SA-BCM baseline                                |           68.20 |         55.48 |           12.71 | —                                     |
| 2D, SA optimisation                                |           92.19 |         75.87 |           16.32 | −2.1 % (vs SA)                        |
| 2D, SA-BCM optimisation                            |           46.13 |         33.03 |           13.10 | −32.4 %                               |
| 2D, SA optimum re-optimised with SA-BCM            |           46.39 |         37.47 |            8.92 | −32.0 %                               |
| 3D AR5 wing (943 104 cells), SA baseline → optimum | 324.47 → 313.45 | 87.14 → 85.99 | 237.34 → 227.46 | −3.40 %                               |
| 3D, SA-BCM baseline (untrimmed) → optimum          | 244.57 → 240.53 | 32.34 → 31.74 | 212.23 → 208.79 | −1.65 % (about −1.4 % after trimming); smooth form |

The 3D SA-BCM optimum above is a smooth-form result. The same optimisation with the original form gave no
reduction (244.39 → 244.80 counts, job 1907421); files in `showcase/optimisation/3d_sabcm_ar5/`.

**Retraction.** Later tests (campaigns `27_bcm_conv_derivs`, `28_article1_bcm_revision`) withdrew the
central claims of that article. They must not be repeated on the strength of this repository.

| Article claim                                                           | Status after re-examination                                                                                                                        |
| ----------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| The original SA-BCM form does not converge                              | False. Different solver recipes had been compared; with the same recipe both forms converge to 1e-14 on the test wing.                             |
| Derivatives of the original form are wrong (median 12 %, maximum 320 %) | Wrong attribution. The primal was stalled (relative residual about 3e-7, a different solution) and the complex-step reference was under-converged. |
| Smoothing makes the derivatives more reliable                           | Not supported.                                                                                                                                     |
| The 2D SA-BCM optimum shows the failure of the original form            | Artefact of the dissipation coefficient: 46.13 counts with `vis4` 0.0156, 84.64 with 0.01, in SA-BCM itself.                                       |
| (new) The smooth form converges faster on the production mesh           | Supported at the optimisation tolerance 1e-6 only: 163 against 375 iterations on the baseline wing.                                                |
| (new) Optimisation with the smooth form advances further                | Trend only (−1.52 % against +0.17 % in 11 major iterations each).                                                                                  |

The latest record (2026-11-03) states that on the 943 104-cell wing neither form converges with a Newton
recipe and the original is equal or better in four jobs. That record was cited second-hand in the
manifest; the underlying files were not opened when the manifest was compiled.

## 10. Cost relative to SA

| Measurement                                                          |   SA | SA-BCM |   GR | Source campaign      |
| -------------------------------------------------------------------- | ---: | -----: | ---: | -------------------- |
| Time per iteration, NLF(1)-0416, 98 736 cells, 16 ranks, matrix [ms] | 24.3 |   25.0 | 62.9 | `17_cost_ab_clean`   |
| Time to 8 orders, same case [s]                                      |  529 |    554 |  699 | `17_cost_ab_clean`   |
| Memory per rank, same case [MB]                                      |  262 |    261 |  333 | `17_cost_ab_clean`   |
| Flow solve relative to SA (to 7 orders, scalar)                      | 1.00 |   1.00 | 1.84 | `18_paper_cond_cost` |
| Adjoint solve relative to SA                                         | 1.00 |   1.70 | 2.27 | `18_paper_cond_cost` |

GR scaling on one node: flow solve 1.0 / 1.9 / 3.2× faster on 32 / 64 / 128 ranks; adjoint 1.0 / 1.64 /
2.96× (campaign `23_scalability`). GR memory in 3D with NK: about 28e-6 GB per cell.

## 11. What is missing

| Gap                                                                                                        |
| ---------------------------------------------------------------------------------------------------------- |
| Converged solutions recomputed with the cleaned branch (none yet).                                         |
| A written per-α difference table and transition locations for the 2D matrix series.                        |
| Sickle wing to 8 orders; station-wise x_tr against the paper; explanation of the in-house-mesh sign of CL. |
| Mesh-converged AR5 forces; 3D optimisation finished and verified from a cold start.                        |
| Final SA-γ calibration and all SA-γ results; final SA-BCM matrix optima.                                   |
| An independent judge (not GR) for the cross-model comparison of optima.                                    |
| Any rotating case.                                                                                         |
