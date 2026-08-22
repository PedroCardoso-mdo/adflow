# A Revised One-Equation Transitional Model for External Aerodynamics

> **Reference transcription for code implementation.** Equations, constants and boundary
> conditions have been transcribed verbatim from the source PDF. Figures are validation
> results (not part of the model definition) and are given as captions plus a short
> description of what each plot contains. Editorial remarks added during transcription are
> clearly marked **[Transcriber note]** and are *not* part of the original paper.
> See the companion file for AIAA 2020-2714 and the shared **Implementer's reconciliation notes**
> at the bottom of that file — several symbols in this paper are defined differently there.

| | |
|---|---|
| **Title** | A Revised One-Equation Transitional Model for External Aerodynamics |
| **Authors** | Samet Çaka Çakmakçıoğlu (TAI – Turkish Aerospace Industries, Ankara, Turkey); Onur Baş (TED University, Ankara, Turkey); Riccardo Mura (Università degli Studi di Cagliari, Cagliari, Italy); Ünver Kaynak (Yıldırım Beyazıt University, Ankara, Turkey) |
| **Venue** | AIAA Aviation Forum, June 15–19, 2020, Virtual Event — AIAA Aviation 2020 Forum |
| **DOI** | 10.2514/6.2020-2706 |
| **Model** | BC → **BCM** (Bas–Cakmakcioglu with Modifications), coupled to Spalart–Allmaras (SA) |
| **Solver** | SU2 (Stanford University Unstructured), open source |

---

## Abstract

A modified version of the BC zero-equation (algebraic) transition model is proposed in order
to eliminate some deficiencies in the original model reported by the experts. The newly
formulated BCM algebraic transition model still relies on local flow information and coupled
with the SA turbulence model in a way that the intermittency factor is multiplied with the
production term of the turbulence model. In the BCM transition model, first, the lack of
Galilean invariance in the original model due to inclusion of the local velocity magnitude in
the equations is removed. Second, the appearance of the Re number in formulation of the
original model is eliminated, which before caused some problems in a way that it brings a
possibly non-unique reference length into play. After fixing these two deficiencies, the new
BCM model is calibrated against the well-known natural transition flat plate test data of
Schubauer and Klebanoff. After the calibration, the BCM model is tested against the
zero-pressure gradient flat plate test cases of Savill and the results are compared with the
original model. It is observed from the results that the original model and the BCM model give
nearly the same results. Finally, the new formulation will be tested against an airfoil, the
Eppler E387, a turbine cascade case, the T106, and a 3-D test case of flow over a 6:1 prolate
spheroid.

## Nomenclature

| Symbol | Definition |
|---|---|
| $C_f$ | skin friction coefficient |
| $Tu_\infty$ | freestream turbulence intensity (%) |
| $U_\infty$ | inlet / free stream velocity |
| $Re$ | Reynolds number |
| $Re_\theta$ | momentum thickness Reynolds number |
| $Re_{\theta c}$ | transition onset momentum thickness Reynolds number ($C_f$ starts to grow) |
| $Re_v$ | vorticity Reynolds number |
| $\gamma_{BC}$ | BC model intermittency function |
| $\nu_T$ | eddy viscosity |
| $\nu_{L}$ | laminar viscosity |
| $\Omega$ | absolute value of the vorticity |

> **[Transcriber note on notation]** The nomenclature lists $\nu_T$ as "eddy viscosity", but in
> the transport equation (Eq. 1/2) $\nu_T$ is used as the SA *transport (working) variable*
> — i.e. what the SA literature and the companion paper 2020-2714 call $\tilde{\nu}$. In the
> Term2 definition (Eqs. 7–8) the symbol $\nu_t$ (lower-case *t*) denotes the *actual eddy
> viscosity* $\mu_t/\rho$. Watch this $\nu_T$ vs $\nu_t$ distinction when coding.

---

## I. Introduction

In most industrial CFD applications, due to usage of fully turbulent models such as the
one-equation Spalart-Allmaras (SA) and the two-equation k-ω SST, laminar to turbulent
transition characteristics of the flow are neglected. Among several reasons for ignoring
transitional flow characteristics, the most significant one is the difficulties arising in
incorporating different transition mechanisms (natural, bypass, separation induced, cross-flow,
etc.) into the unified models. Yet, in the last couple of decades, many successful transition
models coupled with the Reynolds-Averaged Navier-Stokes (RANS) equations have been proposed.
Among these models, the most famed one is the $\gamma - Re_\theta$ model based on the Local
Correlation-Based Transition Modeling (LCTM) concept proposed by Menter et al., which instead of
attempting to actually model the different transition processes, formulate a set of
CFD-compatible transport equations which allow combining experimental correlations in a local
fashion with the underlying turbulence model, the k-ω SST. The model is built on the
experimental observation that the maximum of the vorticity Reynolds number profile is
proportional to the momentum thickness Reynolds number in a Blasius boundary layer. This
observation provided an opportunity to use the vorticity Reynolds number $Re_v$ (a locally
calculated parameter) instead of the momentum thickness Reynolds number $Re_\theta$ (an integral
parameter).

The underlying experimental correlations of the original $\gamma - Re_\theta$ model were not
published until 2009, and because of that there were many researchers who attempted to create
alternative formulations for those correlations. After the full publication of the
$\gamma - Re_\theta$ model, many researchers worked on extending its capabilities into the realm
of more physics such as the cross-flow instability effects by Grabe and Krumbein, the so-called
secondary effects such as roughness by Dassler et al., and compressibility by Kaynak. Besides
the correlation-based modeling approach, a number of successful physics-based models also
appeared in the literature such as the laminar kinetic energy model ($k$-$k_L$-$\omega$) of
Walters and Leylek and Walters and Cokljat. Later, Bas and Cakmakcioglu introduced an algebraic,
or a zero-equation model (namely the BC model) that successfully reproduces the results of the
two- and three-equation models. Meanwhile, Menter et al. proposed a one-equation $\gamma$-model
that is a simplified version of the earlier $\gamma - Re_\theta$ model which possesses similar
quality of results as in the original model with one less equation. Finally, Kubacki et al.
proposed a new algebraic transition model, which once more demonstrated that an equivalent level
of predictive capability compared with the higher order models can be obtained as long as
correctly modeled physics is maintained.

In the BC transition model developed by Bas and Cakmakcioglu, similar to Menter et al. model, the
use of non-local information was avoided by using the vorticity Reynolds number to trigger the
transition through the use of a novel intermittency function formulation rather than an
intermittency transport equation. In this way, the intermittency equation is directly fed into
the production term of the one-equation SA turbulence model. Therefore, the model achieves a
similar effect by using two-less equations (including the additional $Re_\theta$ equation beside
the $\gamma$-equation) compared to Menter et al. model. However, there were some deficiencies of
the model reported by the experts, mentioned in the NASA Langley Research Center Turbulence
Modeling Resource website of the model. One of the deficiencies, although does not pose a severe
limitation, is the lack of Galilean invariance due to the local velocity magnitude in the
equations; and the second one is the appearance of the Re number in the formulation, which may be
problematic due to inclusion of non-unique reference length.

In this paper, the aforementioned deficiencies are removed by re-formulating the Term2 of the
original BC model. With the re-formulation, the modified BC transition model (BCM model) has
become Galilean invariant which is desirable for general-purpose CFD codes, and avoided the Re
number appearing in the model equations while very similar accuracy is maintained compared to the
original formulation. The BCM model is first calibrated against the natural transition flat plate
test case of Schubauer and Klebanoff. Then, the zero pressure gradient flat plate test data of
Savill is validated using the new formulation. In the full paper, the model will be tested against
an airfoil case, which is the Eppler E387, and a turbine cascade case, the T106. Finally, the
model will be employed to solve the flow over a 6:1 prolate spheroid geometry.

---

## II. Computational Method

### A. Flow Solver

Stanford University Unstructured (SU2), an open-source multi-physics solver, is used as the flow
solver. The SU2 is capable of solving 2-D and 3-D incompressible/compressible RANS equations using
linear system solver methods such as LU-SGS, BiCGSTAB and GMRES. In the SU2, central or upwind
methods are used for the discretization of the convective terms. For spatial gradients,
Green-Gauss and weighted-least-squares methods; for temporal discretization Runge-Kutta explicit,
Euler explicit and Euler implicit schemes are available. Also, several state-of-the-art numerical
schemes such as Roe and the JST are implemented in the code. Currently, the code includes several
versions of the SA and k-ω SST models.

### B. Formulation of the Transition Model

The present formulation is summarized, emphasizing the differences in the Term2 and the way to
couple it to the baseline turbulence model. The original transport equation of Spalart Allmaras
turbulence model **without the $f_{t2}$ term** is given below:

$$
\frac{\partial \nu_T}{\partial t} + \frac{\partial}{\partial x_j}\!\left(\nu_T u_j\right)
= C_{b1} S\, \nu_T - C_{w1} f_w \left(\frac{\nu_T}{d}\right)^{2}
+ \frac{1}{\sigma}\left\{ \frac{\partial}{\partial x_j}\!\left[(\nu_L + \nu_T)\frac{\partial \nu_T}{\partial x_j}\right]
+ C_{b2}\,\frac{\partial \nu_T}{\partial x_j}\frac{\partial \nu_T}{\partial x_j} \right\}
\tag{1}
$$

The BC transition model (and the BCM model) relies on modification of the production term to
provide transition to turbulence. Here, the production term is multiplied by the intermittency
distribution function $\gamma_{BC}$:

$$
\frac{\partial \nu_T}{\partial t} + \frac{\partial}{\partial x_j}\!\left(\nu_T u_j\right)
= \boldsymbol{\gamma_{BC}}\, C_{b1} S\, \nu_T - C_{w1} f_w \left(\frac{\nu_T}{d}\right)^{2}
+ \frac{1}{\sigma}\left\{ \frac{\partial}{\partial x_j}\!\left[(\nu_L + \nu_T)\frac{\partial \nu_T}{\partial x_j}\right]
+ C_{b2}\,\frac{\partial \nu_T}{\partial x_j}\frac{\partial \nu_T}{\partial x_j} \right\}
\tag{2}
$$

where the intermittency function is given as:

$$
\gamma_{BC} = 1 - \exp\!\left(-\sqrt{Term1} - \sqrt{Term2}\right)
\tag{3}
$$

First, the Term1, which is the same in both the original and the modified model, is defined as:

$$
Term1 = \frac{\max\!\left(Re_\theta - Re_{\theta c},\; 0.0\right)}{\chi_1 Re_{\theta c}}
\tag{4}
$$

where the value of the calibration constant $\chi_1$ is **0.002** and the momentum thickness
Reynolds number along with the vorticity Reynolds number is given as:

$$
Re_\theta = \frac{Re_{v,\max}}{2.193}
\qquad\text{and}\qquad
Re_v = \frac{\rho\, d_w^{2}}{\mu}\,\Omega
\tag{5}
$$

$Re_{\theta c}$ is defined as the critical momentum thickness Reynolds number, which is a
correlation that is based on the data gathered from transition experiments and given as:

$$
Re_{\theta c} = 803.73\,(Tu_\infty + 0.6067)^{-1.027}
\tag{6}
$$

Second, the Term2 of the **original** BC model is:

$$
Term2 = \frac{\max\!\left(\nu_{BC} - \chi_2,\; 0.0\right)}{\chi_2}
\qquad\text{where}\quad
\nu_{BC} = \frac{\nu_t}{U d_w}
\quad\text{and}\quad
\chi_2 = \frac{5}{Re}
\tag{7}
$$

The local velocity magnitude $U$ appeared in the $\nu_{BC}$ term and the Re number appearing in
the calibration constant $\chi_2$ in Eq. 7 are removed in the **modified BCM model** as following:

$$
Term2 = \max\!\left(\chi_2 \frac{\nu_t}{\nu},\; 0.0\right)
\qquad\text{where}\quad \chi_2 = 50
\tag{8}
$$

The value of the calibration constant $\chi_2$ is found through numerical experimentation, and
the results of different $\chi_2$ constants are presented in the next section. A value of **50**
seems to fit the original model results more faithfully, again shown in the next section. The
constants given in Eq. 4 and Eq. 8 are calibrated against the Schubauer and Klebanoff data, but
obviously a range of data would be required to reach a wider application.

> **[Transcriber note — critical for code]** In Eq. 8 the constant $\chi_2 = 50$ appears **as a
> multiplier** of $\nu_t/\nu$, so the modified Term2 evaluates to $50\,(\nu_t/\nu)$. The companion
> paper (2020-2714) writes the *same* final term as $Term2 = \max\!\left(\tfrac{1}{\chi_2}\tfrac{\mu_T}{\mu},0\right)$
> with $\chi_2 = 0.02$, i.e. $\tfrac{1}{0.02} = 50$. **Both papers therefore agree that the
> implemented value is $50 \times$ (eddy-viscosity ratio).** Do not confuse this with the
> $\chi_2 = 5.0$ that appears in the *derivation* of the original model in 2020-2714. See the
> reconciliation notes in the 2020-2714 file.

---

## III. Validation

### I. Flat Plate Cases

In order to validate the modified BCM model, zero pressure gradient flat plate test cases are
carried out. Under this group, Schubauer & Klebanoff test case for natural transition and Savill's
T3A and T3B test cases for bypass transition were simulated. The free stream conditions are given
in Table 1.

**Table 1. Freestream conditions for the flat plate test cases**

| Case | $U_\infty$ (m/s) | $Re_\infty$ | $Tu$ (%) |
|---|---|---|---|
| Schubauer & Klebanoff | 50.1 | $3.4\times10^{6}$ | 0.18 |
| T3A | 5.4 | $3.7\times10^{5}$ | 3.0 |
| T3B | 9.4 | $6.3\times10^{5}$ | 6.0 |

Schubauer and Klebanoff case is selected for calibration as it represents a natural transition
process. The calibration is done by setting the model coefficients such that the transition onset
point matches with the experimental data and the original BC model results. The results are
presented for different values of $\chi_2$ constant in Figure 1 where excellent match with the
data and the original model is achieved for a value of 50.

> **Figure 1. Calibration case (Schubauer & Klebanoff) results.**
> $C_f$ vs $Re_x$ (0 to $4.8\times10^{6}$). Curves: BC Model; BCM Model $\chi_2 = 1$;
> BCM Model $\chi_2 = 50$; BCM Model $\chi_2 = 10000$; and S&K experiment (diamonds). The
> $\chi_2 = 50$ curve overlays the original BC model and the transition onset (rise in $C_f$ near
> $Re_x \approx 3\times10^{6}$) matches the experiment; $\chi_2 = 1$ and $\chi_2 = 10000$ shift the
> onset location.

Figure 2 compares the results obtained by the BC model and the BCM model with the experimental
data. The figure also includes the numerical results of the two-equation $\gamma - Re_\theta$
model of Menter et al., the most recent one-equation $\gamma$ model of Menter et al., and the
SA – $\gamma - Re_\theta$ model of Medida. Although there is an abrupt rise in the skin friction
and some overshoot in the transition region, same kind of abrupt rise behavior is observed in the
calculations of Menter et al. That is because, as discussed in Ref. 19, the BC and the BCM models
are abrupt models because of not having incorporated a transition length correlation into the
models, due to the observation that the effects of such correlations seem to be minimal.

> **Figure 2. Comparison of the skin friction coefficients for the S&K test case.**
> $C_f$ vs $Re_x$. Curves: BCM Model, BC Model, Menter $\gamma$-$Re_\theta$ Model, Menter
> $\gamma$ Model, Medida SA-$\gamma$-$Re_\theta$ Model, and S&K experiment. BCM and BC overlay each
> other; all models show transition near $Re_x \approx 3\times10^{6}$.

The second case in this group is selected to be the T3A flat plate test case. In Figure 3, it can
be seen that there is a perfect match between the results of the original BC model and the BCM
model. The figure also includes numerous results obtained by different transition models. In
general, the agreement between the numerical and the experimentally measured skin friction
coefficients is good. For this case, the BCM model and Menter et al. $\gamma$ model predict
somehow late transition onset, whereas Menter et al. $\gamma - Re_\theta$ model agrees better with
the experimental data.

> **Figure 3. Comparison of numerical and experimental data for T3A flat plate case.**
> $C_f$ vs $Re_x$ (0 to $6\times10^{5}$). Curves: BCM, BC, Menter $\gamma$-$Re_\theta$, Menter
> $\gamma$, Medida SA-$\gamma$-$Re_\theta$, Suzen & Huang, Walters $k$-$k_L$-$\omega$, plus T3A
> experiment (diamonds). Transition onset near $Re_x \approx 2\times10^{5}$.

The last flat plate case under this group is selected to be T3B which has the highest freestream
turbulence level. As shown in Figure 4, the transition onset is predicted somewhat late with the
BC and the BCM models whereas both Menter et al.'s models predict early onset of transition, but
could not capture the deepest point in the skin friction variation.

> **Figure 4. Comparison of numerical and experimental data for T3B flat plate case.**
> $C_f$ vs $Re_x$ (0 to $9\times10^{5}$). Same model set as Fig. 3 plus T3B experiment.

All in all, the newly formulated BCM model showed excellent agreement with the original BC model.
Also, the results obtained with the BCM model, along with the other models, are quite close to the
experimental data although there is some discrepancy in the overshoot regions for T3A and T3B
cases.

### II. 2-D Test Cases

#### a. Eppler E387 Airfoil Case

The Eppler E387 airfoil was tested by McGhee et al. at Langley Low-Turbulence Pressure Tunnel
(LTPT). In this study, numerical results are obtained with BCM transition model and SA turbulence
model for comparison with the experimental data at a low Re of **200,000** and freestream $Tu$ of
**0.1%**.

Figure 5 compares the experimental lift and drag coefficients with numerical results. As seen in
Fig. 5, drag coefficient prediction of BCM transition model is in good agreement with the
experiment, whereas the fully turbulent solution is rather off the data.

> **Figure 5. Comparison of experimental and numerical data obtained for the E387 airfoil at
> Re = 200,000.** Two panels: (left) $C_L$ vs AoA (−5° to 20°); (right) $C_L$ vs $C_D$ (drag
> polar). Curves: E387 Experiment (diamonds), BCM Model, Fully Turbulent SA Model. BCM tracks the
> experimental polar closely; fully-turbulent SA over-predicts drag.

The pressure coefficient distribution over the E387 airfoil for different angles of attack is also
compared against the experimental data and given in Fig. 6. BCM model seems to capture the
separation bubble remarkably well for all cases, whereas the fully turbulent solution, as
expected, could not do so.

> **Figure 6. Comparison of pressure coefficient distribution for different angles of attack for
> the E387 airfoil.** Three panels ($-C_p$ vs $x/c$) at 0° AoA, 2° AoA, 4° AoA. Curves: BCM Model,
> S-A Fully Turbulent, Experiment. BCM reproduces the laminar-separation-bubble plateau on the
> suction side.

#### b. T106 Turbine Cascade Case

For the second 2-D test case, the experimental study of Stieger et al. has been selected. In the
experiments, five-blade cascade of T106 profile was placed downstream of a moving bar wake
generator in order to investigate the interaction of a convected wake and a separation bubble on
the suction side of a low-pressure turbine blade. For this test case, the geometrical details are
given in Table 2, and the flow conditions correspond to a Re number of **91,000** (based on chord
length of the T106) with a freestream $Tu$ of **0.1%**.

**Table 2. T106 cascade geometrical details**

| Parameter | Value |
|---|---|
| Blade chord | 198 mm |
| Inlet flow angle | 37.7° |
| Blade stagger | 59.3° |
| Cascade pitch | 158 mm |
| Exit flow angle | 63.2° |
| Wake generator bar diameter | 2.05 mm |
| Axial distance from bars to leading edge of the T106 | 70 mm |

Figure 7 shows the experimental pressure coefficients over the T106 blade compared with the
numerical results of the BCM model and the Menter et al. $\gamma - Re_\theta$ model. Looking at
Figure 7, it can be said that both numerical models predict the separation bubble location
successfully; however, it is obvious that the predicted bubble size is slightly larger for the BCM
model. It is interpreted that this discrepancy might be caused by the differences between the
underlying turbulence models.

> **Figure 7. Pressure coefficient distribution over the T106 profile.**
> $C_p$ vs $s/s_{\max}$ (0 to 1). Curves: BCM Model, Menter $\gamma$-$Re_\theta$, T106A experiment
> (diamonds). Both models reproduce the suction-side separation bubble; BCM predicts a slightly
> larger bubble.

### III. 6:1 Prolate Spheroid Case

The last test case is selected to be the experimental work of Kreplin et al., which investigates
incompressible flow over a 6:1 prolate spheroid. The 6:1 prolate spheroid in the experiments have
a minor and major axis lengths of **0.4 m** and **2.4 m**, respectively.

In this study, the test case with a flow velocity of **45 m/s**, an angle of attack of **5
degrees**, and a freestream $Tu$ of **0.1%** is considered. The results are obtained using a coarse
grid of 100×100 nodes over the prolate spheroid and a fine grid of 300×300 nodes. A view of the
coarse grid is shown in Figure 8.

> **Figure 8. Surface grid of the 6:1 prolate spheroid test case (100×100 nodes coarse grid).**
> 3-D rendering of the structured surface mesh on the elongated spheroid.

Figure 9 compares the experimental skin friction coefficients along the top X-Z cutting plane of
the geometry to the numerically obtained values with the BCM model and the Menter et al.
$\gamma - Re_\theta$ model. It is observed from the results that the BCM model predicts an early
transition for the coarse grid whereas the $\gamma - Re_\theta$ model shows a somewhat late
transition. However, the fine grid improved the results of the BCM Model; prediction of the
transition location is very close to that of the experiment. The comparison of the skin friction
contours of the BCM model and the $\gamma - Re_\theta$ model shown in Figure 10 seem quite similar
except that the solution using the coarse mesh with the BCM model.

> **Figure 9. Skin friction coefficients along the top X-Z cutting plane over the 6:1 prolate
> spheroid.** $C_f$ vs $x$ (m) (0.6 to 1.8). Curves: BCM Model Fine Grid, BCM Model Coarse Grid,
> Menter $\gamma$-$Re_\theta$, Experiment (diamonds). Coarse-grid BCM transitions early
> ($x \approx 0.9$ m); fine-grid BCM matches the experimental onset ($x \approx 1.2$ m).

> **Figure 10. Comparison of the skin friction coefficient contours over the 6:1 prolate
> spheroid.** Three surface contour plots of $C_f$ (colorbar 0 to 0.0042): $\gamma - Re_\theta$
> Model; BCM Model Coarse Grid; BCM Model Fine Grid.

---

## IV. Conclusion

A revised version of the algebraic BC transition model is proposed. In the BCM transition model,
first, the lack of Galilean invariance in the original model due to inclusion of the local
velocity magnitude in the equations is removed. Second, the appearance of the Re number in
formulation of the original model is eliminated, which before caused some problems in a way that
it brings a possibly non-unique reference length into play. The BCM model is employed to solve
zero pressure gradient flat plate test cases, an airfoil test case and a turbine cascade case with
success. Finally, the flow over a 6:1 prolate spheroid is investigated which shows that the BCM
model is quite reliable for the 3-D flow cases even if it solves less equations than its
counterparts.

---

## References

1. Spalart, P. R. and Allmaras, S. R., "A One-Equation Turbulence Model for Aerodynamic Flows," AIAA-92-0439, 1992, Reno, NV.
2. Menter, F. R., "Two-equation eddy-viscosity turbulence model for engineering applications," *AIAA Journal*, Vol. 32, 1994, pp. 1598–1605.
3. Menter, F. R., Langtry, R. B., Likki, S. R., Suzen, Y. B., Huang, P. G., Volker, S., "A Correlation based Transition Model using Local Variables Part 1 - Model Formulation," ASME Paper No. GT-2004-53452.
4. Menter, F. R., Esch, T., Kubacki, S., "Transition Modelling based on Local Variables," In: Proc. 5th Int. Sym. on Engineering Turbulence Modelling and Measurements, 2002, Mallorca, Spain.
5. Menter, F. R., Langtry, R. B., Volker, S., "Transition Modelling for General Purpose CFD Codes," *Flow, Turbulence and Combustion*, Vol. 77, 2006, pp. 277-303.
6. Langtry, R. B., Menter, F. R., "Correlation-Based Transition Modeling for Unstructured Parallelized Computational Fluid Dynamics Codes," *AIAA Journal*, Vol. 47, No. 12, 2009, pp. 2984-2906.
7. Wilcox D. C., "Turbulence modeling for CFD," 1st ed. La Canada, CA: DCW Industries, 1993.
8. Content, C., Houdeville, R., "Application of the γ-Reθ Laminar-Turbulent Transition Model in Navier-Stokes Computations," AIAA-2010-4445, 2010.
9. Suluksna, K., Juntasaro, V., Juntasaro, E., "Capability Assessment of Intermittency Transport Equations for Modeling Flow Transition," In: Proc. 19 Conference of Mechanical Engineering Network of Thailand, 2005, Phuket, Thailand.
10. Malan, P., Suluksna, K., Juntasaro, E., "Calibrating the γ-Reθ Transition Model for Commercial CFD," AIAA 2009-1142, 2009.
11. Misaka, T., Obayashi, S., "Application of Local Correlation–Based Transition Model to Flows around Wings," AIAA Paper 2006-918.
12. Piotrowski, W., Elsner, W., Drobniak, S., "Transition Prediction on Turbine Blade Profile with Intermittency Transport Equation," ASME J. of Turbomach., Vol. 132(1), 2009.
13. Grabe, C., Krumbein, A., "Extension of the γ-Reθ model for prediction of crossflow transition," AIAA-2014-1269.
14. Dassler, P., Kozulovic, D., Fiala, A., "An approach for modelling the roughness-induced boundary layer transition using transport equations," 6th European congress on computational methods in applied sciences and engineering (ECCOMAS), 2014.
15. Kaynak, U., "Supersonic boundary-layer transition prediction under the effect of compressibility using a correlation-based model," Proc IMechE, Part G: J Aerospace Engineering, Vol. 226, 2011, pp. 722–739.
16. Walters, D. K., Leylek, J. H., "A New Model for Boundary-Layer Transition Using a Single-Point RANS Approach," ASME J. of Turbomach., Vol. 126(1), 2004, pp. 193–202.
17. Walters, D. K., Cokljat, D., "A Three-Equation Eddy-Viscosity Model for Reynolds-Averaged Navier-Stokes Simulations of Transitional Flows," J. of Fluids Eng., Vol. 130, 2008.
18. Bas, O., Cakmakcioglu, S. C., Kaynak, U., "A novel intermittency distribution based transition model for low-Re number airfoils," AIAA-2013-2531, 2013.
19. Cakmakcioglu, S. C., Bas, O., Kaynak, U., "A correlation-based algebraic transition model," *Proc IMechE, Part C J Mechanical Engineering Science*, Vol. 232(21), 2018, pp. 3915-3929.
20. Menter, F. R., Smirnov, P. E., Liu, T., Avancha, R., "One-equation local correlation-based transition model," *Flow, Turbulence and Combustion*, Vol. 95(4), 2015, pp. 583-619.
21. Kubacki, S., Dick, E., "An algebraic model for bypass transition in turbomachinery boundary layer flows," *Int J Heat Fluid Flow*, Vol. 58, 2016, pp. 68–83.
22. Schubauer, G. B., Klebanoff, P. S., "Contribution on the mechanics of boundary layer transition," NACA Technical Note No. TN-3489, 1955.
23. Savill, A.M., "Some Recent Progress in Turbulence Modeling of By-pass Transition," In: So, R.M.C., Speziale, C.G. and Launder, B.E. (eds) Near-Wall Turbulent Flows, Elsevier, 1993, pp.829-848.
24. McGhee, R. J., Walker, B. S., Millard, B. F., "Experimental results for the Eppler 387 airfoil at low Reynolds numbers in the Langley low-turbulence pressure tunnel," NASA Technical Memorandum No. 4062 1988.
25. Stieger, R., Hollis, D., Hodson, H., "Unsteady surface pressures due to wake induced transition in a laminar separation bubble on a LP turbine cascade," ASME paper no. GT2003-38303, 2003.
26. Kreplin, H. P., Vollmers, H., Meier, H., "Wall shear stress measurements on an inclined prolate spheroid in the DFVLR 3m x 3m low speed wind tunnel," Technical report, DFVLR-AVA.
27. Palacios, F., Alonso, J., Duraisamy, K., Colonno, M., Hicken, J., Aranake, A., Campos, A., Copeland, S., Economon, T., Lonkar, A., Lukaczyk, T., Taylor, T., "Stanford University Unstructured (SU2): An Open-Source Integrated Computation Environment for Multi-Physics Simulation and Design," AIAA Paper No. 2013-0287.
28. Jameson, A., and Seokkwan, Y., "Lower-Upper Implicit Schemes with Multiple grids for the Euler Equations," *AIAA Journal*, Vol. 25(7), 1987, pp. 929-935.
29. Van der Vorst, H. A., "Bi-CGSTAB: A Fast and Smoothly Converging Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems," *SIAM J. Sci. and Stat. Comput.*, Vol. 13(2), 1992, pp.631-644.
30. Saad, Y., Schultz, M. H., "GMRES: A Generalized Minimal Residual Algorithm for Solving Nonsymmetric Linear Systems," *SIAM J. Sci. and Stat. Comput.*, Vol. 7(3), 1986, pp. 856-869.
31. Roe, P. L., "Approximate Riemann Solvers, Parameter Vectors and Difference Schemes," *Journal of Computational Physics*, Vol. 43(2), 1981, pp. 357-372.
32. Jameson, A., Schmidt, W., Turkel, E., "Numerical Solution of the Euler Equations by Finite Volume Methods Using Runge-Kutta Time Stepping Schemes," AIAA Paper No. 81-1259.
33. Medida, S., "Correlation-based transition modeling for external aerodynamic flows," PhD Thesis, University of Maryland, USA, 2014.
34. Suzen, Y. B., Huang, P. G., "Modeling of Flow Transition Using an Intermittency Transport Equation," *ASME J. Fluid Eng.*, Vol. 122(2), 2000, pp. 273-284.
