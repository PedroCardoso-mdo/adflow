# A Revised One-Equation Transitional Model for External Aerodynamics — Part I: Theory, Validation and Base Cases

> **Reference transcription for code implementation.** Equations, constants, boundary
> conditions and numerical settings are transcribed verbatim from the source PDF. Figures are
> validation results and are given as captions plus a short description. Editorial remarks added
> during transcription are marked **[Transcriber note]** and are *not* part of the paper.
> **The definitive, implementation-ready formulation is the Appendix** (complete SA-BCM including
> the negative continuation and the full constant set). A cross-paper **Implementer's
> reconciliation notes** section is appended at the very end.

| | |
|---|---|
| **Title** | A Revised One-Equation Transitional Model for External Aerodynamics – Part I: Theory, Validation and Base Cases |
| **Authors** | Riccardo Mura (Università degli Studi di Cagliari, Cagliari, Italy — PhD student in Industrial Engineering); Samet C. Cakmakcioglu (Turkish Aerospace, Ankara, Turkey — Flight Mechanics and Performance Chief Engineer) |
| **Venue** | AIAA Aviation Forum, June 15–19, 2020, Virtual Event — AIAA Aviation 2020 Forum |
| **DOI** | 10.2514/6.2020-2714 |
| **Model** | SA-BC → **SA-BCM** (revised); baseline = Spalart–Allmaras with Oliver's negative continuation (SA-neg) |
| **Solver** | ANSYS Fluent 17 (via UDF) |

---

## Abstract

The SA-BC (Spalart-Allmaras Bas-Cakmakcioglu) one-equation transitional model has been revised.
Weaknesses of the original formulation are overcome removing the lack of Galilean invariance and
the need for the inclusion of an ambiguous length scale in a term inside the intermittency
function, which is responsible for the growing of the intermittency and then of the turbulence
production inside a laminar boundary layer. After choosing a robust implementation for the baseline
Spalart-Allmaras model, the resulting transitional model has been implemented in the commercial
solver ANSYS Fluent. Calibration against flat-plate test cases makes the revised model behave as
the original one without loss of accuracy. Computations on airfoils at low Reynolds numbers show
fair agreement with experimental data and the results are comparable to those obtained with a more
complex three-equations transitional model.

---

## I. Introduction

An algebraic transition model coupled to the Spalart-Allmaras (SA) turbulence closure was firstly
presented and then, in its final formulation, in Cakmakcioglu et al. Although the good results
obtained with the model, some confusion has arisen about the formulation, leading to a series of
reports and warnings on the well known TMR (Turbulence Modeling Resource) web page: the formulation
there presented matched the equations, but in consequence of some typos they were both different
from the correct sources. The latter were however publicly available because of the "official"
implementation of the model in the open-source code SU2. Beside the difficulties in replicating the
results presented in the aforementioned references, some critical issues were also highlighted on
the TMR web page, mainly regarding the lack of Galilean invariance — i.e., the equation depends on
the reference frame due to the inclusion of the local fluid velocity —, and the presence of a
non-unique length scale to be used for computing the freestream Reynolds number. Moreover, non-local
quantities were employed in the formulation, what is contrary to principles of the Local
Correlation-Based Transition Modeling (LCTM) framework, to which the transitional model belongs.
The present work demonstrates that the two shortcomings can be easily removed by simplifying the
problematic term, without altering the performance offered by the original model significantly. The
revised model is reasonably successful in predicting low Reynolds number flows pertaining to
external aerodynamics cases. The results obtained with a particularly robust implementation in
ANSYS Fluent 17 are presented for three basic test-cases after discussing the revision, while a
companion paper deals with more advanced cases and demonstrates that only minimal modifications are
needed for the original SU2 implementation in order to obtain consistent results.

> Footnote: In the context of RANS turbulence closures, the term *transition model* usually denotes
> a separate model that deals with transition and applies some modification to the baseline full
> turbulent model, either with or without a set of additional equations. This is not always the
> case, and some turbulence closures were generalized so to natively incorporate modeling of
> transitional as well as fully turbulent flows. In this paper, the term *transitional model* is
> employed to denote any complete closure capable of dealing with transition, regardless of whether
> the transition model is "appended" to some existent fully turbulent model or incorporated into a
> set of original equations.

---

## II. Revision of the *Term2* term

Generalizing the production ($P_{\tilde{\nu}}$) and destruction ($D_{\tilde{\nu}}$) terms, and the
diffusivity coefficient ($\Gamma_{\tilde{\nu}}$), the Spalart-Allmaras Bas-Cakmakcioglu (SA-BC)
transitional model has the following formulation:

$$
\begin{cases}
\mu_T = \rho\,\tilde{\nu}\,f_{v1} \\[4pt]
\dfrac{\partial(\rho\tilde{\nu})}{\partial t} + \dfrac{\partial(\rho\tilde{\nu}U_j)}{\partial x_j}
= \gamma\,\rho\,P_{\tilde{\nu}} - \rho\,D_{\tilde{\nu}}
+ \dfrac{1}{\sigma}\left[\dfrac{\partial}{\partial x_j}\!\left(\Gamma_{\tilde{\nu}}\dfrac{\partial\tilde{\nu}}{\partial x_j}\right)
+ \rho\,C_{b2}\dfrac{\partial\tilde{\nu}}{\partial x_j}\dfrac{\partial\tilde{\nu}}{\partial x_j}\right] \\[10pt]
\gamma = 1 - e^{-\left(\sqrt{Term1} + \sqrt{Term2}\right)}
\end{cases}
\tag{1}
$$

$\mu_T$ is the turbulent (or eddy) viscosity used for the closure of the Reynolds-Averaged
Navier-Stokes (RANS) equation. $\tilde{\nu}$ is the working variable transported by the PDE that
constitutes the baseline SA model, called modified turbulent viscosity. $U_j$ are the components of
the local mean flow velocity. $\rho$ is the density of the fluid. $C_{b2}$ is a calibration
constant. The function $f_{v1}$ is given by

$$
f_{v1} = \frac{\chi^3}{\chi^3 + C_{v1}^3}
\tag{2}
$$

where $\chi = \tilde{\nu}/\nu$ is the ratio between the modified turbulent viscosity and the
kinematic viscosity $\nu$ of the fluid. *Term1* is mainly responsible for the estimation of the
so-called intermittency $\gamma$, which determines the fraction of time the flow is turbulent at a
given place, in a statistical sense. The term is expressed as

$$
Term1 = \frac{\max\!\left(\dfrac{Re_v}{2.193} - Re_{\theta c},\; 0\right)}{\chi_1 Re_{\theta c}}
\tag{3}
$$

where $Re_{\theta c}$ is the critical value of the momentum thickness Reynolds number $Re_\theta$,
while $Re_v$ is the vorticity Reynolds number computed as

$$
Re_v = \frac{\rho\, d_w^2}{\mu}\,S
\tag{4}
$$

$d_w$ is the distance from the nearest solid wall. $\mu$ is the dynamic viscosity of the fluid. The
quantity $S$ **must not be confused with the mean strain rate tensor magnitude**: it refers to the
mean vorticity tensor $\Omega_{ij} = 0.5(\partial_j U_i - \partial_i U_j)$, whose magnitude is
$S = \Omega = \sqrt{2\,\Omega_{ij}\Omega_{ij}}$. The empirical correlation relating the critical
momentum thickness Reynolds number to the freestream turbulent intensity $Tu_\infty$ is given by

$$
Re_{\theta c} = 803.73\,(Tu_\infty + 0.6067)^{-1.027}
\tag{5}
$$

Until now, *Term2* was known to be

$$
Term2 = \frac{\max(\nu_{BC} - \chi_2,\; 0)}{\chi_2}
\tag{6}
$$

where:

$$
\begin{cases}
\chi_2 = \dfrac{5.0}{Re} \\[6pt]
\nu_{BC} = \dfrac{\nu_T}{U d_w}
\end{cases}
\tag{7}
$$

$Re$ is the reference Reynolds number of the mean flow. $\nu_T$ is the turbulent viscosity per unit
density. $U$ is the mean velocity magnitude. Consistently with the SU2 implementation, Eqs. 6 and 7
should be rather expressed in an equivalent form, given by

$$
\begin{cases}
Term2 = \dfrac{\max(\nu_{BC} - \nu_{cr},\; 0)}{\nu_{cr}} \\[6pt]
\nu_{cr} = \dfrac{\chi_2}{Re} \\[6pt]
\chi_2 = 5.0
\end{cases}
\tag{8}
$$

Inclusion of *Term2* inside the negative exponential in the intermittency function (last of Eqs. 1)
aims at allowing the growth of intermittency and then the turbulence production inside a laminar
boundary layer: without including *Term2*, once $\gamma$ has gone to zero it cannot grow anymore,
regardless of the values assumed by *Term1*, which becomes therefore ineffective. According to
Cakmakcioglu et al., $\chi_2$ should be interpreted as the critical value of the ratio $\mu_T/\mu$
to be exceeded in order to activate the turbulence production. However, expansion of the *Term2*
expression rather reveals its nature as a critical value of the ratio $Re/Re_T$, where $Re_T$ can be
seen an additional Reynolds number given by the inverse of the $\nu_{BC}$ quantity:

$$
\nu_{BC} = \frac{\nu_T}{U d_w} = \frac{1}{Re_T}
\tag{9}
$$

*Term2* can therefore be written as

$$
Term2 = \frac{\max(1/Re_T - \chi_2/Re,\; 0)}{\chi_2/Re}
= \max\!\left(\frac{1}{\chi_2}\frac{Re}{Re_T} - 1,\; 0\right)
\tag{10}
$$

When $Re/Re_T$ is greater than $\chi_2$, *Term2* assumes positive values, causing the exponential
function to decay and the intermittency value to grow. This expression of the model should not
depend on *a priori* specifications of the Reynolds number. This expression presents at least the
following problematic points:

1. The Reynolds number $Re$ accounts for the scales of the freestream flow, whose relation with the
   local turbulence behavior is not so simple that a comparison of $Re$ and $Re_T$ ($\chi_2 Re_T$)
   could be done based on unambiguous arguments.
2. While a value for the fluid viscosity is certainly known and a velocity scale is ever
   identifiable (the freestream value), a unique length scale cannot always be defined, as it is
   easy to demonstrate referring to the flat plate case, for which $Re$ varies linearly with the
   distance from the leading edge. Several test computations (not reported here) have shown that the
   model behavior is not significantly affected by the choice of the length scale $L$. Nevertheless,
   if the quantity has not a real impact neither should it be included, since it only generates
   additional ambiguities.
3. Including freestream or however chosen non-local quantities inside the formulation is a violation
   of the LCTM essentials. Results provided by the model should not depend on *a priori*
   specifications of the Reynolds number.
4. The inclusion of the local velocity $U$ violates the Galilean invariance, limiting the model
   generality. Although this shortcoming is also present, for instance, in the $\gamma$-$Re_\theta$-SST
   model, it should be stressed that Galilean invariance actually represents a requirement for a
   well-founded turbulence model, and such a fundamental requirement should not be relaxed if not
   strictly needed. It can also be noted that when the $\gamma$-SST model was formulated in an
   attempt to simplify the $\gamma$-$Re_\theta$-SST formulation, attention was paid to restore the
   Galilean invariance.

First of all, the unitary bounding value on the RHS of Eq. 10 can be removed, since also when
turbulence levels are low ($\nu_T < \nu$), $U_\infty L \gg U d_w$ is expected to hold, ensuring
$Re/(\chi_2 Re_T) > 1$. It has been said that inclusion of the mean flow length scale $L$ basically
gives raise to ambiguities when specifying the inputs for the model, without any advantages in terms
of prediction capabilities. In order to keep dimensional consistency, distance from the nearest
wall, $d_w$, is also removed. The latter choice casts a doubt regarding the possibility of a
solution degradation in consequence of neglecting the effects of wall proximity: the results will
show that no significant differences are detected. The ratio $U_\infty/U$ is also ignored, since its
inclusion would be, at this point — at the very least — quite arbitrary. The term is finally reduced
to

$$
Term2 = \max\!\left(\frac{1}{\chi_2}\frac{\nu_T}{\nu},\; 0\right)
= \max\!\left(\frac{1}{\chi_2}\frac{\mu_T}{\mu},\; 0\right)
\tag{11}
$$

> **[Transcriber note — CRITICAL for code]** In Eq. 11 the calibration constant $\chi_2$ appears as
> $1/\chi_2$. The value $\chi_2 = 5.0$ quoted in Eqs. 7–8 belongs to the *original* model
> derivation (where $\chi_2$ is divided by $Re$). The **calibrated value for the revised Eq. 11 is
> $\chi_2 = 0.02$** — stated in §III.D and, definitively, in the **Appendix constant table**. Hence
> $1/\chi_2 = 1/0.02 = 50$ and $Term2 = \max\!\left(50\,\frac{\mu_T}{\mu},\,0\right)$, which matches
> Eq. 8 ($\chi_2 = 50$, multiplier form) of the companion paper 2020-2706. **Using $\chi_2 = 5.0$ in
> Eq. 11 would make Term2 too small by a factor of 250** and break the calibration.

Bounding the term to assume positive values is unnecessary, physically speaking, because the eddy
viscosity cannot be negative. However, the bound is retained for the sake of generality. Standard
(positive) implementations of the SA model can generate ramp solutions, for instance at the edge of
boundary layers and wakes. More severe ramp solutions are likely to arise if turbulence is damped
out to very low levels in the outer region. Such ramp solutions are responsible for the undershoots.
The problem is of particular concern with the SA-BC model, since the intermittency is always zero
far from solid boundaries and outside of wakes. The use of the revised model is intended to be
extended to more challenging computations than those discussed in the remainder of this paper. It
should be noted that the bound in Eq. 11 is not enough in order to mitigate the problem, since
erroneous solutions and numerical instabilities can arise due to the erratic behavior of the
transport equation for $\tilde{\nu}$. The question motivated the choice of a proper implementation in
order to overcome these difficulties. The valuable information provided by Allmaras et al. was used
as a starting point for finding a good strategy, which is one of the subjects of the following
section.

---

## III. Implementation and calibration

The SA model was formulated without any real transition capabilities. However, it was adjusted in
order to provide some flexibility for the treatment of laminar boundary layers and their transition
to the turbulent state: a source term, namely the *trip term*, allowed to directly specify a
transition point in the domain. Two quantities altering the production and destruction terms
(laminar suppression terms) were employed so to preserve the laminar boundary layer until the
prescribed transition point was reached. Such terms are usually not included in modern
implementations, which rather employ the model as a fully turbulent one. The older implementation of
the SA-BC model included the laminar suppression ($f_{t2}$) terms (without trip term), but they were
removed in the following revision, which also presents different calibration constants. Note that
the laminar suppression terms can constrain the turbulence production - forcing the boundary layer to
remain laminar - even when intermittency equals one, altering not only the transition behavior but
also the fully turbulent portions of the solution. The problem is well known and arises when low
$\tilde{\nu}/\nu$ ratios are specified at the inlet of the domain. When dealing with transition,
such values are usually low to very low. It should be stressed that within the LCTM paradigm to
avoid prediction of a laminar flow by the baseline model is a fundamental requirement in order to
ensure proper behavior of the coupled transition model. For the implementation of the revised
transitional model, employment of the laminar suppression terms was carefully avoided, while other
possible issues where removed beforehand following with the choices outlined below.

### A. Positive SA

For a positive, fully turbulent implementation of the SA model, the source term in the second of
Eqs. 1 simply goes to zero, while the production and destruction terms reduce to

$$
\begin{cases}
P_{\tilde{\nu}} = C_{b1}\tilde{S}\tilde{\nu} \\[4pt]
D_{\tilde{\nu}} = C_{w1}f_w\left(\dfrac{\tilde{\nu}}{d_w}\right)^2
\end{cases}
\tag{12}
$$

where $C_{b1}$ and $C_{w1}$ are model constants. In order to avoid numerical problems, the modified
vorticity $\tilde{S}$ must always be positive. Following a modification proposed in Allmaras et al.
and directed to the standard SA model, the quantity can be computed as

$$
\tilde{S} =
\begin{cases}
S + \overline{S} & \overline{S} \geq -C_{v2}S \\[6pt]
S + \dfrac{S\left(C_{v2}^2 S + C_{v3}\overline{S}\right)}{(C_{v3} - 2C_{v2})S - \overline{S}} & \overline{S} < -C_{v2}S
\end{cases}
\tag{13}
$$

where:

$$
\overline{S} = \frac{\tilde{\nu}}{\kappa^2 d_w^2}\,f_{v2}
\tag{14}
$$

$\kappa$ is the Von Kármán constant. Note that Eq. 13 does not belong to the standard SA
implementation: it was proposed in Allmaras et al. — and previously referred to it —, with the
purpose of avoiding negative values of the modified vorticity $\tilde{S}$ arising from a range of
$\chi = \tilde{\nu}/\nu$ values. The modification, also documented in the TMR web page section
dedicated to the SA model, ensures positive $\tilde{S}$ values for nonzero $S$ values and makes the
quantity $C^1$ continuous.

The function $f_w$ has the following expression:

$$
f_w = g\left(\frac{1 + C_{w3}^6}{g^6 + C_{w3}^6}\right)^{1/6}
\tag{15}
$$

where

$$
\begin{cases}
g = r + C_{w2}(r^6 - r) \\[4pt]
r = \dfrac{\tilde{\nu}}{\tilde{S}\kappa^2 d_w^2}
\end{cases}
\tag{16}
$$

and $C_{w2}$ is another calibration constant. The function $f_{v2}$ is defined as:

$$
f_{v2} = 1 - \frac{\chi}{1 + \chi f_{v1}}
\tag{17}
$$

### B. Negative continuation

The SA model is known to be robust. However, it can suffer for locally negative (non-physical)
solutions of the discretized equations. Many implementations rely on explicit clipping, i.e., the
negative $\tilde{\nu}$ values are removed from the solution so ensure that the eddy viscosity is
non-negative everywhere in the domain. Since this strategy can negatively affect convergence of the
discretized residuals, some negative continuation of the model should be formulated instead. A
negative prosecution implies a reformulation of the PDE to be used when the transported variable
assumes negative values. After highlighting the deficiencies of alternative approaches, some
negative continuations are formulated for a SA negative continuation in Allmaras et al., and a
proposal is made, correspondent to the SA-neg model variant. Although the proposal respects all of
the prescribed requirements, it is tied to the inclusion of the $f_{t2}$ terms, which neglect the
coupling with an LCTM transition model. An early variant, without the laminar suppression terms, was
discussed by Oliver in his PhD thesis. Oliver's continuation — as it will be referred to from now on
— was immediately considered for supporting the implementation of the transitional model. It is
respectful of the aforementioned requirements, which are reported here for convenience:

1. the SA formulation must remain unchanged for $\tilde{\nu} \geq 0$;
2. $\mu_T = 0$ for $\tilde{\nu} < 0$;
3. $C^1$ continuity must be ensured at $\tilde{\nu} = 0$ for all the functions in the PDE;
4. the negative continuation must be energy stable, meaning that $\frac{1}{2}\tilde{\nu}^2 < 0$ for
   $\tilde{\nu} < 0$ (the energy associated with the variable cannot increase when the latter is
   negative);
5. non-negative analytic solutions for non-negative boundary conditions.

According to the Oliver's continuation, the production term becomes

$$
P_{\tilde{\nu}} =
\begin{cases}
C_{b1}\tilde{S}\tilde{\nu} & \tilde{\nu} \geq 0 \\[4pt]
C_{b1}S\,g_n\,\tilde{\nu} & \tilde{\nu} < 0
\end{cases}
\tag{18}
$$

where the function $g_n$ is defined by

$$
g_n = 1 - \frac{10^3\chi^2}{1 + \chi^2}
\tag{19}
$$

Accordingly, the destruction term is given by

$$
D_{\tilde{\nu}} =
\begin{cases}
C_{w1}f_w\left(\dfrac{\tilde{\nu}}{d_w}\right)^2 & \tilde{\nu} \geq 0 \\[8pt]
-C_{w1}\left(\dfrac{\tilde{\nu}}{d_w}\right)^2 & \tilde{\nu} < 0
\end{cases}
\tag{20}
$$

Finally, the diffusivity is defined as

$$
\Gamma_{\tilde{\nu}} =
\begin{cases}
\mu + \rho\tilde{\nu} & \tilde{\nu} \geq 0 \\[6pt]
\mu + \rho\tilde{\nu} + \dfrac{1}{2}\dfrac{\rho^2\tilde{\nu}^2}{\mu} & \tilde{\nu} < 0
\end{cases}
\tag{21}
$$

### C. Solver and numerical settings

The SA-BC model was implemented in the cell-centered finite volume commercial CFD solver ANSYS
Fluent 17 via UDF. In addition to the modifications described above, the implementation of the
baseline SA model further differs from the one shipping with the solver for adopting a different
value of the Von Kármán constant — **0.41 as in the original SA model references, against the
solver's default value of 0.4187**. The model that will be employed for comparison, i.e., the
$\gamma$-SST models, follows the standard implementations provided by the solver, as well as the
default settings. In particular, it is worth noting that two different production limiters for the
turbulent kinetic energy are active by default, including the one proposed by Kato and Launder. No
rotation or curvature corrections were employed. All the results presented in this work are obtained
using the same numerical settings. The pressure-based approach with full coupling of momentum and
pressure-correction equations was chosen since only incompressible computations were performed: the
choice was verified to be effective in reducing the overall computational time compared to the
solution of the equations in a segregated fashion. The nominally second-order upwind scheme based on
linear reconstruction of the gradients was used for the convective terms of all equations. Gradients
were discretized with the Least Squares method. In Menter et al. it is emphasized that the transition
onset location is sensitive to the advection scheme used for transporting the turbulent quantities,
and the use of the same upwind scheme is reported to be successful in achieving grid-independent
solutions. Although the advice was the for $\gamma$-$Re_\theta$-SST model, a similar behavior is
expected also for other LCTM models. A second-order central differencing scheme was employed for
pressure interpolation. Iterative convergence was judged on the base of relative changes in the force
coefficients, and a target threshold of $10^{-6}$ for the scaled absolute residuals was found to be
effective in order to achieve full convergence with the employed grids.

> **Figure 1. Flat plate mesh (detail). The entire mesh is extend up to $x/L = 2$.**
> Structured mesh in the $(x/L, x/L)$ plane near the leading edge, refined toward the wall.

### D. Calibration

In order to assign to $\chi_2$ a proper value, a new calibration was needed, and some flat plate test
cases were used for the purpose. Beside the Schubauer & Klebanoff test case, the ERCOFTAC T3A, T3A-
and T3B tests, often employed to calibrate or benchmark transition models, were employed. The aim was
to match the original SA-BC model results rather than improving them: the revised formulation
presently aims at overcoming the main deficiencies of the model, while improvements in terms of
predictive capabilities are left to future works.

The flat plate case based on the Schubauer & Klebanoff test was chosen by the SA-BC model's authors
for its calibration. The case deals with natural transition at the relatively large Reynolds number
of $3.34 \times 10^6$, under a freestream turbulence intensity equal to 0.18%. The computational grid
essentially follows those used in the references, but the stream-wise cell density along the plate
was adjusted in order to achieve grid-independent results. From the leading edge, constant spacing is
assigned to cells so to accommodate different transition locations for the other flat-plate test
cases: this results in a number of cells greater that the one strictly needed for each case, but on
the other hand it removes the need to use different grids. An expansion ratio of 1.05 together with a
geometric law was adopted in the wall-normal direction for optimal cell growth, while the height of
the first cell was adjusted to ensure a non-dimensional wall distance $y^+ < 1$: these settings will
be used for all the grids in the following, therefore the details about the cell distribution normally
to the walls will not be repeated. Since the Schubauer & Klebanoff test is the one involving the
highest Reynolds numbers, the constraint $y^+ < 1$ is satisfied *a fortiori* for all other tests. The
grid is shown in Fig. 1.

A constant velocity field is applied at the inlet, where the ratio $\nu' = \tilde{\nu}/\nu = 0.02$ is
also specified, choosing the central value in the range $[0.015 - 0.025]$ suggested in the references.
The ratio is there expressed as $\tilde{\nu}/\nu_\infty$, which is the correct specification in the
general (compressible) case: here $\nu$ assumes a prescribed value for the entire domain in
consequence of the incompressible assumption, and then the subscript $\infty$ is dropped. At the
outlet, an outflow boundary condition, i.e., zero order extrapolation, is used. For the lower boundary
before the plate (negative $x$ values), a symmetry conditions is imposed, while for the upper boundary
a free-slip boundary condition is assigned. The canonical no-slip wall condition is used for the plate
surface ($x/L \in [0, 2]$). **$\chi_2 = 0.02$ was found to be an optimal value** for minimizing the
difference between the SA-BC and the revised SA-BCM (Spalart-Allmaras Bas-Cakmakcioglu with
Modifications) models. The results are shown in Figs. 2, 3, 4 and 5. They show almost a perfect
agreement between the old the revised models. Note that bad behavior for the T3A- test was already
shown in the references, and it is not specific to the revised version, nor it depends on the mesh or
boundary conditions.

> **[Transcriber note]** This confirms the implemented calibration constant $\chi_2 = 0.02$ used in
> Eq. 11 / Appendix, distinct from the $\chi_2 = 5.0$ of Eqs. 7–8. Inlet condition for the SA
> working variable: $\tilde{\nu}/\nu = 0.02$ (range 0.015–0.025).

> **Figure 2. Skin Friction coefficient distribution for the Schubauer & Klebanoff flat plate
> test.** $C_f$ vs $Re_x$ (0 to $4.8\times10^6$). Curves: SA ($\nu' = 5$); SA–BC
> ($Tu = 0.18\%$, $\nu' = 0.02$); SA–BCM ($Tu = 0.18\%$, $\nu' = 0.02$); Experiments (squares).
> SA–BC and SA–BCM overlay; transition onset near $Re_x \approx 2.9\times10^6$.

In the following section, the SA-BCM model will be compared to the $\gamma$-SST model, considered here
the most appropriate competitor. Three cases dealing with low Reynolds number airfoils will be used
for the purpose. The results computed with the respective fully turbulent baseline models will also be
included in order to asses the capabilities of the LCTM models in this category of flows. The results
computed with the SA-BC model will be included as well, so to validate the SA-BCM model calibration.

> **Figure 3. Skin Friction coefficient distribution for the ERCOFTAC T3A flat plate test.**
> $C_f$ vs $Re_x$ (0 to 600,000). Curves: SA ($\nu' = 5$); SA–BC ($Tu = 3\%$, $\nu' = 0.02$);
> SA–BCM ($Tu = 3\%$, $\nu' = 0.02$); Experiments. Transition onset near $Re_x \approx 200{,}000$.

> **Figure 4. Skin Friction coefficient distribution for the ERCOFTAC T3A- flat plate test.**
> $C_f$ vs $Re_x$ (0 to $2.5\times10^6$). Curves: SA ($\nu' = 5$); SA–BC ($Tu = 0.9\%$, $\nu' = 0.02$);
> SA–BCM ($Tu = 0.9\%$, $\nu' = 0.02$); Experiments. The known poor T3A- behavior (late onset near
> $Re_x \approx 1\times10^6$ vs. experimental early rise) is visible.

> **Figure 5. Skin Friction coefficient distribution for the ERCOFTAC T3B flat plate test.**
> $C_f$ vs $Re_x$ (0 to $1\times10^6$). Curves: SA ($\nu' = 5$); SA–BC ($Tu = 6\%$, $\nu' = 0.02$);
> SA–BCM ($Tu = 6\%$, $\nu' = 0.02$); Experiments.

---

## IV. Results for low Reynolds number airfoils

### A. NACA-0021 airfoil  *(paper heading reads "NACA-0012" — see note)*

> **[Transcriber note]** The section heading and the Fig. 6 caption in the source read
> **"NACA-0012"**, but the body text, flow conditions and every other figure (Figs. 7–13) refer to the
> **NACA-0021** airfoil (the Menter et al. $\gamma$-SST validation case). Treat "NACA-0012" in the
> heading and Fig. 6 caption as a typo for **NACA-0021**.

The NACA-0021 airfoil at $Re = 2.7 \times 10^5$ case was used by Menter et al. to show the ability of
the $\gamma$-SST model in predicting the stall angle and the consequent steep Lift decrease with 2D
computations when laminar separations bubbles (LSBs) dominate the airfoil behavior at different AoAs
(Angles of Attack). The simulations have been performed in the range $0° \leq \text{AoA} \leq 28°$.
The mesh used for the computations is shown in Fig. 6. Lift coefficient values computed by the
transition models are in fair agreement with experimental data, and provide a good estimate of the
stall angle (Fig. 7). The SA-BC and SA-BCM models provide almost identical results. They perform
better than the $\gamma$-SST model in predicting AoA around 10° and also for the higher angles, where
the underlying SST model under-predicts Lift: it should be remembered that for $\gamma \to 1$ a
transition model reverts to the fully turbulent model it is based on, and therefore any deficiency of
the latter necessarily affects the results computed by means of the former. The better predictions of
the SA model give an advantage to the SA-based transition models, but it can be shown (Figs. 8, 9, 10,
11, 12, 13) that the SA-BCM and the $\gamma$-SST models do not significantly differ in how they predict
the separation and reattachment locations: the laminar separation bubbles have similar extensions.

Data reported in Swalwell are not very clear, thus it was not possible to accurately extract Drag
coefficient values and pressure coefficient profiles for extending the discussion.

> **Figure 6. O-grid mesh for the NACA-0021 airfoil** *(caption reads "NACA-0012")*. 721×242 points,
> 173520 cells. Boundaries are located at approximately 50 times the chord length from the airfoil.

> **Figure 7. Lift coefficient at different angles of attack for the NACA-0021 airfoil.** Experimental
> points are the same used in Menter et al. and taken from Swalwell. $C_L$ vs AoA [deg] (0 to 30).
> Curves: SA ($\nu' = 5$); SA–BC ($Tu = 0.6\%$, $\nu' = 0.02$); SA–BCM ($Tu = 0.6\%$, $\nu' = 0.02$);
> SST ($Tu = 0.924\%$, $\mu/\mu_t = 100$); $\gamma$-SST ($Tu = 0.924\%$, $\mu/\mu_t = 100$);
> Experiments. SA–BC and SA–BCM overlay and capture stall near AoA ≈ 15°.

> **Figure 8. Contours of non-dimensional x-velocity for the NACA-0021 airfoil at AoA = 0°.** Four
> panels: Upper left SA; Upper right SA-BCM; Lower left SST; Lower right $\gamma$-SST.
> ($U/U_\text{inf}$ colorbar −0.1 to 1.2.)

> **Figure 9. Contours of non-dimensional x-velocity for the NACA-0021 airfoil at AoA = 5°.** Four
> panels (SA / SA-BCM / SST / $\gamma$-SST). ($U/U_\text{inf}$ colorbar −0.1 to 1.4.)

> **Figure 10. Contours of non-dimensional x-velocity for the NACA-0021 airfoil at AoA = 10°.** Four
> panels (SA / SA-BCM / SST / $\gamma$-SST). ($U/U_\text{inf}$ colorbar −0.05 to 1.75.)

> **Figure 11. Contours of non-dimensional x-velocity for the NACA-0021 airfoil at AoA = 15°.** Four
> panels (SA / SA-BCM / SST / $\gamma$-SST). ($U/U_\text{inf}$ colorbar −0.15 to 1.80.)

> **Figure 12. Contours of non-dimensional x-velocity for the NACA-0021 airfoil at AoA = 20°.** Four
> panels (SA / SA-BCM / SST / $\gamma$-SST). ($U/U_\text{inf}$ colorbar −0.15 to 1.15; large
> suction-side separation.)

> **Figure 13. Contours of non-dimensional x-velocity for the NACA-0021 airfoil at AoA = 25°.** Four
> panels (SA / SA-BCM / SST / $\gamma$-SST). ($U/U_\text{inf}$ colorbar −0.15 to 1.15.)

### B. FFA-W3-301 airfoil

The FFA-W3-301 airfoil at $Re = 3 \times 10^6$ was used in Corson et al. to compare the behavior of
the $\gamma$ and $\gamma - Re_\theta$ models using both the SA and the SST models as a base, and
employing the implementations available in the commercial solver AcuSolve. A series of proprietary
airfoils for state-of-the-art horizontal axis wind turbines was also used for the same purpose,
demonstrating the reliability of such models. The computations were limited to AoAs before stall:
according to Corson et al., the inability of RANS simulations to predict stall is not related to the
ability in predicting transition phenomena. The statement seems to be contradicted by the results
previously shown for the NACA-0021 airfoil; it is, however, nothing but meaningless (see Sec. IV.C).
At $Re = 3 \times 10^6$, the flow over a FFA-W3-301 airfoil does not exhibit significant recirculation
regions, as it can be deduced from the Pressure coefficient profiles (Figs. 17, 18, 19, 20). The Lift
values are still better predicted by the transitional models (Fig. 15), but differences from the fully
turbulent results are less pronounced than in the previous case. What results to be really improved is
the Drag coefficient prediction, and consequently its ratio to the Lift coefficient (Fig. 16): the
fact demonstrates that transitional models do not only make sense when a transitional behavior has to
be reproduced, since their ability to damp the turbulence inside the boundary layer is sufficient to
make them performing significantly better than traditional RANS models in low Reynolds number flows.
In this case the pressure profile over the airfoil plays only a secondary role, while a correct
prediction of the viscous Drag allows for better estimates of the airfoil performance at angles prior
to the stall.

Comparison of the results here presented with the ones reported in Corson et al. is encouraged.
Firstly, it is interesting to note that the SA-BC and SA-BCM results almost match those obtained in
the cited reference with the $\gamma$-SA and $\gamma$-$Re_\theta$-SA models, which in turn only
marginally differ from the ones computed with the respective SST-based models. Drag is however better
predicted by the present computations for both fully turbulent and transitional models, whereas a
problematic $C_p$ profile is obtained in Corson et al. at AoA = −4.2° with the Spalart-Allmaras model.
It would be hard to guess where exactly the differences originate from, but it should be noted that
both the meshing strategy - structured O-grid against unstructured - and the numerical methods -
cell-centered FVM against FEM -, could be responsible for a significant impact on the solution.

> **Figure 14. O-grid mesh for the FFA-W3-301 airfoil.** 1079×290 points, 311542 cells. Boundaries at
> approximately 50 times the chord length from the airfoil.

> **Figure 15. Lift coefficient at different AoAs for the FFA-W3-301 airfoil.** $C_L$ vs AoA [deg]
> (−6 to 12). Curves: SA ($\nu' = 5$); SA–BC ($Tu = 0.1\%$, $\nu' = 0.02$); SA–BCM
> ($Tu = 0.1\%$, $\nu' = 0.02$); SST ($Tu = 0.157\%$, $\mu/\mu_t = 30$); $\gamma$-SST
> ($Tu = 0.157\%$, $\mu/\mu_t = 30$); Experiments (from Corson et al.). All curves nearly collapse
> onto the experimental line.

> **Figure 16. Ratio of Drag coefficient on Lift coefficient for the FFA-W3-301 airfoil.** $C_D$ vs
> $C_L$ (drag polar, −0.5 to 2). Same model set as Fig. 15. Transitional models (SA–BC/SA–BCM,
> $\gamma$-SST) track experiment; fully turbulent SA and SST over-predict $C_D$ substantially.

> **Figures 17–20. Pressure coefficient distribution for the FFA-W3-301 airfoil** at AoA = −4°
> (exp. −4.2°), 2° (2.0°), 6° (6.1°) and 12° (12.3°) respectively. $C_p$ vs $x/C$. Same model set.
> All models agree closely with experiment; no significant separation.

> **Figures 21–24. Contours of Eddy Viscosity Ratio (EVR) for the FFA-W3-301 airfoil** at AoA = −4°,
> 2°, 6° and 12° respectively. Four panels each (SA / SA-BCM / SST / $\gamma$-SST). Transitional models
> show delayed/reduced eddy-viscosity growth over the forward chord relative to fully turbulent
> baselines.

### C. NLF(1)-0416 airfoil

While NACA-0021 and FFA-W3-301 airfoils are used for blades in vertical and horizontal axis wind
turbines, respectively, the NLF(1)-0416 airfoil has been designed and tested at the NASA-Langley
Research Center for employment in the natural laminar flow design of wings. The case (at
$Re = 10^6$) is used here to demonstrate that an accurate stall prediction is not always achievable
with the LCTM transitional models. The flow over the airfoil is simulated at the lower Reynolds
number for which experimental data are available. The flow regime ensures the presence of laminar
separation bubbles, expected not to show at Reynolds numbers equal to or higher than $4 \times 10^6$
for this airfoil. Fig. 26 clearly reveals the inadequacy of the three transitional models to predict
stall conditions at the expected AoA, while Lift is slightly over-predicted for lower angles, as it
was found with the FFA-W3-301 airfoil. Differences between the SA-BC model - or the SA-BCM model,
giving the same results - and the $\gamma$-SST model probably depend only on the different baseline
RANS model. Again, the Drag coefficient values fairly resemble the experimental data before stall,
and the advantage over fully turbulent models is not to be questioned.

Differently from what was observed with the FFA-W3-301 airfoil, the $C_p$ profiles exhibits some
notable differences between SA-BCM and $\gamma$-SST models, and at AoA = 5° also between the SA-BCM and
the original SA-BC model. The discrepancies are evidently due to different LSB length estimates (cf.
Figs. 33 and 34). Since complete $C_p$ profiles are reported only for $Re = 4 \times 10^6$ in the
reference, it was not possible to determine which model provides the better predictions. The results
were however reported for completeness.

> **Figure 25. O-grid mesh for the NLF(1)-0416 airfoil.** 1081×272 points, 292680 cells. Boundaries at
> approximately 50 times the chord length from the airfoil.

> **Figure 26. Lift coefficient at different AoAs for the NLF(1)-0416 airfoil.** $C_L$ vs AoA [deg]
> (−5 to 20). Curves: SA ($\nu' = 5$); SA–BC ($Tu = 0.1\%$, $\nu' = 0.02$); SA–BCM
> ($Tu = 0.1\%$, $\nu' = 0.02$); SST ($Tu = 0.157\%$, $\mu/\mu_t = 10$); $\gamma$-SST
> ($Tu = 0.157\%$, $\mu/\mu_t = 10$); Experiments (from Somers). All models over-predict $C_L$ past
> AoA ≈ 12° (stall not captured).

> **Figure 27. Ratio of Drag coefficient on Lift coefficient for the NLF(1)-0416 airfoil.** $C_D$ vs
> $C_L$ (−0.5 to 2). Same model set. Transitional models match the experimental polar; fully turbulent
> models over-predict drag.

> **Figures 28–31. Pressure coefficient distribution for the NLF(1)-0416 airfoil** at AoA = −5°, 0°,
> 5° and 10° respectively. $C_p$ vs $x/C$. Same model set. Sharp $C_p$ kinks mark the LSB
> transition location; model differences appear in bubble length.

> **Figure 32. Contours of Eddy Viscosity Ratio for the NLF(1)-0416 airfoil at AoA = −5°.** Four panels
> (SA / SA-BCM / SST / $\gamma$-SST). EVR colorbar 10 to 110.

> **Figures 33–34. Contours of non-dimensional x-velocity for the NLF(1)-0416 airfoil** at AoA = 0° and
> 5° respectively. Four panels each (SA / SA-BCM / SST / $\gamma$-SST). Zoomed on the LSB region
> (x/C ≈ 0.4–0.7); reveal differing bubble lengths between models.

> **Figure 35. Contours of Eddy Viscosity Ratio for the NLF(1)-0416 airfoil at AoA = 10°.** Four panels
> (SA / SA-BCM / SST / $\gamma$-SST). EVR colorbar 25 to 300.

---

## V. Conclusion

The SA-BCM transitional model is successful in fixing the major deficiencies of the SA-BC formulation
without altering its results in low Reynolds numbers flows significantly. A robust implementation was
chosen so to make the model as independent from the specific solver as possible, but it is not
strictly required in order to exploit the advantages over the old formulation. Some basic cases were
employed for showing how the model can compete with a more complex - but also more general - one at a
reduced cost (one equation instead of three), what could be valuable in contexts where minimization of
the computational cost is mandatory, e.g., design optimization. The model still suffers from some
weaknesses - e.g., the inability to evolve the turbulence intensity levels along the domain -, but
these are not necessarily a negative factor for the external flows it is aimed to - specifying a
turbulence intensity once and for all avoids the need for the user to carefully tune the boundary
conditions as needed for the SST-based models. The most important lesson learned is that a complex
formulation is not always the best answer to practical needs, even in the context of transition
modeling. Without moving out from the LCTM framework, some categories of flows can be dealt with by
very simple models, given that they are implemented the right way so ensure mathematical and physical
consistency and numerical stability. The task is left to a companion paper to extend the validation of
the model with more challenging test cases, not necessarily involving external aerodynamics. Although
the model is not as general as other currently available alternatives, and not equally well-founded
through theoretical arguments, it could be still considered a valuable alternative to more demanding
formulations, or used as a base to develop a new generation of transitional models having the
computational efficiency as their main objective. Nevertheless, improvements and generalization will
hopefully be achieved in the near future so to make the model more competitive. Coupling of the BCM
transition model with other fully turbulent RANS closures would also be interesting to explore.

---

## Appendix — Complete formulation as implemented in ANSYS Fluent

> The complete formulation of the revised model as implemented in ANSYS Fluent, including the negative
> continuation of the baseline Spalart-Allmaras model, is reported below.
> **This block is the authoritative, copy-for-code version of the SA-BCM model.**

**Transport equation**

$$
\frac{\partial(\rho\tilde{\nu})}{\partial t} + \frac{\partial(\rho\tilde{\nu}U_j)}{\partial x_j}
= \rho\,\gamma\,P_{\tilde{\nu}} - \rho\,D_{\tilde{\nu}}
+ \frac{1}{\sigma}\left[\frac{\partial}{\partial x_j}\!\left(\Gamma_{\tilde{\nu}}\frac{\partial\tilde{\nu}}{\partial x_j}\right)
+ \rho\,C_{b2}\frac{\partial\tilde{\nu}}{\partial x_j}\frac{\partial\tilde{\nu}}{\partial x_j}\right]
$$

**Eddy viscosity**

$$
\mu_T =
\begin{cases}
\rho\tilde{\nu}f_{v1} & \tilde{\nu} > 0 \\
0 & \tilde{\nu} \leq 0
\end{cases}
\qquad
f_{v1} = \frac{\chi^3}{\chi^3 + C_{v1}^3}
\qquad
\chi = \frac{\tilde{\nu}}{\nu}
$$

**Production** (with negative continuation)

$$
P_{\tilde{\nu}} =
\begin{cases}
C_{b1}\tilde{S}\tilde{\nu} & \tilde{\nu} \geq 0 \\
C_{b1}S\,g_n\,\tilde{\nu} & \tilde{\nu} < 0
\end{cases}
\qquad
g_n = 1 - \frac{10^3\chi^2}{1 + \chi^2}
$$

**Modified vorticity**

$$
\overline{S} = \frac{\tilde{\nu}}{\kappa^2 d_w^2}\,f_{v2}
\qquad
f_{v2} = 1 - \frac{\chi}{1 + \chi f_{v1}}
$$

$$
\tilde{S} =
\begin{cases}
S + \overline{S} & \overline{S} \geq -C_{v2}S \\[6pt]
S + \dfrac{S\left(C_{v2}^2 S + C_{v3}\overline{S}\right)}{(C_{v3} - 2C_{v2})S - \overline{S}} & \overline{S} < -C_{v2}S
\end{cases}
$$

$$
S = \Omega = \sqrt{2\,\Omega_{ij}\Omega_{ij}}
\qquad
\Omega_{ij} = \frac{1}{2}\left(\frac{\partial U_j}{\partial x_i} - \frac{\partial U_i}{\partial x_j}\right)
$$

**Destruction** (with negative continuation)

$$
D_{\tilde{\nu}} =
\begin{cases}
C_{w1}f_w\left(\dfrac{\tilde{\nu}}{d_w}\right)^2 & \tilde{\nu} \geq 0 \\[8pt]
-C_{w1}\left(\dfrac{\tilde{\nu}}{d_w}\right)^2 & \tilde{\nu} < 0
\end{cases}
\qquad
f_w = g\left(\frac{1 + C_{w3}^6}{g^6 + C_{w3}^6}\right)^{1/6}
$$

$$
g = r + C_{w2}(r^6 - r)
\qquad
r = \frac{\tilde{\nu}}{\tilde{S}\kappa^2 d_w^2}
$$

**Diffusivity** (with negative continuation)

$$
\Gamma_{\tilde{\nu}} =
\begin{cases}
\mu + \rho\tilde{\nu} & \tilde{\nu} \geq 0 \\[6pt]
\mu + \rho\tilde{\nu} + \dfrac{1}{2}\dfrac{\rho^2\tilde{\nu}^2}{\mu} & \tilde{\nu} < 0
\end{cases}
$$

**Intermittency / transition (BCM)**

$$
\gamma = 1 - e^{-\left(\sqrt{Term1} + \sqrt{Term2}\right)}
$$

$$
Term1 = \frac{\max\!\left(\dfrac{Re_v}{2.193} - Re_{\theta c},\; 0\right)}{\chi_1 Re_{\theta c}}
$$

$$
Re_v = \frac{\rho\, d_w^2}{\mu}\,S
\qquad
Re_{\theta c} = 803.73\,(Tu_\infty + 0.6067)^{-1.027}
$$

$$
Term2 = \max\!\left(\frac{1}{\chi_2}\frac{\mu_T}{\mu},\; 0\right)
$$

**Constants**

| Constant | Value | | Constant | Value | | Constant | Value |
|---|---|---|---|---|---|---|---|
| $\sigma$ | $2/3$ | | $C_{b1}$ | $0.1355$ | | $C_{b2}$ | $0.622$ |
| $C_{v1}$ | $7.1$ | | $C_{v2}$ | $0.7$ | | $C_{v3}$ | $0.9$ |
| $C_{w1}$ | $\dfrac{C_{b1}}{\kappa^2} + \dfrac{1 + C_{b2}}{\sigma}$ | | $C_{w2}$ | $0.3$ | | $C_{w3}$ | $2.0$ |
| $\kappa$ | $0.41$ | | $\chi_1$ | $0.002$ | | $\chi_2$ | $0.02$ |

> $C_{w1} = \dfrac{C_{b1}}{\kappa^2} + \dfrac{1 + C_{b2}}{\sigma} = \dfrac{0.1355}{0.41^2} + \dfrac{1 + 0.622}{2/3} \approx 0.8062 + 2.433 = 3.239$

**Acknowledgments.** The authors gratefully thank Dr. Christopher L. Rumsey at the NASA Langley
Research Center for having encouraged their mutual cooperation.

---

## References

1. Spalart, P. R., and Allmaras, S. R., "A One-Equation Turbulence Model for Aerodynamic Flows," *La Recherche Aérospatiale*, Vol. 1, 1994, pp. 5–21.
2. Bas, O., Cakmakcioglu, S. C., and Kaynak, U., "A Novel Intermittency Distribution Based Transition Model For Low-Re Number Airfoils," *31st AIAA Applied Aerodynamics Conference*, AIAA, 2013. https://doi.org/10.2514/6.2013-2531.
3. Cakmakcioglu, S. C., Bas, O., and Kaynak, U., "A correlation-based algebraic transition model," *Proceedings of the Institution of Mechanical Engineers, Part C: Journal of Mechanical Engineering Science*, 2018. https://doi.org/10.1177/0954406217743537.
4. Rumsey, C. L., "The Spalart-Allmaras 1-equation BC Transitional Model," 2019. URL https://turbmodels.larc.nasa.gov/sa-bc_1eqn.html, last Updated: 03/17/2019.
5. Rumsey, C. L., Smith, B., and Huang, G., "Description of a website resource for turbulence modeling verification and validation," *40th Fluid Dynamics Conference and Exhibit*, 2010, p. 4742. https://doi.org/10.2514/6.2010-4742.
6. Palacios, F., Alonso, J. J., Duraisamy, K., Colonno, M., Hicken, J., Aranake, A., Campos, A., Copeland, S. S., Economon, T. D., Lonkar, A. K., and ant T. Taylor, T. W. L., "Stanford University Unstructured (SU2): An open-source integrated computational environment for multi-physics simulation and design," *51st AIAA Aerospace Sciences Meeting including the New Horizons Forum and Aerospace Exposition*, 2013. https://doi.org/10.2514/6.2013-287.
7. Palacios, F., Economon, T. D., Aranake, A., Copeland, S. R., Lonkar, A. K., Lukaczyk, T. W., Monosalvas, D. E., Naik, K. R., Padron, S., Tracey, B., A. Variyar, and Alonso, J. J., "Stanford University Unstructured (SU2): Analysis and Design Technology for Turbulent Flows," *52nd Aerospace Sciences Meeting*, 2014. https://doi.org/10.2514/6.2014-0243.
8. Menter, F. R., Esch, T., and Kubacki, S., "Transition modelling based on local variables," *Engineering Turbulence Modelling and Experiments 5*, Elsevier, 2002, pp. 555–564. https://doi.org/10.1016/B978-008044114-6/50053-3.
9. Langtry, R. B., "A correlation-based transition model using local variables for unstructured parallelized CFD codes," Ph.D. thesis, 2006.
10. Menter, F. R., Langtry, R. B., and Völker, S., "Transition Modelling for General Purpose CFD Codes," *Flow, Turbulence and Combustion*, 2006. https://doi.org/10.1007/s10494-006-9047-1.
11. Langtry, R. B., and Menter, F. R., "Correlation-Based Transition Modeling for Unstructured Parallelized Computational Fluid Dynamics Codes," *AIAA journal*, Vol. 47, No. 12, 2009, pp. 2894–2906. https://doi.org/10.2514/1.42362.
12. Menter, F. R., P., S. E., Liu, T., and Avancha, R., "A One-Equation Local Correlation-Based Transition Model," *Flow, Turbulence and Combustion*, 2015. https://doi.org/10.1007/s10494-015-9622-4.
13. Walters, D. K., and Cokljat, D., "A Three-Equation Eddy-Viscosity Model for Reynolds-Averaged Navier–Stokes Simulations of Transitional Flow," *Journal of fluids engineering*, Vol. 130, No. 12, 2008, p. 121401. https://doi.org/10.1115/1.2979230.
14. ANSYS, *Fluent Theory Guide*, 2016. Release 17.0.
15. Langtry, R. B., and Menter, F. R., "Transition Modeling for General CFD Applications in Aeronautics," *43rd AIAA Aerospace Sciences Meeting and Exhibit*, 2005. https://doi.org/10.2514/6.2005-522.
16. Menter, F. R., Langtry, R. B., Blair, M., Likki, S. R., Suzen, Y. B., Huang, P. G., and Völker, S., "A Correlation-Based Transition Model Using Local Variables Part I: Model Formulation," *Journal of Turbomachinery*, 2006. https://doi.org/10.1115/1.2184352.
17. Wilcox, D. C., *Turbulence modeling for CFD*, DCW Industries, Inc, 1993.
18. Allmaras, S. R., Johnson, F. T., and Spalart, P. R., "Modifications and clarifications for the implementation of the Spalart-Allmaras turbulence model," *7th international conference on computational fluid dynamics (ICCFD7)*, 2012, pp. 1–11.
19. Rumsey, C. L., "Apparent transition behavior of widely-used turbulence models," *International Journal of Heat and Fluid Flow*, Vol. 28, 2007. https://doi.org/10.1016/j.ijheatfluidflow.2007.04.003.
20. Rumsey, C. L., and Spalart, P. R., "Turbulence Model Behavior in Low Reynolds Number Regions of Aerodynamic Flowfields," *AIAA Journal*, 2009. https://doi.org/10.2514/1.39947.
21. Oliver, T. A., "A high-order, adaptive, discontinuous Galerkin finite element method for the Reynolds-Averaged Navier-Stokes equations," Ph.D. thesis, Massachusetts Institute of Technology, 2008.
22. Rumsey, C. L., "The Spalart-Allmaras Turbulence Model," 2018. URL https://turbmodels.larc.nasa.gov/spalart.html, last Updated: 09/10/2019.
23. ANSYS, *Fluent UDF Manual*, 2016. Release 17.0.
24. Kato, M., and Launder, B., "The modelling of turbulent flow around stationary and vibrating square cylinders," *9th Symp. on Turbulent Shear Flows*, Vol. 1, 1993, pp. 10–4.
25. Barth, T., and Jespersen, D., "The design and application of upwind schemes on unstructured meshes," *27th Aerospace sciences meeting*, 1989. https://doi.org/10.2514/6.1989-366.
26. Schubauer, G. B., and Klebanoff, P. S., "Contributions on the mechanics of boundary-layer transition," Tech. rep., NACA, 09 1955. 3498.
27. Swalwell, K. E., "The effect of turbulence on stall of horizontal axis wind turbines," Ph.D. thesis, Monash University, 2005.
28. Corson, D. A., Zamora, A., and Medida, S., "A Comparative Assessment of Correlation-based Transition Models for Wind Power Applications," *34th AIAA Applied Aerodynamics Conference*, 2016, p. 3129. https://doi.org/10.2514/6.2016-3129.
29. Somers, D. M., "Design and experimental results for a natural-laminar-flow airfoil for general aviation applications," 1981.

---

## Implementer's reconciliation notes (2020-2706 ↔ 2020-2714)

> **[Transcriber-authored — not part of either paper.]** These are the points where the two papers'
> notation or values differ in ways that will silently break a code port if missed. They are the
> reason to keep both files side by side.

**1. The `Term2` calibration constant — the single most dangerous trap.**
The final, implemented modified term is the same in both papers:

$$
Term2 = \max\!\left(50 \cdot \frac{\mu_T}{\mu},\; 0\right)
$$

- **2020-2706, Eq. 8:** written in *multiplier* form, $Term2 = \max(\chi_2\,\nu_t/\nu, 0)$ with
  $\chi_2 = 50$.
- **2020-2714, Eq. 11 + Appendix:** written in *reciprocal* form, $Term2 = \max\!\left(\tfrac{1}{\chi_2}\,\mu_T/\mu, 0\right)$
  with $\chi_2 = 0.02$ → $1/\chi_2 = 50$.
- **Do NOT** plug the $\chi_2 = 5.0$ from 2020-2714 Eqs. 7–8 into Eq. 11. That 5.0 belongs to the
  *original* SA-BC model (where $\chi_2$ is divided by $Re$). Using it in Eq. 11 gives
  $Term2 = 0.2\,\mu_T/\mu$ — wrong by a factor of **250**.

**2. Vorticity Reynolds number: local vs. maximum.**
- **2020-2706, Eq. 5:** $Re_\theta = Re_{v,\max}/2.193$ — uses the **maximum of $Re_v$ across the
  boundary layer** (a non-local search).
- **2020-2714, Eq. 3 + Appendix:** $Term1$ uses $Re_v/2.193$ — the **local** $Re_v$, no max.
- For a strictly local (LCTM-consistent) implementation, follow 2020-2714 and use the local $Re_v$
  pointwise. The published/official BC-BCM implementations (SU2, NASA TMR) compute $Re_v$ locally;
  the "max" in 2020-2706 Eq. 5 restates Menter's boundary-layer *observation*, not a per-cell
  operation. Verify against your reference implementation before shipping.

**3. Working-variable notation.**
- **2020-2706** uses $\nu_T$ for the SA *transport variable* (Eqs. 1–2) and $\nu_t$ for the actual
  *eddy viscosity* in Term2 (Eqs. 7–8). Same-looking symbols, different meaning.
- **2020-2714** uses the standard SA notation: $\tilde{\nu}$ for the transport variable,
  $\mu_T = \rho\tilde{\nu}f_{v1}$ (or $\nu_T = \mu_T/\rho$) for the eddy viscosity.
- In Term2, the ratio that matters is **(eddy viscosity)/(molecular viscosity)** =
  $\mu_T/\mu = \nu_T/\nu = \nu_t/\nu$ — i.e. $f_{v1}\chi$. It is **not** the raw transport variable
  ratio $\chi = \tilde{\nu}/\nu$ (they differ by the $f_{v1}$ factor, which matters in the viscous
  sublayer where $f_{v1}\to 0$).

**4. Baseline SA robustness — use 2020-2714's Appendix.**
2020-2706 references a positive SA baseline (no $f_{t2}$) but gives no details. 2020-2714 provides
the full **SA-neg (Oliver's negative continuation)**: modified vorticity limiter $\tilde{S}$
(Eq. 13), negative-branch production via $g_n$ (Eqs. 18–19), negative-branch destruction (Eq. 20),
and negative-branch diffusivity (Eq. 21). For a numerically robust SA-BCM, implement the Appendix,
not just Eqs. 12/15–17.

**5. Von Kármán constant.** 2020-2714 uses $\kappa = 0.41$ (matching the original SA references and
its own Appendix), explicitly overriding ANSYS Fluent's default $0.4187$. Keep $\kappa = 0.41$ so
that $C_{w1} = C_{b1}/\kappa^2 + (1+C_{b2})/\sigma$ evaluates consistently ($\approx 3.239$).

**6. `S` is vorticity, not strain.** Both papers define $S = \Omega = \sqrt{2\,\Omega_{ij}\Omega_{ij}}$
(mean vorticity magnitude), as in standard SA — explicitly *not* the strain-rate magnitude. The
$\Omega_{ij}$ sign convention differs between 2020-2714's body ($0.5(\partial_j U_i - \partial_i U_j)$)
and its Appendix ($\tfrac12(\partial_i U_j - \partial_j U_i)$), but the magnitude is identical, so
it does not affect the code.

**7. Constant set (shared).** $\chi_1 = 0.002$, $Re_{\theta c} = 803.73(Tu_\infty + 0.6067)^{-1.027}$,
and the divisor $2.193$ are identical in both papers. Standard SA constants: $\sigma = 2/3$,
$C_{b1} = 0.1355$, $C_{b2} = 0.622$, $C_{v1} = 7.1$, $C_{v2} = 0.7$, $C_{v3} = 0.9$, $C_{w2} = 0.3$,
$C_{w3} = 2.0$.

**8. Naming.** "BC" = Bas–Cakmakcioglu (original); "BCM" / "SA-BCM" = the revised model of these two
papers; the baseline is "SA" / "SA-neg". $\gamma$-SST and $\gamma$-$Re_\theta$-SST are Menter's
competitor models used only for comparison.

**9. Coupling.** Transition enters exactly one place: the SA **production** term is multiplied by
$\gamma$ (a.k.a. $\gamma_{BC}$). Destruction, diffusion and $C_{b2}$ cross-diffusion are untouched.
$\gamma \in (0,1]$; $\gamma \to 1$ recovers the fully turbulent SA model.
