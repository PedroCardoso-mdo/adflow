# Piotrowski & Zingg (2022, Special Session) — Numerical Behaviour of a Smooth Local Correlation-based Transition Model in a Newton-Krylov Flow Solver

> **Full transcript** — verbatim text extracted from the source PDF with
> `pdftotext -layout`, page-delimited, nothing removed except bare page-number
> lines. Layout is preserved, so tables and multi-column material read as they
> appear in the PDF. **Equations are raw PDF extractions**: superscripts,
> subscripts and math symbols may be flattened or split across lines — for any
> equation you intend to implement, confirm against the PDF itself (path below)
> or the digest companion file.
>
> Source PDF: `/home/mdo/Desktop/ARTIGOS SMOTH/PiotrowskiandZingg2022SpecialSession (3).pdf`

---

<!-- ===== PDF page 1 ===== -->

    Numerical Behaviour of a Smooth Local Correlation-based
       Transition Model in a Newton-Krylov Flow Solver

                                                Michael G. H. Piotrowski∗
                                   NASA Ames Research Center, Moffett Field, CA 94035
                           Institute for Aerospace Studies, University of Toronto, ON M3H 5T6

                                                     David W. Zingg†
                           Institute for Aerospace Studies, University of Toronto, ON M3H 5T6

             The numerical behaviour of transport-equation-based transition models, including both
         iterative and grid convergence, is influenced by the source terms. Transition models contain
         source terms that are large and highly nonlinear, and can be destabilizing in a strong implicit
         solver. Linearization strategies with varying levels of coupling are evaluated in conjunction
         with a source-term time step restriction to determine best-practices for solving the SA-sLM2015
         smooth local correlation-based transition model in an implicit Newton-Krylov flow solver.
         Achieving deep iterative convergence facilitates a detailed investigation of the grid convergence
         of these free-transition simulations, which are evaluated relative to fully-turbulent simulations
         performed using the Spalart-Allmaras turbulence model. Simulations of the NLF0416 general
         aviation airfoil, VA-2 supercritical airfoil, and NASA CRM-NLF wing-body geometry are
         performed over a range of grid levels. The results demonstrate that both a fully-coupled
         linearization strategy and a source-term time step restriction improve nonlinear convergence
         as the complexity of the free-transition simulations increases. In general, additional grid
         resolution is required for free-transition simulations relative to fully-turbulent simulations
         in order to achieve a similar level of accuracy, with the grid convergence of free-transition
         simulations sensitive to the streamwise grid spacings in the transition regions.

                                                      I. Introduction
    Local transport-equation-based transition models contain large and highly nonlinear source terms [1–4]. These
source terms introduce numerical stiffness through time scales that can vary significantly from those introduced by
the convection and diffusion terms. Unless suitable steps are taken to modify the basic algorithm, the implementation
of transition model equations in a strong implicit solver can lead to numerical instability. The development of
best-practices for solving these transport-equation-based transition models in implicit solvers can facilitate the deep
iterative convergence necessary for investigating their grid convergence [5]. Eça et al. [6] demonstrated that iterative
error must be reduced two to three orders of magnitude below the discretization error to prevent contamination, while
the elimination of numerical error, including both iterative and discretization error, is necessary in order to properly
evaluate the modelling error of these RANS-based transition models [7]. Recent progress has been made to evaluate
the grid convergence of transport-equation-based transition models [8], including the grid convergence of the original
Langtry-Menter transition model [1]. However, the smooth-variant coupled with the Spalart-Allmaras (SA) turbulence
model [9], SA-sLM2015 [10, 11], has not been thoroughly evaluated, where stiff and non-smooth functions in the
original transition model source terms could impact grid convergence.
    The question of how to best solve transport equations with strong source terms was emphasized with the introduction
of two-equation turbulence models, such as those based on the popular k- and k-ω models. Moryossef and Levy [12]
provide a review of solution strategies for the stable treatment of source terms in these models. The methods focus on
positivity preservation of the independent working variables though the implicit treatment of sinks and the explicit
treatment of sources. However, clipping is still introduced in order to prevent excessive solution growth. Work by
Lian et al. [13] demonstrated the effect of explicit and implicit treatments of source terms through stability-analysis
    ∗ NASA Pathways Ph.D. Candidate, University of Toronto, NASA Computational Aeroscience Branch, AIAA Student Member,

m.piotrowski@mail.utoronto.ca
    † University of Toronto Distinguished Professor of Computational Aerodynamics and Sustainable Aviation, Associate Fellow AIAA,

dwz@oddjob.utias.utoronto.ca

<!-- ===== PDF page 2 ===== -->

of a linear convection equation with a source term. The results demonstrate that the stability of the system can be
represented using the non-dimensional source-term time step, where a large positive source-term time step can lead to
strong solution growth and subsequent flow solver divergence. While implicit methods provide more stable updates than
explicit methods for the treatment of sinks, i.e. negative sources, the implicit treatment of sources can produce unstable
solution updates for large source-term time steps [13]. Specifically, implicit Euler time-marching with a source-term
time step larger than unity can produce solution updates with an unphysical sign change. This behaviour can help to
explain the common problem of positivity preservation of turbulence and transition model working variables [14].
    Strategies for the stable treatment of turbulence model source terms have recently been applied to transition model
equations as well. Similar to the work on positivity preservation of two-equation turbulence models by Ilinca and
Pelletier [15], Coder [16] reformulated the intermittency equation in the amplification factor transport (AFT) transition
model in order to solve for the logarithmic variable, γ̃ = log(γ). Although this strategy ensures a positive solution, it
introduces difficulties with the treatment of homogenous boundary conditions at the wall, i.e. γ̃ = 0, and can require
additional clipping to prevent excessive variable growth [12]. Mosahebi and Laurendeau [17, 18] investigated different
coupling and time-marching strategies for the Langtry-Menter LM2009 [1] transition model and the shear-stress
transport (SST) turbulence model [19]. Specifically, they investigated a fully-segregated treatment of the transition and
turbulence models, including a modified segregated treatment where the coupling between the two models is under
relaxed, and coupled treatments where the two transition and turbulence models are each solved independently, resulting
in two 2x2 Jacobians, and fully-coupled, resulting in a 4x4 Jacobian. The results demonstrate that of the four strategies
considered, the under-relaxed segregated approach provided the most robust convergence for the subsonic airfoil test
cases investigated [17, 18].
    Newton-Krylov solvers provide the potential for greatly accelerated nonlinear convergence [20–23]. The rate of
convergence of Newton’s method is closely linked to the the accuracy of the Jacobian. Neglecting terms or otherwise
manipulating the Jacobian can prevent quadratic convergence [22]. However, Blanco and Zingg [21] demonstrated the
potential benefits of using a loosely-coupled strategy for solving the SA turbulence model equations in an unstructured
Newton-Krylov flow solver. Specifically, they demonstrated that the loosely-coupled approach produced similar or
better nonlinear convergence relative to fully-coupled Newton-Krylov solvers, with the loosely-coupled approach also
benefiting from a 30% reduction in memory requirements. More recently, Yildrim et al. [24] investigated coupling
strategies for the SA turbulence model and mean-flow equations in an approximate Newton-Krylov flow solver. They
demonstrated that solving the SA turbulence model loosely-coupled to the mean-flow equations in the globalization
phase could improve robustness and efficiency, with the turbulence model solved fully-coupled in the inexact-Newton
phase to facilitate nonlinear convergence.
    As previously discussed, developing best-practices for achieving deep iterative convergence of these transport-
equation-based transition models facilitates a discussion of their grid convergence. Lopes et al. [8] investigated the grid
convergence of three popular transport-equation-based transition models [1, 3, 4] coupled with the SST turbulence
model [19]. They performed a thorough analysis of the effects of discretization strategy for the turbulence and transition
model convective terms on solution accuracy. The results demonstrate that the discretization strategy for the transition
model equations has a small influence on the solution accuracy, but the SST turbulence model discretization strategy
has a strong impact. However, this may be due to a difference in the freestream decay of the turbulent kinetic energy and
therefore the turbulence intensity, which has a significant impact on the transition location [25].
    The goal of the current work is to investigate the iterative and grid convergence of the SA-sLM2015 smooth
local correlation-based transition model [10, 11] in an implicit Newton-Krylov flow solver. To achieve this goal, the
numerical behaviour of free-transition simulations performed using two- and three-dimensional subsonic and transonic
transition test cases are evaluated relative to fully-turbulent simulations performed using the SA turbulence model. For
the free-transition simulations, fully-, loosely-coupled, and decoupled linearization strategies for the turbulence and
transition model equations are investigated along with a source-term time step restriction to determine best-practices for
solving transport equations with large and highly non-linear source terms in a strong implicit solver. Grid-refinement
studies are performed to investigate the grid convergence of the smooth transition model relative to fully-turbulent
simulations. Although the focus is on the SA-sLM2015 transition model, we expect some of our conclusions to be
applicable to other transport-equation-based transition models as well.
    The remainder of the paper is organized as follows: Section II provides an overview of the Newton-Krylov solution
algorithm, linearization strategies, and source-term time step restriction investigated in the current work, as well as
the SA-sLM2015 smooth local correlation-based transition model, Section III presents the numerical behaviour of the
fully-turbulent and free-transition simulations using a series of two- and three-dimensional transition test cases, with
Section IV discussing conclusions and future work.

<!-- ===== PDF page 3 ===== -->

                                                  II. Methodology

A. Flow Solution Algorithm
    The flow solution algorithm used in the current work is a three-dimensional structured multi-block finite-difference
solver [23, 26]. The governing equations are spatially discretized using summation-by-parts operators, with simultaneous
approximation terms applied to enforce boundary conditions and inter-block coupling. The mean-flow equations are
discretized using second-order summation-by-parts operators with artificial dissipation provided using the matrix-based
dissipation model developed by Swanson and Turkel [27]. First-order upwinding is applied to the convection terms in the
turbulence and transition model equations. The computational domain is decomposed into multiple blocks, resulting in
multi-block structured grids, which allow for efficient parallel computations. The flow solver was extensively validated
with fully-turbulent flow using the SA turbulence model as part of the Fifth AIAA Drag Prediction Workshop [28].
    A Newton-Krylov-Schur solution algorithm making use of a pseudo-transient continuation strategy is applied to
the set of discretized equations to drive the residual to a converged steady-state solution. Globalization is achieved
using an approximate-Newton phase, where an implicit Euler time marching strategy with local time linearization is
applied using an approximate, analytically derived Jacobian. After the residual drops several orders of magnitude,
an inexact-Newton method using a full second-order Jacobian is used to converge the system to a residual norm of
machine zero. The Jacobian can be formed explicitly or a Jacobian-free approach can be used. The large system of
linear equations generated at each iteration in both phases is solved using the GMRES linear solver [29], preconditioned
using an approximate-Schur parallel preconditioner [30].

B. SA-sLM2015 Transition Model
    In previous work, the authors developed a smooth, three-equation transition model, based on the LM2009 [1]
empirical correlations for bypass transition and transition due to Tollmien-Schlichting instabilities as well as the LM2015
[31] helicity-based empirical correlations for stationary crossflow instabilities. A novel coupling was introduced with the
one-equation SA turbulence model [9] that maintains the fully-local approach spearheaded by Langtry and Menter [1].
Smooth approximations were developed to replace stiff and non-differentiable functions in the LM2009 and LM2015
source terms to improve the numerical behaviour of the model, while maintaining the predictive capability of the original
correlations. Details of the transition model, designated SA-sLM2015, including validation with subsonic two- and
three-dimensional transition test cases, can be found in [10]. Recent work has focused on developing compressibility
corrections to extend the SA-sLM2015 transition model to transonic flow conditions, with the new extension designated
SA-sLM2015cc [11].

C. Linearization Strategies
    Four linearization strategies for the turbulence and transition model equations with varying levels of coupling
are considered: a fully-coupled approach, a loosely-coupled approach with the turbulence and transition model
equations coupled, a loosely-coupled approach with the transition model equations coupled with an independent
linearization for the turbulence model, and a decoupled approach with an independent linearization for each of the
turbulence and transition model equations. The four linearization strategies are investigated each with and without the
source-term time step restriction, presented below in Section II.D, with the coupling of the source-term Jacobian used to
evaluate the source-term eigenvalues mirroring the solver linearization strategy. Other than the source-term time step
restriction, no special considerations are applied to the turbulence and transition model production and destruction
terms. Loosely-coupled and decoupled linearization strategies for the turbulence and transition model equations are
implemented by neglecting off-diagonal components of the fully-coupled 8x8 Jacobian. The block-diagonal Jacobians
for the four linearization strategies are presented in Figure 1.
    The source terms in the turbulence and transition models are most destabilizing during the early stages of convergence.
To recover the nonlinear convergence properties of Newton’s method, the loosely-coupled and decoupled linearization
strategies are only adopted in the approximate-Newton phase of the solver, with the mean-flow, turbulence, and transition
model equations fully coupled in the inexact-Newton phase. The emphasis of this study is on the impact of decoupling
the turbulence and transition model equations on nonlinear convergence relative to a fully-coupled linearization strategy,
in particular the impact on robustness. Although using a sequential segregated algorithm can reduce wall-clock times for
these loosely-coupled and decoupled approaches in the approximate-Newton phase, the current implementation isolates
the effects of varying coupling strategies and avoids differences introduced by using a different solution algorithm, such
as the introduction of different preconditioners. A physicality check, solution update damping algorithm, and unsteady

<!-- ===== PDF page 4 ===== -->

        (a) fully coupled             (b) turb-trans coupled           (c) trans coupled        (d) decoupled

Fig. 1 Block-diagonal Jacobians produced by the four linearization strategies for the turbulence and transition
model equations with varying levels of coupling, with the grey elements representing filled Jacobian entries.

                                             analytical
                                      6      implicit
                                             explicit

                              u/u o

                                      -2

                                      -4

                                       -3       -2        -1       0        1        2     3
                                                               b t
Fig. 2 Normalized solution update over a range of time steps produced by implicit and explicit Euler time
marching applied to a linear convection equation with a scalar source term, b. Reproduced from [13].

residual line-search algorithm with backtracking are used to maintain stable solution updates, and an equation and
variable scaling approach is used to improve linear solver performance, with these modifications to the Newton-Krylov
algorithm presented in [10].

D. Source-Term Time Step Restriction
    Lian et al. [13] demonstrated the effects of treating source terms with implicit and explicit Euler time marching by
investigating a linear convection equation with a source term. The effect of the pseudo time-marching strategy on the
normalized solution update for this scalar equation is illustrated in Figure 2. An implicit Euler time marching strategy
applied to a source term produces large solution updates, which asymptote at a source-term time step of unity, after
which it produces solution updates with the opposite sign. These results can help to explain the common problem of
positivity preservation for turbulence and transition model dependent variables [14, 32]. Maintaining a source-term
time step less than unity allows for rapid solution growth while preventing unbound and unstable solution updates.
The source-term time step is determined by the product of the largest positive eigenvalue of the source-term Jacobian
matrix and the local time step. A source-term time step restriction of 0.8 was used for the results presented in the
current work. Because the turbulence and transition model source terms can remain active in a fully-converged solution,
the source-term time step restriction is disabled in the inexact-Newton phase after five successive iterations for which
neither the physicality check nor the unsteady residual line-search algorithm is active, and is re-introduced if these
measures become active. More details of the source-term time step restriction can be found in [10].

<!-- ===== PDF page 5 ===== -->

E. Estimating Grid Convergence
    The grid convergence of the fully-turbulent and free-transition simulations is evaluated using Richardson extrapolation
[5]. There are several important requirements for using Richardson extrapolation to report grid-converged values; for
example, the grids must be sufficiently fine to be in the asymptotic range of convergence in order to have full confidence
in the results [5]. Previous work has demonstrated the difficulty in achieving asymptotic convergence for simple test
cases with fully-turbulent RANS solutions [33–35]. The presence of discontinuities, such as shocks, and singularities
can complicate grid-convergence studies [5]. In the current work, Richardson extrapolation is used to evaluate the
relative grid convergence of the fully-turbulent and free-transition simulations. The grid-refinement results are evaluated
using the apparent order of convergence, p, the actual fractional error (AFE), and grid-convergence index (GCI), each
calculated for the finest grid level. These grid-convergence metrics are calculated according to ASME standards [36]
with some minor modifications∗ as follows,
                                                   ln |( f 3 − f 2 )/( f 2 − f 1 )| + q(p)
                                                p=                                          ,                                 (1)
                                                                    ln(r 21 )
                                                    rp − s                         f −f 
                                                                                         3      2
                                          q(p) = ln 21 p         ; s = 1 · sign                   ,                           (2)
                                                     r 32 − s                          f 2 −  f 1
                                                                             p
                                      f exact − f 1                      r 21 f 1 − f 2               f1 − f2
                                AFE =               ;        f exact =           p        ;   EFE =           ,               (3)
                                          f exact                            r 21 − 1                    f1

where f 1 , f 3 , and f exact represent the values on the finest and coarsest of the three grid levels being investigated, and the
grid-converged value estimated using Richardson extrapolation, respectively, r is the refinement factor, and EFE is the
estimated fractional error.
     For the results presented, the total number of grid points is approximately doubled between refinement levels,
with the refinement factors calculated using the total number of grid nodes at each grid level: r 21 = (N1 /N2 ) 1/2 and
r 21 = (N1 /N2 ) 1/3 for two- and three-dimensional cases, respectively. Using a consistent family of grids is important for
evaluating grid convergence. This can be done by removing every second node in each coordinate direction repeatedly
from a fine mesh. However, it can be preferable to use a refinement factor less than two. Grid sequences for the
NLF0416 airfoil and CRM-NLF geometry are produced using an automatic grid refinement algorithm that uses an
analytical representation of the grid based on B-spline volumes [37]. The VA-2 airfoil grid family is obtained by
generating a sequence of meshes with suitably varying grid spacings and growth ratios. Eça and Hoekstra introduced
modifications to the GCI calculation to prevent under- and over-estimating discretization error due to scatter in the data
and the possibility of the data being outside the asymptotic range [6]. The following modification is used in the current
work according to the Turbulence Modelling Resources website guidelines∗ ;
     • if p > 3.05,
                                                                1.25 · EFE 1.25∆ 
                                                                                    M
                                                GCI = max            3 −1
                                                                           ,          ,                                       (4)
                                                                   r 21       | f1|
                                            ∆ M = max(| f 2 − f 1 |, | f 3 − f 2 |, | f 3 − f 1 |),                           (5)

    • else,
                                                                    1.25 · EFE
                                                          GCI =         p      .                                              (6)
                                                                      r 21 − 1

                                                              III. Results
    Fully-turbulent and free-transition grid-refinement studies are presented for the NLF0416 general aviation airfoil,
VA-2 supercritical airfoil, and NASA CRM-NLF wing-body geometry. The grid-refinement studies are performed using
the fully-coupled linearization strategy presented in Section II.C and the source-term time step restriction presented
in Section II.D. Linearization strategies are investigated using free-transition simulations both with and without the
source-term time step restriction to determine best-practices for solving the transition model in an implicit solver.
Convergence histories are illustrated by the L2-norm of the total residual, which includes the RANS, turbulence, and
    ∗ https://turbmodels.larc.nasa.gov/uncertainty_summary.pdf, accessed November 2021

<!-- ===== PDF page 6 ===== -->

transition model equations, including the row (equation) scaling measure presented in [10], normalized by the total
residual at the first iteration,
                                                               ||R(n) ||2
                                                     Rd(n) =              ,                                                (7)
                                                               ||R(0) ||2
versus the computational cost, which is measured in units of equivalent residual evaluations. This is determined as the
wall-clock time normalized by the average cost of computing the total residual for each simulation. It is important to
note that the addition of the two transition model transport equations approximately doubles the computational cost of
each nonlinear iteration, which is not illustrated when using units of equivalent residual evaluations. In addition, the
increased computational cost associated with increasing the number of grid points is not explicitly illustrated. Instead,
plotting against equivalent residual evaluations illustrates the relative difficulty of solving the nonlinear system. The
grid-refinement results are plotted versus the grid spacing squared, h2 , while the grid-convergence metrics are calculated
according to Equations 1 - 6.
    Limiters are provided in the matrix-based dissipation model to prevent zero artificial viscosity at stagnation points
and sonic lines, where the linear and nonlinear eigenvalue scaling factors approach zero, respectively [27]. The linear
eigenvalue limiter parameter Vl , is set to zero, as it was found to be overly dissipative in the laminar boundary layer, while
the nonlinear eigenvalue limiter parameter, Vn , is set to 0.25, which is between the upper and lower bound of 0.20 and
0.30 recommended by Swanson and Turkel [27]. For the subsonic NLF0416 airfoil transition test case, a fourth-difference
dissipation coefficient of 0.04 is used, while for the VA-2 supercritical airfoil second- and fourth-difference coefficients
of 2.00 and 0.04 are used, respectively. Due to the increased complexity of the NASA CRM-NLF wing-body geometry,
the dissipation coefficients are increased to 3.00 and 0.06 for the second- and fourth-difference coefficients, respectively,
and to 0.30 for the nonlinear eigenvalue limiter, Vn .
    Residual drop tolerances, Rd(n) (Equation 7), of 10−4 and 10−5 are used for the switch from the approximate- to the
inexact-Newton phase for the two- and three-dimensional simulations, respectively. Relative and absolute convergence
tolerances for the total residual are specified as 10−15 and 10−10 , respectively. However, for the majority of the results
presented the absolute residual tolerance is met first. The same solver settings are used for the fully-turbulent and
free-transition simulations. Free-transition simulations are performed using the SA-sLM2015 transition model [10],
with compressibility corrections for Tollmien-Schlichting and stationary crossflow instabilities used for the transonic
VA-2 and CRM-NLF transition test cases [11]. The SA-noft2-neg variant of the SA turbulence model [9, 14] is used for
both the fully-turbulent and free-transition simulations, and the QCR2000 correction [38] is used for simulations of the
NASA CRM-NLF.

A. NLF0416 General Aviation Airfoil
    The NLF0416 airfoil is a natural-laminar-flow general aviation airfoil designed by Somers [39]. Experimental
results of the NLF0416 airfoil were obtained in the low turbulence pressure tunnel (LTPT) at NASA Langley Research
Center [39]. Transition locations were measured for several angles of attack at a Mach number of 0.10 and a Reynolds
number of 4.0 × 106 . The transition mechanisms for these conditions are natural and separation-induced transition.
The turbulence intensity in the wind tunnel used for the experimental study was not specified; however a value of
Tu = 0.15% is assumed in order to be consistent with results obtained using the Transition Modelling and Predictive
Capabilities Seminar (TCMPS) guidelines [40]. The SA-sLM2015 transition model was validated using this test case in
previous work by the current authors [10]. Structured C-topology multi-block grids with characteristics similar to the
NLF0416 grids provided by the TMPCS are used in the current work [40], with leading- and trailing-edge aspect ratios
of approximately 200 used for the initial grid-refinement study. The dimensions of the structured C-grids are presented
in Table 1.
    Grid-refinement studies are presented at zero angle of attack using the grids in Table 1 for fully-turbulent flow and
with free transition. Low-Mach-number preconditioning is not used for the results presented; however, simulations at a
higher Mach number of 0.30 were performed with the results demonstrating similar iterative and grid convergence.
The residual convergence histories for the fully-turbulent and free-transition grid-refinement studies are presented in
Figure 3. The total residual converges approximately 12-13 orders of magnitude for each case, with the solver exiting
after hitting the absolute residual tolerance of 10−10 . Although the simulations performed on the L5 grid level satisfy
the same convergence tolerance, the results are not included in Figure 3. The additional stiffness from the small grid
spacings significantly affects the linear solver performance, requiring larger Krylov subspace sizes and significantly
more equivalent residual evaluations to converge. The free-transition simulations converge with approximately 1.2 to
1.6 times the number of equivalent residual evaluations required for the fully-turbulent simulations.

<!-- ===== PDF page 7 ===== -->

                                         Table 1      NLF0416 structured multi-block C-grid dimensions.

                                         grid level     chord × off-wall nodes                      ∆s (chord)    avg/max y+
                                            L0                443 × 79                              5.00 × 10−6    0.61 / 1.65
                                            L1                613 × 109                             3.54 × 10−6    0.42 / 1.03
                                            L2                885 × 157                             2.50 × 10−6    0.30 / 0.63
                                            L3               1225 × 217                             1.77 × 10−6    0.21 / 0.42
                                            L4               1769 × 313                             1.25 × 10−6    0.15 / 0.28
                                            L5               2449 × 433                             0.89 × 10−6    0.13 / 0.21

                        103                                                                          103
                                                                     L0                                                                    L0
                        101                                          L1                              101                                   L1
Relative Residual, d

                                                                            Relative Residual, d
                        10−1                                         L2                             10−1                                   L2
                                                                     L3                                                                    L3
                        10−3                                         L4                             10−3                                   L4
                        10−5                                                                        10−5
                        10−7                                                                        10−7
                        10−9                                                                        10−9
                   10−11                                                                       10−11
                   10−130.0        0.5        1.0          1.5        2.0                      10−130.0           0.5        1.0    1.5    2.0
                               Equivalent Residual Evaluations 1e5                                          Equivalent Residual Evaluations 1e5
                                         (a) SA                                                                    (b) SA-sLM2015

  Fig. 3 Residual convergence histories for the NLF0416 airfoil simulated fully turbulent and with free transition
  at Ma = 0.10, Re = 4.0 × 106 , Tu = 0.15%, and α = 0◦ using the grids presented in Table 1.

       A linearization study for the free-transition simulations on the L4 grid level is presented, with the residual convergence
  histories illustrated in Figure 4. With the source-term time step restriction, the fully-coupled solution strategy requires
  more equivalent residual evaluations to converge relative to the loosely-coupled and decoupled approaches. Although
  the increased coupling can provide a more accurate solution update, the linearization of a large source term can require
  a smaller local time step to prevent unstable solution updates, as demonstrated in Figure 2. This behaviour is illustrated
  by comparing Figures 4a and 4b, where the fully-coupled solution without the source-term time step restriction fails
  to converge. Counter to the behaviour of the fully-coupled simulations, the loosely-coupled and decoupled strategies
  require more equivalent residual evaluations to converge with the source-term time step restriction. This behaviour can
  likely be attributed to the relative simplicity of the problem. As the flow complexity increases, the more accurate update
  provided by the fully-coupled approach improves convergence, as demonstrated in Sections III.B and III.C.
       Diskin et al. [35] demonstrate that the grid aspect ratio at geometric singularities, such as a sharp trailing edge,
  is the primary factor affecting accuracy and grid convergence of fully-turbulent simulations. To investigate the grid
  convergence of free-transition simulations, grid-refinement studies were performed with trailing-edge grid aspect ratios
  of 200, 20, and 2. The grids in the current work are generated with streamwise cluster locations only specified at the
  leading and trailing edges. Therefore, as the trailing-edge aspect ratio decreases, the streamwise spacings on the upper
  and lower surfaces of the airfoil increase away from the trailing edge. While this does not have a noticeable impact on
  fully-turbulent results, the transition onset locations move downstream as the streamwise grid spacing decreases in the
  vicinity of the transition points with increasing trailing-edge aspect ratio. These small changes in the transition onset
  location can significantly affect both lift and drag. A trailing-edge aspect ratio of approximately 20 is used for the grid
  convergence results presented, as it provides sufficient streamwise grid resolution at the transition locations and at the
  trailing edge. The leading-edge aspect ratio is maintained at approximately 200. The fully-turbulent and free-transition
  grid-convergence results produced at zero angle of attack are presented in Figure 5.

<!-- ===== PDF page 8 ===== -->

                        103                                                                                103
                                                          fully coupled                                                                       fully coupled
                        101                               turb-trans coupled                               101                                turb-trans coupled
Relative Residual, d

                                                                                   Relative Residual, d
                        10−1                              trans coupled                                    10−1                               trans coupled
                                                          decoupled                                                                           decoupled
                        10−3                                                                               10−3
                        10−5                                                                               10−5
                        10−7                                                                               10−7
                        10−9                                                                               10−9
                   10−11                                                                              10−11
                   10−130.0              0.5           1.0        1.5        2.0                      10−130.0               0.5            1.0       1.5         2.0
                                 Equivalent Residual Evaluations 1e5                                                 Equivalent Residual Evaluations 1e5
                                                (a) w/o STTS                                                                       (b) w/ STTS

 Fig. 4 Linearization study residual convergence histories for the NLF0416 airfoil simulated on the L4 grid
 (Table 1) using the SA-sLM2015 transition model at Ma = 0.10, Re = 4.0 × 106 , Tu = 0.15%, and α = 0◦ .

          .0097                                                    SA                                                                                 SA
                                pL3-L5 = 1.07                                                                      pL3-L5 = 1.01
                                                                   SA-sLM2015            0.488                                                        SA-sLM2015

          .0096
                                                                                         0.487
                              AFE = 0.37%                                                                        AFE = 0.09%
                              GCI = 0.47%                                                                        GCI = 0.11%
          .0095
        Cd

                         //                                                                                 //
                                                                                       Cl

          .0053                                                                          0.452
                              AFE = 0.20%
                              GCI = 0.25%                                                                          pL3-L5 = 1.01

          .0052                                                                          0.451

                                pL3-L5 = 1.20                                                                    AFE = 0.03%
                                                                                                                 GCI = 0.04%
          .0051                                                                          0.450
               0                         5E-06      N-1   1E-05         1.5E-05               0                             5E-06     N-1     1E-05         1.5E-05

  Fig. 5 Grid-convergence results for the NLF0416 airfoil simulated fully turbulent and with free transition at
  Ma = 0.10, Re = 4.0 × 106 , Tu = 0.15%, and α = 0◦ .

      The fully-turbulent and free-transition lift and drag coefficient grid-convergence metrics for the finest grid level are
  calculated using the methods described in Section II.E, with the results presented in Figure 5. The free-transition grid
  convergence for lift is similar to the fully-turbulent results, but with higher AFE and GCI values on the finest grid level.
  The drag coefficient decreases with increasing grid refinement for the free-transition simulations, and increases with
  refinement for the simulations with fully-turbulent flow, with the free-transition results producing smaller AFE and GCI
  values. Figure 6 presents the fully-turbulent and free-transition pressure and skin friction coefficient profiles produced
  on the L1 to L5 grid levels. The free-transition simulations produce pressure coefficient profiles that agree better with
  the pressure profile from experiment [39]. As previously discussed, for the free-transition simulations the transition
  onset locations move downstream with increasing grid refinement. The difference in the transition onset locations
  between the L3 and L5 grids is less than 2% chord and is approximately 1% chord between the L4 and L5 grids.
      The NLF0416 sharp trailing-edge transition test case highlights the impact of streamwise grid spacing at the transition
  locations in addition to the trailing-edge geometric singularity on the grid convergence of free-transition simulations.

<!-- ===== PDF page 9 ===== -->

                                            LTPT Experiment                                                                LTPT Experiment
                                            SAL1                                                                           SA-sLM2015L1
                                            SAL2                  0.015                                                    SA-sLM2015L2          0.015
     -1                                     SAL3                                   -1                                      SA-sLM2015L3
                                            SAL4                                                                           SA-sLM2015L4
                                            SAL5                                                                           SA-sLM2015L5

   -0.5                                                                          -0.5
                                                                  0.01                                                                           0.01

       0                                                                              0

                                                                      Cf,x

                                                                                                                                                     Cf,x
  Cp

                                                                                 Cp
                     upper

                                                                  0.005                                                                          0.005
    0.5                                                                           0.5
                                                                                                    upper

       1                                                                              1
                                                                                                                          lower
                                            lower
                                                                  0                                                                              0
           0   0.2       0.4          0.6       0.8           1                           0   0.2       0.4         0.6           0.8        1
                               x/c                                                                            x/c
                             (a) SA                                                                   (b) SA-sLM2015

Fig. 6 Pressure and skin friction coefficient profiles for the NLF0416 airfoil simulated fully turbulent and with
free transition at Ma = 0.10, Re = 4.0 × 106 , Tu = 0.15%, and α = 0◦ compared with the pressure profile from
the experiment [39].

                               Table 2      VA-2 structured multi-block O-grid dimensions.

                             grid level     chord × off-wall nodes                ∆s (chord)        avg/max y+
                                L0                541 × 121                       1.01 × 10−6        0.25 / 0.49
                                L1                761 × 171                       0.71 × 10−6        0.18 / 0.35
                                L2               1081 × 241                       0.50 × 10−6        0.13 / 0.25
                                L3               1521 × 341                       0.36 × 10−6        0.09 / 0.19

The results also demonstrate the potential benefits of using using adaptive mesh refinement and redistribution methods
for free-transition simulations. The blunt trailing-edge VA-2 and CRM-NLF cases presented in the following sections
are simulated using O-topology grids with trailing-edge aspect ratios on the order of 100. The need for small grid aspect
ratio at the trailing edge is not as acute for these cases.

B. VA-2 Supercritical Airfoil
    The VA-2 supercritical airfoil was recently investigated in the DNW-TWG wind tunnel by Costantini et al. [41] at a
Mach number of 0.72, Reynolds number of 10 × 106 , and angles of attack of -0.4 to 2.0 degrees. The upper surface
skin friction coefficient distributions were determined using a global luminescent oil-film skin friction field estimation
method (GLOFSFE), which measures the development of the thickness of an oil film, the distribution of which can be
used to calculate the skin friction, based on its luminescent intensity. The pressure distributions used for comparison in
the current work were obtained before applying the oil film. Mach number and angle of attack shifts of -0.01 Mach
and -0.30 degrees, respectively, were identified by Hebler et al. for experimental studies of the CAST10-2 airfoil in the
DNW-TWG wind tunnel [42], and are applied for the VA-2 simulations in the current work. The turbulence intensity
for the VA-2 experimental investigations in the DNW-TWG was not provided. Therefore, the turbulence intensity is
assumed to be consistent with previously estimated values for the DNW-TWG wind tunnel [43], with a value of 0.25%
used. This test case was previously investigated by the current authors using the SA-sLM2015 transition model both
with and without the Tollmien-Schlichting compressibility correction [11]. Structured O-type topology multi-block
grids were used to mesh the blunt trailing-edge VA-2 airfoil, with the grid dimensions presented in Table 2.
    Grid-refinement studies are presented for fully-turbulent flow and with free transition at two angles of attack,

<!-- ===== PDF page 10 ===== -->

                        103                                                                           103
                                                                       L0                                                                           L0
                        101                                            L1                             101                                           L1
Relative Residual, d

                                                                              Relative Residual, d
                        10−1                                           L2                             10−1                                          L2
                                                                       L3                                                                           L3
                        10−3                                                                          10−3
                        10−5                                                                          10−5
                        10−7                                                                          10−7
                        10−9                                                                          10−9
                   10−11                                                                         10−11
                   10−13
                      0.00                   0.25               0.50                             10−13
                                                                                                    0.00                    0.25             0.50
                               Equivalent Residual Evaluations 1e5                                           Equivalent Residual Evaluations 1e5
                                    (a) α = −0.40◦ : SA                                                      (b) α = −0.40◦ : SA-sLM2015cc
                        103                                                                           103
                                                                       L0                                                                           L0
                        101                                            L1                             101                                           L1
Relative Residual, d

                                                                              Relative Residual, d
                        10−1                                           L2                             10−1                                          L2
                                                                       L3                                                                           L3
                        10−3                                                                          10−3
                        10−5                                                                          10−5
                        10−7                                                                          10−7
                        10−9                                                                          10−9
                   10−11                                                                         10−11
                   10−130.0            0.5                1.0          1.5                       10−130.0             0.5           1.0             1.5
                               Equivalent Residual Evaluations 1e5                                           Equivalent Residual Evaluations 1e5
                                     (c) α = 1.80◦ : SA                                                       (d) α = 1.80◦ : SA-sLM2015cc

  Fig. 7 Residual convergence histories for the VA-2 airfoil simulated fully turbulent and with free transition at
  M = 0.71, Re = 10 × 106 , Tu = 0.25%, and α = −0.40◦, 1.80◦ .

-0.40 and 1.80 degrees. A shock forms on the upper surface of the airfoil at the higher angle of attack. The residual
convergence histories for the total residual for the fully-turbulent and free-transition grid-refinement studies are illustrated
in Figure 7. Similar to the NLF0416 results presented in Section III.A, all cases exit after hitting the absolute residual
convergence tolerance of 10−10 . The free-transition results demonstrate a significant increase in the computational cost
for the 1.80 degree angle of attack case, which features shock-induced transition on the upper surface of the airfoil at all
grid levels, relative to the -0.40 degrees angle of attack case that has a near-zero pressure gradient on the upper surface.
     A linearization study for the free-transition simulations at each angle of attack, -0.40 and 1.80 degrees, on the L3
grid is presented, with the total residual convergence histories illustrated in Figure 8. For the -0.40 degree case, with the
source-term time step restriction all four linearization strategies produce similar convergence behaviour. The number of
equivalent residual evaluations required to converge increases for each linearization strategy without the source-term
time step restriction, demonstrating that the effect of the source terms is more significant for this case relative to the
NLF0416 simulations. At 1.80 degrees angle of attack, with the source-term time step restriction the fully-coupled
linearization strategy provides improved nonlinear convergence relative to the loosely-coupled and decoupled approaches.
The fully-coupled linearization strategy provides a more complete Jacobian and produces a more accurate solution
update, which significantly improves convergence for this more challenging case. All solutions fail to converge within
the allocated compute time at 1.80 degrees angle of attack without the source-term time step restriction.
     The grid-convergence results for the fully-turbulent and free-transition simulations are presented in Figure 9. For

<!-- ===== PDF page 11 ===== -->

                        103                                                                        103
                                                    fully coupled                                                              fully coupled
                        101                         turb-trans coupled                             101                         turb-trans coupled
Relative Residual, d

                                                                           Relative Residual, d
                        10−1                        trans coupled                                  10−1                        trans coupled
                                                    decoupled                                                                  decoupled
                        10−3                                                                       10−3
                        10−5                                                                       10−5
                        10−7                                                                       10−7
                        10−9                                                                       10−9
                   10−11                                                                      10−11
                   10−130.0                   0.5                  1.0                        10−130.0                   0.5                  1.0
                               Equivalent Residual Evaluations 1e5                                        Equivalent Residual Evaluations 1e5
                                 (a) α = −0.40◦ : w/o STTS                                                   (b) α = −0.40◦ : w/ STTS
                        103                                                                        103
                                                    fully coupled                                                              fully coupled
                        101                         turb-trans coupled                             101                         turb-trans coupled
Relative Residual, d

                                                                           Relative Residual, d
                        10−1                        trans coupled                                  10−1                        trans coupled
                                                    decoupled                                                                  decoupled
                        10−3                                                                       10−3
                        10−5                                                                       10−5
                        10−7                                                                       10−7
                        10−9                                                                       10−9
                   10−11                                                                      10−11
                   10−130.0         0.5        1.0           1.5    2.0                       10−130.0         0.5        1.0           1.5    2.0
                               Equivalent Residual Evaluations 1e5                                        Equivalent Residual Evaluations 1e5
                                  (c) α = 1.80◦ : w/o STTS                                                   (d) α = 1.80◦ : w/ STTS

  Fig. 8 Linearization study residual convergence histories for the VA-2 airfoil simulated on the L3 grid (Table
  2) using the SA-sLM2015cc transition model at M = 0.71, Re = 10 × 106 , Tu = 0.25%, and α = −0.40◦, 1.80◦ .

 both angles of attack, the fully-turbulent and free-transition simulations produce similar lift coefficient grid convergence.
 At -0.40 degrees angle of attack, the free-transition results produce smaller AFE and GCI values and a higher order
 of convergence relative to the fully-turbulent results; however, at 1.80 degrees the results are reversed, with the
 free-transition results producing larger AFE and GCI values. At both angles of attack, the free-transition drag coefficient
 grid-convergence results demonstrate an increase in the discretization error produced on the finest grid level relative
 to the fully-turbulent results. The results demonstrate that additional mesh resolution is required for free-transition
 simulations in order to achieve a similar level of accuracy in the drag coefficient to the fully-turbulent simulations.

       The fully-turbulent and free-transition pressure and upper surface skin friction coefficient profiles are presented in
  Figure 10. The free-transition simulations more accurately reproduce the pressure coefficient profiles from experiment
  [41]. The increased extent of the laminar boundary pushes the location of the shock at 1.80 degrees angle of attack aft
  to better agree with experiment, with the shock location moving downstream with increasing grid refinement. For both
  angles of attack, the upper surface transition location moves downstream as the grid is refined. The difference in the
  upper surface transition onset locations between the L0 and L3 grids for both angles of attack is less than 3% chord, and
  is less than 1% chord between the L1 and L3 grids.

<!-- ===== PDF page 12 ===== -->

               AFE = 0.05%                     SA                                                                    SA
               GCI = 0.07%                     SA-sLM2015cc             0.360                                        SA-sLM2015cc
  .0089                                                                                 pL1-L3 = 1.35

                                                                        0.350
  .0088
                  pL1-L3 = 2.26
                                                                        0.340        AFE = 1.17%
  .0087                                                                              GCI = 1.48%
 Cd

          //                                                                    //

                                                                        Cl
               AFE = 0.18%
  .0060        GCI = 2.07%                                              0.320           pL1-L3 = 1.02

  .0059                                                                 0.310

  .0058           pL1-L3 = 3.44                                         0.300        AFE = 1.50%
                                                                                     GCI = 1.91%
                     5E-06            -1   1.5E-05                                         5E-06            -1   1.5E-05
                                  N                                                                     N
                                                            (a) α = −0.40◦

  .0104                                                                 0.810
               AFE = 0.12%                     SA                                                                    SA
               GCI = 0.14%                     SA-sLM2015cc                             pL1-L3 = 1.12                SA-sLM2015cc
  .0103                                                                 0.800

  .0102                                                                 0.790
                  pL1-L3 = 1.91                                                      AFE = 0.88%
                                                                                     GCI = 1.10%
  .0101                                                                 0.780
 Cd

          //                                                                    //
                                                                        Cl

                                                                        0.750
               AFE = 0.24%                                                              pL1-L3 = 1.10
  .0076        GCI = 0.30%
                                                                        0.740
  .0075
                                                                        0.730
                  pL1-L3 = 2.25                                                      AFE = 0.71%
  .0074                                                                              GCI = 0.90%
                                      -1                                0.720                               -1
                     5E-06        N        1.5E-05                                         5E-06        N        1.5E-05

                                                             (b) α = 1.80◦

Fig. 9 Grid-convergence results for the VA-2 airfoil simulated fully turbulent and with free transition at
M = 0.71, Re = 10 × 106 , Tu = 0.25%, and α = −0.40◦, 1.80◦ .

C. NASA CRM-NLF
    The NASA CRM-NLF was recently investigated in the NASA Langley National Transonic Facility (NTF) wind
tunnel [44, 45] and was the primary test case for the First AIAA Transition Modelling and Prediction Workshop† .
Transition visualizations, pressure coefficient profiles, and integrated forces were provided from the experiment over a
range of angles of attack [44, 45]. For the current work, the CRM-NLF is simulated at conditions similar to the 2524
run condition. Specifically, the wing-body geometry is investigated at a Mach number of 0.86, Reynolds number of
15 × 106 , 2◦ angle of attack, turbulence intensity of 0.24%, and with a surface roughness of 1.0µin. Four mesh levels are
investigated using structured multi-block grids following the gridding guidelines provided by the workshop committee‡ ,
with the grid details listed in Table 3.
   † https://transitionmodeling.larc.nasa.gov/wp-content/uploads/sites/109/2020/02/TransitionMPW_CaseDescriptions.pdf, accessed June 2021
   ‡ https://transitionmodeling.larc.nasa.gov/wp-content/uploads/sites/109/2020/02/CRM-NLF_GridGuidelines.pdf, accessed June 2021

<!-- ===== PDF page 13 ===== -->

                                      Costantini et al (2021)                                                           Costantini et al (2021)
   -1.5                               SAL0                          0.015         -1.5                                  Costantini et al (2021)       0.015
                                      SAL1                                                                              SA-sLM2015ccL0
                                      SAL2                                                                              SA-sLM2015ccL1
                                      SAL3                                                                              SA-sLM2015ccL2
       -1                                                                              -1                               SA-sLM2015ccL3

   -0.5                                                             0.01          -0.5                                                                0.01

                                                                        Cf,x

                                                                                                                                                          Cf,x
 Cp

                                                                                 Cp
       0                                                                               0

                                                                    0.005                                                                             0.005
      0.5                                                                             0.5

       1                                                                               1

                                                                    0                                                                                 0
      1.5                                                                             1.5
            0   0.2      0.4         0.6        0.8             1                           0   0.2        0.4         0.6        0.8             1
                               x/c                                                                               x/c
                      (a) α = −0.40◦ : SA                                                       (b) α = −0.40◦ : SA-sLM2015cc

                                      Costantini et al (2021)                                                           Costantini et al (2021)
   -1.5                               SAL0                          0.015         -1.5                                  Costantini et al (2021)       0.015
                                      SAL1                                                                              SA-sLM2015ccL0
                                      SAL2                                                                              SA-sLM2015ccL1
                                      SAL3                                                                              SA-sLM2015ccL2
       -1                                                                              -1                               SA-sLM2015ccL3

   -0.5                                                             0.01          -0.5                                                                0.01
                                                                        Cf,x

                                                                                                                                                          Cf,x
 Cp

                                                                                 Cp

       0                                                                               0

                                                                    0.005                                                                             0.005
      0.5                                                                             0.5

       1                                                                               1

                                                                    0                                                                                 0
      1.5                                                                             1.5
            0   0.2      0.4         0.6        0.8             1                           0   0.2        0.4         0.6        0.8             1
                               x/c                                                                               x/c
                      (c) α = 1.80◦ : SA                                                        (d) α = 1.80◦ : SA-sLM2015cc

Fig. 10 Pressure and upper surface skin friction coefficient profiles for the VA-2 airfoil simulated fully turbulent
and with free transition at M = 0.71, Re = 10 × 106 , Tu = 0.25%, and α = −0.40◦, 1.80◦ compared with the
pressure and skin friction coefficient profiles from the experiment [41].

                        Table 3      NASA CRM-NLF structured multi-block grid dimensions.

                               grid level     # of nodes                avg/max ∆s (chord)            avg/max y+
                                  L0          8,893,456                  1.00/2.24 × 10−6              0.36 / 1.50
                                  L1          16,691,200                 0.78/1.75 × 10−6              0.28 / 1.06
                                  L2          32,787,200                 0.60/1.35 × 10−6              0.22 / 0.78
                                  L3          64,330,000                 0.47/1.06 × 10−6              0.18 / 0.64

<!-- ===== PDF page 14 ===== -->

                        101                                                                    101
                                                                L0                                                                      L0
                        10−1                                    L1                             10−1                                     L1
Relative Residual, d

                                                                       Relative Residual, d
                        10−3                                    L2                             10−3                                     L2
                                                                L3                                                                      L3
                        10−5                                                                   10−5
                        10−7                                                                   10−7
                        10−9                                                                   10−9
                   10−11                                                                  10−11
                   10−13                                                                  10−13
                   10−150.0       0.5       1.0     1.5   2.0   2.5                       10−150.0       0.5    1.0     1.5       2.0   2.5
                               Equivalent Residual Evaluations 1e5                                    Equivalent Residual Evaluations 1e5
                                        (a) SA-QCR2000                                                 (b) SA-QCR2000-sLM2015cc

  Fig. 11 Residual convergence histories for the NASA CRM-NLF simulated fully turbulent and with free
  transition at approximately M = 0.86, Re = 15 × 106 , Tu = 0.24%, and α = 2.0◦ .

       Grid-refinement studies are presented for fully-turbulent flow and with free transition, with the residual convergence
  histories for the total residuals illustrated in Figure 11. All cases successfully meet either the absolute or relative residual
  convergence tolerances.
       A linearization study for the free-transition simulations on the L2 grid is presented, with the residual convergence
  histories illustrated in Figure 12. As demonstrated in Figure 12a, without the source-term time step restriction all
  linearization strategies fail to converge and instead stall after a 4 order of magnitude drop in the total residual. The
  fully-coupled linearization strategy without the source-term time step restriction exits after approximately 5 × 104
  equivalent residual evaluations after the solver is unable to prevent an unphysical solution update. However, with
  the source-term time step restriction (Figure 12b) the fully-coupled linearization strategy provides the best nonlinear
  convergence. The turb-trans coupled and trans coupled solutions stall after the total residual drops 5 orders of magnitude.
  The decoupled solution strategy appears to be converging, but exits due to a wall-clock constraint. The results in Figure
  12 demonstrate that the fully-coupled linearization strategy provides a more accurate solution update that improves
  nonlinear convergence; however, the source-term time step restriction is required in order to prevent unstable solution
  updates.
       The fully-turbulent and free-transition grid-convergence results are presented in Figure 13. While the fully-turbulent
  and free-transition simulations produce similar drag coefficient grid convergence, the free-transition lift coefficient grid
  convergence produces a larger GCI value. The high free-transition lift coefficient order of convergence demonstrates
  that the results are not in the asymptotic range of convergence [6], and indicates that additional grid-refinement levels
  are required in order to draw further conclusions.
       The pressure and skin friction coefficient profiles for the fully-turbulent and free-transition simulations are extracted
  at nine spanwise stations across the wing and compared with the experimental pressure profiles in Figures 14 and 15,
  respectively. For both the fully-turbulent and free-transition simulations, the shocks become better resolved as the
  mesh is refined. The local reduction in the upper surface skin friction coefficient produced by the upper surface shocks
  can be used to help differentiate the upper and lower surface skin friction coefficient profiles. For the free-transition
  results, the lower surface transition locations, which are predominantly produced by stationary crossflow instabilities,
  appear most sensitive to grid refinement at the η = 0.163 and η = 0.910 stations. The results demonstrate that the
  streamwise vorticity, which influences the growth of stationary crossflow instabilities, is better resolved with increased
  grid refinement. At the η = 0.163 station, the lower surface transition location moves upstream with refinement,
  indicating streamwise vorticity increases, and moves significantly downstream at the η = 0.910 station, demonstrating
  that the streamwise vorticity decreases. At all other stations, the difference between the lower surface transition locations
  between the L1 and L3 grids is less than 1%. The upper surface transition locations are most sensitive to grid refinement
  at the η = 0.730 and η = 0.910 stations, where the transition locations move as the shock is better resolved. Elsewhere,
  the difference between the upper surface transition locations between the L1 and L3 grids is less than 1%. The results
  again highlight the potential benefits of using adaptive mesh refinement and redistribution methods for these simulations

<!-- ===== PDF page 15 ===== -->

                        103                                                                             103
                                                          fully coupled                                 101
                                                                                                                                       fully coupled
                                                          turb-trans coupled                                                           turb-trans coupled
Relative Residual, d

                                                                                Relative Residual, d
                        10−1                              trans coupled                                 10−1                           trans coupled
                        10−3                              decoupled                                     10−3                           decoupled
                        10−5                                                                            10−5
                        10−7                                                                            10−7
                        10−9                                                                            10−9
                   10−11                                                                           10−11
                   10−13                                                                           10−13
                   10−150.0                 0.5               1.0        1.5                       10−150.0              0.5              1.0         1.5
                                 Equivalent Residual Evaluations 1e5                                           Equivalent Residual Evaluations 1e5
                                            (a) w/o STTS                                                                 (b) w/ STTS

 Fig. 12 Linearization study residual convergence histories for the NASA CRM-NLF simulated on the L2 grid
 (Table 3) using the SA-QCR2000-sLM2015cc transition model at approximately M = 0.86, Re = 15 × 106 ,
 Tu = 0.24%, and α = 2.0◦ .

                                                      SA-QCR2000                      0.430                                        SA-QCR2000
                                                      SA-QCR2000-sLM2015cc                                                         SA-QCR2000-sLM2015cc

                                                                                                                   pL1-L3 = 4.15
                                      pL1-L3 = 1.71                                   0.420
          .0187

                              AFE = 0.78%                                             0.410                AFE = 0.18%
                              GCI = 0.99%                                                                  GCI = 1.21%
            CD

                                                                                        CL

                                      pL1-L3 = 1.64                                   0.400                        pL1-L3 = 2.81
          .0182

                                                                                      0.390
                              AFE = 0.65%                                                                  AFE = 0.39%
                              GCI = 0.82%                                                                  GCI = 0.48%
          .0177                                                                       0.380
            0.0E+00                   1.0E-05             2.0E-05    3.0E-05            0.0E+00                    1.0E-05           2.0E-05     3.0E-05
                                                  N-2/3                                                                      N-2/3
  Fig. 13 Grid-convergence results for the NASA CRM-NLF simulated fully turbulent and with free transition
  at approximately M = 0.86, Re = 15 × 106 , Tu = 0.24%, and α = 2.0◦ .

  in order to improve grid resolution in these transition regions.

                                                                      IV. Conclusions
      The numerical behaviour of the SA-sLM2015 smooth local correlation-based transition model in an implicit
  Newton-Krylov flow solver is investigated through convergence and grid-refinement studies of the NLF0416 general
  aviation airfoil, VA-2 supercritical airfoil, and NASA CRM-NLF wing-body geometry. The results are evaluated
  relative to fully-turbulent simulations performed using the SA turbulence model. Four linearization strategies for the
  globalization phase of the free-transition simulations are investigated with varying levels of coupling for the turbulence
  and transition model equations. The linearization strategies are evaluated in conjunction with a source-term time step
  restriction in order to determine best-practices for solving the smooth transport-equation-based transition model in a

<!-- ===== PDF page 16 ===== -->

                                                                                    SA-QCR2000L0
                                                                                    SA-QCR2000L1
                                                                                    SA-QCR2000L2
                                                                                    SA-QCR2000L3
                                                                                    Lynde et al (2019)

                                                   0.009                                                       0.009                                                       0.009

  -0.5                                                        -0.5                                                        -0.5

                                                   0.006                                                       0.006                                                       0.006
       0                                               Cf,x        0                                                           0

                                                                                                                   Cf,x

                                                                                                                                                                               Cf,x
  Cp

                                                              Cp

                                                                                                                          Cp
   0.5                                                         0.5                                                         0.5
                                                   0.003                                                       0.003                                                       0.003

       1                                                           1                                                           1

                                                   0                                                           0                                                           0
           0   0.2     0.4         0.6   0.8   1                       0   0.2     0.4         0.6   0.8   1                       0   0.2     0.4         0.6   0.8   1
                             x/c                                                         x/c                                                         x/c

                     (a) η = 0.163                                               (b) η = 0.252                                               (c) η = 0.370

                                                   0.009                                                       0.009                                                       0.009

  -0.5                                                        -0.5                                                        -0.5

                                                   0.006                                                       0.006                                                       0.006
       0                                                           0                                                           0
                                                       Cf,x

                                                                                                                   Cf,x

                                                                                                                                                                               Cf,x
  Cp

                                                              Cp

   0.5                                                         0.5
                                                                                                                          Cp
                                                                                                                           0.5
                                                   0.003                                                       0.003                                                       0.003

       1                                                           1                                                           1

                                                   0                                                           0                                                           0
           0   0.2     0.4         0.6   0.8   1                       0   0.2     0.4         0.6   0.8   1                       0   0.2     0.4         0.6   0.8   1
                             x/c                                                         x/c                                                         x/c

                     (d) η = 0.460                                               (e) η = 0.550                                               (f) η = 0.640

                                                   0.009                                                       0.009                                                       0.009

  -0.5                                                        -0.5                                                        -0.5

                                                   0.006                                                       0.006                                                       0.006
       0                                                           0                                                           0
                                                       Cf,x

                                                                                                                   Cf,x

                                                                                                                                                                               Cf,x
  Cp

                                                              Cp

                                                                                                                          Cp

   0.5                                                         0.5                                                         0.5
                                                   0.003                                                       0.003                                                       0.003

       1                                                           1                                                           1

                                                   0                                                           0                                                           0
           0   0.2     0.4         0.6   0.8   1                       0   0.2     0.4         0.6   0.8   1                       0   0.2     0.4         0.6   0.8   1
                             x/c                                                         x/c                                                         x/c

                     (g) η = 0.730                                               (h) η = 0.820                                               (i) η = 0.910

Fig. 14 Pressure and skin friction coefficient profiles for the CRM-NLF simulated using the SA-QCR2000
turbulence model at approximately M = 0.86, Re = 15 × 106 , Tu = 0.24%, and α = 2.0◦ compared with the
pressure profiles from the experiment [44].

<!-- ===== PDF page 17 ===== -->

                                                                                 SA-QCR2000-sLM2015ccL0
                                                                                 SA-QCR2000-sLM2015ccL1
                                                                                 SA-QCR2000-sLM2015ccL2
                                                                                 SA-QCR2000-sLM2015ccL3
                                                                                 Lynde et al (2019)

                                                   0.009                                                       0.009                                                       0.009

  -0.5                                                        -0.5                                                        -0.5

                                                   0.006                                                       0.006                                                       0.006
       0                                               Cf,x        0                                                           0

                                                                                                                   Cf,x

                                                                                                                                                                               Cf,x
  Cp

                                                              Cp

                                                                                                                          Cp
   0.5                                                         0.5                                                         0.5
                                                   0.003                                                       0.003                                                       0.003

       1                                                           1                                                           1

                                                   0                                                           0                                                           0
           0   0.2     0.4         0.6   0.8   1                       0   0.2     0.4         0.6   0.8   1                       0   0.2     0.4         0.6   0.8   1
                             x/c                                                         x/c                                                         x/c

                     (a) η = 0.163                                               (b) η = 0.252                                               (c) η = 0.370

                                                   0.009                                                       0.009                                                       0.009

  -0.5                                                        -0.5                                                        -0.5

                                                   0.006                                                       0.006                                                       0.006
       0                                                           0                                                           0
                                                       Cf,x

                                                                                                                   Cf,x

                                                                                                                                                                               Cf,x
  Cp

                                                              Cp

   0.5                                                         0.5
                                                                                                                          Cp
                                                                                                                           0.5
                                                   0.003                                                       0.003                                                       0.003

       1                                                           1                                                           1

                                                   0                                                           0                                                           0
           0   0.2     0.4         0.6   0.8   1                       0   0.2     0.4         0.6   0.8   1                       0   0.2     0.4         0.6   0.8   1
                             x/c                                                         x/c                                                         x/c

                     (d) η = 0.460                                               (e) η = 0.550                                               (f) η = 0.640

                                                   0.009                                                       0.009                                                       0.009

  -0.5                                                        -0.5                                                        -0.5

                                                   0.006                                                       0.006                                                       0.006
       0                                                           0                                                           0
                                                       Cf,x

                                                                                                                   Cf,x

                                                                                                                                                                               Cf,x
  Cp

                                                              Cp

                                                                                                                          Cp

   0.5                                                         0.5                                                         0.5
                                                   0.003                                                       0.003                                                       0.003

       1                                                           1                                                           1

                                                   0                                                           0                                                           0
           0   0.2     0.4         0.6   0.8   1                       0   0.2     0.4         0.6   0.8   1                       0   0.2     0.4         0.6   0.8   1
                             x/c                                                         x/c                                                         x/c

                     (g) η = 0.730                                               (h) η = 0.820                                               (i) η = 0.910

Fig. 15 Pressure and skin friction coefficient profiles for the CRM-NLF simulated using the SA-QCR2000-
sLM2015cc transition model at approximately M = 0.86, Re = 15 × 106 , Tu = 0.24%, and α = 2.0◦ compared
with the pressure profiles from the experiment [44].

<!-- ===== PDF page 18 ===== -->

strong implicit solver.
    The linearization study residual convergence histories demonstrate that both the fully-coupled linearization strategy
and the source-term time step restriction improve nonlinear convergence as the flow complexity increases. Without the
source-term time step restriction, the fully-coupled linearization strategy produces unstable solution updates. This is
particularly relevant for the test cases featuring shock-induced transition. The results from the grid-convergence studies
demonstrate that while the free-transition simulations produce similar grid-convergence behaviour to fully-turbulent
simulations, the free-transition simulations in general require additional streamwise grid resolution in order to achieve a
similar level of accuracy, primarily driven by the need for sufficient streamwise resolution in the transition region.

                                                      Acknowledgments
    This work was partially funded by the NASA Transformational Tools and Technologies (TTT) and Advanced Air
Transport Technology (AATT) projects, the Natural Sciences and Engineering Research Council (NSERC), and the
University of Toronto. Computations were performed on the Niagara supercomputer at the SciNet HPC Consortium, a
part of Compute Canada. The authors gratefully acknowledge the input and support provided by Dr. Ercan Dumlupinar,
Dr. Cetin Kiris, and the rest of the LAVA team.

                                                           References
 [1] Langtry, R. B., and Menter, F. R., “Correlation-based transition modeling for unstructured parallelized computational fluid
     dynamics codes,” AIAA Journal, Vol. 47, No. 12, 2009, pp. 2894–2906. doi:https://doi.org/10.2514/1.42362.
 [2] Medida, S., and Baeder, J. D., “Application of the Correlation-based γ- Re        ˜ θt Transition Model to the Spalart-Allmaras
     Turbulence Model,” 20th AIAA Computational Fluid Dynamics Conference, AIAA Paper 2011-3979, Honolulu, HI, 2011.
     doi:https://doi.org/10.2514/6.2011-3979.
 [3] Coder, J. G., and Maughmer, M. D., “Computational fluid dynamics compatible transition modeling using an amplification
     factor transport equation,” AIAA Journal, Vol. 52, No. 11, 2014, pp. 2506–2512. doi:https://doi.org/10.2514/1.J052905.
 [4] Menter, F. R., Smirnov, P. E., Liu, T., and Avancha, R., “A One-Equation Local Correlation-Based Transition Model,” Flow,
     Turbulence and Combustion, Vol. 95, 2015, pp. 583–619. doi:https://doi.org/10.1007/s10494-015-9622-4.
 [5] Roache, P. J., “Perspective: A Method for Uniform Reporting of Grid Refinement Studies,” Journal of Fluids Engineering, Vol.
     116, No. 3, 1994, pp. 405–413. doi:https://doi.org/10.1115/1.2910291.
 [6] Eça, L., and Hoekstra, M., “Evaluation of numerical error estimation based on grid refinement studies with
     the method of the manufactured solutions,” Computers & Fluids, Vol. 38, No. 8, 2009, pp. 1580–1591.
     doi:https://doi.org/10.1016/j.compfluid.2009.01.003.
 [7] Zingg, D. W., and Godin, P., “A perspective on turbulence models for aerodynamic flows,” International Journal of Computational
     Fluid Dynamics, Vol. 23, No. 4, 2009, pp. 327–335. doi:https://doi.org/10.1080/10618560902776802.
 [8] Lopes, R., Eça, L., Vaz, G., and Kerkvliet, M., “Assessing Numerical Aspects of Transitional Flow Simulations Using
     the RANS Equations,” International Journal of Computational Fluid Dynamics, Published Online, January 14, 2021.
     doi:https://doi.org/10.1080/10618562.2020.1870962.
 [9] Spalart, P. R., and Allmaras, S. R., “A One-Equation Turbulence Model for Aerodynamic Flows,” 30th AIAA Aerospace Sciences
     Meeting and Exhibit, AIAA Paper 092-0439, Reno, Nevada, United States 1992. doi:https://doi.org/10.2514/6.1992-439.
[10] Piotrowski, M. G. H., and Zingg, D. W., “Smooth Local Correlation-based Transition Model for the Spalart-Allmaras Turbulence
     Model,” AIAA Journal, Vol. 59, No. 2, 2021, pp. 474–492. doi:https://doi.org/10.2514/1.J059784.
[11] Piotrowski, M. G. H., and Zingg, D. W., “Compressibility Corrections to Extend a Smooth Local Correlation-Based Transition
     Model to Transonic Flows,” AIAA Journal, Submitted July, 2021.
[12] Moryossef, Y., and Levy, Y., “Unconditionally positive implicit procedure for two-equation turbulence models:
     Application to k-ω turbulence models,” Journal of Computational Physics, Vol. 220, 2006, pp. 88–108.
     doi:https://doi.org/10.1016/j.jcp.2006.05.001.
[13] Lian, C., Xia, G., and Merkle, C. L., “Impact of source terms on reliability of CFD algorithms,” Computers and Fluids Journal,
     Vol. 39, 2010, pp. 1909–1922. doi:https://doi.org/10.1016/j.compfluid.2010.06.021.
[14] Allmaras, S. R., and Johnson, F. T., “Modifications and clarifications for the implementation of the Spalart-Allmaras turbulence
     model,” Seventh International Conference on Computational Fluid Dynamics, ICCFD7 Paper 1902, 2012, pp. 1–11.
[15] Ilinca, F., and Pelletier, D., “Positivity preservation and adaptive solution of two-equation models of turbulence,” International
     Journal of Thermal Sciences, Vol. 38, No. 7, 1999, pp. 560–571. doi:https://doi.org/10.1016/S0035-3159(99)80036-1.
[16] Coder, J. G., “Further Development of the Amplification Factor Transport Transition Model for Aerodynamic Flows,” 2019
     AIAA Aerospace Sciences Meeting, AIAA Paper 2019-0039, 2019. doi:https://doi.org/10.2514/6.2019-0039.
[17] Mosahebi, A., and Laurendeau, E., “Convergence Characteristics of Fully and Loosely Coupled Numerical Approaches for
     Transition Models,” AIAA Journal, Vol. 53, No. 5, 2015, pp. 1399–1404. doi:https://doi.org/10.2514/1.J053722.
[18] Mosahebi, A., and Laurendeau, E., “Introduction of a modified segregated numerical approach for efficient simulation

<!-- ===== PDF page 19 ===== -->

     of γ- Re
            ˜ θt transition model,” International Journal of Computational Fluid Dynamics, Vol. 29, 2015, pp. 357–375.
     doi:https://doi.org/10.1080/10618562.2015.1093624.
[19] Menter, F. R., “Two-Equation Eddy-Viscosity Turbulence Models for Engineering Applications,” AIAA Journal, Vol. 32, No. 8,
     1994, pp. 1598–1605. doi:https://doi.org/10.2514/3.12149.
[20] Knoll, D., and Keyes, D., “Jacobian-free Newton-Krylov methods: a survey of approaches and applications,” Journal of
     Computational Physics, Vol. 193, No. 2, 2004, pp. 357–397. doi:https://doi.org/10.1016/j.jcp.2003.08.010.
[21] Blanco, M., and Zingg, D. W., “Newton-Krylov Algorithm with a Loosely Coupled Turbulence Model for Aerodynamic Flows,”
     AIAA Journal, Vol. 45, No. 5, 2007, pp. 980–987. doi:https://doi.org/10.2514/1.22972.
[22] Chisholm, T. T., and Zingg, D. W., “A Jacobian-free Newton-Krylov algorithm for compressible turbulent fluid flows,” Journal
     of Computational Physics, Vol. 228, No. 9, 2009, pp. 3490–3507. doi:https://doi.org/10.1016/j.jcp.2009.02.004.
[23] Osusky, M., and Zingg, D. W., “Parallel Newton-Krylov-Schur Flow Solver for the Navier-Stokes Equations,” AIAA Journal,
     Vol. 51, No. 12, 2013, pp. 2833–2851. doi:https://doi.org/10.2514/1.J052487.
[24] Yildirim, A., Kenway, G. K. W., Mader, C. A., and Martins, J. R. R. A., “A Jacobian-free approximate Newton-Krylov startup strat-
     egy for RANS simulations,” Journal of Computational Physics, Vol. 397, 2019. doi:https://doi.org/10.1016/j.jcp.2019.06.018.
[25] Lopes, R., Eça, L., and Vaz, G., “On the Numerical Behavior of RANS-Based Transition Models,” Journal of Fluids Engineering,
     Vol. 142, No. 5, 2020. doi:https://doi.org/10.1115/1.4045576.
[26] Hicken, J. E., and Zingg, D. W., “Parallel Newton-Krylov solver for the Euler equations discretized using simultaneous
     approximation terms,” AIAA Journal, Vol. 46, No. 11, 2008, pp. 2773–2786. doi:https://doi.org/10.2514/1.34810.
[27] Swanson, R. C., and Turkel, E., “On central-difference and upwind schemes,” Journal of Computational Physics, Vol. 101,
     1992, pp. 292–306. doi:https://doi.org/10.1016/0021-9991(92)90007-L.
[28] Osusky, M., Boom, P. D., and Zingg, D. W., “Results from the Fifth AIAA Drag Prediction Workshop obtained with a
     parallel Newton-Krylov-Schur flow solver discretized using summation-by-parts operators,” 31st AIAA Applied Aerodynamics
     Conference, AIAA 2013-2511, 2013. doi:10.2514/6.2013-2511.
[29] Saad, Y., and Schultz, M. H., “GMRES: A generalized minimal residual algorithm for solving nonsymmetric linear systems,”
     SIAM Journal on Scientific and Statistical Computing, Vol. 7, No. 3, 1986, pp. 856–869. doi:https://doi.org/10.1137/0907058.
[30] Saad, Y., and Sosonkina, M., “Distributed Schur complement techniques for general sparse linear systems,” SIAM Journal of
     Scientific Computing, Vol. 21, 1999, pp. 1337–1357. doi:https://doi.org/10.1137/S1064827597328996.
[31] Langtry, R. B., Sengupta, K., Yeh, D. T., and Dorgan, A. J., “Extending the γ- Re     ˜ θt Correlation based Transition Model for
     Crossflow Effects,” 45th AIAA Fluid Dynamics Conference, AIAA Paper 2015-2474, 2015. doi:https://doi.org/10.2514/6.2015-
     2474.
[32] Modisette, J. M., “An automated reliable method for two-dimensional Reynolds-Averaged Navier-Stokes simulations,” Ph.D.
     thesis, Massachusetts Institute of Technology, 2011.
[33] Zingg, D. W., “Viscous airfoil computations using Richardson extrapolation,” 10th Computational Fluid Dynamics Conference,
     AIAA 1991-1559, 1991. doi:https://doi.org/10.2514/6.1991-1559.
[34] Zingg, D. W., “Grid studies for thin-layer Navier-Stokes computations of airfoil flowfields,” AIAA Journal, Vol. 30, No. 10,
     1992, pp. 2561–2564. doi:https://doi.org/10.2514/3.11265.
[35] Diskin, B., Thomas, J. L., Rumsey, C. L., and Schwöppe, A., “Grid Convergence for Turbulent Flows (Invited),” 53rd AIAA
     Aerospace Sciences Meeting, AIAA 2015-1746, 2015. doi:https://doi.org/10.2514/6.2015-1746.
[36] Celik, I. B., Ghia, U., Roache, P. J., Freitas, C. J., Coleman, H., and Raad, P. E., “Procedure for Estimation and Reporting
     of Uncertainty Due to Discretization in CFD Applications,” Journal of Fluids Engineering, Vol. 130, No. 7, July 2008.
     doi:https://doi.org/10.1115/1.2960953.
[37] Osusky, L., Buckley, H. P., Reist, T. A., and Zingg, D. W., “Drag Minimization Based on the Navier-Stokes Equations Using a
     Newton-Krylov Approach,” AIAA Journal, Vol. 53, No. 6, 2015, pp. 1555–1577. doi:https://doi.org/10.2514/1.J053457.
[38] Spalart, P. R., “Strategies for Turbulence Modelling and Simulation,” International Journal of Heat and Fluid Flow, Vol. 21,
     No. 3, 2000, pp. 252–263. doi:https://doi.org/10.1016/S0142-727X(00)00007-2.
[39] Somers, D. M., “Design and experimental results for a natural-laminar-flow airfoil for general aviation applications,” NASA
     TP-1861, 1981.
[40] Coder, J. G., “Standard Test Cases for Transition Model Verification and Validation in Computational Fluid Dynamics,” 2018
     AIAA Aerospace Sciences Meeting, AIAA Paper 2018-0029, 2018. doi:https://doi.org/10.2514/6.2018-0029.
[41] Costantini, M., Lee, T., Nonomura, T., Asai, K., and Klein, C., “Feasibility of skin-friction field measurements in a transonic wind
     tunnel using a global luminescent oil film,” Experiments in Fluids, Vol. 62, No. 1, 2021, p. 21. doi:https://doi.org/10.1007/s00348-
     020-03109-z.
[42] Hebler, A., Schojda, L., and Mai, H., “Experimental Investigation of the Aeroelastic Behavior of a Laminar Airfoil in Transonic
     Flow,” Proc. IFASD 2013, IFASD, Bristol, 2013. doi:https://elib.dlr.de/83897/.
[43] Fehrs, M., van Rooij, A. C., and Nitzsche, J., “Influence of Boundary Layer Transition on the Flutter Behavior of a Supercritical
     Airfoil,” CEAS Aeronautical Journal, Vol. 6, No. 2, 2015, pp. 291–303. doi:https://doi.org/10.1007/s13272-014-0147-7.
[44] Lynde, M. N., Campbell, R. L., Rivers, M. B., Viken, S. A., Chan, D. T., Watkins, N. A., and Goodliff, S. L., “Preliminary
     Results from an Experimental Assessment of a Natural Laminar Flow Design Method,” 2019 AIAA Aerospace Sciences Meeting,
     AIAA Paper 2019-2298, 2019. doi:https://doi.org/10.2514/6.2019-2298.

<!-- ===== PDF page 20 ===== -->

[45] Lynde, M. N., Campbell, R. L., and Viken, S. A., “Additional Findings from the Common Research Model Natural Laminar
     Flow Wind Tunnel Test,” AIAA Aviation 2019 Forum, AIAA Paper 2019-3292, 2019. doi:https://doi.org/10.2514/6.2019-3292.
