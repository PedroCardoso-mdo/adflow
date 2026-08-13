# Piotrowski & Zingg (2022) — Investigation of a Smooth Local Correlation-based Transition Model in a Discrete-Adjoint Aerodynamic Shape Optimization Algorithm

> **Full transcript** — verbatim text extracted from the source PDF with
> `pdftotext -layout`, page-delimited, nothing removed except bare page-number
> lines. Layout is preserved, so tables and multi-column material read as they
> appear in the PDF. **Equations are raw PDF extractions**: superscripts,
> subscripts and math symbols may be flattened or split across lines — for any
> equation you intend to implement, confirm against the PDF itself (path below)
> or the digest companion file.
>
> Source PDF: `/home/mdo/Desktop/ARTIGOS SMOTH/PiotrowskiandZingg2022 (2).pdf`

---

<!-- ===== PDF page 1 ===== -->

 Investigation of a Smooth Local Correlation-based Transition
 Model in a Discrete-Adjoint Aerodynamic Shape Optimization
                          Algorithm

                                                Michael G. H. Piotrowski∗
                                   NASA Ames Research Center, Moffett Field, CA 94035
                           Institute for Aerospace Studies, University of Toronto, ON M3H 5T6

                                                     David W. Zingg†
                           Institute for Aerospace Studies, University of Toronto, ON M3H 5T6

             A smooth local correlation-based transition model is fully coupled to a RANS-based Newton-
         Krylov flow solver and discrete-adjoint gradient-based optimization algorithm. The free-
         transition optimization framework is evaluated using lift-constrained drag minimizations of
         airfoils at design conditions ranging from light to single-aisle aircraft and an infinite swept wing
         at design conditions representative of a transonic strut-braced wing aircraft. The impact of
         the streamwise grid resolution on the ability of the optimization algorithm to delay boundary-
         layer transition is investigated, with the results demonstrating that streamwise grid resolution
         requirements increase as the transition length decreases with increasing Reynolds number. The
         optimization problem at the light aircraft design conditions is demonstrated to be multi-modal,
         with the optimization algorithm producing two distinct designs: one with a thin, reflexed
         trailing edge and steep pressure recovery regions, the other with increased aft loading, with
         the latter design outperforming the former. A drag minimization of an airfoil at transonic
         design conditions demonstrates that the optimization algorithm successfully trades a decrease
         in viscous drag by delaying boundary-layer transition with an increase in wave drag, while the
         drag minimization of an infinite swept wing demonstrates the capability of the optimization
         algorithm to delay both Tollmien-Schlichting and stationary crossflow instabilities.

                                                      I. Introduction
    Laminar-flow wing designs produce an increased laminar extent of the boundary layer, resulting in a decrease
in viscous drag [1]. These drag savings have the potential to significantly reduce the fuel burn of transport aircraft,
where viscous drag constitutes approximately 50% of the total drag [2]. The potential of laminar-flow designs has been
demonstrated by various studies which suggest the application of laminar flow control to large commercial aircraft can
reduce aerodynamic drag by approximately 10% [3].
    In the near term, natural laminar flow (NLF) is being investigated for winglet, tail, and nacelle designs, with the
Boeing 787-8 and 777X representing the first commercial applications of NLF nacelles to large transport aircraft, and
the 737 Max aircraft exploiting laminar flow for both the nacelle and winglet designs [4]. In the future, the design
of wings which exploit significant regions of laminar flow can result in more significant drag reduction. The Airbus
BLADE project is investigating this, while experimental studies of the NASA common research model NLF (CRM-NLF)
demonstrate that it is possible to design transport wings with high sweep and Reynolds numbers with significant regions
of laminar flow [5–7].
    Efficient high-fidelity numerical optimization algorithms provide the designers of next generation aircraft with
powerful tools for the development of more fuel-efficient designs [8–11]. The integration of transition prediction
methodologies into aerodynamic shape optimization algorithms enables the study of various trade-offs in the design
of novel NLF wing configurations. Although there are several examples in the literature of researchers combining
transition prediction with aerodynamic optimization using lower-fidelity models [12–17], or with gradient-free or costly
    ∗ NASA Pathways Ph.D. Candidate, University of Toronto, NASA Computational Aeroscience Branch, AIAA Student Member,

m.piotrowski@mail.utoronto.ca
    † University of Toronto Distinguished Professor of Computational Aerodynamics and Sustainable Aviation, Associate Fellow AIAA,

dwz@oddjob.utias.utoronto.ca

<!-- ===== PDF page 2 ===== -->

finite-difference gradient approximations [18–24], the focus of the current review will be limited to methods using
RANS-based analysis with adjoint-based shape optimization algorithms.
     The work by Driver and Zingg [25] represents an early example of the use of RANS-based optimization to investigate
the design of NLF airfoils. The coupled inviscid-viscous boundary-layer solver MSES [26] was used to calculate
transition locations, which were enforced using the trip terms in the Spalart-Allmaras (SA) turbulence model [27].
Through the application of the resulting algorithm to single-point optimizations, Driver and Zingg [25] produced NLF
airfoil designs with laminar-flow design characteristics similar to those designed by Liebeck [28] and Zingg [29].
     Lee and Jameson [30] coupled a boundary-layer solver and database e N method, developed by Kroo and Sturdza
[16] for the design of supersonic laminar flow wings, to the Baldwin-Lomax algebraic turbulence model [31] in a
RANS-based solver and discrete-adjoint gradient-based optimization algorithm. Although their work demonstrated
the significance of laminar flow on NLF wing design, the gradients did not include transition prediction. Instead, the
optimizations were focused on reducing wave drag in order to increase aerodynamic performance.
     Rashad and Zingg [32] performed single- and multi-point optimizations at subsonic and transonic design conditions
using the simplified e N envelope method developed by Drela and Giles [26] coupled to the SA turbulence model [27].
Rashad conducted a thorough investigation of coupling strategies for the RANS-e N flow solver and demonstrated the
importance of providing a smooth ramp-up of the eddy-viscosity, as well as using a tight tolerance for the transition
residual, in order to ensure a smooth design space [32, 33]. The coupled-adjoint system is solved using preconditioned
GMRES with a non-iterative solution strategy. Single-point optimizations demonstrate that the resulting algorithm is
capable of extending the laminar boundary layer significantly, with multi-point optimizations producing NLF airfoil
designs robust to changes in the flow conditions and disturbance environments. The drag minimizations performed
by Rashad and Zingg [32] represent suitable two-dimensional reference benchmarks for evaluating NLF optimization
frameworks.
     Shi et al. [34] used a simplified e N method based on a database of linearized stability theory (LST) results in a
Jacobian-free discrete-adjoint optimization algorithm. Similar to the work by Rashad [33] and Rashad and Zingg
[32], Shi et al. [34] emphasized the importance of using a smooth intermittency function to couple the simplified e N
method with the RANS solutions, as well as a tight transition residual tolerance. Shi et al. investigated single-point
airfoil optimizations at the Cessna 172R design conditions initially proposed by Rashad and Zingg [32], producing an
optimized NLF design similar to that developed by Rashad and Zingg, and multi-point optimizations of an airfoil at
design conditions representative of the HondaJet aircraft [35]. Their work was recently extended to infinite swept wing
optimizations using the C1 criterion [36].
     Zhu and Qin [37] applied the simplified e N envelope method by Drela and Giles [26], extended with the crossflow
criterion developed by Kroo and Sturdza [16], and using Poll’s criterion [38] as a constraint to prevent attachment-line
transition, to the optimization of infinite swept wings using a discrete-adjoint optimization algorithm. Similar to Lee
and Jameson [30], the transition locations from stability analysis were enforced using the Baldwin-Lomax algebraic
turbulence model [31]. Drag-minimizations of infinite swept wings at transonic flight conditions were performed
both with and without a shock-control bump. The shock-control bump reduces wave drag, enabling the optimization
algorithm to reduce sweep, which attenuates the crossflow instabilities, moving the transition front downstream. A
systems-level benefit assessment was performed demonstrating that a low-sweep NLF wing design with shock-control
bumps applied to a narrow-body commercial aircraft can reduce fuel consumption by 11.1% [39].
     Methods based on stability analysis require non-local boundary-layer information and a large infrastructure to apply
the code, usually consisting of a boundary-layer solver and an e N method. While previous work has successfully
demonstrated that methods for automatic transition prediction using stability analysis can be developed for complex
three-dimensional configurations [40, 41], methods based on local transition criteria, such as the Langtry-Menter
transition model (LM2009) [42], can be more easily integrated into a RANS-based flow solver and gradient-based
optimization algorithm.
     Khayatzadeh and Nadarajah [43, 44] coupled the Langtry-Menter LM2009 transition model [42] with a RANS-based
flow solver and discrete-adjoint optimization algorithm. Through the application of the resulting algorithm to subsonic
lift-constrained drag-minimization and lift-to-drag ratio maximization optimizations, they demonstrated the importance
of fully coupling the transition model equations to the adjoint system, and successfully produced designs that delayed
the onset of transition.
     Halila et al. [45] integrated a smooth version of the amplification factor transport (AFT) transition model [46],
coupled to the SA turbulence model, with an approximate Newton-Krylov flow solver and discrete-adjoint optimization
algorithm. The resulting algorithm was applied to single- and multi-point subsonic and single-point transonic airfoil
drag minimizations.

<!-- ===== PDF page 3 ===== -->

    Yang and Mavriplis [47] integrated the AFT transition model [46] with a Newton-Krylov flow solver and discrete-
adjoint framework and applied the algorithm to two-dimensional subsonic optimizations. Following this work, Mavripilis
et al. [48] applied the one-equation intermittency transition model developed by Menter et al. [49], coupled to the SA
turbulence model [27], in a discrete-adjoint optimization algorithm to a slotted transonic truss-braced wing (TTBW)
geometry.
    In general, optimized designs resulting from algorithms combining gradient-based optimization with local transport-
equation-based transition models tend to underperform NLF designs produced using stability-analysis approaches
coupled with gradient-based optimization algorithms. Specifically, gradient-based optimization algorithms with
transport-equation-based transition models struggle to produce designs with significantly increased regions of laminar
flow [43, 44, 47, 50], especially at higher Mach and Reynolds numbers [45]. This is particularly evident when comparing
with the designs produced by Rashad and Zingg [32] that achieve significant regions of laminar flow at high Mach
numbers and flight-scale Reynolds numbers.
    The goal of this work it to evaluate the capabilities of a discrete-adjoint gradient-based optimization algorithm
coupled with a smooth local correlation-based transition model. To achieve this goal, the SA-sLM2015 transition model
[51–53] is fully coupled into a RANS-based Newton-Krylov flow solver [51, 53] and discrete-adjoint gradient-based
optimization algorithm, with details and verification of the gradients of the resulting algorithm presented in Section II.
The free-transition optimization framework is investigated using airfoil drag minimizations at design conditions ranging
from subsonic, representative of light aircraft, up to transonic conditions representative of single-aisle aircraft, and
an infinite swept wing drag minimization at design conditions representative of a transonic strut-braced wing aircraft.
The influence of the streamwise resolution of the grids and constraints applied to the design space on the behaviour of
the optimization algorithm is investigated, with the results presented in Section III. Conclusions and future work are
presented in Section IV.

                                                  II. Methodology
    The high-fidelity discrete-adjoint gradient-based aerodynamic shape optimization framework used in the current
work, Jetstream, consists of three primary components: a geometry parameterization and control scheme [54–56], a
parallel structured multi-block Newton-Krylov-Schur RANS-based flow solver [51, 53, 57, 58], and a discrete-adjoint
gradient-based optimization algorithm [54, 59]. Jetstream was recently cross-validated with an industry RANS-based
discrete-adjoint optimization algorithm, with the two algorithms producing very similar geometries and performance
improvements [60]. The flow solver, Diablo, was extensively validated with fully-turbulent flow as part of the Fifth
AIAA Drag Prediction Workshop [61], and has recently been extended to perform free-transition simulations using the
SA-sLM2015 transition model in subsonic and transonic flow regimes [51–53]. This section will provide an overview of
these components with an emphasis on the flow solver and gradient evaluation, which were modified for the introduction
of the SA-sLM2015 transition model equations.

A. Flow Solution Algorithm
    The flow solution algorithm used in the current work is a three-dimensional structured multi-block finite-difference
solver [57, 58] which solves the RANS equations fully-coupled to the three-equation SA-sLM2015 transition model
[51–53]. The governing equations are spatially discretized using summation-by-parts operators, with simultaneous
approximation terms applied to enforce boundary conditions and inter-block coupling. The mean-flow equations are
discretized using second-order summation-by-parts operators with artificial dissipation provided using the matrix-based
dissipation model developed by Swanson and Turkel [62]. First-order upwinding is applied to the turbulence and
transition model equations. The computational domain is decomposed into multiple blocks, resulting in multi-block
structured grids, which allow for efficient parallel computations. A fully-coupled, implicit Newton-Krylov-Schur
solution algorithm making use of a pseudo-transient continuation strategy is applied to the set of discretized equations
to drive the residual to a converged steady-state solution. Details of the modifications made to the Newton-Krylov-Schur
algorithm to integrate the transition model equations can be found in [51, 53].
    The work by Rashad [33] and Rashad and Zingg [32] emphasizes the importance of ensuring a smooth design
space for gradient-based optimization. However, many local transport-based transition models [42, 46, 49, 63]
contain discontinuous and non-differentiable functions. In previous work, the current authors developed the smooth,
three-equation SA-sLM2015 transition model [51], which is based on the LM2009 [42] empirical correlations for bypass
transition and transition due to Tollmien-Schlichting instabilities as well as the LM2015 [64] helicity-based empirical
correlations for stationary crossflow instabilities. The empirical correlations in the SA-sLM2015 transition model were

<!-- ===== PDF page 4 ===== -->

recently extended to transonic flow regimes [52].

B. Gradient-based Optimization

1. Optimization Algorithm
    The optimization algorithm consists of an integrated geometry parametrization and mesh movement scheme
developed by Hicken and Zingg [54] that represents the initial geometry by a set of B-spline surface patches. Shape
control is achieved using free-form deformation (FFD) B-spline volumes with an axial deformation geometry control
system [55, 56]. Optimization is performed using a gradient-based optimization strategy where the gradients are
calculated using the discrete-adjoint approach [54, 59]. The flow adjoint system is solved using a simplified and flexible
variant of GCROT [65], while the mesh adjoint system is solved using a preconditioned conjugate-gradient method.
The design variables are updated using the sparse quadratic programming algorithm SNOPT [66], which is capable
of handling both linear and non-linear constraints and has been demonstrated to be efficient for problems with large
numbers of design variables.

2. Gradient Evaluation
    The discrete-adjoint Lagrangian function for the PDE-constrained optimization problem containing the flow, ψ f ,
and mesh adjoint variables, ψ m , is given by,

                           L(Q, C, X, ψ f, ψ m ) = J (Q, C, X) − ψ Tf (R(Q, C, X)) − ψ Tm M(C, X),                     (1)

where Q, C, and X represent the array of conserved flow variables, B-spline control points, and design variables, and R
and M represent the total residual, which includes the RANS, turbulence, and transition model equations, and the mesh
movement equations, respectively. Setting the first derivatives of the Lagrangian with respect to Q, C, X, ψ f , and ψ m to
zero, we recover the following first-order optimality conditions,
                                           ∂L
                                                = 0 = M(C, X),                                                         (2)
                                          ∂ψ m
                                           ∂L
                                                = 0 = R(Q, C, X),                                                      (3)
                                           ∂ψ f
                                            ∂L        ∂J        ∂R
                                                =0=      − ψ Tf    ,                                                   (4)
                                           ∂Q         ∂Q        ∂Q
                                            ∂L        ∂J        ∂R         ∂M
                                                =0=      − ψ Tf     − ψ Tm    ,                                        (5)
                                            ∂C        ∂C        ∂C         ∂C
                                            ∂L        ∂J        ∂R         ∂M
                                                =0=      − ψ Tf     − ψ Tm    ,                                        (6)
                                            ∂X        ∂X        ∂X         ∂X
where in the current work the partial derivatives are formed analytically using hand linearization with some entries
determined using the complex-step method.
     To perform free-transition optimization, the partial derivatives of the total residual must be updated. Specifically,
                      ∂R                              ∂R
the flow Jacobian, ∂Q    , and metric linearization, ∂C  , used for the flow and mesh adjoint systems (Equations 4 and 5),
respectively, are updated to include the modifications to the SA turbulence model and the transition model equations,
enabling free-transition gradient-based aerodynamic shape optimization. Because these partial derivatives are formed
analytically it is important to verify their implementation. The analytical flow Jacobian and metric linearization are each
verified through comparing with the complex-step method, while the analytical gradient is verified using a second-order
finite-difference approximation. It is important to note that the axial and FFD geometric design variables are formed
using the B-spline control points, C, [55] and that the residuals depend on the B-spline control points through the
grid metrics and the off-wall spacing, the latter of which is used in the turbulence and transition model source terms.
Furthermore, partial derivatives involving the angle of attack and sideslip design variables are determined using the
                                                     ∂R                                ∂R
complex-step method. Therefore, the analytical ∂X        is verified concurrently with ∂C .
     It is important that the verification tests include all relevant flow features. To ensure this, an infinite swept wing
geometry consisting of a blunt trailing-edge RAE2822 airfoil extruded with a 25 degree sweep and periodic boundary
conditions is simulated using the smooth transition model with the local helicity-based crossflow correlations and
compressibility corrections, SA-sLM2015cc, presented in [52]. The infinite swept wing is simulated at a Mach number

<!-- ===== PDF page 5 ===== -->

                                                                 Fig. 1                       Infinite swept wing geometry with FFD design variables in blue.

                                                                                        Flow Jacobian                                                                                                                                                                                                     Drag
                                                                                                                                                                                            Metric Linearization                              0
                                                                                                                                                                                                                                                                                                          Lift
              1                                                                                                                       1
        10                                                                                                                      10
                                                                                                                                                                                                                                         10-1
        10-1                                                                                                                    10-1
              -3                                                                                                                      -3
                                                                                                                                                                                                                                         10-2
        10                                                                                                                      10
                                                                                                                                                                                                                                         10-3
        10-5                                                                                                                    10-5
Error

                                                                                                                        Error

                                                                                                                                                                                                                                 Error
                                                                                                                                                                                                                                              -4
        10-7                                                                                                                    10-7
                                                                                                                                                                                                                                              -5
        10-9                                                                                                                    10-9                                                                                                     10

             -11
                                                                                                                                     -11                                                                                                 10-6

             -13
                                                                                                                                     -13
                                                                                                                                                                                                                                         10-7
             -15                                                                                                                     -15
        10         -2        -4        -6        -8        -10        -12        -14        -16        -18        -20
                                                                                                                                10              -6        -8        -10        -12        -14        -16        -18        -20
                                                                                                                                                                                                                                         10-8 -2             -4        -6            -8        -10        -12
              10        10        10        10        10         10         10         10         10         10                            10        10        10         10         10         10         10         10                   10           10        10            10        10         10
                                                             ε                                                                                                                  ε                                                                                           ε
                                            ∂R
                         (a) Flow Jacobian: ∂Q                                                                                                                       ∂R
                                                                                                                                           (b) Metric Linearization: ∂C                                                                            (c) Directional Derivative: Dz J

Fig. 2 Verification for the analytical flow Jacobian, metric linearization, and directional derivative including
the SA-sLM2015cc transition model equations using both complex-step and finite-difference approximations.

of 0.785, Reynolds number of 20.3 × 106 , lift coefficient of 0.56, turbulence intensity of 0.07%, and with an assumed
surface roughness of 1.8µm. The infinite swept wing O-grid consists of 561, 121, and 11 nodes in the streamwise,
off-wall, and spanwise directions, respectively, with average and maximum y+ values of 0.25 and 0.50. The geometry is
parameterized using B-spline surface patches that are embedded in an FFD volume with 6 FFD design variables on each
of the upper and lower surfaces, as presented in Figure 1. The design variables are constrained to be constant along
the span and the leading- and trailing-edge design variables are constrained to move symmetrically, resulting in 10
effective geometric design variables plus angle of attack. Although the current test does not involve axial control system
design variables, the two-level axial and FFD control system was not modified by the inclusion of the transition model
equations. The changes to the partial derivatives that were made to enable free-transition optimization can be verified by
considering only the FFD geometric design variables and angle of attack.
    The verification cases are performed using a fully-converged SA-sLM2015cc simulation on the infinite swept wing
to ensure that the turbulence and transition model source terms are active. The results for the verification study are
presented in Figure 2. Machine-zero agreement is demonstrated for both the analytical flow Jacobian and metric
linearization with approximations formed using the complex-step method. A directional derivative is used to verify the
analytical gradient components simultaneously [67, 68], defined as,
                                                                                                                                                                                     ∂J
                                                                                                                                                               Dz J =                   z,                                                                                                                  (7)
                                                                                                                                                                                     ∂X
and the second-order finite-difference approximation is,
                                                                                                                                            J (X + z) − J (X − z)
                                                                                                                        Dz J =                                      + O( 2 ),                                                                                                                              (8)
                                                                                                                                                      2

<!-- ===== PDF page 6 ===== -->

Table 1 Design conditions for the two-dimensional lift-constrained drag minimizations [32]. For each case
turbulence intensity, Tu, is specified as 0.07%.

                  design condition            lift coefficient       Mach number     Reynolds number (MAC)
                   Cessna 172R                      0.30                 0.19               5.6 × 106
              De Havilland Dash8-Q400               0.42                 0.60              15.7 × 106
                  Boeing 737-800                    0.50              0.71 (corr.)         20.3 × 106

where  is the perturbation parameter. The difference between the values produced by the analytical directional derivative
and the second-order finite-difference approximation is presented in Figure 2c. Good agreement is demonstrated
to approximately 6 and 7 orders of magnitude with the objective function set to lift and drag, respectively. Figure
2c demonstrates increased round-off error with large step sizes, which is a consequence of using a finite-difference
approximation.

                                                       III. Results
    In this section, the free-transition optimization framework presented in Section II is applied to single-point lift-
constrained drag minimizations of airfoils and an infinite swept wing geometry. Airfoil optimizations are performed
using an initial geometry based on a blunt trailing-edge RAE2822 airfoil. For the infinite swept wing optimization, the
airfoil is extruded one chord length with 11 nodes in the spanwise direction and with periodic boundary conditions
applied at the wing root and tip. The airfoil and wing surfaces are parameterized using B-spline surface patches, which
are controlled using FFD design variables. For the results presented, 6 streamwise FFD design variables are applied
to each of the upper and lower surfaces, with the infinite swept wing root and tip design variables constrained to be
equal and the leading- and trailing-edge design variables constrained to move symmetrically, for a total of 10 effective
geometric design variables plus angle of attack.

A. Airfoil Optimization
    For the two-dimensional airfoil optimizations, three design conditions are considered that are representative of
aircraft classes ranging from light aircraft up to single-aisle aircraft based on the work of Rashad and Zingg [32], as
presented in Table 1. The Cessna 172R and De Havilland Dash8-Q400 aircraft have nominally zero sweep; however, the
Boeing 737-800 design Mach number is corrected from a cruise Mach number of 0.785 to an effective Mach number of
0.71 based on a wing sweep of 25 degrees. Although the two-dimensional Boeing 737-800 design condition does not
include the effects of crossflow instabilities, and therefore is not representative of the disturbance environment for a
single-aisle aircraft wing, drag minimizations at these design conditions are valuable for investigating the streamwise
grid resolution requirements for NLF optimizations at large Mach and Reynolds numbers, as well as for comparing the
optimized designs with designs previously developed by Rashad and Zingg [32]. A drag minimization of an infinite
swept wing that includes the effects of both Tollmien-Schlichting and stationary crossflow instabilities is investigated in
Section III.B.

1. Cessna 172R Skyhawk
    Preliminary optimizations revealed that the capability of the optimization algorithm to delay the onset of boundary-
layer transition when applied to free-transition design problems was sensitive to the streamwise resolution of the
grid. To investigate this further, three structured multi-block O-type topology grids with varying levels of streamwise
grid resolution were considered for the Cessna 172R drag minimizations, with the grid details presented in Table
2. Lift-constrained drag minimizations are performed at each grid level with the cross-sectional area of the airfoil
constrained so that it cannot decrease and with loose thickness-to-chord ratio constraints enforced at each FFD design
variable pair. The optimization problem, which will be referred to as the baseline optimization problem, is given by,

                                                 min Cd (Q, X),                                                        (9)
                                                   X
                                                  s.t. Cl = Cl∗,                                                      (10)
                                                         A ≥ Ainit,                                                   (11)

<!-- ===== PDF page 7 ===== -->

Table 2 Structured multi-block O-grid dimensions for the two-dimensional optimizations at the Cessna 172R
design conditions (see Table 1).

                                                         grid level             chord x off-wall nodes                                                                   avg/max ∆s × 10−6 (chord)                                                                           avg/max y+
                                                            L0                        281x121                                                                                   1.09 / 1.19                                                                                   0.17 / 0.54
                                                            L1                        561x121                                                                                        "                                                                                             "
                                                            L2                        1121x121                                                                                       "                                                                                             "

                                                  Optimality                                                                                    2
                                                                                                                                                                   Optimality                                                                                       2
                                                                                                                                                                                                                                                                                         Optimality
                          10                      Feasibility                                            .0060                             10                      Feasibility                                               .0060                             10                        Feasibility                                             .0060
                                                                                    Merit function

                                                                                                                                                                                                        Merit function

                                                                                                                                                                                                                                                                                                                            Merit function
                                                  Merit function                                                                                                   Merit function                                                                                                        Merit function
                                                                                                         .0055                                                                                                               .0055                                                                                                               .0055
Optimality, Feasibility

                                                                                                                 Optimality, Feasibility

                                                                                                                                                                                                                                     Optimality, Feasibility
                          100                                                                                                              100                                                                                                                 100
                                                                                                         .0050                                                                                                               .0050                                                                                                               .0050

                               -2                                                                        .0045                             10
                                                                                                                                                -2                                                                           .0045                             10
                                                                                                                                                                                                                                                                    -2                                                                           .0045

                                                                                                         .0040                                                                                                               .0040                                                                                                               .0040
                          10-4                                                                           .0035                             10-4                                                                              .0035                             10-4                                                                              .0035

                               -6
                                                                                                         .0030                                  -6
                                                                                                                                                                                                                             .0030                                  -6
                                                                                                                                                                                                                                                                                                                                                 .0030
                          10                                                                                                               10                                                                                                                  10
                                                                                                         .0025                                                                                                               .0025                                                                                                               .0025
                               -8                                                                                                               -8                                                                                                                  -8
                          10                                                                             .0020                             10                                                                                .0020                             10                                                                                .0020
                                    0         5                10              15                                                                    0   20       40         60          80      100                                                                     0    20       40          60     80      100       120
                                              Design iteration                                                                                                   Design iteration                                                                                                      Design iteration

                                                                        Optimized                                                                                                             Optimized                                                                                                           Optimized
                                                                        Initial                          .015                                                                                 Initial                        .015                                                                                 Initial                        .015

                      -0.6                                                                                                             -0.6                                                                                                                -0.6
                                                                                                         .012                                                                                                                .012                                                                                                                .012

                          0.0                                                                                                              0.0                                                                                                                 0.0
                                                                                                         .009                                                                                                                .009                                                                                                                .009
Cp

                                                                                                                 Cp

                                                                                                                                                                                                                                     Cp
                          0.6                                                                                                              0.6                                                                                                                 0.6
                                                                                                           Cf

                                                                                                                                                                                                                               Cf

                                                                                                                                                                                                                                                                                                                                                   Cf
                                                                                                         .006                                                                                                                .006                                                                                                                .006
                          1.2                                                                                                              1.2                                                                                                                 1.2
                                         upper                                                                                                                 upper                                                                                                                 upper
                                                                                                         .003                                                                                                                .003                                                                                                                .003
                          1.8                       upper                                                                                  1.8                           lower                                                                                 1.8                             lower
                                                     lower                                                                                                                                      lower                                                                                                               lower
                                                                lower                                                                                                               upper                                                                                                                 upper
                          2.4                                                                            .000                              2.4                                                                               .000                              2.4                                                                               .000
                                    0   0.2        0.4         0.6       0.8                         1                                               0   0.2           0.4         0.6         0.8                       1                                               0     0.2           0.4         0.6       0.8                       1
                                                         x/c                                                                                                                 x/c                                                                                                                   x/c

                                                     (a) L0                                                                                                             (b) L1                                                                                                                (c) L2

Fig. 3 Optimization convergence histories and cross-sectional profiles of the initial and optimized designs
produced by drag minimizations with the baseline optimization problem (Equations 9-12) at the Cessna 172R
conditions with varying streamwise grid resolution.

                                                                                                                                                                  t/c ≥ 0.15t/cinit,                                                                                                                                                         (12)
where the vectors of solution and design variables are represented by Q and X, respectively, and Cl∗ represents the design
lift coefficient presented in Table 1.
     The results from the baseline drag minimizations are presented in Figure 3. The results suggest that the large
gradients in the transition region can be under-resolved due to poor streamwise grid resolution. Instead of delaying
transition, the results in Figure 3a demonstrate that for the L0 grid level the optimization algorithm generates significant
adverse pressure gradients downstream of the transition locations by reducing the trailing-edge thickness, where the
flow deceleration reduces the turbulent skin friction coefficient and skin friction drag. This behaviour is also produced
in the Dash8-Q400 optimizations presented in the following section (Figure 5a). The optimization algorithm produces
more efficient designs with larger extents of laminar flow as the streamwise resolution of the grid is increased. Figure 3
illustrates that at the L1 grid level the optimization algorithm is able to delay the onset of boundary-layer transition on
both surfaces significantly. Diminishing returns are encountered by doubling the streamwise resolution going from the
L1 to the L2 grid level, which provides a modest reduction in drag but significantly increases the computational cost of
the optimizations.
     The optimized designs in Figures 3b and 3c feature concave pressure recovery regions on the airfoil upper surface

<!-- ===== PDF page 8 ===== -->

Table 3 Aerodynamic performance of the initial and optimized designs at the Cessna 172R conditions on the
L1 grid level.

                                                              Cd (cnts.)       Cl      Cm       L/D       aoa (deg.)
                           Initial                             52.04          0.30   -0.060    57.647       0.971
               Baseline Optimization Problem                   29.04          0.30   0.034    103.294       2.542
                   Additional Constraints                      28.88          0.30   -0.060   103.872       0.551
       Re-optimized w/ Baseline Optimization Problem           27.75          0.30   -0.058   108.127       0.600

and extended regions of favourable pressure gradient on the lower surface. Instead of a flat pressure plateau on the
upper surface, the optimization algorithm has produced a flow deceleration section followed by flow acceleration.
The optimization algorithm has achieved these designs by increasing the angle of attack, decreasing the trailing-edge
thickness, and moving the aerodynamic loading upstream. However, from a manufacturability standpoint the thin
trailing edge is undesirable. Furthermore, the design features a reflexed trailing edge that reduces lift in order to
satisfy the target lift coefficient, and steep adverse pressure gradients that could reduce the performance of the airfoil at
higher angles of attack. Rashad and Zingg [32] produced a more optimal design at the same design conditions using a
discrete-adjoint optimization algorithm coupled with a stability-analysis framework that features more aft-loading with
increased trailing-edge thickness. To investigate whether the optimization algorithm can recover a geometry similar to
that of Rashad and Zingg [32] and to address the undesirable characteristics mentioned above, the optimization problem
was modified to include a pitching-moment constraint and more conservative thickness-to-chord ratio constraints. The
optimization problem with the additional constraints is presented below,

                                                 min Cd (Q, X),                                                         (13)
                                                   X
                                                   s.t. Cl = Cl∗,                                                       (14)
                                                        Cm = Cm,init,                                                   (15)
                                                        A ≥ Ainit,                                                      (16)
                                                        t/c ≥ 0.85t/cinit .                                             (17)

     Based on the results from the streamwise grid-resolution study presented in Figure 3, the drag minimization with
the optimization problem with additional constraints is performed on the L1 grid level. The results from the drag
minimizations with the baseline and more constrained optimization problems are presented in Figure 4. The optimization
with the more constrained optimization problem (Equations 13-17) produces a design with extended regions of near-zero
pressure gradient on both the upper and lower surfaces of the airfoil that delay boundary-layer transition and increase the
laminar extent of the boundary layer. The optimization algorithm produces this design by maintaining more aft-loading
relative to the baseline optimized design, with the new design also featuring increased trailing-edge thickness. The
resulting design closely resembles the design produced by the lift-constrained drag minimization performed by Rashad
and Zingg [32]. Moreover, the optimization history demonstrates that the design produced by the optimization with
additional constraints produces lower drag than the design produced by the baseline optimization problem, despite
lying in the design space of the latter. The results suggest that the Cessna 172R design space is multi-modal, with the
optimization algorithm producing (at least) two distinct local minima. Details of the aerodynamic performance of the
initial blunt trailing-edge RAE2822 airfoil and the two optimized designs are presented in Table 3.
     To confirm that the optimized design produced using the additional constraints is a local minimum in the baseline
optimization problem design space, the more constrained optimized design was used as an initial geometry with the
relaxed constraints of the baseline optimization problem. The results in Figure 4c and Table 3 demonstrate that the
optimization algorithm is able to reduce drag by an additional count by decreasing the trailing-edge thickness to create
stronger adverse pressure gradients near the trailing edge, which reduce the skin friction in the turbulent boundary
layer. However, the new design maintains a similar profile to the initial design produced with the additional constraints.
The results in Figure 4c and Table 3 confirm that the Cessna 172R design space has at least two distinct local minima,
one with a higher design angle of attack and extended regions of laminar flow on the lower surface of the airfoil and
one with more aft-loading and more balanced extents of laminar flow on the upper and lower surfaces, with the latter
outperforming the former. The results demonstrate the significant risk multimodality presents at these design conditions,

<!-- ===== PDF page 9 ===== -->

                                                  Optimality                                                                                       2
                                                                                                                                                                          Optimality                                                                                       2
                                                                                                                                                                                                                                                                                             Optimality
                          10                      Feasibility                                               .0060                             10                          Feasibility                                               .0060                             10                     Feasibility                                            .0060

                                                                                       Merit function

                                                                                                                                                                                                               Merit function

                                                                                                                                                                                                                                                                                                                               Merit function
                                                  Merit function                                                                                                          Merit function                                                                                                     Merit function
                                                                                                            .0055                                                                                                                   .0055                                                                                                           .0055
Optimality, Feasibility

                                                                                                                    Optimality, Feasibility

                                                                                                                                                                                                                                            Optimality, Feasibility
                          100                                                                                                                 100                                                                                                                     100
                                                                                                            .0050                                                                                                                   .0050                                                                                                           .0050

                               -2                                                                           .0045                             10
                                                                                                                                                   -2                                                                               .0045                             10
                                                                                                                                                                                                                                                                           -2                                                                       .0045

                                                                                                            .0040                                                                                                                   .0040                                                                                                           .0040
                          10-4                                                                              .0035                             10-4                                                                                  .0035                             10-4                                                                          .0035

                               -6
                                                                                                            .0030                                  -6
                                                                                                                                                                                                                                    .0030                                  -6
                                                                                                                                                                                                                                                                                                                                                    .0030
                          10                                                                                                                  10                                                                                                                      10
                                                                                                            .0025                                                                                                                   .0025                                                                                                           .0025
                               -8                                                                                                                  -8                                                                                                                      -8
                          10                                                                                .0020                             10                                                                                    .0020                             10                                                                            .0020
                                    0   20       40         60          80      100                                                                     0            20             40                 60                                                                       0             10                     20
                                                Design iteration                                                                                                       Design iteration                                                                                                     Design iteration

                                                                             Optimized                                                                                                             Optimized                                                                                                        Optimized
                                                                             Initial                        .015                                                                                   Initial                          .015                                                                            Initial                         .015

                      -0.6                                                                                                                -0.6                                                                                                                    -0.6
                                                                                                            .012                                                                                                                    .012                                                                                                            .012

                          0.0                                                                                                                 0.0                                                                                                                     0.0
                                                                                                            .009                                                                                                                    .009                                                                                                            .009
Cp

                                                                                                                    Cp

                                                                                                                                                                                                                                            Cp
                          0.6                                                                                                                 0.6                                                                                                                     0.6
                                                                                                              Cf

                                                                                                                                                                                                                                      Cf

                                                                                                                                                                                                                                                                                                                                                      Cf
                                                                                                            .006                                                                                                                    .006                                                                                                            .006
                          1.2                                                                                                                 1.2                                                                                                                     1.2
                                              upper                                                                                                                  upper
                                                                                                            .003                                                                                                                    .003                                                                                                            .003
                          1.8                           lower                                                                                 1.8                              lower                                                                                  1.8                                        upper
                                                                               lower                                                                                                           upper                                                                                                             upper lower
                                                                   upper                                                                                                                               lower                                                                                                           lower
                          2.4                                                                               .000                              2.4                                                                                   .000                              2.4                                                                           .000
                                    0   0.2           0.4         0.6         0.8                       1                                               0      0.2           0.4         0.6           0.8                      1                                               0     0.2      0.4         0.6         0.8                      1
                                                            x/c                                                                                                                    x/c                                                                                                               x/c

                          (a) Baseline Optimization Problem                                                                                                 (b) Additional Constraints                                                      (c) Re-optimized w/ Baseline Optimiza-
                                         filler                                                                                                                        filler                                                               tion Problem

Fig. 4 Optimization convergence histories and cross-sectional profiles of the initial and optimized designs
produced by drag minimizations with the baseline (Equations 9-12) and more constrained optimization problem
(Equations 13-17) at the Cessna 172R conditions on the L1 grid level.

Table 4 Structured multi-block O-grid dimensions for the two-dimensional optimizations at the Dash8-Q400
design conditions (see Table 1).

                                                            grid level                 chord x off-wall nodes                                                                  avg/max ∆s × 10−6 (chord)                                                                            avg/max y+
                                                               L1                            561x121                                                                                  0.65 / 0.74                                                                                    0.26 / 0.53
                                                               L2                            1121x121                                                                                      "                                                                                              "

as the baseline optimization problem produces an inferior design, and highlights the importance of efficient global
optimization techniques, such as gradient-based multistart methods [69, 70].

2. De Havilland Dash8-Q400
    Based on the results from the Cessna 172R drag minimizations, grids with streamwise resolutions representative of
the L1 and L2 grid levels are investigated for the drag minimizations at the higher Mach and Reynolds numbers of the
Dash8-Q400 design conditions. The Dash8-Q400 design conditions are presented in Table 1, with the grid details for
the L1 and L2 grid levels presented in Table 4. Similar to the Cessna 172R case, a streamwise grid-resolution study is
performed using the L1 and L2 grids and the baseline optimization problem (Equations 9-12) in order to evaluate the
capabilities of the optimization algorithm to explore the Dash8-Q400 design space. The optimization histories and
cross-sectional profiles for the initial and optimized designs are presented in Figure 5.
    The results demonstrate that there is a significant difference in both airfoil profile shape and performance of the
designs produced by the optimization algorithm on the L1 and L2 grid levels. Where the L1 grid level provides adequate
resolution for the lower Reynolds number drag minimizations at the Cessna 172R design conditions, the optimizations at
the higher Mach and Reynolds numbers of the Dash8-Q400 design conditions require finer streamwise grid resolutions.

<!-- ===== PDF page 10 ===== -->

                                                                           Optimality                                                                                 2
                                                                                                                                                                                              Optimality
                                              10                           Feasibility                                        .0060                              10                           Feasibility                                         .0060

                                                                                                         Merit function

                                                                                                                                                                                                                             Merit function
                                                                           Merit function                                                                                                     Merit function
                                                                                                                              .0055                                                                                                               .0055

                    Optimality, Feasibility

                                                                                                                                       Optimality, Feasibility
                                              100                                                                                                                100

                                                                                                                              .0050                                                                                                               .0050
                                                   -2                                                                                                                 -2
                                              10                                                                                                                 10
                                                                                                                              .0045                                                                                                               .0045
                                              10-4                                                                                                               10-4
                                                                                                                              .0040                                                                                                               .0040

                                                   -6                                                                                                                 -6
                                              10                                                                                                                 10
                                                                                                                              .0035                                                                                                               .0035

                                              10-8                                                                            .0030                              10-8                                                                             .0030
                                                        0           5       10         15         20     25                                                                0           5       10       15            20     25
                                                                          Design iteration                                                                                                   Design iteration

                                                                                                  Optimized                                                                                                           Optimized
                                                                                                  Initial                     .015                                                                                    Initial                     .015
                                          -1.2                                                                                                               -1.2

                                          -0.6                                                                                .012                           -0.6                                                                                 .012

                                              0.0                                                                                                                0.0
                                                                                                                              .009                                                                                                                .009
                                              0.6                                                                                                                0.6
                    Cp

                                                                                                                                       Cp
                                                                                                                                Cf

                                                                                                                                                                                                                                                    Cf
                                              1.2                                                                             .006                               1.2                                                                              .006

                                              1.8                                                                                                                1.8
                                                            upper                                                                                                              upper
                                                                          upper
                                                                                                                              .003                                                                                                                .003
                                                                                                                                                                                              upper
                                              2.4                           lower                                                                                2.4                                  lower
                                                                                    lower                                                                                                                     lower
                                              3.0                                                                             .000                               3.0                                                                              .000
                                                        0           0.2       0.4           0.6    0.8                    1                                                0           0.2      0.4         0.6        0.8                    1
                                                                                    x/c                                                                                                               x/c

                                                                               (a) L1                                                                                                            (b) L2

Fig. 5 Optimization convergence histories and cross-sectional profiles of the initial and optimized designs
produced by drag minimizations with the baseline optimization problem (Equations 9-12) at the De Havilland
Dash8-Q400 conditions with varying streamwise grid resolution.

There is a decrease in the transition length predicted by the transition model at the Dash8-Q400 design conditions,
illustrated in Figure 5, relative to the transition length produced at the Cessna 172R conditions, which are presented in
Figure 3. This decrease in transition length due to the increase in Reynolds number can help to explain the increased
sensitivity of the optimizations to streamwise grid resolution, as finer streamwise grid spacing is required in order to
maintain a similar resolution in the transition region to the Cessna 172R optimizations. Instead of delaying transition
on the L1 grid, the optimization algorithm prioritizes developing strong adverse pressure gradients downstream of
the transition location, which reduce the turbulent skin friction coefficient. This is similar to the behaviour of the
drag minimization at the Cessna 172R design conditions on the L0 grid level (Figure 3a). The results suggest that the
streamwise resolution of the grid must be increased as the Reynolds number increases in order to accurately resolve the
large gradients in the transition region.
    Similar to the designs in Figure 3, the designs in Figure 5 optimized using the baseline optimization problem produce
excessively thin trailing edges. In order to improve manufacturability and to further explore the design space, drag
minimizations are performed on the L2 grid level using the optimization problem with additional constraints introduced
in Section III.A.1 (Equations 13-17). The optimized designs produced using the two optimization problems on the
L2 grid level are presented in Figure 6. The results demonstrate that the optimization with the additional constraints
produces a more aft-loaded design with increased trailing-edge thickness and a more gradual concave pressure recovery.
However, opposite of the behaviour observed at the Cessna 172R design conditions (Figure 4), the more constrained
design produces increased drag relative to the design optimized with the baseline optimization problem. Although the
two optimized designs both delay transition to approximately 45% chord on the upper surface of the airfoil, the design
optimized with the baseline optimization problem maintains a larger region of laminar flow on the lower surface. In
addition, the steep pressure recovery produced by the baseline optimization problem decelerates the flow, producing
decreased skin friction in the turbulent boundary layer, resulting in a decrease in skin friction drag relative to the more
constrained design. Details of the aerodynamic performance of the initial blunt trailing-edge RAE2822 airfoil and two
optimized designs on the L2 grid level are presented in Table 5. The results demonstrate that the optimization algorithm

<!-- ===== PDF page 11 ===== -->

                                                                           Optimality                                                                                  2
                                                                                                                                                                                                   Optimality
                                              10                           Feasibility                                         .0060                              10                               Feasibility                                        .0060

                                                                                                          Merit function

                                                                                                                                                                                                                                 Merit function
                                                                           Merit function                                                                                                          Merit function
                                                                                                                               .0055                                                                                                                  .0055

                    Optimality, Feasibility

                                                                                                                                        Optimality, Feasibility
                                              100                                                                                                                 100

                                                                                                                               .0050                                                                                                                  .0050
                                                   -2                                                                                                                  -2
                                              10                                                                                                                  10
                                                                                                                               .0045                                                                                                                  .0045
                                              10-4                                                                                                                10-4
                                                                                                                               .0040                                                                                                                  .0040

                                                   -6                                                                                                                  -6
                                              10                                                                                                                  10
                                                                                                                               .0035                                                                                                                  .0035

                                              10-8                                                                             .0030                              10-8                                                                                .0030
                                                        0           5       10       15            20     25                                                                0                 10                  20        30
                                                                          Design iteration                                                                                                    Design iteration

                                                                                                   Optimized                                                                                                            Optimized
                                                                                                   Initial                     .015                                                                                     Initial                       .015
                                          -1.2                                                                                                                -1.2

                                          -0.6                                                                                 .012                           -0.6                                                                                    .012

                                              0.0                                                                                                                 0.0
                                                                                                                               .009                                                                                                                   .009
                                              0.6                                                                                                                 0.6
                    Cp

                                                                                                                                        Cp
                                                                                                                                 Cf

                                                                                                                                                                                                                                                        Cf
                                              1.2                                                                              .006                               1.2                                                                                 .006

                                              1.8                                                                                                                 1.8
                                                            upper                                                                                                               upper
                                                                                                                               .003                                                                                                                   .003
                                                                           upper                                                                                                                    lower
                                              2.4                                  lower                                                                          2.4                                       upper
                                                                                           lower                                                                                                            lower

                                              3.0                                                                              .000                               3.0                                                                                 .000
                                                        0           0.2      0.4         0.6        0.8                    1                                                0           0.2          0.4          0.6    0.8                      1
                                                                                   x/c                                                                                                                      x/c

                                              (a) Baseline Optimization Problem                                                                                                 (b) Additional Constraints

Fig. 6 Optimization convergence histories and cross-sectional profiles of the initial and optimized designs
produced by drag minimizations with the baseline (Equations 9-12) and more constrained optimization problem
(Equations 13-17) at the De Havilland Dash8-Q400 conditions on the L2 grid level.

Table 5 Aerodynamic performance of the initial and optimized designs at the Dash8-Q400 conditions on the
L2 grid level.

                                                                                                          Cd (cnts.)                                                    Cl                  Cm                      L/D            aoa (deg.)
                            Initial                                                                        54.85                                                       0.42               -0.074                   76.438            1.329
                Baseline Optimization Problem                                                              35.85                                                       0.42               -0.017                  117.124            1.998
                    Additional Constraints                                                                 40.18                                                       0.42               -0.074                  104.521            0.722

was able to reduce drag by approximately 19 and 15 drag counts using the baseline and more constrained optimization
problems, respectively. Therefore, adding a pitching moment constraint and more conservative thickness-to-chord ratio
constraints incurs a drag penalty of 4 drag counts. Further investigation is required to determine if the design produced
using the optimization problem with additional constraints is a local minima of the baseline optimization problem
design space.

3. Boeing 737-800
     To explore the performance of the free-transition optimization framework when applied to airfoil optimizations at
design conditions similar to that of a single-aisle aircraft, lift-constrained drag minimizations were performed at the
sweep-corrected Boeing 737-800 design conditions presented in Table 1. It is important to note again that because these
optimizations are two-dimensional the effects of crossflow instabilities are not included. However, the purpose of this
design problem is to evaluate the capabilities of the optimization algorithm to reduce Tollmien-Schlichting instabilities
at the high Mach and Reynolds numbers typical of a single-aisle transport aircraft. Drag minimizations are performed on
grids with streamwise grid resolutions representative of the L1 and L2 grid levels previously investigated, with the grid

<!-- ===== PDF page 12 ===== -->

Table 6 Structured multi-block O-grid dimensions for the two-dimensional optimizations at the sweep-
corrected Boeing 737-800 design conditions (see Table 1).

                   grid level                                        chord x off-wall nodes                                              avg/max ∆s × 10−6 (chord)                                                              avg/max y+
                      L1                                                   561x121                                                              0.54 / 0.60                                                                      0.25 / 0.50
                      L2                                                   1121x121                                                                  "                                                                                "

                                                                        Optimality                                                                                         2
                                                                                                                                                                                                    Optimality
                                               10                       Feasibility                                               .0060                               10                            Feasibility                                             .0060

                                                                                                             Merit function

                                                                                                                                                                                                                                       Merit function
                                                                        Merit function                                                                                                              Merit function
                                                                                                                                  .0055                                                                                                                     .0055
                     Optimality, Feasibility

                                                                                                                                            Optimality, Feasibility
                                               100                                                                                                                    100

                                                                                                                                  .0050                                                                                                                     .0050
                                                    -2                                                                                                                     -2
                                               10                                                                                                                     10
                                                                                                                                  .0045                                                                                                                     .0045
                                               10-4                                                                                                                   10-4
                                                                                                                                  .0040                                                                                                                     .0040

                                                    -6                                                                                                                     -6
                                               10                                                                                                                     10
                                                                                                                                  .0035                                                                                                                     .0035

                                               10-8                                                                               .0030                               10-8                                                                                .0030
                                                         0     5        10       15         20          25                                                                      0   2           4        6         8       10     12                    14
                                                                      Design iteration                                                                                                          Design iteration

                                                                                                 Optimized                                                                                                                  Optimized
                                                                                                 Initial                          .015                                                                                      Initial                         .015
                                           -1.2                                                                                                                   -1.2

                                           -0.6                                                                                   .012                            -0.6                                                                                      .012

                                               0.0                                                                                                                    0.0
                                                                                                                                  .009                                                                                                                      .009
                                               0.6                                                                                                                    0.6
                     Cp

                                                                                                                                            Cp
                                                                                                                                    Cf

                                                                                                                                                                                                                                                              Cf
                                               1.2                                                                                .006                                1.2                                                                                   .006

                                               1.8                                                                                                                    1.8
                                                             upper                                                                                                                      upper
                                                                      upper                                                       .003                                                                                                                      .003
                                                                                                                                                                                                     lower
                                               2.4                            lower                                                                                   2.4                           lower
                                                                              lower                                                                                                                                upper

                                               3.0                                                                                .000                                3.0                                                                                   .000
                                                         0     0.2        0.4         0.6         0.8                         1                                                 0       0.2           0.4          0.6          0.8                     1
                                                                                x/c                                                                                                                          x/c

                                                                              (a) L1                                                                                                                   (b) L2

Fig. 7 Optimization convergence histories and cross-sectional profiles of the initial and optimized designs
produced by drag minimizations with the baseline optimization problem (Equations 9-12) at the sweep-corrected
Boeing 737-800 conditions with varying streamwise grid resolution.

details presented in Table 6. Optimizations are performed using the baseline optimization problem (Equations 9-12).
    The results, which are presented in Figure 7, demonstrate that the grid-resolution requirements identified in the
Cessna 172R and De Havilland Dash8-Q400 drag minimizations are amplified at the higher Mach and Reynolds numbers
of the Boeing 737-800 conditions. Specifically, as the Reynolds number increases the transition region decreases further.
Therefore, the drag minimization on the L1 grid (Figure 7a) fails to significantly move the transition fronts aft as the
effective grid resolution in the transition region decreases. The optimization algorithm appears to encounter difficulty
trading a decrease in viscous drag associated with delaying boundary-layer transition with the increase in wave drag
associated with the development of a favourable pressure gradient to delay transition. This is also demonstrated by the
optimization history in Figure 7a, which illustrates that convergence of the merit function has stalled. The performance
of the optimization algorithm improves significantly going from the L1 to the L2 grid level. The optimization algorithm
successfully produces a favourable pressure gradient on the upper surface of the airfoil to delay transition and in the
process produces a shock wave at approximately 50% chord. The results in Figure 7b demonstrate that provided enough
grid resolution in the transition region the optimization algorithm is able to successfully trade a decrease in viscous drag
produced by delaying boundary-layer transition with an increase in wave drag.
    The aerodynamic performance metrics for the initial and optimized design on the L2 grid level are presented in
Table 7, including drag breakdowns for the initial and optimized designs. The results demonstrate that the optimized

<!-- ===== PDF page 13 ===== -->

Table 7 Aerodynamic performance of the initial and optimized designs at the sweep-corrected Boeing 737-800
conditions on the L2 grid level.

                                Cd (cnts.)    Cd, p   Cd, f      Cl       Cm       L/D      aoa (deg.)
                   Initial       54.33        18.62   35.71     0.50    -0.086    92.018      1.291
                  Optimized      42.49        13.30   29.19     0.50    -0.088   117.519      1.133

design with a larger extent of laminar flow on the upper surface of the airfoil produces decreased skin friction drag as
well as a net decrease in pressure drag, as the reduction in the pressure component of viscous drag is larger than the
increase in wave drag.

B. Infinite Swept Wing Optimization
    An infinite swept wing optimization was performed at design conditions similar to a transonic strut-braced wing
aircraft [11]. The goal of the drag minimization at this design condition is to evaluate the ability of the optimization
algorithm to delay both Tollmien-Schlichting and stationary crossflow instabilities. Specifically, a lift-constrained drag
minimization was conducted using an initial geometry consisting of an RAE2822 airfoil extruded with a 30 degree
sweep. The infinite swept wing O-grid consists of 857, 103, and 11 nodes in the streamwise, off-wall, and spanwise
directions, respectively, with an average y+ value of 0.50. Design conditions were specified as a Mach number of 0.80,
Reynolds number of 12.12 × 106 , turbulence intensity of 0.10%, surface roughness of 1µinch, with a target cruise lift
coefficient of 0.56. An optimization problem similar to the baseline optimization problem (Equations 9-12) was used
for the infinite swept wing optimization; however, to expand the design space the cross-sectional area constraint was
removed and more conservative thickness bounds were placed on the FFD design variable pairs. The infinite swept
wing optimization problem is defined as follows,
                                                min CD (Q, X),                                                       (18)
                                                  X
                                                 s.t. CL = 0.56,                                                     (19)
                                                       t/c ≥ 0.85t/cinit .                                           (20)
Similar to the airfoil optimizations, 6 streamwise FFD design variables are applied to each of the upper and lower
surfaces, with the infinite swept wing root and tip design variables constrained to be equal and the leading- and
trailing-edge design variables constrained to move symmetrically, for a total of 10 effective geometric design variables
plus angle of attack.
     The results from the drag minimization are presented in Figure 8. The optimization history illustrated in Figure 8a
demonstrates that the optimization algorithm successfully reduces drag by approximately 12 counts. The optimized
cross-sectional profiles in Figure 8b demonstrate that the optimization algorithm produces near-zero pressure gradients
on the upper and lower surfaces of the wing, which delay transition to 60% and 40% chord, respectively. While the
initial design already featured a near-zero upper surface pressure gradient, it is important to note that the optimization
algorithm has removed the favourable pressure gradient on the lower surface of the wing to prevent the growth of
crossflow instabilities and therefore delay transition on the lower surface. The optimization algorithm produced these
near-zero pressure gradients on the upper and lower surfaces by reducing the wing thickness and moving the aerodynamic
loading aft.
     The integrated forces for the initial and optimized infinite swept wing designs are presented in Table 8. The results
demonstrate that the pressure and skin friction drag components account for approximately 60% and 40% of the total
drag reduction, respectively. The relatively larger reduction in pressure drag can be explained by Figure 8b, which
illustrates that the optimized design simultaneously reduces the strength of the shock while increasing the laminar
extent of the boundary layer, which the optimizer achieved by reducing the maximum thickness-to-chord ratio and
the cross-sectional area. Future optimizations will investigate this optimization problem with a more conservative
cross-sectional area constraint, similar to that used in the airfoil optimizations.

                                                  IV. Conclusions
   The SA-sLM2015 smooth local correlation-based transition model is fully coupled with a Newton-Krylov RANS-
based flow solver and discrete-adjoint gradient-based optimization algorithm, with the analytical partial derivatives and

<!-- ===== PDF page 14 ===== -->

                                                                     Optimality                                                                                               Optimized
                                              10                     Feasibility                               .0070                                                          Initial         .015

                                                                                              Merit function
                                                                                                                            -1.2
                                                                     Merit function
                                                                                                               .0065        -0.6

                    Optimality, Feasibility
                                              100                                                                                                                                             .012

                                                                                                               .0060         0.0
                                                   -2
                                              10                                                                                                                                              .009
                                                                                                                             0.6

                                                                                                                        Cp

                                                                                                                                                                                                Cf
                                                                                                               .0055
                                              10-4                                                                           1.2                                                              .006
                                                                                                               .0050
                                                                                                                             1.8                                      upper   upper
                                                   -6
                                              10                                                                                                lower                                         .003
                                                                                                               .0045                                    lower
                                                                                                                             2.4

                                              10-8                                                             .0040         3.0                                                              .000
                                                        0       50            100       150                                        0      0.2           0.4           0.6        0.8      1
                                                                 Design iteration                                                                               x/c

                                                            (a) Optimization history                                                   (b) Cross-sectional profiles

 Fig. 8   Infinite swept wing drag minimization at the transonic strut-braced wing aircraft design conditions.

Table 8 Aerodynamic performance of the initial and optimized infinite swept wing designs at the transonic
strut-braced wing aircraft design conditions.

                                                                 CD (cnts.)           CD, p                    CD, f           CL           CM                    L/D                 aoa (deg.)
                    Initial                                       58.10               25.97                    32.13          0.56         -0.23                 96.38                  2.10
                   Optimized                                      46.12               18.99                    27.13          0.56         -0.28                 121.42                 0.94

gradients verified using finite-difference and complex-step approximations. The free-transition optimization framework
is applied to lift-constrained drag minimizations of airfoils and an infinite swept wing geometry. Airfoil optimizations
are performed at design conditions ranging from a light aircraft to a single-aisle aircraft, with the infinite swept wing
optimized at design conditions similar to a transonic strut-braced wing aircraft.
     The results demonstrate that increasing the streamwise resolution, and therefore better resolving the large gradients
in the boundary-layer transition region, improves the capability of the optimization algorithm to delay boundary-layer
transition. Airfoil optimizations at the Cessna 172R light aircraft design conditions demonstrate that the design space is
multi-modal, with the optimization algorithm producing (at least) two distinct local minima: one with a thin, reflexed
trailing edge and steep pressure recovery regions, the other similar to the design produced by Rashad and Zingg [32] at
the same design conditions featuring more aft loading. The latter design outperforms the former, demonstrating the
importance of addressing multimodality for such design problems.
     As the transition length decreases with increasing Reynolds number, finer streamwise mesh spacings are required to
maintain sufficient mesh resolution to adequately resolve the streamwise gradients in the transition region. Despite
increasing the streamwise grid resolution, drag minimizations at the De Havilland Dash8-Q400 and sweep-corrected
Boeing 737-800 design conditions fail to delay boundary-layer transition as far aft as the stability-analysis based
optimization framework developed by Rashad and Zingg [32]. However, optimizations at the sweep-corrected Boeing
737-800 design conditions demonstrate that the optimization algorithm is able to successfully trade a decrease in viscous
drag associated with delaying boundary-layer transition with an increase in wave drag. The decrease in the pressure
component of viscous drag due to the increased laminar extent of the boundary layer is more significant than the increase
in wave drag, resulting in a net decrease in both pressure and skin friction drag. The drag minimization of the infinite
swept wing geometry at the transonic strut-braced wing aircraft design conditions demonstrates that the optimization
algorithm successfully delays both Tollmien-Schlichting and stationary crossflow instabilities.
     Future work will investigate methods for reducing the streamwise grid requirements of the free-transition gradient-
based optimization framework, and the impact of increasing the number of design variables on the ability of the
optimization algorithm to delay boundary-layer transition.

                                                                                       Acknowledgments
   This work was partially funded by the NASA Transformational Tools and Technologies (TTT) and Advanced Air
Transport Technology (AATT) projects, the Natural Sciences and Engineering Research Council (NSERC), and the

<!-- ===== PDF page 15 ===== -->

University of Toronto. Computations were performed on the Niagara supercomputer at the SciNet HPC Consortium, a
part of Compute Canada. The authors gratefully acknowledge the input and support provided by Timothy Chau, Dr.
Cetin Kiris, and the LAVA team.

                                                          References
 [1] Green, J. E., “Civil aviation and the environment: the next frontier for the aerodynamicist,” Aeronautical Journal, Vol. 110, No.
     1110, 2006, pp. 469–486. doi:https://doi.org/10.1017/S0001924000001378.
 [2] Bushnell, D. M., “Overview of aircraft drag reduction technology,” AGARD Report 786 (Special Course on Skin Friction Drag
     Reduction), 1992.
 [3] Malik, M. R., Crouch, J. D., Saric, W. S., Lin, J. C., and Whalen, E. A., “Application of Drag Reduction Tech-
     niques to Transport Aircraft,” In Encyclopedia of Aerospace Engineering (eds R. Blockley and W. Shyy), 2015.
     doi:https://doi.org/10.1002/9780470686652.eae1013.
 [4] Crouch, J. D., “Boundary-Layer Transition Prediction for Laminar Flow Control,” 45th AIAA Fluid Dynamics Conference,
     AIAA Paper 2015-2472, Dallas, TX, 2015. doi:https://doi.org/10.2514/6.2015-2472.
 [5] Campbell, R. L., and Lynde, M. N., “Natural Laminar Flow Design for Wings with Moderate Sweep,” 34th AIAA Applied
     Aerodynamics Conference, AIAA Paper 2016-4326, Washington D.C., 2016. doi:https://doi.org/10.2514/6.2016-4326.
 [6] Lynde, M. N., Campbell, R. L., Rivers, M. B., Viken, S. A., Chan, D. T., Watkins, N. A., and Goodliff, S. L., “Preliminary
     Results from an Experimental Assessment of a Natural Laminar Flow Design Method,” 2019 AIAA Aerospace Sciences Meeting,
     AIAA Paper 2019-2298, 2019. doi:https://doi.org/10.2514/6.2019-2298.
 [7] Lynde, M. N., Campbell, R. L., and Viken, S. A., “Additional Findings from the Common Research Model Natural Laminar
     Flow Wind Tunnel Test,” AIAA Aviation 2019 Forum, AIAA Paper 2019-3292, 2019. doi:https://doi.org/10.2514/6.2019-3292.
 [8] Kenway, G. K. W., and Martins, J. R. R. A., “Multipoint High-Fidelity Aerostructural Optimization of a Transport Aircraft
     Configuration,” Journal of Aircraft, Vol. 51, No. 1, 2014, pp. 144–160. doi:https://doi.org/10.2514/1.C032150.
 [9] Reist, T. A., and Zingg, D. W., “High-Fidelity Aerodynamic Shape Optimization of a Lifting-Fuselage Concept for Regional
     Aircraft,” Journal of Aircraft, Vol. 54, No. 3, 2016, pp. 1085–1097. doi:https://doi.org/10.2514/1.C033798.
[10] Reist, T. A., Zingg, D. W., Rakowitz, M., Potter, G., and Banerjee, S., “Multifidelity Optimization of Hybrid
     Wing-Body Aircraft with Stability and Control Requirements,” Journal of Aircraft, Vol. 56, No. 2, 2018, pp. 442–456.
     doi:https://doi.org/10.2514/1.C034703.
[11] Chau, T., and Zingg, D. W., “Aerodynamic Design Optimization of a Transonic Strut-Braced-Wing Regional Aircraft,” Journal
     of Aircraft, Vol. 0, No. 0, 2021, pp. 1–19. doi:https://doi.org/10.2514/1.C036389.
[12] Drela, M., “Design and optimization method for multi-element airfoils,” Aerospace Design Conference, 1993.
     doi:https://doi.org/10.2514/6.1993-969.
[13] Dodbele, S. S., “Design optimization of natural laminar flow bodies in compressible flow,” Journal of aircraft, Vol. 29, No. 3,
     1992, pp. 343–347. doi:https://doi.org/10.2514/3.46167.
[14] Green, B. E., Whitesides, J. L., Campbell, R. L., and Mineck, R. E., “Method for the Constrained Design of Natural Laminar
     Flow Airfoils,” Journal of Aircraft, Vol. 34, No. 6, 1997, pp. 706–712. doi:https://doi.org/10.2514/2.2248.
[15] Pralits, J. O., “Optimal design of natural and hybrid laminar flow control on wings,” Ph.D. thesis, Department of Mechanics,
     Royal Institute of Technology, Stockholm, Sweden, 2003.
[16] Kroo, I., and Sturdza, P., “Design-Oriented Aerodynamic Analysis for Supersonic Laminar Flow Wings,” 41st Aerospace
     Sciences Meeting and Exhibit, 2003. doi:https://doi.org/10.2514/6.2003-774.
[17] Amoignon, O., Pralits, J., Hanifi, A., Berggren, M., and Henningson, D., “Shape Optimization for Delay of Laminar-Turbulent
     Transition,” AIAA Journal, Vol. 44, No. 5, 2006, pp. 1009–1024. doi:https://doi.org/10.2514/1.12431.
[18] Zhao, K., Gao, Z., and Huang, J., “Robust Design of Natural Laminar Flow Supercritical Airfoil by Multi-Objective Evolution
     Method,” Applied Mathematics and Mechanics, Vol. 35, No. 2, 2014, pp. 191–202. doi:https://doi.org/10.1007/s10483-014-
     1783-6.
[19] Robitaille, M., Mosahebi, A., and Laurendeau, E., “Design of adaptive transonic laminar airfoils using the γ- Re  ˜ θt transition
     model,” Aerospace Science and Technology, Vol. 46, 2015, pp. 60–71. doi:https://doi.org/10.1016/j.ast.2015.06.027.
[20] Zhang, Y., Fang, X., Chen, H., Fu, S., Duan, Z., and Zhang, Y., “Supercritical natural laminar flow airfoil
     optimization for regional aircraft wing design,” Aerospace Science and Technology, Vol. 43, 2015, pp. 152–164.
     doi:https://doi.org/10.1016/j.ast.2015.02.024.
[21] Han, Z.-H., Chen, J., Zhang, K.-S., Xu, Z.-M., Zhu, Z., and Song, W.-P., “Aerodynamic Shape Optimization of
     Natural-Laminar-Flow Wing Using Surrogate-Based Approach,” AIAA Journal, Vol. 56, No. 7, 2018, pp. 2579–2593.
     doi:https://doi.org/10.2514/1.J056661.
[22] Tang, Z., Chen, Y., and Zhang, L., “Natural laminar flow shape optimization in transonic regime with competitive Nash game
     strategy,” Applied Mathematical Modelling, Vol. 48, 2017, pp. 534–547. doi:https://doi.org/10.1016/j.apm.2017.04.012.
[23] Hollom, J., and Qin, N., “Uncertainty analysis and robust shape optimisation for laminar flow aerofoils,” The Aeronautical
     Journal, Vol. 125, No. 1284, 2021, pp. 365–388. doi:https://doi.org/10.1017/aer.2020.63.
[24] Sabater, C., Bekemeyer, P., and Görtz, S., “Robust Design of Transonic Natural Laminar Flow Wings under Environmental and

<!-- ===== PDF page 16 ===== -->

     Operational Uncertainties,” AIAA Scitech 2021 Forum, 2021. doi:10.2514/6.2021-0071.
[25] Driver, J., and Zingg, D. W., “Numerical Aerodynamic Optimization Incorporating Laminar-Turbulent Transition Prediction,”
     AIAA Journal, Vol. 45, No. 8, 2007, pp. 1810–1818. doi:https://doi.org/10.2514/1.23569.
[26] Drela, M., and Giles, M. B., “Viscous-inviscid analysis of transonic and low Reynolds number airfoils,” AIAA Journal, Vol. 25,
     No. 10, 1987, pp. 1347–1355. doi:https://doi.org/10.2514/3.9789.
[27] Spalart, P. R., and Allmaras, S. R., “A One-Equation Turbulence Model for Aerodynamic Flows,” 30th AIAA Aerospace Sciences
     Meeting and Exhibit, AIAA Paper 092-0439, Reno, Nevada, United States 1992. doi:https://doi.org/10.2514/6.1992-439.
[28] Liebeck, R. H., “A Class of Airfoils Designed for High Lift in Incompressible Flow,” Journal of Aircraft, Vol. 10, No. 10, 1973,
     pp. 610–617. doi:https://doi.org/10.2514/3.60268.
[29] Zingg, D. W., “An Approach to the Design of Airfoils with High Lift to Drag Ratios,” Master’s thesis, University of Toronto,
     Toronto, Ontario, Canada, 1981.
[30] Lee, J., and Jameson, A., “Natural-Laminar-Flow Airfoil and Wing Design by Adjoint Method and Automatic Transition
     Prediction,” 47th AIAA Aerospace Sciences Meeting and Exhibit, AIAA Paper 2009-897, Orlando, Florida, January, 2009.
     doi:https://doi.org/10.2514/6.2009-897.
[31] Baldwin, B., and Lomax, H., “Thin-layer approximation and algebraic model for separated turbulent flows,” 16th Aerospace
     Sciences Meeting, 1978. doi:https://doi.org/10.2514/6.1978-257.
[32] Rashad, R., and Zingg, D. W., “Aerodynamic shape optimization for natural laminar flow using a discrete-adjoint approach,”
     AIAA Journal, Vol. 54, No. 11, 2016, pp. 3321–3337. doi:https://doi.org/10.2514/1.J054940.
[33] Rashad, R., “High-Fidelity Aerodynamic Shape Optimization for Natural Laminar Flow,” Ph.D. thesis, University of Toronto,
     Toronto, Ontario, Canada, 2016.
[34] Shi, Y., Mader, C. A., He, S., Halila, G. L. O., and Martins, J. R. R. A., “Natural Laminar-Flow Airfoil Optimization Design Using
     a Discrete Adjoint Approach,” AIAA Journal, Vol. 58, No. 11, 2020, pp. 4702–4722. doi:https://doi.org/10.2514/1.J058944.
[35] Fujino, M., Yoshizaki, Y., and Kawamura, Y., “Natural-Laminar-Flow Airfoil Development for a Lightweight Business Jet,”
     Journal of Aircraft, Vol. 40, No. 4, 2003, pp. 609–615. doi:https://doi.org/10.2514/2.3145.
[36] Shi, Y., Mader, C. A., and Martins, J. R. R. A., “Natural laminar flow wing optimization using a discrete adjoint approach,”
     Structural and Multidisciplinary Optimization, Vol. 64, No. 2, 2021, pp. 541–562. doi:https://doi.org/10.1007/s00158-021-
     02936-w.
[37] Zhu, M., and Qin, N., “Balancing Laminar Extension and Wave Drag for Transonic Swept Wings,” AIAA Journal, Vol. 59,
     No. 5, 2021, pp. 1660–1672. doi:https://doi.org/10.2514/1.J059914.
[38] Poll, D., “Transition in the Infinite Swept Attachment Line Boundary Layer,” Aeronautical Quarterly, Vol. 30, No. 4, 1979, pp.
     607–629. doi:https://doi.org/10.1017/S0001925900008763.
[39] Fan, Z., Yu, X., Qin, N., and Zhu, M., “Benefit Assessment of Low-Sweep Transonic Natural Laminar Flow Wing for
     Commercial Aircraft,” Journal of Aircraft, Vol. 0, No. 0, 0, pp. 1–8. doi:https://doi.org/10.2514/1.C036138.
[40] Krumbein, A., Krimmelbein, N., and Schrauf, G., “Automatic transition prediction in hybrid flow solver, part 1: Methodology
     and sensitivities,” Journal of Aircraft, Vol. 46, No. 4, 2009, pp. 1176–1190. doi:https://doi.org/10.2514/1.39736.
[41] Krumbein, A., Krimmelbein, N., and Schrauf, G., “Automatic transition prediction in hybrid flow solver, part 2: Practical
     application,” Journal of Aircraft, Vol. 46, No. 4, 2009, pp. 1191–1199. doi:https://doi.org/10.2514/1.39738.
[42] Langtry, R. B., and Menter, F. R., “Correlation-based transition modeling for unstructured parallelized computational fluid
     dynamics codes,” AIAA Journal, Vol. 47, No. 12, 2009, pp. 2894–2906. doi:https://doi.org/10.2514/1.42362.
[43] Khayatzadeh, P., and Nadarajah, S. K., “Aerodynamic Shape Optimization via Discrete Viscous Adjoint Equations for the k-ω
     SST Turbulence and γ- Re  ˜ θt Transition Models,” 49th AIAA Aerospace Sciences Meeting including the New Horizons Forum
     and Aerospace Exposition, AIAA Paper 2011-1247, Orlando, Florida, January, 2011. doi:https://doi.org/10.2514/6.2011-1247.
[44] Khayatzadeh, P., and Nadarajah, S. K., “Aerodynamic Shape Optimization of Natural Laminar Flow (NLF) Airfoils,” 50th AIAA
     Aerospace Sciences Meeting including the New Horizons Forum and Aerospace Exposition, AIAA Paper 2012-0061, Nashville,
     Tennessee, January, 2012. doi:https://doi.org/10.2514/6.2012-61.
[45] Halila, G. L. O., Martins, J. R. R. A., and Fidkowski, K. J., “Adjoint-Based Aerodynamic Shape Optimiza-
     tion Including Transition to Turbulence Effects,” Aerospace Science and Technology, Vol. 107, 2020, pp. 1–15.
     doi:https://doi.org/10.1016/j.ast.2020.106243.
[46] Coder, J. G., and Maughmer, M. D., “Computational fluid dynamics compatible transition modeling using an amplification
     factor transport equation,” AIAA Journal, Vol. 52, No. 11, 2014, pp. 2506–2512. doi:https://doi.org/10.2514/1.J052905.
[47] Yang, Z., and Mavriplis, D. J., “Discrete Adjoint Formulation for Turbulent Flow Problems with Transition Modelling on
     Unstructured Meshes,” AIAA Scitech 2019 Forum, 2019. doi:10.2514/6.2019-0294.
[48] Mavriplis, D. J., Yang, Z., and Anderson, E. M., “Adjoint Based Optimization of a Slotted Natural Laminar Flow Wing for
     Ultra Efficient Flight,” AIAA Scitech 2020 Forum, 2020. doi:https://doi.org/10.2514/6.2020-1292.
[49] Menter, F. R., Smirnov, P. E., Liu, T., and Avancha, R., “A One-Equation Local Correlation-Based Transition Model,” Flow,
     Turbulence and Combustion, Vol. 95, 2015, pp. 583–619. doi:https://doi.org/10.1007/s10494-015-9622-4.
[50] Kaya, H., and Tuncer, I. H., “Discrete Adjoint-Based Aerodynamic Shape Optimization Framework for Natural Laminar Flows,”
     AIAA Journal, 2021. doi:https://doi.org/10.2514/1.J059923.
[51] Piotrowski, M. G. H., and Zingg, D. W., “Smooth Local Correlation-based Transition Model for the Spalart-Allmaras Turbulence

<!-- ===== PDF page 17 ===== -->

     Model,” AIAA Journal, Vol. 59, No. 2, 2021, pp. 474–492. doi:https://doi.org/10.2514/1.J059784.
[52] Piotrowski, M. G. H., and Zingg, D. W., “Compressibility Corrections to Extend a Smooth Local Correlation-Based Transition
     Model to Transonic Flows,” AIAA Journal, Submitted July, 2021.
[53] Piotrowski, M. G. H., and Zingg, D. W., “Numerical Behaviour of a Smooth Local Correlation-based Transition Model in a
     Newton-Krylov Flow Solver,” AIAA Scitech 2022 Forum, 2022.
[54] Hicken, J. E., and Zingg, D. W., “Aerodynamic Optimization Algorithm with Integrated Geometry Parameterization and Mesh
     Movement,” AIAA Journal, Vol. 48, No. 2, 2010, pp. 400–413. doi:https://doi.org/10.2514/1.44033.
[55] Gagnon, H., and Zingg, D. W., “Two-Level Free-Form and Axial Deformation for Exploratory Aerodynamic Shape Optimization,”
     AIAA Journal, Vol. 53, No. 7, 2015, pp. 2015–2026. doi:10.2514/1.J053575.
[56] Lee, C., Koo, D., and Zingg, D. W., “Comparison of B-spline Surface and Free-form Deformation Geometry Control Methods
     for Aerodynamic Shape Optimization,” AIAA Journal, Vol. 55, No. 1, 2017, pp. 228–240. doi:https://doi.org/10.2514/1.J055102.
[57] Hicken, J. E., and Zingg, D. W., “Parallel Newton-Krylov solver for the Euler equations discretized using simultaneous
     approximation terms,” AIAA Journal, Vol. 46, No. 11, 2008, pp. 2773–2786. doi:https://doi.org/10.2514/1.34810.
[58] Osusky, M., and Zingg, D. W., “Parallel Newton-Krylov-Schur Flow Solver for the Navier-Stokes Equations,” AIAA Journal,
     Vol. 51, No. 12, 2013, pp. 2833–2851. doi:https://doi.org/10.2514/1.J052487.
[59] Osusky, L., Buckley, H. P., Reist, T. A., and Zingg, D. W., “Drag Minimization Based on the Navier-Stokes Equations Using a
     Newton-Krylov Approach,” AIAA Journal, Vol. 53, No. 6, 2015, pp. 1555–1577. doi:https://doi.org/10.2514/1.J053457.
[60] Reist, T. A., Koo, D., Zingg, D. W., Bochud, P., Castonguay, P., and Leblond, D., “Cross Validation of Aerodynamic Shape
     Optimization Methodologies for Aircraft Wing-Body Optimization,” AIAA Journal, Vol. 58, No. 6, 2020, pp. 2581–2595.
     doi:https://doi.org/10.2514/1.J059091.
[61] Osusky, M., Boom, P. D., and Zingg, D. W., “Results from the Fifth AIAA Drag Prediction Workshop obtained with a
     parallel Newton-Krylov-Schur flow solver discretized using summation-by-parts operators,” 31st AIAA Applied Aerodynamics
     Conference, AIAA 2013-2511, 2013. doi:10.2514/6.2013-2511.
[62] Swanson, R. C., and Turkel, E., “On central-difference and upwind schemes,” Journal of Computational Physics, Vol. 101,
     1992, pp. 292–306. doi:https://doi.org/10.1016/0021-9991(92)90007-L.
[63] Medida, S., and Baeder, J. D., “Application of the Correlation-based γ- Re       ˜ θt Transition Model to the Spalart-Allmaras
     Turbulence Model,” 20th AIAA Computational Fluid Dynamics Conference, AIAA Paper 2011-3979, Honolulu, HI, 2011.
     doi:https://doi.org/10.2514/6.2011-3979.
[64] Langtry, R. B., Sengupta, K., Yeh, D. T., and Dorgan, A. J., “Extending the γ- Re    ˜ θt Correlation based Transition Model for
     Crossflow Effects,” 45th AIAA Fluid Dynamics Conference, AIAA Paper 2015-2474, 2015. doi:https://doi.org/10.2514/6.2015-
     2474.
[65] Hicken, J. E., and Zingg, D. W., “A Simplified and flexible variant of GCROT for solving nonsymmetric linear systems,” SIAM
     Journal on Scientific Computing, Vol. 32, No. 3, 2010, pp. 1672–1694. doi:https://doi.org/10.1137/090754674.
[66] Gill, P. E., Murray, W., and Saunders, M. A., “SNOPT: an SQP algorithm for large-scale Constrained Optimization,” SIAM
     Journal on Optimization, Vol. 12, No. 4, 2002, pp. 979–1006. doi:https://doi.org/10.1137/S0036144504446096.
[67] Hicken, J. E., “Efficient algorithms for future aircraft design: contributions to aerodynamic shape optimization,” Ph.D. thesis,
     Graduate Department of Aerospace Science and Engineering, University of Toronto, 2009.
[68] Osusky, L., “A Numerical Methodology for Aerodynamic Shape Optimization in Turbulent Flow Enabling Large Geometric
     Variation,” Ph.D. thesis, Graduate Department of Aerospace Science and Engineering, University of Toronto, 2014.
[69] Chernukhin, O., and Zingg, D. W., “Multimodality and Global Optimization in Aerodynamic Design,” AIAA Journal, Vol. 51,
     No. 6, 2013, pp. 1342–1354. doi:https://doi.org/10.2514/1.J051835.
[70] Streuber, G. M., and Zingg, D. W., “Evaluating the Risk of Local Optima in Aerodynamic Shape Optimization,” AIAA Journal,
     Vol. 59, No. 1, 2021, pp. 75–87. doi:https://doi.org/10.2514/1.J059826.
