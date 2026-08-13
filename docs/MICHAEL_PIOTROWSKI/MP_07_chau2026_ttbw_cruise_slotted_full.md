# Chau, Piotrowski & Duensing (2026, Journal of Aircraft) — Aerodynamic Optimization of a Cruise-Slotted Transonic Truss-Braced-Wing Aircraft Configuration

> **Full transcript** — verbatim text extracted from the source PDF with
> `pdftotext -layout`, page-delimited, nothing removed except bare page-number
> lines. Layout is preserved, so tables and multi-column material read as they
> appear in the PDF. **Equations are raw PDF extractions**: superscripts,
> subscripts and math symbols may be flattened or split across lines — for any
> equation you intend to implement, confirm against the PDF itself (path below)
> or the digest companion file.
>
> Source PDF: `/home/mdo/Desktop/ARTIGOS SMOTH/TChau_MPiotrowski_JoA2026 (1).pdf`

---

<!-- ===== PDF page 1 ===== -->

      Aerodynamic Optimization of a Cruise-Slotted Transonic
            Truss-Braced-Wing Aircraft Configuration

                                                     Timothy Chau∗
                              Analytical Mechanics Associates, Moffett Field, CA, 94035, USA

                                      Michael G. Piotrowski† and Jared C. Duensing‡
                                 NASA Ames Research Center, Moffett Field, CA, 94035, USA

               The transonic truss-braced-wing configuration (TTBW) and cruise-slotted wing are two

          advanced technologies with the potential to significantly improve the aerodynamic efficiency

          of next-generation transport aircraft. The TTBW achieves reduced induced drag through

          the use of an ultra-high-aspect-ratio wing supported by a truss, while the cruise-slotted wing

          employs slotted airfoils to mitigate shock formation and boundary layer separation at transonic

          flight conditions. However, it remains to be seen whether these technologies can be combined

          while retaining their individual benefits. To address this question, this paper investigates the

          aerodynamic design and performance of a cruise-slotted TTBW configuration through the

          application of aerodynamic shape optimization based on the Reynolds-averaged Navier-Stokes

          equations. Aerodynamic shape optimization is also applied to the aerodynamic design of a

          non-slotted TTBW reference aircraft. Results indicate that aerodynamic shape optimization

          is successful in achieving an aerodynamically efficient cruise-slotted TTBW design for fully

          turbulent flow, with similar features to those of the optimized non-slotted TTBW. In terms

          of aerodynamic performance, the optimized cruise-slotted TTBW experiences 3.0% more

          drag at Mach 0.80 than the non-slotted variant but benefits from a more gradual drag rise

          as the operating Mach number increases. This trend persists even when both configurations

          are optimized for higher design Mach numbers. However, the low wave drag of both sets of

          optimized designs mitigates the airfoil technology factor advantage of the cruise-slotted wing,

          resulting in lower maximum aerodynamic range efficiency compared to the non-slotted wing.

                                                          Nomenclature

𝛼        =    Angle of attack

    ∗ Senior Research Scientist/Engineer, Computational Aerosciences Branch, Building N258, AIAA Member, timothy.chau@nasa.gov
   † Previously Aerospace Engineer, Computational Aerosciences Branch, Building N258; currently Engineering Professional, Bombardier, AIAA

Member, michael.piotrowski@aero.bombardier.com
   ‡ Branch Chief, Computational Aerosciences Branch, Building N258, AIAA Member, jared.c.duensing@nasa.gov

<!-- ===== PDF page 2 ===== -->

𝐷         =   Drag

𝐶𝐷        =   Drag coefficient

𝐶𝐿        =   Lift coefficient

𝐶𝑃        =   Pressure coefficient

Cgeo      =   Geometric equality constraints

Ggeo      =   Geometric inequality constraints

𝐿         = Lift

𝐿/𝐷       = Lift-to-drag ratio

𝑀         = Mach number

x         = Design variables

𝑥, 𝑦, 𝑧   = Cartesian coordinates

𝑦+        = Nondimensional off-wall distance

                                                     Acronyms

FFD                  = Free-Form Deformation

LAVA                 = Launch, Ascent, and Vehicle Aerodynamics

LE                   = Leading Edge

MFW                  = Maximum Fuel Weight

MTOW                 = Maximum Takeoff Weight

MZFW                 = Maximum Zero Fuel Weight

OEW                  = Operating Empty Weight

RANS                 = Reynolds-Averaged Navier-Stokes

SA                   = Spalart-Allmaras

SNOPT                = Sparse Nonlinear Optimizer

SOC                  = Start of Cruise

TE                   = Trailing Edge

TTBW                 = Transonic Truss-Braced Wing

TTBW-150             = Single-Aisle Transonic Truss-Braced Wing

TTBW-CS-150          = Single-Aisle Cruise-Slotted Transonic Truss-Braced Wing

<!-- ===== PDF page 3 ===== -->

                                                   I. Introduction
    According to the Aviation Climate Action Plan released by the United States in 2021, the U.S. aviation sector is

projected to double its CO2 emissions per year by 2050 if no further action is taken due to the steady growth in the demand

for aviation [1]. NASA’s response is the Aeronautics Research Mission Directorate (ARMD) Strategic Implementation

Plan, where one of the main strategic thrusts is to achieve ultra-efficient subsonic transport by accelerating performance

improvements in drag, weight, and propulsion to reduce aviation energy consumption and environmental impact [2].

With regard to improvements in drag, two advanced air transport technologies that can potentially contribute significantly

are the transonic truss-braced wing (TTBW) aircraft configuration and the cruise-slotted wing.

    The TTBW [3] is an unconventional aircraft configuration that has the main advantage of a much larger wing span,

typically 25-50% greater than that of an equivalent conventional tube-and-wing aircraft. This significantly reduces

induced drag, improving aerodynamic efficiency and hence energy efficiency. To support the increased wing-bending

loads that come from such a large wing span, the wing is structurally braced by a main strut, with the option of

one or more jury struts. Further, the structural efficiency of the truss-braced wing topology can also be leveraged

for other design trades, for example, for reductions in wing weight, reductions in wing thickness-to-chord ratio, and

reductions in wing sweep. Reducing wing sweep can also reduce wing weight, which can be further traded for thinner

wings to counteract any penalties in wave drag. This can enhance opportunities for introducing natural-laminar-flow

wing technology by combining reduced cross-flow instabilities with the already low Reynolds numbers of the TTBW

configuration. Some examples of major efforts focused on the TTBW concept include Boeing’s Mach 0.70 SUGAR

High [4] and SUGAR Volt [5], as well as Boeing’s Mach 0.745 [6] and Mach 0.80 [7, 8] TTBW vision vehicles, all

developed for the single-aisle aircraft class. These examples feature wing spans of up to 170 ft—requiring folding wings

to accommodate the 118 ft gate limits that currently constrain the incumbent Boeing B737-8 and Airbus A320neo

aircraft. Until recently, Boeing was focused on the Mach 0.80 variant and had been engaged with NASA to develop the

Sustainable Flight Demonstrator, or X-66 [9], a research aircraft under NASA ARMD’s Sustainable Flight National

Partnership [10]. However, this effort was paused to concentrate near-term efforts on thin wing technology [11].

    The cruise-slotted wing is another advanced air transport technology, originally proposed by Whitcomb and Clark

in 1965 [12] to improve transonic aerodynamic performance. This wing configuration features a multi-element or

slotted airfoil, as shown in Figure 1, where airflow through the slot allows for an increase in aerodynamic loading over

the main element without boundary layer separation. Overall, this can enable a wing design with increased loading,

featuring a very strong pressure recovery toward the trailing edge of the aft element. The slotted airfoil can also enable

operation at higher Mach numbers given its ability to mitigate boundary layer separation encountered with shock

formation, i.e. shock-induced boundary layer separation, on the main element [12]. Conceptually, this translates to

higher airfoil technology factors, which reduce compressibility drag through the Korn equation [13]. Studies performed

by Boeing suggest that a properly designed cruise-slotted wing can allow transport aircraft to operate at increased

<!-- ===== PDF page 4 ===== -->

                                         Fig. 1   Notional slotted airfoil profile.

cruise Mach numbers for a similar maximum aerodynamic range efficiency, i.e. 𝑀 (𝐿/𝐷), where 𝑀 is the Mach

number and 𝐿/𝐷 is the lift-to-drag ratio [14]. This involved wind tunnel tests comparing a partial-span slotted wing

configuration—including the wing, fuselage, horizontal and vertical tails, nacelles, and pylons—with an equivalent

conventional wing configuration, accounting for system-level changes [14].

    The gradual pressure recoveries over the main and aft elements can also provide opportunities for natural-laminar-flow

wing technology. For example, Hiller et al. [15] applied a knowledge-based aerodynamic design framework based on

the Reynolds-averaged Navier-Stokes (RANS) equations to investigate the design of a transonic slotted airfoil and found

that their Aft-Laminar Multielement Airfoil can provide a 3% drag reduction compared to an equivalent conventional

transonic airfoil. Hiller et al. [16] then extended this to the design of a partially slotted wing based on a Mach 0.80 variant

of the NASA Common Research Model [17]. Although a similar drag count was obtained between the cruise-slotted

and conventional wing designs, a net drag benefit was achieved when considering Mach numbers exceeding the design

point, resulting from the more gradual drag rise behavior of the slotted airfoils. This was found to be the case under both

natural laminar flow and turbulent flow conditions. Other researchers have also investigated the design of slotted airfoils

for natural laminar flow, albeit for reduced transonic Mach numbers. See for example work on the Slotted Natural

Laminar Flow airfoil by Somers [18, 19], and Coder and Somers [20]. This airfoil was also implemented on the wing of

a Mach 0.70 TTBW and analyzed with free transition by Perkins et al. [21]. Although the truss was not included, results

from these simulations suggest the potential for extended regions of laminar flow over the lower and upper surfaces of

the wing.

    The design of a cruise-slotted wing is not without its challenges, however. For example, one must consider the

potential penalties in skin-friction drag that comes from the higher wetted area of the slotted airfoil. Moreover, Hiller et

al. [15, 16] showed that the high-speed flow over the aft element also causes it to exhibit relatively high skin-friction

drag coefficients that can mitigate the slotted airfoil’s overall drag benefit. Several challenges at the aircraft system level

must also be addressed, such as increased nose-down pitching moments from the increased wing aft loading, and the

structural integration of the main and aft wing elements.

    Notwithstanding these system-level challenges, a properly designed cruise-slotted wing can reduce drag relative

to an equivalent conventional wing [15, 16], while providing the benefits of operating at higher Mach numbers for a

<!-- ===== PDF page 5 ===== -->

similar aerodynamic range efficiency [14]. Therefore, the question that follows is whether the TTBW and cruise-slotted

wing technologies can be combined to improve overall vehicle aerodynamic performance. The first major challenge will

be addressing transonic interference effects arising from the narrow—and sometimes enclosed—regions formed by the

wing and strut(s), particularly when incorporating the cruise-slotted wing topology. For TTBWs, this can cause shock

formation and boundary layer separation, which can significantly diminish the overall advantages of the configuration.

    Recent advances in computational tools and capabilities for high-fidelity aerodynamic shape optimization have

proven to be invaluable in this regard. For example, Secco and Martins [22] addressed the transonic interference effects

of a Mach 0.72 single-aisle strut-braced wing through the application of aerodynamic shape optimization based on

the RANS equations. Chau and Zingg then demonstrated that these adverse effects can also be eliminated for Mach

0.78 regional jet [23, 24] and single-aisle [25] strut-braced-wing aircraft operating at the relatively high lift coefficients

required to obtain optimum aerodynamic performance. Single- and multi-point RANS-based aerodynamic shape

optimization was also applied, with results from systems analysis indicating that block fuel savings on the order of 7

to 8% can be achieved in both aircraft classes relative to conventional tube-and-wing aircraft based on the E190-E2

and A320neo. Once it has been demonstrated that the transonic interference effects of the cruise-slotted TTBW can

be mitigated, the question that remains is whether a net aircraft-level benefit can be achieved—either in the form of

reduced cruise drag, or operation at higher cruise Mach numbers without significantly compromising aerodynamic

range efficiency.

    The present paper investigates the aerodynamic design and performance of the cruise-slotted TTBW configuration

under fully turbulent flow through the application of aerodynamic shape optimization based on the RANS equations.

Specifically, a high-sweep, Mach 0.80, single-aisle TTBW configuration with a full-span slotted wing is considered to

serve as a proof of concept. The main objectives are to (i) determine whether the unique aerodynamic design challenges

of the TTBW configuration can be addressed when including the slotted airfoil, and (ii) evaluate the aerodynamic

performance of the cruise-slotted TTBW relative to a similarly-optimized non-slotted TTBW. Investigations of the

cruise-slotted TTBW’s aerodynamic performance will also consider higher Mach numbers via Mach number sensitivity

studies. Figure 2 shows concepts of the two aircraft configurations considered in this work.

    The paper is organized as follows. Section II presents the computational design and analysis tools and capabilities

used in the present paper, including details on the high-fidelity aerodynamic shape optimization framework. Section III

presents the aircraft concepts, their characteristics, and their operating conditions. Section IV then presents the

aerodynamic shape optimization problems and results. Mach number sensitivity studies are included in Section V, with

conclusions and future work presented in Section VI.

<!-- ===== PDF page 6 ===== -->

                    (a) Non-Slotted TTBW                                       (b) Cruise-Slotted TTBW

                Fig. 2   Concepts of the non-slotted and cruise-slotted TTBW aircraft configurations.

                     II. Computational Design and Analysis Tools and Capabilities

A. High-Fidelity Aerodynamic Shape Optimization Framework

   The aerodynamic design of the non-slotted and cruise-slotted TTBWs is performed through the Jetstream high-fidelity

aerodynamic shape optimization framework developed by the Computational Aerodynamics Lab at the University of

Toronto Institute for Aerospace Studies [26]. It consists of (i) an integrated mesh parameterization and deformation

method, (ii) a free-form and axial deformation geometry control system, (iii) a Newton-Krylov-Schur flow solver

for the RANS equations, (iv) the discrete-adjoint method for gradient evaluation, and (v) SNOPT for gradient-based

optimization.

   Components (i) and (ii) comprise a two-level geometry and mesh parameterization and deformation method [27]

where a baseline grid is first parameterized with B-spline curves, surfaces, and volumes. The control points of these

B-splines provide the means for mesh and hence geometry control, and form the basis of the mesh deformation

scheme [28]. Although the B-spline control points that define the aerodynamic surfaces can be used to directly control

the geometry, they are instead embedded within free-form deformation (FFD) volumes [29], which decouple the

geometry parameterization and geometry control schemes. This allows the user to control the resolution of the geometry

parameterization separately from the resolution of the geometry control system, which drives the number of optimization

design variables.

   In this method, the FFD volumes are themselves defined as B-spline volumes whose control points can be manipulated

in 𝑥, 𝑦, and 𝑧 to smoothly deform the embedded objects. In this work, the number of degrees of freedom is further

reduced by translating the FFD-volume control point coordinates to more meaningful aerodynamic design variables

as will be described in Section IV.C. As will also be described in Section IV.C, the FFD volumes are also attached

to axial curves, which aid in defining the aerodynamic design variables while also allowing for additional geometric

freedom. Once the B-spline control points defining the aerodynamic surfaces have been deformed by the geometry

control system, the perturbations are propagated across the B-spline control grid through a linear-elasticity-based model.

The deformed computational grid is then obtained by re-evaluating the parametric locations of the computational nodes

<!-- ===== PDF page 7 ===== -->

on the deformed B-spline control grid.

   Component (iii) provides the means for simulating aerodynamics and computing the aerodynamic functionals. This

parallel implicit structured multiblock Newton-Krylov-Schur algorithm [30] is used to solve the steady-state RANS

equations fully coupled with the Spalart-Allmaras (SA) [31] turbulence model. Specifically, SA-neg [32] is used with

the quadratic constitutive relation, QCR2000 [33], which has been shown to aid in more accurately detecting regions

of boundary layer separation [34], relevant to the design of wing-fuselage junctions and more importantly, wing-strut

junctions. This flow solver has been shown to accurately predict transonic aerodynamic flows when compared against

experimental data; see, for example, Osusky and Zingg [30] and Osusky et al. [35]. In the present work, boundary

layers are assumed to be fully turbulent, consistent with the focus on turbulent designs.

   Spatial discretization of the RANS equations is achieved through second-order centered summation-by-parts operators,

with first-order upwinding applied to the convective terms of the turbulence model. Simultaneous approximation

terms are also included to weakly couple block-to-block interfaces and enforce boundary conditions [36]. In addition,

a pressure-based shock sensor with scalar dissipation is included to stabilize the flow solution around shock waves.

The numerical solution is driven to steady-state via an inexact Newton method, with globalization achieved using

pseudo-transient continuation. The linear system at each Newton iteration is preconditioned with an approximate Schur

complement [37], following an incomplete LU factorization, and solved using a flexible variant of the general minimum

residual (GMRES) method [38, 39].

   Numerical optimization is provided by components (iv) and (v), which include the gradient-based optimizer

SNOPT [40] and the discrete-adjoint method [41, 42] for the computation of objective function and constraint gradients

that depend on the flow solution or the mesh deformation. These provide a computationally efficient means for solving

the aerodynamic design optimization problems considered in this work, which include on the order of hundreds of design

variables. For all other sensitivities, gradients are calculated analytically or through the complex-step method [43].

B. The Launch, Ascent, and Vehicle Aerodynamics (LAVA) Curvilinear Flow Solver

   In order to reduce computational cost, gradient-based optimizations are performed on coarse grids that are still

capable of detecting adverse flow features for the optimizer to address. Post-optimization grid-convergence studies are

then performed using the Launch, Ascent, and Vehicle Aerodynamics (LAVA) Curvilinear flow solver [44]. This flow

solver is also based on the RANS equations and is similar in accuracy and efficiency to that of Jetstream. However, the

LAVA framework enables the use of structured overset meshes, which was found to be advantageous for maintaining

the quality and resolution of the grid in and around regions of high geometric complexity such as the slot of the

cruise-slotted wing around the wing tip. These regions were found to be problematic when utilizing the automatic grid

refinement capabilities within Jetstream, leading to reduced flow solver robustness when considering very fine grids.

Fortunately, the flow features around these problem areas do not largely affect vehicle aerodynamic performance and

<!-- ===== PDF page 8 ===== -->

were found to have minimal impact on the design optimizations. As such, the coarser, but robust, structured multiblock

meshes used for optimization were still deemed sufficient.

    The version of the LAVA Curvilinear flow solver used in the present paper is a second-order finite-difference code for

the RANS equations coupled with the one-equation SA turbulence model [31]. Rotation curvature (RC) corrections [45]

and QCR2000 [33] are included, resulting in the SA-RC-QCR2000 turbulence model. Boundary layers are also assumed

to be fully turbulent. Convective fluxes are calculated using a modified Roe scheme implemented with third-order

upwind flux reconstruction and the Koren flux limiter [46, 47]. For the turbulence model’s convective terms, a first-order

upwind approach is applied. Steady-state solutions are obtained with a parallel implicit Newton-Krylov method, and

uses pseudo-transient continuation for improved robustness. Each Newton step involves solving a linear system using

the GMRES algorithm [48], preconditioned with incomplete LU factorization.

    This flow solver has also been extensively validated against experimental data. Representative examples include

transonic aerodynamic analyses of the Mach 0.745 and 0.80 Boeing TTBW configurations reported by Maldonaldo et

al. [6, 8].

                        III. Aircraft Characteristics and Operating Conditions
    The non-slotted and cruise-slotted TTBWs are developed based on a common aircraft concept, which assumes

a similar overall aircraft sizing and operating envelope between the two configurations. This aircraft concept is

shown in Figure 3 and is based on the high-sweep, Mach 0.800 variant of the Boeing SUGAR High from Bradley et

al. [4], both in terms of its overall dimensions and its top-level aircraft requirements. This concept is a 154-passenger

single-aisle transport aircraft and is developed through the application of a conceptual multidisciplinary design

optimization framework [23], which has been used to develop concepts for several conventional and unconventional

aircraft configurations based on first-order disciplinary trades (see for example References [23, 25, 49]). Table 1

summarizes the key characteristics, design requirements, and operating conditions for the conceptual model. This

includes a design payload and range of 30,800 lb and 3,500 nmi, respectively, and an economy payload and range of

30,800 lb and 900 nmi, respectively.

    For both aircraft configurations, aerodynamic shape optimization is applied at the start of cruise (SOC) conditions of

the economy range mission. This mission represents a typically-flown mission and is used to define the design operating

conditions, namely, cruise at Mach 0.800, an altitude of 43,000 ft, and a lift coefficient of 0.560. For the remainder of

the present paper, the 154-passenger non-slotted and cruise-slotted TTBW configurations are herein also referred to as

the TTBW-150 and TTBW-CS-150, respectively.

<!-- ===== PDF page 9 ===== -->

                                                               Table 1 Aircraft characteristics, design requirements,
                                                               and operating conditions for the reference TTBW
                                                               aircraft concept.

                                                                    Parameter                               Value
                                                                    Mean aerodynamic chord (ref.), ft        9.17
                                                                    Wing area (gross, ref.), ft2            1,467
                                                                    Wing span, ft                           170.0
                                                                    Wing aspect ratio                       19.70
                                                                    Wing sweep (mid chord), deg              30.0
                                                                    Wing loading, lb/ft2                    110.0
                                                                    Thrust-to-weight ratio                  0.320
                                                                    MTOW, lb                             156,010
                                                                    MZFW, lb                             132,295
                                                                    OEW, lb                                86,300
                                                                    MFW, lb                                33,430
                                                                    Maximum payload, lb                    46,000
                                                                    Design payload, lb                     30,800
     Fig. 3    Reference TTBW aircraft concept.                     Design range, nmi                       3,500
                                                                    Economy range, nmi                        900
                                                                    Mach number (SOC, economy)              0.800
                                                                    Altitude (SOC, economy), ft            43,000
                                                                    Lift coefficient (SOC, economy)         0.560

                          IV. High-Fidelity Aerodynamic Shape Optimization

A. Baseline Geometries

   Baseline geometries are developed for the TTBW-150 and TTBW-CS-150 configurations based on the aircraft

concept presented in Section III. These models, shown in Figure 4, include the wing, strut, and fuselage, which are

common between the two configurations except for the wing sections. For the baseline TTBW-150, the wing is defined by

RAE-2822 airfoils, which provide a reasonable starting point for transonic wing design. For the baseline TTBW-CS-150,

the wing is defined by slotted airfoils based on the original profiles proposed by Whitcomb and Clark [12], but modified

with an increased aft-element leading-edge (LE) radius following the designs presented in References [14, 15]. This is

done to improve aerodynamic performance robustness over a range of operating conditions. As will be described in

Section IV.C, the current work considers the design of the wing sections via a unified approach that views the wing

sections as integrated airfoils. Therefore, such detailed modifications must be introduced a priori. In addition, bracket

attachments, which connect the main and aft elements, are not modeled in the current work given the conceptual nature

of the present study.

   The strut geometries consist only of a main strut, which ends with a small vertical segment that intersects the wing

<!-- ===== PDF page 10 ===== -->

                          (a) TTBW-150                                                 (b) TTBW-CS-150

                    Fig. 4     Baseline geometries for aerodynamic shape optimization and analysis.

Table 2 Grid information for the structured multiblock grids used with the Jetstream framework for aerodynamic
shape optimization.

  Configuration         Grid Level       Number of Blocks         Number of Nodes   Avg. Off-Wall Spacing (a)   Avg. 𝑦 +(b)
  TTBW-150                   L0                   898               18.25×106              1.45×10−6              0.494
  TTBW-CS-150                L0                  1,406              22.65×106              1.42×10−6              0.511
  (a) Off-wall spacings are in units of mean aerodynamic chord.
  (b) Based on flow solutions at Mach 0.800.

normal to the surface, similar to the configurations studied by Secco and Martins [22] and Chau and Zingg [23–25].

Although such a design is less desirable in terms of structural efficiency, it readily enables geometric control of the

wing-strut junction using the free-form and axial deformation approach, making it more suitable for the exploratory

studies presented herein. Symmetric supercritical SC(2)-0012 airfoils are used to define the baseline struts.

   Lastly, the fuselage geometries are based on the Boeing SUGAR High [4], while the wing and strut fairings are

modeled after the publicly-released strut-braced-wing geometry from the PADRI workshop [50].

B. Computational Grids

   For aerodynamic shape optimization, structured multiblock grids are developed for use with the Jetstream framework.

These three-dimensional grids have O-O blocking topologies with no bi-directional boundary layer meshes, and are

developed based on the gridding guidelines of the Fourth Drag Prediction Workshop (DPW) [51]. Specifically, the

grids are developed based on the medium grid requirements and consist of 18.25 and 22.65 million grid nodes for the

TTBW-150 and TTBW-CS-150 configurations, respectively. Table 2 provides a summary of the grid information for

these coarser, optimization level grids. Note that the difference in the number of blocks is due to the higher geometric

complexity of the cruise-slotted TTBW configuration. Figure 5 shows the surface meshes and patch topologies for the

TTBW-150 and TTBW-CS-150 models.

<!-- ===== PDF page 11 ===== -->

                   (a) TTBW-150 (L0)                                                   (b) TTBW-CS-150 (L0)

Fig. 5 Surface meshes and patch topologies of the structured multiblock grids used with the Jetstream framework
for aerodynamic shape optimization.

Table 3 Grid information for the structured overset grids used with the LAVA Curvilinear flow solver for
post-optimization analysis.

                      Grid Level       Number of Nodes             Avg. Off-Wall Spacing (a)   Avg. 𝑦 +(b)
                                                            TTBW-150
                          L0               33.10×106                      2.14×10−5              0.937
                          L1               64.78×106                      1.15×10−5              0.684
                          L2              130.54×106                      5.14×10−6              0.560
                                                          TTBW-CS-150
                          L0               35.67×106                      2.08×10−5              0.938
                          L1               69.58×106                      1.11×10−5              0.676
                          L2              140.33×106                      5.01×10−6              0.568
                      (a) Off-wall spacings are in units of mean aerodynamic chord.
                      (b) Based on flow solutions at Mach 0.800.

   Post-optimization flow analysis for computing aerodynamic performance is performed with the LAVA Curvilinear

flow solver, which uses structured overset meshes in the current work. These grids are generated for the optimized

TTBW-150 and TTBW-CS-150 designs and include coarse, medium, and fine grid levels for performing grid-convergence

studies. These grid families provide the means for estimating grid-converged aerodynamic functionals, enabling more

direct comparisons between the two aircraft configurations, which differ in geometric complexity. Grid families are also

developed for the additional optimized designs presented in Section V, which are a part of the Mach number sensitivity

studies. Grid information for the grid families developed for the nominal operating conditions is included in Table 3.

The corresponding surface meshes and patch topologies for each grid family are shown in Figure 6.

<!-- ===== PDF page 12 ===== -->

                  (a) TTBW-150 (L0)                                     (d) TTBW-CS-150 (L0)

                  (b) TTBW-150 (L1)                                     (e) TTBW-CS-150 (L1)

                  (c) TTBW-150 (L2)                                     (f) TTBW-CS-150 (L2)

Fig. 6 Surface meshes and patch topologies of the structured overset grids used with the LAVA Curvilinear flow
solver for post-optimization analysis.

<!-- ===== PDF page 13 ===== -->

C. Optimization Problems

   The aerodynamic design of the TTBW-150 and TTBW-CS-150 configurations are formulated as optimization

problems of the form:

                                             minimize       𝐶𝐷

                                        with respect to     x

                                             subject to     𝐿 = 𝑊,

                                                            Cgeo = 0,    Ggeo ≤ 0

where 𝐶𝐷 is the cruise drag coefficient at the design point and x = [𝛼, xTgeo ] T is the design variable vector, which

includes angle of attack, 𝛼, and the vector of geometric design variables, xgeo . Cgeo and Ggeo are the geometric equality

and inequality constraints, respectively.

   Figure 7 shows the free-form and axial deformation geometry control systems, which define the geometric design

variables. For both aircraft configurations, the B-spline control points that parameterize the wing and strut are embedded

within FFD volumes, with the FFD volume control points arrayed in 2 × 11 cross-sections, which are used to define

twist and section shape design variables. The twist design variable is defined as a rotation of a given FFD volume

cross-section about the local origin, while the section shape design variables are vertical displacements along the

local 𝑧-axis. The local coordinate systems are defined by the axial curves positioned at the 1/4-chord of the wing and

strut. These axial curves are defined as B-spline curves whose control points can also be manipulated to drive the

attached FFD-volume cross-sections, providing more global degrees of freedom such as sweep, span, and dihedral

design variables. However, these degrees of freedom are omitted in the current work to minimize the impact of the

aerodynamic design on other disciplines.

   Note that, to reduce optimization bias toward any specific area, the FFD-volume cross-sections are distributed as

uniformly as possible along the wing and strut. However, exceptions are made at the wing root and tip, and around

the wing-strut junction, where additional degrees of freedom are necessary for addressing adverse effects in high

gradient regions of the flow, including transonic interference. In the present paper, the main and aft elements of the

TTBW-CS-150 configuration are also only embedded within one FFD volume to reduce geometric complexity. This

provides an integrated approach to design, which views the multi-element topology as a unified wing section, and is

more suitable to the implicit geometry control of FFD methods.

   In order to minimize drag under cruise conditions, a constant lift constraint is included based on the weight of the

aircraft. Another important consideration for cruise flight is the pitching moment of the aircraft, which can be trimmed

<!-- ===== PDF page 14 ===== -->

                       (a) TTBW-150                                               (b) TTBW-CS-150

Fig. 7 Axial and free-form deformation geometry control systems for aerodynamic shape optimization with
axial curves (magenta lines) and FFD volumes (blue points and surfaces).

with a zero pitching moment constraint. However, since the computational aircraft models do not include the tail systems,

such a constraint is not included. This reduces not only the computational cost of a given flow solution, but also reduces

the computational cost of the optimization problem by not requiring an additional flow-based adjoint solution. More

importantly, it decouples the present initial aerodynamic assessment of the cruise-slotted TTBW configuration from

uncertainties in the conceptual design, including center-of-gravity location and horizontal tail sizing. This allows the

study to focus on the relative aerodynamic performance of the two configurations under comparable cruise conditions.

   Geometric constraints are also included to maintain the feasibility of the aircraft concepts. This includes a minimum

wing and strut volume constraint based on the volume of the fuel tanks estimated from systems analysis, as well as

minimum (𝑡/𝑐)max constraints, which aid in maintaining a minimum structural depth, also based on systems analysis. As

shown in Figure 8, the minimum (𝑡/𝑐)max constraints are enforced at 30 wing and strut sections and are calculated at each

station based on the maximum of 22 thickness probes distributed along the chord. These constraints are formulated with

a smooth maximum based on the work of Kreisselmeier and Steinhauser [52]. To reduce geometric complexity, linear

twist constraints are also included over the strut FFD-volume cross-sections. Specifically, a linear twist interpolation

is imposed on the FFD-volume cross-sections over the transition strut, and a twist linking constraint is imposed over

the vertical strut. A linear twist interpolation constraint is also applied over the wing across the wing-strut junction to

prevent highly local design features, which were encountered during earlier aerodynamic design optimizations.

   Lastly, minimum 𝑡/𝑐 scaling and fixed LE and trailing-edge (TE) constraints are also included. The former enforces

a lower and upper bound on the local vertical distance between a given pair of FFD-volume control points. The latter

prevents shear twist and translations of the wing systems by enforcing equal and opposite deformations of the upper and

lower control points for each pair of leading- and trailing-edge control points. Tables 4 and 5 summarize the design

variables and constraints, respectively, for the TTBW-150 and TTBW-CS-150 optimizations.

<!-- ===== PDF page 15 ===== -->

                      (a) TTBW-150                                                     (b) TTBW-CS-150

 Fig. 8      Minimum (𝒕/𝒄)max constraint probes (teal points and lines) for aerodynamic shape optimization.

                          Table 4        Aerodynamic shape optimization design variables.

                              Design variable      Quantity     Lower bound        Upper bound
                              Angle of attack          1           −2.0 deg          +2.0 deg
                              Twist                   35           −10.0 deg         +10.0 deg
                                 Wing                 17           −10.0 deg         +10.0 deg
                                 Strut                18           −10.0 deg         +10.0 deg
                              Section shape           770             0.5                2.0
                                 Wing                 374             0.5                2.0
                                 Strut                396             0.5                2.0
                              Total                   806            ——                ——

                              Table 5      Aerodynamic shape optimization constraints.

Constraint            Quantity     Description
Lift                      1        Constrains the aircraft lift to equal the weight at the start of cruise (nonlinear)
Minimum volume            1        Constrains the minimum wing (and strut) OML volume based on fuel storage requirements
                                   (nonlinear)
Minimum (𝑡/𝑐)max         30        Minimum maximum thickness-to-chord ratio constraints based on structural requirements
                                   (nonlinear)
Minimum 𝑡/𝑐
                        385        Constrains the local vertical separation between each pair of FFD-volume control points
scaling
                                   to be greater or equal to 50% of its initial value (linear)
Fixed LE/TE              70        Constrains section shape design variables at the LE and TE to be equal and opposite
                                   between the lower and upper FFD-volume control points (linear)
Linear junction
                          1        Interpolates the twist design variables across the 9th and 11th wing FFD-volume cross-
wing twist
                                   sections (linear)
Linked vertical
                          1        Links the twist design variables of the vertical strut segment (linear)
strut twist
Linear transition
                          1        Interpolates the twist design variables across the transition strut segment (linear)
strut twist

<!-- ===== PDF page 16 ===== -->

D. Optimization Results

   Optimization of the TTBW-150 and TTBW-CS-150 configurations was completed in approximately 100 major

iterations, with the objective functions converged to within one-tenth of a drag count and constraint violations reduced

to 10−5 or less. This required approximately 200 function evaluations, each consisting of one RANS solution and two

adjoint solutions—one for the objective function and one for the lift constraint—at roughly half the computational cost

of a flow solution each. Figures 9 and 10 provide overviews of the initial and optimized design and flow features for the

TTBW-150 and TTBW-CS-150 configurations, respectively. From Figures 9a and 10a, it can be seen that the initial

shock waves on the wing upper surface and within the wing-strut junction of both configurations have been largely

eliminated. Although RANS-based aerodynamic shape optimization tools have previously demonstrated success in

addressing the transonic interference effects within the wing-strut junction at high transonic Mach numbers typical of

current single-aisle transport aircraft [23–25], this has not yet been shown to be possible for a cruise-slotted variant of

the TTBW. However, this should not come as a surprise since the unique features of the cruise-slotted wing topology are

positioned downstream of the location where the strut intersects with the wing, reducing its influence on the junction

flow.

   Figures 9c and 10c show the initial and optimized airfoils and pressure distributions. Note that the airfoil profiles are

not rotated with respect to the angle of attack. Following optimization, the steep pressure gradients on the initial wing

upper surface of the TTBW-150 and TTBW-CS-150 designs have been smoothed to achieve gradual pressure recoveries,

indicating that shock formation has been largely mitigated. For the optimized TTBW-CS-150, the pressure distributions

indicate higher aft loading over the main element, with the high-speed slot flow from the wing lower surface appearing

to help prevent boundary layer separation. The aft element is also highly loaded compared to a similar length of chord

of the non-slotted wing. This leads to a more uniformly distributed chordwise loading with a fairly modest pressure

recovery. However, the trade-off is that the swept wing now produces more of a nose-down pitching moment, which

results in trim drag penalties not accounted for in the present study.

   Figures 9b and 10b show the initial and optimized surface pressure contours over the upper wing and strut surfaces.

For the optimized TTBW-150, the pressure contours vary smoothly with well-aligned isobars. This is consistent with

previous single-point optimization results [23, 25], which indicate low wave drag and delayed drag divergence. On the

other hand, although the surface pressure contours of the optimized TTBW-CS-150 approach a similar pattern over

the main element, its contours are disrupted by non-uniform spanwise pressure distributions over the aft element. For

example, the aft element experiences more load near the wing root, which tapers toward the intersection of the wing and

strut. Nonetheless, such non-uniformities may be the result of a lack of geometric freedom available to the optimizer for

more carefully tailoring the shape of the main and aft elements independently. Multipoint optimization could also aid

in producing a design more similar to the optimized TTBW-150 when considering aerodynamic performance over a

range of operating conditions. However, multipoint optimization is not considered in the current work for reasons of

<!-- ===== PDF page 17 ===== -->

                         (a) Shock surfaces

                   (b) Surface pressure contours

                                                                          (c) Airfoil profiles and pressure distributions

Fig. 9 TTBW-150: Initial and optimized design and flow features. Flow features are computed on the L0 grid
level using Jetstream at 𝑴 = 0.800 and 𝑪𝑳 = 0.560.

computational cost.

   With the optimization problems formulated as lift-constrained drag minimizations without constraints on pitching

moment, the aerodynamic shape optimizations result in net spanwise lift distributions that are near elliptical. This is shown

in Figures 11a and 11b for the TTBW-150 and TTBW-CS-150 configurations, respectively. For the optimized TTBW-150,

the spanwise lift distribution over the wing features a local increase near the wing-strut junction, compensating for the

negative lift over the strut and achieving a relatively smooth elliptical distribution across the intersection. This negative

strut lift was shown previously [22–25, 53] to be necessary for mitigating shock formation and, in some cases, boundary

layer separation, which arises from the transonic channel effect. For the TTBW-CS-150 configuration, the optimizer is

also capable of achieving a net spanwise lift distribution that is smooth and close to elliptical, despite the disruptions

in the upper surface pressure contours observed in Figure 10c. The optimizer also introduces a similar positive and

<!-- ===== PDF page 18 ===== -->

                         (a) Shock surfaces

                   (b) Surface pressure contours

                                                                           (c) Airfoil profiles and pressure distributions

Fig. 10 TTBW-CS-150: Initial and optimized design and flow features. Flow features are computed on the L0
grid level using Jetstream at 𝑴 = 0.800 and 𝑪𝑳 = 0.560.

negative shift over the wing and strut, respectively, near the wing-strut junction. Interestingly, the lift compensation

occurs only on the main element of the cruise-slotted wing, while the spanwise lift distribution on the aft element

remains smooth and nearly monotonic from root to tip. Although the aft element loading is expected to be only weakly

sensitive to the effective angle of attack, which is essentially set by the main element [19], it is still notable that changes

to the main-element airfoils immediately upstream of the slot along the span have minimal impact on the aft-element

loading.

    Figures 12 and 13 show the initial and optimized wing and strut airfoils and pressure distributions near the wing-strut

junction, providing a closer look at the local flow field. Note that sections A, B, C, and D correspond to 47%, 48%,

49%, and 50% semispan, respectively. Wing and strut surface pressure contours are also shown, with inner (inboard)

and outer (outboard) views of the wing-strut junction. Here, it can be seen more clearly that, in the case of the initial

<!-- ===== PDF page 19 ===== -->

                       (a) TTBW-150                                               (b) TTBW-CS-150

Fig. 11 Optimized spanwise lift distributions computed on the L0 grid level using Jetstream at 𝑴 = 0.800 and
𝑪𝑳 = 0.560.

geometries, the wing-strut junctions behave as transonic channels, accelerating the flow and causing a shock wave to

develop. This comes in the form of sharp pressure drops on the lower wing surface of both aircraft configurations.

More gradual pressure drops are also seen on the upper strut surfaces. By introducing an inset and bump on the lower

wing surfaces, the optimizer is able to increase the pressure coefficient and mitigate this effect. This coincides with an

outwards force distribution over the strut sections. Although these design features were also observed by Chau and

Zingg [23–25], what is remarkable is that the optimizer opts to introduce a similar set of features for the cruise-slotted

wing topology. Another important observation is that the optimizer does not appear to change the shape of the aft

element much when reviewing Figures 10c and 13. However, this may be due to the lack of geometric freedom that

comes as a consequence of designing the cruise-slotted wing as an integrated, single-element wing.

    Table 6 presents the aerodynamic performance of the initial and optimized TTBW-150 and TTBW-CS-150 designs.

Based on results obtained from Jetstream on the optimization-level (L0) grids, the optimized designs reduce drag by

17.5% and 19.3%, respectively, relative to the initial designs, while maintaining the design lift coefficient. Table 6 also

includes estimates of the grid-converged aerodynamic functionals (L∞), obtained through Richardson extrapolation. As

described in Section II.B, post-optimization grid-convergence studies are performed using the LAVA Curvilinear flow

solver to reduce uncertainties arising from differences in grid resolution and topology between the computational grids

of each aircraft configuration.

    Figure 14 shows a comparison of pressure distributions between the Jetstream flow solutions computed on the

structured multiblock grids and those of the LAVA Curvilinear flow solver computed on the structured overset grids.

Good agreement is observed between the two solvers, even on the coarsest (L0) grid levels, with the exception of

small pressure disturbances on the upper wing surfaces that grow toward the root. These discrepancies arise from the

geometry extraction process, in which the optimized wing and strut designs are re-installed on the baseline fuselage

<!-- ===== PDF page 20 ===== -->

           (a) Surface pressure contours with inner (left) and outer (right) views of the wing-strut junction.

                                     (b) Airfoil profiles and pressure distributions.

Fig. 12 TTBW-150: Initial and optimized (a) surface pressure contours and (b) pressure distributions computed
on the L0 grid level using Jetstream at 𝑴 = 0.800 and 𝑪𝑳 = 0.560.

Table 6 Aerodynamic performance computed on the optimization level grids (L0) using Jetstream and
extrapolated (L∞) based on LAVA Curvilinear flow solutions.

               Parameter                Initial (L0)    Optimized (L0)       Δ (L0)     Optimized (L∞)
                                                       TTBW-150
               Angle of attack, deg         1.78              2.00           ——              ——
               Lift coefficient            0.560             0.560           +0.0%           0.560
               Drag coefficient           0.0285            0.0235          −17.5%          0.0216
               Lift-to-drag ratio           19.6              23.8          +21.4%            25.9
                                                    TTBW-CS-150
               Angle of attack, deg         1.21              2.00           ——              ——
               Lift coefficient            0.560             0.560           +0.0%           0.560
               Drag coefficient           0.0311            0.0251          −19.3%          0.0223
               Lift-to-drag ratio           18.0              22.3          +23.9%            25.2

<!-- ===== PDF page 21 ===== -->

             (a) Surface pressure contours with inner (left) and outer (right) views of the wing-strut junction.

                                       (b) Airfoil profiles and pressure distributions

Fig. 13 TTBW-CS-150: Initial and optimized (a) surface pressure contours and (b) pressure distributions
computed on the L0 grid level using Jetstream at 𝑴 = 0.800 and 𝑪𝑳 = 0.560.

and fairing geometries prior to structured overset mesh generation, and are not expected to significantly influence the

overall aerodynamic trends. Figure 15 shows the results of the grid-convergence studies, which demonstrate monotonic

convergence of the drag coefficient at constant lift for both configurations.

   Based on an assumed second-order rate of convergence, Richardson extrapolation indicates that the optimized

TTBW-CS-150 experiences slightly more cruise drag than the optimized TTBW-150—specifically, a 3.0% increase,

corresponding to a 3.0% deficit in lift-to-drag ratio at cruise. This is likely due to the increased wetted area of the

multi-element wing, as well as the presence of supervelocities in and around the slot and aft element, which increase

skin-friction drag.

   Indeed, Hiller et al. [16] showed that the sectional skin-friction drag of a partial-span cruise-slotted wing is higher

over the region of the slot than an equivalent conventional wing under the conditions of fully turbulent flow. Figure 16a

shows a similar trend between the optimized TTBW-150 and TTBW-CS-150 designs, where the sectional skin-friction

drag of the cruise-slotted wing is on average 18% greater than the non-slotted wing. For the cruise-slotted wing, the

<!-- ===== PDF page 22 ===== -->

                        (a) TTBW-150                                                 (b) TTBW-CS-150

Fig. 14 Comparison of optimized pressure distributions between Jetstream and the LAVA Curvilinear flow
solver on the L0 grid levels. Results are computed at 𝑴 = 0.800 and 𝑪𝑳 = 0.560.

Fig. 15 Grid-convergence studies performed using the LAVA Curvilinear flow solver at 𝑴 = 0.800 and
𝑪𝑳 = 0.560. Richardson extrapolation is performed assuming second-order convergence.

average sectional skin-friction drag of the main wing element alone is 91% that of the non-slotted wing, suggesting

that the majority of the drag increase comes from increased skin-friction drag coefficients rather than increased wetted

area. For reference, the average sectional skin-friction drag breakdown between the main and aft wing elements of

the cruise-slotted wing is 77% and 23%, respectively, based on the results obtained on the optimization level grids.

Figure 16b indicates that the optimized cruise-slotted wing also experiences similar sectional pressure drag compared to

the non-slotted wing. This is consistent with the results of Hiller et al. [16], who found that the partial-span cruise-slotted

wing had a similar sectional pressure drag distribution to that of the equivalently-designed conventional wing across the

entire span in turbulent flow.

<!-- ===== PDF page 23 ===== -->

               (a) Sectional skin-friction drag                                (b) Sectional pressure drag

Fig. 16 Optimized sectional wing skin-friction and pressure drag coefficients computed on the L0 grid level
using Jetstream at 𝑴 = 0.800 and 𝑪𝑳 = 0.560.

   Figure 17 compares the skin-friction drag distributions at various spanwise stations of the optimized wings. The

skin-friction drag coefficient is higher overall for the slotted airfoils, particularly in the regions surrounding the slot

and aft element. This increase is due in large part to the formation of a second thin, newly developed boundary layer

originating from the stagnation point of the aft element. In addition, the skin-friction drag coefficient increases along the

lower surface of the main element toward its trailing edge due to the influence of accelerated slot flow, which may also

be contributing to elevated skin-friction levels over the upper surface of the aft element in the mid-chord region. It is

also interesting to note that just inboard of the wing-strut junction at 45% semispan, both the optimized TTBW-150 and

TTBW-CS-150 wings experience an increase in the skin-friction drag coefficient upstream of where the shock waves

initially form. This results from the optimizer’s solution to the transonic channel effect, which involves shifting the

region of flow acceleration upstream to achieve an increase in local pressure.

   Figure 18 shows velocity magnitude contours at 15% semispan for both aircraft configurations. The contours

suggest that the increased skin-friction drag on the aft element arises in part from accelerated flow over its suction side.

Additionally, higher-speed flow is observed near the TE of the cruise-slotted wing’s main element. These flow features

are consistent with those at other spanwise stations across both wings.

   Figure 19 shows the pressure and skin-friction drag breakdown for the optimized TTBW-150 and TTBW-CS-150

designs. These values are estimated based on Richardson extrapolation. Here, it can be seen that the proportion of

skin-friction drag relative to the total drag is greater for the cruise-slotted TTBW relative to the non-slotted TTBW.

Comparing the drag components directly, the optimized TTBW-CS-150 experiences 7.6% more skin-friction drag and

2.6% less pressure drag than the optimized TTBW-150. Note, however, that the pressure drag computed from the RANS

solutions inherently includes viscous pressure drag arising from boundary-layer growth. While this contribution may be

higher for slotted airfoils, it represents an important component of the aerodynamic trade space, as the slotted airfoil

<!-- ===== PDF page 24 ===== -->

Fig. 17 Optimized wing skin-friction drag distributions computed on the L0 grid level using Jetstream at
𝑴 = 0.800 and 𝑪𝑳 = 0.560.

                      (a) TTBW-150                                             (b) TTBW-CS-150

Fig. 18 Optimized velocity magnitude contours at 15% semispan, normalized by freestream conditions.
Solutions are computed on the L0 grid level using Jetstream at 𝑴 = 0.800 and 𝑪𝑳 = 0.560.

topology enables operation at higher loading or Mach numbers through delayed boundary-layer separation. This effect

is investigated through Mach number sensitivity studies presented in the following section.

                                   V. Mach Number Sensitivity Studies
   As described in Section I, one of the main benefits of the cruise-slotted wing is its potential to achieve a similar

maximum aerodynamic range efficiency, namely, 𝑀 (𝐿/𝐷) as that of a non-slotted wing, but at higher Mach numbers [14].

<!-- ===== PDF page 25 ===== -->

                               Fig. 19    Pressure and skin-friction drag breakdowns.

In principle, this would enable an improvement to airline productivity without compromising aircraft fuel efficiency. In

this section, Mach number sensitivity studies are performed to investigate this hypothesis for the cruise-slotted variant of

the TTBW. Specifically, the optimized TTBW-150 and TTBW-CS-150 designs presented in Section IV.D are analyzed

over a range of Mach numbers, from Mach 0.740 to 0.830, trimmed with respect to angle of attack to maintain the

design lift coefficient, and compared with one another.

   Results indicate that the cruise-slotted variant of the TTBW experiences a more gradual drag rise than the non-slotted

TTBW, but the maximum 𝑀 (𝐿/𝐷) remains lower for the TTBW-CS-150 configuration due to the drag penalties

described in Section IV.D. Because the single-point optimizations at Mach 0.800 led to the elimination of most, if not

all, of the wave drag encountered by the initial TTBW-150 and TTBW-CS-150 designs, it was posited that the benefit of

the cruise-slotted wing’s increased airfoil technology factor could not manifest as an aerodynamic trade advantage over

the non-slotted wing. This is especially the case at the design point, as discussed in Section IV.D, where the increased

skin-friction drag experienced by the optimized TTBW-CS-150 led to higher overall drag compared to the optimized

TTBW-150 since both configurations experienced similar pressure drag and minimal wave drag. In practice, transonic

wing designs typically maintain some wave drag as a component of the aerodynamic trades, whether for reasons of

aerodynamic efficiency, performance robustness, or airline productivity. This is supported in part by investigations

involving multipoint optimization, where weak to moderate shocks often remain on the upper surfaces of optimized

airfoils and wings. See for example, Nemec et al. [54], Kenway and Martins [55], and Chau and Zingg [24, 25].

   This motivates the need for multipoint optimizations of the TTBW-150 and TTBW-CS-150 configurations. However,

to reduce computational cost, the present study instead repeats the single-point optimizations but at higher design Mach

numbers. Specifically, design Mach numbers of 0.815 and 0.830 are considered for the additional designs under the

assumption that the optimizer will make compromises in the form of wave drag when asked to push the envelope of the

aerodynamic designs. These additional design optimizations consider the same problem formulations as those presented

<!-- ===== PDF page 26 ===== -->

                   (a) TTBW-150 M-0.800                                       (d) TTBW-CS-150 M-0.800

                   (b) TTBW-150 M-0.815                                        (e) TTBW-CS-150 M-0.815

                   (c) TTBW-150 M-0.830                                        (f) TTBW-CS-150 M-0.830

Fig. 20 Initial and optimized design-point shock surfaces computed on the L0 grid level using Jetstream at
𝑪𝑳 = 0.560.

in Section IV.C and are converged to similar optimization tolerances. All optimized designs are also based on the same

aircraft concept described in Section III.

   Figure 20 shows the initial and optimized shock surfaces for all six of the optimized designs, while optimized

pressure distributions are included in Figure 21. From Figure 20, it can be seen that as the initial designs are operated at

successively higher Mach numbers, the shock waves grow in strength as expected. However, what is surprising is that

the optimizer is capable of nearly completely mitigating shock formation across all of the wing and strut designs, even

when the design Mach number is increased to 0.830. The wing-strut junctions also remain shock-free for all optimized

designs. This is consistent with the pressure distributions shown in Figure 21, which remain relatively smooth across all

design Mach numbers.

   Although the unexpected success of the single-point optimizations reinforce the need for investigations involving

multipoint optimization, the Mach number sensitivity studies are nonetheless performed for all six configurations. In

order to estimate the grid-converged aerodynamic functionals at each operating point, as done previously in Section IV.D,

grid-convergence studies are first performed. These are performed on the grid families presented in Section IV.B for

the Mach 0.800 designs, as well as grid families generated for the TTBW-150 and TTBW-CS-150 designs optimized

at Mach 0.815 and 0.830, which are developed with similar grid properties. From Figure 22, it can be seen that

monotonic convergence of the drag coefficient is obtained over each grid level, improving the reliability of the Richardson

extrapolations.

   Figure 23 presents the drag rise and aerodynamic range efficiency curves for the Mach 0.815 and 0.830 TTBW-150

and TTBW-CS-150 designs, as well as those of the Mach 0.800 designs previously discussed in this section. From

Figure 23a, it can be seen that for each configuration, as the design Mach number is increased, the drag rise becomes

<!-- ===== PDF page 27 ===== -->

                                                (a) TTBW-150

                                              (b) TTBW-CS-150

Fig. 21 Optimized pressure distributions for the Mach 0.800, 0.815, and 0.830 designs operating at the design
Mach numbers and 𝑪𝑳 = 0.560. Solutions are computed on the L0 grid level using Jetstream.

<!-- ===== PDF page 28 ===== -->

                                               (a) Mach 0.800 designs

                  (b) Mach 0.815 designs                                     (c) Mach 0.830 designs

Fig. 22 Grid-convergence studies performed using the LAVA Curvilinear flow solver to estimate grid-converged
aerodynamic functionals at each operating condition.

more gradual within the Mach 0.800 to 0.830 range. This comes at the cost of increased drag at Mach numbers below

0.800. However, when comparing the optimized TTBW-150 and TTBW-CS-150 configurations developed for a given

design Mach number, it can be seen that the cruise-slotted wing topology results in a more gradual drag rise than the

non-slotted variant as the operating Mach number is increased. This trend comes with diminishing returns, with the

advantage reducing most significantly when the design Mach number is increased to 0.830.

   Figure 23b shows that for a given configuration, as the design Mach number is increased, the location of maximum

𝑀 (𝐿/𝐷) shifts to higher Mach numbers. The maximum 𝑀 (𝐿/𝐷) also increases with the design Mach number for the

optimized TTBW-CS-150 but not for the optimized TTBW-150. However, the expected cross-over in the 𝑀 (𝐿/𝐷)

curves as the design Mach number is increased does not occur, leaving the optimized TTBW-CS-150 at a performance

deficit in terms of maximum 𝑀 (𝐿/𝐷) relative to the optimized TTBW-150 for all design Mach numbers. Overall, these

trends are consistent with the results of Vassberg et al. [14], although the penalty in maximum 𝑀 (𝐿/𝐷) is larger than

expected. It is also important to note that Figure 23b presents 𝑀 (𝐿/𝐷) across a range of Mach numbers at constant

<!-- ===== PDF page 29 ===== -->

                  (a) 𝑪𝑫 vs. Mach number                                   (b) 𝑴 (𝑳/𝑫) vs. Mach number

Fig. 23 Drag rise and aerodynamic range efficiency curves for Mach 0.800, 0.815, and 0.830 optimized designs.
Solutions are computed using the LAVA Curvilinear flow solver at 𝑪𝑳 = 0.560 and 𝑹𝒆 = 12.12 × 106 .

lift coefficient and Reynolds number, while Vassberg et al. [14] shows the maximum 𝑀 (𝐿/𝐷) as a function of Mach

number.

                                                  VI. Conclusions
   This paper investigates the aerodynamic design and performance of a Mach 0.800, single-aisle, cruise-slotted TTBW

configuration developed for turbulent flow through the application of aerodynamic shape optimization based on the

RANS equations. Gradient-based optimization is used to minimize vehicle drag at cruise with respect to wing and strut

twist and section shape design variables, while subject to constant lift and minimum wing and strut volume constraints,

and several other geometric constraints introduced to maintain concept feasibility. Optimization is performed on

geometric models consisting of a wing, strut, and fuselage, and include the cruise-slotted TTBW or TTBW-CS-150

configuration, as well as a non-slotted TTBW or TTBW-150 configuration, the latter serving as a reference aircraft.

   Results indicate that at the design point, the optimizer is capable of achieving aerodynamically efficient designs for

both the TTBW-150 and TTBW-CS-150 configurations. This includes smooth wing and strut pressure distributions,

near elliptical spanwise lift distributions, and minimal wave drag. With regard to wave drag, the optimizer was also

found to be successful in mitigating the transonic channel effect unique to truss-braced wings operating at high transonic

Mach numbers typical of single-aisle aircraft. Moreover, mitigation of shock formation in and around the wing-strut

junction of the TTBW-CS-150 configuration involved remarkably similar design features to those of the optimized

TTBW-150, including an outwards force distribution, and novel bumps and insets to control local pressure drops due

to sudden accelerations of the flow. However, comparisons of on-design aerodynamic performance between the two

configurations indicate that the optimized TTBW-CS-150 experiences 3.0% more drag than the optimized TTBW-150,

which translates to a lift-to-drag ratio deficit of 3.0% for the same lift. This was found to be primarily due to higher

skin-friction drag across the wing span, likely arising from the slotted region.

<!-- ===== PDF page 30 ===== -->

    Mach number sensitivity studies are also performed for both the TTBW-150 and TTBW-CS-150 configurations.

Results show that the TTBW-CS-150 configuration does not provide a maximum 𝑀 (𝐿/𝐷) benefit over the TTBW-150,

even at higher operating Mach numbers. This remains to be the case when comparing configurations optimized for

higher design Mach numbers. However, this may be due to the unexpected success of the single-point optimizations in

eliminating nearly all sources of wave drag over the wing and struts, leaving the airfoil technology factor advantage of

the cruise-slotted wing with minimal impact. Nonetheless, for all optimized TTBW-CS-150 designs, the cruise-slotted

wing topology was found to lead to more gradual rises in drag as Mach number is increased compared to the optimized

TTBW-150 configurations. This suggests that the potential advantage of the cruise-slotted wing can manifest itself

if the wing system is designed under more practical considerations, such as performance robustness via multipoint

optimization. This would likely lead to the return of wave drag trades, in which the cruise-slotted wing is well-positioned.

    Overall, aerodynamic shape optimization based on the RANS equations is successful in demonstrating that an

efficient and viable cruise-slotted TTBW configuration can be developed for the single-aisle class, with optimal design

features very much similar to an equivalent non-slotted TTBW configuration. However, more work is required to

minimize the impact of operating at higher cruise Mach numbers on maximum aerodynamic range efficiency when

considering fully turbulent flow conditions. For instance, although the current geometry control approach enables

assessment of the cruise-slotted wing as an integrated system, it is limited in its ability to directly control the shape of

the slot and aft element. More explicit geometry control schemes, such as those employed by Hiller et al. [15, 16],

which directly manipulate the shape, location, and orientation of the aft element, may therefore be required to realize the

performance benefits indicated by Vassberg et al. [14]. Natural laminar flow variants of the cruise-slotted wing are also

of interest due to their potential to further improve cruise aerodynamic efficiency, as demonstrated by Somers [18, 19],

Coder and Somers [20], and by Hiller et al. [15, 16]. Transition prediction can also be incorporated into gradient-based

aerodynamic shape optimization, as demonstrated for transonic swept wings by Husain et al. [56] and for transonic

strut-braced wings by Saadeh et al. [57].

                                                 Acknowledgments
    This work is funded by the NASA Advanced Air Transport Technology (AATT) Project, which is part of the

Advanced Air Vehicles Program (AAVP) within the NASA Aeronautics Research Mission Directorate (ARMD). The

authors gratefully acknowledge early technical contributions from Ercan Dumlupinar and Elisha V. Makarevich, as well

as leadership, guidance, and support from William E. Milholen, Cetin C. Kiris, Jeffrey A. Housman, and David A. Craig

Penner. Appreciation is also extended to David W. Zingg for access to the Jetstream aerodynamic shape optimization

framework, and Brandon M. Lowe and Joshua L. Anibal for internal review. Computational resources were provided by

the NASA High-End Computing Capability (HECC) Program through the NASA Advanced Supercomputing (NAS)

Division at Ames Research Center. The artificial intelligence tool ChatGPT (OpenAI) was used for minor language

<!-- ===== PDF page 31 ===== -->

editing of the manuscript.

                                                        References
 [1] “United States Aviation Climate Action Plan 2021,” United States, retrieved on 15 August 2023. URL https://www.faa.

     gov/sites/faa.gov/files/2021-11/Aviation_Climate_Action_Plan.pdf.

 [2] “NASA ARMD Strategic Implementation Plan 2023,” NASA, retrieved on 15 August 2023. URL https://www.nasa.gov/

     directorates/armd/armd-strategic-implementation-plan/.

 [3] Pfenninger, W., “Design Considerations of Large Subsonic Long Range Transport Airplanes with Low Drag Boundary Layer

     Suction,” Tech. rep., Northrop Aircraft Incorporated, November 1954. NAI-54-800 (BLC-67).

 [4] Bradley, M. K., Droney, C. K., and Allen, T. J., “Subsonic Ultra Green Aircraft Research: Phase II - Volume I - Truss Braced

     Wing Design Exploration,” Tech. rep., Boeing Research and Technology, April 2015. NASA/CR-2015-218704/Volume I.

 [5] Bradley, M. K., and Droney, C. K., “Subsonic Ultra Green Aircraft Research: Phase II - Volume II - Hybrid Electric Design

     Exploration,” Tech. rep., Boeing Research and Technology, April 2015. NASA/CR-2015-218704/Volume II.

 [6] Maldonado, D., Housman, J. A., Duensing, J. C., Jensen, J. C., Kiris, C. C., Viken, S. A., Hunter, C. A., Frink, N. T., and

     McMillin, S. N., “Computational Simulations of a Mach 0.745 Transonic Truss-Braced Wing Design,” AIAA Scitech Forum,

     AIAA 2020-1649, Orlando, Florida, January 2020. https://doi.org/10.2514/6.2020-1649.

 [7] Harrison, N. A., Beyar, M. D., Dickey, E. D., Hoffman, K., Gatlin, G. M., and Viken, S. A., “Development of an Efficient

     Mach = 0.80 Transonic Truss-Braced Wing Aircraft,” AIAA Scitech Forum, AIAA 2020-0011, Orlando, Florida, January 2020.

     https://doi.org/10.2514/6.2020-0011.

 [8] Maldonado, D., Housman, J. A., Piotrowski, M. G. H., Kiris, C. C., Hunter, C. A., Viken, S. A., McMillin, S. N., and

     Milholen, W. E., “Improvements in Simulating a Mach 0.80 Transonic Truss-Braced Wing Configuration using the Spalart-

     Allmaras and k-omega SST Turbulence Models,” AIAA Scitech Forum, AIAA 2021-1531, Virtual Event, January 2021.

     https://doi.org/10.2514/6.2021-1531.

 [9] “Sustainable Flight Demonstrator Project,” NASA, retrieved on 4 February 2023.            URL https://www.nasa.gov/

     aeroresearch/programs/iasp/sfd/description/.

[10] “Sustainable Flight National Partnership,” NASA, retrieved on 4 February 2023.            URL https://www.nasa.gov/

     directorates/armd/sfnp/.

[11] “Boeing Puts X-66 On Ice But Will Continue Thin Wing Studies,” Aviation Week, retrieved on 24 April

     2025.    URL https://aviationweek.com/aerospace/emerging-technologies/boeing-puts-x-66-ice-will-

     continue-thin-wing-studies.

<!-- ===== PDF page 32 ===== -->

[12] Whitcomb, R. T., and Clark, L. R., “An Airfoil Shape for Efficient Flight at Supercritical Mach Numbers,” Tech. rep., NASA,

     July 1965. NASA TM X-1109.

[13] Malone, B., and Mason, W. H., “Multidisciplinary Optimization in Aircraft Design Using Analytic Technology Models,” Aircraft

     Design and Operations Meeting, AIAA 1991-3187, Baltimore, MD, September 1991. https://doi.org/10.2514/6.1991-

     3187.

[14] Vassberg, J. C., Gea, L. M., McLean, J. D., Witowski, D. P., Krist, S. E., and Campbell, R. L., “Slotted Aircraft Wing,” , U.S.

     Patent No. 7,048,228. May 2006.

[15] Hiller, B. R., Campbell, R. L., Banchy, M. N., and Boyett, T. K., “Design Exploration of a Transonic Cruise Slotted Airfoil,”

     AIAA Aviation Forum, AIAA 2021-2525, Virtual, August 2021. https://doi.org/10.2514/6.2021-2525.

[16] Hiller, B. R., Campbell, R. L., and Banchy, M. N., “Transonic Cruise Slotted Wing Design for Commercial Transport Aircraft using

     CDISC,” AIAA SciTech Forum, AIAA 2024-0677, Orlando, FL, January 2024. https://doi.org/10.2514/6.2024-0677.

[17] Vassberg, J., Dehaan, M., Rivers, M., and Wahls, R., “Development of a Common Research Model for Applied CFD

     Validations Studies,” 26th AIAA Applied Aerodynamics Conference, AIAA 2008-6919, Honolulu, HI, August 2008.

     https://doi.org/10.2514/6.2008-6919.

[18] Somers, D. M., “Design of a Slotted, Natural-Laminar-Flow Airfoil for Business-Jet Applications,” Tech. rep., NASA, July

     2012. NASA/CR 2012-217559.

[19] Somers, D. M., “Design of a Slotted, Natural-Laminar-Flow Airfoil for a Transport Aircraft,” Tech. rep., NASA, September

     2019. NASA/CR-2019-220403.

[20] Coder, J. G., and Somers, D. M., “Design of a Slotted, Natural-Laminar-Flow Airfoil for Commercial Transport Applications,”

     Aerospace Science and Technolog, Vol. 106, No. 1, 2020, pp. 106–216. https://doi.org/10.1016/j.ast.2020.106217.

[21] Perkins, C. L., Yang, Z., Mavriplis, D. J., Coder, J. G., Shoemake, L., and Axten, C. J., “Effect of Transition Modeling for

     Analysis of a Slotted, Natural-Laminar-Flow Transonic Truss-Braced Wing Aircraft Configuration,” AIAA SciTech Forum,

     AIAA 2023-2454, National Harbor, MD, January 2023. https://doi.org/10.2514/6.2023-2454.

[22] Secco, N. R., and Martins, J. R. R. A., “RANS-Based Aerodynamic Shape Optimization of a Strut-Braced Wing with Overset

     Meshes,” Journal of Aircraft, Vol. 56, No. 1, 2019, pp. 217–227. https://doi.org/10.2514/1.C034934.

[23] Chau, T., and Zingg, D. W., “Aerodynamic Design Optimization of a Transonic Strut-Braced-Wing Regional Aircraft,” Journal

     of Aircraft, Vol. 59, No. 1, 2022, pp. 253–271. https://doi.org/10.2514/1.C036389.

[24] Chau, T., and Zingg, D. W., “Fuel Burn Evaluation of a Transonic Strut-Braced-Wing Regional Aircraft through Multipoint

     Aerodynamic Optimisation,” The Aeronautical Journal, Vol. 127, No. 1308, 2022, pp. 305–329. https://doi.org/10.

     1017/aer.2022.64.

<!-- ===== PDF page 33 ===== -->

[25] Chau, T., and Zingg, D. W., “Aerodynamic Optimization and Fuel Burn Evaluation of a Transonic Strut-Braced-Wing

     Single-Aisle Aircraft,” Journal of Aircraft, Vol. 60, No. 5, 2023, pp. 1638–1658. https://doi.org/10.2514/1.C037158.

[26] Osusky, L., Buckley, H. P., Reist, T. A., and Zingg, D. W., “Drag Minimization Based on the Navier-Stokes Equations Using a

     Newton-Krylov Approach,” AIAA Journal, Vol. 53, No. 6, 2015, pp. 1555–1577. https://doi.org/10.2514/1.J053457.

[27] Gagnon, H., and Zingg, D. W., “Two-Level Free-Form and Axial Deformation for Exploratory Aerodynamic Shape Optimization,”

     AIAA Journal, Vol. 53, No. 7, 2015, pp. 2015–2026. https://doi.org/10.2514/1.J053575.

[28] Hicken, J. E., and Zingg, D. W., “Aerodynamic Optimization Algorithm with Integrated Geometry Parameterization and Mesh

     Movement,” AIAA Journal, Vol. 48, No. 2, 2010, pp. 400–413. https://doi.org/10.2514/1.44033.

[29] Sederberg, T. W., and Parry, S. R., “Free-Form Deformation of Solid Geometric Models,” ACM SIGGRAPH Computer Graphics,

     Vol. 20, No. 4, 1986, pp. 151–160. https://doi.org/10.1145/15886.15903.

[30] Osusky, M., and Zingg, D. W., “Parallel Newton-Krylov-Schur Solver for the Navier-Stokes Equations Discretized Using

     Summation-By-Parts Operators,” AIAA Journal, Vol. 51, No. 12, 2013, pp. 2833–2851. https://doi.org/10.2514/1.

     J052487.

[31] Spalart, P. R., and Allmaras, S. R., “A One-Equation Turbulence Model for Aerodynamic Flows,” 30th AIAA Aerospace Sciences

     Meeting and Exhibit, AIAA 92-0439, Reno, Nevada, January 1992. https://doi.org/10.2514/6.1992-439.

[32] Allmaras, S. R., Johnson, F. T., and Spalart, P. R., “Modifications and Clarifications for the Implementation of the Spalart-

     Allmaras Turbulence Model,” 7th International Conference on Computational Fluid Dynamics, ICCFD7-1902, Big Island,

     Hawaii, July 2012.

[33] Spalart, P. R., “Strategies for Turbulence Modelling and Simulations,” International Journal of Heat and Fluid Flow, Vol. 21,

     No. 3, 2000, pp. 252–263. https://doi.org/10.1016/S0142-727X(00)00007-2.

[34] Tinoco, E. N., Brodersen, O. P., Keye, S., Laflin, K. R., Feltrop, E., Vassberg, J. C., Mani, M., Rider, B., Wahls, R. A., Morrison,

     J. H., Hue, D., Roy, C. J., Mavriplis, D. J., and Murayama, M., “Summary Data from the Sixth AIAA CFD Drag Prediction

     Workshop: CRM Cases,” Journal of Aircraft, Vol. 55, No. 4, 2018, pp. 1352–1379. https://doi.org/10.2514/1.C034409.

[35] Osusky, M., Boom, P. D., and Zingg, D. W., “Results from the Fifth AIAA Drag Prediction Workshop Obtained with a Parallel

     Newton-Krylov-Schur Flow Solver Discretized using Summation-By-Parts Operators,” 31st AIAA Applied Aerodynamics

     Conference, AIAA 2013-2511, San Diego, CA, June 2013. https://doi.org/10.2514/6.2013-2511.

[36] Del Rey Fernández, D. C., Hicken, J. E., and Zingg, D. W., “Review of Summation-by-Parts Operators with Simultaneous

     Approximation Terms for the Numerical Solution of Partial Differential Equations,” Computers & Fluids, Vol. 95, No. 22, 2014,

     pp. 171–196. https://doi.org/10.1016/j.compfluid.2014.02.016.

<!-- ===== PDF page 34 ===== -->

[37] Saad, Y., and Sosonkina, M., “Distributed Schur Complement Techniques for General Sparse Linear Systems,” SIAM

     Journal on Scientific and Statistical Computing, Vol. 21, No. 4, 1999, pp. 1337–1357. https://doi.org/10.1137/

     S1064827597328996.

[38] Saad, Y., and Schultz, M. H., “GMRES: A Generalized Minimal Residual Algorithm for Solving Nonsymmetric Linear

     Systems,” SIAM Journal on Scientific and Statistical Computing, Vol. 7, No. 3, 1986, pp. 856–869. https://doi.org/10.

     1137/0907058.

[39] Hicken, J. E., and Zingg, D. W., “A Simplified and Flexible Variant of GCROT for Solving Nonsymmetric Linear Systems,”

     SIAM Journal of Scientific Computing, Vol. 32, No. 3, 2010, pp. 1672–1694. https://doi.org/10.1137/090754674.

[40] Gill, P. E., Murray, W., and Saunders, M. A., “SNOPT: An SQP Algorithm for Large-Scale Constrained Optimization,” SIAM

     Journal on Optimization, Vol. 14, No. 4, 2002, pp. 979–1006. https://doi.org/10.1137/S0036144504446096.

[41] Pironneau, O., “On Optimum Design in Fluid Mechanics,” Journal of Fluid Mechanics, Vol. 64, No. 1, 1974, pp. 97–110.

     https://doi.org/10.1017/S0022112074002023.

[42] Jameson, A., “Aerodynamic Design via Control Theory,” Journal of Scientific Computing, Vol. 3, No. 3, 1998, pp. 223–260.

     https://doi.org/10.1007/BF01061285.

[43] Squire, W., and Trapp, G., “Using Complex Variables to Estimate Derivatives of Real Functions,” SIAM Review, Vol. 40, No. 1,

     1998, pp. 110–112. https://doi.org/10.1137/S003614459631241X.

[44] Kiris, C. C., Housman, J. A., Barad, M. F., Brehm, C., Sozer, E., and Moini-Yekta, S., “Computational Framework for

     Launch, Ascent, and Vehicle Aerodynamics (LAVA),” Aerospace Science and Technology, Vol. 55, No. 1, 2016, pp. 189–219.

     https://doi.org/10.1016/j.ast.2016.05.008.

[45] Shur, M. L., Strelets, M. K., Travin, A. K., and Spalart, P. R., “Turbulence Modeling in Rotating and Curved Channels: Assessing

     the Spalart-Shur Correction,” AIAA Journal, Vol. 38, No. 5, 2000, pp. 784–792. https://doi.org/10.2514/2.1058.

[46] Housman, J. A., Kiris, C. C., and Hafez, M. M., “Time-Derivative Preconditioning Methods for Multicomponent Flows - Part I:

     Riemann Problems,” Journal of Applied Mechanics, Vol. 76, No. 2, 2009, pp. 1–13. https://doi.org/10.1115/1.3072905.

[47] Housman, J. A., Kiris, C. C., and Hafez, M. M., “Time-Derivative Preconditioning Methods for Multicomponent Flows - Part

     II: Two-Dimensional Applications,” Journal of Applied Mechanics, Vol. 76, No. 3, 2009, pp. 1–12. https://doi.org/10.

     1115/1.3086592.

[48] Saad, Y., “A Flexible Inner-Outer Preconditioned GMRES Algorithm,” SIAM Journal on Scientific and Statistical Computing,

     Vol. 14, No. 2, 1993, pp. 461–469. https://doi.org/10.1137/0914028.

[49] Chau, T., Anibal, J. L., Lowe, B. M., Machado, L. M., and Duensing, J. C., “High-Fidelity Aeropropulsive Assessment

     of Distributed Electric Propulsion and Boundary Layer Ingestion for the SUSAN Electrofan,” AIAA SciTech Forum, AIAA

     2025-1275, Orlando, FL, January 2025. https://doi.org/10.2514/6.2025-1275.

<!-- ===== PDF page 35 ===== -->

[50] “Strategic Research & Innovation Agenda,” Advisory Council for Aviation Research and Innovation in Europe, retrieved on 20

     October 2020. URL https://www.acare4europe.org/sria.

[51] Vassberg, J. C., Tinoco, E. N., Mani, M., Rider, B., Zickuhr, T., Levy, D. W., Brodersen, O. P., Eisfeld, B., Crippa, S., Wahls,

     R. A., Morrison, J. H., Mavriplis, D. J., and Murayama, M., “Summary of the Fourth AIAA Computational Fluid Dynamics Drag

     Prediction Workshop,” Journal of Aircraft, Vol. 51, No. 4, 2014, pp. 1070–1089. https://doi.org/10.2514/1.C032418.

[52] Kreisselmeier, G., and Steinhauser, R., “Systematic Control Design by Optimizing a Vector Performance Index,” International

     Federation of Active Controls Symposium on Computer-Aided Design of Control Systems, Zurich, Switzerland, August 1979.

     https://doi.org/10.1016/S1474-6670(17)65584-8.

[53] Gagnon, H., and Zingg, D. W., “Euler-Equation-Based Drag Minimization of Unconventional Aircraft,” Journal of Aircraft,

     Vol. 53, No. 5, 2016, pp. 1361–1371. https://doi.org/10.2514/1.C033591.

[54] Nemec, M., Zingg, D. W., and Pulliam, T. H., “Multipoint and Multi-Objective Aerodynamic Shape Optimization,” AIAA

     Journal, Vol. 42, No. 6, 2004, pp. 1057–1065. https://doi.org/10.2514/1.10415.

[55] Kenway, G. K. W., and Martins, J. R. R. A., “Multipoint Aerodynamic Shape Optimization Investigations of the Common

     Research Model Wing,” AIAA Journal, Vol. 54, No. 1, 2016, pp. 113–128. https://doi.org/10.2514/1.J054154.

[56] Husain, F., Simmons, I., and Zingg, D. W., “Application of Aerodynamic Shape Optimization to Swept Natural Laminar Flow

     Wings,” AIAA Aviation Forum and ASCEND, AIAA 2024-3507, Las Vegas, NV, July 2024. https://doi.org/10.2514/6.

     2024-3507.

[57] Saadeh, N., Chau, T., and Zingg, D. W., “Aerodynamic Shape Optimization of a Transonic Strut-Braced-Wing Regional Jet

     With Natural Laminar Flow,” AIAA SciTech Forum, AIAA 2026-0868, Orlando, FL, January 2026. https://doi.org/10.

     2514/6.2026-0868.
