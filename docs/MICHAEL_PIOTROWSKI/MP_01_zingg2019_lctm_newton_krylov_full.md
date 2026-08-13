# Piotrowski & Zingg (2019) — Investigation of a Local Correlation-based Transition Model in a Newton-Krylov Algorithm

> **Full transcript** — verbatim text extracted from the source PDF with
> `pdftotext -layout`, page-delimited, nothing removed except bare page-number
> lines. Layout is preserved, so tables and multi-column material read as they
> appear in the PDF. **Equations are raw PDF extractions**: superscripts,
> subscripts and math symbols may be flattened or split across lines — for any
> equation you intend to implement, confirm against the PDF itself (path below)
> or the digest companion file.
>
> Source PDF: `/home/mdo/Desktop/ARTIGOS SMOTH/PiotrowskiandZingg2019 (2).pdf`

---

<!-- ===== PDF page 1 ===== -->

 Investigation of a Local Correlation-based Transition Model in
                   a Newton-Krylov Algorithm

                                                    Michael G. H. Piotrowski∗
                               Institute for Aerospace Studies, University of Toronto, ON M3H 5T6
                                       NASA Ames Research Center, Moffett Field, CA 94035

                                                         David W. Zingg†
                               Institute for Aerospace Studies, University of Toronto, ON M3H 5T6

              A fully-coupled, implicit Newton-Krylov algorithm with laminar-turbulent boundary-layer
          transition prediction capabilities for three-dimensional external aerodynamic flows has been
          developed. The γ- Re
                             ˜ θt Local Correlation-Based Transition Model with helicity-based crossflow
          correlations is modified and fully coupled to the Spalart-Allmaras turbulence model in a
          parallel Newton-Krylov-Schur flow solver. Deep convergence is obtained using an efficient
          pseudo-transient continuation time marching strategy, consisting of an approximate-Newton
          phase followed by an inexact-Newton phase, with source-term time stepping and exponential
          penalty functions used to address the strong sources and non-smooth functions introduced
          by the transition model. A modification to the Flength function is introduced to improve the
          numerical convergence behaviour of the model. The current implementation is validated
          using two- and three-dimensional transition test cases, with the results obtained using the
          modified transition model comparing well with experiment. The elimination of discontinuities
          in the transition model makes it a suitable candidate for gradient-based shape optimization
          algorithms.

                                                           Nomenclature
~
U           unit velocity vector
U           local velocity magnitude
~
Ω           vorticity vector
Ω           vorticity magnitude
Ωstreamwise streamwise vorticity, Helicity
Hcrossflow crossflow strength parameter
γ           intermittency
 ˜ θt
Re          transported quantity of the local transition onset momentum thickness Reynolds number
Reθt        transition onset momentum thickness Reynolds number
Rev         vorticity Reynolds number
d           minimum off-wall distance
ρ           density
µ           molecular viscosity
µt          eddy viscosity
ui          Cartesian velocity component
Tu          freestream turbulence intensity
λθ          pressure gradient parameter
h           surface roughness height
δ           boundary layer thickness
θt          transition onset momentum thickness
S           strain rate magnitude
γsource     largest positive source-term Jacobian eigenvalue
    ∗ NASA Pathways Ph.D. Student, Computational Aeroscience Branch, AIAA Student member, michael.g.piotrowski@nasa.gov
    † Professor, University of Toronto Distinguished Professor of Computational Aerodynamics and Sustainable Aviation, Associate Fellow AIAA,

dwz@oddjob.utias.utoronto.ca

<!-- ===== PDF page 2 ===== -->

                                                    I. Introduction
     Natural-laminar-flow designs produce an increased laminar extent of the boundary layer, which results in a decrease
in viscous drag [1]. This reduction in viscous drag has the potential to significantly reduce the fuel burn of transport
aircraft, where it has been demonstrated to constitute approximately 50% of the total drag budget [2]. This potential
has been demonstrated by NASA-supported studies which suggest the application of laminar flow control to large
commercial aircraft can reduce aerodynamic drag by approximately 10% [3].
     In the near term, natural laminar flow is being investigated for winglet, tail, and nacelle designs, with the Boeing
787-8 and 777X representing the first commercial applications of natural-laminar-flow nacelles to large transport
aircraft, and the 737 Max aircraft exploiting laminar flow for both the nacelle and winglet designs [4]. However, in
the future, the design of wing sections which exploit significant regions of laminar flow can result in more significant
drag reduction. The Boeing 757 EcoDemonstrator is investigating the effect of anti-insect devices and surface coatings
for the application of laminar flow on wing leading edge surfaces [5], the development of which would improve the
feasibility of such designs.
     Efficient high-fidelity numerical shape optimization algorithms provide the designers of next generation aircraft with
a powerful tool for the development of more fuel-efficient designs [6–9]. The integration of a flow solver with transition
prediction capabilities into an aerodynamic shape optimization algorithm will allow for the study of various trade-offs in
the design of novel natural-laminar-flow wing configurations. In the future, the design of aircraft configurations using
design tools capable of incorporating and exploiting natural laminar flow, which requires the ability to efficiently predict
laminar-turbulent transition, could play a key role in reducing the environmental impact of aviation.
     There are several mechanisms for boundary-layer transition. The method of transition depends greatly on the
application being studied, such as flow past wings, airframes, rotor blades, and wind turbine blades, as well as the
characteristics of the flow, such as Reynolds number, angle of attack, pressure gradient, and turbulence intensity.
For external aerodynamic flows the two dominant mechanisms for transition result from Tollmien-Schlichting wave
growth leading to natural transition and crossflow instabilities resulting from highly swept wings [10]. The effects of
laminar-turbulent transition can be included in numerical simulations using various transition prediction and modeling
techniques. These methods can consist of either local or non-local operations, which either attempt to predict or model
the transition process depending on the fidelity of the method [11, 12].
     Flow solvers based on the Reynolds-Averaged Navier-Stokes equations provide an excellent balance between
accuracy, robustness, and efficiency, which make them a suitable choice for practical engineering applications that
require the solution of large and complex geometries. However, turbulence models used in RANS solvers do not have
the stand-alone capability to predict boundary-layer transition. In order to predict transition, one must apply a transition
criterion. Reviews by Arnal et al. [11, 12] and Pasquale and Rona [13] describe the advantages and disadvantages
of several approaches for coupling transition prediction criteria into RANS solvers. The most widely implemented
transition prediction strategies for practical engineering applications are based on either the e N criterion or on transition
onset functions.
     The e N method for predicting transition is a semi-empirical approach, with the N-factor depending on free-stream
flow conditions, based on linear stability theory. To apply the e N criterion one must first approximate the N-factor
curves which represent the amplification ratios of the unstable frequencies of the disturbances in the boundary-layer.
This can be achieved using stability analysis methods such as the linearized or Parabolized Stability Equations (PSE)
[14], or through the use of simplified e N envelope methods, developed by Drela and Giles [15], and used by Rashad and
Zingg [16] and Mayda [17], depending on the level of fidelity required in a simulation.
     A boundary-layer code is often used to increase the accuracy in determining the boundary-layer edge, which is
required to satisfy the requirements of the stability analysis [18]. However, recent work conducted by Rashad and Zingg
[16] has demonstrated that with sufficient grid resolution, consisting of approximately 100 nodes in the boundary layer,
it is possible to define the boundary-layer edge and extract the boundary-layer properties based on the RANS solution
accurately, without using a boundary-layer solver. Transition models based on the e N method are able to accurately
predict natural-, bypass-, and separation-induced transition. In addition, correlations for crossflow instabilities (such as
the C1 criterion) have been successfully combined with these approaches, with experimental validation demonstrating
accurate transition prediction on transonic swept wings in three dimensions [19–21]. However, e N methods require
non-local boundary-layer quantities that are not directly accessible in RANS-based CFD codes [22], and require a large
infrastructure to apply the code.
     In the context of the current work, the requirement that the transition model be formulated locally is significant.
Modern CFD codes utilize domain decomposition and parallel computation to increase the efficiency of simulations.
Non-local operations, such as integrating quantities over the boundary-layer, introduce difficulties to this process. Of

<!-- ===== PDF page 3 ===== -->

the many transition models developed for RANS-based turbulence models, two which satisfy the requirement of a
fully local formulation are the empirical correlation-based γ- Re˜ θt Local Correlation-based Transition Model (LCTM)
developed by Langtry and Menter [23, 24], and the Amplification Factor Transport (AFT) transition model developed
by Coder and Maughmer [25, 26] and later revised by Coder [27–29].
    The AFT transition model is based on the approximate envelope method for Tollmien-Schlichting wave growth,
developed by Drela and Giles [15], and is able to predict natural and separation-induced transition in subsonic and
transonic flow regimes. The model consists of a single transport equation for the envelope amplification factor ñ with
a locally-derived approximation for the boundary-layer shape factor. The transition model was originally coupled
directly to the one-equation Spalart-Allmaras (SA) eddy-viscosity model [30], resulting in a two-equation coupled
system [25–28]. The most recent model [29] is coupled to an additional, intermittency equation, inspired by the work of
Langtry and Menter [23], resulting in a three-equation coupled system. Simulations of airfoils using the AFT model and
the γ- Re˜ θt LCTM coupled to the two-equation k-ω SST turbulence model demonstrate that the AFT transition model
produces results which more accurately agree with experiment than the LCTM [25, 31]. However, at the time of this
writing correlations for crossflow instabilities have not been published for the AFT transition model.
    The transition model formulated by Langtry and Menter [23] is composed of two transport equations. The first is an
equation for the intermittency, γ, which is used to trigger the turbulence transition process by controlling the level of
turbulence in the boundary layer. The second is an equation for the transition momentum thickness Reynolds number,
 ˜ θt , which in the free-stream is equal to the value produced by the empirical correlations, Reθt , and is diffused into
Re
the boundary layer using a standard diffusion term. Transition due to crossflow instabilities represents the dominant
transition mechanism for swept wings. Of the various crossflow transition models developed for the γ- Re    ˜ θt transition
model, two maintain a fully-local formulation, methods based on the local C1 criterion [32–34] and local helicity
approaches [24, 35]. The local C1 approach uses the solutions of the Falkner-Skan and Cooke equations and is therefore,
unlike methods based on the local helicity, only expected to work for wing-like geometries with high aspect ratios [33].
Helicity, alternatively known as the streamwise vorticity Ωstreamwise , is the magnitude of the scalar product of the local
                                                            ~
                                    ~ and vorticity vector, Ω,
values of the unit velocity vector, U,

                                 ~ = √
                                               u              v                w       
                                 U                       ,√             ,√                ,                             (1)
                                           u2 + v 2 + w 2 u2 + v 2 + w 2 u2 + v 2 + w 2
                                 ~ = ∂w − ∂v , ∂u − ∂w , ∂v − ∂u ,
                                                                        
                                 Ω                                                                                      (2)
                                         ∂ y ∂z ∂z ∂ x ∂ x ∂ y
and can be represented using the non-dimensional crossflow strength Hcrossflow as,

                                                Ωstreamwise = |U   ~
                                                               ~ · Ω|,                                                  (3)
                                                               yΩstreamwise
                                                 Hcrossflow =               .                                           (4)
                                                                   U
    Crossflow instability correlations based on the local helicity were initially developed by Muller and Herbst [35], with
more recent correlations, which include a framework for introducing roughness effects, developed by Langtry et al. [24].
Grabe and Krumbein [33] conducted simulations using correlations based on both the C1 and helicity approaches. The
C1 approach yielded more accurate results in terms of the comparison of predicted and measured transition locations
on infinite swept wing flows, while the local helicity approach produced more accurate transition locations for non
wing-like geometries [33].
    The goal of this work is to develop a RANS-based Newton-Krylov flow solver with laminar-turbulent boundary-
layer transition prediction capabilities, to be used in a high-fidelity aerodynamic shape optimization algorithm for
three-dimensional external aerodynamic flows. The resulting algorithm will be used to investigate the use of natural-
laminar-flow in both conventional, and unconventional aircraft configurations. Therefore, the ability to accurately predict
transition on both conventional and unconventional aircraft designs, where the lifting surfaces can differ drastically from
conventional wing-like geometries, is paramount.
    The local helicity-based correlations developed by Muller and Herbst [35] and Langtry et al. [24] consist of an
additional source term in the Re˜ θt equation, which acts as a sink to reduce the value of Re ˜ θt in regions of significant
crossflow and therefore triggers transition further upstream. While the crossflow correlations of Muller and Herbst [35]
and Langtry et al. [24] are similar, the latter represents a more thoroughly validated model. Furthermore, the ability of
the Langtry model to include roughness effects expands the applicability of the model. For these reasons, the transition
model developed by Langtry et al. [24] was implemented in the current work.

<!-- ===== PDF page 4 ===== -->

                                                   II. Methodology

A. Solution Algorithm
    An efficient, robust, and accurate flow solver is crucial to the performance of an optimization algorithm, as flow
analysis is conducted many times over the course of the optimization process, with a variety of geometries expected. The
flow solver used in the current work is a three-dimensional multi-block structured finite-difference solver developed by
Hicken and Zingg [36] for the solution of the Euler equations, and extended to solve the RANS equations coupled to the
one-equation SA turbulence model by Osusky and Zingg [37]. The governing equations are spatially discretized using
summation-by-parts operators, with simultaneous approximation terms applied to enforce boundary conditions and
inter-block coupling. To decrease the computational time required to complete a flow solution, the computational domain
is decomposed into multiple blocks, resulting in multi-block structured grids, which allow for parallel computations.
    A fully-coupled, implicit Newton-Krylov-Schur solution algorithm making use of a pseudo-transient continuation
(PTC) strategy is applied to the set of nonlinear algebraic equations produced by the spatial discretization of the governing
equations to drive the solution from an initial guess to a converged steady-state solution. The PTC strategy consists of
two phases, the approximate-Newton and inexact-Newton. Globalization is achieved using the approximate-Newton
phase, where an implicit Euler time marching strategy with local time linearization is applied using an approximate,
analytically derived Jacobian. After the residual drops several orders of magnitude, Newton’s method using the full
explicit Jacobian is used to converge the system to a steady-state flow solution. The large system of linear equations
generated at each iteration of the PTC strategy is solved using the GMRES linear solver [38], preconditioned using an
approximate-Schur parallel preconditioner [39].

B. Modified Local Correlation-based Transition Model
    The original γ- Re
                     ˜ θt formulation developed by Langtry and Menter [23] was coupled with the two-equation SST
turbulence model [40]. Subsequent extensions of this model to include crossflow instabilities with correlations based
on the local helicity by Grabe and Krumbein [33], Muller and Herbst [35], and Langtry et al. [24] also included this
coupling. Medida and Baeder [41] modified the Langtry-Menter empirical Local Correlation-based Transition Model
to be compatible with the SA turbulence model, and developed a new set of correlations validated using several test
cases. Recent work demonstrated the γ- Re  ˜ θt -SA model coupled with the crossflow correlations developed by Muller
and Herbst [35] was able to accurately predict transition for several three-dimensional test cases, including a swept
wing, requiring only a slight modification [42]. However, their current formulation contains a function, Gonset , which
requires non-local flow information.
    Recent work by Shengjun et al. [43] investigated the performance of the Local Correlation-based Transition Model
coupled with the SA turbulence model. Their work demonstrated that coupling the γ- Re      ˜ θt transition model with the
helicity-based crossflow correlations developed by Langtry et al. [24] to the SA turbulence model produced accurate
results for two- and three-dimensional flows. However, they did not provide details of the modifications necessary
for this coupling. In the current work, the γ- Re˜ θt transition model is coupled to the SA turbulence model [30] using
the crossflow correlations developed by Langtry et al. [24]. For compatibility, the modifications to the RT and Fwake
functions developed by Medida and Baeder [41] and Schucker [44], respectively, were implemented. In addition,
adjustments to the Fonset2 and Fonset3 terms, inspired by the work of Medida and Beder [41], were applied. The modified
Local Correlation-based Transition Model with helicity-based crossflow correlations coupled to the SA turbulence
model, designated LCTM-SA, is presented in the following three sections.

1. Momentum-thickness Reynolds-number transport equation
   The transport equation for the momentum-thickness Reynolds number is defined below. The convection terms for
both transition model equations have been modified to be consistent with the SA model formulation.

                              ∂ Re
                                 ˜ θt      ∂ Re
                                              ˜ θt                   ∂                ∂ Re
                                                                                          ˜ θt 
                                      + uj         = Pθt + Dscf +        σθt (µ + µt )                                   (5)
                                ∂t          ∂xj                     ∂xj                 ∂xj
                                                ρ
                                      Pθt = cθt (Reθt − Re ˜ θt )(1 − Fθt )                                              (6)
                                                 t
                                                               d 4
                                                                          γ − 1/c  2  
                                                                                   e2
                                      Fθt = min max Fwake e−( δ ) , 1 −                  ,1                              (7)
                                                                            1 − 1/ce2

<!-- ===== PDF page 5 ===== -->

                                               50dΩ               15                   ˜ θt µ
                                                                                      Re
                                          δ=        δ BL ;       δ BL =
                                                                     θ BL ; θ BL =                                   (8)
                                                U                  2                    ρU
                                                                        q
                                                                                dU
                                                    −R e U
                                                           i, j        ρ dU   i    j 2
                                                                           dx j dxi d
                                          Fwake = e 6E +4 ; ReUi, j =                                                (9)
                                                                              µ
                                          cθt = 0.03; σθt = 2.0                                                     (10)
                                              500µ
                                          t=        .                                                               (11)
                                               ρU 2
The natural and bypass transition empirical correlations are defined as follows:
                                                                             
                            ρUθ 
                                        1173.51 − 589.428Tu       +  0.2196
                                                                                F (λ θ ) Tu ≤ 1.3
                   Reθ t =         =                                   T u2
                                     
                                                                                                                    (12)
                             µ        331.50[Tu − 0.5658]−0.671 F (λ θ )
                                     
                                                                                         Tu > 1.3
                                     
                                                                                             T u 1.5
                                      1 − [−12.986λ θ − 123.66λ 2θ − 405.689λ 3θ ]e−[ 1.5 ]            λθ ≤ 0
                          F (λ θ ) = 
                                     
                                     
                                                       [−35λ θ ] ]e−[ T0.5u ]
                                                                                                                    (13)
                                      1 + 0.275[1 − e                                                  λθ > 0

                                                ρθ 2 dU                           1
                                          λθ =           ; U = (u2 + v 2 + w 2 ) 2                                  (14)
                                                 µ ds
                                          dU  u  dU  v  dU  w  dU 
                                             =             +             +                                          (15)
                                          ds      U dx         U dy         U dz
                                          dU                       1
                                                                       du     dv      dw 
                                             = (u2 + v 2 + w 2 ) − 2 · u    +v      +w                              (16)
                                          dx                             dx    dx      dx
                                          dU                 2 − 21
                                                                       du     dv      dw 
                                             = (u + v + w ) · u
                                                 2     2
                                                                            +v      +w                              (17)
                                          dy                             dy    dy      dy
                                          dU                       1
                                                                       du     dv      dw 
                                             = (u2 + v 2 + w 2 ) − 2 · u    +v      +w      .                       (18)
                                          dz                             dz    dz      dz
The local turbulence intensity is not available when using the SA turbulence model. However, Suluksana et al. [45]
demonstrated it can be redundant to use the local value of the turbulence intensity, Tu, as well as the pressure gradient
parameter λ θ . The free-stream turbulence intensity value, Tu∞ , is used throughout the domain in the current work,
which is consistent with the approach taken by Medida and Baeder [41].
   The helicity-based stationary crossflow correlations developed by Langtry et al. are defined below [24]:
                                                            ρ
                                                  Dscf = cθt ccrossflow min(Rescf − Re
                                                                                     ˜ θt , 0)(Fθt2 )               (19)
                                                            t                 
                                                                        d 4
                                                  Fθt2 = min Fwake e−( δ ) , 1                                      (20)
                                              ccrossflow = 0.6                                                      (21)
                                     
                             ρ 0.82
                                U
                                    θt                    h
                   Rescf =                 = −35.088 ln       + 319.51 + f (+∆Hcrossflow ) − f (−∆Hcrossflow )      (22)
                                  µ                        θt
                                                                                     
                                           ∆Hcrossflow = Hcrossflow 1.0 + min RT , 0.4                              (23)
                                          +∆Hcrossflow = max(0.1066 − ∆Hcrossflow, 0)                               (24)
                                      f (+∆Hcrossflow ) = 6200(+∆Hcrossflow ) + 50000(+∆Hcrossflow )    2
                                                                                                                    (25)
                                          −∆Hcrossflow = max(−(0.1066 − ∆Hcrossflow ), 0)                           (26)
                                                                  −∆H           
                                                                       crossflow
                                      f (−∆Hcrossflow ) = 75tanh                   .                                (27)
                                                                     0.0125

<!-- ===== PDF page 6 ===== -->

2. Intermittency transport equation
    The transport equation for intermittency is:

                                      ∂γ      ∂γ                 ∂          µt  ∂γ 
                                         + uj      = Pγ − Eγ +          µ+                                                (28)
                                      ∂t      ∂xj               ∂xj          σf ∂x j
                                         Pγ = Flength ρca1 S γFonset (1 − ce1 γ)
                                                            p
                                                                                                                          (29)
                                          Fonset = max(Fonset2 − Fonset3, 0)                                              (30)
                                                             R 3 
                                                                 T             µt
                                          Fonset3 = max 2 −         , 0 ; RT =                                            (31)
                                                               2.5             µ
                                          Fonset2 = min(max(Fonset1, Fonset1
                                                                             ), 4)                                        (32)
                                                      Reν                         ρd 2
                                          Fonset1 =          ; Reν =      S                                               (33)
                                                   2.193Reθc           µ
                                          ce1 = 1.0; ca1 = 2.0; σ f = 1.0                                                 (34)
                                          Eγ = ca2 ρΩγFturb (ce2 γ − 1)                                                   (35)
                                                        R
                                                      −( 4T ) 4
                                          Fturb = e                                                                       (36)
                                          ce2 = 50;        ca2 = 0.06.                                                    (37)

The separation-induced transition modification is given by:
                                                          Reν                    
                                γsep = min 2.0 · max 0,               − 1 Freattach, 2 Fθt                                (38)
                                                          3.235Reθc
                                               RT 4
                               Freattach = e−( 20 )                                                                       (39)
                                   γeff = max(γ, γsep ).                                                                  (40)

The Flength and Reθc correlations are as follows:

          
          
           398.189 · 10−1 + (−119.270 · 10−4 ) Re  ˜ θt + (−132.567 · 10−6 ) Re˜ 2θt                       ˜ θt < 400
                                                                                                           Re
          
           263.404 + (−123.939 · 10−2 ) Re                            ˜ 2θt + (−101.695 · 10−8 ) Re
                                             ˜ θt + (194.548 · 10−5 ) Re                           ˜ 3θt           ˜ θt < 596
          
                                                                                                           400 ≤ Re
          
          
Flength = 
          
          
                    ˜ θt − 596.0) · 3 · 10−4 )
            0.5 − ( Re                                                                                     596 ≤ Re˜ θt < 1200
          
                                                                                                                    ˜ θt
          
           0.3188                                                                                         1200 ≤ Re
          
          
                                                                                                                          (41)
        
        
          ˜ θt − (396.035 · 10−2 + (−120.656 · 10−4 ) Re
          Re                                                                       ˜ 2θt
                                                         ˜ θt + (868.230 · 10−6 ) Re
                                3                          4
 Reθc =  +(−696.506 · 10−9 ) Re
                               ˜ θt + (174.105 · 10−12 ) Re
                                                          ˜ θt )                                ˜ θt ≤ 1870
        
                                                                                              Re                         (42)
         θt − (593.11 + 0.482( Reθt − 1870.0))
           ˜                        ˜                                                           ˜ θt > 1870.
        
         Re                                                                                   Re
        

   Initial simulations of two- and three-dimensional test cases revealed the existing correlations, when coupled with the
SA turbulence model, consistently predicted transition too far upstream. To remedy this, the Fonset2 and Fonset3 functions
were modified from their original values. The results presented in Section 4 demonstrate that the current model agrees
well with experimental values using this modification.

3. Coupling to SA turbulence model
    The transition model is coupled to the SA turbulence model production term, Pν̃ , through the following:

                                                              P̃ν̃ = γeff Pν̃ .                                           (43)

                                               III. Numerical Behaviour
   The solution algorithm used in the current work consists of a fully-coupled, implicit Newton-Krylov-Schur algorithm.
The implementation of the LCTM-SA transition model equations in such an algorithm produces substantial numerical

<!-- ===== PDF page 7 ===== -->

instabilities. The source terms in the γ- Re
                                           ˜ θt equations are large and highly non-linear. This leads to numerical stiffness
through the presence of restrictive source-term time scales [46, 47]. Several authors [46, 48, 49] have suggested
that fully-implicit approaches can lead to numerical instabilities due to the formation of Jacobians with unfavourable
characteristics. More specifically, it has been demonstrated that the implicit treatment of sources can lead to unbounded
solution updates, and that a modified explicit treatment of sources can improve solver performance [47]. Additionally,
recent work has suggested that fully-coupled approaches may lead to reduced performance due to increased numerical
stiffness [50–53]. It has been shown that relaxing the coupling between the turbulence and mean-flow, or transition and
turbulence equations, can reduce this stiffness, leading to reduced computational cost [53], and improved convergence
[51, 52].
     However, the rate of convergence of Newton’s method is closely linked to the the accuracy of the Jacobian. Neglecting
terms or otherwise manipulating the Jacobian can prevent quadratic convergence [54]. An alternative approach would
be to use a segregated approach in the approximate-Jacobian phase of the Newton-Krylov solver, with the linearization
of the source terms neglected, and use the full Jacobian for the inexact-Newton phase. However, while this has been
demonstrated to be effective for fully-turbulent simulations, the transition model equations are highly sensitive to
changes in the Jacobian. Currently, the switch to the inexact-Newton phase using the fully-coupled approach results in a
slight jump in the transition model residuals as we switch to the full Jacobian. Switching from a segragated approach
to an exact Jacobian has been shown to destabilize the solution further. Rather than modifying the Jacobian matrix,
Lian et al. [47] demonstrated that restricting the local time step based on the non-dimensional source-term number can
improve solver convergence. This smooths the jump from the approximate to the inexact-Newton phase, and allows the
sources to be effectively restrained while using the full Jacobian. A source-term time step restriction, based on the ideas
developed by Lian et al. [47], has been implemented in the current work in order to address these instabilities.
     The LCTM-SA model contains several non-smooth functions in the form of min/max operators, which negatively
affect the numerical behaviour of the model. These non-smooth functions, referred to as simple kinks [55, 56], produce
discontinuities in the Jacobian. This presents a problem for both the Newton-Krylov solution algorithm and the
discrete-adjoint gradient-based optimization algorithm, and can lead to solution stalling. One such function, the Fonset
function, introduces simple kinks through the presence of several layers of min/max operators. Previous studies using the
Langtry-Menter transition model encountered solution stalling as a result of this function and attributed to it additional
stiffness introduced through the coupling of all three turbulence/transition equations, ν̃, γ, Re˜θt [51, 52]. They addressed
this stiffness through relaxing the coupling of the transition and turbulence model. However, in the current work this
stalling has been addressed through the application of smooth approximations.

A. Source-term time stepping
    In their work, Lian et al. [47] demonstrated that an implicit scheme becomes unstable with source-term time step
values greater than one. Maintaining a source-term time step less than unity can prevent unbounded solution updates,
where the source-term time step is determined by the product of the eigenvalues of the source-term Jacobian with
the local time step, λ source ∆t. For the current work the eigenvalues of the block diagonal source-term Jacobian are
determined using a QR algorithm. For nodes where the source-term time step exceeds 0.5, the local time step is
restricted as follows,
                                               if     λ source ∆t ≥ 0.5,                                                (44)
                                             then            ∆t = 0.5/λ source,                                         (45)
where λ source represents the largest positive eigenvalue of the Jacobian of the source-term vector, in this case consisting
of the SA and transition model source terms.

B. Smooth min/max approximation by exponential penalty functions
    Several examples of smooth approximations of non-smooth functions can be found in the literature. Rivest [57]
suggested using a generalized mean-value operator, such as the p-mean,
                                                          P             1
                                                             n
                                                               (g  (x)) p p
                                                            i    i
                                               φ p (x) ,                    .                                      (46)
                                                                   n
As p approaches positive and negative infinity, we recover the maximum and minimum operators, respectively. However,
this approximation only holds if g = (g1, ..., gn ) is a vector of positive real numbers, which limits its applicability.
Another approach is to use an exponential penalty function [58–60], defined as

<!-- ===== PDF page 8 ===== -->

                                                                            max                      0.03
                                                                             (x)                                                            max
                                                                            diff                 0.025                                       (x)
                                                                                                                                            diff
          1.5
                                                                                                     0.02

                                                                                                 0.015
Fonset3

                                                                                       Fonset3
                                                                                                     0.01

                                                                                                 0.005

          0.5                                                                                          0

                                                                                                 -0.005

           0                                                                                      -0.01
                0   0.5   1        1.5    2    2.5          3         3.5          4                  3.14           3.145   3.15   3.155      3.16
                                         RT                                                                                  RT

                          Fig. 1     Exponential penalty function applied to Fonset3 using a p value of 300

                                                                              P                             
                                                                                   n
                                                                        log        i exp(pgi (x))
                                                          φ p (x) ,                                              .                           (47)
                                                                                                 p
Again, as p approaches positive infinity, the function produces the maximum of gi , while p approaching negative
infinity produces the minimum. The advantage of this approach is g = (g1, ..., gn ) can be a vector of any real numbers.
However, a drawback is the presence of log and exponential functions, which while producing a smooth function can
increase the computational cost. This is especially true as gi approaches values close to zero. This can introduce
denormalized floating point numbers which significantly slow the algorithm. To prevent this, a switch has been derived
to activate the exponential penalty function only in the region near the discontinuity. The smooth approximation of a
minimum/maximum operator for two functions, g1 (x) and g2 (x), including this switch is given as follows;

                                                                               log(|p| · pmach )
                                                     if    |g1 (x) − g2 (x)| > −                                                             (48)
                                                                                       |p|
                                               then        φ p (x) = min / max(g1 (x), g2 (x))                                               (49)
                                                                        P                 
                                                                     log in exp(pgi (x))
                                               else        φ p (x) =                         ,                                               (50)
                                                                              p

where pmach represents a value close to machine zero, ie. 10−15 . The Fonset3 function calculated using the exponential
penalty function with a p value of 300, and with the switch active, is illustrated in Figure 1.

1. Flength approximation
    In addition to the minimum/maximum operators, the original Flength formulation introduces simple kinks. In the
current work this is avoided by approximating Flength with a Gaussian model, and using the exponential penalty functions
defined above. The new Flength function is defined below, and presented in Figure 2,

                                                                 X                ˜ θt − bi  2 
                                                                               Re
                                                          t1 =         ai exp −                                                              (51)
                                                                  i
                                                                                     ci
                                                          t 2 = 0.5 − ( Re
                                                                         ˜ θt − 596.0) · 3 · 10−4                                            (52)
                                                          t 3 = φ300 (t 1, t 2 )                                                             (53)
                                                Flength = φ300 (t 3, 0.3188),                                                                (54)

where the Gaussian model constants are defined by,

<!-- ===== PDF page 9 ===== -->

          40                                                                     2
                                                   Flength                                                                   Flength
          35                                       Flength mod                                                               Flength mod
                                                    (x)                                                                          (x)
          30                                                                    1.5

Flength

                                                                      Flength
          20                                                                     1

          10                                                                    0.5

           0                                                                     0
               0      500               1000                  1500               540         560        580          600   620          640
                              Re t                                                                            Re t

                               Fig. 2     Smooth Flength function using a p value of 300

                             Table 1      Flow conditions used for the convergence study

                                               M             Re                  α      Tu (%)       h (µm)
                              subsonic      0.30     5.00 × 106                 1.00      0.20        3.30
                              transonic     0.75     2.00 × 107                 2.00      0.20        3.30

                                        a1 = 38.85;          b1 = −13.94;              c1 = 273.1                                      (55)
                                        a2 = 13.00;          b2 = 262.50;              c2 = 159.5                                      (56)
                                        a3 = 3.699;          b3 = 375.90;              c3 = 81.92.                                     (57)

    Convergence plots of the LCTM-SA model with these modifications, presented in the following section, demonstrate
the improved numerical behaviour of the smooth Flength and min/max functions.

C. Convergence Study
    Simulations of the ONERA M6 geometry [61] were conducted in order to demonstrate the effect of source-term
time stepping and exponential smoothing, including the modified Flength function, on the convergence characteristics of
the LCTM-SA transition model. The ONERA M6 is a single element semispan swept, tapered wing with an aspect ratio
of AR = 3.8 and a taper ratio of λ = 0.562. The leading and trailing edge sweep angles are 30◦ and 15.8◦ , respectively,
resulting in a 26.7◦ sweep angle at the quarter chord. The cross-sectional elements of the wing are based on a symmetric
airfoil using the ONERA D section perpendicular to the 40% line. Two flight conditions, which are representative of
subsonic and transonic flow conditions were chosen, as presented in Table 1. Computational results were obtained using
the LCTM-SA transition model on a 512 block structured mesh containing five million nodes, with a minimum off-wall
spacing of 2 × 10−7 chords, resulting in y + values of 0.07 and 0.23 for the subsonic and transonic cases, respectively.
The unscaled intermittency convergence histories are plotted versus the outer, non-linear iteration count in Figure 3.
    For the subsonic cases without source-term time stepping, large turbulence and transition model updates destabilize
the solution resulting in flow solver divergence. In addition, as the solution advanced, large off-diagonal components
introduced in the LCTM-SA linearization led to poor conditioning of the preconditioner, which is based on an
approximate Jacobian matrix, causing the linear solver to fail. Although the transonic cases did not diverge, negative
intermittency updates caused the solution to stall due to excessive variable clipping to keep intermittency positive. Lian
et al. demonstrated that for source-term time step values greater than unity an implicit method produces solutions with
the opposite sign [47], which helps to explain this behaviour. Enabling source-term time stepping eliminated these

<!-- ===== PDF page 10 ===== -->

                                        1010                                                                                                   1010
                                                                               stts = 0, smth = 0                                                                                      stts = 0, smth = 0
L2 norm of the intermittency residual

                                                                                                       L2 norm of the intermittency residual
                                                                               stts = 0, smth = 1                                                                                      stts = 0, smth = 1
                                                                               stts = 1, smth = 0                                                                                      stts = 1, smth = 0
                                         105                                   stts = 1, smth = 1                                               105                                    stts = 1, smth = 1

                                         100                                                                                                    100

                                        10-5                                                                                                   10-5

                                        10-10                                                                                                  10-10

                                        10-15                                                                                                  10-15
                                                0   200     400          600       800          1000                                                   0   200     400           600       800          1000
                                                           outer iterations                                                                                        outer iterations
                                                          (a) subsonic                                                                                           (b) transonic

Fig. 3 Intermittency variable convergence histories with combinations of source-term time stepping (stts) and
smoothing methods (smth) active (1) and inactive (0)

negative intermittency updates, and improved the conditioning of the preconditioner through a reduction in the local
time step.
    The subsonic case without exponential smoothing failed to converge below a relative residual drop of 10−5 , while
the transonic case converged only three orders of magnitude. Although lift and drag converged to acceptable tolerances,
for numerical shape optimization algorithms deep numerical convergence is necessary. Implementing the exponential
smoothing functions, in addition to the modified Flength function, eliminated the discontinuities in the LCTM-SA model,
resulting in convergence to machine precision. The results presented in the following section use source-term time
stepping and the aforementioned smoothing techniques.

                                                                                            IV. Results
    To validate the current implementation, including the coupling of the LCTM transition model with the SA turbulence
model and the smoothing techniques described above, a range of two- and three-dimensional transition test cases were
simulated. These test cases include natural and separation-induced transition and transition due to crossflow instabilities,
with the results compared to experimental values.

A. Two-Dimensional Cases

1. NACA0012
    Experimental transition locations for the NACA0012 airfoil were obtained by Gregory and O’Reilly [62], with
computational results produced by Johansen and Sorensen [63] and Mosahebi and Laurendeau [52] using an e N method
and the Local Correlation-based Transition Model, respectively. Experimental results were obtained at a Reynolds
number of 2.88 × 106 and Mach number of 0.16, over a range of angles of attack. A free-stream turbulence intensity
value was not given for the wind tunnel. Mosahebi and Laurendeau determined the turbulence intensity as the value
which produced transition locations at zero angle of attack which best matched experimental data [52], resulting in a
turbulence intensity value of 0.50%. However, a Mach number of 0.10 was used in their work. A turbulence intensity of
0.20% was selected for the current work using this approach with a Mach number of 0.16 as described in experiment.
Computational results obtained using the LCTM-SA model are compared with experimental results [62] in Figure 4.
    The computational mesh consists of 60,320 nodes with an off-wall spacing of 1.69 × 10−6 chords, resulting in an
average y + value of 0.18. For angles of attack less than 8◦ the primary transition mechanism is natural transition, with
a laminar separation bubble forming and triggering transition on the upper surface of the airfoil at angles of attack
of 9◦ and 10◦ [62]. The results demonstrate that the computed transition locations generated using the LCTM-SA
model match well with experimental, and previously simulated results, over all angles of attack. However, the drag

<!-- ===== PDF page 11 ===== -->

                            1                                                                                                      10-3
                           0.9
                           0.8
transition location, x/c

                           0.7

                           0.6                                                 Experiment, upper
                                                                               Experiment, lower
                           0.5

                                                                                                        Cd
                                                                               LCTM-SA, upper                              9
                                                                               LCTM-SA, lower
                           0.4
                           0.3
                           0.2
                                                                                                                                                                                    Experiment
                           0.1                                                                                             6                                                        LCTM-SST [52]
                                                                                                                                                                                    LCTM-SA
                            0                                                                                              5
                                 0               2          4           6             8            10                          0          0.2           0.4         0.6      0.8         1          1.2
                                                           angle of attack                                                                                          Cl
                                         (a) Transition locations versus angle of attack                                                                (b) Drag polar
                             10      5                                                                                    102
                                                                                                                                                                                               Cl
                                                                                            SA                                                                                                 Cd
                                                                                            Re t
L2 norm of the residual

                                                                                                                         10-2
                                                                                                        relative error

                            10-5                                                                                         10-4

                                                                                                                         10-6
                                 -10
                                                                                                                         10-8

                           10-15                                                                                         10-10
                                         0           100          200           300              400                               0      50      100     150       200    250     300   350        400
                                                            outer iterations                                                                                  outer iterations
                                               (c) Residual convergence history                                                                (d) Force convergence history

Fig. 4 Experimental and computational results for transitional flow over the NACA0012 airfoil [62], where lift
and drag coefficient error represents the difference between the value at each non-linear iteration and the final,
converged value. Convergence behaviour is presented for the 3◦ angle of attack case

polar illustrates the results under-predict drag relative to experiment. This behaviour was also identified in the work of
Mosahebi and Laurendeau [52] using the original LCTM model.
    The convergence behaviour for the 3◦ angle of attack case is presented in Figures 4 (c) and (d), which illustrate the
L 2 norm of the residuals and the error in the lift and drag coefficients relative to the converged value, respectively. Deep
convergence is obtained for all angles of attack using the local source-term time stepping strategy with exponential
smoothing in the fully-coupled fully-implicit algorithm.

2. NLF0416 General Aviation Airfoil
    The NLF0416 airfoil is a natural-laminar-flow general aviation airfoil designed by Somers [64]. Experimental
results of the NLF0416 airfoil were obtained in the low turbulence pressure tunnel at NASA Langley Research Centre
[64]. Transition locations were measured for several angles of attack at a Mach number of 0.1, and a Reynolds number
of 4.0 × 106 . Unfortunately the turbulence intensity in the wind tunnel used for the experimental study was not specified.
A value of Tu = 0.15% was assumed in order to be consistent with results obtained using the Transition Modelling and
Predictive Capabilities Seminar (TCMPS) guidelines [65].

<!-- ===== PDF page 12 ===== -->

                          1.5                                                                                              1.5

                               1                                                                                             1

                          0.5                                                                                              0.5
                                                                                                                                                                                 Experiment

                                                                                                         Cl
Cl

                                                                                                                                                                                 LCTM-SA

                               0                                                                                             0

                                           Experiment, upper                                                               -0.5
                          -0.5
                                           Experiment, lower
                                           LCTM-SA, upper
                                           LCTM-SA, lower
                                                                                                                            -1
                           -1                                                                                               0.005                           0.01                          0.015
                                   0       0.1    0.2     0.3    0.4    0.5    0.6    0.7          0.8
                                                                                                                                                            Cd
                                                         transition location
                                       (a) Transition locations versus lift coefficient                                                             (b) Drag polar
                           105                                                                                             102
                                                                                                                                                                                          Cl
                                                                                            SA                                                                                            Cd
                                                                                            Re t
L2 norm of the residual

                                                                                                                          10-2
                                                                                                         relative error

                          10-5                                                                                            10-4

                                                                                                                          10-6
                               -10
                                                                                                                          10-8

                          10-15                                                                                           10-10
                                     0             100           200           300               400                              0   50      100     150   200      250   300      350        400
                                                           outer iterations                                                                           outer iterations
                                             (c) Residual convergence history                                                              (d) Force convergence history

Fig. 5 Experimental and computational results for the NLF0416 airfoil, where the convergence history is
presented for the 4◦ angle of attack case [64]

    Computational results were obtained using a mesh consisting of 137,760 nodes with an off-wall spacing of
1.38 × 10−6 chords, resulting in an average y + value of 0.15. Results obtained using the LCTM-SA model are compared
with experimental values in Figure 5. Transition locations obtained using the LCTM-SA model compare well with
experimental values. The drag polar indicates that for negative angles of attack the computed values match experiment;
however the LCTM-SA model begins to over-predict drag for larger angles of attack. Further validation will require a
mesh-refinement study.

3. S809 Wind Turbine Airfoil
    The S809 airfoil is a natural-laminar-flow wind turbine airfoil designed by Somers [66] and studied in the low-
turbulence wind tunnel of the Delft University of Technology Low Speed Laboratory. Experimental results for the S809
airfoil were obtained at a Reynolds number of 2.0 × 106 at angles of attack ranging from 0◦ to 9◦ . The free-stream
turbulence intensity in the wind tunnel varied between Tu = 0.02% and Tu = 0.04%. For a Reynolds number of
2.0 × 106 and an angle of attack of 0◦ , Somers noted the primary transition mechanism on both the upper and lower
surfaces of the airfoil was separation-induced transition due to the presence of laminar separation bubbles [66]. As the
angle of attack increased, the separation bubble on the upper surface decreased in length, disappearing completely at an
angle of attack of 5.13◦ . Numerical results were obtained using a Mach number of 0.10, Reynolds number of 2.0 × 106 ,

<!-- ===== PDF page 13 ===== -->

                           0.6
                                                                                                                        1.2

                           0.5                                                                                           1
transition location, x/c

                           0.4                                                                                          0.8

                                                                                                      Cl
                           0.3                                                                                          0.6

                           0.2                                                                                          0.4

                                             Experiment, upper                                                          0.2
                           0.1               Experiment, lower                                                                                                                     Experiment
                                             LCTM-SA, upper                                                                                                                        LCTM-SA
                                             LCTM-SA, lower                                                               0
                            0                                                                                            0.006       0.008         0.01     0.012    0.014     0.016        0.018
                                 1           2       3       4     5      6     7      8          9
                                                                                                                                                             Cd
                                                            angle of attack
                                         (a) Transition locations versus angle of attack                                                            (b) Drag polar
                                     5                                                                                       2
                             10                                                                                         10
                                                                                                                                                                                            Cl
                                                                                           SA                                                                                               Cd
                                                                                           Re t
L2 norm of the residual

                                                                                                                       10-2
                                                                                                      relative error

                            10-5                                                                                       10-4

                                                                                                                       10-6
                           10-10
                                                                                                                       10-8

                           10-15                                                                                       10-10
                                         0            100          200          300             400                              0   50      100      150    200     250     300      350        400
                                                             outer iterations                                                                         outer iterations
                                                 (c) Residual convergence history                                                         (d) Force convergence history

Fig. 6 Experimental and computational results for the S809 airfoil, where the convergence history is presented
for the 4◦ angle of attack case [66]

and a free-stream turbulence intensity of 0.07%, as recommended by the TCMPS guidelines [65]. .
    Computational results were obtained using a mesh consisting of 137,760 nodes with an off-wall spacing of 1.38×10−6
chords, resulting in an average y + value of 0.10. Results obtained using the LCTM-SA model are compared with
experimental results in Figure 6. As with the other results presented in this section, all cases converged to machine
precision, with the exception of the 9◦ angle of attack case, where a large separated region on the upper leading edge
surface of the airfoil caused the solution to oscillate. The plotted forces were taken as the average value of this oscillation.
In their work, Denison and Pulliam [31] also encountered this large separated region when using the original LCTM
model.

B. NASA NLF2-0415 Infinite Swept Wing
    The NASA NLF2-0415 Infinite Swept Wing consists of an NLF2-0415 airfoil extruded with a 45◦ sweep angle.
Experimental measurements were obtained at an angle of attack of −4◦ , over a range of Mach and Reynolds numbers [67],
and including a variety of surface roughness levels [68]. Transition locations were determined using the naphthalene
visualization technique. Transition locations for a roughness level of 3.3µm, representing a painted surface, were
obtained using the LCTM-SA model and are compared with experimental values [68] in Figure 7. The Reynolds number

<!-- ===== PDF page 14 ===== -->

                           0.9
                                                                               Experiment                           0.01
                                                                               LCTM-SA                                                                             Re = 1.8E6
                           0.8                                                                                                                                     Re = 2.0E6
                                                                                                                                                                   Re = 2.2E6
transition location, x/c

                                                                                                                                                                   Re = 2.5E6
                           0.7                                                                                     0.008
                                                                                                                                                                   Re = 3.0E6
                                                                                                                                                                   Re = 3.5E6
                           0.6
                                                                                                                   0.006

                                                                                                                   Cf
                           0.5

                           0.4                                                                                     0.004

                           0.3

                                                                                                                   0.002
                           0.2
                             1.5               2             2.5          3                 3.5
                                                                                                                                  0.2      0.4         0.6         0.8          1
                                                    Reynolds number                    106                                                       x/c
(a) Upper surface transition locations versus Reynolds number
                                                                                                                               (b) Upper surface skin friction plots
                             105                                                                                    102
                                                                                                                                                                                    CL
                                                                                     SA                                                                                             CD
                                                                                     Re t
L2 norm of the residual

                                                                                                                    10-2
                                                                                                  relative error

                            10-5                                                                                    10-4

                                                                                                                    10-6
                                 -10
                                                                                                                    10-8

                           10-15                                                                                   10-10
                                       0     100       200         300        400         500                              0      100        200             300         400             500
                                                      outer iterations                                                                      outer iterations
                                           (c) Residual convergence history                                                       (d) Force convergence history

Fig. 7 Experimental and computational results for the NLF2-0415 infinite swept wing [68], where the conver-
gence histories are presented for a Reynolds number of 2.2 million

based on free-stream velocity magnitude and chord length was varied from 1.8 × 106 to 3.5 × 106 . The free-stream
turbulence intensity for the wind tunnel was measured as 0.20%. Transition locations are presented only for the upper
surface of the wing, as experimental results are not available for the lower surface.
    The computational mesh used to obtain these results consists of 9.1 million nodes split into 160 blocks, with 81,608
nodes per cross-section, and an off-wall spacing of 0.96 × 10−6 chords, resulting in an average y + value of 0.06. Periodic
boundary conditions were used to simulate the infinite swept wing. The results demonstrate that the LCTM-SA model is
able to closely predict the location of transition for all cases except for the 2 × 106 Reynolds number case, where there is
an error of approximately 14% between the experimental and predicted values, respectively. This error is also present in
the work of Jung and Baeder [42], where they coupled a semi-local correlation-based transition model formulation
with the SA model using the crossflow correlations developed by Muller and Herbst [35], and in the original model
developed by Langtry et al. [24]. In the work produced by Radeztsky et al. [68], the results demonstrate that there is
significant experimental error associated with this case, which helps to explain this discrepancy.
    It is important to note that the choice of numerical dissipation strongly affects the accuracy of the crossflow transition
model. Both the scalar dissipation model developed by Jameson et al. [69] and later refined by Pulliam [70], and the

<!-- ===== PDF page 15 ===== -->

matrix-based dissipation model developed by Swanson and Turkel [71] were investigated for this case. The scalar
dissipation model drastically under predicted the crossflow strength and helicity values. This was determined to be due
to the increased dissipation of the scalar dissipation model in and above the laminar boundary layer, which affected
the vorticity profile within the boundary layer. For turbulent boundary layers this increased dissipation did not have a
significant effect on the vorticity within the boundary layer. For the results presented in the current work, matrix-based
dissipation was used with Vl and Vn values of 0, and a fourth-difference dissipation coefficient of 0.04.

                                                        V. Conclusions
    The γ- Re
            ˜ θt Local Correlation-based Transition Model with helicity-based crossflow correlations developed by
Langtry et al. [24] has been modified to be compatible with the one-equation SA turbulence model, resulting in a
new model designated LCTM-SA. The modifications maintain the fully local formulation pioneered by Langtry and
Menter [23]. Deep convergence is achieved through the use of source-term time stepping and exponential smoothing
functions in a fully-coupled, fully-implicit Newton-Krylov algorithm. Two- and three-dimensional test cases demonstrate
that the LCTM-SA model is able to accurately predict transition due to natural, separation-, and crossflow-induced
transition. The resulting model has been integrated into an aerodynamic shape optimization framework, with future
work investigating the performance of the algorithm in the design of natural-laminar-flow wings.

                                                     Acknowledgments
    This work was partially funded by the NASA Transformational Tools and Technologies (TTT) and Advanced
Air Transport Technology (AATT) projects, the Natural Sciences and Engineering Research Council (NSERC), the
Ontario Graduate Scholarship program, and the University of Toronto. Computations were performed on the Niagara
supercomputer at the SciNet HPC Consortium, a part of Compute Canada. The authors gratefully acknowledge the
input and support provided by Dr. Masayuki Yano, Dr. Thomas Reist, Dr. Gaetan Kenway, Dr. Thomas Pulliam, and Dr.
Cetin Kiris.

                                                          References
 [1] Green, J. E., “Civil aviation and the environment: the next frontier for the aerodynamicist,” Aeronautical Journal, Vol. 110, No.
     1110, 2006, pp. 469–486.
 [2] Bushnell, D. M., “Overview of aircraft drag reduction technology,” AGARD Report 786 (Special Course on Skin Friction Drag
     Reduction), 1992.
 [3] Malik, M. R., Crouch, J. D., Saric, W. S., Lin, J. C., and Whalen, E. A., “Application of Drag Reduction Techniques to Transport
     Aircraft,” Encyclopedia of Aerospace Engineering, 2015.
 [4] Crouch, J. D., “Boundary-Layer Transition Prediction for Laminar Flow Control,” 45th AIAA Fluid Dynamics Conference, No.
     AIAA-2015-2472, Dallas, TX, 2015.
 [5] Week, A., “757 EcoDemo Focuses on Laminar and Active Flow,” http://aviationweek.com/technology/757-ecodemo-
     focuses-laminar-and-active-flow, 2015. Accessed: 2017-05-25.
 [6] Osusky, L., Buckley, H. P., Reist, T. A., and Zingg, D. W., “Drag Minimization Based on the Navier-Stokes Equations Using a
     Newton-Krylov Approach,” AIAA Journal, Vol. 53, No. 6, 2015, pp. 1555–1577.
 [7] Kenway, G. K. W., and Martins, J. R. R. A., “Multipoint High-Fidelity Aerostructural Optimization of a Transport Aircraft
     Configuration,” Journal of Aircraft, Vol. 51, No. 1, 2014, pp. 144–160.
 [8] Reist, T. A., and Zingg, D. W., “High-Fidelity Aerodynamic Shape Optimization of a Lifting-Fuselage Concept for Regional
     Aircraft,” Journal of Aircraft, Vol. 54, No. 3, 2016, pp. 1085–1097.
 [9] Reist, T., Zingg, D. W., Rakowitz, M., Potter, G., and Banerjee, S., “Multifidelity Optimization of Hybrid Wing-Body Aircraft
     with Stability and Control Requirements,” Journal of Aircraft, Advance online publication. DOI: 10.2514/1.C034703, 2018.
[10] Reed, H. L., and Saric, W. S., “Transition mechanisms for transport aircraft,” 38th Fluid Dynamics Conference and Exhibit, No.
     AIAA-2008-3743, Seattle, Washington, 2008.
[11] Arnal, D., Casalis, G., and Houdeville, R., “Practical Transition Prediction Methods: Subsonic and Transonic Flows,” VKI
     Lecture Series: Advances in Laminar-Turbulent Transition Modeling, 2009.
[12] Aupoix, B., Arnal, D., Bezard, H., Chaouat, B., and Chedevergne, F., “Transition and Turbulence Modeling,” AerospaceLab,
     hal-01181225, 2011, pp. 1–13.
[13] Pasquale, D., Rona, A., and Garrett, S. J., “A selective review of CFD transition models,” 39th AIAA Fluid Dynamics Conference,
     No. AIAA-2009-3812, San Antonio, Texas, 2009.
[14] Savill, A. M., “By-pass transition using conventional closures,” Closure strategies for turbulent and transitional flows,
     Cambridge University Press, Cambridge, Vol. 17, 2002, pp. 464–492.

<!-- ===== PDF page 16 ===== -->

[15] Drela, M., and Giles, M. B., “Viscous-inviscid analysis of transonic and low Reynolds number airfoils,” AIAA Journal, Vol. 25,
     No. 10, 1987, pp. 1347–1355.
[16] Rashad, R., and Zingg, D. W., “Aerodynamic shape optimization for natural laminar flow using a discrete-adjoint approach,”
     AIAA Journal, Vol. 54, No. 11, 2016, pp. 3321–3337.
[17] Mayda, E., “Boundary Layer Transition Prediction for Reynolds-Averaged Navier-Stokes Methods,” Ph.D. thesis, University of
     California, Davis, 2007.
[18] Stock, H. W., “e N Transition Prediction in Three-Dimensional Boundary Layers on Inclined Prolate Spheroids,” AIAA Journal,
     Vol. 44, No. 1, 2006, pp. 108–118.
[19] Moens, F., Perraud, J., Krumbein, A., Toulorge, T., Iannelli, P., and Hanifi, A., “Transition prediction and impact on a
     three-dimensional high-lift-wing configuration,” Journal of Aircraft, Vol. 45, No. 5, 2008, pp. 1751–1766.
[20] Streit, T., Horstmann, H., Shraut, G., Hein, S., Fey, U., Egmai, Y., Perraud, J., Salah El Din, I., Cella, U., and Quest, J.,
     “Complementary numerical and experimental data analysis of the ETW Telfona Pathfinder wing transition tests,” 49th Aerospace
     Sciences Meeting and Exhibit, No. AIAA-2011-881, Orlando, Florida, January, 2011.
[21] Shi, Y., Gross, R., Mader, C. A., and Martins, J., “Transition Prediction Based on Linear Stability Theory with the RANS
     Solver for Three-Dimensional Configurations,” 2018 AIAA Aerospace Sciences Meeting, No. AIAA-2018-0819, 2018.
[22] Seyfert, C., and Krumbein, A., “Evaluation of a correlation-Based Transition Model and Comparison with the e N Method,”
     Journal of Aircraft, Vol. 49, No. 6, 2012, pp. 1765–1773.
[23] Langtry, R. B., and Menter, F. R., “Correlation-based transition modeling for unstructured parallelized computational fluid
     dynamics codes,” AIAA Journal, Vol. 47, No. 12, 2009, pp. 2894–2906.
[24] Langtry, R. B., Sengupta, K., Yeh, D. T., and Dorgan, A. J., “Extending the γ- Re ˜ θt Correlation based Transition Model for
     Crossflow Effects,” 45th AIAA Fluid Dynamics Conference, No. AIAA-2015-2474, 2015.
[25] Coder, J. G., and Maughmer, M. D., “A CFD-Compatible Transition Model using an Amplification Factor Transport Equation,”
     51st AIAA Aerospace Sciences Meeting including the New Horizons Forum and Aerospace Exposition, No. AIAA-2013-0253,
     2013.
[26] Coder, J. G., and Maughmer, M. D., “Computational fluid dynamics compatible transition modeling using an amplification
     factor transport equation,” AIAA Journal, Vol. 52, No. 11, 2014, pp. 2506–2512.
[27] Coder, J. G., “Enhancement of the Amplification Factor Transport Transition Modeling Framework,” 55th AIAA Aerospace
     Sciences Meeting, No. AIAA-2017-1709, 2017.
[28] Stefanski, D., Glasby, R., Erwin, J. T., and Coder, J. G., “Development of a Predictive Capability for Laminar-Turbulent
     Transition in HPCMP CREATETM-AV Kestrel Component COFFE using the Amplification Factor Transport Model,” 2018
     AIAA Aerospace Sciences Meeting, No. AIAA-2018-1041, 2018.
[29] Coder, J. G., Pulliam, T. H., and Jensen, J. C., “Contributions to HiLiftPW-3 Using Structured, Overset Grid Methods,” 2018
     AIAA Aerospace Sciences Meeting, No. AIAA-2018-1039, 2018.
[30] Spalart, P. R., and Allmaras, S. R., “A one equation turbulence model for aerodynamic flows,” 30th AIAA Aerospace Sciences
     Meeting and Exhibit, No. AIAA-092-0439, Reno, Nevada, United States 1992.
[31] Denison, M., and Pulliam, T. H., “Implementation and Assessment of the Amplification Factor Transport Laminar-Turbulent
     Transition Model,” 2018 AIAA Aerospace Sciences Meeting, No. AIAA-2018-3382, 2018.
[32] Choi, J. H., and Kwon, O. J., “Enhancement of a Correlation-Based Transition Turbulence Model for Simulating Crossflow
     Instability,” AIAA Journal, Vol. 53, No. 10, 2015, pp. 3063–3072.
[33] Grabe, C., Shengyang, N., and Krumbein, A., “Transition Transport Modeling for the Prediction of Crossflow Transition,” 34th
     AIAA Applied Aerodynamics Conference, No. AIAA-2016-3572, Washington D.C., June, 2016.
[34] Xu, J. K., Bai, J. Q., Qiao, L., and Zhang, Y., “Correlation-Based Transition Transport Modeling for Simulating Crossflow
     Instabilities,” Journal of Applied Fluid Mechanics, Vol. 9, No. 5, 2016, pp. 2435–2442.
[35] Muller, C., and Herbst, F., “Modelling of Crossflow-Induced Transition Based on Local Variables,” 6th European Conference
     on Computational Fluid Dynamics (ECFD), 2014.
[36] Hicken, J. E., and Zingg, D. W., “Parallel Newton-Krylov solver for the Euler equations discretized using simultaneous
     approximation terms,” AIAA Journal, Vol. 46, No. 11, 2008, pp. 2773–2786.
[37] Osusky, M., and Zingg, D. W., “Parallel Newton-Krylov-Schur Flow Solver for the Navier-Stokes Equations,” AIAA Journal,
     Vol. 51, No. 12, 2013, pp. 2833–2851.
[38] Saad, Y., and Schultz, M. H., “GMRES: A generalized minimal residual algorithm for solving nonsymmetric linear systems,”
     SIAM Journal on Scientific and Statistical Computing, Vol. 7, No. 3, 1986, pp. 856–869.
[39] Saad, Y., and Sosonkina, M., “Distributed Schur complement techniques for general sparse linear systems,” SIAM Journal of
     Scientific Computing, Vol. 21, 1999, pp. 1337–1357.
[40] Menter, F. R., “Two-Equation Eddy-Viscosity Turbulence Models for Engineering Applications,” AIAA Journal, Vol. 32, No. 8,
     1994, pp. 1598–1605.
[41] Medida, S., and Baeder, J. D., “Application of the Correlation-based γ- Re    ˜ θt Transition Model to the Spalart-Allmaras
     Turbulence Model,” 20th AIAA Computational Fluid Dynamics Conference, No. AIAA-2011-3979, Honolulu, HI, 2011.
                                          ˜ θt -SA with Crossflow Transition Model using Hamiltonian-Strand Approach,” 2018 AIAA
[42] Jung, Y. S., and Baeder, J. D., “γ- Re
     Aerospace Sciences Meeting, No. AIAA-2018-1040, 2018.

<!-- ===== PDF page 17 ===== -->

[43] Shengjun, J., Chao, Y., and Zhifei, Y., “Integrating Stanford University Unstructured Code with Transition Model,” 7th
     International Conference on Mechanical and Aerospace Engineering, 2016, pp. 573–577.
[44] Schucker, J., “Development of a Three-Equation γ- Re    ˜ θt -Spalart-Allmaras Turbulence-Transition Model,” Master’s thesis,
     DLR, 2012.
[45] Suluksna, K., Dechaumphai, P., and Juntasaro, E., “Correlations for modeling transitional boundary layers under influences of
     freestream turbulence and pressure gradient,” International Journal of Heat and Fluid Flow, Vol. 30, No. 1, 2009, pp. 66–75.
[46] Moryossef, Y., and Levy, Y., “Unconditionally positive implicit procedure for two-equation turbulence models: Application to
     k-ω turbulence models,” Journal of Computational Physics, Vol. 220, 2006, pp. 88–108.
[47] Lian, C., Xia, G., and Merkle, C. L., “Impact of source terms on reliability of CFD algorithms,” Computers & Fluids, Vol. 39,
     2010, pp. 1909–1922.
[48] Moryossef, Y., and Levy, Y., “The unconditionally positive-convergent implicit time integration scheme for two-equation
     turbulence models: Revisited,” Computers & Fluids, Vol. 38, 2009, pp. 1984–1994.
[49] Moryossef, Y., “Unconditionally stable time marching scheme for Reynolds stress models,” Journal of Computational Physics,
     Vol. 276, 2014, pp. 635–664.
[50] Lee, J., and Jameson, A., “Natural-Laminar-Flow Airfoil and Wing Design by Adjoint Method and Automatic Transition
     Prediction,” 47th AIAA Aerospace Sciences Meeting and Exhibit, No. AIAA-2009-897, Orlando, Florida, January, 2009.
[51] Mosahebi, A., and Laurendeau, E., “Convergence Characteristics of Fully and Loosely Coupled Numerical Approaches for
     Transition Models,” AIAA Journal, Vol. 53, No. 5, 2015, pp. 1399–1404.
[52] Mosahebi, A., and Laurendeau, E., “Introduction of a modified segregated numerical approach for efficient simulation of
     γ- Re
         ˜ θt transition model,” International Journal of Computational Fluid Dynamics, Vol. 29, 2015, pp. 357–375.
[53] Candler, G. V., Subbareddy, P. K., and Nompelis, I., “Decoupled Implicit Method for Aerothermodynamics and Reacting
     Flows,” AIAA Journal, Vol. 51, No. 5, 2013, pp. 1245–1254.
[54] Chisholm, T. T., and Zingg, D. W., “A Jacobian-free Newton-Krylov algorithm for compressible turbulent fluid flows,” Journal
     of Computational Physics, Vol. 228, No. 9, 2009, pp. 3490–3507.
[55] Bertsekas, D. P., “Nondifferentiable optimization via approximation,” Mathematical Programming Study 3, 1975, pp. 1–25.
[56] Chakraborty, A., Sinha Roy, A., and Dasgupta, B., “Non-parametric Smoothing for Gradient Methods in Non-differentiable
     Optimization Problems,” 2016 IEEE International Conference on Systems, Man, and Cybernetics, 2016.
[57] Rivest, R. L., “Game Tree Searching by Min/Max Approximation,” Artificial Intelligence, Vol. 34, No. 0004-3702, 1988, pp.
     77–96.
[58] Pee, E. Y., “On Solving Large-Scale Finite Minimax Problems using Exponential Smoothing,” Journal of Optimization Theory
     and Applications, Vol. 148, No. 2, 2011, pp. 390–421.
[59] Polak, E., Royset, J. O., and Womersley, R. S., “Algorithms with Adaptive Smoothing for Finite Minimax Problems,” Journal of
     Optimization Theory and Applications, Vol. 119, No. 3, 2003, pp. 459–484.
[60] Li, X., “An Entropy-Based Aggregate Method for Minimax Optimization,” Engineering Optimization, Vol. 18, No. 0004-3702,
     1997, pp. 277–285.
[61] Schmitt, V., Monneris, B., Dorey, G., and Capelier, C., “Etude De La Couche Limite Tridimensionnelle Sur Une Aile En
     Fleche,” ONERA, Vol. Rept. 14/1713-AN, 1975.
[62] Gregory, N., and O’Reilly, C. L., “Low-Speed Aerodynamic Characteristics of NACA0012 Airfoil Section Including the Effects
     of Upper Surface Roughness Simulating Hoar Frost,” TR, NPL AERO Rept. 1308, Middlesex, England, UK, 1970.
[63] Johansen, J., and Sorensen, J., “Prediction of LaminarTurbulent Transition in Airfoil Flows,” Journal of Aircraft, Vol. 36, No. 4,
     1999, pp. 731–734.
[64] Somers, D. M., “Design and experimental results for a natural-laminar-flow airfoil for general aviation applications,” NASA
     Langley Research Center, Hampton, VA, United States, , No. Technical Paper 1861, 1981.
[65] Coder, J. G., “Standard Test Cases for Transition Model Verification and Validation in Computational Fluid Dynamics,” 2018
     AIAA Aerospace Sciences Meeting, No. AIAA-2018-0029, 2018.
[66] Somers, D. M., “Design and experimental results for the S809 airfoil,” National Renewable Energy Laboratory, Golden,
     Colorado, NRELSR-440-6918, 1989.
[67] Dagenhart, J., and Saric, W., “Crossflow stability and transition experiments in swept-wing flow,” NASA Langley Technical
     Report Server, Vol. NASA/TP-1999-209344, 1999.
[68] Radeztsky, R. H., Reibert, M. S., and Saric, W. S., “Effect of Micron-Sized Roughness on Transition in Swept-Wing Flows,”
     No. AIAA-93-0076, 1993.
[69] Jameson, A., Schmidt, W., and Turkel, E., “Numerical solution of the Euler equations by finite volume methods using
     Runge-Kutta time-stepping schemes,” 14th Fluid and Plasma Dynamics Conference, No. AIAA-81-1259, 1981.
[70] Pulliam, T. H., “Efficient solution methods for the Navier-Stokes equations,” Tech. rep., Lecture Notes for the von Karman Inst.
     for Fluid Dynamics Lecture Series: Numerical Techniques for Viscous Flow Computation in Turbomachinery Bladings, 1986.
[71] Swanson, R. C., and Turkel, E., “On central-difference and upwind schemes,” Journal of Computational Physics, Vol. 101,
     1992, pp. 292–306.
