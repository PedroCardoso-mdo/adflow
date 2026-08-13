# Piotrowski & Zingg (2020) — Smooth Local Correlation-based Transition Model for the Spalart-Allmaras Turbulence Model

> **Full transcript** — verbatim text extracted from the source PDF with
> `pdftotext -layout`, page-delimited, nothing removed except bare page-number
> lines. Layout is preserved, so tables and multi-column material read as they
> appear in the PDF. **Equations are raw PDF extractions**: superscripts,
> subscripts and math symbols may be flattened or split across lines — for any
> equation you intend to implement, confirm against the PDF itself (path below)
> or the digest companion file.
>
> Source PDF: `/home/mdo/Desktop/ARTIGOS SMOTH/PiotrowskiandZingg2020 (7).pdf`

---

<!-- ===== PDF page 1 ===== -->

      Smooth Local Correlation-based Transition Model for the
               Spalart-Allmaras Turbulence Model

                                                 Michael G. H. Piotrowski∗
                                    NASA Ames Research Center, Moffett Field, CA 94035
                            Institute for Aerospace Studies, University of Toronto, ON M3H 5T6

                                                      David W. Zingg†
                            Institute for Aerospace Studies, University of Toronto, ON M3H 5T6

                A smooth version of the γ- Re
                                            ˜ θt local correlation-based transition model (LM2015) for

             the Spalart-Allmaras (SA) turbulence model is presented. The LM2015 model with helicity-

             based crossflow correlations is modified and coupled to the SA turbulence model, designated

             SA-LM2015. The LM2015 and SA-LM2015 transition models include source terms that con-

             tain stiff and non-smooth functions. Approximations to these functions are introduced to

             eliminate discontinuities and improve the numerical behaviour of the model, with the smooth

             model designated SA-sLM2015. Deep convergence is achieved using a fully-coupled, implicit

             Newton-Krylov algorithm globalized using an efficient pseudo-transient continuation strat-

             egy. Modifications to the Newton-Krylov algorithm are introduced, including a source term

             time stepping strategy, to address the large sources introduced by the turbulence and tran-

             sition model equations. Two- and three-dimensional transition test cases demonstrate that

             both models, SA-LM2015 and SA-sLM2015, are able to predict transition due to a variety of

             mechanisms accurately, with the smooth variant displaying significantly improved numerical

             behaviour. Several of the strategies presented, including the smoothing techniques and source

             term time stepping, could be useful in the context of other transition models as well.

                                                      Nomenclature
d              minimum off-wall distance

Hcrossflow     crossflow strength parameter

h              surface roughness height (rms)

M              Mach number

p              smooth maximum/minimum function parameter

pswitch        proximity switch cutoff value, 10−15
    ∗ NASA Pathways Ph.D. Candidate, University of Toronto, NASA Computational Aeroscience Branch, AIAA Student Member,

michael.g.piotrowski@nasa.gov, m.piotrowski@mail.utoronto.ca
    † University of Toronto Distinguished Professor of Computational Aerodynamics and Sustainable Aviation, Associate Fellow AIAA,

dwz@oddjob.utias.utoronto.ca

<!-- ===== PDF page 2 ===== -->

Re          Reynolds number
 ˜ θt
Re          transported quantity of the local transition onset momentum-thickness Reynolds number

Reθt        transition onset momentum-thickness Reynolds number

ReS         strain-rate magnitude Reynolds number

ReΩ         vorticity Reynolds number

RT          eddy viscosity ratio

S           strain rate magnitude

Tu          turbulence intensity
~
U           unit velocity vector

U           local velocity magnitude

ui          Cartesian velocity component
~
Ω           vorticity vector

Ω           vorticity magnitude

Ωstreamwise streamwise vorticity, helicity

γ           intermittency

ρ           density

µ           molecular viscosity

µt          eddy viscosity

λθ          pressure gradient parameter

δ           boundary layer thickness

θt          transition onset momentum thickness

γsource     largest positive source-term Jacobian eigenvalue

φp          smooth maximum/minimum function

φswitch     proximity switch for the smooth maximum/minimum function

                                                  I. Introduction

     Laminar-flow wing designs produce an increased laminar extent of the boundary layer, resulting in a decrease

in viscous drag [1]. These drag savings have the potential to significantly reduce the fuel burn of transport aircraft,

where viscous drag constitutes approximately 50% of the total drag [2]. The potential of laminar-flow designs has been

demonstrated by various studies which suggest the application of laminar flow control to large commercial aircraft can

reduce aerodynamic drag by approximately 10% [3].

<!-- ===== PDF page 3 ===== -->

    In the near term, natural laminar flow (NLF) is being investigated for winglet, tail, and nacelle designs, with the

Boeing 787-8 and 777X representing the first commercial applications of NLF nacelles to large transport aircraft, and

the 737 Max aircraft exploiting laminar flow for both the nacelle and winglet designs [4]. In the future, the design

of wings which exploit significant regions of laminar flow can result in more significant drag reduction. The Airbus

BLADE project is investigating this, while experimental studies of the NASA common research model NLF (CRM-NFL)

demonstrate that it is possible to design transport wings with high sweep and Reynolds numbers with significant regions

of laminar flow [5–7]. However, there are several operational challenges that need to be addressed. For example, the

Boeing 757 EcoDemonstrator is investigating the effect of anti-insect devices and surface coatings for the application of

laminar flow on wing leading edge surfaces.

    Efficient high-fidelity numerical optimization algorithms provide the designers of next generation aircraft with

powerful tools for the development of more fuel-efficient designs [8–10]. The integration of a flow solver with transition

prediction capability into an aerodynamic shape optimization algorithm will allow for the study of various trade-offs

in the design of novel NLF wing configurations. The design of aircraft configurations using design tools capable of

incorporating and exploiting NLF, which requires the ability to predict laminar-turbulent transition efficiently, could

play a key role in reducing the environmental impact of aviation.

    There are several mechanisms for boundary-layer transition. The method of transition depends on the application

being studied, such as flow past wings, airframes, rotor blades, and wind turbine blades, as well as the characteristics

of the flow, such as Reynolds number, angle of attack, pressure gradient, and turbulence intensity. For a swept

transport-aircraft wing the dominant mechanisms for transition result from Tollmien-Schlichting wave growth leading

to natural transition, crossflow instabilities resulting from highly swept wings, concave curvature producing Görtler

instabilities, and attachment-line instabilities as a result of large leading edge radius and sweep at the root of the wing [11].

While the last two mechanisms can be prevented by appropriate profile design, balancing natural and crossflow-induced

transition can be difficult as favourable pressure gradients used to stabilize stream-wise instabilities destabilize crossflow

instabilities [11]. The effects of laminar-turbulent transition can be included in numerical simulations using various

transition prediction and modeling techniques, consisting of either local or non-local operations [12–14].

    Flow solvers based on the Reynolds-Averaged Navier-Stokes (RANS) equations provide an excellent balance between

accuracy, robustness, and efficiency, making them a suitable choice for practical engineering applications that require the

solution of large and complex geometries. However, the turbulence models typically used in RANS solvers do not have

the stand-alone capability to predict boundary-layer transition. In order to predict transition, one must apply a transition

criterion. Reviews by Arnal et al. [12, 13] and Pasquale and Rona [14] describe the advantages and disadvantages

of several approaches for coupling transition prediction criteria into RANS solvers. The most widely implemented

transition prediction strategies for practical engineering applications are based on either the e N criterion or on transition

onset functions.

<!-- ===== PDF page 4 ===== -->

   The e N method for predicting transition is a semi-empirical approach, with the N-factor depending on free-stream

flow conditions, based on linear stability theory. To apply the e N criterion one must first approximate the N-factor

curves which represent the amplification ratios of the unstable frequencies of the disturbances in the boundary-layer.

This can be achieved using stability analysis methods such as the linearized or Parabolized Stability Equations (PSE)

[15, 16], or through the use of simplified e N envelope methods, developed by Drela and Giles [17], and used by Rashad

and Zingg [18] and Mayda [19], depending on the level of fidelity required in a simulation.

   A boundary-layer code can be used to increase the accuracy in determining the boundary-layer edge, which is

required to satisfy the requirements of the stability analysis [20, 21]. However, with sufficient grid resolution it is

possible to define the boundary-layer edge and extract the boundary-layer properties based on the RANS solution

accurately, without using a boundary-layer solver [18, 21]. The influence of the grid resolution on the stability analysis

was studied in detail by Krimmelbein and Radespiel [21].

   Transition models based on the e N method are able to accurately predict natural- and separation-induced transition.

In addition, correlations for crossflow instabilities (such as the C1 criterion [22]) have been successfully combined with

these approaches, with experimental validation demonstrating accurate transition prediction on transonic swept wings

in three dimensions [23–27]. However, e N methods require non-local boundary-layer quantities that are not directly

accessible in RANS-based CFD codes [28].

   In the context of the current work, the requirement that the transition model be formulated locally is significant.

Modern CFD codes utilize domain decomposition and parallel computation to increase the efficiency of simulations.

Non-local operations, such as integrating quantities over the boundary-layer, introduce difficulties to this process. Of the

many transition models developed for RANS-based turbulence models, two which satisfy the requirement of a fully

local formulation are the empirical correlation-based γ- Re
                                                          ˜ θt local correlation-based transition model developed by

Langtry and Menter (LM2009) [29] and the Amplification Factor Transport (AFT) transition model developed by Coder

and Maughmer [30, 31] and later revised by Coder [32–35].

   The AFT transition model is based on the approximate envelope method for Tollmien-Schlichting wave growth,

developed by Drela and Giles [17], and is able to predict natural and separation-induced transition in subsonic and

transonic flow regimes. The model consists of a single transport equation for the envelope amplification factor ñ with

a locally-derived approximation for the boundary-layer shape factor. The transition model was originally coupled

directly to the one-equation Spalart-Allmaras (SA) eddy-viscosity model [36], resulting in a two-equation coupled

system [30–33]. The most recent model [34] is coupled to an additional, intermittency equation, inspired by the work of

Langtry and Menter [29], resulting in a three-equation coupled system. Recent work by Xu et al. [37] extended the

AFT model to include crossflow effects, using the Falkner-Skan and Cooke equations to develop correlations for local

approximations of nonlocal boundary-layer quantities, with promising results. This work involved coupling the model

to the SST turbulence model, resulting in a four-equation system.

<!-- ===== PDF page 5 ===== -->

    The transition model formulated by Langtry and Menter [29] is composed of two transport equations. The first is an

equation for the intermittency, γ, which is used to trigger the turbulence transition process by controlling the level of

turbulence in the boundary layer. The second is an equation for the transition momentum-thickness Reynolds number,
 ˜ θt , which in the free-stream is equal to the value produced by the empirical correlations, Reθt , and is diffused into the
Re

boundary layer using a standard diffusion term. Of the various crossflow transition models developed for the γ- Re
                                                                                                                 ˜ θt

transition model, two maintain a fully-local formulation, methods based on the local C1 criterion [38–41] and local

helicity approaches [42, 43]. The local C1 approach uses the solutions of the Falkner-Skan and Cooke equations and is

therefore only expected to work for wing-like geometries with high aspect ratios [40], unlike methods based on the local

helicity.

    Crossflow instability correlations based on the local helicity were initially developed by Muller and Herbst [42],

with more recent correlations, which include a framework for introducing roughness effects, developed by Langtry et al.

(LM2015) [43]. Grabe and Krumbein [40] conducted simulations using correlations based on both the C1 and helicity

approaches. The C1 approach yielded more accurate results in terms of the comparison of predicted and measured

transition locations on infinite swept wing flows, while the local helicity approach produced more accurate transition

locations for non wing-like geometries [40].

    The local helicity-based correlations developed by Muller and Herbst [42] and Langtry et al. [43] both consist of
                                  ˜ θt equation that acts as a sink to reduce the value of Re
an additional source term in the Re                                                         ˜ θt in regions of significant

crossflow, triggering transition further upstream. While the crossflow correlations of Muller and Herbst [42] and Langtry

et al. [43] are similar, the latter represents a more thoroughly validated model [44]. Furthermore, the ability of the

Langtry model to include roughness effects expands the applicability of the model.

    The original γ- Re
                     ˜ θt LM2009 model developed by Langtry and Menter [29] was coupled with the two-equation SST

turbulence model [45]. Subsequent extensions of this model to include crossflow instabilities with correlations based

on the local helicity by Grabe and Krumbein [40], Muller and Herbst [42], and Langtry et al. [43] also included this

coupling, while Nie et al. [46] investigated the performance of the Langtry-Menter model with C1 and helicity-based

crossflow correlations coupled with a Reynolds stress model.

    Medida and Baeder [47] modified the γ- Re
                                            ˜ θt model to be compatible with the SA turbulence model, and developed

a new set of correlations validated using several test cases, designating this model the γ- Re
                                                                                             ˜ θt -SA transition model.

Recent work demonstrated the γ- Re
                                 ˜ θt -SA model coupled with the crossflow correlations developed by Muller and

Herbst [42] was able to predict transition for several three-dimensional test cases, including a swept wing, requiring

only a slight modification [48]. However, their current formulation contains a function, Gonset , which requires non-local

flow information. More recently, Carnes and Coder [49] successfully coupled the algebraic helicity-based crossflow

correlations developed by Langtry et al. [43] to the AFT model, with the results comparing well to experiment.

    The goal of this work is to develop a RANS-based Newton-Krylov flow solver with laminar-turbulent boundary-layer

<!-- ===== PDF page 6 ===== -->

transition prediction capabilities, to be used in a high-fidelity discrete-adjoint gradient-based aerodynamic shape

optimization algorithm. To achieve this goal, the γ- Re
                                                      ˜ θt LM2015 transition model with helicity-based crossflow

correlations [43] is modified and coupled to the one-equation SA turbulence model, while maintaining a fully-local

formulation, reducing the computational cost of the original four-equation system (Section II). This model will be

referred to as the SA-LM2015 transition model. Smooth approximations of non-differentiable functions are developed

to eliminate discontinuities in the model. Additionally, modifications to numerically stiff functions are produced in

order to improve numerical convergence. This new, smooth model is designated the SA-sLM2015 transition model

(Section III). Numerical methods for addressing the strong source terms present in the SA-LM2015 and SA-sLM2015

models, including a source term time stepping strategy, are developed, and their implementation in the Newton-Krylov

framework is discussed (Section IV). A range of two- and three-dimensional transition test cases are simulated, with

the results demonstrating that the SA-LM2015 and SA-sLM2015 transition models are able to accurately predict

transition due to a variety of transition mechanisms, with the SA-sLM2015 model demonstrating significantly improved

convergence characteristics (Section V).

                                       II. SA-LM2015 Transition Model

   The γ- Re
           ˜ θt LM2015 model with helicity-based crossflow correlations [43] is modified and coupled to the SA

turbulence model [36]. The LM2015 model is decoupled from the SST turbulence model using modifications to the

RT and Fwake functions. Schücker [50] suggests replacing the specific turbulence dissipation rate, ω, in the original

Fwake function with velocity gradients. This approach is adopted in the current work using the strain-rate magnitude.

A simplification of the Fθt function is introduced, which reduces the amount of non-smooth minimum/maximum

operators present in the model. An adjustment to the Fonset1 term was applied to recalibrate the model, and the strain-rate

magnitude in the original intermittency production term is replaced by vorticity magnitude to make the model more

stable near laminar separation bubbles. The separation-induced transition modification developed by Langtry [51] was

neglected, as the delay in the growth of turbulence after laminar separation bubbles produced by the SST turbulence

model is not reproduced by the SA turbulence model, which is demonstrated by the results in Section V.

   The remaining terms and coefficients remain unchanged from the original model. The modified LM2015 model

with helicity-based crossflow correlations coupled to the SA turbulence model, designated SA-LM2015, is presented as

follows.

A. Momentum-thickness Reynolds number transport equation

   The transport equation for the momentum-thickness Reynolds number, non-dimensionalized by freestream quantities

and written in non-conservative form to be consistent with the SA turbulence model, is defined as:

<!-- ===== PDF page 7 ===== -->

                               ∂ Re
                                  ˜ θt      ∂ Re
                                               ˜ θt                  1 ∂                   ∂ Re
                                                                                               ˜ θt 
                                       + uj         = Pθt + Dscf +            σθt (ν + νt )           ,                 (1)
                                 ∂t          ∂xj                     Re ∂ x j                ∂xj
                                              cθt
                                       Pθt =               ˜ θt )(1 − Fθt ),
                                                  (Reθt − Re                                                            (2)
                                               t
                                                           d 4
                                       Fθt = Fwake e−( δ ) ,                                                            (3)
                                             50dΩ                      15                ˜ θt µ 1
                                                                                        Re
                                       δ=         δ BL ;         δ BL =   θ BL ; θ BL =           ,                     (4)
                                              U                         2                 ρU Re
                                                   −Re S               ρd 2 S
                                       Fwake = e 1E +6 ;         ReS =        Re,                                       (5)
                                                                         µ

                                       cθt = 0.03;         σθt = 2.0,                                                   (6)
                                            500µ 1
                                       t=            ,                                                                  (7)
                                             ρU 2 Re

where Re = ρ∞µa∞∞ l , ‘a’ is the sound speed, ‘l’ is the characteristic (reference) length, and ‘∞’ denotes a freestream

quantity. The non-dimensionalization procedure for the transition model equations is consistent with the procedure for

the Navier-Stokes equations, which is described in [52] and is used, for example, in the well-known OVERFLOW flow

solver∗ . For the cases presented in the current work the reference length is specified as the root chord. The natural and

bypass transition empirical correlations are defined as follows:

                                                                      
                                         1173.51 − 589.428Tu + 0.2196
                                     
                            ρUθ
                                     
                                     
                                                                T u 2   F (λ θ )               Tu ≤ 1.3,
                    Reθ t =     Re = 
                                     
                                                                                                                        (8)
                             µ       
                                      331.50[Tu − 0.5658]−0.671 F (λ θ )                       Tu > 1.3,
                                     
                                     
                                     
                                                                                            T u 1.5
                                            1 − [−12.986λ θ − 123.66λ 2θ − 405.689λ 3θ ]e−[ 1.5 ]           λ θ ≤ 0,
                                          
                                          
                                          
                                          
                               F (λ θ ) = 
                                          
                                                                                                                       (9)
                                           1 + 0.275[1 − e[−35λθ ] ]e−[ T0.5u ]
                                          
                                                                                                           λ θ > 0,
                                          

where the local pressure gradient parameter is given by:

                                              ρθ 2 dU                              1
                                          λθ =        Re; U = (u2 + v 2 + w 2 ) 2 ,                                    (10)
                                               µ ds
                                        dU  u  dU  v  dU  w  dU 
                                           =            +              +         ,                                     (11)
                                        ds      U dx         U dy         U dz
                                        dU                       1
                                                                     du     dv      dw 
                                           = (u2 + v 2 + w 2 ) − 2 · u    +v    +w        ,                            (12)
                                        dx                             dx    dx      dx
                                        dU                       1
                                                                      du    dv      dw 
                                           = (u2 + v 2 + w 2 ) − 2 · u    +v    +w        ,                            (13)
                                        dy                             dy    dy      dy
                                        dU                       1
                                                                     du     dv      dw 
                                           = (u2 + v 2 + w 2 ) − 2 · u    +v    +w        .                            (14)
                                        dz                             dz    dz      dz

The local turbulence intensity is not directly available when using the SA turbulence model. However, Suluksana et al.

[53] demonstrated it can be redundant to use the local value of the turbulence intensity, Tu, as well as the pressure
   ∗ https://overflow.larc.nasa.gov/home/users-manual-for-overflow-2-2/, accessed August 2020

<!-- ===== PDF page 8 ===== -->

gradient parameter λ θ , since the effects of pressure gradient are reduced by the decay of the local turbulence intensity.

A favourable pressure gradient accelerates the flow leading to a decrease in the turbulence intensity, and therefore

momentum-thickness Reynolds number, while an adverse pressure gradient leads to a reduction in turbulence intensity.

The free-stream turbulence intensity value, Tu∞ , is used throughout the domain in the current work, which is consistent

with the approach taken by Medida and Baeder [47].
   The helicity-based stationary crossflow correlations developed by Langtry et al. are defined below [43]:

                                                             cθt
                                                    Dscf =                               ˜ θt , 0)(Fθt ),
                                                                 ccrossflow min(Rescf − Re                            (15)
                                                              t
                                                ccrossflow = 0.6,                                                     (16)
                                   
                           ρ 0.82
                              U
                                  θt                         h
                                                                                    +
                 Rescf =                 Re = −35.088 ln            + 319.51 + f (∆Hcrossflow          −
                                                                                              ) − f (∆Hcrossflow ),   (17)
                                µ                            θt
                                                                                                        µt
                                             ∆Hcrossflow = Hcrossflow 1.0 + min RT , 0.4 ;          RT =       ,      (18)
                                                                                                            µ
                                               +
                                             ∆Hcrossflow = max(0.1066 − ∆Hcrossflow, 0),                              (19)
                                               +                    +                     +
                                          f (∆Hcrossflow ) = 6200(∆Hcrossflow ) + 50000(∆Hcrossflow ) 2,              (20)
                                                 −
                                              ∆Hcrossflow  = max(−(0.1066 − ∆Hcrossflow ), 0),                        (21)
                                                                     ∆H −         
                                                                         crossflow
                                                −
                                          f (∆Hcrossflow ) = 75tanh                  .                                (22)
                                                                       0.0125

Helicity, alternatively known as the streamwise vorticity Ωstreamwise , is the magnitude of the scalar product of the local
                                                            ~
                                    ~ and vorticity vector, Ω:
values of the unit velocity vector, U,

                                       
                                    ~ = √       u              v             w          
                                    U                    ,√             ,√                ,                           (23)
                                           u2 + v 2 + w 2 u2 + v 2 + w 2 u2 + v 2 + w 2
                                    ~ = ∂w − ∂v , ∂u − ∂w , ∂v − ∂u ,
                                                                       
                                    Ω                                                                                 (24)
                                         ∂ y ∂z ∂z ∂ x ∂ x ∂ y

and can be represented using the non-dimensional crossflow strength Hcrossflow as,

                                                    Ωstreamwise = |U   ~
                                                                   ~ · Ω|,                                            (25)
                                                                     dΩstreamwise
                                                      Hcrossflow =                .                                   (26)
                                                                         U

B. Intermittency transport equation
   The transport equation for intermittency, also non-dimensionalized and written in non-conservative form, is:

                                        ∂γ      ∂γ              1 ∂       νt  ∂γ 
                                           + uj     = Pγ − Eγ +          ν+           ,                               (27)
                                        ∂t      ∂xj             Re ∂ x j    σf ∂x j

<!-- ===== PDF page 9 ===== -->

                                                                 √
                                        Pγ = ca1 Flength Fonset Ω γ(1 − ce1 γ),                                                (28)

                                        Fonset = max( Fonset2 − Fonset3, 0),
                                                      p
                                                                                                                               (29)
                                                             R 3 
                                                                 T
                                        Fonset3 = max 1 −          ,0 ,                                                        (30)
                                                               2.5
                                        Fonset2 = min(max(Fonset1, Fonset1
                                                                           ), 2),                                              (31)
                                                      ReS
                                        Fonset1 =           ,                                                                  (32)
                                                    2.6Reθc
                                        ce1 = 1.0;     ca1 = 2.0;      σ f = 1.0,                                              (33)

                                        Eγ = ca2 Fturb Ωγ(ce2 γ − 1),                                                          (34)
                                                     RT 4
                                        Fturb = e−( 4 ) ,                                                                      (35)

                                        ce2 = 50;     ca2 = 0.06.                                                              (36)

The Flength and Reθc correlations are as follows:

            398.189 · 10−1 + (−119.270 · 10−4 ) Re  ˜ θt + (−132.567 · 10−6 ) Re˜ 2θt                          ˜ θt < 400,
          
          
          
          
          
                                                                                                             Re
          
          
          
                                                                       ˜ 2θt + (−101.695 · 10−8 ) Re
                                                                                                   ˜ 3θt
          
            263.404 + (−123.939 · 10−2 ) Re  ˜ θt + (194.548 · 10−5 ) Re                                             ˜ θt < 596,
          
          
          
                                                                                                             400 ≤ Re
Flength = 
          
          
                     ˜ θt − 596.0) · 3 · 10−4                                                                        ˜ θt < 1200,
          
          
          
          
          
           0.5 − ( Re                                                                                        596 ≤ Re
          
          
          
          
                                                                                                                      ˜ θt ,
          
           0.3188
                                                                                                             1200 ≤ Re
          
                                                                                                                               (37)

           ˜ θt − (396.035 · 10−2 + (−120.656 · 10−4 ) Re                           ˜ 2θt
                                                          ˜ θt + (868.230 · 10−6 ) Re
        
        
        
        
        
         Re
        
        
        
 Reθc =                       ˜ 3θt + (174.105 · 10−12 ) Re
                                                           ˜ 4θt )
        
        
         +(−696.506 · 10−9 ) Re                                                                    ˜ θt ≤ 1870,
                                                                                                   Re                          (38)
        
        
        
        
           ˜ − (593.11 + 0.482( Re   ˜ θt − 1870.0))                                                ˜ θt > 1870.
        
         Re                                                                                       Re
         θt
        

C. Coupling to SA turbulence model
    The SA-neg-noft2 variant of the SA turbulence model [36, 54] was used in the current work and is presented in

non-conservative, non-dimensional form below:

                           ∂ ν̃      ∂ ν̃                  1  ∂              ∂ ν̃          ∂ ν̃ ∂ ν̃ 
                                + uj      = Pν̃ − Dν̃ +               (ν + ν̃)          + cb2             ,                    (39)
                           ∂t        ∂xj                 σRe ∂ x j             ∂xj            ∂ xi ∂ xi
                                                      cb1               cw1 f w  ν̃  2
                                               Pν̃ =      S̃ ν̃; Dν̃ =                   .                                     (40)
                                                      Re                  Re d

Osusky [55] presents the non-dimensionalization procedure for the SA turbulence model, which is consistent with that

used for the Navier-Stokes equations. The transition model is coupled to the production term, Pν̃ , through the following:

<!-- ===== PDF page 10 ===== -->

                                                         P̃ν̃ = γPν̃ .                                              (41)

                                      III. SA-sLM2015 Transition Model
   The LM2015 and SA-LM2015 models contain several non-smooth functions in the form of min/max operators

and conditional statements. These non-smooth functions, referred to as simple kinks [56, 57], produce discontinuities

in the Jacobian. This presents a problem for both the Newton-Krylov solution algorithm and the discrete-adjoint

gradient-based optimization algorithm. In the current work, this has been addressed through the application of smooth

approximations to min/max operators and conditional statements. In addition, modifications are introduced, such as new

Fonset and Fturb functions, in order to reduce numerical stiffness.
   The modifications below describe the simplification and smoothing of the SA-LM2015 transition model, presented

in Section II. For the remainder of this work it will be designated SA-sLM2015.

A. Minimum/maximum operators
   Several examples of smooth approximations to minimum and maximum operators can be found in the literature.

Rivest [58] suggested using a generalized mean-value operator, such as the p-mean:
                                                      P              p1
                                                         n         p
                                                        i (gi (x))
                                            φ p (x) ,                     .                                         (42)
                                                              n

As p approaches positive and negative infinity, we recover the maximum and minimum operators, respectively. However,

this approximation only holds if g = (g1, ..., gn ) is a vector of positive real numbers, which limits its applicability.

Another approach is to use an exponential penalty function
                                                       P [59–61], defined
                                                                        as:
                                                          n
                                                    log i exp(pgi (x))
                                          φ p (x) ,                     .                                           (43)
                                                             p

Again, as p approaches positive infinity, the function produces the maximum of gi , while p approaching negative infinity

produces the minimum.
   The advantage of this approach is that g = (g1, ..., gn ) can be a vector of any real numbers. A drawback is the

presence of log and exponential functions, which can increase the computational cost. This is especially true as gi

approaches values close to zero. This can introduce denormalized floating point numbers which significantly slow the

algorithm. To prevent this, a proximity switch has been developed (Equation 43) to activate the exponential penalty

function only in the region near the discontinuity. In addition, when considering a two-variable array the exponential

penalty function can be rewritten as:                                                 
                                                       log 1 + exp(p(g2 (x) − g1 (x)))
                                  φ p (x) , g1 (x) +                                       .                        (44)
                                                                         p

In combination with the proximity switch, this form of the exponential penalty function prevents both underflow and

overflow. The proximity switch is given by:

<!-- ===== PDF page 11 ===== -->

Algorithm 1 Smooth maximum/minimum function of two variables
 1: function φ p (g1 (x), g2 (x))
 2:    pswitch ← 10−15                                                                           . Cutoff value for proximity check
 3:    a ← max(g1 (x), g2 (x))                                                                                   . Save larger value
 4:    b ← min(g1 (x), g2 (x))                                                                                 . Save smaller value
 5:    if p > 0 then                                                                                        . Determine maximum
                              log( |p | ·p   )
 6:           if |a − b| > −     |p |
                                      switch
                                             then                                                                . Proximity switch
 7:                φ p (x) ← a                                                                                  . Return maximum
 8:                return
 9:           else                                  
                                  log(1+exp(p(b−a)))
10:                φ p (x) ← a +             p                                             . Calculate exponential penalty function
11:       else                                                                                              . Determine minimum
                             log( |p | ·pswitch )
12:           if |a − b| > −        |p |          then                                                         . Proximity switch
13:                φ p (x) ← b                                                                                 . Return minimum
14:                return
15:           else                                     
                                     log(1+exp(p(a−b)))
16:                φ p (x) ← b +                  p                                        . Calculate exponential penalty function

                                                                   log(|p| · pswitch )
                                                     φswitch , −                       ,                                       (45)
                                                                          |p|

where pswitch represents a cutoff value for the proximity check, which for the current work is set to a value of 10−15 . If

the magnitude of the difference between the two variables, |g1 (x) − g2 (x)|, is greater than φswitch , demonstrating we are

far from the discontinuity, a simple min/max operator is used. As we approach the kink, if the difference between the

two variables is less than φswitch , the exponential penalty function is used.
      Algorithm 1 presents the complete method for approximating a minimum/maximum operator, based on the value of

the smoothing parameter, p, for two variables, g1 (x) and g2 (x), using the exponential penalty functions including the

proximity switch. Figure 1 illustrates the effect of this smoothing by approximating the maximum of a variable and a

constant using various smoothing values, p. For large values of p the region of smoothing is small, and the sharp region

is closely modelled. As p is reduced, the domain where the exponential penalty functions are active increases, increasing

computational cost and producing a smoother approximation of the kink region. Although adaptively updating the value

of p has been shown to help prevent ill-conditioned systems [60], all min/max operators present in the SA-sLM2015

model were approximated using this algorithm with fixed magnitudes of p of 300, φ±300 , which was found through

numerical experimentation to provide a good balance between convergence properties and solution accuracy.

B. Intermittency Source Terms

1. Fonset
      The Fonset function (Equations 29-32) couples the SA, intermittency, and momentum-thickness Reynolds number

equations. It was developed using several layers of minimum and maximum operators to activate the production of

intermittency in the boundary layer, triggering the transition process. Solving the mean flow, turbulence, and transition

<!-- ===== PDF page 12 ===== -->

                                            0.8
                                                    max
                                            0.7      60    0.06
                                            0.6            0.04

                              (0,x-1)

                                            0.5            0.02

                                    p
                              max(0,x-1),
                                            0.4                0
                                                                        0.95        1   1.05
                                            0.3

                                            0.2

                                            0.1

                                             0.5                                1              1.5
                                                                                x

                       Fig. 1 max(0, x − 1) approximated using φ p with various values for p.

equations in a fully-coupled manner results in the Fonset function being particularly problematic, and often the bottleneck

in flow solver convergence. Rashad [62] investigated the design space of several transition region models and noted

that the smoothness of the ramp-up of the eddy-viscosity was found to be particularly important to the smoothness

of the design space. He also demonstrated the importance of achieving deep numerical convergence of the transition

model residual. Therefore, special care was taken in ensuring that the Fonset function be as smooth as possible, while

maintaining the same characteristics of the original function.
    The transition process is triggered when Fonset1 (Equation 32) exceeds unity, and is sustained as the eddy-viscosity

ratio increases. To smooth and simplify this process the Fonset function was reformulated using a hyperbolic tangent

function as follows:
                                                               s
                                                                        ReS  2   2
                                                   Fonset1 =                    + RT ,                                (46)
                                                                       2.6Reθc
                                                               tanh(6(Fonset1 − 1.35)) + 1
                                                   Fonset =                                .                          (47)

Similar to the original, the new Fonset function is activated when Fonset1 exceeds unity; however, once a turbulent

boundary layer has formed Fonset1 remains active. This increases the stability of the model around laminar separation

bubbles, where the solution can oscillate as the bubble grows and shrinks, but prevents the model from being able to

predict relaminarization. For flight conditions typical of transport aircraft we do not expect relaminarization to occur;

however, the reader should be aware of this limitation.

2. Fturb
    In the SA-LM2015 model the Fturb function (Equation 35) remains active until a threshold eddy viscosity ratio is

reached, approximately a value of 4. Until this threshold value is reached, the destruction term remains active and

<!-- ===== PDF page 13 ===== -->

                                                                                                   LM2015
                                         35                                                        sLM2015

                                                                  1.5

                               Flength
                                                                  0.5

                                         15                           0
                                                                      540    560   580      600   620   640

                                              0                 500                   1000                1500

                                          Fig. 2        Original and smooth Flength functions.

interferes with the production term, slowing convergence and often producing oscillations. To remedy this, a new

function was developed which deactivates the destruction term when Fonset becomes active, and rapidly deactivates Fturb

as turbulent eddy viscosity, RT , is produced. It is given as follows,
                                                         Fturb = (1 − Fonset ) exp(−RT ).                             (48)

3. Flength
    The original Flength function (Equation 37) consists of conditional statements and is non-differentiable at the

boundaries. In the current work, this is avoided by approximating Flength using a modified, flexible version of a sigmoid

function known as a generalized logistic function [63]. The new Flength function is defined below and plotted in Figure 2:

                                         Flength1 = exp(−3 · 10−2 ( Re
                                                                     ˜ θt − 460)),                                    (49)
                                                                                       ˜ θt − 596))
                                                              44 − (0.50 − 3 · 10−4 ( Re
                                          Flength = 44 −                                     1
                                                                                                              .       (50)
                                                                            (1 + Flength1 ) 6

4. Reθc
    In previous work approximations of the critical momentum-thickness Reynolds number, Reθc , as linear functions of
 ˜ θt are introduced [47, 48, 50]. In the current work a sinusoidal approximation is developed in order to better represent
Re
           ˜ θt behaviour of the original correlation (Equation 38). The critical momentum-thickness Reynolds number
the lower Re

is simplified from the original version as follows:
                                                                                  ˜
                                                                                Re
                                                                                      θt
                                                                                                
                                                  Reθc = 0.67 Re
                                                               ˜ θt + 24 sin               + 0.5 + 14.                (51)

                                                                                     ˜ θt the new correlation produces
The original and smooth Reθc functions are plotted in Figure 3. For large values of Re

smaller Reθc values relative to the original model; however this difference has been accounted for in the model

<!-- ===== PDF page 14 ===== -->

                              Re c
                                     1000             480         500        520

                                                                                          LM2015
                                                                                          sLM2015
                                            0               500         1000       1500    2000

                                       Fig. 3         Original and smooth Reθc functions.

calibration.

5. Pγ and Eγ

   The maximum value of vorticity magnitude in the intermittency production and destruction terms, Pγ (Equation 28)

and Eγ (Equation 34), was limited to improve stability. Vorticity magnitude is included to obtain the correct dimensions

for the source terms; however, it was found that large values of vorticity, especially near laminar separation bubbles, can

lead to the solution oscillating. Limiting the maximum vorticity value in the intermittency source terms restricts the

relative magnitude of the sources, improving the numerical behaviour of the model. The limiting is a function of the

freestream Mach number in order to account for the non-dimensionalization of velocity by the speed of sound, a∞ , and

the freestream Reynolds number to account for the scaling of vorticity in laminar boundary layers with respect to the

square root of the Reynolds number. The Reynolds number is multiplied by the Mach number to recover the Reynolds

number based on the freestream velocity. It is important to note that while the presence of the Mach number is the

result of the particular non-dimensionalization used, the square root of the freestream Reynolds number based on the

freestream velocity is a physical scaling that is independent of the non-dimensionalization. Furthermore, the vorticity

limiting is introduced only to improve the numerical behaviour of the model and does not aim to improve its predictive

capability. The new source terms are given as follows,

                                                            M √ M Re  √
                               Pγ = ca1 Flength Fonset φ−300 Ω,               γ(1 − ce1 γ),                           (52)
                                                      M √ M Re  
                               Eγ = ca2 Fturb φ−300 Ω,               γ(ce2 γ − 1),                                    (53)

and are illustrated in Figures 4 (a-d) alongside the original SA-LM2015 source terms. The total source terms,

Sγ = Pγ − Eγ , are plotted with fixed Flength and vorticity magnitude, Ω, values of unity versus the transition criterion

<!-- ===== PDF page 15 ===== -->

                         LM2015                                                       sLM2015

                            (a)                                                           (b)

                             (c)                                                          (d)

Fig. 4 Original and smooth intermittency source terms, Sγ = Pγ − Eγ , plotted versus the transition criterion
(Equation 32), intermittency, γ, and the eddy viscosity ratio, RT , with Flength = Ω = 1.

(Equation 32), intermittency, γ, and the eddy viscosity ratio, RT . The reformulated intermittency source terms provide a

smoother design space, while the results in Section V demonstrate that the predicted transition locations with these

modifications match well with the SA-LM2015 model and with experiment.

C. Momentum-thickness Reynolds number Empirical Correlations

1. F (λ θ )

    Similar to the Flength function, the F (λ θ ) correlations (Equation 9) are conditional, depending on the sign of the

pressure gradient parameter, λ θ . Here F (λ θ ) is approximated using exponential penalty functions:

                                                                      Tu
                            F (λ θ )1 = 1 + 0.275[1 − e[−35λθ ] ]e−[ 0.5 ],                                         (54)

<!-- ===== PDF page 16 ===== -->

                                   1.5
                                             LM2015
                                     1       sLM2015

                                   0.5

                                   -0.5                    1.01

                              )
                              F(
                                    -1
                                   -1.5

                                    -2                     0.99

                                   -2.5                                 -1       0         1
                                                                                               10-3
                                    -3
                                     -0.3        -0.2          -0.1          0       0.1                0.2

                Fig. 5   Original and smooth F (λ θ ) functions for a turbulence intensity of 0.05%.

                            F (λ θ )2 = φ300 (F (λ θ )1, 1),                                                         (55)
                                                                                                      T u 1.5
                            F (λ θ )3 = 1 − [−12.986λ θ − 123.66λ 2θ − 405.689λ 3θ ]e−[ 1.5 ] ,                      (56)

                             F (λ θ ) = φ−300 (F (λ θ )2, F (λ θ )3 ).                                               (57)

The original and smooth F (λ θ ) functions are plotted in Figure 5.

                                                IV. Solution Algorithm

A. Parallel Newton-Krylov-Schur flow solver

   An efficient, robust, and accurate flow solver is crucial to the performance of an optimization algorithm, as flow

analysis is conducted many times over the course of the optimization process, with a variety of geometries expected.

The flow solver used in the current work is a three-dimensional multi-block structured finite-difference solver developed

by Hicken and Zingg [64] for the solution of the Euler equations, and extended to the RANS equations coupled to

the one-equation SA turbulence model by Osusky and Zingg [65]. The governing equations are spatially discretized

using summation-by-parts operators, with simultaneous approximation terms applied to enforce boundary conditions

and inter-block coupling. To decrease the computational time required to complete a flow solution, the computational

domain is decomposed into multiple blocks, resulting in multi-block structured grids, which allow for efficient parallel

computations.

   A fully-coupled, implicit Newton-Krylov-Schur solution algorithm making use of a pseudo-transient continuation

strategy is applied to the set of nonlinear algebraic equations produced by the spatial discretization of the governing

equations to drive the solution from an initial guess to a converged steady-state solution. The solution strategy consists

of two phases, an approximate-Newton phase followed by an inexact-Newton phase. Globalization is achieved using the

<!-- ===== PDF page 17 ===== -->

approximate-Newton phase, where an implicit Euler time marching strategy with local time linearization is applied

using an approximate, analytically derived Jacobian. After the residual drops several orders of magnitude, Newton’s

method using a full second-order Jacobian, formed either analytically or using a forward difference, is used to converge

the system to a residual norm of machine zero. The large system of linear equations generated at each iteration in both

phases is solved using the GMRES linear solver [66], preconditioned using an approximate-Schur parallel preconditioner

[67].

B. Considerations for SA-LM2015 and SA-sLM2015 transition models

   Unless suitable steps are taken to modify the basic algorithm, the implementation of the SA-LM2015 and SA-

sLM2015 transition model equations in a Newton-Krylov solver can lead to numerical instability. The source terms in

the γ- Re
        ˜ θt equations are large and highly non-linear. This leads to numerical stiffness through the presence of restrictive

source-term time scales [68, 69]. Several authors [68, 70, 71] have suggested that fully-implicit approaches can lead to

numerical instabilities due to the formation of Jacobians with unfavourable characteristics. More specifically, it has been

demonstrated that the implicit treatment of sources can lead to unbounded solution updates, and that a modified explicit

treatment of sources can improve solver performance [69]. Additionally, recent work has suggested that fully-coupled

approaches may lead to reduced performance due to increased numerical stiffness [72–74]. It has been shown that

relaxing the coupling between the turbulence and mean-flow, or transition and turbulence equations, can reduce this

stiffness, leading to reduced computational cost [74], and improved convergence [72, 73].

   However, the rate of convergence of Newton’s method is closely linked to the the accuracy of the Jacobian. Neglecting

terms or otherwise manipulating the Jacobian can prevent quadratic convergence [75]. Yildirim et al. [76] demonstrated

that segregating the turbulence model in an approximate-Newton phase to globalize Newton’s method can be effective for

fully-turbulent simulations. Similar approaches have been investigated for the SA-LM2015 and SA-sLM2015 models;

however, this can result in the solution stalling in the approximate-Newton phase. The transition model equations

require a tight coupling in order to achieve deep convergence, often requiring stricter linear solver tolerances relative to

fully-turbulent simulations. Moreover, the transition model source terms remain active throughout the inexact-Newton

phase. Rather than modifying the Jacobian matrix, Lian et al. [69] demonstrated that restricting the local time step based

on the non-dimensional source-term number, defined as the product of the largest positive eigenvalue of the source-term

Jacobian matrix and the local time step, can be an effective way of stabilizing large sources in a fully-coupled, implicit

algorithm. A source-term time step restriction based on the ideas developed by Lian et al. [69] has been implemented in

the current work.

   Additional modifications to the Newton-Krylov algorithm were implemented to improve robustness. These include

equation and variable scaling measures, a solution update damping routine, and an unsteady line-search backtracking

algorithm.

<!-- ===== PDF page 18 ===== -->

1. Equation and variable scaling

   The turbulence and transition model variables can vary significantly in magnitude from the mean-flow variables

and do not contain inherent geometric scaling [75]. In addition, large off-diagonal entries are introduced through the

linearization of the turbulence and transition model equations with respect to the mean flow variables. If not accounted

for, these factors can result in improper scaling of the linear system and poor performance of the linear solver. To

address this, row and column scaling measures, including residual auto scaling, are used based on the work by Chisholm

and Zingg [75] and Osusky and Zingg [55, 65]. Instead of solving the linear system of equations directly, a scaled

version of the equations is solved, given by:

                                                        I         
                                                Sa Sr        + A(n) Sc Sc−1 ∆Q(n) = −Sa Sr R(n),                                                                      (58)
                                                        ∆t

where Sa represents the residual auto-scaling matrix and

              23                                                                                                                                  
               Ji                                                                                  1                                               
                           ...                                                                                      ...
                                                                                                                                                    
                                                                                                                                                     
                                2                                                                                                                     
                                 Ji                                                                                       1
      Sri =                                                                           ,    Sci =                                                         ,
                                                                                                                                                         
                                            −1
                                    ν̃max
                                        −1 J3
                                            i                                                                             ν̃max                            
                                                                                                          
                                                             −1                                                                                               
                    
                                                 γmax
                                                    −1 J
                                                        i
                                                                                                           
                                                                                                                                    γmax                       
                                                                              −1                                                                                 
                                                                     −1
                                                                   ˜ θt max J 3                                                             ˜ θt max 
                                                                  Re                                                                         Re
                                                                                                             
                                                                            i                                                                                  

are the i th blocks of the diagonal row and column scaling matrices, respectively. Through a thorough algorithm

optimization procedure, Osusky and Zingg [65] determined an optimal normalization value for ν̃max of 103 , which

accounts for the maximum turbulence value that is likely to be encountered in the flow solve. Following a similar

approach, γmax and Re
                    ˜ θt max values of 10 and 104 , respectively, were selected. The partially scaled residual, Sr R(n) , is

used when monitoring convergence. Full details on the auto, row, and column scaling measures can be found in Osusky

and Zingg [55, 65].

2. Solution update

   The transition model equations can become unstable if the intermittency and transition onset momentum-thickness

Reynolds number are not bounded. Solution-limited time stepping methods, such as those developed by Chisholm

and Zingg [75], and Lian et al. [69, 77], were investigated. However these methods can slow convergence for cases

with stable updates and lead to restrictive time steps for variables approaching zero. In order to ensure stable physical

transition model updates, a solution update damping procedure was implemented, as described in Algorithm 2. Clipping

<!-- ===== PDF page 19 ===== -->

Algorithm 2 Transition model solution update damping
 1: θ fac ← 0.99                                                                                                      . Damping factor
 2: for j = 1, nnodes do                                                                                      . Loop over local nodes
 3:      for l = 1, nvar do                                                                                     . Loop over variables
 4:          Q(n+1)
               j,l
                     = Q(n)
                          j,l
                               + ∆Q(n)
                                     j,l
                                                                                                                     . Update solution
 5:          m=0                                                                                            . Initialize count variable
 6:          while Q(n+1)
                      j,7     > 2 or Q(n+1)
                                         j,7 < 10−10 do                                                          . Limit intermittency
 7:              m = m+1                                                                         . Count number of damping iterations
 8:              Q(n+1)
                  j,7   = Q(n)           m  (n)
                           j,7 + (θ fac ) ∆Q j,7                                                                        . Damp update
 9:          m=0                                                                                            . Initialize count variable
10:          while Q(n+1)
                      j,8   < 20 do                                                                                       . Limit Re
                                                                                                                                   ˜ θt
11:             m = m+1                                                                          . Count number of damping iterations
12:             Q(n+1)
                  j,8     = Q(n)           m  (n)
                             j,8 + (θ fac ) ∆Q j,8                                                                      . Damp update

the transition model variables at these bounds was also investigated; however placing a hard limit on variable bounds

can potentially lead to stalling of the nonlinear solver. It is important to note that with the implementation of the

source-term time stepping strategy, which prevents unstable solution updates, in combination with a more dissipative

first-order upwind discretization strategy for the SA, SA-LM2015, and SA-sLM2015 equations, the upper and lower

bounds specified in Algorithm 2 are rarely reached. Furthermore, steady-state solutions for intermittency respect the

bounds of 0.02 and 1 implicitly defined in the source terms. In addition to the solution update damping, a backtracking

line-search method was implemented based on the work of Modisette [78]

3. Source-term time stepping

      Lian et al. [69] demonstrated that an implicit scheme becomes unstable with source-term time step values greater

than unity. Maintaining a source-term time step less than unity can prevent unbounded solution updates, where the

source-term time step is determined by the product of the largest positive eigenvalue of the source-term Jacobian matrix,

Asource , and the local time step, λ source ∆t (n)
                                               j,k,m
                                                     . For the current work the eigenvalues of the block diagonal source-term

Jacobian are determined using a QR algorithm. The local time step is restricted to ensure the following source-term

time step condition is met:
                                                       ∂Sν̃                ∂Sν̃    ∂Sν̃ 
                                                                                             
                                                        ∂ ν̃               ∂γ     ∂ Re
                                                                                       ˜ θt 
                                                                                               
                                                               ∂Sγ           ∂Sγ     ∂Sγ 
                                                        
                                         Asource =                                           ,
                                                          ∂ ν̃             ∂γ     ∂ Re
                                                                                       ˜ θt 
                                                                                                
                                                            ∂SRe       ∂SRe       ∂SRe
                                                          
                                                                 ˜ θt       ˜ θt        ˜ θ t 
                                                             ∂ ν̃          ∂γ     ∂ Re
                                                                                       ˜ θt 
                                                                                                 

                                                         λ source ∆t (n)
                                                                     j,k,m
                                                                           ≤ 0.9,                                                 (59)

where S represents the local source-term vector, in this case consisting of the SA-LM2015, SA-sLM2015, and SA

turbulence model source terms. Several source-term time step values were investigated, however 0.9 was determined

<!-- ===== PDF page 20 ===== -->

to be an effective balance between speed and robustness. The spatially varying time step and the time step ramping

procedure are described in [65].
   As demonstrated by Lian et al. [69], source-term eigenvalues can remain large throughout convergence. This can

cause the local time step to remain excessively small throughout the solution process, delaying convergence. However,

for most cases the large sources are only destabilizing in the early stages of convergence. To accelerate convergence, a

switch was developed for deactivating the source-term time stepping in the inexact-Newton phase. The source-term time

stepping is active throughout the approximate-Newton phase, and is deactivated after five successive inexact-Newton

iterations for which backtracking is not active. Source-term time stepping is reactivated if backtracking is required, or if

the relative residual drop (Equation 60) increases beyond the residual drop tolerance for the switch from the approximate-

to the inexact-Newton phase, which is typically a value of 10−5 .

                                                       V. Results
   Both transition models, SA-LM2015 and SA-sLM2015, are validated in the following section using a range of two-

and three-dimensional transition test cases. These test cases include natural and separation-induced transition as well as

transition due to crossflow instabilities, with the results compared to experimental values. Convergence histories are

illustrated by the L2-norm of the fully-coupled residual, including the row (equation) scaling measure presented in

Section IV, normalized by the initial residual,
                                                              ||R(n) ||2
                                                    Rd(n) =              ,                                             (60)
                                                              ||R(0) ||2

versus the computational time in seconds. All simulations were conducted on Intel Skylake processors with a core

clock speed of 2.4GHz and one MPI thread per computational block. All results were obtained using a second-order

discretization for the mean flow equations with matrix-based dissipation using Vl and Vn values of 0 and a fourth-difference

dissipation coefficient of 0.04 [79]. First-order upwinding was used for the SA, SA-LM2015, and SA-sLM2015 models.

The same solution algorithm parameters were used for the 2D and 3D cases, respectively, demonstrating the robustness

of the solver. These settings are similar to the values recommended by Osusky and Zingg [55], with the settings for the

source term time stepping strategy and scaling measures described previously in Section IV.

A. NLF0416 General Aviation Airfoil
   The NLF0416 airfoil is a natural-laminar-flow general aviation airfoil designed by Somers [80]. Experimental

results of the NLF0416 airfoil were obtained in the low turbulence pressure tunnel (LTPT) at NASA Langley Research

Center [80]. Transition locations were measured for several angles of attack at a Mach number of 0.1 and a Reynolds

number of 4.0 × 106 . The transition mechanisms for these conditions are natural and separation-induced transition.

The turbulence intensity in the wind tunnel used for the experimental study was not specified, however a value of Tu

= 0.15% was assumed in order to be consistent with results obtained using the Transition Modelling and Predictive

<!-- ===== PDF page 21 ===== -->

                           Table 1      TMPCS NLF0416 structured C-grid dimensions [81]

                            grid level     dimensions   ∆s (chord)      average/maximum y+
                             coarse         529 × 73    4.47 × 10−6           0.57/1.19
                            medium          705 × 97    3.25 × 10−6           0.51/0.89
                               fine        1057 × 145   2.10 × 10−6           0.26/0.45

Capabilities Seminar (TCMPS) guidelines [81].

   The coarse, medium, and fine NLF0416 family of grids provided by the TMPCS were used in the current work [81].

The dimensions of the structured C-grids are presented in Table 1. A grid-refinement study was conducted for the 0◦

angle of attack case using both transition models and is presented in Figure 6. The refinement study demonstrates a clear

second-order convergence trend, as the inverse of the number of grid points, N −1 , is proportional to the square of the

mesh spacing, although they converge to slightly different lift and drag values due to a difference in predicted transition

location. Figure 6 also illustrates the computational results obtained using the SA-LM2015 and SA-sLM2015 transition

models on the fine grid level compared to experimental values. The results demonstrate that both the SA-LM2015

model and the SA-sLM2015 model do a reasonable job of matching the forces and transition locations from experiment.

Figure 6(e) illustrates that both models under-predict drag at low angles of attack; however this is consistent with results

obtained from previous studies using variations of the γ- Re
                                                           ˜ θt model [48, 82, 83] and does not appear to be caused by

incorrect transition point locations.

   The residual convergence histories for these cases are presented in Figure 7. For these simulations the grids were

split into 8 computational blocks. The convergence histories demonstrate the improved numerical behaviour of the

SA-sLM2015 model. At angles of attack of −6◦ and −4◦ a large separation bubble forms on the lower surface of the

airfoil, increasing the stiffness of the linear system. This is evident in the decreased slope of the residual curve for the

SA-sLM2015 model, as each outer inexact-Newton iteration required more inner GMRES linear solver iterations. It

should also be noted that small trailing edge-spacings on the C-grids provided by the TMPCS propagated to the far-field

of the domain with little splaying, further increasing the stiffness of the linear system and increasing the computational

cost of the simulations.

B. S809 Wind Turbine Airfoil

   The S809 airfoil is a natural-laminar-flow wind turbine airfoil designed by Somers [84] and studied in the low-

turbulence wind tunnel of the Delft University of Technology Low Speed Laboratory. Experimental results for the S809

airfoil were obtained at a Reynolds number of 2.0 × 106 at angles of attack ranging from 0◦ to 9◦ . The free-stream

turbulence intensity in the wind tunnel varied between Tu = 0.02% and Tu = 0.04%. For a Reynolds number of

2.0 × 106 and an angle of attack of 0◦ , Somers noted that the primary transition mechanism on both the upper and lower

<!-- ===== PDF page 22 ===== -->

                                 SA-LM2015fine                                                                SA-LM2015cl
         0.008                   SA-LM2015medium                                .0055                         SA-LM2015cd                 .4900
                                 SA-LM2015coarse                                                              SA-sLM2015cl
                                 SA-sLM2015fine                                                               SA-sLM2015cd
         0.006                   SA-sLM2015medium                               .0054                                                     .4890
                                 SA-sLM2015coarse

                                                                                .0053                                                     .4880
         0.004

                                                                            Cd

                                                                                                                                             Cl
                             upper
        Cf

                                                                                .0052                                                     .4870
         0.002

                                                       lower                    .0051                                                     .4860

                                                                                .0050
                   0       0.2       0.4         0.6         0.8       1                    1E-05          1.5E-05   2E-05        2.5E-05
                                           x/c                                                                 N-1
                       (a) Skin friction grid convergence                               (b) Grid convergence for lift and drag

                                                         SA-LM2015fine                                     SA-LM2015fine
                                                         SA-sLM2015fine           1.5                      SA-sLM2015fine
              2                                          LTPT Experiment                                   LTPT Experiment

             1.5                                                                   1

                                                                                  0.5
        Cl

                                                                            Cl

             0.5

                                                                                 -0.5
                        .0050      .0100         .0150         .0200                            -5              0            5
                                           Cd                                                                 aoa

                                  (c) Drag polar                                                       (d) Lift curve

                                                         SA-LM2015fine                                                  SA-LM2015fine
                                                         SA-sLM2015fine                                                 SA-sLM2015fine
         .0090                                           LTPT Experiment          1.5                                   LTPT Experiment

         .0080                                                                     1
                                                                                           upper                                  lower
        Cd

                                                                                  0.5
                                                                            Cl

         .0070

         .0060

                                                                                 -0.5
         .0050

                            -5             0             5                                 0         0.2       0.4    0.6        0.8
                                           aoa                                                 transition location (x/c)
                                  (e) Drag curve                                               (f) Transition locations

Fig. 6 NLF0416 grid convergence at 0◦ angle of attack, and angle of attack sweep simulated using the SA-
LM2015 and SA-sLM2015 transition models at M = 0.1, Re = 4.0 × 106 , and Tu = 0.15%. The ‘fine’ grid level
was used for the sweep, with the experimental results obtained at the LTPT wind tunnel [80].

<!-- ===== PDF page 23 ===== -->

                         101                                                                101                                                                101
                                                      SA-LM2015                                                          SA-LM2015                                                          SA-LM2015
                        10−1                          SA-sLM2015                           10−1                          SA-sLM2015                           10−1                          SA-sLM2015
Relative Residual, d

                                                                   Relative Residual, d

                                                                                                                                      Relative Residual, d
                        10−3                                                               10−3                                                               10−3
                        10−5                                                               10−5                                                               10−5
                        10−7                                                               10−7                                                               10−7
                        10−9                                                               10−9                                                               10−9
                   10−11                                                              10−11                                                              10−11
                   10−13                                                              10−13                                                              10−13
                               0   3000   6000   9000 12000                                       0   3000   6000   9000 12000                                       0   3000   6000   9000 12000
                                            Time(s)                                                           Time(s)                                                            Time(s)
                                    (a) aoa = −6◦                                                      (b) aoa = −4◦                                                      (c) aoa = −2◦
                         101                                                                101                                                                101
                                                      SA-LM2015                                                          SA-LM2015                                                          SA-LM2015
                        10−1                          SA-sLM2015                           10−1                          SA-sLM2015                           10−1                          SA-sLM2015
Relative Residual, d

                                                                   Relative Residual, d

                                                                                                                                      Relative Residual, d
                        10−3                                                               10−3                                                               10−3
                        10−5                                                               10−5                                                               10−5
                        10−7                                                               10−7                                                               10−7
                        10−9                                                               10−9                                                               10−9
                   10−11                                                              10−11                                                              10−11
                   10−13                                                              10−13                                                              10−13
                               0   3000   6000   9000 12000                                       0   3000   6000   9000 12000                                       0   3000   6000   9000 12000
                                           Time(s)                                                             Time(s)                                                            Time(s)
                                     (d) aoa = 0◦                                                       (e) aoa = 2◦                                                       (f) aoa = 4◦
                                                                                                                         SA-LM2015
                                                                                           10−1                          SA-sLM2015
                                                                   Relative Residual, d

                                                                                           10−3
                                                                                           10−5
                                                                                           10−7
                                                                                           10−9
                                                                                      10−11
                                                                                      10−13
                                                                                                  0   3000   6000   9000 12000
                                                                                                              Time(s)
                                                                                                        (g) aoa = 6◦

Fig. 7 NLF0416 residual convergence histories for the SA-LM2015 and SA-sLM2015 transition models on
the ‘fine’ grid level split into 8 computational blocks. The airfoil was simulated at M = 0.1, Re = 4.0 × 106 ,
Tu = 0.15% and over a range of angles of attack.

surfaces of the airfoil was separation-induced transition due to the presence of large laminar separation bubbles [84]. As

the angle of attack increased, the separation bubble on the upper surface decreased in length, disappearing completely at

an angle of attack of 5.13◦ , with the transition location moving upstream rapidly on the upper surface and the lower

surface separation bubble fixing transition near mid-chord over the range of angles of attack. Numerical results were

obtained using a Mach number of 0.10, Reynolds number of 2.0 × 106 , and a free-stream turbulence intensity of 0.07%,

as recommended by the TCMPS guidelines [81].

                        Computational results were obtained using the S809 grid families supplied through the TCMPS, which were

generated using the same conformal mapping method as for the NLF0416 [81], with the coarse, medium, and fine grid

dimensions presented in Table 2. A grid-refinement study was conducted for the S809 airfoil at 6◦ angle of attack

and is plotted in Figure 8. This condition was chosen as it represents the point beyond the upper corner of the drag

<!-- ===== PDF page 24 ===== -->

                            Table 2     TMPCS S809 structured C-grid dimensions [81]

                           grid level    dimensions    ∆s (chord)     average/maximum y+
                            coarse        529 × 73     8.90 × 10−6          0.60/1.73
                           medium         705 × 97     6.47 × 10−6          0.52/0.81
                              fine       1057 × 145    4.19 × 10−6          0.28/0.54

bucket where the transition location is rapidly moving, which is is expected to amplify the effects of grid refinement on

the predictions [81]. Figure 8 also illustrates the forces and transition locations obtained using the SA-LM2015 and

SA-sLM2015 transition models on the fine grid level. Although the transition locations match experiment well, both

models under-predict drag before and after the upper corner of the drag bucket. Again, this trend has been identified by

previous studies using variations of the γ- Re
                                             ˜ θt model [48, 82, 83], as well as studies using the AFT transition model

[35, 82], and does not appear to be caused by inaccurate transition prediction.

   The residual convergence histories for these cases are presented in Figure 9 for computations using 8 processors.

Similar to the NLF0416 airfoil, the large separation bubbles at low angles of attack led to stiff linear systems which

slowed convergence. For all angles of attack the SA-sLM2015 model displays substantially improved convergence. In

particular, for the 9◦ angle of attack case the SA-LM2015 model was unable to converge due to the unsteady interaction

of the transition model with the laminar separation bubble on the upper surface of the airfoil. This behaviour was also

identified by Denison and Pulliam using the LM2015 model [82]. Limiting vorticity in the transition model source

terms prevented this oscillation without affecting the accuracy of the solution.

C. NASA NLF2-0415 Infinite Swept Wing

   The NASA NLF2-0415 Infinite Swept Wing consists of an NLF2-0415 airfoil extruded with a 45◦ sweep angle.

Experimental measurements were obtained at an angle of attack of −4◦ , over a range of Mach and Reynolds numbers [85],

and including a variety of surface roughness levels [86]. Transition locations were determined using the naphthalene

visualization technique. Transition locations for a roughness level, h, of 3.3µm, representing a painted surface, were

obtained using the SA-LM2015 and SA-sLM2015 models and are compared with experimental values [86] in Figure 10.

The Reynolds number based on free-stream velocity magnitude and chord length was varied from 1.8 × 106 to 3.5 × 106 .

The free-stream turbulence intensity for the wind tunnel was measured as 0.20%. Transition locations are presented

only for the upper surface of the wing, as experimental results are not available for the lower surface.

   The computational grid used to obtain these results consists of 270,336 nodes split into 128 blocks, with 90,112

nodes per cross-section, and an off-wall spacing of 2.5 × 10−6 chords, resulting in average and maximum y + values of

0.21 and 0.46, respectively. The analytical full Jacobian was used in the inexact-Newton phase for this case and the Tu

Braunschweig Sickle Wing. Periodic boundary conditions were used to simulate the infinite swept wing. The results in

<!-- ===== PDF page 25 ===== -->

                                                       SA-LM2015fine                                             SA-LM2015cl
         0.012                                         SA-LM2015medium          .0094                            SA-LM2015cd
                                                       SA-LM2015coarse                                           SA-sLM2015cl
          0.01                                         SA-sLM2015fine                                            SA-sLM2015cd                .8640
                                                                                .0092
                                                       SA-sLM2015medium
         0.008                                         SA-sLM2015coarse
                                                                                .0090
                       upper                                                                                                                 .8600
         0.006
                                                                                .0088
        Cf

                                                                                                                                                 Cl
                                                                            Cd
         0.004
                                                                                .0086                                                        .8560
         0.002
                                                                                .0084
                                                lower
              0                                                                                                                              .8520
                                                                                .0082
        -0.002
                   0       0.2     0.4         0.6        0.8      1                            1E-05       1.5E-05          2E-05     2.5E-05
                                         x/c                                                                     N-1
                       (a) Skin friction grid convergence                               (b) Grid convergence for lift and drag

                                                     SA-LM2015fine                                        SA-LM2015fine
                                                     SA-sLM2015fine                                       SA-sLM2015fine
                                                     TU Delft Experiment                                  TU Delft Experiment
             1.2                                                                  1.2

             0.8
                                                                                  0.8
        Cl

                                                                            Cl

             0.6
                                                                                  0.6
             0.4

             0.2                                                                  0.4

              0                                                                   0.2
                                 .0100                   .0150                          0          2         4         6         8      10
                                         Cd                                                                      aoa

                                 (c) Drag polar                                                           (d) Lift curve

                                                     SA-LM2015fine                                                           SA-LM2015fine
                                                     SA-sLM2015fine               12                                         SA-sLM2015fine
                                                     TU Delft Experiment                                                     TU Delft Experiment
         .0160                                                                    10

                                                                                                        upper                  lower
         .0140                                                                     8
                                                                            aoa

         .0120
        Cd

         .0100

         .0080
         .0060
                           2       4           6           8                                0       0.2       0.4      0.6       0.8    1
                                         aoa                                                      transition location (x/c)
                                 (e) Drag curve                                                    (f) Transition locations

Fig. 8 S809 grid convergence at 6◦ angle of attack, and angle of attack sweep simulated using the SA-LM2015
and SA-sLM2015 transition models at M = 0.1, Re = 2.0 × 106 , and Tu = 0.07%. The ‘fine’ grid level was used
for the sweep, with the experimental results obtained at the TU Delft wind tunnel [84].

<!-- ===== PDF page 26 ===== -->

                         101                                                                 101                                                                 101
                                                       SA-LM2015                                                           SA-LM2015                                                          SA-LM2015
                        10−1                           SA-sLM2015                           10−1                           SA-sLM2015                           10−1                          SA-sLM2015
Relative Residual, d

                                                                    Relative Residual, d

                                                                                                                                        Relative Residual, d
                        10−3                                                                10−3                                                                10−3
                        10−5                                                                10−5                                                                10−5
                        10−7                                                                10−7                                                                10−7
                        10−9                                                                10−9                                                                10−9
                   10−11                                                               10−11                                                               10−11
                   10−13                                                               10−13                                                               10−13
                               0   3000   6000   9000 12000                                        0   3000   6000   9000 12000                                        0   3000   6000   9000 12000
                                           Time(s)                                                              Time(s)                                                             Time(s)
                                     (a) aoa = 1◦                                                        (b) aoa = 2◦                                                        (c) aoa = 3◦
                         101                                                                 101                                                                 101
                                                      SA-LM2015                                                           SA-LM2015                                                           SA-LM2015
                        10−1                          SA-sLM2015                            10−1                          SA-sLM2015                            10−1                          SA-sLM2015
Relative Residual, d

                                                                    Relative Residual, d

                                                                                                                                        Relative Residual, d
                        10−3                                                                10−3                                                                10−3
                        10−5                                                                10−5                                                                10−5
                        10−7                                                                10−7                                                                10−7
                        10−9                                                                10−9                                                                10−9
                   10−11                                                               10−11                                                               10−11
                   10−13                                                               10−13                                                               10−13
                               0   3000 6000 9000 12000                                            0   3000   6000   9000 12000                                        0   3000   6000   9000 12000
                                           Time(s)                                                             Time(s)                                                              Time(s)
                                     (d) aoa = 4◦                                                        (e) aoa = 5◦                                                        (f) aoa = 6◦
                         101                                                                 101                                                                 101
                                                       SA-LM2015                                                           SA-LM2015                                                          SA-LM2015
                        10−1                           SA-sLM2015                           10−1                           SA-sLM2015                           10−1                          SA-sLM2015
Relative Residual, d

                                                                    Relative Residual, d

                                                                                                                                        Relative Residual, d
                        10−3                                                                10−3                                                                10−3
                        10−5                                                                10−5                                                                10−5
                        10−7                                                                10−7                                                                10−7
                        10−9                                                                10−9                                                                10−9
                   10−11                                                               10−11                                                               10−11
                   10−13                                                               10−13                                                               10−13
                               0   3000   6000   9000 12000                                        0   3000   6000   9000 12000                                        0   3000   6000   9000 12000
                                            Time(s)                                                             Time(s)                                                            Time(s)
                                     (g) aoa = 7◦                                                        (h) aoa = 8◦                                                        (i) aoa = 9◦

  Fig. 9 S809 residual convergence histories for the SA-LM2015 and SA-sLM2015 transition models on the ‘fine’
  grid level split into 8 computational blocks. The airfoil was simulated at M = 0.1, Re = 2.0 × 106 , Tu = 0.07%
  and over a range of angles of attack.

  Figure 10 demonstrate that both SA-LM2015 and SA-sLM2015 transition models predict transition accurately, while

  plots of the residual convergence histories in Figure 11 demonstrate that the SA-sLM2015 model provides significantly

  improved convergence to steady state.

                        The results demonstrate that the SA-LM2015 model is able to closely predict the location of transition for all

  Reynolds numbers except for 2 × 106 , where there is an error of approximately 14% between the experimental and

  predicted values. This error is also present in the work of Jung and Baeder [48], where they coupled a semi-local

  correlation-based transition model formulation with the SA model using the crossflow correlations developed by Muller

  and Herbst [42], and in the original model developed by Langtry et al. [43]. In the work produced by Radeztsky et al.

  [86], the results demonstrate that there is significant experimental error associated with this case, which may help to

  explain this discrepancy.

<!-- ===== PDF page 27 ===== -->

                                                                                                    SA-LM2015                                                                                        SA-LM2015
                                                                                                    SA-sLM2015            0.006                                                                      SA-sLM2015
                                                                                                    Experiment
                              0.8
  transition location (x/c)

                              0.6                                                                                         0.004

                                                                                                                          Cf
                              0.4
                                                                                                                          0.002

                                           2E+06                        3E+06                                                     0                                           0.5                              1
                                                    Reynolds number                                                                                                           x/c
                                                    (a) Transition locations                                                                   (b) Skin friction profiles

Fig. 10 Experimental and computational results for the NLF2-0415 infinite swept wing [86] simulated using
the SA-LM2015 and SA-sLM2015 transition models at M = 0.15, α = −4◦ , Tu = 0.20%, and over a range of
Reynolds numbers from 1.8 − 3.5 × 106 .

                                101                                                              101                                                                101
                                                               SA-LM2015                                                          SA-LM2015                                                          SA-LM2015
                               10−1                                                             10−1                                                               10−1
                                                               SA-sLM2015                                                         SA-sLM2015                                                         SA-sLM2015
Relative Residual, d

                                                                            Relative Residual, d

                                                                                                                                               Relative Residual, d
                               10−3                                                             10−3                                                               10−3
                               10−5                                                             10−5                                                               10−5
                               10−7                                                             10−7                                                               10−7
                               10−9                                                             10−9                                                               10−9
                              10−11                                                            10−11                                                              10−11
                              10−13                                                            10−13                                                              10−13
                              10−15                                                            10−15                                                              10−15
                                      0      1000       2000     3000                                  0      1000        2000        3000                                0         1000     2000       3000
                                                     Time(s)                                                          Time(s)                                                              Time(s)
                                          (a) Re = 1.8 × 106                                               (b) Re = 2.0 × 106                                                  (c) Re = 2.2 × 106
                                101                                                              101                                                                101
                                                               SA-LM2015                                                          SA-LM2015                                                           SA-LM2015
                               10−1                                                             10−1                                                               10−1
                                                               SA-sLM2015                                                         SA-sLM2015                                                          SA-sLM2015
Relative Residual, d

                                                                            Relative Residual, d

                                                                                                                                               Relative Residual, d

                               10−3                                                             10−3                                                               10−3
                               10−5                                                             10−5                                                               10−5
                               10−7                                                             10−7                                                               10−7
                               10−9                                                             10−9                                                               10−9
                              10−11                                                            10−11                                                              10−11
                              10−13                                                            10−13                                                              10−13
                              10−15                                                            10−15                                                              10−15
                                      0      1000       2000     3000                                  0      1000        2000        3000                                0         1000      2000       3000
                                                     Time(s)                                                          Time(s)                                                              Time(s)
                                          (d) Re = 2.5 × 106                                               (e) Re = 3.0 × 106                                                  (f) Re = 3.5 × 106

Fig. 11 NLF2-0415 residual convergence histories for the SA-LM2015 and SA-sLM2015 transition models,
with the grid split into 128 computational blocks. The infinite swept wing was simulated at M = 0.15, α = −4◦ ,
Tu = 0.20%, and over a range of Reynolds numbers.

                              It is important to note that the choice of numerical dissipation strongly affects the accuracy of the crossflow transition

model. Both the scalar dissipation model developed by Jameson et al. [87] and later refined by Pulliam [88], and

<!-- ===== PDF page 28 ===== -->

                              Table 3    Sickle Wing structured O-O grid dimensions

            grid level   total node count, N    N-span     N-chord       ∆s (chord)    average/maximum y+
             coarse           2.32 × 106         171         111         7.81 × 10−6         0.70/2.64
            medium            5.19 × 106         217         139         6.02 × 10−6         0.54/1.55
               fine          10.19 × 106         273         175         4.29 × 10−6         0.39/1.13

the matrix-based dissipation model developed by Swanson and Turkel [79] were investigated for this case. On the

mesh used, the scalar dissipation model drastically under predicts the crossflow strength and helicity values. This was

determined to be due to the increased dissipation of the scalar dissipation model in and above the laminar boundary

layer, which affected the vorticity profile within the boundary layer.

D. TU Braunschweig Sickle Wing

   The TU Braunschweig Sickle Wing [89] was developed as a test case for crossflow transition on a complex geometry.

The model consists of five main sections designed to create three-dimensional boundary layers with increasing crossflow

towards the wing tip. The first section raises the model away from the turbulent boundary layer of the wind tunnel wall,

with the three main wing sections swept at 30, 45, and 55◦ with distinct kinks in the leading edge sweep. As a result,

large span-wise gradients are present and the assumptions of linear local stability theory are challenged [89]. The wing

is capped with a rounded wing tip. At the experimental flow conditions, the upper surface of the wing has extended

regions of strong favourable pressure gradients with crossflow instabilities developing at the section kinks. This results

in mixed mode transition scenarios with Tollmien-Schlichting and crossflow instabilities triggering transition, except on

the upper surface near the root of the wing, where a laminar separation bubble triggers transition at approximately 80%

chord.

   The Sickle Wing was simulated at a Reynolds number of 2.75 × 106 , Mach number of 0.16, and −2.6◦ angle of attack.

The rms surface roughness of the experimental model was measured to be 1.47µm with a five-point average of the

maximum peak to peak height of 9.78µm [89]. In order to be consistent with results obtained using the LM2015 model

[43, 44] the surface roughness for the simulations was set to the peak to peak height value of 9.78µm with a free-stream

turbulence intensity of 0.17%, which is consistent with the turbulence intensity measured in the wind tunnel [89].

   Three levels of O-O topology structured multi-block grids were simulated with the grid details presented in Table 3,

including the number of spanwise and chordwise nodes on each of the upper and lower surfaces. Each grid level was

split into 592 computational blocks. Figure 12 illustrates the skin friction profiles on the upper and lower surfaces of

the Sickle Wing using both the SA-LM2015 and SA-sLM2015 transition models for all three grid levels overlaid with

experimental transition locations [89]. Both models predict the transition fronts reasonably well.

   As the grid is refined the transition front on the lower surface, which is dominated by natural transition, moves

<!-- ===== PDF page 29 ===== -->

                coarse                                 medium                                     fine

                                                                                                                 Cf
                                                                                                                0.006
                                                                                                                0.005
                                                                                                                0.004
                                                                                                                0.003
                                                                                                                0.002
                                                                                                                0.001

     (a) SA-LM2015 upper surface             (b) SA-LM2015 upper surface             (c) SA-LM2015 upper surface

                                                                                                                 Cf
                                                                                                                0.006
                                                                                                                0.005
                                                                                                                0.004
                                                                                                                0.003
                                                                                                                0.002
                                                                                                                0.001

    (d) SA-sLM2015 upper surface             (e) SA-sLM2015 upper surface            (f) SA-sLM2015 upper surface

                                                                                                                 Cf
                                                                                                                0.006
                                                                                                                0.005
                                                                                                                0.004
                                                                                                                0.003
                                                                                                                0.002
                                                                                                                0.001

     (g) SA-LM2015 lower surface             (h) SA-LM2015 lower surface              (i) SA-LM2015 lower surface

                                                                                                                 Cf
                                                                                                                0.006
                                                                                                                0.005
                                                                                                                0.004
                                                                                                                0.003
                                                                                                                0.002
                                                                                                                0.001

     (j) SA-sLM2015 lower surface            (k) SA-sLM2015 lower surface            (l) SA-sLM2015 lower surface

Fig. 12 Skin friction profiles on the upper and lower surfaces of the Sickle Wing for all grid levels overlaid with
experimental transition locations [89]. Results were obtained using the SA-LM2015 and SA-sLM2015 transition
models at M = 0.16, α = −2.6◦ , Tu = 0.17%, h = 9.78µm, and a Reynolds number of 2.75 × 106 .

slightly downstream and agrees well with experiment. The transition front on the upper surface, where crossflow-induced

transition is the primary mechanism, appears more sensitive to the grid refinement. The transition locations at the kink

regions move upstream as the stationary crossflow vortices are better resolved. The lift and drag grid convergence for

both models is presented in Figure 13(b).

   The pressure coefficient profiles are extracted at the midspan of each sweep section for the fine grid, corresponding

to spanwise sections at 25%, 55%, and 85% span, and are compared with the experimental profiles presented by Kruse

<!-- ===== PDF page 30 ===== -->

                        lower                                                                SA-LM2015cl
                                                                                             SA-LM2015cd
     -0.4                                                           .0063                    SA-sLM2015cl
                                                 upper                                       SA-sLM2015cd
                                                                                                                     -.0640
     -0.2

                                                                                                                      CL
                                                                    .0062
         0                SA-LM201525%

                                                                    CD
   Cp

                          SA-LM201555%
                                                                                                                     -.0650
        0.2               SA-LM201585%
                          SA-sLM201525%                             .0061
                          SA-sLM201555%
        0.4               SA-sLM201585%
                          Experiment25%
                          Experiment55%                             .0060                                            -.0660
        0.6
                          Experiment85%
              0   0.2      0.4         0.6      0.8      1
                                                                                 3E-05       4E-05        5E-05
                                 x/c
                                                                                             N-2/3
  (a) Pressure coefficient profiles obtained on the fine grid at
                                                                            (b) Grid convergence for lift and drag
  three spanwise sections compared with experiment [90]

Fig. 13 Experimental and fine grid pressure coefficient profiles, and grid convergence for the SA-LM2015 and
SA-sLM2015 transition models. The Sickle Wing was simulated at M = 0.16, α = −2.6◦ , Tu = 0.17%, and a
Reynolds number of 2.75 × 106 [90].

et al. [90] in Figure 13(a). The pressure coefficient profiles from the simulations lie on top of one another and the results

demonstrate the ‘wiggles’ observed in experiment at approximately x/c = 0.35, corresponding to natural-transition on

the lower surface of the Sickle Wing, are accurately reproduced[90]. There is a slight deviation between the simulations

and experiment at 25% span on the upper surface of the Sickle Wing where both transition models, SA-LM2015 and

SA-sLM2015, predict crossflow-induced transition upstream of the transition location in the experiment, which is

triggered at x/c = 0.80 by a laminar separation bubble [89, 90]. This deviation is also illustrated in Figures 12(a-f) and

is consistent with results obtained using the LM2015 transition model [44].

   The relative residual and drag coefficient convergence histories are presented in Figure 14. The flow complexity

in this case slows convergence for both models, however the smooth-variant displays improved residual convergence

behaviour in Figures 14(a-c). The drag convergence, the magnitude of the difference between the current drag coefficient

and the converged drag coefficient, plotted in Figures 14(d-f) show that the SA-sLM2015 model leads to faster

convergence of the drag coefficient.

                                             VI. Conclusions and Future Work
   A novel coupling of the γ- Re
                               ˜ θt LM2015 transition model with helicity-based crossflow correlations [43] and the

one-equation SA turbulence model has been developed and designated SA-LM2015. The modifications maintain the

fully local formulation pioneered by Langtry and Menter [29]. Smooth approximations to non-differentiable and stiff

functions in the SA-LM2015 model are introduced in order to eliminate discontinuities and improve numerical behaviour,

<!-- ===== PDF page 31 ===== -->

                                           coarse                                                                  medium                                                                 fine
                           101                                                                      101                                                                  101
                                                            SA-LM2015                                                            SA-LM2015                                                            SA-LM2015
                          10−1                                                                     10−1                                                                 10−1
                                                            SA-sLM2015                                                           SA-sLM2015                                                           SA-sLM2015
Relative Residual, d

                                                                         Relative Residual, d

                                                                                                                                              Relative Residual, d
                          10−3                                                                     10−3                                                                 10−3
                          10−5                                                                     10−5                                                                 10−5
                          10−7                                                                     10−7                                                                 10−7
                          10−9                                                                     10−9                                                                 10−9
                         10−11                                                                    10−11                                                                10−11
                         10−13                                                                    10−13                                                                10−13
                         10−15                                                                    10−15                                                                10−15
                                    0   2000         4000    6000                                            0   5000 10000 15000 20000                                           0   20000 40000 60000 80000
                                                Time(s)                                                                Time(s)                                                              Time(s)
                                               (a)                                                                    (b)                                                                  (c)
                              104                           SA-LM2015                                  104                       SA-LM2015                                  104                       SA-LM2015
  Drag Convergence(counts)

                                                                           Drag Convergence(counts)

                                                                                                                                                Drag Convergence(counts)
                              102                           SA-sLM2015                                 102                       SA-sLM2015                                 102                       SA-sLM2015
                              100                                                                      100                                                                  100
                             10−2                                                                     10−2                                                                 10−2
                             10−4                                                                     10−4                                                                 10−4
                             10−6                                                                     10−6                                                                 10−6
                             10−8                                                                     10−8                                                                 10−8
                        10−10                                                                    10−10                                                                10−10
                                    0   2000     4000        6000                                            0   5000 10000 15000 20000                                           0   20000 40000 60000 80000
                                                Time(s)                                                                Time(s)                                                              Time(s)
                                               (d)                                                                    (e)                                                                  (f)

Fig. 14 Relative residual and drag coefficient convergence histories for the SA-LM2015 and SA-sLM2015
transition models on all grid levels, with the grids split into 592 computational blocks. The Sickle Wing was
simulated at M = 0.16, α = −2.6◦ , Tu = 0.17%, h = 9.78µm, and a Reynolds number of 2.75 × 106 .

with the smooth model designated SA-sLM2015. Deep convergence is achieved through the use of a fully-coupled,

fully-implicit Newton-Krylov algorithm with modifications to the solution algorithm, including a source-term time

stepping strategy, introduced to address the large source terms in the turbulence and transition model equations. Several

of the strategies presented, including the smoothing techniques and source term time stepping, could be useful in the

context of other transition models, specifically transition models which use similar intermittency transport equations,

such as the following [35, 39, 40, 48, 91].

                             Simulations of the S809 and NLF0416 natural-laminar-flow airfoil test cases demonstrate that both models,

SA-LM2015 and SA-sLM2015, are able to predict natural and separation-induced transition accurately, with the

SA-sLM2015 model displaying improved numerical behaviour. Specifically, the new smooth formulations improve

convergence, and combined with the new source term formulation, including limiting vorticity magnitude, oscillations

near laminar separation bubbles produced by the LM2015 and SA-LM2015 models are eliminated. Three-dimensional

simulations of the NLF2-0415 infinite swept wing and Sickle Wing test cases demonstrate that the model is able to

capture the effects of crossflow accurately.

                             The SA-LM2015 and SA-sLM2015 transition models are validated in the current work using popular transition

test cases under subsonic conditions. Future work will involve applying both transition models, SA-LM2015 and

SA-sLM2015, to transonic transition test cases in order to confirm that the trends produced in subsonic flow extend to

<!-- ===== PDF page 32 ===== -->

higher Mach numbers.

                                                     Acknowledgments
    This work was partially funded by the NASA Transformational Tools and Technologies (TTT) and Advanced Air

Transport Technology (AATT) projects, the Natural Sciences and Engineering Research Council (NSERC), and the

University of Toronto. Computations were performed on the Niagara supercomputer at the SciNet HPC Consortium, a

part of Compute Canada. The authors gratefully acknowledge the input and support provided by Dr. Masayuki Yano, Dr.

Thomas Reist, Dr. Gaetan Kenway, Dr. Thomas Pulliam, Dr. James Coder, and Dr. Cetin Kiris.

                                                          References
 [1] Green, J. E., “Civil aviation and the environment: the next frontier for the aerodynamicist,” Aeronautical Journal, Vol. 110, No.

     1110, 2006, pp. 469–486. doi:https://doi.org/10.1017/S0001924000001378.

 [2] Bushnell, D. M., “Overview of aircraft drag reduction technology,” AGARD Report 786 (Special Course on Skin Friction Drag

     Reduction), 1992.

 [3] Malik, M. R., Crouch, J. D., Saric, W. S., Lin, J. C., and Whalen, E. A., “Application of Drag Reduction Techniques to Transport

     Aircraft,” Encyclopedia of Aerospace Engineering, 2015. doi:https://doi.org/10.1002/9780470686652.eae1013.

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

[11] Reed, H. L., and Saric, W. S., “Transition mechanisms for transport aircraft,” 38th Fluid Dynamics Conference and Exhibit,

     AIAA Paper 2008-3743, Seattle, Washington, 2008. doi:https://doi.org/10.2514/6.2008-3743.

[12] Arnal, D., Casalis, G., and Houdeville, R., “Practical Transition Prediction Methods: Subsonic and Transonic Flows,” VKI

<!-- ===== PDF page 33 ===== -->

     Lecture Series: Advances in Laminar-Turbulent Transition Modeling, 2009.

[13] Aupoix, B., Arnal, D., Bezard, H., Chaouat, B., and Chedevergne, F., “Transition and Turbulence Modeling,” AerospaceLab,

     2011, pp. 1–13.

[14] Pasquale, D., Rona, A., and Garrett, S. J., “A selective review of CFD transition models,” 39th AIAA Fluid Dynamics Conference,

     AIAA Paper 2009-3812, San Antonio, Texas, 2009. doi:https://doi.org/10.2514/6.2009-3812.

[15] Savill, A. M., “By-pass transition using conventional closures,” Closure strategies for turbulent and transitional flows,

     Cambridge University Press, Cambridge, Vol. 17, 2002, pp. 464–492. doi:https://doi.org/10.1017/CBO9780511755385.019.

[16] Halila, G. L. O., Chen, G., Shi, Y., Fidkowksi, K., Martins, J. R. R. A., and Mendonca, M. T., “High-Reynolds number

     transitional flow simulation via parabolized stability equations with an adaptive RANS solver,” Aerospace Science and

     Technology, Vol. 91, 2019, pp. 321–336. doi:https://doi.org/10.1016/j.ast.2019.05.018.

[17] Drela, M., and Giles, M. B., “Viscous-inviscid analysis of transonic and low Reynolds number airfoils,” AIAA Journal, Vol. 25,

     No. 10, 1987, pp. 1347–1355. doi:https://doi.org/10.2514/3.9789.

[18] Rashad, R., and Zingg, D. W., “Aerodynamic shape optimization for natural laminar flow using a discrete-adjoint approach,”

     AIAA Journal, Vol. 54, No. 11, 2016, pp. 3321–3337. doi:https://doi.org/10.2514/1.J054940.

[19] Mayda, E., “Boundary Layer Transition Prediction for Reynolds-Averaged Navier-Stokes Methods,” Ph.D. thesis, University of

     California, Davis, 2007.

[20] Stock, H. W., “e N Transition Prediction in Three-Dimensional Boundary Layers on Inclined Prolate Spheroids,” AIAA Journal,

     Vol. 44, No. 1, 2006, pp. 108–118. doi:https://doi.org/10.2514/1.16026.

[21] Krimmelbein, N., and Radespiel, R., “Transition prediction for three-dimensional flows using parallel computation,” Computers

     and Fluids, Vol. 38, No. 1, 2009, pp. 121–136. doi:https://doi.org/10.1016/j.compfluid.2008.01.004.

[22] Arnal, D., Habiballah, M., and Coutols, E., “Laminar instability theory and transition criteria in two and three-dimensional

     flow,” La Recherche Aérospatiale, , No. 2, 1984, pp. 125–143.

[23] Krumbein, A., Krimmelbein, N., and Schrauf, G., “Automatic transition prediction in hybrid flow solver, part 1: Methodology

     and sensitivities,” Journal of Aircraft, Vol. 46, No. 4, 2009, pp. 1176–1190. doi:https://doi.org/10.2514/1.39736.

[24] Krumbein, A., Krimmelbein, N., and Schrauf, G., “Automatic transition prediction in hybrid flow solver, part 2: Practical

     application,” Journal of Aircraft, Vol. 46, No. 4, 2009, pp. 1191–1199. doi:https://doi.org/10.2514/1.39738.

[25] Moens, F., Perraud, J., Krumbein, A., Toulorge, T., Iannelli, P., and Hanifi, A., “Transition prediction and im-

     pact on a three-dimensional high-lift-wing configuration,” Journal of Aircraft, Vol. 45, No. 5, 2008, pp. 1751–1766.

     doi:https://doi.org/10.2514/1.36238.

[26] Streit, T., Horstmann, H., Shraut, G., Hein, S., Fey, U., Egmai, Y., Perraud, J., Salah El Din, I., Cella, U., and Quest, J.,

     “Complementary numerical and experimental data analysis of the ETW Telfona Pathfinder wing transition tests,” 49th Aerospace

     Sciences Meeting and Exhibit, AIAA Paper 2011-881, Orlando, Florida, January, 2011. doi:https://doi.org/10.2514/6.2011-881.

[27] Shi, Y., Gross, R., Mader, C. A., and Martins, J., “Transition Prediction Based on Linear Stability Theory with the RANS

     Solver for Three-Dimensional Configurations,” 2018 AIAA Aerospace Sciences Meeting, AIAA Paper 2018-0819, 2018.

     doi:https://doi.org/10.2514/6.2018-0819.

<!-- ===== PDF page 34 ===== -->

[28] Seyfert, C., and Krumbein, A., “Evaluation of a Correlation-Based Transition Model and Comparison with the e N Method,”

     Journal of Aircraft, Vol. 49, No. 6, 2012, pp. 1765–1773. doi:https://doi.org/10.2514/1.C031448.

[29] Langtry, R. B., and Menter, F. R., “Correlation-based transition modeling for unstructured parallelized computational fluid

     dynamics codes,” AIAA Journal, Vol. 47, No. 12, 2009, pp. 2894–2906. doi:https://doi.org/10.2514/1.42362.

[30] Coder, J. G., and Maughmer, M. D., “A CFD-Compatible Transition Model using an Amplification Factor Transport Equation,”

     51st AIAA Aerospace Sciences Meeting including the New Horizons Forum and Aerospace Exposition, AIAA Paper 2013-0253,

     2013. doi:https://doi.org/10.2514/6.2013-253.

[31] Coder, J. G., and Maughmer, M. D., “Computational fluid dynamics compatible transition modeling using an amplification

     factor transport equation,” AIAA Journal, Vol. 52, No. 11, 2014, pp. 2506–2512. doi:https://doi.org/10.2514/1.J052905.

[32] Coder, J. G., “Enhancement of the Amplification Factor Transport Transition Modeling Framework,” 55th AIAA Aerospace

     Sciences Meeting, AIAA Paper 2017-1709, 2017. doi:https://doi.org/10.2514/6.2017-1709.

[33] Stefanski, D., Glasby, R., Erwin, J. T., and Coder, J. G., “Development of a Predictive Capability for Laminar-Turbulent

     Transition in HPCMP CREATETM-AV Kestrel Component COFFE using the Amplification Factor Transport Model,” 2018

     AIAA Aerospace Sciences Meeting, AIAA Paper 2018-1041, 2018. doi:https://doi.org/10.2514/6.2018-1041.

[34] Coder, J. G., Pulliam, T. H., and Jensen, J. C., “Contributions to HiLiftPW-3 Using Structured, Overset Grid Methods,” 2018

     AIAA Aerospace Sciences Meeting, AIAA Paper 2018-1039, 2018. doi:https://doi.org/10.2514/6.2018-1039.

[35] Coder, J. G., “Further Development of the Amplification Factor Transport Transition Model for Aerodynamic Flows,” 2019

     AIAA Aerospace Sciences Meeting, AIAA Paper 2019-0039, 2019. doi:https://doi.org/10.2514/6.2019-0039.

[36] Spalart, P. R., and Allmaras, S. R., “A one equation turbulence model for aerodynamic flows,” 30th AIAA Aerospace Sciences

     Meeting and Exhibit, AIAA Paper 092-0439, Reno, Nevada, United States 1992. doi:https://doi.org/10.2514/6.1992-439.

[37] Xu, J., Han, X., Qiao, L., Bai, J., and Zhang, Y., “Fully Local Amplification Factor Transport Equation for Stationary Crossflow

     Instabilities,” AIAA Journal, Vol. 57, No. 7, 2019, pp. 2682–2693. doi:https://doi.org/10.2514/1.J057502.

[38] Choi, J. H., and Kwon, O. J., “Enhancement of a Correlation-Based Transition Turbulence Model for Simulating Crossflow

     Instability,” AIAA Journal, Vol. 53, No. 10, 2015, pp. 3063–3072. doi:https://doi.org/10.2514/1.J053887.

[39] Grabe, C., and Krumbein, A., “Correlation-Based Transition Modeling for Three-Dimensional Aerodynamic Configurations,”

     Journal of Aircraft, Vol. 50, No. 5, 2013, pp. 1533–1539. doi:https://doi.org/10.2514/1.C032063.

[40] Grabe, C., Shengyang, N., and Krumbein, A., “Transition Transport Modeling for the Prediction of Crossflow

     Transition,” 34th AIAA Applied Aerodynamics Conference, AIAA Paper 2016-3572, Washington D.C., June, 2016.

     doi:https://doi.org/10.2514/6.2016-3572.

[41] Xu, J. K., Bai, J. Q., Qiao, L., and Zhang, Y., “Correlation-Based Transition Transport Modeling for Sim-

     ulating Crossflow Instabilities,” Journal of Applied Fluid Mechanics, Vol. 9, No. 5, 2016, pp. 2435–2442.

     doi:https://doi.org/10.18869/acadpub.jafm.68.236.25356.

[42] Muller, C., and Herbst, F., “Modelling of Crossflow-Induced Transition Based on Local Variables,” 6th European Conference

     on Computational Fluid Dynamics (ECFD), 2014.

[43] Langtry, R. B., Sengupta, K., Yeh, D. T., and Dorgan, A. J., “Extending the γ- Re
                                                                                     ˜ θt Correlation based Transition Model for

<!-- ===== PDF page 35 ===== -->

     Crossflow Effects,” 45th AIAA Fluid Dynamics Conference, AIAA Paper 2015-2474, 2015. doi:https://doi.org/10.2514/6.2015-

     2474.

[44] Venkatachari, B. S., Paredes, P., Derlaga, J. M., Buning, P., Choudhari, M., Li, F., and Chang, C. L., “Assessment of Transition

     Modeling Capability in OVERFLOW with Emphasis on Swept-Wing Configurations,” 2020 AIAA Aerospace Sciences Meeting,

     AIAA Paper 2020-1034, 2020. doi:https://doi.org/10.2514/6.2020-1034.

[45] Menter, F. R., “Two-Equation Eddy-Viscosity Turbulence Models for Engineering Applications,” AIAA Journal, Vol. 32, No. 8,

     1994, pp. 1598–1605. doi:https://doi.org/10.2514/3.12149.

[46] Nie, S., Krimmelbein, N., Krumbein, A., and Grabe, C., “Extension of a Reynolds-stress-based Transition Transport Model for

     Crossflow Transition,” Journal of Aircraft, Vol. 55, No. 4, 2018, pp. 1641–1654. doi:https://doi.org/10.2514/1.C034586.

[47] Medida, S., and Baeder, J. D., “Application of the Correlation-based γ- Re
                                                                              ˜ θt Transition Model to the Spalart-Allmaras

     Turbulence Model,” 20th AIAA Computational Fluid Dynamics Conference, AIAA Paper 2011-3979, Honolulu, HI, 2011.

     doi:https://doi.org/10.2514/6.2011-3979.
                                          ˜ θt -SA with Crossflow Transition Model using Hamiltonian-Strand Approach,” 2018 AIAA
[48] Jung, Y. S., and Baeder, J. D., “γ- Re

     Aerospace Sciences Meeting, AIAA Paper 2018-1040, 2018. doi:https://doi.org/10.2514/6.2018-1040.

[49] Carnes, J. A., and Coder, J. G., “Effect of Crossflow on the S-76 and PSP Rotors in Hover,” 2020 AIAA Aerospace Sciences

     Meeting, AIAA Paper 2020-0773, 2020. doi:https://doi.org/10.2514/6.2020-0773.

[50] Schucker, J., “Development of a Three-Equation γ- Re
                                                        ˜ θt -Spalart-Allmaras Turbulence-Transition Model,” Master’s thesis,

     DLR, 2012.

[51] Langtry, R. B., “A correlation-based transition model using local variables for unstructured parallelized CFD codes,” Ph.D.

     thesis, DLR, 2006.

[52] Pulliam, T. H., and Zingg, D. W., Fundamental algorithms in computational fluid dynamics, Springer International Publishing,

     2014.

[53] Suluksna, K., Dechaumphai, P., and Juntasaro, E., “Correlations for modeling transitional boundary layers under influences of

     freestream turbulence and pressure gradient,” International Journal of Heat and Fluid Flow, Vol. 30, No. 1, 2009, pp. 66–75.

     doi:https://doi.org/10.1016/j.ijheatfluidflow.2008.09.004.

[54] Allmaras, S. R., and Johnson, F. T., “Modifications and clarifications for the implementation of the Spalart-Allmaras turbulence

     model,” Seventh International Conference on Computational Fluid Dynamics, ICCFD7 Paper 1902, 2012, pp. 1–11.

[55] Osusky, M., “A Parallel Newton-Krylov-Schur Algorithm for the Reynolds-Averaged Navier-Stokes Equations,” Ph.D. thesis,

     Graduate Department of Aerospace Science and Engineering, University of Toronto, 2013.

[56] Bertsekas, D. P., “Nondifferentiable optimization via approximation,” Mathematical Programming Study 3, 1975, pp. 1–25.

     doi:https://doi.org/10.1007/BFb0120696.

[57] Chakraborty, A., Sinha Roy, A., and Dasgupta, B., “Non-parametric Smoothing for Gradient Methods in

     Non-differentiable Optimization Problems,” 2016 IEEE International Conference on Systems, Man, and Cybernetics, 2016.

     doi:https://doi.org/10.1109/SMC.2016.7844819.

[58] Rivest, R. L., “Game Tree Searching by Min/Max Approximation,” Artificial Intelligence, Vol. 34, No. 1, 1988, pp. 77–96.

<!-- ===== PDF page 36 ===== -->

     doi:https://doi.org/10.1016/0004-3702(87)90004-X.

[59] Pee, E. Y., “On Solving Large-Scale Finite Minimax Problems using Exponential Smoothing,” Journal of Optimization Theory

     and Applications, Vol. 148, No. 2, 2011, pp. 390–421. doi:https://doi.org/10.1007/s10957-010-9759-1.

[60] Polak, E., Royset, J. O., and Womersley, R. S., “Algorithms with Adaptive Smoothing for Finite Mini-

     max Problems,” Journal of Optimization Theory and Applications, Vol. 119, No. 3, 2003, pp. 459–484.

     doi:https://doi.org/10.1023/B:JOTA.0000006685.60019.3e.

[61] Li, X., “An Entropy-Based Aggregate Method for Minimax Optimization,” Engineering Optimization, Vol. 18, No. 4, 1997, pp.

     277–285. doi:https://doi.org/10.1080/03052159208941026.

[62] Rashad, R., “High-Fidelity Aerodynamic Shape Optimization for Natural Laminar Flow,” Ph.D. thesis, University of Toronto,

     2016.

[63] Richards, F. J., “A Flexible Growth Function for Empirical Use,” Journal of Experimental Botany, Vol. 10, No. 2, 1959, pp.

     290–301. doi:https://doi.org/10.1093/jxb/10.2.290.

[64] Hicken, J. E., and Zingg, D. W., “Parallel Newton-Krylov solver for the Euler equations discretized using simultaneous

     approximation terms,” AIAA Journal, Vol. 46, No. 11, 2008, pp. 2773–2786. doi:https://doi.org/10.2514/1.34810.

[65] Osusky, M., and Zingg, D. W., “Parallel Newton-Krylov-Schur Flow Solver for the Navier-Stokes Equations,” AIAA Journal,

     Vol. 51, No. 12, 2013, pp. 2833–2851. doi:https://doi.org/10.2514/1.J052487.

[66] Saad, Y., and Schultz, M. H., “GMRES: A generalized minimal residual algorithm for solving nonsymmetric linear systems,”

     SIAM Journal on Scientific and Statistical Computing, Vol. 7, No. 3, 1986, pp. 856–869. doi:https://doi.org/10.1137/0907058.

[67] Saad, Y., and Sosonkina, M., “Distributed Schur complement techniques for general sparse linear systems,” SIAM Journal of

     Scientific Computing, Vol. 21, 1999, pp. 1337–1357. doi:https://doi.org/10.1137/S1064827597328996.

[68] Moryossef, Y., and Levy, Y., “Unconditionally positive implicit procedure for two-equation turbulence models:

     Application to k-ω turbulence models,” Journal of Computational Physics, Vol. 220, 2006, pp. 88–108.

     doi:https://doi.org/10.1016/j.jcp.2006.05.001.

[69] Lian, C., Xia, G., and Merkle, C. L., “Impact of source terms on reliability of CFD algorithms,” Computers & Fluids, Vol. 39,

     2010, pp. 1909–1922. doi:https://doi.org/10.1016/j.compfluid.2010.06.021.

[70] Moryossef,    Y.,   and Levy,     Y.,   “The unconditionally positive-convergent implicit time integration scheme

     for two-equation turbulence models: Revisited,”           Computers & Fluids,         Vol. 38,     2009,   pp. 1984–1994.

     doi:https://doi.org/10.1016/j.compfluid.2009.06.005.

[71] Moryossef, Y., “Unconditionally stable time marching scheme for Reynolds stress models,” Journal of Computational Physics,

     Vol. 276, 2014, pp. 635–664. doi:https://doi.org/10.1016/j.jcp.2014.07.047.

[72] Mosahebi, A., and Laurendeau, E., “Convergence Characteristics of Fully and Loosely Coupled Numerical Approaches for

     Transition Models,” AIAA Journal, Vol. 53, No. 5, 2015, pp. 1399–1404. doi:https://doi.org/10.2514/1.J053722.

[73] Mosahebi, A., and Laurendeau, E., “Introduction of a modified segregated numerical approach for efficient simulation

     of γ- Re
            ˜ θt transition model,” International Journal of Computational Fluid Dynamics, Vol. 29, 2015, pp. 357–375.

     doi:https://doi.org/10.1080/10618562.2015.1093624.

<!-- ===== PDF page 37 ===== -->

[74] Candler, G. V., Subbareddy, P. K., and Nompelis, I., “Decoupled Implicit Method for Aerothermodynamics and Reacting

     Flows,” AIAA Journal, Vol. 51, No. 5, 2013, pp. 1245–1254. doi:https://doi.org/10.2514/1.J052070.

[75] Chisholm, T. T., and Zingg, D. W., “A Jacobian-free Newton-Krylov algorithm for compressible turbulent fluid flows,” Journal

     of Computational Physics, Vol. 228, No. 9, 2009, pp. 3490–3507. doi:https://doi.org/10.1016/j.jcp.2009.02.004.

[76] Yildirim, A., Kenway, G. K. W., Mader, C. A., and Martins, J. R. R. A., “A Jacobian-free approximate Newton-Krylov startup

     strategy for RANS simulations,” Journal of Computational Physics, 2019. doi:https://doi.org/10.1016/j.jcp.2019.06.018.

[77] Lian, C., Xia, G., and Merkle, C. L., “Solution-limited time stepping to enhance reliability in CFD applications,” Journal of

     Computational Physics, Vol. 228, 2009, pp. 4836–4857. doi:https://doi.org/10.1016/j.jcp.2009.03.040.

[78] Modisette, J. M., “An automated reliable method for two-dimensional Reynolds-Averaged Navier-Stokes simulations,” Ph.D.

     thesis, Massachusetts Institute of Technology, 2011.

[79] Swanson, R. C., and Turkel, E., “On central-difference and upwind schemes,” Journal of Computational Physics, Vol. 101,

     1992, pp. 292–306. doi:https://doi.org/10.1016/0021-9991(92)90007-L.

[80] Somers, D. M., “Design and experimental results for a natural-laminar-flow airfoil for general aviation applications,” NASA

     TP-1861, 1981.

[81] Coder, J. G., “Standard Test Cases for Transition Model Verification and Validation in Computational Fluid Dynamics,” 2018

     AIAA Aerospace Sciences Meeting, AIAA Paper 2018-0029, 2018. doi:https://doi.org/10.2514/6.2018-0029.

[82] Denison, M., and Pulliam, T. H., “Implementation and Assessment of the Amplification Factor Transport Laminar-Turbulent

     Transition Model,” 2018 AIAA Aerospace Sciences Meeting, AIAA Paper 2018-3382, 2018. doi:https://doi.org/10.2514/6.2018-

     3382.

[83] Khayatzadeh, P., and Nadarajah, S. K., “Laminar-Turbulent Flow Simulation for Wind Turbine Profiles using the γ- Re
                                                                                                                       ˜ θt

     Transition Model,” Wind Energy, Vol. 17, No. 6, 2014, pp. 901–918. doi:https://doi.org/10.1002/we.1606.

[84] Somers, D. M., “Design and experimental results for the S809 airfoil,” National Renewable Energy Laboratory, Golden,

     Colorado, NRELSR-440-6918, 1989. doi:https://doi.org/10.2172/437668.

[85] Dagenhart, J., and Saric, W., “Crossflow stability and transition experiments in swept-wing flow,” NASA Langley Technical

     Report Server, Vol. NASA/TP-1999-209344, 1999.

[86] Radeztsky, R. H., Reibert, M. S., and Saric, W. S., “Effect of Micron-Sized Roughness on Transition in Swept-Wing Flows,”

     31st Aerospace Sciences Meeting and Exhibit, AIAA Paper 93-0076, 1993. doi:https://doi.org/10.2514/6.1993-76.

[87] Jameson, A., Schmidt, W., and Turkel, E., “Numerical solution of the Euler equations by finite volume methods us-

     ing Runge-Kutta time-stepping schemes,” 14th Fluid and Plasma Dynamics Conference, AIAA Paper 81-1259, 1981.

     doi:https://doi.org/10.2514/6.1981-1259.

[88] Pulliam, T. H., “Efficient solution methods for the Navier-Stokes equations,” Tech. rep., Lecture Notes for the von Karman Inst.

     for Fluid Dynamics Lecture Series: Numerical Techniques for Viscous Flow Computation in Turbomachinery Bladings, 1986.

[89] Petzold, R., and Radespiel, R., “Transition on a Wing with Spanwise Varying Crossflow and Linear Stability Analysis,” AIAA

     Journal, Vol. 53, No. 2, 2015, pp. 321–335. doi:https://doi.org/10.2514/1.J053127.

[90] Kruse, M., Munoz, F., and Radespiel, R., “Transition Prediction Results for Sickle Wing and NLF(1)-0416 Test Cases,” 2018

<!-- ===== PDF page 38 ===== -->

     AIAA Aerospace Sciences Meeting, AIAA Paper 2018-0537, 2018. doi:https://doi.org/10.2514/6.2018-0537.

[91] Menter, F. R., Smirnov, P. E., Liu, T., and Avancha, R., “A One-Equation Local Correlation-Based Transition Model,” Flow,

     Turbulence and Combustion, Vol. 95, 2015, pp. 583–619. doi:https://doi.org/10.1007/s10494-015-9622-4.
