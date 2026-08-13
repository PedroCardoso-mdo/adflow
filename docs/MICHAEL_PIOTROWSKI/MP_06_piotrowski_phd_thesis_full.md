# Piotrowski (PhD thesis, University of Toronto) — Development of a Transition Prediction Methodology Suitable for Aerodynamic Shape Optimization

> **Full transcript** — verbatim text extracted from the source PDF with
> `pdftotext -layout`, page-delimited, nothing removed except bare page-number
> lines. Layout is preserved, so tables and multi-column material read as they
> appear in the PDF. **Equations are raw PDF extractions**: superscripts,
> subscripts and math symbols may be flattened or split across lines — for any
> equation you intend to implement, confirm against the PDF itself (path below)
> or the digest companion file.
>
> Source PDF: `/home/mdo/Desktop/ARTIGOS SMOTH/Piotrowski_Michael_GH_PHD_thesis (1).pdf`

---

<!-- ===== PDF page 1 ===== -->

Development of a Transition Prediction Methodology Suitable for
                Aerodynamic Shape Optimization

                                   by

                        Michael G. H. Piotrowski

          A thesis submitted in conformity with the requirements
                   for the degree of Doctor of Philosophy
        Graduate Department of Aerospace Science and Engineering
                         University of Toronto

              c Copyright 2022 by Michael G. H. Piotrowski

<!-- ===== PDF page 3 ===== -->

  Development of a Transition Prediction Methodology Suitable for Aerodynamic Shape
                                     Optimization

                                  Michael G. H. Piotrowski
                                   Doctor of Philosophy
               Graduate Department of Aerospace Science and Engineering
                                University of Toronto

                                         Abstract

Aircraft configurations that exploit significant regions of natural laminar flow could play a
key role in reducing the environmental impact of aviation. The design of natural-laminar-
flow configurations can be accelerated using computational tools capable of accurately and
efficiently predicting boundary-layer transition. Toward the development of a suitable com-
putational design tool, methods for the efficient prediction of boundary-layer transition in
a Reynolds-averaged Navier-Stokes-based flow solver and integration in a discrete-adjoint
gradient-based optimization algorithm are presented. A local correlation-based transition
model is modified and coupled to the Spalart-Allmaras turbulence model and integrated in a
Newton-Krylov-Schur flow solver. Modifications to the solution strategy are introduced, in-
cluding a source-term time step restriction, in order to prevent unstable solution updates for
the fully coupled, fully implicit solver. A smooth-variant of the transition model is developed
with approximations to discontinuous and stiff source-term functions. Both transition mod-
els are validated using two- and three-dimensional subsonic transition test cases, with the
new, smooth model producing significantly improved iterative convergence. Compressibility
corrections are developed and applied to extend the transition model empirical correlations
to transonic flow regimes typical of commercial transport aircraft, with the resulting model
investigated using two- and three-dimensional transonic transition test cases. The results
demonstrate that the compressibility corrections produce substantially improved agreement
with the experimental transition locations, particularly for higher Reynolds number applica-
tions. Finally, the smooth transition model is integrated in a discrete-adjoint gradient-based
optimization algorithm and applied to two- and three-dimensional drag-minimization stud-
ies across a range of design conditions. The results demonstrate that the capability of the
current framework to explore the natural-laminar-flow configuration design space is sensitive
to the streamwise grid resolution in the transition regions, with grid requirements increas-
ing as the transition length decreases with increasing Reynolds number. The light aircraft
design space is found to be multi-modal, with the optimization framework producing two
distinct local minima. Transonic drag minimizations demonstrate that the optimization
framework can successfully trade a decrease in viscous drag with an increase in wave drag,
and an infinite swept wing optimization demonstrates that the framework can delay both
Tollmien-Schlichting and stationary crossflow instabilities.

                                              ii

<!-- ===== PDF page 4 ===== -->

iii

<!-- ===== PDF page 5 ===== -->

                                   Acknowledgments

   It has been a privilege to work with Professor Zingg over the last several years. Beyond
his knowledge of and passion for CFD and aircraft design is his immense dedication and
commitment to his students. I have benefited greatly from his leadership, curiosity, and
drive, both on and off the basketball court, and I could not have asked for a better mentor
to have guided me throughout this journey, thank you.
   Thank you to my committee members, Drs. Yano and Lavoie, for the thoughtful questions
and discussions over the years. Since our first discussion on numerical robustness, Dr. Yano
has always found the key, illuminating questions to ask, while Dr. Lavoie’s guidance has
helped to keep the big picture in sight. I thank Drs. Allmaras and Sullivan for their involve-
ment in the final oral examination. In particular, I greatly appreciate the insight provided
by Dr. Allmaras in writing the appraisal letter. It was a privilege to discuss aspects of the
thesis with someone so deeply connected to the field.
   To Dr. Kiris, thank you for providing me with the opportunity to work in the LAVA group,
and for your unwavering support, mentorship, and guidance. Being exposed to such high-
calibre research has improved the quality of this work immeasurably. Thank you to everyone
in the LAVA group, whose dedication and work ethic I strive to match. Specifically, to
Drs. Housman and Kenway, thank you for your unfailing patience in answering my many
questions. To everyone in the Computational Sciences Branch, thank you for welcoming
me in with open arms. Thank you Dr. Pulliam for opening the door, both metaphorically
and physically on the first day, and for your generous support and guidance, as well as
the occasional hiking advice. A special thank you to Coach Jeff and the running group for
helping me survive my first triathlon, the beach volleyball group, and the hiking group for
keeping me sane during the pandemic.
   Thank you Dr. Coder for the conversations and support over the years. To Dr. Costantini,
thank you for providing extensive data for the VA-2 experiment, as well as for the several
insightful conversations on the intricacies of wind-tunnel experimental campaigns.
   Thank you to everyone at UTIAS for creating such an enjoyable, hidden, academic oasis.
To Drs. Reist and Apponsah, thank you for helping me learn all the intricacies and surprises
that Jetstream has to offer. To everyone in the CFD lab past and present, thank you for
the comradery and the many basketball, ultimate frisbee, and softball games. To Gregg and
Tim, thank you for your friendship — we finally made it.
   Finally, to my Family. Thank you to my parents, and to Jilly Willy, for providing me
with the opportunities and support to push myself further than I ever thought possible, and
for keeping me fed in the process. To my sister, Alicia, thank you for always providing me
with a role model to look up to. I’ll continue to strive to follow in your footsteps. Lastly, to
Taryn, for being my steadfast proofreader and partner in crime, thank you.
   Financial support from the Natural Sciences and Engineering Research Council of Canada,

                                               iv

<!-- ===== PDF page 6 ===== -->

the University of Toronto, and the NASA Advanced Air Transport Technology and Transfor-
mational Tools and Technologies projects is gratefully acknowledged, as well as the compu-
tational resources and support provided by the SciNet High-Performance Computing Con-
sortium.

                                                                  Michael Piotrowski

University of Toronto Institute for Aerospace Studies
May, 2022

                                             v

<!-- ===== PDF page 7 ===== -->

Contents

List of Symbols and Abbreviations                                                            xviii

1 Introduction                                                                                  1
  1.1 Motivation . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .      1
  1.2 Boundary-Layer Transition . . . . . . . . . . . . . . . . . . . . . . . . . . . .         2
       1.2.1 Streamwise Instabilities . . . . . . . . . . . . . . . . . . . . . . . . .         3
       1.2.2 Spanwise Instabilities . . . . . . . . . . . . . . . . . . . . . . . . . . .       5
  1.3 Transition Prediction and Modelling Methods . . . . . . . . . . . . . . . . .             7
       1.3.1 eN Criterion . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .       8
       1.3.2 Transport-Equation-Based Transition Models . . . . . . . . . . . . .               8
       1.3.3 Effects of Flow Compressibility . . . . . . . . . . . . . . . . . . . . .         10
  1.4 Shape Optimization with Boundary-Layer Transition . . . . . . . . . . . . .              11
  1.5 Objectives . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .     14
  1.6 Thesis Outline . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .     14

2 Transition Model Equations                                                                   15
  2.1 SA-LM2015 Transition Model . . . . . . . . . . . . . . . . . . . . . . . . . .           15
      2.1.1 Momentum-Thickness Reynolds Number Transport Equation . . . . .                    16
      2.1.2 Intermittency Transport Equation . . . . . . . . . . . . . . . . . . . .           18
      2.1.3 Coupling to SA Turbulence Model . . . . . . . . . . . . . . . . . . . .            19
      2.1.4 Initial and Boundary Conditions . . . . . . . . . . . . . . . . . . . . .          19
  2.2 SA-sLM2015 Transition Model . . . . . . . . . . . . . . . . . . . . . . . . . .          20
      2.2.1 Minimum/Maximum Operators . . . . . . . . . . . . . . . . . . . . .                20
      2.2.2 Intermittency Source Terms . . . . . . . . . . . . . . . . . . . . . . .           22
      2.2.3 Momentum-Thickness Reynolds Number Empirical Correlations . . .                    26
  2.3 Compressibility Corrections . . . . . . . . . . . . . . . . . . . . . . . . . . .        27
      2.3.1 Tollmien-Schlichting Instabilities . . . . . . . . . . . . . . . . . . . .         28
      2.3.2 Stationary Crossflow Instabilities . . . . . . . . . . . . . . . . . . . .         33

3 Optimization Framework                                                                       37
  3.1 Parallel Implicit Newton-Krylov-Schur Flow Solver . . . . . . . . . . . . . .            37

                                               vi

<!-- ===== PDF page 8 ===== -->

         3.1.1 Considerations for the Transition Model Equations . . . . . . . . . .        39
         3.1.2 Equation and Variable Scaling . . . . . . . . . . . . . . . . . . . . . .    41
         3.1.3 Solution Update Damping . . . . . . . . . . . . . . . . . . . . . . . .      42
         3.1.4 Source-Term Time Step Restriction . . . . . . . . . . . . . . . . . . .      44
   3.2   Gradient-based Optimization . . . . . . . . . . . . . . . . . . . . . . . . . .    47
         3.2.1 Optimization Algorithm . . . . . . . . . . . . . . . . . . . . . . . . .     47
         3.2.2 Gradient Evaluation . . . . . . . . . . . . . . . . . . . . . . . . . . .    47

4 Results: Analysis                                                                         51
  4.1 Subsonic Transition Test Cases . . . . . . . . . . . . . . . . . . . . . . . . .      52
      4.1.1 NLF0416 General Aviation Airfoil . . . . . . . . . . . . . . . . . . . .        53
      4.1.2 S809 Wind Turbine Airfoil . . . . . . . . . . . . . . . . . . . . . . . .       55
      4.1.3 NASA NLF2-0415 Infinite Swept Wing . . . . . . . . . . . . . . . . .            58
      4.1.4 TU Braunschweig Sickle Wing . . . . . . . . . . . . . . . . . . . . . .         60
  4.2 Transonic Transition Test Cases . . . . . . . . . . . . . . . . . . . . . . . . .     61
      4.2.1 CAST10-2 Airfoil . . . . . . . . . . . . . . . . . . . . . . . . . . . . .      64
      4.2.2 VA-2 Supercritical Airfoil . . . . . . . . . . . . . . . . . . . . . . . .      68
      4.2.3 NASA CRM-NLF Wing-Body Geometry . . . . . . . . . . . . . . . .                 73

5 Results: Optimization                                                                     85
  5.1 Airfoil Optimization . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .    85
      5.1.1 Cessna 172R Skyhawk . . . . . . . . . . . . . . . . . . . . . . . . . .         86
      5.1.2 De Havilland Dash8-Q400 . . . . . . . . . . . . . . . . . . . . . . . .         90
      5.1.3 Boeing 737-800 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .      93
  5.2 Infinite Swept Wing Optimization . . . . . . . . . . . . . . . . . . . . . . . .      95

6 Contributions, Conclusions, and Recommendations                                           97
  6.1 Contributions and Conclusions . . . . . . . . . . . . . . . . . . . . . . . . . .     97
  6.2 Recommendations . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .     99

A Fonset,scf Source-Term Validation                                                   115
  A.1 TU Braunschweig Sickle Wing . . . . . . . . . . . . . . . . . . . . . . . . . . 115

B Further Evaluation of Iterative and Grid Convergence                                      117
  B.1 Methodology . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 118
      B.1.1 Linearization Strategies . . . . . . . . . . . . . . . . . . . . . . . . . 118
      B.1.2 Estimating Grid Convergence . . . . . . . . . . . . . . . . . . . . . . 119
  B.2 Results . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 120
      B.2.1 NLF0416 General Aviation Airfoil . . . . . . . . . . . . . . . . . . . . 120
      B.2.2 VA-2 Supercritical Airfoil . . . . . . . . . . . . . . . . . . . . . . . . 123
      B.2.3 NASA CRM-NLF Wing-Body Geometry . . . . . . . . . . . . . . . . 125

                                              vii

<!-- ===== PDF page 9 ===== -->

C Transition Length Modification                                                         129
  C.1 Modified Intermittency Source Terms . . . . . . . . . . . . . . . . . . . . . . 129
  C.2 Airfoil Optimization . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 132
      C.2.1 Cessna 172R Skyhawk . . . . . . . . . . . . . . . . . . . . . . . . . . 132
      C.2.2 De Havilland Dash8-Q400 . . . . . . . . . . . . . . . . . . . . . . . . 133
      C.2.3 Boeing 737-800 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 134
  C.3 Summary . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 135

                                            viii

<!-- ===== PDF page 10 ===== -->

ix

<!-- ===== PDF page 11 ===== -->

List of Tables

 4.1   TMPCS NLF0416 structured C-grid dimensions [27], where ∆s represents the
       first off-wall grid spacing. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .   53
 4.2   TMPCS S809 structured C-grid dimensions [27]. . . . . . . . . . . . . . . . .            56
 4.3   Sickle Wing structured O-O grid dimensions. . . . . . . . . . . . . . . . . . .          61
 4.4   CAST10-2 structured O-grid dimensions. . . . . . . . . . . . . . . . . . . . .           65
 4.5   VA-2 structured O-grid dimensions. . . . . . . . . . . . . . . . . . . . . . . .         69
 4.6   CRM-NLF wind tunnel test conditions . . . . . . . . . . . . . . . . . . . . .            74
 4.7   CRM-NLF structured grid characteristics. . . . . . . . . . . . . . . . . . . .           74

 5.1   Design conditions for the two-dimensional lift-constrained drag minimiza-
       tions [150]. For each case turbulence intensity, T u, is specified as 0.07%. . . .       86
 5.2   Structured multi-block O-grid dimensions for the two-dimensional optimiza-
       tions at the Cessna 172R design conditions (see Table 5.1). . . . . . . . . . .          86
 5.3   Aerodynamic performance of the initial and optimized designs at the Cessna
       172R conditions on the L1 grid level. . . . . . . . . . . . . . . . . . . . . . .        88
 5.4   Structured multi-block O-grid dimensions for the two-dimensional optimiza-
       tions at the Dash8-Q400 design conditions (see Table 5.1). . . . . . . . . . .           90
 5.5   Aerodynamic performance of the initial and optimized designs at the Dash8-
       Q400 conditions on the L2 grid level. . . . . . . . . . . . . . . . . . . . . . .        92
 5.6   Structured multi-block O-grid dimensions for the two-dimensional optimiza-
       tions at the sweep-corrected Boeing 737-800 design conditions (see Table 5.1).           93
 5.7   Aerodynamic performance of the initial and optimized designs at the sweep-
       corrected Boeing 737-800 conditions on the L2 grid level. . . . . . . . . . . .          94
 5.8   Aerodynamic performance of the initial and optimized infinite swept wing
       designs at the transonic strut-braced wing aircraft design conditions. . . . .           96

 B.1 NLF0416 structured multi-block C-grid dimensions. . . . . . . . . . . . . . . 120

 C.1 Design conditions for the extruded RAE2822 airfoil simulated at the swept
     and sweep-corrected Boeing 737-800 design conditions. For each configuration
     turbulence intensity, T u, is specified as 0.07%. . . . . . . . . . . . . . . . . . 132

                                               x

<!-- ===== PDF page 12 ===== -->

xi

<!-- ===== PDF page 13 ===== -->

List of Figures

 1.1   Historical barrier for NLF designs based on experiments with varying combi-
       nations of Reynolds number and wing leading-edge sweep angle, along with
       the primary transition mechanisms. Reproduced from [63] with CRM-NLF
       data provided by [104]. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .    2
 1.2   Natural transition process over a flat plate. Reproduced from [191]. . . . . .         3

 2.1   max(0, x − 1) approximated using φp with a range of values for the smoothing
       parameter, p. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .   22
 2.2   Original and smooth Flength functions. . . . . . . . . . . . . . . . . . . . . . .    23
 2.3   Original and smooth Reθc functions. . . . . . . . . . . . . . . . . . . . . . . .     24
 2.4   Original and smooth intermittency source terms, Sγ = Pγ − Eγ , plotted ver-
       sus the transition criterion (Equation 2.35), intermittency, γ, and the eddy
       viscosity ratio, RT , with Flength = Ω = 1. . . . . . . . . . . . . . . . . . . . .   25
 2.5   Original and smooth F (λθ ) functions for a turbulence intensity of 0.05%. . .        26
                                ±
 2.6   Original and smooth f (∆Hcrossflow ) functions. . . . . . . . . . . . . . . . . . .   27
 2.7   Sensitivity of the LM2009 empirical correlation [92], the stability-based model,
       and stability analysis [136] to pressure gradient and Mach number. A higher
       momentum-thickness Reynolds number delays boundary-layer transition. . .              29
 2.8   Sensitivity of the normalized stability-based model (blue), stability analy-
       sis [136], and the initial compressibility correction, ψinit (grey) (Equations 2.67–
       2.69), to Mach number with varying turbulence intensity and pressure gradient. 30
 2.9   Sensitivity of the initial, ψinit (grey) (Equations 2.67–2.69), and calibrated, ψ
       (red) (Equations 2.70–2.72), compressibility corrections to Mach number with
       varying turbulence intensity and pressure gradient. . . . . . . . . . . . . . .       32
 2.10 Effects of Mach number on the TS and stationary crossflow instability com-
      pressibility corrections, ψ and ψscf , respectively. . . . . . . . . . . . . . . . .   36

 3.1   Normalized solution update over a range of time steps produced by implicit
       and explicit Euler time marching applied to a linear convection equation with
       a source term. Reproduced from [100]. . . . . . . . . . . . . . . . . . . . . .       45
 3.2   Infinite swept wing geometry with FFD design variables in blue. . . . . . . .         48

                                             xii

<!-- ===== PDF page 14 ===== -->

3.3   Verification for the analytical flow Jacobian, metric linearization, and direc-
      tional derivative including the SA-sLM2015cc transition model equations us-
      ing both complex-step and finite-difference approximations. . . . . . . . . . .         49

4.1  NLF0416 grid convergence at 0◦ angle of attack, and angle-of-attack sweep
     simulated using the SA-LM2015 and SA-sLM2015 transition models at M =
     0.1, Re = 4.0 × 106 , and T u = 0.15%. The ‘fine’ grid level is used for the
     sweep, with the experimental results obtained at the LTPT wind tunnel [171].             54
4.2 NLF0416 residual convergence histories for the SA-LM2015 and SA-sLM2015
     transition models on the ‘fine’ grid level. The airfoil is simulated at M = 0.1,
     Re = 4.0 × 106 , T u = 0.15% and over a range of angles of attack. . . . . . .           55
4.3 S809 grid convergence at 6◦ angle of attack, and angle-of-attack sweep simu-
     lated using the SA-LM2015 and SA-sLM2015 transition models at M = 0.1,
     Re = 2.0 × 106 , and T u = 0.07%. The ‘fine’ grid level is used for the sweep,
     with the experimental results obtained at the TU Delft wind tunnel [172]. . .            57
4.4 S809 residual convergence histories for the SA-LM2015 and SA-sLM2015 tran-
     sition models on the ‘fine’ grid level. The airfoil is simulated at M = 0.1,
     Re = 2.0 × 106 , T u = 0.07% and over a range of angles of attack. . . . . . .           58
4.5 Experimental and computational results for the NLF2-0415 infinite swept
     wing [148] simulated using the SA-LM2015 and SA-sLM2015 transition mod-
     els at M = 0.15, α = −4◦ , T u = 0.20%, and over a range of Reynolds numbers
     from 1.8–3.5 × 106 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .     59
4.6 NLF2-0415 residual convergence histories for the SA-LM2015 and SA-sLM2015
     transition models. The infinite swept wing is simulated at M = 0.15, α = −4◦ ,
     T u = 0.20%, and over a range of Reynolds numbers. . . . . . . . . . . . . . .           60
4.7 Skin friction profiles on the upper and lower surfaces of the Sickle Wing for
     all grid levels overlaid with experimental transition locations [138]. Results
     are obtained using the SA-LM2015 and SA-sLM2015 transition models at
     M = 0.16, α = −2.6◦ , T u = 0.17%, h = 9.78µm, and a Reynolds number of
     2.75 × 106 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .   62
4.8 Experimental and ‘fine’ grid pressure coefficient profiles, and grid convergence
     for the SA-LM2015 and SA-sLM2015 transition models. The Sickle Wing is
     simulated at M = 0.16, α = −2.6◦ , T u = 0.17%, and a Reynolds number of
     2.75 × 106 [90]. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .     63
4.9 Relative residual and drag coefficient convergence histories for the SA-LM2015
     and SA-sLM2015 transition models on all grid levels. The Sickle Wing is
     simulated at M = 0.16, α = −2.6◦ , T u = 0.17%, h = 9.78µm, and a Reynolds
     number of 2.75 × 106 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .     63
4.10 CAST10-2 grid-convergence results produced at −0.39 and 0.82 degrees angle
     of attack and M = 0.74, Re = 2 × 106 , and T u = 0.25%. . . . . . . . . . . .            65

                                            xiii

<!-- ===== PDF page 15 ===== -->

4.11 CAST10-2 pressure and upper-surface skin friction coefficient profiles pro-
     duced at M = 0.74, Re = 2 × 106 , T u = 0.25%, and at three angles of attack
     overlaid with the pressure profiles from the experiment [66]. . . . . . . . . .        66
4.12 CAST10-2 upper-surface transition locations (dashed) and lift curve (solid)
     produced on the L3 grid at M = 0.74, Re = 2 × 106 , T u = 0.25%, and over a
     range of angles of attack compared with the results from the experiment [66].          67
4.13 VA-2 grid-convergence results produced at −0.40 and 1.80 degrees angle of
     attack and M = 0.71, Re = 10 × 106 , and T u = 0.25%. . . . . . . . . . . . .          69
4.14 VA-2 pressure and upper-surface skin friction coefficient profiles produced at
     M = 0.71, Re = 10 × 106 , T u = 0.25%, and over a range of angles of attack
     overlaid with the results from the experiment [34]. . . . . . . . . . . . . . . .      71
4.15 VA-2 upper-surface transition locations produced on the L3 grid at M = 0.71,
     Re = 10 × 106 , T u = 0.25%, and over a range of angles of attack compared
     with results from the experiments [34]. . . . . . . . . . . . . . . . . . . . . .      73
4.16 CRM-NLF grid-convergence results at the 2524 test conditions (α ≈ 2.0◦ ). .            75
4.17 CRM-NLF grid-refinement study residual convergence histories at the 2524
     test conditions (α ≈ 2.0◦ ). . . . . . . . . . . . . . . . . . . . . . . . . . . . .   75
4.18 CRM-NLF grid-refinement study pressure and skin friction coefficient profiles
     produced at the 2524 test conditions (α ≈ 2.0◦ ) compared with the pressure
     profiles from the experiment [104] at varying spanwise stations, η. . . . . . .        77
4.19 CRM-NLF grid-refinement study upper-surface skin friction coefficient profiles
     at the 2524 test conditions (α ≈ 2.0◦ ) overlaid with the estimated transition
     front from the experiment [104]. . . . . . . . . . . . . . . . . . . . . . . . . .     79
4.20 CRM-NLF angle-of-attack sweep upper-surface skin friction coefficient profiles
     obtained on the L1 grid overlaid with the estimated transition front from the
     experiment [104]. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .    80
4.21 CRM-NLF angle-of-attack sweep force curves for lift, drag, and pitching mo-
     ment coefficient produced on the L1 grid compared with results from the
     experiment [104]. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .    81
4.22 CRM-NLF angle-of-attack sweep residual convergence histories produced on
     the L1 grid. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .   81
4.23 CRM-NLF pressure and skin friction coefficient profiles produced at the 2523
     test conditions (α ≈ 1.5◦ ) on the L1 grid compared with the pressure profiles
     from the experiment [104] at varying spanwise stations, η. . . . . . . . . . .         82
4.24 CRM-NLF pressure and skin friction coefficient profiles produced at the 2525
     test conditions (α ≈ 2.5◦ ) on the L1 grid compared with the pressure profiles
     from the experiment [104] at varying spanwise stations, η. . . . . . . . . . .         83

                                            xiv

<!-- ===== PDF page 16 ===== -->

4.25 CRM-NLF pressure and skin friction coefficient profiles produced at the 2526
     test conditions (α ≈ 3.0◦ ) on the L1 grid compared with the pressure profiles
     from the experiment [104] at varying spanwise stations, η. . . . . . . . . . .          84

5.1   Optimization convergence histories and cross-sectional profiles of the initial
      and optimized designs produced by drag minimizations with the baseline op-
      timization problem (Equations 5.1–5.4) at the Cessna 172R conditions with
      varying streamwise grid resolution. . . . . . . . . . . . . . . . . . . . . . . .      87
5.2   Optimization convergence histories and cross-sectional profiles of the initial
      and optimized designs produced by drag minimizations with the baseline
      (Equations 5.1–5.4) and more constrained optimization problem (Equations 5.5–
      5.9) at the Cessna 172R conditions on the L1 grid level. . . . . . . . . . . . . 89
5.3   Optimization convergence histories and cross-sectional profiles of the initial
      and optimized designs produced by drag minimizations with the baseline op-
      timization problem (Equations 5.1–5.4) at the De Havilland Dash8-Q400 con-
      ditions with varying streamwise grid resolution. . . . . . . . . . . . . . . . .       90
5.4   Optimization convergence histories and cross-sectional profiles of the initial
      and optimized designs produced by drag minimizations with the baseline
      (Equations 5.1–5.4) and more constrained optimization problem (Equations 5.5–
      5.9) at the De Havilland Dash8-Q400 conditions on the L2 grid level. . . . . 92
5.5   Optimization convergence histories and cross-sectional profiles of the initial
      and optimized designs produced by drag minimizations with the baseline opti-
      mization problem (Equations 5.1–5.4) at the sweep-corrected Boeing 737-800
      conditions with varying streamwise grid resolution. . . . . . . . . . . . . . .        94
5.6   Infinite swept wing drag minimization at the transonic strut-braced wing air-
      craft design conditions. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .   96

A.1 Sickle Wing grid-refinement study upper and lower surface skin friction coef-
    ficient profiles produced at M = 0.16, α = −2.6◦ , T u = 0.17%, h = 9.78µm,
    and Re = 2.75×106 overlaid with transition locations from the experiment [138].116

B.1 Block-diagonal Jacobians produced by the four linearization strategies for the
    turbulence and transition model equations with varying levels of coupling,
    with the grey elements representing filled Jacobian entries. . . . . . . . . . . 118
B.2 Residual convergence histories for the NLF0416 airfoil simulated fully turbu-
    lent and with free transition at M a = 0.10, Re = 4.0 × 106 , T u = 0.15%, and
    α = 0◦ using the grids presented in Table B.1. . . . . . . . . . . . . . . . . . 121
B.3 Linearization study residual convergence histories for the NLF0416 airfoil sim-
    ulated on the L4 grid (Table B.1) using the SA-sLM2015 transition model at
    M a = 0.10, Re = 4.0 × 106 , T u = 0.15%, and α = 0◦ . . . . . . . . . . . . . . 122

                                             xv

<!-- ===== PDF page 17 ===== -->

B.4 Grid-convergence results for the NLF0416 airfoil simulated fully turbulent and
    with free transition at M a = 0.10, Re = 4.0 × 106 , T u = 0.15%, and α = 0◦ .           122
B.5 Residual convergence histories for the VA-2 airfoil simulated fully turbulent
    and with free transition at M = 0.71, Re = 10 × 106 , T u = 0.25%, and
    α = −0.40◦ , 1.80◦ . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .   123
B.6 Linearization study residual convergence histories for the VA-2 airfoil simu-
    lated on the L3 grid (Table 4.5) using the SA-sLM2015cc transition model at
    M = 0.71, Re = 10 × 106 , T u = 0.25%, and α = −0.40◦ , 1.80◦ . . . . . . . . .          124
B.7 Grid-convergence results for the VA-2 airfoil simulated fully turbulent and
    with free transition at M = 0.71, Re = 10 × 106 , T u = 0.25%, and α =
    −0.40◦ , 1.80◦ . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .   125
B.8 Linearization study residual convergence histories for the NASA CRM-NLF
    simulated on the L2 grid (Table 4.7) using the SA-QCR2000-sLM2015cc tran-
    sition model at approximately M = 0.86, Re = 15 × 106 , T u = 0.24%, and
    α = 2.0◦ . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .   126
B.9 Grid-convergence results for the NASA CRM-NLF simulated fully turbulent
    and with free transition at approximately M = 0.86, Re = 15 × 106 , T u =
    0.24%, and α = 2.0◦ . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .      127

C.1 Sensitivity of the different Fonset formulations to the transition criterion (Equa-
    tion 2.35) and the eddy viscosity ratio. . . . . . . . . . . . . . . . . . . . . . 130
C.2 Pressure and upper-surface skin friction coefficient profiles produced by simu-
    lations of the VA-2 airfoil on the L0 grid level (Table 4.5) with varying Mach
    and Reynolds number; α = 0.00, T u = 0.25%. . . . . . . . . . . . . . . . . . 131
C.3 Pressure and skin friction coefficient profiles produced by simulations of a
    swept and sweep-corrected RAE2822 airfoil at the Boeing 737-800 design con-
    ditions (Table C.1) with the new intermittency source-term formulation. . . 132
C.4 Optimization convergence histories and cross-sectional profiles of the initial
    and optimized designs produced by drag minimizations with the baseline
    (Equations 5.1–5.4) and more constrained optimization problem (Equations 5.5–
    5.9) at the Cessna 172R conditions (Table 5.1) on the L1 grid level (Table 5.2).133
C.5 Optimization convergence histories and cross-sectional profiles of the initial
    and optimized designs produced by drag minimizations with the baseline
    (Equations 5.1–5.4) and more constrained optimization problem (Equations 5.5–
    5.9) at the De Havilland Dash8-Q400 conditions (Table 5.1) on the L1 grid
    level (Table 5.4). . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 134
C.6 Drag minimization with the baseline optimization problem (Equations 5.1–
    5.4) at the sweep-corrected Boeing 737-800 conditions (Table 5.1) on the L1
    grid level (Table 5.6). . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 135

                                            xvi

<!-- ===== PDF page 18 ===== -->

xvii

<!-- ===== PDF page 19 ===== -->

List of Symbols and Abbreviations

Symbols

a           sound speed
d           distance to nearest solid wall
Hcrossflow crossflow strength parameter
h           surface roughness height (rms)
l           characteristic length
M           Mach number
N           number of grid nodes
p           smooth maximum/minimum function parameter
pswitch     proximity switch cutoff value, 10−15
P           static pressure
Re          Reynolds number based on the characteristic length
Re˜ θt      transported local transition onset momentum-thickness Reynolds number
Reθc        critical momentum-thickness Reynolds number
Reθt        transition onset momentum-thickness Reynolds number
ReS         strain-rate magnitude Reynolds number
ReΩ         vorticity Reynolds number
RT          eddy viscosity ratio
s           streamwise coordinate
∆s          first off-wall grid spacing
S           strain-rate magnitude
S̃          modified strain-rate magnitude
Tu          turbulence intensity
~
U           unit velocity vector
U           local velocity magnitude
ui          Cartesian velocity component
y+          non-dimensional wall distance
~
Ω           vorticity vector
Ω           vorticity magnitude
Ωstreamwise streamwise vorticity, helicity

                                          xviii

<!-- ===== PDF page 20 ===== -->

γ          intermittency
ρ          density
µ          molecular viscosity
µt         eddy viscosity
ν̃         modified eddy viscosity
λθ         pressure gradient parameter
δ          boundary-layer thickness
θ          momentum thickness
θt         transition onset momentum thickness
γsource    largest positive source-term Jacobian eigenvalue
φp         smooth maximum/minimum function
φswitch    proximity switch for the smooth maximum/minimum function
ψ          compressibility correction for Tollmien-Schlichting instability correlations
ψscf       compressibility correction for stationary crossflow instability correlations
κ          heat capacity ratio
Subscripts
e          boundary-layer edge
comp       compressible
∞          far field
scf        stationary crossflow

Abreviations

AFT           amplification factor transport
AHD           Arnal-Habiballah-Delcourt
CFD           computational fluid dynamics
CRM           common research model
CRM-NLF       common research model natural-laminar-flow variant
DNS           direct numerical simulation
FFD           free-form deformation
LES           large-eddy simulation
LST           linear stability theory
NLF           natural laminar flow
PSE           parabolized stability equations
RANS          Reynolds-averaged Navier-Stokes
SA            Spalart-Allmaras
SST           shear stress transport
TS            Tollmien-Schlichting
TTBW          transonic truss-braced wing

                                             xix

<!-- ===== PDF page 21 ===== -->

Chapter 1

Introduction

1.1    Motivation

Laminar-flow wing designs produce an increased laminar extent of the boundary layer, result-
ing in a decrease in viscous drag [63]. These drag savings have the potential to significantly
reduce the fuel burn of transport aircraft, where viscous drag constitutes approximately 50%
of the total drag [15]. The potential of laminar-flow designs has been demonstrated by vari-
ous studies which suggest the application of laminar flow control to large commercial aircraft
can reduce aerodynamic drag by approximately 10% [109].
   In the near term, the application of natural-laminar-flow (NLF) designs has been limited
to low sweep and/or low Reynolds number configurations, such as commercial transport
aircraft winglet, tail, and nacelle design, and fuselage and wing design for lightweight business
aircraft [54]. The Boeing 787-8 and 777X represent the first commercial applications of NLF
nacelles to large transport aircraft, while the 737 Max aircraft exploits laminar flow for both
the nacelle and winglet designs [37]. In the future, the design of commercial transport aircraft
wings which exploit significant regions of laminar flow can result in more significant drag
reduction. The Airbus BLADE project is investigating this, while experimental studies of
the NASA common research model NLF variant (CRM-NLF) demonstrate that it is possible
to design transport wings with high sweep and Reynolds numbers with significant regions
of laminar flow [17, 104, 105]. There are several operational challenges that still need to be
addressed; however, the Boeing 757 EcoDemonstrator is investigating the effect of anti-insect
devices and surface coatings which are key enablers for the application of laminar flow on
wing surfaces.
   The reason for the limited application of NLF designs to date is presented in Figure 1.1,
which illustrates the historical barrier for NLF designs based on previous experiments with
varying combinations of Reynolds number and wing leading-edge sweep angle [63]. Signifi-
cant regions of laminar flow can be obtained for low sweep designs up to flight scale Reynolds
numbers, where the dominant transition mechanism is natural transition due to Tollmien-
Schlichting (TS) instabilities. As the wing leading-edge sweep angle increases, crossflow and

<!-- ===== PDF page 22 ===== -->

2                                                                                                   CHAPTER 1. INTRODUCTION

                                                      NLF limit

                                            degrees
                                                                      NASA CRM-NLF

                       Leading edge sweep
                                                                         crossflow

                                                                                     Tollmien-Schlichting

                                                      10          20          30          40         x10
                                                                  Reynolds number
Figure 1.1: Historical barrier for NLF designs based on experiments with varying combinations of Reynolds
number and wing leading-edge sweep angle, along with the primary transition mechanisms. Reproduced
from [63] with CRM-NLF data provided by [104].

attachment-line instabilities become dominant, which reduces the transition Reynolds num-
ber. However, through redesigning the cross-sections of the NASA common research model
(CRM) geometry, a team at NASA Langley has demonstrated that this historical NLF barrier
can be extended [17, 104]. Section shapes were designed using a crossflow-attenuated natural
laminar flow (CATNLF) method to produce pressure gradients that attenuate both crossflow
and TS instabilities, with geometric constraints enforced to limit Görtler and attachment-line
instabilities.
    Efficient high-fidelity numerical optimization algorithms provide the designers of next-
generation aircraft with powerful tools for the development of more fuel-efficient designs [79,
153, 154, 150, 22]. The integration of transition prediction methodologies into aerodynamic
shape optimization algorithms enables the study of various trade-offs in the design of novel
NLF wing configurations. The design of aircraft configurations using design tools capable
of incorporating and exploiting NLF, which requires the ability to predict laminar-turbulent
transition efficiently, could play a key role in reducing the environmental impact of aviation.

1.2     Boundary-Layer Transition

There are several mechanisms for boundary-layer transition. The mechanism of transition
depends on the characteristics of the flow, such as Reynolds number, angle of attack, pres-
sure gradient, and turbulence intensity, which depend on the application being studied, such
as flow past wings, airframes, rotor blades, and wind turbine blades. For a swept transport-
aircraft wing the dominant mechanisms for transition result from TS wave growth leading
to natural transition, crossflow instabilities resulting from highly swept wings, concave cur-
vature producing Görtler instabilities, and attachment-line instabilities as a result of large
leading-edge radius and sweep at the root of the wing [151]. While the last two mechanisms

<!-- ===== PDF page 23 ===== -->

1.2. BOUNDARY-LAYER TRANSITION                                                                   3

                                         spanwise vorticity           turbulent spots

                                     TS waves         3D vortex breakdown

                       U

                                         Recrit                             Retran

                       U

                               laminar                 transitional                  turbulent

            Figure 1.2: Natural transition process over a flat plate. Reproduced from [191].

can be prevented by appropriate profile design, balancing natural and crossflow-induced
transition can be difficult, as favourable pressure gradients used to stabilize streamwise in-
stabilities destabilize crossflow instabilities [151].
   The transition process can be broadly broken down into three phases: receptivity, instabil-
ity growth or decay, and nonlinear breakdown. For low disturbance environments, boundary-
layer transition follows all three steps in what is known as natural transition (Figure 1.2).
However, for high disturbance environments, produced by freestream turbulence, pressure
gradients, or roughness elements, parts of the natural transition process can be skipped. The
transition from a laminar to a turbulent boundary layer results in an increase in skin friction
drag due to an increased off-wall velocity gradient and in pressure drag due to an increased
boundary-layer displacement thickness. The location of transition onset can be identified as
the location where the skin friction coefficient first begins to increase [114]. This section will
discuss each transition mechanism and its impact on the transition process.

1.2.1   Streamwise Instabilities
Natural Transition

Natural transition is characterized by the exponential growth of two-dimensional TS insta-
bilities, leading to nonlinear disturbances and eventually fully turbulent flow. Figure 1.2
illustrates the process of natural transition over a flat plate.
    In the first region of transition, the laminar boundary layer is destabilized by the presence
of viscous instability waves, TS waves [163]. Instabilities enter the boundary layer in a process
known as receptivity [19]. When the freestream turbulence level is low, i.e. below 1% [113],
a laminar boundary layer becomes linearly unstable beyond a critical Reynolds number at
which TS waves start to grow [91]. This instability is due to viscosity effects which destabilize
the TS waves, causing them to slowly grow throughout the transition process, leading to a
breakdown of the laminar boundary layer. This process might not be complete until a

<!-- ===== PDF page 24 ===== -->

4                                                                     CHAPTER 1. INTRODUCTION

streamwise distance as large as 20 times farther downstream from the leading edge than
the initial starting position of linear instability [192]. Inviscid mechanisms produce three-
dimensional disturbances once the waves have become nonlinear in the form of turbulent
spots. These spots grow in the surrounding laminar flow until a fully turbulent boundary
layer is produced.
   This transition process can be accelerated in the presence of an adverse pressure gradi-
ent, which forms an inflection point in the boundary-layer velocity profile. According to
Rayleigh’s criterion, an inflection point in the boundary-layer velocity profile is a necessary
condition for an inviscid instability to occur [19]. In such a flow environment, transition to a
turbulent boundary layer occurs before the formation of turbulent spots can be observed [72].

Bypass Transition

In contrast to natural transition, bypass transition occurs at high freestream turbulence
levels, i.e. greater than 1%. It is usually associated with turbomachinery applications where
upstream blade rows generate large levels of turbulence. For bypass transition, transition
occurs through the diffusion of freestream turbulence into the boundary layer. This results
in the instability growth or decay phase of natural transition being bypassed. Currently,
researchers have been unable to detect TS waves when the freestream turbulence level is
greater than 1% [113]. Nonlinear instabilities form directly in the laminar boundary layer,
leading to the formation of turbulent spots. Alternatively, surface roughness can create
disturbances at the wall, which can also lead to bypass transition. In the case of bypass
transition, the flow before the transition location is fully laminar, as opposed to natural
transition where flow disturbances, such as TS waves, are evident before this point. As
with natural transition, the development of these turbulent spots eventually leads to the
formation of fully turbulent flow.

Separation-Induced Transition

Another important mechanism for transition is separation-induced transition. This oc-
curs when a laminar boundary layer separates in the presence of a strong adverse pres-
sure gradient. Transition may occur in the separated shear layer, developing a laminar-
separation/turbulent-reattachment bubble [113]. This transition process can include many
of the steps involved in natural transition. The length of this bubble is dependent on the
Reynolds number, angle of attack, surface roughness, and turbulence intensity. Separation
bubbles are usually classified as being either long or short. Short separation bubbles af-
fect the local pressure distribution and can be useful for tripping the boundary layer to a
fully turbulent profile. Turbulent boundary-layer profiles allow for larger adverse pressure
gradients before the aerodynamic body encounters stalling, as turbulent boundary layers
are able to sustain an adverse pressure gradient and prevent flow separation better than an
equivalent laminar boundary layer. Long separation bubbles influence the overall pressure

<!-- ===== PDF page 25 ===== -->

1.2. BOUNDARY-LAYER TRANSITION                                                                    5

distribution and can generate significant pressure losses. Depending on flow conditions, long
separation bubbles can develop from short bubbles and vice versa, this process is known as
bursting [114]. The sudden change from a short to a long separation bubble can result in a
drastic loss of lift and cause the airfoil to stall if the separated flow does not reattach.

Görtler Instabilities

While TS instabilities result from small disturbances superimposed to the flow, Görtler
instabilities are initiated by centrifugal instability resulting from the curvature of concave
surfaces [19]. Görtler instabilities consist of vortices with axes oriented in the streamwise
direction. The strength of a Görtler instability can be evaluated using a modified form of
the Taylor number,
                                                     r
                                                         θ
                                         T a = Reθ          ,                                (1.1)
                                                         r0
where θ is the momentum thickness, Reθ is the Reynolds number based on the momentum
thickness, and r0 is the radius of curvature. Dryden demonstrated that unstable Taylor
numbers range between 6 and 9 depending on the freestream turbulence intensity of the
flow [47]. Forest [53] developed a correlation that can be used to relate the Taylor number
to the freestream turbulence intensity, T u, as follows:

                                     T a = 9 exp(−17.3T u).                                  (1.2)

Görtler instabilities play a significant role in the design of supercritical laminar airfoils that
often have concave regions near the lower surface leading and trailing edges.

1.2.2    Spanwise Instabilities

Two additional transition mechanisms are introduced when the flow becomes three-
dimensional due to spanwise velocity gradients, such as those introduced by wing-
sweep: crossflow instabilities leading to crossflow-induced transition and leading-edge
contamination/attachment-line transition.

Crossflow-Induced Transition

The natural transition process can be accelerated in the case of swept wings, where a signif-
icant velocity component in the sweep direction can form crossflow instabilities. Crossflow
instabilities are produced by streamwise vortices roughly aligned with the external velocity
streamlines. These vortices deform the streamwise velocity profile and produce an inflection
point in the boundary-layer velocity profile, which can result in transition occurring further
upstream than for an equivalent zero-sweep design [19, 72]. This inflectional instability can

<!-- ===== PDF page 26 ===== -->

6                                                                      CHAPTER 1. INTRODUCTION

trigger transition in regions of favourable pressure gradient. The receptivity process for cross-
flow instabilities is different from viscous TS instabilities, which are sensitive to freestream
disturbances such as turbulence intensity. There are two forms of crossflow instability, sta-
tionary and traveling [36]. Stationary instabilities are produced by steady micron-sized sur-
face roughness, while unsteady sources, such as freestream turbulence growth are responsible
for traveling crossflow instabilities. Although traveling crossflow instabilities can be highly
amplified, boundary-layer transition is determined by the characteristics of the stationary
instabilities [148]. Crossflow instabilities are heavily dependent on the pressure gradient.
Favourable pressure gradients, which dampen the growth of TS instabilities, amplify the
growth of crossflow instabilities as the flow acceleration increases the crossflow velocity.

Leading-Edge Contamination and Attachment-Line Transition

Primary attachment lines form on the leading edge of swept wings where the inviscid and
boundary-layer flow diverge symmetrically with respect to the upper and lower surface of
the wing. In the spanwise direction there exists a non-zero velocity and therefore a boundary
layer forms along the attachment line. Similar to the upper and lower surfaces of the wing,
this boundary layer can be laminar, transitional, or turbulent. This process of transition
is referred to as attachment-line transition [72]. A fully turbulent boundary layer that
forms upstream of the wing, such as on the fuselage or at the fuselage/wing junction, can
contaminate the leading edge of the wing. This results in a fully turbulent boundary layer
forming at the attachment line and is known as leading-edge contamination [19, 72]. Leading-
edge contamination and attachment-line transition can be evaluated using the spanwise
Reynolds number based on momentum thickness, which is defined by,

                                         we θz
                                 Reθz =        ,                                            (1.3)
                                           ν
                                         Z δ          
                                              w     w
                                    θz =         1−     ∂y,                                 (1.4)
                                          0 we      we

where ν is the kinematic viscosity, δ is the boundary-layer thickness, w is the spanwise
velocity, θz is the momentum thickness along the attachment line, and the subscript ‘e’
denotes a value at the boundary-layer edge. According to experiment, all disturbances will
damp out for a value of Reθz less than 100, small disturbances damp out for a value less
than 230, and laminar flow is unattainable for values larger than 240 [19, 72].
   Poll developed an approximate method for evaluating the spanwise momentum-thickness
Reynolds number, Reθz , at the leading edge of an infinitely swept cylinder [144]. This
criterion was modified by Campbell and Lynde [16] to use the leading-edge radius, resulting
in the following:

                                           p
                               Reθz = 0.404 (Re · rle · tan Λ),                             (1.5)

<!-- ===== PDF page 27 ===== -->

1.3. TRANSITION PREDICTION AND MODELLING METHODS                                               7

where Re is the Reynolds number based on the local chord, rle is the local leading-edge
radius, and Λ is the local wing leading-edge sweep angle.

1.3    Transition Prediction and Modelling Methods

The highest-fidelity approach for predicting the onset of boundary-layer transition involves
solving the Navier-Stokes equations directly using direct numerical simulation (DNS). DNS
solves the time-dependent Navier-Stokes equations directly and is able to produce all of the
mechanisms of boundary-layer transition, such as laminar flow breakdown, the development
of turbulent spots, and the transition to fully turbulent flow [200]. However, the number
of floating-point operations required to complete a simulation, which is proportional to the
number of mesh points and time steps, scales significantly with the Reynolds number, which
can quickly lead to impractical computational costs at higher Reynolds numbers. Large-
eddy simulation (LES) was developed in order to reduce the computational cost associated
with DNS simulations. Similar to DNS, LES computations solve the Navier-Stokes equations
directly. However, while the large-scale turbulent eddies are fully resolved, the small-scale
eddies are modelled using an eddy viscosity approach [169]. Although more computationally
efficient, the location of transition has been shown to be dependent on the Smagorinsky
constant specified to calibrate the subgrid eddy viscosity [56]. Germano implemented a
method to determine this constant locally [56]; however, work performed by Michelassi et
al. [119] demonstrated noticeable differences in the quantitative comparison of DNS and
LES simulations for transitional flows. As with DNS simulations, LES simulations incur
substantial computational costs for higher Reynolds number applications.
    Flow solvers based on the Reynolds-averaged Navier-Stokes (RANS) equations provide an
excellent balance between accuracy, robustness, and efficiency, making them a suitable choice
for practical engineering applications involving high Reynolds number flow with limited
amounts of separated flow. However, the turbulence models typically used in RANS solvers
do not have the stand-alone capability to predict boundary-layer transition. In order to
predict transition, one must apply a transition criterion. Reviews by Arnal et al. [8, 11] and
Pasquale and Rona [134] describe the advantages and disadvantages of several approaches
for coupling transition prediction criteria into RANS solvers. The most widely implemented
transition prediction strategies for practical engineering applications are based on either
the eN criterion, which is based on stability analysis, empirical transition onset criteria,
or on transport-equation-based transition models. The current review will focus on the
first and the last of these three approaches. Although there has been considerable work
applied to the development of empirical transition onset criteria [61, 1, 113, 183, 9, 10, 137],
the focus of recent work in the literature has been on the integration of these criteria in
transport-equation-based transition models to simplify their implementation in modern CFD
codes [92, 29, 133, 180].

<!-- ===== PDF page 28 ===== -->

8                                                                      CHAPTER 1. INTRODUCTION

1.3.1   eN Criterion

Approximations can be introduced to the Navier-Stokes equations to recover boundary-layer
stability equations, such as the parabolized stability equations (PSE) and the linear stability
equations derived from linear stability theory (LST), including the Orr-Sommerfeld equa-
tion, which provide significantly reduced computational cost relative to DNS. These stability
equations can be coupled with the eN method [170, 187], which tracks the most amplified
disturbances in the boundary layer from the stability equations, and, using empirical correla-
tions to evaluate the critical N-factor, can be used to identify the location of transition. The
N-factor curves produced by stability analysis can also be approximated using simplified eN
envelope methods, developed by Drela and Giles [44], and used by Rashad and Zingg [150]
and Mayda [112], depending on the level of fidelity required in a simulation.
   A boundary-layer code can be used to increase the accuracy in determining the boundary-
layer edge, which is required to satisfy the requirements of the stability analysis [177, 86].
However, with sufficient grid resolution, it is possible to define the boundary-layer edge and
extract the boundary-layer properties based on the RANS solution accurately without using
a boundary-layer solver [86, 150]. The influence of the grid resolution on the stability analysis
was studied in detail by Krimmelbein and Radespiel [86].
   Transition models based on the eN method are able to accurately predict natural- and
separation-induced transition. In addition, correlations for crossflow instabilities (such as the
C1 criterion [9]) have been successfully combined with these approaches, with experimental
validation demonstrating accurate transition prediction on transonic swept wings in three
dimensions [88, 89, 121, 178, 166]. However, eN methods require non-local boundary-layer
quantities that are not directly accessible in RANS-based computational fluid dynamics
(CFD) codes [165]. In the context of the current work, the requirement that the transition
model is formulated locally is significant. Modern CFD codes utilize domain decomposition
and parallel computation to increase the efficiency of simulations. Non-local operations, such
as integrating quantities over the boundary layer, introduce difficulties to this process. In
addition, non-local operations complicate the implementation in a discrete-adjoint gradient-
based optimization algorithm.

1.3.2   Transport-Equation-Based Transition Models

Of the many transport-equation-based transition models developed for RANS-based turbu-
                                                                                       ˜ θt
lence models, two which satisfy the requirement of a fully local formulation are the γ-Re
local correlation-based transition model developed by Langtry and Menter (LM2009) [92]
and the amplification factor transport (AFT) transition model developed by Coder and
Maughmer [30, 29] and later revised by Coder [26, 176, 31, 28].
   The AFT transition model is based on the approximate envelope method for TS wave
growth, developed by Drela and Giles [44], and is able to predict natural and separation-

<!-- ===== PDF page 29 ===== -->

1.3. TRANSITION PREDICTION AND MODELLING METHODS                                                9

induced transition. The model consists of a single transport equation for the envelope ampli-
fication factor, ñ, with a locally derived approximation for the boundary-layer shape factor.
The transition model was originally coupled directly to the one-equation Spalart-Allmaras
(SA) eddy viscosity model [174], resulting in a two-equation coupled system [30, 29, 26, 176].
The most recent model [31] is coupled to an additional, intermittency equation, inspired by
the work of Langtry and Menter [92], resulting in a three-equation coupled system. Re-
cent work by Xu et al. [193] extended the AFT model to include crossflow effects, using
the Falkner-Skan and Cooke equations to develop correlations for local approximations of
non-local boundary-layer quantities, with promising results. This work involved coupling the
model to the shear stress transport (SST) turbulence model [117], resulting in a four-equation
coupled system.
   The transition model formulated by Langtry and Menter [92] is composed of two transport
equations. The first is an equation for the intermittency, γ, which is used to trigger the
turbulence transition process by controlling the level of turbulence in the boundary layer.
The second is an equation for the transition momentum-thickness Reynolds number, Re      ˜ θt ,
which in the freestream is equal to the value produced by the empirical correlations, Reθt ,
and is diffused into the boundary layer using a standard diffusion term. A one-equation
variant of the model was recently developed by Menter et al. [118]. Of the various crossflow
transition models developed for the γ-Re ˜ θt transition model, two maintain a fully local
formulation, methods based on the local C1 criterion [25, 58, 59, 195] and local helicity
approaches [125, 93, 59]. The local C1 approach [9] was developed using solutions to the
Falkner-Skan and Cooke equations and is therefore only expected to work for wing-like
geometries with high aspect ratios [58, 59], unlike methods based on the local helicity.
    Crossflow instability correlations based on the local helicity were initially developed by
Müller and Herbst [125], with more recent correlations, which include a framework for in-
troducing roughness effects, developed by Langtry et al. (LM2015) [93]. Grabe and Krum-
bein [59] conducted simulations using correlations based on both the C1 and helicity ap-
proaches. The C1 approach yielded more accurate results in terms of the comparison of
predicted and measured transition locations on infinite swept wing flows, while the local he-
licity approach produced more accurate transition locations for non-wing-like geometries [59].
   The local helicity-based correlations developed by Müller and Herbst [125] and Langtry
et al. [93] both consist of an additional source term in the Re ˜ θt equation that acts as a sink
to reduce the value of Re ˜ θt in regions of significant crossflow, triggering transition further
upstream. While the crossflow correlations of Müller and Herbst [125] and Langtry et al. [93]
are similar, the latter represents a more thoroughly validated model [188]. Furthermore, the
ability of the Langtry model to include roughness effects expands the applicability of the
model.
   The original γ-Re ˜ θt LM2009 model developed by Langtry and Menter [92] was coupled
with the two-equation SST turbulence model [117]. Subsequent extensions of this model

<!-- ===== PDF page 30 ===== -->

10                                                                    CHAPTER 1. INTRODUCTION

to include crossflow instabilities with correlations based on the local helicity by Grabe and
Krumbein [59], Müller and Herbst [125], and Langtry et al. [93] were also coupled to the SST
model, while Nie et al. [126] investigated the performance of the Langtry-Menter model with
C1 and helicity-based crossflow correlations coupled with a Reynolds stress model. Relative
to multi-equation turbulence models, the one-equation SA turbulence model [174] provides
reduced computational cost, particularly for implicit solvers where the cost of solving the
linear system grows significantly with the number of equations. In addition, the model
constants in the SA turbulence model [174] are tuned for external aerodynamic flows, and
the model is shown to give good performance in boundary layers with adverse pressure
gradients [190].
   Medida and Baeder [116] modified the γ-Re     ˜ θt model to be compatible with the SA tur-
bulence model, and developed a new set of correlations validated using several test cases,
designating this model the γ-Re  ˜ θt -SA transition model. Recent work demonstrated the γ-
 ˜ θt -SA model coupled with the crossflow correlations developed by Müller and Herbst [125]
Re
was able to predict transition for several three-dimensional test cases, including a swept
wing, requiring only a slight modification [76]. However, their current formulation contains
a function, Gonset , which requires non-local flow information. More recently, Carnes and
Coder [18] successfully coupled the algebraic helicity-based crossflow correlations developed
by Langtry et al. [93] to the AFT model, with the results comparing well to experiment.
   Two transport-equation-based transition models, each consisting of four transport equa-
tions that act to create a local framework for evaluating the non-local Arnal-Habiballah-
Delcourt (AHD) transition criterion [9], have recently been developed and validated for low
turbulence intensity flows over a large range of Reynolds and Mach numbers [133, 180, 181].
Both models are coupled to the two-equation SST turbulence model [117], resulting in six-
equation coupled systems.

1.3.3   Effects of Flow Compressibility

The location of boundary-layer transition from laminar to turbulent flow can have a signif-
icant impact on the aerodynamic performance of a transonic wing. An upstream laminar
boundary layer directly influences shock location and strength, while laminar separation
bubbles can lead to severe drag penalties, stability and control issues, and adverse stall char-
acteristics. In order to design wings optimized for low drag for commercial transport aircraft,
which typically fly at transonic Mach numbers and high Reynolds numbers, it is necessary
to be able to accurately and efficiently predict boundary-layer transition in transonic flow
regimes.
   One of the first investigations into the effects of flow compressibility on the stability
of laminar boundary layers was performed by Lees and Lin who developed results for the
generalized inflection point and for the effects of wall heating and cooling in subsonic and
transonic flow [97]. More recently, the stabilizing effect of flow compressibility on laminar

<!-- ===== PDF page 31 ===== -->

1.4. SHAPE OPTIMIZATION WITH BOUNDARY-LAYER TRANSITION                                        11

boundary layers was investigated first by Mack [106] and then by Arnal [5, 7]. Due to the
lack of experimental data, these effects are best demonstrated using stability theory [186].
By applying the eN method with linear stability theory to adiabatic flat plates, Mack [106]
and Arnal [5] observed that compressibility has a strong stabilizing effect in transonic flow,
before becoming destabilizing from Mach 2 to 3.5, and stabilizing again for hypersonic flow.
   The empirical correlations of Langtry and Menter (LM2009) [92] and Langtry et al.
(LM2015) [93] were developed based on results from subsonic experiments and stability anal-
ysis to accurately predict two- and three-dimensional transition mechanisms, respectively.
However, recent studies have demonstrated that models based on the LM2009 and LM2015
correlations severely under-predict the extent of laminar flow when applied to transonic test
cases [189, 50], such as the NASA CRM-NLF [17, 104, 105]. This behaviour was foreseen by
Arnal who stated, “simple transition criteria developed for low-speed flows cannot be used
with confidence in configurations where compressibility effects become significant” [7].
   The non-local AHD transition criterion [9], which was developed as a simpler alternative
to the eN method [170, 187] based on stability analysis of Falkner-Skan attached self-similar
boundary-layer velocity profiles, has been extended to include the effects of compressibility on
TS instabilities. The new criterion was developed based on stability analysis of compressible
local similarity solutions first up to Mach 1.6 [10], and more recently up to Mach 4 [137].
Both transport-equation-based transition models based on the AHD criterion include this
correction [133, 181]. In addition, a compressibility correction for a transport-equation-based
transition model based on the Langtry and Menter model [92] was developed by Fers [50] in
order to incorporate the effects of compressibility on laminar boundary-layer velocity profiles
and therefore the vorticity Reynolds number in the transition criterion.
   Crossflow instabilities are dominated by the properties of the inflection point, which
are not affected by flow compressibility to the same degree as viscous instabilities, such as
TS waves [39, 6]. Malik et al. [108] developed a compressibility correction for a crossflow
Reynolds number criterion. This correction was used by Kroo and Sturdza [87] to develop
a compressible crossflow Reynolds number criterion for the design of laminar supersonic
swept wings, which was subsequently used by Lee and Jameson [96] to perform aerodynamic
shape optimization of transonic NLF wings. More recently, the C1 criterion [9], a crossflow
criterion developed based on stability analysis of the solutions of the Falkner-Skan and Cooke
equations, was extended to take into account compressibility effects using this correction [10].
A similar correction was also used by Xu et al. [194] to extend a local model for stationary
crossflow instabilities to transonic flows.

1.4    Shape Optimization with Boundary-Layer Transition

Although there are several examples in the literature of researchers combining transition
prediction with aerodynamic optimization using lower-fidelity models [43, 42, 62, 145, 87, 4],

<!-- ===== PDF page 32 ===== -->

12                                                                    CHAPTER 1. INTRODUCTION

or with gradient-free or costly finite-difference gradient approximations [199, 159, 198, 65,
185, 73, 162], the focus of the current review will be limited to methods using RANS-based
analysis with adjoint-based shape optimization algorithms.
   The work by Driver and Zingg [46] represents an early example of the use of RANS-
based optimization to investigate the design of NLF airfoils. The coupled inviscid-viscous
boundary-layer solver MSES [44] was used to calculate transition locations, which were
enforced using the trip terms in the SA turbulence model [174]. Through the application
of the resulting algorithm to single-point optimizations, Driver and Zingg [46] produced
NLF airfoil designs with laminar-flow design characteristics similar to those designed by
Liebeck [101] and Zingg [202].
   Lee and Jameson [96] coupled a boundary-layer solver and database eN method, developed
by Kroo and Sturdza [87] for the design of supersonic laminar flow wings, to the Baldwin-
Lomax algebraic turbulence model [12] in a RANS-based solver and discrete-adjoint gradient-
based optimization algorithm. Although their work demonstrated the significance of laminar
flow on NLF wing design, the gradients did not include transition prediction. Instead,
the optimizations were focused on reducing wave drag in order to increase aerodynamic
performance.
   Rashad and Zingg [150] performed single- and multi-point optimizations at subsonic and
transonic flight conditions using the simplified eN envelope method developed by Drela and
Giles [44] coupled to the SA turbulence model [174]. Rashad conducted a thorough investi-
gation of coupling strategies for the RANS-eN flow solver and demonstrated the importance
of providing a smooth ramp-up of the eddy viscosity, as well as using a tight tolerance for the
transition residual, in order to ensure a smooth design space [149, 150]. The coupled-adjoint
system is solved using preconditioned GMRES. Single-point optimizations demonstrate that
the resulting algorithm is capable of extending the laminar boundary layer significantly, with
multi-point optimizations producing NLF airfoil designs robust to changes in the flow condi-
tions and disturbance environments, where the latter are represented using a range of critical
N-factors. The drag minimizations performed by Rashad and Zingg [150] represent suitable
two-dimensional reference benchmarks for evaluating NLF optimization frameworks.
   Shi et al. [167] used a simplified eN method based on a database of LST results in a
Jacobian-free discrete-adjoint optimization algorithm. Similar to the work by Rashad [149]
and Rashad and Zingg [150], Shi et al. emphasized the importance of using a smooth in-
termittency function to couple the simplified eN method with the RANS solutions, as well
as a tight transition residual tolerance. Shi et al. investigated single-point airfoil optimiza-
tions at the Cessna 172R design conditions initially proposed by Rashad and Zingg [150],
producing an optimized NLF design similar to that developed by Rashad and Zingg, and
multi-point optimizations of an airfoil at flight conditions representative of the HondaJet
aircraft [54]. Recently, their work was extended to infinite swept wing optimizations using
the C1 criterion [168].

<!-- ===== PDF page 33 ===== -->

1.4. SHAPE OPTIMIZATION WITH BOUNDARY-LAYER TRANSITION                                     13

   Zhu and Qin [201] applied the simplified eN envelope method by Drela and Giles [44],
extended with the crossflow criterion developed by Kroo and Sturdza [87], and using Poll’s
criterion [144] as a constraint to prevent attachment-line transition, to the optimization
of infinite swept wings using a discrete-adjoint optimization algorithm. Similar to Lee
and Jameson [96], the transition locations from stability analysis were enforced using the
Baldwin-Lomax algebraic turbulence model [12]. Drag-minimizations of infinite swept wings
at transonic flight conditions were performed both with and without a shock-control bump.
The shock-control bump reduces wave drag, enabling the optimization algorithm to reduce
sweep, which attenuates the crossflow instabilities, moving the transition front downstream.
A systems-level benefit assessment was performed demonstrating that a low-sweep NLF wing
design with shock-control bumps applied to a narrow-body commercial aircraft can reduce
fuel consumption by 11.1% [49].
   As previously discussed, methods based on stability analysis require non-local boundary-
layer information and a large infrastructure to apply the code, usually consisting of a
boundary-layer solver and a module for solving the stability equations, which can in-
crease computational cost and complicate their implementation in highly parallel solvers.
Methods based on local transition criteria, such as the Langtry-Menter transition model
(LM2009) [92], can be more easily integrated into a RANS-based flow solver and gradient-
based optimization algorithm.
   Khayatzadeh and Nadarajah [81, 82] coupled the Langtry-Menter LM2009 transition
model [92] with a RANS-based flow solver and discrete-adjoint optimization algorithm.
Through the application of the resulting algorithm to subsonic lift-constrained drag-
minimization and lift-to-drag maximization optimization studies, they demonstrated the
importance of fully coupling the transition model equations to the adjoint system, and suc-
cessfully produced designs that delayed the onset of transition.
   Halila et al. [64] integrated a smooth version of the AFT transition model [29], coupled
to the SA turbulence model, with an approximate Newton-Krylov flow solver and discrete-
adjoint optimization algorithm. The resulting algorithm was applied to single- and multi-
point subsonic and single-point transonic airfoil drag-minimization studies.
   Yang and Mavriplis [196] integrated the AFT transition model [29] with a Newton-Krylov
flow solver and discrete-adjoint framework and applied the algorithm to two-dimensional
subsonic optimizations. Following this work, Mavripilis et al. [111] applied the one-equation
intermittency transition model developed by Menter et al. [118], coupled to the SA turbulence
model [174], in a discrete-adjoint optimization algorithm to a slotted transonic truss-braced
wing (TTBW) geometry.
   In general, optimized designs resulting from algorithms combining gradient-based op-
timization with local transport-equation-based transition models tend to underperform
NLF designs produced using stability-analysis approaches coupled with gradient-based op-
timization algorithms. Specifically, gradient-based optimization algorithms with transport-

<!-- ===== PDF page 34 ===== -->

14                                                                    CHAPTER 1. INTRODUCTION

equation-based transition models struggle to produce designs with significantly increased
regions of laminar flow [81, 82, 196, 77], especially at higher Mach and Reynolds num-
bers [64]. This is particularly evident when comparing with the designs produced by Rashad
and Zingg [150] that achieve significant regions of laminar flow at high Mach numbers and
flight-scale Reynolds numbers.

1.5       Objectives

The objective of this work is to develop and investigate a high-fidelity aerodynamic shape
optimization framework based on the RANS equations that includes the effects of boundary-
layer transition which can be applied to the design of subsonic and transonic NLF wings. To
ensure a realistic design, the relevant two- and three-dimensional transition mechanisms must
be included. For a swept wing, these include transition due to TS and crossflow instabilities,
which are sensitive to Reynolds number, turbulence intensity, pressure gradient, surface
roughness, and Mach number. Although relevant, Görtler and attachment-line instabilities
can be prevented through effectively constraining the design space, which will be the focus of
future work. To enable efficient design space exploration, the transition prediction strategy
must be compatible with the discrete-adjoint method. To facilitate this, methods for robust
and efficient free-transition function and gradient evaluation must be developed.
   The following subgoals are necessary in order to achieve the primary objective:

     1. Extend a three-dimensional RANS solver with a local transition prediction model.

     2. Develop a framework for achieving deep and efficient iterative convergence.

     3. Validate the model for a range of subsonic and transonic flow regimes.

     4. Integrate the model in a discrete-adjoint gradient-based optimization algorithm.

     5. Apply the framework to drag-minimization studies over a range of design conditions.

1.6       Thesis Outline

The remainder of the thesis is organized as follows: Chapter 2 presents the turbulence
and transition model equations. Chapter 3 presents the flow solution strategy, including
modifications introduced to stabilize the implicit treatment of source terms in the turbulence
and transition model equations, and the integration in the discrete-adjoint gradient-based
optimization algorithm. Validation studies using two- and three-dimensional transition test
cases in both subsonic and transonic flow regimes are presented in Chapter 4. Free-transition
drag-minimization studies across a range of design conditions are presented in Chapter 5. A
summary and conclusions are presented in Chapter 6.

<!-- ===== PDF page 35 ===== -->

Chapter 2

Transition Model Equations

This chapter presents the turbulence and transition model equations. These include the
SA turbulence model, the SA-LM2015 and SA-sLM2015 transition models developed in the
current work [139], as well as the compressibility corrections developed and applied to extend
the LM2009 [92] and LM2015 [93] empirical correlations to transonic flow regimes [142]. The
curvilinear coordinate transformation for the RANS equations and the SA turbulence model
can be found in the work produced by Osusky and Zingg [129], which provides the framework
for the current work.

2.1    SA-LM2015 Transition Model

The γ-Re˜ θt LM2015 model with helicity-based crossflow correlations [93] is modified and cou-
pled to the SA turbulence model [174], while maintaining a fully local formulation, reducing
the computational cost of the original four-equation system based on the SST turbulence
model. The LM2015 model is decoupled from the SST turbulence model using modifications
to the eddy viscosity ratio, RT , and the Fwake function. Schücker [164] suggests replacing
the specific turbulence dissipation rate, ω, in the original Fwake function with velocity gra-
dients. This approach is adopted in the current work using the strain-rate magnitude. A
simplification of the Fθt function is introduced, which reduces the amount of non-smooth
minimum/maximum operators present in the model. An adjustment to the Fonset1 term is
applied to recalibrate the model, and the strain-rate magnitude in the original intermittency
production term is replaced by vorticity magnitude to make the model more stable near
laminar separation bubbles. The separation-induced transition modification developed by
Langtry [91] is neglected, as the delay in the growth of turbulence after laminar separation
bubbles produced by the SST turbulence model is not reproduced by the SA turbulence
model, which is demonstrated by the results in Chapter 4.
   The remaining terms and coefficients remain unchanged from the original model. The
modified LM2015 model with helicity-based crossflow correlations coupled to the SA turbu-
lence model, designated SA-LM2015, is presented as follows.

<!-- ===== PDF page 36 ===== -->

16                                                                       CHAPTER 2. TRANSITION MODEL EQUATIONS

2.1.1       Momentum-Thickness Reynolds Number Transport Equation

The transport equation for the momentum-thickness Reynolds number, non-dimensionalized
by freestream quantities and written in non-conservative form to be consistent with the SA
turbulence model, is defined as,

                        ˜ θt
                      ∂ Re          ˜ θt
                                  ∂ Re                    1 ∂
                                                                                     ˜ θt 
                                                                                    ∂ Re
                             + uj        = Pθt + Dscf +               σθt (ν + νt )          ,           (2.1)
                        ∂t         ∂xj                   Re∞ ∂xj                     ∂xj
                                            cθt         ˜ θt )(1 − Fθt ),
                                      Pθt =     (Reθt − Re                                               (2.2)
                                             t
                                                             d 4
                                             Fθt = Fwake e−( δ ) ,                                       (2.3)
                               50dΩ                 15                   ˜ θt µ 1
                                                                        Re
                           δ=        δBL , δBL = θBL , θBL =                          ,                  (2.4)
                                  U                  2                    ρU Re∞
                                           −ReS ×10−6              ρd2 S
                                  Fwake = e           , ReS =            Re∞ ,                           (2.5)
                                                                    µ
                                                cθt = 0.03, σθt = 2.0,                                   (2.6)
                                                         500µ 1
                                                    t=            ,                                      (2.7)
                                                         ρU 2 Re∞

where Re∞ = ρ∞µa∞∞ l , ‘a’ is the sound speed, ‘l’ is the characteristic (reference) length, and
‘∞’ denotes a freestream quantity. The non-dimensionalization procedure for the transition
model equations is consistent with the procedure for the Navier-Stokes equations, which is
described in [147] and is used, for example, in the well-known OVERFLOW flow solver1 .
 ˜ θt refers to the transported quantity of the momentum-thickness Reynolds number, while
Re
Reθt is the value calculated using the empirical correlations.
      The natural and bypass transition empirical correlations are defined as follows:
                                                                 
                            1173.51 − 589.428T u +       0.2196
                                                                    F (λθ ) T u ≤ 1.3,
                           
              ρU θ                                          T u2
       Reθt =      Re∞ =                                                                                 (2.8)
               µ           331.50[T u − 0.5658]−0.671 F (λ )                 T u > 1.3,
                           
                                                                 θ
                           
                           1 − [−12.986λ − 123.66λ2 − 405.689λ3 ]e−[ T1.5u ]1.5 λ ≤ 0,
                                          θ                θ                θ             θ
                 F (λθ ) =                               Tu
                                                                                                         (2.9)
                           1 + 0.275[1 − e[−35λθ ] ]e−[ 0.5 ]                           λθ > 0,

where the local pressure gradient parameter is given by,

                                     ρθ2 dU                            1
                                 λθ =       Re∞ , U = (u2 + v 2 + w2 ) 2 ,                              (2.10)
                                      µ ds
                                                          
                                dU       u dU       v dU      w dU
                                   =           +           +             ,                              (2.11)
                                ds       U dx       U dy      U dz
                                                                            
                                dU     2    2    2 − 21   du    dv       dw
                                   = (u + v + w ) · u + v          +w         ,                         (2.12)
                                dx                        dx    dx        dx
     1 https://overflow.larc.nasa.gov/home/users-manual-for-overflow-2-3/, accessed August 2020

<!-- ===== PDF page 37 ===== -->

2.1. SA-LM2015 TRANSITION MODEL                                                                     17
                                                         
                         dU     2   2   2 − 21   du dv dw
                            = (u + v + w ) · u + v + w      ,                                   (2.13)
                         dy                      dy dy dy
                                                         
                         dU     2   2   2 − 21   du dv dw
                            = (u + v + w ) · u + v + w     .                                    (2.14)
                         dz                      dz dz dz

This formulation of the local pressure gradient parameter, which uses the streamwise velocity
gradient formed using Cartesian velocity gradients, is used in the original model by Langtry
and Menter [92] in combination with a blending function, Fθt (Equation 2.3), to disable it in
the boundary layer where it is no longer accurate. An alternative formulation was investi-
gated that uses the local Cartesian pressure gradients (Equations 2.75–2.79). However, the
original correlations by Langtry and Menter were calibrated using the velocity-based pres-
sure gradient. Although the two formulations produce similar results for non-zero pressure
gradients, the velocity-based pressure gradient is used in the current work to evaluate the
empirical correlations in order to be consistent with Langtry and Menter [92]. It is important
to note that this parameter is not Galilean invariant due to the use of the velocity magni-
tude. While this is not an issue for the current work, which is focused on fixed reference
frames with no moving walls, a Galilean invariant pressure gradient parameter similar to
that introduced by Menter et al. [118] could be adopted in the future.
    The local turbulence intensity is not directly available when using the SA turbulence
model. However, Suluksana et al. [182] demonstrated that it can be redundant to use the
local value of the turbulence intensity, T u, as well as the pressure gradient parameter λθ , since
the effects of pressure gradient are reduced by the decay of the local turbulence intensity.
A favourable pressure gradient accelerates the flow, leading to a decrease in the turbulence
intensity, and therefore an increase in the momentum-thickness Reynolds number, while
an adverse pressure gradient leads to an increase in turbulence intensity. The freestream
turbulence intensity value, T u∞ , is used throughout the domain in the current work, which
is consistent with the approach taken by Medida and Baeder [116].
  The helicity-based stationary crossflow momentum-thickness Reynolds number, Reθt,scf ,
empirical correlation developed by Langtry et al. is defined below [93]:

                       cθt                           ˜ θt , 0)(Fθt ), ccrossflow = 0.6,
                Dscf =     ccrossflow min(Reθt,scf − Re                                         (2.15)
                        t
                       ρ( U )θt
           Reθt,scf = 0.82 Re∞
                            µ
                                      
                                       h                        +                 −
                     = −35.088 ln          + 319.51 + f (∆Hcrossflow    ) − f (∆Hcrossflow ),   (2.16)
                                       θt
                                                                     µt
         ∆Hcrossflow = Hcrossflow (1.0 + min(RT , 0.4)), RT = ,                                 (2.17)
                                                                     µ
           +
         ∆Hcrossflow = max(0.1066 − ∆Hcrossflow , 0),                                           (2.18)
          +                    +                     +
     f (∆Hcrossflow ) = 6200(∆Hcrossflow ) + 50000(∆Hcrossflow )2 ,                             (2.19)
           −
         ∆Hcrossflow = max(−(0.1066 − ∆Hcrossflow ), 0),                                        (2.20)

<!-- ===== PDF page 38 ===== -->

18                                                            CHAPTER 2. TRANSITION MODEL EQUATIONS

                                       −        
          −                          ∆Hcrossflow
     f (∆Hcrossflow ) = 75tanh                    .                                          (2.21)
                                      0.0125

Helicity, alternatively known as the streamwise vorticity, Ωstreamwise , is the magnitude of the
                                                                ~ , and vorticity vector, Ω:
scalar product of the local values of the unit velocity vector, U                         ~
                                                                            
                   ~ = √        u                v                w
                   U                     ,√               ,√                   ,          (2.22)
                           u2 + v 2 + w2    u2 + v 2 + w2    u2 + v 2 + w2
                                                       
                   ~      ∂w ∂v ∂u ∂w ∂v ∂u
                   Ω=        −     ,    −     ,    −      ,                               (2.23)
                          ∂y    ∂z ∂z      ∂x ∂x ∂y

and can be represented using the non-dimensional crossflow strength, Hcrossflow , as,

                                                       ~ · Ω|,
                                        Ωstreamwise = |U   ~                                 (2.24)
                                                      dΩstreamwise
                                         Hcrossflow =              .                         (2.25)
                                                           U
The empirical correlations for Reθt (Equation 2.8) and Reθt,scf (Equation 2.16) are implicit
functions of the momentum thickness, θ, and the transition onset momentum thickness, θt ,
which are solved for in the current work using Newton’s method. It is important to note
that the use of the velocity magnitude in the local helicity formulation once again violates
Galilean invariance. The development of a Galilean invariant crossflow correlation could be
a focus of future work.

2.1.2   Intermittency Transport Equation

The transport equation for intermittency, also non-dimensionalized and written in non-
conservative form, is,
                                                                     
                    ∂γ      ∂γ                    1 ∂            νt ∂γ
                       + uj     = Pγ − Eγ +                ν+             ,      (2.26)
                    ∂t      ∂xj                 Re∞ ∂xj          σf ∂xj
                                                     √
                            Pγ = ca1 Flength Fonset Ω γ(1 − ce1 γ),              (2.27)
                                 ce1 = 1.0,     ca1 = 2.0,        σf = 1.0,                  (2.28)
                                      Eγ = ca2 Fturb Ωγ(ce2 γ − 1),                          (2.29)
                                                         R
                                                       −( 4T )4
                                             Fturb = e            ,                          (2.30)
                                         ce2 = 50,    ca2 = 0.06,                            (2.31)

where the Fonset function, which activates the transition process, is given by,
                                          p
                             Fonset = max( Fonset2 − Fonset3 , 0),                           (2.32)
                                              3 
                                                RT
                            Fonset3 = max 1 −         ,0 ,                                   (2.33)
                                                2.5
                            Fonset2 = min(max(Fonset1 , Fonset1 ), 2),                       (2.34)

<!-- ===== PDF page 39 ===== -->

2.1. SA-LM2015 TRANSITION MODEL                                                                  19

                                            ReS
                              Fonset1 =           ,                                          (2.35)
                                          2.6Reθc
and the Flength and Reθc empirical correlations are as follows:
              
              
               398.189 · 10−1 + (−119.270 · 10−4 )Re  ˜ θt +
              
                                −6 ˜ 2
              
              
              
              
               (−132.567   · 10  )Reθt                                ˜ θt < 400,
                                                                      Re
              
              
              263.404 + (−123.939 · 10−2 )Re
                                               ˜ θt +
    Flength =                                                                                (2.36)
              (194.548 · 10−5 )Re
                                 ˜ 2θt + (−101.695 · 10−8 )Re˜ 3θt 400 ≤ Re  ˜ θt < 596,
              
              
                        ˜ θt − 596.0) · 3 · 10−4                              ˜ θt < 1200,
              
                0.5 − (Re                                             596 ≤ Re
              
              
              
              
              
                                                                                ˜ θt ,
              
                0.3188                                                1200 ≤ Re
              
              
              
                ˜ θt − (396.035 · 10−2 + (−120.656 · 10−4 )Re
                Re                                              ˜ θt +
              
                                  ˜ 2θt + (−696.506 · 10−9 )Re˜ 3θt +
              
              (868.230 · 10−6 )Re
              
     Reθc =                                                                                  (2.37)
                                  ˜4)
                (174.105 · 10−12 )Re                                       ˜ θt ≤ 1870,
                                                                          Re
              
              
                                     θt
              ˜
              
                                           ˜ θt − 1870.0))                 ˜ θt > 1870.
                Reθt − (593.11 + 0.482(Re                                 Re

2.1.3   Coupling to SA Turbulence Model

The SA-neg-noft2 variant of the SA turbulence model [174, 3] is used in the current work
and is presented in non-conservative, non-dimensional form below:
                                                                                     
            ∂ ν̃      ∂ ν̃                    1    ∂             ∂ ν̃         ∂ ν̃ ∂ ν̃
                 + uj      = P̃ν̃ − Dν̃ +               (ν + ν̃)        + cb2             , (2.38)
             ∂t       ∂xj                 σRe∞ ∂xj               ∂xj          ∂xi ∂xi
                                                              2
                                      cb1              cw1 fw ν̃
                              Pν̃ =       S̃ ν̃, Dν̃ =                .                     (2.39)
                                     Re∞               Re∞ d

Osusky [129] presents the non-dimensionalization procedure for the SA turbulence model,
which is consistent with that used for the Navier-Stokes equations. The transition model is
coupled to the production term, Pν̃ , through the following:

                                               P̃ν̃ = γPν̃ .                                 (2.40)

2.1.4   Initial and Boundary Conditions

The initial conditions for the turbulence and transition model equations for a free-transition
simulation are defined by,

                    ν̃initial = 0.1,   γinitial = 1,   ˜ θt,initial = Reθt (λθ = 0),
                                                       Re                                    (2.41)

while the boundary conditions are:

                      farfield:        ν̃ = 0.1,   γ = 1,      ˜ θt = Reθt (λθ = 0).
                                                               Re                            (2.42)

<!-- ===== PDF page 40 ===== -->

20                                                        CHAPTER 2. TRANSITION MODEL EQUATIONS

                                   ∂ ν̃         ∂γ           ˜ θt
                                                          ∂ Re
                   symmetry:            = 0,       = 0,           = 0.                     (2.43)
                                   ∂n           ∂n          ∂n
                                               ∂γ          ˜ θt
                                                        ∂ Re
                         wall:     ν̃ = 0,        = 0,          = 0.                       (2.44)
                                               ∂n         ∂n

2.2     SA-sLM2015 Transition Model

The LM2015 and SA-LM2015 models contain several non-smooth functions in the form of
min/max operators and conditional statements. These non-smooth functions, referred to
as simple kinks [13, 21], produce discontinuities in the Jacobian. This presents a problem
for both the Newton-Krylov solution algorithm and the discrete-adjoint gradient-based op-
timization algorithm. For example, these discontinuities can cause the nonlinear solver for
the flow solution to stall, preventing the deep convergence necessary for the discrete-adjoint
method. In the current work, this has been addressed through the application of smooth ap-
proximations to min/max operators and conditional statements. In addition, modifications
are introduced, such as new Fonset and Fturb functions, in order to reduce numerical stiffness.
   The modifications below describe the simplification and smoothing of the SA-LM2015
transition model, presented in the previous section. For the remainder of this work the
smooth model will be designated SA-sLM2015.

2.2.1   Minimum/Maximum Operators

Several examples of smooth approximations to minimum and maximum operators can be
found in the literature. Rivest [157] suggested using a generalized mean-value operator, such
as the p-mean:
                                             X n            p1
                                             1            p
                                   φp (x) ,        gi (x)        .                         (2.45)
                                             n i=1

As p approaches positive and negative infinity, we recover the maximum and minimum
operators, respectively. However, this approximation only holds if g = (g1 , ..., gn ) is a vector
of positive real numbers, which limits its applicability.
   Another approach is to use an exponential penalty function [135, 143, 98], which is similar
to the Kreisselmeier-Steinhauser function for constraint aggregation [84, 110, 78], defined as,

                                        log( ni=1 exp(pgi (x)))
                                            P
                              φp (x) ,                          .                           (2.46)
                                                   p

Again, as p approaches positive infinity, the function produces the maximum of gi , while
p approaching negative infinity produces the minimum. The advantage of this approach is
that g = (g1 , ..., gn ) can be a vector of any real numbers. A drawback is the presence of
log and exponential functions, which can increase the computational cost. This is especially
true as gi approaches values close to zero. This can introduce denormalized floating-point

<!-- ===== PDF page 41 ===== -->

2.2. SA-SLM2015 TRANSITION MODEL                                                                                  21

Algorithm 1 Smooth maximum/minimum function of two variables
 1: function φp (g1 (x), g2 (x))
 2:    pswitch ← 10−15                                                           . Cutoff value for proximity check
 3:    a ← max(g1 (x), g2 (x))                                                                   . Save larger value
 4:    b ← min(g1 (x), g2 (x))                                                                 . Save smaller value
 5:    if p > 0 then                                                                        . Determine maximum
 6:        if |a − b| > − log(|p|·p
                                 |p|
                                    switch )
                                             then                                                . Proximity switch
 7:            φp (x) ← a                                                                       . Return maximum
 8:            return
 9:        else
10:            φp (x) ← a + log(1+exp(p(b−a)))
                                         p                                . Calculate exponential penalty function
11:    else                                                                                 . Determine minimum
12:       if |a − b| > − log(|p|·p
                                |p|
                                   switch )
                                            then                                               . Proximity switch
13:           φp (x) ← b                                                                       . Return minimum
14:           return
15:       else
16:           φp (x) ← b + log(1+exp(p(a−b)))
                                        p                                 . Calculate exponential penalty function

numbers which significantly slow the algorithm. To prevent this, a proximity switch has been
developed (Equation 2.48) to activate the exponential penalty function only in the region
near the kink. In addition, when considering a two-variable array the exponential penalty
function can be rewritten as,

                                                    log(1 + exp(p(g2 (x) − g1 (x))))
                            φp (x) , g1 (x) +                                        .                      (2.47)
                                                                  p

In combination with the proximity switch, this form of the exponential penalty function
prevents both underflow and overflow. The proximity switch is given by,

                                                         log(|p| · pswitch )
                                          φswitch , −                        ,                              (2.48)
                                                                |p|

where pswitch represents a cutoff value for the proximity check, which for the current work
is set to a value of 10−15 . If the magnitude of the difference between the two variables,
|g1 (x) − g2 (x)|, is greater than φswitch , demonstrating that we are far from the discontinuity,
a simple min/max operator is used. As we approach the kink, if the difference between the
two variables is less than φswitch , the exponential penalty function is used.
    Algorithm 1 presents the complete method for approximating a minimum/maximum op-
erator, based on the value of the smoothing parameter, p, for two variables, g1 (x) and g2 (x),
using the exponential penalty functions including the proximity switch. Figure 2.1 illustrates
the effect of this smoothing by approximating the maximum of a variable and a constant
using various smoothing values, p. For large values of p, the region of smoothing is small,
and the sharp region is closely modelled. As p is reduced, the domain where the exponen-
tial penalty functions are active increases, increasing computational cost and producing a
smoother approximation of the kink region. Although adaptively updating the value of p has
been shown to help prevent ill-conditioned systems [143], all min/max operators present in

<!-- ===== PDF page 42 ===== -->

22                                                                 CHAPTER 2. TRANSITION MODEL EQUATIONS

                            0.8
                                             0.06
                                             0.04
                            0.6
                                             0.02
                            0.4
                                                    0.95       1   1.05

                            0.2

                              0.5                          1                    1.5

 Figure 2.1: max(0, x − 1) approximated using φp with a range of values for the smoothing parameter, p.

the SA-sLM2015 model are approximated using this algorithm with fixed magnitudes of p of
300, φ±300 , which was found through numerical experimentation to provide a good balance
between convergence properties and solution accuracy.

2.2.2   Intermittency Source Terms

Fonset Function

The Fonset function (Equations 2.32–2.35) couples the SA, intermittency, and momentum-
thickness Reynolds number equations. It was developed using several layers of minimum
and maximum operators to activate the production of intermittency in the boundary layer,
triggering the transition process. Solving the mean-flow, turbulence, and transition model
equations in a fully coupled manner results in the Fonset function being particularly prob-
lematic, and often the bottleneck in flow solver convergence. Rashad [149] investigated the
design space of several transition region models and noted that the smoothness of the ramp-
up of the eddy viscosity was found to be particularly important to the smoothness of the
design space. He also demonstrated the importance of achieving deep numerical convergence
of the transition model residual. Therefore, special care is taken in ensuring that the Fonset
function be as smooth as possible while maintaining the same characteristics of the original
function.
   The transition process is triggered when Fonset1 (Equation 2.35) exceeds unity, and is
sustained as the eddy viscosity ratio increases. To smooth and simplify this process, the
Fonset function is reformulated using a hyperbolic tangent function as follows:

                                        tanh(6(Fonset1 − 1.35)) + 1
                                  Fonset =                          ,                             (2.49)
                                        s           2
                                                    2  2
                                              ReS
                              Fonset1 =                + RT .                                     (2.50)
                                            2.6Reθc

<!-- ===== PDF page 43 ===== -->

2.2. SA-SLM2015 TRANSITION MODEL                                                                  23

                                               1.5

                           15                  0.5

                           10                      0
                                                   540   560   580   600   620   640

                                0            500                1000               1500

                         Figure 2.2: Original and smooth Flength functions.

Similar to the original, the new Fonset function is activated when Fonset1 exceeds unity; how-
ever, once a turbulent boundary layer has formed Fonset1 remains active. This increases the
stability of the model around laminar separation bubbles, where the solution can oscillate
as the bubble grows and shrinks, but prevents the smooth model from being able to predict
relaminarization, which is a feature of the original model. For flight conditions typical of
transport aircraft, we do not expect relaminarization to occur; however, the reader should
be aware of this limitation.

Fturb Function

In the SA-LM2015 model the Fturb function (Equation 2.30) remains active until a threshold
eddy viscosity ratio is reached, approximately a value of 4. Until this threshold value is
reached, the destruction term remains active and interferes with the production term, slowing
convergence and often producing oscillations. To remedy this, a new function is introduced
that deactivates the destruction term when Fonset becomes active, and rapidly deactivates
Fturb as turbulent eddy viscosity is produced. It is given as follows:

                                    Fturb = (1 − Fonset ) exp(−RT ).                          (2.51)

Flength Correlation

The original Flength function (Equation 2.36) consists of conditional statements and is non-
differentiable at the boundaries. In the current work, this is avoided by approximating
Flength using a modified, flexible version of a sigmoid function known as a generalized logistic
function [155]. The new Flength function is defined as follows and plotted in Figure 2.2:

                                                                ˜ θt − 596))
                                         44 − (0.50 − 3 · 10−4 (Re
                      Flength = 44 −                                       1              ,   (2.52)
                                                (1 + Flength1 ) 6
                                                ˜ θt − 460)).
                      Flength1 = exp(−3 · 10−2 (Re                                            (2.53)

<!-- ===== PDF page 44 ===== -->

24                                                                   CHAPTER 2. TRANSITION MODEL EQUATIONS

                                           480         500     520

                                 0               500         1000     1500    2000

                          Figure 2.3: Original and smooth Reθc functions.

Reθc Correlation

In previous work, approximations of the critical momentum-thickness Reynolds number,
Reθc , as linear functions of Re˜ θt are introduced [116, 164, 76]. In the current work, a
                                                                             ˜ θt behaviour of
sinusoidal approximation is developed in order to better represent the lower Re
the original correlation (Equation 2.37). The critical momentum-thickness Reynolds number
is simplified from the original version as follows:
                                                 ˜         
                                   ˜ θt + 24 sin Re θt
                        Reθc = 0.67Re                  + 0.5 + 14.                                  (2.54)

The original and smooth Reθc functions are plotted in Figure 2.3. For large values of Re ˜ θt ,
the new correlation produces smaller Reθc values relative to the original model; however, this
difference has been accounted for in the model calibration. Emphasis was placed on closely
approximating the original function for 300 ≤ Re  ˜ θt ≤ 1600, which is the range of values
typically encountered in a converged flow solution.

Pγ and Eγ Source Terms

The maximum value of vorticity magnitude in the intermittency production and destruction
terms, Pγ (Equation 2.27) and Eγ (Equation 2.29), is limited to improve stability. Vorticity
magnitude is included to obtain the correct dimensions for the source terms; however, it
was found that large values of vorticity, especially near laminar separation bubbles, can
lead to the solution oscillating with the transition location moving up and downstream in a
periodic manner. Limiting the maximum vorticity value in the intermittency source terms
restricts the relative magnitude of the sources, improving the numerical behaviour of the
model. The limiting is a function of the freestream Mach number in order to account for the
non-dimensionalization of velocity by the speed of sound, a∞ , and the freestream Reynolds
number to account for the scaling of vorticity in laminar boundary layers with respect to

<!-- ===== PDF page 45 ===== -->

2.2. SA-SLM2015 TRANSITION MODEL                                                                      25

                   (a) SA-LM2015                                      (b) SA-sLM2015

Figure 2.4: Original and smooth intermittency source terms, Sγ = Pγ − Eγ , plotted versus the transition
criterion (Equation 2.35), intermittency, γ, and the eddy viscosity ratio, RT , with Flength = Ω = 1.

the square root of the Reynolds number. The Reynolds number is multiplied by the Mach
number to recover the Reynolds number based on the freestream velocity. It is important
to note that while the presence of the Mach number is the result of the particular non-
dimensionalization used, the square root of the freestream Reynolds number based on the
freestream velocity is a physical scaling that is independent of the non-dimensionalization.
Furthermore, the vorticity limiting is introduced only to improve the numerical behaviour
of the model and does not aim to improve its predictive capability. The new source terms
are given as follows:
                                                  √     
                                                 M∞ M∞ Re∞ √
                Pγ = ca1 Flength Fonset φ−300 Ω,               γ(1 − ce1 γ),                     (2.55)
                                                √      
                                               M∞ M∞ Re∞
                    Eγ = ca2 Fturb φ−300 Ω,                γ(ce2 γ − 1),                         (2.56)

and are illustrated in Figure 2.4 alongside the original SA-LM2015 source terms. The total
source terms, Sγ = Pγ −Eγ , are plotted with fixed Flength and vorticity magnitude, Ω, values of

<!-- ===== PDF page 46 ===== -->

26                                                               CHAPTER 2. TRANSITION MODEL EQUATIONS

                               1.2

                               0.8                       1.01

                               0.6
                                                         0.99

                               0.4                              -1   0      1
                                                                                 -3
                                 -0.1        -0.05         0         0.05             0.1

          Figure 2.5: Original and smooth F (λθ ) functions for a turbulence intensity of 0.05%.

unity versus the transition criterion (Equation 2.35), intermittency, γ, and the eddy viscosity
ratio, RT . The reformulated intermittency source terms provide a smoother contour map,
while the results in Chapter 4 demonstrate that the predicted transition locations with these
modifications match well with the SA-LM2015 model and with experimental data for the
cases considered.

2.2.3    Momentum-Thickness Reynolds Number Empirical Correlations

F (λθ ) Correlation

Similar to the Flength function, the F (λθ ) correlations (Equation 2.9), which predict the effects
of pressure gradient on the transition onset momentum-thickness Reynolds number, Reθt ,
are conditional, depending on the sign of the pressure gradient parameter, λθ . Here F (λθ )
is reformulated using minimum and maximum operators, which can then be approximated
using exponential penalty functions:

                      F (λθ ) = φ−300 (F (λθ )2 , F (λθ )3 ),                                           (2.57)
                                                                                       −[ T u 1.5
                  F (λθ )3 = 1 − [−12.986λθ − 123.66λ2θ − 405.689λ3θ ]e                   1.5
                                                                                              ]
                                                                                                    ,   (2.58)
                  F (λθ )2 = φ300 (F (λθ )1 , 1),                                                       (2.59)
                                                                Tu
                  F (λθ )1 = 1 + 0.275[1 − e[−35λθ ] ]e−[ 0.5 ] .                                       (2.60)

The original and smooth F (λθ ) functions are plotted in Figure 2.5.

     ±
f (∆Hcrossflow ) Functions

          ±
The f (∆Hcrossflow ) functions either increase or decrease the stationary crossflow Reynolds
number, Reθt,scf , empirical correlation depending on the strength of the local helicity [93].
The original functions presented in Equations 2.18–2.21 contain two maximum operators.

<!-- ===== PDF page 47 ===== -->

2.3. COMPRESSIBILITY CORRECTIONS                                                            27

                           800                   50
                           600                       0

                           400                  -50
                                                     0.09 0.1 0.11 0.12 0.13

                          -200
                                 0    0.05     0.1        0.15       0.2       0.25

                                                          ±
                     Figure 2.6: Original and smooth f (∆Hcrossflow ) functions.

While these could be approximated using the exponential penalty functions presented pre-
viously, an approximation is introduced using a three-term Gaussian model to reduce com-
putational cost,
                                     3                            2 
                      ±
                                     X             ∆Hcrossflow − bi
                 f (∆Hcrossflow )=       ai exp −                        − 75,         (2.61)
                                     i=1
                                                         ci

where the Gaussian model constants are defined by,

                        a1 = 309.2,      b1 = 0.04824,           c1 = 0.03661,         (2.62)
                        a2 = 1482,       b2 = −0.0236,           c2 = 0.05637,         (2.63)
                        a3 = 85.97,      b3 = 0.08573,           c3 = 0.02410.         (2.64)

The original and smooth F (λθ ) functions are plotted in Figure 2.6. This modification was
found to have a negligible impact on the predicted transition locations.

2.3    Compressibility Corrections

Two compressibility corrections are applied, ψ and ψscf , in order to extend the LM2009
(Equations 2.8–2.9) [92] and LM2015 (Equations 2.16–2.21) [93] empirical correlations in the
SA-sLM2015 transition model to include the stabilizing effects of Mach number on TS and
stationary crossflow instabilities, respectively. The empirical correlations in these models
were originally developed through calibration with experimental results. However, there
is limited availability of experimental results at transonic flight conditions, especially at
higher Reynolds numbers where the stabilizing effects of compressibility are more significant.
This can be attributed to the relative significance of surface imperfections on thin wind-
tunnel model boundary layers at high Reynolds numbers and high wind-tunnel freestream
turbulence intensities at high Mach numbers, as an increase in the wind-tunnel freestream

<!-- ===== PDF page 48 ===== -->

28                                                         CHAPTER 2. TRANSITION MODEL EQUATIONS

Mach number is often associated with an increase in the freestream turbulence intensity [156].
Therefore, the approach taken in the current work is to build off of existing and well-validated
correlations in the literature.
   The TS instability compressibility correction, ψ, is developed based on results from the
compressible extension to the AHD criterion [10, 136] and is calibrated with the smooth
transition model using the test cases presented in Section 4.2. In order to facilitate the
integration of the stationary crossflow instability compressibility correction developed by
Malik et al. [108], ψscf , with the SA-sLM2015 transition model, a new crossflow instability
source term, Fonset,scf , is introduced in the intermittency transport equation to replace the
original Dscf source term (Equation 2.15) in the momentum-thickness Reynolds number
transport equation. The smooth transition model with the compressibility corrections is
designated SA-sLM2015cc.
   It is important to note that these compressibility corrections do not attempt to model
the complex flow physics of compressible boundary-layer transition. Consistent with the
underlying transition model, they provide a framework for including the first-order effects
of flow compressibility on the location of boundary-layer transition. The corrections ex-
tend the domain of applicability for these commonly used empirical correlations and do not
affect the predictive capabilities of the models in the incompressible flow regime, as both
compressibility corrections, ψ and ψscf , approach unity for low Mach numbers.

2.3.1   Tollmien-Schlichting Instabilities

In an effort to compare the behaviour of the LM2009 empirical correlation [92] to stability
analysis, Perraud et al. [136] developed a simplified model based on the database of stability
analysis results that were used to develop the compressible extension to the AHD criterion [9,
10, 137]. The model was developed for Mach numbers up to 1.1 and low turbulence intensities
(T u < 1%). Their model, which approximates the sensitivity of the linear stability equations
to varying turbulence intensities, T u, pressure gradients, λθ,e , and flow compressibility, Me ,
is defined as,

                                          ρθ2 dUe
     Reθt = f (T u, λθ , Me ),   λθ,e =           Re∞ ,   T u = T u∞ ,                    (2.65)
                                           µ ds
     Reθt = −(177Me2 − 22Me + 210) ln((7Me + 4.8)T u/100) exp((5Me + 27)λθ,e ),           (2.66)

where Re∞ is again introduced as part of the non-dimensionalization procedure.
   A comparison of the sensitivities of the LM2009 empirical correlation and the stability-
based model to pressure gradient and Mach number is illustrated in Figure 2.7. The LM2009
empirical correlation [92] appears to under-predict the effect of pressure gradient relative to
the stability-based model. The LM2009 empirical correlation [92] is based on the empirical
correlations by Abu-Ghannam et al. [1], which were developed based on experimental data
with a focus on higher turbulence intensity conditions (T u > 0.3%). It is also important

<!-- ===== PDF page 49 ===== -->

2.3. COMPRESSIBILITY CORRECTIONS                                                                         29

  4000                                                  2200

  1500                                                  1600

     0                                                  1200
     -0.1       -0.05        0          0.05      0.1          0   0.2       0.4    0.6     0.8      1

                (a) Pressure gradient                                    (b) Mach number

Figure 2.7: Sensitivity of the LM2009 empirical correlation [92], the stability-based model, and stability
analysis [136] to pressure gradient and Mach number. A higher momentum-thickness Reynolds number
delays boundary-layer transition.

to note that the AHD criterion is calculated using an averaged parameter, λ̄θ,e , which is
computed by integrating the local pressure gradient parameter, λθ,e , from the critical point
to the location being evaluated, while the LM2009 empirical correlation is calculated using
the local pressure gradient, λθ . Together these factors can help to explain the differences
in Figure 2.7a. However, as Figure 2.7b demonstrates, the LM2009 empirical correlation
under-predicts the momentum-thickness Reynolds number at higher Mach numbers by not
including the stabilizing effects of flow compressibility.
   To isolate the effects of flow compressibility predicted by the stability-based model, the
momentum-thickness Reynolds number is investigated over a range of boundary-layer edge
Mach numbers, Me , normalized by the value at a boundary-layer edge Mach number of zero.
The effects of flow compressibility at various turbulence intensities and pressure gradients
predicted by the normalized stability-based model are compared in Figure 2.8. The trends
predicted by the stability-based model are similar over the range of turbulence intensities,
Figure 2.8a, and pressure gradients, Figure 2.8b, investigated. An increasing Mach number
produces a more significant stabilizing effect in more favourable disturbance environments,
which are represented by low turbulence intensities, and in regions of more favourable pres-
sure gradient (here corresponding to a more positive λθ,e ). This behaviour is consistent with
the analysis by Ströer et al. [181], who investigated the compressible AHD criterion on a flat
plate with varying pressure gradients. As the Reynolds number increases, transition often
occurs in regions of more favourable pressure gradient. The stabilizing effect of flow com-
pressibility is therefore most relevant at these higher Reynolds numbers, which are typical
of commercial transport aircraft.
   In order to include the effects of flow compressibility on the LM2009 empirical correlation

<!-- ===== PDF page 50 ===== -->

30                                                                      CHAPTER 2. TRANSITION MODEL EQUATIONS

      2                                                         2

     1.8                                                      1.8

     1.6                                                      1.6

     1.4                                                      1.4

     1.2                                                      1.2

      1                                                         1

           0       0.2     0.4       0.6        0.8   1             0       0.2        0.4      0.6       0.8    1

                     (a) Turbulence intensity                                     (b) Pressure gradient

Figure 2.8: Sensitivity of the normalized stability-based model (blue), stability analysis [136], and the initial
compressibility correction, ψinit (grey) (Equations 2.67–2.69), to Mach number with varying turbulence
intensity and pressure gradient.

for TS instabilities, an initial compressibility correction, ψinit , was developed to reproduce the
trends predicted by the normalized stability-based model. The TS compressibility correction
includes the effects of flow compressibility with respect to the Mach number at the edge of
the boundary layer, Me , including the effects of freestream turbulence intensity, T u∞ , and
a local approximation to the boundary-layer edge pressure gradient, λθ,e . The correction is
given by,

                                            Reθc,comp = ψinit Reθc ,                                            (2.67)
                                                                         p     
               ψinit = (a1 Me + a2 Me + a3 ) exp[(b1 λθ,e )Me ] exp c1 + c2 T u∞ Me ,                           (2.68)

               a1 = 0.44, a2 = −0.38, a3 = 1.00;          b1 = 5.00;       c1 = 0.41, c2 = −0.27,               (2.69)

and is compared with the normalized stability-based model in Figure 2.8.
   The correction is integrated in the SA-sLM2015 smooth local correlation-based transition
model by scaling the incompressible critical momentum-thickness Reynolds number, Reθc ,
which is calculated using the smooth empirical correlation (Equation 2.54). This approach is
similar to scaling the local value of the transported momentum-thickness Reynolds number
 ˜ θt , as Reθc in the smooth transition model is a near-linear function of Re
Re                                                                            ˜ θt . Applying the
compressibility correction to the local momentum-thickness Reynolds number, Reθt , which is
then convected and diffused into the boundary layer by the Re  ˜ θt transport equation, was also
investigated. However, this resulted in the correction being overly stabilizing, particularly
for the transonic CAST10-2 airfoil test case (presented in Section 4.2) at near-zero angles
of attack. The rapid acceleration at the leading edge of an airfoil produced large upstream
                                                        ˜ θt . Due to the delay in the influence
values of the compressibility correction and therefore Re

<!-- ===== PDF page 51 ===== -->

2.3. COMPRESSIBILITY CORRECTIONS                                                                 31

of Reθt in the boundary layer, where the Pθt source term (Equation 2.2) is not active, this
resulted in the inability of the correction to sufficiently reduce the stabilizing effects of
compressibility and therefore Re˜ θt in regions of adverse pressure gradient, and consistently
produced predicted transition locations downstream of experiment. A consequence of the
current approach is that the Flength empirical correlation (Equation 2.52) is not scaled by
the compressibility correction. However, the stability-based model and the AHD criterion
do not predict the effects of compressibility on the transition region length. For example,
the AHD-based transition models developed by Pascal et al. [133] and Ströer et al. [180, 181]
rely on user-specified inputs for the transition region length. Further work is required to
develop a correlation for transition length in compressible flows.
    Initial simulations using the transonic transition test cases presented in Section 4.2 re-
vealed that the smooth transition model with the compressibility correction presented in
Equations 2.67–2.69 often over-predicted the extent of the laminar boundary layer, partic-
ularly for the transonic VA-2 airfoil test case at lower angles of attack, and appeared to
over-estimate the impact of favourable and adverse pressure gradients on the stabilizing
effect of compressibility. This behaviour may be due to the aforementioned inherent dif-
ferences between the stability-based approaches, including the AHD criterion, and the fully
local LM2009 empirical correlation. For example, as the local pressure gradient formulation
does not take into account pressure gradient history effects, the compressibility correction
in a fully local framework is more sensitive to local pressure gradients.
   Therefore, the compressibility correction was calibrated with the smooth transition model.
Rather than calibrating for a single test case, the correction was calibrated using the cases
in Section 4.2 in order to improve the predictive capability of the model across the full range
of test cases and flow conditions investigated (0.71 ≤ M ≤ 0.856, 2 × 106 ≤ Re ≤ 15 × 106 ).
These test cases involve flow conditions where the stabilizing effect of flow compressibility is
expected to be both small (CAST10-2 airfoil) and significant (VA-2 airfoil and CRM-NLF
geometry). During this calibration, the compressibility correction was made less sensitive to
the effects of Mach number and pressure gradient. Furthermore, due to the limited availabil-
ity of experimental data at transonic flow conditions, the sensitivity of the compressibility
correction to turbulence intensity is unmodified from the normalized stability-based model.
The calibrated compressibility correction, ψ, is given by,

                                       Reθc,comp = ψReθc ,                                  (2.70)
                                                                   p     
             ψ = (a1 Me + a2 Me + a3 ) exp[(b1 λθ,e )Me ] exp c1 + c2 T u∞ Me ,             (2.71)

       a1 = 0.34, a2 = −0.38, a3 = 1.00;       b1 = 3.00;    c1 = 0.41, c2 = −0.27,         (2.72)

and is compared with the initial correction, ψinit (Equations 2.67–2.69), in Figure 2.9. The
a1 and b1 constants were reduced relative to the initial correction, ψinit , in order to reduce the
sensitivity to Mach number and pressure gradient, respectively. The results in Section 4.2

<!-- ===== PDF page 52 ===== -->

32                                                                CHAPTER 2. TRANSITION MODEL EQUATIONS

      2                                                   2

     1.8                                                1.8

     1.6                                                1.6

     1.4                                                1.4

     1.2                                                1.2

      1                                                   1

           0   0.2     0.4       0.6        0.8   1           0       0.2        0.4      0.6       0.8   1

                 (a) Turbulence intensity                                   (b) Pressure gradient

Figure 2.9: Sensitivity of the initial, ψinit (grey) (Equations 2.67–2.69), and calibrated, ψ (red) (Equa-
tions 2.70–2.72), compressibility corrections to Mach number with varying turbulence intensity and pressure
gradient.

demonstrate that the calibrated compressibility correction coupled with the smooth transi-
tion model produces results that match the transonic experiments well; however, additional
experimental data is required to perform a more detailed validation study of the resulting
model at the high Reynolds number transonic flow conditions where the stabilizing effect of
flow compressibility is expected to be significant. It is also important to note that recent
work [180] suggests that the Langtry-Menter transition model over- and under-predicts the
extent of the laminar boundary layer at low and high Reynolds numbers, respectively, due
to the reduced sensitivity to pressure gradient illustrated in Figure 2.7a, which may impact
the results.
   The compressibility corrections require the Mach number at the edge of the boundary
layer, Me , which can be approximated locally using isentropic relations and the compressible
Bernoulli equation as follows [29]:
                                           s
                                                        1− κ1
                                Ue            κP∞ P
                          Me = , ae =                           ,                      (2.73)
                                ae             ρ ∞ P∞
                                s
                                                             1− κ1 
                                    2 +
                                          2κ   P ∞         P
                          Ue = U∞                    1−                 ,              (2.74)
                                         κ − 1 ρ∞         P∞

where κ is the heat capacity ratio. As previously mentioned, the pressure gradient parameter,
λθ (Equations 2.10–2.14), used in the LM2009 empirical correlation is based on the local
streamwise velocity gradient formed using Cartesian velocity gradients, which is not valid in
the boundary layer. To avoid problems with this formulation, the LM2009 transition model
introduces the Fθt function (Equation 2.3) in the Pθt source term (Equation 2.2) to disable the
TS empirical correlation inside the boundary layer. However, the compressibility corrections

<!-- ===== PDF page 53 ===== -->

2.3. COMPRESSIBILITY CORRECTIONS                                                               33

are most active in the middle of the boundary layer where the transition onset process
begins. In order to address this, the pressure gradient parameter used for the compressibility
correction in the current work is calculated using a local approximation to the boundary-
layer edge pressure gradient, λθ,e , adopted from Grabe et al. [60, 59], which is approximated
using Cartesian pressure gradients and isentropic relations as follows:

                                       ρθ2 dUe
                                λθ,e =         Re∞ ,                                      (2.75)
                                        µ ds
                                                           
                        dUe     u dUe        v dUe          w dUe
                            =           +             +             ,                     (2.76)
                         ds     U dx         U dy           U dz
                                                    − κ1
                                dUe        1     P         dP
                                     =−                       ,                           (2.77)
                                dx       ρ∞ Ue P∞          dx
                                                    − κ1
                                dUe        1     P         dP
                                     =−                       ,                           (2.78)
                                dy       ρ∞ Ue P∞          dy
                                                    − κ1
                                dUe        1     P         dP
                                     =−                       .                           (2.79)
                                 dz      ρ∞ Ue P∞          dz

Similar to the original pressure gradient parameter, λθ (Equations 2.10–2.14), λθ,e is not
Galilean invariant.
    It is important to re-emphasize that the boundary-layer edge pressure gradient parameter,
λθ,e , defined above is not used to calculate the Reθt empirical correlation (Equations 2.8–2.9).
Although the two pressure gradient formulations produce similar results for favourable and
adverse pressure gradients, the original, velocity-gradient-based pressure gradient parameter,
λθ , is used to be consistent with the Langtry-Menter model [92], as it provides better agree-
ment with experimental data for zero-pressure-gradient flat plate cases. This is likely due
to the calibration of the model. For numerical robustness, λθ,e is limited using the bounds
introduced by Langtry and Menter [92] for λθ (−0.1 ≤ λθ,e ≤ 0.1).

2.3.2   Stationary Crossflow Instabilities

In the LM2015 [93] and sLM2015 [139] transition models, the stationary crossflow
momentum-thickness Reynolds number correlation, Reθt,scf (Equation 2.16), is included using
                                  ˜ θt transport equation, Dscf (Equation 2.15). The natural
an additional source term in the Re
and bypass transition correlations influence Re˜ θt in the freestream, which is then convected
and diffused into the boundary layer where the stationary crossflow source term is active.
The combined effect of these correlations is realized through the Fonset (Equation 2.32) func-
tion in the intermittency source terms through the critical momentum-thickness Reynolds
number, Reθc .
   In order to include the stabilizing effects of compressibility on each set of correlations,
and subsequent Reynolds numbers, Reθt and Reθt,scf , separately, Fonset is separated into two
functions, one for two-dimensional mechanisms and one for stationary crossflow instabilities.

<!-- ===== PDF page 54 ===== -->

34                                                        CHAPTER 2. TRANSITION MODEL EQUATIONS

A similar approach to that developed by Grabe et al. [59] and Carnes and Coder [18] is
adopted. A new onset function, Fonset,scf , is introduced for the intermittency transport equa-
tion replacing the Dscf source term in the Re  ˜ θt transport equation. In the original LM2015
model [93], the effects of Reθt,scf are realized through the Fonset function by reducing Re ˜ θt
and therefore Reθc in the boundary layer. To be consistent with Langtry et al. [93], the
Fonset,scf function was developed using the same approach as the original Fonset function,
and therefore the smooth approximation (Equation 2.49). The stationary crossflow critical
momentum-thickness Reynolds number, Reθc,scf , correlation in the new Fonset,scf function was
developed using the NLF2-0415 Infinite Swept Wing [38, 148] transition test case, leading
to the following:

                                        tanh(6(Fonset,scf,1 − 1.35)) + 1
                           Fonset,scf =                                  ,                 (2.80)
                           s                          2
                                           2  2
                                   ReS
            Fonset,scf,1 =                    + RT , Reθc,scf = 0.6231Reθt,scf ,           (2.81)
                              2.60Reθc,scf

where ReS and RT represent the strain-rate magnitude Reynolds number and eddy viscosity
ratio, respectively. The transition process is triggered when Fonset,scf,1 exceeds unity, and is
sustained as the eddy viscosity ratio, RT , increases. The Reθc,scf correlation consists of a
value of 0.6231 applied to Reθt,scf in order to initiate the transition process upstream of the
predicted transition location. This value closely resembles a linear fit of the smooth Reθc
function (Equation 2.54) and is similar to the value of 0.62 used by Medida [115].
    The Fonset,scf function is incorporated into the intermittency source term equations as
follows:
                                                             √        
                                                             M∞ M∞ Re∞ √
     Pγ = ca1 Flength [φ300 (Fonset , Fonset,scf )] φ−300 Ω,                  γ(1 − ce1 γ), (2.82)
                                                        √       
                                                      M∞ M∞ Re∞
                       Eγ = ca2 Fturb φ−300 Ω,                      γ(ce2 γ − 1).           (2.83)

The Fturb function deactivates the destruction term when Fonset becomes active or until a
threshold eddy viscosity ratio, RT , is reached and is modified to include the new crossflow
onset function. The new Fturb function is defined as,

                        Fturb = (1 − [φ300 (Fonset , Fonset,scf )] exp(−RT ).              (2.84)

   In order to validate the new Fonset,scf crossflow source term, the SA-sLM2015 transition
model is applied to the TU Braunschweig Sickle Wing [138] subsonic transition test case both
with the new source term and with the original Dscf crossflow source term, with the results
presented in Appendix A. The transition fronts on both the upper and lower surfaces of the
wing are overlaid with the transition fronts from experiment [138]. The results demonstrate
that the SA-sLM2015 transition model with the new Fonset,scf crossflow source term does a

<!-- ===== PDF page 55 ===== -->

2.3. COMPRESSIBILITY CORRECTIONS                                                                  35

reasonable job of reproducing the transition fronts from the experiment, and matches well
the transition front produced by the transition model with the original Dscf crossflow source
term. Validation studies of the SA-sLM2015 transition model with the original Dscf crossflow
source term with the NASA NLF2-0415 [38, 148] and TU Braunschweig Sickle Wing [138]
transition test cases are presented in Section 4.1 [139].
   Simulations of transonic swept wings revealed that both crossflow source-term approaches,
Dscf and Fonset,scf , produce a transition front that is dependent on the initialization of the
flow field. A solution with the crossflow correlation active initialized with far-field conditions
converges to a transition front upstream of the converged transition front produced with the
crossflow correlation active initialized with a solution with the crossflow correlation inactive.
It is important to note that this behaviour does not appear on subsonic transition test
cases with crossflow, such as the NASA NLF2-0415 [38, 148] and TU Braunschweig Sickle
Wing [138]. To prevent this behaviour, for transonic cases the crossflow source term is
activated after the total residual drops several orders of magnitude without the crossflow
source term, Fonset,scf , active. A relative residual drop tolerance of five orders of magnitude for
the total residual was found to be sufficient to prevent these inaccurate, upstream transition
fronts and is the tolerance used for the results presented in the current work.
   Crossflow instabilities, and their associated inflection points, are less influenced by flow
compressibility than viscous instabilities, such as TS waves [39, 6]. Malik et al. developed
a compressibility correction for a crossflow Reynolds number [108] which is used in the cur-
rent work by scaling the incompressible stationary crossflow momentum-thickness Reynolds
number, Reθt,scf , in the Fonset,scf function as follows:

                                    Reθt,scf,comp = ψscf Reθt,scf ,                          (2.85)
                                                   κ−1 2
                                     ψscf = 1 +           Me .                               (2.86)
The two compressibility corrections, ψ and ψscf , for TS and stationary crossflow instabilities,
respectively, are presented in Figure 2.10. As expected, the compressibility correction for TS
instabilities produces a stronger stabilizing effect than the crossflow instability compressibil-
ity correction, especially in a favourable flow environment (T u∞ = 0.05% and λθ,e = 0.04).

<!-- ===== PDF page 56 ===== -->

                             1.8

                             1.6

                             1.4

                             1.2

                                   0   0.2      0.4        0.6   0.8      1

Figure 2.10: Effects of Mach number on the TS and stationary crossflow instability compressibility correc-
tions, ψ and ψscf , respectively.

<!-- ===== PDF page 57 ===== -->

Chapter 3

Optimization Framework

The high-fidelity discrete-adjoint gradient-based aerodynamic shape optimization framework
used in the current work, Jetstream, consists of three primary components: a geometry pa-
rameterization and control scheme [70, 55, 95], a parallel structured multi-block Newton-
Krylov-Schur RANS-based flow solver [69, 131], and a discrete-adjoint gradient-based op-
timization algorithm [70, 128]. Jetstream was recently cross-validated with an industry
RANS-based discrete-adjoint optimization algorithm, with the two algorithms producing
very similar geometries and performance improvements [152]. The flow solver, Diablo, was
extensively validated with fully turbulent flow as part of the Fifth AIAA Drag Prediction
Workshop [130]. This chapter will provide an overview of these components with an empha-
sis on the flow solver and gradient evaluation, which were modified for the introduction of
the transition model equations [139, 140, 141].

3.1    Parallel Implicit Newton-Krylov-Schur Flow Solver

An efficient, robust, and accurate flow solver is crucial to the performance of an optimization
algorithm, as flow analysis is conducted many times over the course of the optimization
process, with a variety of geometries expected. The flow solver used in the current work
is a three-dimensional structured multi-block finite-difference solver developed by Hicken
and Zingg [69] for the solution of the Euler equations and extended to the RANS equations
coupled to the one-equation SA turbulence model by Osusky and Zingg [131]. The governing
equations are spatially discretized using summation-by-parts operators, with simultaneous
approximation terms applied to enforce boundary conditions and inter-block coupling. To
decrease the computational time required to complete a flow solution, the computational
domain is decomposed into multiple blocks, resulting in multi-block structured grids, which
allow for efficient parallel computations.
   Spatial discretization applied to the governing equations produces the following system
of nonlinear ordinary differential equations:

<!-- ===== PDF page 58 ===== -->

38                                                           CHAPTER 3. OPTIMIZATION FRAMEWORK

                                      dQ
                                          + R(Q) = 0,                                    (3.1)
                                       dt
where Q and R represent the array of conserved variables and the spatially discretized total
residual consisting of the mean-flow, turbulence, and transition model equations. Consider-
ing a steady-state solution produces the following nonlinear system of algebraic equations:

                                         R(Q) = 0.                                       (3.2)

Equation 3.2 can be solved using a range of explicit and implicit methods. RANS-based
simulations require small off-wall spacings to satisfy the grid resolution requirements of the
turbulence model, which produce high-aspect-ratio cells in the boundary layer. These cells
can reduce the performance of explicit methods due to restrictive time steps. An implicit
strategy allows for large time steps, and when combined with a local time linearization, can
provide stable and efficient convergence. Newton’s method is a widely adopted nonlinear
solver, which can provide rapid convergence if an acceptable initial guess is supplied.
    In the current work, a fully coupled, implicit Newton-Krylov-Schur solution algorithm
making use of a pseudo-transient continuation strategy is applied to Equation 3.2 to drive
the solution from an initial guess to a converged steady-state solution. The solution strat-
egy consists of two phases, an approximate-Newton phase followed by an inexact-Newton
phase. Globalization of the inexact-Newton method is achieved using the approximate-
Newton phase, where an implicit Euler time marching strategy with local time linearization
is applied using an approximate, analytically derived Jacobian. The application of implicit
Euler time marching produces the following system of linear equations at each outer iteration,
(n):

                             (T(n) + A(n) )∆Q (n) = −R(Q (n) ),                          (3.3)

where T(n) is a diagonal matrix containing the inverse local time steps and A(n) is the flow
Jacobian given by,

                                       (n)   ∂Ri
                                     Aij =       (Q (n) ).                               (3.4)
                                             ∂Qj

The large system of linear equations (Equation 3.3) generated at each iteration in both
phases is solved inexactly using the GMRES linear solver [160], preconditioned using an
approximate-Schur parallel preconditioner [161]. Details on the formulation of the approx-
imate Jacobian and local time step can be found in [129]. After the residual drops several
orders of magnitude, the inexact-Newton method using the full Jacobian is used to converge
the system to a residual norm of machine zero, where the inexact-Newton method is recov-
                                                        (n)
ered by aggressively ramping the reference time step, ∆tref . The full Jacobian can be formed
explicitly or a Jacobian-free approach can be used.

<!-- ===== PDF page 59 ===== -->

3.1. PARALLEL IMPLICIT NEWTON-KRYLOV-SCHUR FLOW SOLVER                                         39

   The switch from the approximate- to the inexact-Newton phase is controlled using a
residual drop tolerance based on the L2 -norm of the total residual, including the row (equa-
tion) scaling measure presented in Section 3.1.2, normalized by the total residual at the first
iteration,

                                       (n)   ||R(Q (n) )||2
                                     Rd =                   .                               (3.5)
                                             ||R(Q (0) )||2

Residual drop tolerances of 10−4 and 10−5 are used for the two- and three-dimensional
simulations, respectively, presented in the current work.

3.1.1   Considerations for the Transition Model Equations

Local transport-equation-based transition models contain large and highly nonlinear source
terms [92, 116, 29, 118]. These source terms introduce numerical stiffness through time
scales that can vary significantly from those introduced by the convection and diffusion
terms. Unless suitable steps are taken to modify the basic algorithm, the implementation of
transition model equations in a strong implicit solver can lead to numerical instability.
    The question of how to best solve transport equations with strong source terms was em-
phasized with the introduction of two-equation turbulence models, such as those based on the
popular k- and k-ω models. Moryossef and Levy [122] provide a review of solution strategies
for the stable treatment of source terms in these models. The methods focus on positivity
preservation of the independent working variables through the implicit treatment of sinks
and the explicit treatment of sources. However, solution update clipping is still introduced
in order to prevent excessive solution growth. Work by Lian et al. [100] demonstrated the
effect of explicit and implicit treatments of source terms through stability-analysis of a linear
convection equation with a source term. The results demonstrate that the stability of the
system can be represented using the non-dimensional source-term time step, where a large
positive source-term time step can lead to strong solution growth and subsequent flow solver
divergence. While implicit methods provide more stable updates than explicit methods for
the treatment of sinks, i.e. negative sources, the implicit treatment of sources can produce
unstable solution updates for large source-term time steps [100]. Specifically, implicit Euler
time-marching with a source-term time step larger than unity can produce solution updates
with an unphysical sign change. This behaviour can help to explain the common problem
of positivity preservation of turbulence and transition model working variables [120, 3]. In
this regard, Spalart and Allmaras developed a linearization strategy for the SA turbulence
model that guarantees positivity of the modified eddy viscosity when using a first-order up-
wind discretization strategy [174]. Their strategy involves modifying the turbulence model
Jacobian based on the sign of the linearization of the source terms.
    Strategies for the stable treatment of turbulence model source terms have recently been
applied to transition model equations as well. Similar to the work on positivity preservation

<!-- ===== PDF page 60 ===== -->

40                                                         CHAPTER 3. OPTIMIZATION FRAMEWORK

of two-equation turbulence models by Ilinca and Pelletier [74], Coder [28] reformulated the
intermittency equation in the AFT transition model in order to solve for the logarithmic vari-
able, γ̃ = log(γ). Although this strategy ensures a positive solution, it introduces difficulties
with the treatment of homogenous boundary conditions at the wall, i.e. γ̃ = 0, and can
require additional solution update clipping to prevent excessive variable growth [122]. Mosa-
hebi and Laurendeau [123, 124] investigated different coupling and time-marching strategies
for the Langtry-Menter LM2009 [92] transition model and the SST turbulence model [117].
Specifically, four linearization strategies were investigated: a decoupled approach with an
independent linearization for each of the turbulence and transition model equations, an
under-relaxed decoupled approach with a blending function applied to the effective intermit-
tency, a loosely coupled approach with the transition and turbulence models solved sepa-
rately resulting in two 2 × 2 Jacobians, and a loosely coupled approach with the turbulence
and transition models solved coupled resulting in one 4 × 4 Jacobian. For each approach,
the mean-flow equations were solved decoupled from the turbulence and transition model
equations. The results demonstrate that of the four strategies considered, the under-relaxed
decoupled approach provided the most robust convergence for the subsonic airfoil test cases
investigated [123, 124].
    While Newton-Krylov solvers provide the potential for greatly accelerated nonlinear con-
vergence [83, 14, 24, 131], the rate of convergence of Newton’s method is closely linked to the
accuracy of the Jacobian. Neglecting terms or otherwise manipulating the Jacobian can pre-
vent quadratic convergence [24]. However, Blanco and Zingg [14] demonstrated the potential
benefits of using a loosely coupled strategy for solving the SA turbulence model equations
in an unstructured Newton-Krylov flow solver. Specifically, they demonstrated that the
loosely coupled approach produced similar or better nonlinear convergence relative to fully
coupled Newton-Krylov solvers, with the loosely coupled approach also benefiting from a
30% reduction in memory requirements. More recently, Yildrim et al. [197] investigated
coupling strategies for the SA turbulence model and mean-flow equations in an approximate
Newton-Krylov flow solver. They demonstrated that solving the SA turbulence model loosely
coupled to the mean-flow equations in the globalization phase could improve robustness and
efficiency, with the turbulence model solved fully coupled in the inexact-Newton phase to
facilitate nonlinear convergence.
   Similar approaches have been investigated for the SA-sLM2015 transition model with the
results presented in Appendix B. The results demonstrate that a fully coupled approach
improves nonlinear convergence as the complexity of the flow increases. Specifically, for
more complicated problems, the results demonstrate that loosely coupled and decoupled
solution strategies can result in the solution stalling in the approximate-Newton phase, and
that the mean-flow, turbulence, and transition model equations require a tight coupling in
order to achieve deep nonlinear convergence. However, a strategy is required to stabilize the
implicit treatment of the large transition model source terms, where the source terms can

<!-- ===== PDF page 61 ===== -->

3.1. PARALLEL IMPLICIT NEWTON-KRYLOV-SCHUR FLOW SOLVER                                       41

remain active throughout the inexact-Newton phase. Rather than modifying the Jacobian
matrix, Lian et al. [100] demonstrated that restricting the product of the largest positive
eigenvalue of the source-term Jacobian matrix and the local time step can be an effective
way of stabilizing large sources in a fully coupled, implicit algorithm. A source-term time
step restriction based on the ideas developed by Lian et al. [100] has been implemented in
the current work with the details presented in Section 3.1.4, and is investigated in more
detail in Appendix B.

   Additional modifications to the Newton-Krylov algorithm were implemented to improve
robustness. These include equation and variable scaling measures and a solution update
damping algorithm, including a physics-based restriction and an unsteady residual back-
tracking line search [120].

3.1.2   Equation and Variable Scaling

The turbulence and transition model variables can vary significantly in magnitude from the
mean-flow variables and do not contain inherent geometric scaling [24]. In addition, large
off-diagonal entries are introduced through the linearization of the turbulence and transition
model equations with respect to the mean-flow variables. If not accounted for, these factors
can result in improper scaling of the linear system and poor performance of the linear solver.
To address this, row and column scaling measures, including residual auto scaling, are used
based on the work by Chisholm and Zingg [24] and Osusky and Zingg [129, 131]. Instead
of solving the linear system of equations (Equation 3.3) directly, a scaled version of the
equations is solved, given by:

                      Sa Sr (T(n) + A(n) )Sc Sc−1 ∆Q (n) = −Sa Sr R(Q (n) ),              (3.6)

where Sa represents a residual auto-scaling matrix [129, 131] and

                         2                                                          
                         Ji3
                        
                              ...                                                   
                                                                                     
                                                                                    
                                      2                                             
                                    Ji3                                             
                 Sr,i =                           −1
                                                                                     ,   (3.7)
                                             −1
                                           ν̃max Ji 3
                                                                                    
                                                                                    
                                                              −1                    
                                                         −1
                        
                                                       γmax Ji 3                    
                                                                                     
                                                                                −1
                                                                    ˜ −1
                                                                    Re θt,max Ji

<!-- ===== PDF page 62 ===== -->

42                                                             CHAPTER 3. OPTIMIZATION FRAMEWORK
                                                                         
                                  .
                                     ..
                                                                          
                                                                         
                                                                         
                                                                         
                          Sc,i =                                         ,                 (3.8)
                                                                         
                                 
                                         ν̃max                           
                                                                          
                                                γmax
                                                                         
                                                                         
                                                               ˜ θt,max
                                                               Re

are the ith blocks of the diagonal row and column scaling matrices, respectively. Through a
thorough algorithm optimization procedure, Osusky and Zingg [131] determined an optimal
normalization value for ν̃max of 103 , which accounts for the maximum turbulence value that
is likely to be encountered in the flow solve. Following a similar approach, γmax and Re   ˜ θt,max
values of 10 and 104 , respectively, are used. The partially scaled residual, Sr R(Q (n) ), is used
when monitoring convergence (Equation 3.5). Full details on the auto, row, and column
scaling measures can be found in Osusky and Zingg [129, 131].

3.1.3   Solution Update Damping

An upstream laminar boundary layer can increase shock strength, while a laminar boundary
layer is more sensitive to the effects of pressure gradient and is more prone to boundary-
layer separation. If not properly accounted for, these flow features can cause the solution
algorithm to produce unphysical solution updates leading to flow solver divergence. To
maintain stable solution updates while improving the convergence rate of the flow solver, a
solution update damping algorithm was developed based on the transition and mean-flow
variable updates and an unsteady residual backtracking line-search method [120]. As an
alternative to solution update damping, solution-limited time-stepping methods, such as
those developed by Chisholm and Zingg [24], and Lian et al. [99, 100], were investigated.
However, these methods can slow convergence for cases with stable updates and lead to
restrictive time steps for variables approaching zero, which is common for the intermittency
and eddy viscosity-like variables.

Physics-Based Restriction

A physics-based restriction was implemented based on the work of Modisette [120]. To
prevent unphysical solution updates, negative solution updates are limited based on the size
of the density and energy solution updates relative to the current values. The physics-based
restriction is defined as,

                                          (n)          (n)
                                       |∆Qρ,e | ≤ 0.90Qρ,e ,                                 (3.9)

where the solution update limit is set to 90% of the current value, which was found to
provide a good balance between speed and robustness. To enforce this limit, negative solution

<!-- ===== PDF page 63 ===== -->

3.1. PARALLEL IMPLICIT NEWTON-KRYLOV-SCHUR FLOW SOLVER                                                    43

Algorithm 2 Physics-based restriction
 1: δlocal ← 1                                                                   . Initialize damping factor
 2: for j = 1, nnodes do                                                            . Loop over local nodes
 3:     for l = [1, 5] do                                                              . Loop over variables
                   (n)
 4:        if ∆Qj,l < 0 then                                                       . Limit negative updates
                              (n)
                                             −1
                             ∆Qj,l      1
 5:           δlocal = max      (n) , δlocal                         . Obtain local minimum damping factor
                              0.90Qj,l

 6: δphys ← MPIAllreduce(δlocal , min)                         . Determine global minimum damping factor
 7: if δphys < 0.01 then                                                            . Poor solution update
             (n)
 8:    if ∆tref > ∆tref,min then                                         . Reject current update and retry
             (n)              (n)
 9:       ∆tref = max (0.5∆tref , ∆tref,min )                                 . Reduce reference time step
10:       δphys ← 0                                                                        . Reject update
11:       return                                                                           . Retry solution
12:    else
13:       δphys = 0.01                                                          . Clip and continue solution

updates, ∆Q(n) , are damped using a damping factor, δphys , which is defined as,

                                                      0.90Q(n)
                                                           ρ,e
                                            δphys =              .                                   (3.10)
                                                       ∆Q(n)
                                                         ρ,e

The full algorithm for the physics-based solution update damping is presented in Algorithm 2.
For unstable solution updates that require a damping factor of 1% or less to maintain a
physical update, we attempt to reject the solution update and try again with the reference
               (n)
time step, ∆tref , cut in half. This makes the linear system more diagonally dominant due
to the inverse time step in Equation 3.3, which improves convergence of the linear system
and reduces the magnitude of the solution update. However, if the time step is at the user-
specified lower bound, ∆tref,min , we clip the update factor and damp the local solution update
to maintain positivity of density and energy and continue with the solution update.

Transition Model Variable Bounds

The transition model equations can become unstable if the intermittency and transition onset
momentum-thickness Reynolds number are not bounded. In order to ensure stable transition
model updates that respect the transition model variable bounds, a block-based solution
update damping procedure was implemented, as described in Algorithm 3. The solution
updates are determined including the physics-based damping factor presented in Algorithm 2.
Clipping the transition model variables at these bounds was also investigated; however,
placing a hard limit on variable bounds can potentially lead to stalling of the nonlinear
solver. It is important to note that with the implementation of the source-term time step
restriction, which prevents unstable solution updates, in combination with a more dissipative
first-order upwind discretization strategy for the turbulence and transition model equations,
the upper and lower bounds specified in Algorithm 3 are rarely reached. Furthermore,
steady-state solutions for intermittency respect the bounds of 0.02 and 1 implicitly defined

<!-- ===== PDF page 64 ===== -->

44                                                               CHAPTER 3. OPTIMIZATION FRAMEWORK

Algorithm 3 Solution update with transition model update damping
 1: θfac ← 0.99                                                                    . Damping factor
 2: for j = 1, nnodes do                                                    . Loop over local nodes
 3:     for l = 1, nvar do                                                    . Loop over variables
            (n+1)     (n)           (n)
 4:        Qj,l    = Qj,l + δphys ∆Qj,l                                           . Update solution
 5:        m=0                                                            . Initialize count variable
                    (n+1)          (n+1)
 6:        while (Qj,7 > 2 or Qj,7 < 10−10 ) and (θfac )m > 0.01 do           . Limit intermittency
 7:            m=m+1                                           . Count number of damping iterations
                (n+1)     (n)                 (n)
 8:            Qj,7 = Qj,7 + δphys (θfac )m ∆Qj,7                                    . Damp update
 9:        m=0                                                                 . Initialize count variable
                  (n+1)                                                                               ˜ θt
10:        while Qj,8 < 20 and (θfac )m > 0.01 do                                            . Limit Re
11:          m=m+1                                                  . Count number of damping iterations
               (n+1)    (n)                 (n)
12:          Qj,8 = Qj,8 + δphys (θfac )m ∆Qj,8                                           . Damp update

in the source terms.

Unsteady Residual Backtracking Line Search

In addition to the solution update damping, a backtracking line-search method was im-
plemented based on the work of Modisette [120]. Line searches increase reliability of the
nonlinear solver by expanding the global sphere of convergence of Newton’s method [120].
The principle of the unsteady residual backtracking line-search method is to damp the global
solution update in order to ensure a decrease in the unsteady residual at each time step,
which would be the case if Equation 3.1 was solved time accurately. The unsteady residual
is given by,

                             R(Q (n+1) )unst = ||T(n) ∆Q (n) + R(Q (n+1) )||2 .                    (3.11)

The global solution update is damped to ensure a decrease in the unsteady residual at each
iteration. This damping procedure is presented in Algorithm 4. The unsteady residual line
search is only attempted if the current global solution update damping factor is greater than
the lower bound of 1%. Similar to the physics-based restriction, if the damping factor falls
below 1% and if the reference time step is above the user-specified lower bound we attempt
to reject the solution and try again with a lower time step. If the reference time step is below
the lower bound, we continue with the solution with the backtracked residual and solution
variables.

3.1.4    Source-Term Time Step Restriction

Lian et al. [100] demonstrated the effects of treating source terms with implicit and explicit
Euler time marching by investigating the following linear convection equation with a source
term:
                                             ∂u    ∂u
                                                +a    = bu,                                        (3.12)
                                             ∂t    ∂x

<!-- ===== PDF page 65 ===== -->

3.1. PARALLEL IMPLICIT NEWTON-KRYLOV-SCHUR FLOW SOLVER                                                      45

Algorithm 4 Unsteady residual backtracking line search
 1: if δphys > 0.01 then
 2:     δls ← δphys                                       . Initialize damping factor with physics-based factor
           (n+1)
 3:    R(Q        )unst = ||T δls ∆Q + R(Q + δls ∆Q (n) )||2
                           (n)         (n)        (n)
                                                                                 . Calculate unsteady residual
 4:    while R(Q (n+1) )unst ≥ ||R(Q (n) )||2 and δls > 0.01 do                     . Reduce unsteady residual
 5:        δls = 0.90δls                                                              . Reduce damping factor
 6:        R(Q (n+1) )unst = ||T(n) δls ∆Q (n) + R(Q (n) + δls ∆Q (n) )||2       . Calculate unsteady residual
 7:    if δls < 0.01 then                                                               . Poor solution update
                 (n)
 8:        if ∆tref > ∆tref,min then                                         . Reject current update and retry
                  (n)               (n)
 9:            ∆tref = max (0.5∆tref , ∆tref,min )                                . Reduce reference time step
                 (n+1)       (n)
10:            Q        ←Q                                                               . Reset flow variables
11:            return                                                                          . Retry solution
12:        else                                                                       . Continue with solution

                                 -2

                                 -4

                                  -3         -2   -1     0       1       2        3

Figure 3.1: Normalized solution update over a range of time steps produced by implicit and explicit Euler
time marching applied to a linear convection equation with a source term. Reproduced from [100].

where a represents the wave speed and b is the scalar source term. The effect of the pseudo-
time-marching strategy on the normalized solution update for this scalar equation is illus-
trated in Figure 3.1. An implicit Euler time marching strategy applied to a source term
produces large solution updates, which asymptote at a source-term time step of unity, af-
ter which it produces solution updates with the opposite sign. These results can help to
explain the common problem of positivity preservation for turbulence and transition model-
dependent variables [120, 3]. Maintaining a source-term time step less than unity allows
for rapid solution growth while preventing unbounded and unstable solution updates. The
source-term time step is determined by the product of the largest positive eigenvalue of the
source-term Jacobian matrix, Asource , given by,

<!-- ===== PDF page 66 ===== -->

46                                                             CHAPTER 3. OPTIMIZATION FRAMEWORK
                                                                       
                                         ∂Sν̃        ∂Sν̃       ∂Sν̃
                                       ∂ ν̃
                                                    ∂γ          ˜ θt 
                                                              ∂ Re      
                                       ∂S           ∂Sγ        ∂Sγ 
                                            γ
                            A(n)    =                                   ,                 (3.13)
                                                                       
                             source
                                       ∂ ν̃
                                      
                                                     ∂γ          ˜
                                                              ∂ Reθt 
                                       ∂SRe        ∂SRe      ∂SRe
                                                                       
                                           ˜ θt        ˜ θt       ˜ θt 
                                                                      ,
                                         ∂ ν̃        ∂γ          ˜ θt
                                                              ∂ Re

                                                        (n)
and the local time step at node (j, k, m), λsource ∆tj,k,m . For the current work, the eigenvalues
of the block diagonal source-term Jacobian are determined using a QR algorithm. Similar to
the work performed by Allmaras [2], the approach taken to address the instabilities produced
by the linearization of the source terms is to reduce the local time step. However, in the
current work, a limit is placed on the time step based on the source-term strength. Specifi-
cally, the local time step is limited to ensure the following source-term time step restriction
is met:

                                                  (n)
                                      λsource ∆tj,k,m ≤ 0.8,                               (3.14)

where S represents the local source-term vector, in this case consisting of the SA turbulence
model and the SA-LM2015 or SA-sLM2015 transition model source terms. Several source-
term time step values were investigated; a value of 0.8 was determined to be an effective
balance between speed and robustness. More details on the spatially varying time step,
   (n)
∆tj,k,m , and the time step ramping procedure for the approximate- and inexact-Newton
phases can be found in [131].

   As demonstrated by Lian et al. [100], source-term eigenvalues can remain large throughout
convergence. This can cause the local time step to remain excessively small throughout
the solution process, delaying convergence. However, for most cases, the large sources are
only destabilizing in the early stages of convergence. To accelerate convergence, a switch
was developed for deactivating the source-term time step restriction in the inexact-Newton
phase. The source-term time step restriction is active throughout the approximate-Newton
phase, and is deactivated after five successive inexact-Newton iterations for which solution
update damping is not active (δphys = δls = 1). The source-term time step restriction is
reactivated if damping of the solution update is required, or if the relative residual drop for
the total residual (Equation 3.5) increases beyond the residual drop tolerance for the switch
from the approximate- to the inexact-Newton phase.

   It is important to note that in addition to the free-transition simulations investigated
in the current work, the source-term time step restriction has been found to improve the
robustness of fully turbulent simulations of complex geometries using the SA turbulence
model, such as in the design of transonic strut-braced wing and hybrid wing-body aircraft [22,
94].

<!-- ===== PDF page 67 ===== -->

3.2. GRADIENT-BASED OPTIMIZATION                                                             47

3.2     Gradient-based Optimization

3.2.1   Optimization Algorithm

The optimization algorithm consists of an integrated geometry parametrization and mesh
movement scheme developed by Hicken and Zingg [70] that represents the initial geometry
by a set of B-spline surface patches. Shape control is achieved using free-form deformation
(FFD) B-spline volumes with an axial deformation geometry control system [55, 95]. Opti-
mization is performed using a gradient-based optimization strategy where the gradients are
calculated using the discrete-adjoint approach [70, 128]. The flow adjoint system is solved
using a simplified and flexible variant of GCROT [71], while the mesh adjoint system is
solved using a preconditioned conjugate-gradient method. The design variables are updated
using the sparse quadratic programming algorithm SNOPT [57], which is capable of handling
both linear and nonlinear constraints and has been demonstrated to be efficient for problems
with large numbers of design variables.

3.2.2   Gradient Evaluation

The discrete-adjoint Lagrangian function for the PDE-constrained optimization problem
containing the flow, ψ f , and mesh adjoint variables, ψ m , is given by,

        L(Q, C , X , ψ f , ψ m ) = J(Q, C , X ) − ψ Tf R(Q, C , X ) − ψ Tm M (C , X ),   (3.15)

where Q, C , and X represent the array of conserved flow variables, B-spline control points,
and design variables, and R and M represent the total residual and the mesh movement
equations, respectively. Setting the first derivatives of the Lagrangian with respect to Q,
C , X , ψ f , and ψ m to zero, we recover the following first-order optimality conditions:

                              ∂L
                                   = 0 = M (C , X ),                                     (3.16)
                             ∂ψ m
                              ∂L
                                   = 0 = R(Q, C , X ),                                   (3.17)
                              ∂ψ f
                              ∂L         ∂J        ∂R
                                   =0=      − ψ Tf     ,                                 (3.18)
                              ∂Q         ∂Q        ∂Q
                              ∂L         ∂J        ∂R           ∂M
                                   =0=      − ψ Tf      − ψ Tm     ,                     (3.19)
                              ∂C         ∂C        ∂C           ∂C
                              ∂L         ∂J        ∂R           ∂M
                                   =0=      − ψ Tf       − ψ Tm    ,                     (3.20)
                              ∂X         ∂X        ∂X           ∂X
where in the current work the partial derivatives are formed analytically using hand lin-
earization with some entries determined using the complex-step method.
   To perform free-transition optimization, the partial derivatives of the total residual must
                                             ∂R                                   ∂R
be updated. Specifically, the flow Jacobian, ∂Q , and grid metric linearization, ∂C  , used for
the flow and mesh adjoint systems (Equations 3.18 and 3.19), respectively, are updated to

<!-- ===== PDF page 68 ===== -->

48                                                             CHAPTER 3. OPTIMIZATION FRAMEWORK

              Figure 3.2: Infinite swept wing geometry with FFD design variables in blue.

include the modifications to the SA turbulence model and the transition model equations,
enabling free-transition gradient-based aerodynamic shape optimization. Since first-order
upwinding is applied to the turbulence and transition model convective terms, the analytical
linearization of these equations for the approximate Jacobian, which is used in the global-
ization phase of the flow solver and to form the preconditioner, is sufficient. As these partial
derivatives are formed analytically, it is important to verify their implementation. The ana-
lytical flow Jacobian and metric linearization are each verified through comparing with the
complex-step method, while the analytical gradient of the objective function is verified using
a second-order finite-difference approximation. It is important to note that the axial and
FFD geometric design variables are formed using the B-spline control points, C , [55] and
that the residuals depend on the B-spline control points through the grid metrics and the
off-wall spacing, the latter of which is used in the turbulence and transition model source
terms. Furthermore, partial derivatives involving the angle of attack and sideslip design
                                                                                          ∂R
variables are determined using the complex-step method. Therefore, the analytical ∂X          is
                            ∂R
verified concurrently with ∂C .
   It is important that the verification tests include all relevant flow features. To ensure this,
an infinite swept wing geometry consisting of a blunt trailing-edge RAE2822 airfoil extruded
with a 25-degree sweep and periodic boundary conditions is simulated using the smooth
transition model with the local helicity-based crossflow correlations and compressibility cor-
rections, SA-sLM2015cc, presented in Chapter 2. The infinite swept wing is simulated at a
Mach number of 0.785, Reynolds number of 20.3 × 106 , lift coefficient of 0.56, turbulence
intensity of 0.07%, and with an assumed surface roughness of 1.8µm. The infinite swept
wing O-grid consists of 561, 121, and 11 nodes in the streamwise, off-wall, and spanwise di-
rections, respectively, with average and maximum y+ values of 0.25 and 0.50. The geometry
is parameterized using B-spline surface patches that are embedded in an FFD volume with
6 FFD design variables on each of the upper and lower surfaces, as presented in Figure 3.2.
The design variables are constrained to be constant along the span and the leading- and
trailing-edge design variables are constrained to move symmetrically, resulting in 10 effec-

<!-- ===== PDF page 69 ===== -->

3.2. GRADIENT-BASED OPTIMIZATION                                                                                                                                                                                                                                                                           49

                                                                                        Flow Jacobian                                                                                                                                                                                                Drag
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
                                                                                                                                                                                                                                         10-8 -2        -4        -6            -8        -10        -12
              10        10        10        10        10         10         10         10         10         10                            10        10        10         10         10         10         10         10                   10      10        10            10        10         10
                                                             ε                                                                                                                  ε                                                                                      ε
                                            ∂R                                                                                                                       ∂R
                         (a) Flow Jacobian: ∂Q                                                                                             (b) Metric Linearization: ∂C                                                                    (c) Directional Derivative: Dz J

Figure 3.3: Verification for the analytical flow Jacobian, metric linearization, and directional derivative
including the SA-sLM2015cc transition model equations using both complex-step and finite-difference ap-
proximations.

tive geometric design variables plus angle of attack. Although the current test does not
involve axial control system design variables, the two-level axial and FFD control system
was not modified by the inclusion of the transition model equations. The changes to the
partial derivatives that were made to enable free-transition optimization can be verified by
considering only the FFD geometric design variables and angle of attack.
   The verification cases are performed using a fully converged SA-sLM2015cc simulation
on the infinite swept wing to ensure that the turbulence and transition model source terms
are active. The results for the verification study are presented in Figure 3.3. Machine-zero
agreement is demonstrated for both the analytical flow Jacobian and metric linearization
with approximations formed using the complex-step method. A directional derivative is
used to verify the analytical gradient components simultaneously [68, 127], defined as,

                                                                                                                                                                                     ∂J
                                                                                                                                                      Dz J =                            z,                                                                                                      (3.21)
                                                                                                                                                                                     ∂X
and the second-order finite-difference approximation is,

                                                                                                                                 J(X + z) − J(X − z)
                                                                                                        Dz J =                                         + O(2 ),                                                                                                                                (3.22)
                                                                                                                                           2
where  is the perturbation parameter. The difference between the values produced by
the analytical directional derivative and the second-order finite-difference approximation is
presented in Figure 3.3c. Good agreement is demonstrated to approximately 6 and 7 orders
of magnitude with the objective function set to lift and drag, respectively. Figure 3.3c
demonstrates increased round-off error with large step sizes, which is a consequence of using
a finite-difference approximation.

<!-- ===== PDF page 71 ===== -->

Chapter 4

Results: Analysis

This chapter presents results for the SA-LM2015 and SA-sLM2015 transition models ap-
plied to subsonic flow regimes (Section 4.1) [139], and the SA-sLM2015 and SA-sLM2015cc
transition models applied to transonic flow regimes (Section 4.2) [142]. The focus of the
subsonic transition test cases is to validate the SA-LM2015 transition model and the new
smooth variant, SA-sLM2015, in incompressible flow regimes, as well as to demonstrate the
improved numerical behaviour of the latter. The transonic transition test cases are used
to investigate the compressibility corrections for the TS and stationary crossflow instability
empirical correlations applied in the SA-sLM2015cc transition model, and to evaluate the
results relative to the original, incompressible correlations in the smooth transition model,
SA-sLM2015. For the results presented, the simulated transition locations are identified by
the sharp rise in the laminar skin friction coefficient.
    A second-order discretization is used for the mean-flow equations with artificial viscosity
added using the matrix-based dissipation model of Swanson and Turkel [184]. Limiters are
provided in the matrix-based dissipation model to prevent zero artificial viscosity at stagna-
tion points and sonic lines, where the linear and nonlinear eigenvalue scaling factors approach
zero, respectively [184]. The linear eigenvalue limiter parameter Vl , is set to zero, as it was
found to be overly dissipative in the laminar boundary layer, while the nonlinear eigenvalue
limiter parameter, Vn , is set to 0.25, which is between the upper and lower bound of 0.20
and 0.30 recommended by Swanson and Turkel [184]. For the subsonic transition test cases,
a fourth-difference dissipation coefficient of 0.04 is used, while for the transonic airfoil test
cases second- and fourth-difference coefficients of 2.00 and 0.04 are used, respectively. Due
to the increased complexity of the NASA CRM-NLF wing-body geometry, the dissipation
coefficients are increased to 3.00 and 0.06 for the second- and fourth-difference coefficients,
respectively, and to 0.30 for the nonlinear eigenvalue limiter, Vn . In addition, the QCR2000
correction [173] is used for simulations of the NASA CRM-NLF. First-order upwinding is
used for the turbulence and transition model convective terms, with second-order centered
differencing used for the diffusion and source terms.
   Lopes et al. [103] performed a thorough analysis of the effects of discretization strategy for

<!-- ===== PDF page 72 ===== -->

52                                                                 CHAPTER 4. RESULTS: ANALYSIS

the turbulence and transition model convective terms on solution accuracy by investigating
the SST turbulence model coupled to several transport-equation-based transition models.
The results demonstrate that the discretization strategy for the transition model equations
has a small influence on the solution accuracy, but the SST turbulence model discretization
strategy has a strong impact. However, this may be due to a difference in the freestream
decay of the turbulent kinetic energy and therefore the turbulence intensity, which has a sig-
nificant impact on the transition location [102]. This is avoided in the current work using the
SA turbulence model and a far-field turbulence intensity in the empirical correlations [139].
   In addition, Swanson and Turkel recommend fourth-difference dissipation coefficients be-
        1      1
tween 64  and 32 [184], which are smaller than the values of 0.04 and 0.06 used in the current
work. However, through participation at the Fifth AIAA Drag Prediction Workshop Osusky
and Zingg [131, 130] demonstrated excellent grid convergence using a fourth-difference coeffi-
cient of 0.04 and first-order upwinding applied to the SA turbulence model convection term,
with the grid-converged value within 1 drag count of the median obtained by all workshop
participants. Care is taken for each test case to ensure all results are sufficiently iteratively
and grid converged to produce the low levels of numerical error required to investigate the
modelling error and validity of the transition model variants. The numerical behaviour of
the smooth transition model, including both iterative and grid convergence, is investigated
in more detail for the subsonic NLF0416 airfoil, transonic VA-2 airfoil, and transonic CRM-
NLF geometry and evaluated relative to fully turbulent simulations performed using the
SA turbulence model in Appendix B. Grid convergence results are plotted versus N −1 and
N −2/3 for two- and three-dimensional simulations, respectively, where N is the total number
of compute nodes, which is proportional to the grid spacing squared.
                                                                                      (n)
   Residual convergence histories are illustrated by the normalized total residual, Rd (Equa-
tion 3.5), versus the computational cost, which is measured in units of equivalent residual
evaluations. This is determined as the wall-clock time normalized by the average cost of
computing the total residual for each simulation. Residual drop tolerances of 10−4 and 10−5
are used for the switch from the approximate- to the inexact-Newton phase of the Newton-
Krylov solver for the two- and three-dimensional simulations, respectively, while relative
and absolute convergence tolerances are specified as 10−15 and 10−10 . For the majority of
the results presented the absolute residual tolerance is met first. Similar solution algorithm
parameters to the values recommended by Osusky and Zingg [129] are used for the two-
and three-dimensional test cases, with the settings for the source-term time step restriction,
solution update damping, and scaling measures described previously in Chapter 3.

4.1    Subsonic Transition Test Cases

The SA-LM2015 and SA-sLM2015 transition models are validated in this section using a
range of two- and three-dimensional subsonic transition test cases. These test cases include

<!-- ===== PDF page 73 ===== -->

4.1. SUBSONIC TRANSITION TEST CASES                                                                  53

Table 4.1: TMPCS NLF0416 structured C-grid dimensions [27], where ∆s represents the first off-wall grid
spacing.
                grid level   chord × off-wall nodes   ∆s (chord)    average/maximum y+
                 coarse            529 × 73           4.47 × 10−6         0.57/1.19
                 medium            705 × 97           3.25 × 10−6         0.51/0.89
                   fine           1057 × 145          2.10 × 10−6         0.26/0.45

natural and separation-induced transition as well as transition due to stationary crossflow
instabilities, with the results compared to experimental values.

4.1.1    NLF0416 General Aviation Airfoil

The NLF0416 airfoil is a natural-laminar-flow general aviation airfoil designed by
Somers [171]. Experimental results of the NLF0416 airfoil were obtained in the low turbu-
lence pressure tunnel (LTPT) at NASA Langley Research Center [171]. Transition locations
were measured for several angles of attack at a Mach number of 0.1 and a Reynolds number
of 4.0 × 106 . The transition mechanisms for these conditions are natural and separation-
induced transition. The turbulence intensity in the wind tunnel used for the experimental
study was not specified; however a value of T u = 0.15% is assumed in order to be consistent
with results obtained using the Transition Modelling and Predictive Capabilities Seminar
(TCMPS) guidelines [27].
   The coarse, medium, and fine NLF0416 family of grids provided by the TMPCS com-
mittee [27] is used for the results presented. The dimensions of the structured C-grids are
presented in Table 4.1. A grid-refinement study is conducted for the 0◦ angle of attack case
using both transition models with the results presented in Figure 4.1. The results demon-
strate that both models produce similar grid convergence, although the two models appear
to converge to slightly different lift and drag values, which is likely due to the difference
in the predicted transition locations demonstrated in Figure 4.1b. A more thorough in-
vestigation of the grid convergence of the SA-sLM2015 transition model for this test case
is presented in Appendix B. Figure 4.1 also illustrates the computational results obtained
using the SA-LM2015 and SA-sLM2015 transition models on the ‘fine’ grid level compared
to experimental values. The results demonstrate that both the SA-LM2015 model and the
SA-sLM2015 model do a reasonable job of matching the forces and transition locations from
the experiment. Figure 4.1e illustrates that both models under-predict drag at low angles of
attack; however this is consistent with results obtained from previous studies using variations
of the γ-Re˜ θt model [76, 40, 80] and does not appear to be caused by incorrect transition
onset locations.
   Residual convergence histories for these cases are presented in Figure 4.2. The convergence
histories demonstrate the improved numerical behaviour of the SA-sLM2015 model. At
angles of attack of −6◦ and −4◦ , a large separation bubble forms on the lower surface of the
airfoil, increasing the stiffness of the linear system. This is evident in the decreased slope

<!-- ===== PDF page 74 ===== -->

54                                                                                                        CHAPTER 4. RESULTS: ANALYSIS

                                 SA-LM2015fine                                                            SA-LM2015cl
        0.008                    SA-LM2015medium                           .0055                          SA-LM2015cd                  .4910
                                 SA-LM2015coarse                                                          SA-sLM2015cl
                                 SA-sLM2015fine                            .0054                          SA-sLM2015cd                 .4900
        0.006                    SA-sLM2015medium
                                 SA-sLM2015coarse                          .0053                                                       .4890

        0.004                                                              .0052                                                       .4880

                                                                           Cd

                                                                                                                                          Cl
                             upper
        Cf

                                                                           .0051                                                       .4870
        0.002
                                                                           .0050                                                       .4860
                                                       lower
                                                                           .0049                                                       .4850

                                                                           .0048
                   0       0.2       0.4         0.6         0.8       1           0         1E-05              2E-05         3E-05
                                           x/c                                                            N-1
                       (a) Skin friction grid convergence                         (b) Grid convergence for lift and drag

                                                         SA-LM2015fine                                 SA-LM2015fine
                                                         SA-sLM2015fine         1.5                    SA-sLM2015fine
              2                                          LTPT Experiment                               LTPT Experiment

             1.5                                                                 1

                                                                                0.5
        Cl

                                                                           Cl

             0.5

                                                                            -0.5
                         .0050     .0100         .0150         .0200                        -5             0              5
                                           Cd                                                             aoa

                                 (c) Drag polar                                                   (d) Lift curve

                                                         SA-LM2015fine                                                  SA-LM2015fine
                                                         SA-sLM2015fine                                                 SA-sLM2015fine
        .0090                                            LTPT Experiment        1.5                                     LTPT Experiment

        .0080                                                                    1
                                                                                       upper                                   lower
        Cd

                                                                                0.5
                                                                           Cl

        .0070

        .0060

                                                                            -0.5
        .0050

                            -5             0             5                             0         0.2      0.4       0.6       0.8
                                           aoa                                             transition location (x/c)
                                 (e) Drag curve                                            (f) Transition locations

Figure 4.1: NLF0416 grid convergence at 0◦ angle of attack, and angle-of-attack sweep simulated using the
SA-LM2015 and SA-sLM2015 transition models at M = 0.1, Re = 4.0 × 106 , and T u = 0.15%. The ‘fine’
grid level is used for the sweep, with the experimental results obtained at the LTPT wind tunnel [171].

<!-- ===== PDF page 75 ===== -->

4.1. SUBSONIC TRANSITION TEST CASES                                                                                                                                                                                       55

                        101                                                                         101                                                                         101
                                                           SA-LM2015                                                                   SA-LM2015                                                             SA-LM2015
                        10−1                               SA-sLM2015                               10−1                               SA-sLM2015                               10−1                         SA-sLM2015

Relative Residual, d

                                                                            Relative Residual, d

                                                                                                                                                        Relative Residual, d
                        10−3                                                                        10−3                                                                        10−3
                        10−5                                                                        10−5                                                                        10−5
                        10−7                                                                        10−7                                                                        10−7
                        10−9                                                                        10−9                                                                        10−9
                   10−11                                                                       10−11                                                                       10−11
                   10−13                                                                       10−13                                                                       10−13
                   10−15       0.0      0.5   1.0   1.5   2.0   2.5   3.0                      10−15       0.0      0.5   1.0   1.5   2.0   2.5   3.0                      10−15 0.0      0.5   1.0   1.5   2.0   2.5   3.0
                                     Equivalent Residual Evaluations 1e5                                         Equivalent Residual Evaluations 1e5                                   Equivalent Residual Evaluations 1e5
                                        (a) aoa = −6◦                                                               (b) aoa = −4◦                                                          (c) aoa = −2◦
                        101                                                                         101                                                                         101
                                                           SA-LM2015                                                                   SA-LM2015                                                             SA-LM2015
                        10−1                               SA-sLM2015                               10−1                               SA-sLM2015                               10−1                         SA-sLM2015
Relative Residual, d

                                                                            Relative Residual, d

                                                                                                                                                        Relative Residual, d
                        10−3                                                                        10−3                                                                        10−3
                        10−5                                                                        10−5                                                                        10−5
                        10−7                                                                        10−7                                                                        10−7
                        10−9                                                                        10−9                                                                        10−9
                   10−11                                                                       10−11                                                                       10−11
                   10−13                                                                       10−13                                                                       10−13
                   10−15 0.0            0.5   1.0   1.5   2.0   2.5   3.0                      10−15 0.0            0.5   1.0   1.5   2.0   2.5   3.0                      10−15 0.0      0.5   1.0   1.5   2.0   2.5   3.0
                                     Equivalent Residual Evaluations 1e5                                         Equivalent Residual Evaluations 1e5                                   Equivalent Residual Evaluations 1e5
                                          (d) aoa = 0◦                                                                (e) aoa = 2◦                                                          (f) aoa = 4◦
                                                                                                                                       SA-LM2015
                                                                                                    10−1                               SA-sLM2015
                                                                            Relative Residual, d

                                                                                                    10−3
                                                                                                    10−5
                                                                                                    10−7
                                                                                                    10−9
                                                                                               10−11
                                                                                               10−13
                                                                                               10−15 0.0            0.5   1.0   1.5   2.0   2.5   3.0
                                                                                                                 Equivalent Residual Evaluations 1e5
                                                                                                                      (g) aoa = 6◦

Figure 4.2: NLF0416 residual convergence histories for the SA-LM2015 and SA-sLM2015 transition models
on the ‘fine’ grid level. The airfoil is simulated at M = 0.1, Re = 4.0 × 106 , T u = 0.15% and over a range
of angles of attack.

of the residual curve for the SA-sLM2015 model, as each outer inexact-Newton iteration
required more inner GMRES linear solver iterations. It should also be noted that small
trailing-edge spacings in the C-grids provided by the TMPCS committee propagate to the
far field of the domain with little splaying, further increasing the stiffness of the linear system
and the computational cost of the simulations.

4.1.2                                 S809 Wind Turbine Airfoil

The S809 airfoil is a natural-laminar-flow wind turbine airfoil designed by Somers [172] and
studied in the low-turbulence wind tunnel of the Delft University of Technology Low-Speed
Laboratory. Experimental results for the S809 airfoil were obtained at a Reynolds number
of 2.0 × 106 at angles of attack ranging from 0◦ to 9◦ . The freestream turbulence intensity
in the wind tunnel varied between T u = 0.02% and T u = 0.04%. For a Reynolds number of

<!-- ===== PDF page 76 ===== -->

56                                                                     CHAPTER 4. RESULTS: ANALYSIS

                     Table 4.2: TMPCS S809 structured C-grid dimensions [27].
              grid level   chord × off-wall nodes   ∆s (chord)    average/maximum y+
               coarse            529 × 73           8.90 × 10−6         0.60/1.73
               medium            705 × 97           6.47 × 10−6         0.52/0.81
                 fine           1057 × 145          4.19 × 10−6         0.28/0.54

2.0 × 106 and an angle of attack of 0◦ , Somers noted that the primary transition mechanism
on both the upper and lower surfaces of the airfoil was separation-induced transition due
to the presence of large laminar separation bubbles [172]. As the angle of attack increased,
the separation bubble on the upper surface decreased in length, disappearing completely
at an angle of attack of 5.13◦ , with the transition location moving upstream rapidly on
the upper surface and the lower surface separation bubble fixing transition near mid-chord
over the range of angles of attack. Numerical results are obtained using a Mach number
of 0.10, Reynolds number of 2.0 × 106 , and a freestream turbulence intensity of 0.07%, as
recommended by the TCMPS guidelines [27].

   Computational results are obtained using the S809 grid families supplied through the
TCMPS, which were generated using the same conformal mapping method as for the
NLF0416 airfoil [27], with the coarse, medium, and fine grid dimensions presented in Ta-
ble 4.2. A grid-refinement study is performed for the S809 airfoil at 6◦ angle of attack with
the results presented in Figure 4.3. This condition is investigated as it represents the point
beyond the upper corner of the drag bucket where the transition location is rapidly moving,
which is expected to amplify the effects of grid refinement on the predictions [27]. Figure 4.3
also illustrates the forces and transition locations obtained using the SA-LM2015 and SA-
sLM2015 transition models on the ‘fine’ grid level. Although the transition locations match
the experiment well, both models under-predict drag before and after the upper corner of
the drag bucket. Again, this trend has been identified by previous studies using variations
of the γ-Re˜ θt model [76, 40, 80], as well as studies using the AFT transition model [40, 28],
and does not appear to be caused by inaccurate transition prediction. This discrepancy may
be due to differences between the free-air simulations and the wind tunnel environment.

   Residual convergence histories for these cases are presented in Figure 4.4. Similar to
the NLF0416 airfoil, the large separation bubbles at low angles of attack lead to stiff linear
systems which delay convergence. For all angles of attack the SA-sLM2015 model displays
substantially improved convergence. In particular, for the 9◦ angle of attack case the SA-
LM2015 model is unable to converge due to the unsteady interaction of the transition model
with the laminar separation bubble on the upper surface of the airfoil. This behaviour was
also identified by Denison and Pulliam using the LM2015 model [40]. Limiting vorticity in
the transition model source terms prevents this oscillation without affecting the accuracy of
the solution.

<!-- ===== PDF page 77 ===== -->

4.1. SUBSONIC TRANSITION TEST CASES                                                                                                                    57

                                                        SA-LM2015fine                                            SA-LM2015cl
        0.012                                           SA-LM2015medium                                          SA-LM2015cd
                                                        SA-LM2015coarse     .0094                                SA-sLM2015cl
         0.01                                           SA-sLM2015fine                                           SA-sLM2015cd                  .8700
                                                        SA-sLM2015medium    .0092
        0.008                                           SA-sLM2015coarse
                                                                            .0090
                        upper                                                                                                                  .8600
        0.006                                                               .0088

                                                                            Cd
        Cf

                                                                                                                                                  Cl
        0.004                                                               .0086
                                                                                                                                               .8500
                                                                            .0084
        0.002
                                                 lower                      .0082
              0                                                                                                                                .8400
                                                                            .0080
        -0.002
                   0       0.2     0.4          0.6        0.8      1               0                1E-05             2E-05           3E-05
                                          x/c                                                                    N-1
                       (a) Skin friction grid convergence                           (b) Grid convergence for lift and drag

                                                      SA-LM2015fine                                     SA-LM2015fine
                                                      SA-sLM2015fine                                    SA-sLM2015fine
                                                      TU Delft Experiment                               TU Delft Experiment
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
                                  .0100                   .0150                         0        2           4         6         8      10
                                          Cd                                                                     aoa

                                 (c) Drag polar                                                         (d) Lift curve

                                                      SA-LM2015fine                                                          SA-LM2015fine
                                                      SA-sLM2015fine              12                                         SA-sLM2015fine
                                                      TU Delft Experiment                                                    TU Delft Experiment
        .0160                                                                     10

                                                                                                      upper                    lower
        .0140                                                                      8
                                                                            aoa

        .0120
        Cd

        .0100

        .0080
        .0060
                            2      4            6           8                               0     0.2        0.4       0.6       0.8     1
                                         aoa                                                    transition location (x/c)
                                 (e) Drag curve                                                 (f) Transition locations

Figure 4.3: S809 grid convergence at 6◦ angle of attack, and angle-of-attack sweep simulated using the SA-
LM2015 and SA-sLM2015 transition models at M = 0.1, Re = 2.0 × 106 , and T u = 0.07%. The ‘fine’ grid
level is used for the sweep, with the experimental results obtained at the TU Delft wind tunnel [172].

<!-- ===== PDF page 78 ===== -->

58                                                                                                                                                                              CHAPTER 4. RESULTS: ANALYSIS

                        101                                                                         101                                                                         101
                                                           SA-LM2015                                                                   SA-LM2015                                                             SA-LM2015
                        10−1                               SA-sLM2015                               10−1                               SA-sLM2015                               10−1                         SA-sLM2015
Relative Residual, d

                                                                            Relative Residual, d

                                                                                                                                                        Relative Residual, d
                        10−3                                                                        10−3                                                                        10−3
                        10−5                                                                        10−5                                                                        10−5
                        10−7                                                                        10−7                                                                        10−7
                        10−9                                                                        10−9                                                                        10−9
                   10−11                                                                       10−11                                                                       10−11
                   10−13                                                                       10−13                                                                       10−13
                   10−15       0.0      0.5   1.0   1.5   2.0   2.5   3.0                      10−15       0.0      0.5   1.0   1.5   2.0   2.5   3.0                      10−15 0.0      0.5   1.0   1.5   2.0   2.5   3.0
                                     Equivalent Residual Evaluations 1e5                                         Equivalent Residual Evaluations 1e5                                   Equivalent Residual Evaluations 1e5
                                          (a) aoa = 1◦                                                                (b) aoa = 2◦                                                          (c) aoa = 3◦
                        101                                                                         101                                                                         101
                                                           SA-LM2015                                                                   SA-LM2015                                                             SA-LM2015
                        10−1                               SA-sLM2015                               10−1                               SA-sLM2015                               10−1                         SA-sLM2015
Relative Residual, d

                                                                            Relative Residual, d

                                                                                                                                                        Relative Residual, d
                        10−3                                                                        10−3                                                                        10−3
                        10−5                                                                        10−5                                                                        10−5
                        10−7                                                                        10−7                                                                        10−7
                        10−9                                                                        10−9                                                                        10−9
                   10−11                                                                       10−11                                                                       10−11
                   10−13                                                                       10−13                                                                       10−13
                   10−15 0.0            0.5   1.0   1.5   2.0   2.5   3.0                      10−15 0.0            0.5   1.0   1.5   2.0   2.5   3.0                      10−15 0.0      0.5   1.0   1.5   2.0   2.5   3.0
                                     Equivalent Residual Evaluations 1e5                                         Equivalent Residual Evaluations 1e5                                   Equivalent Residual Evaluations 1e5
                                          (d) aoa = 4◦                                                                (e) aoa = 5◦                                                          (f) aoa = 6◦
                        101                                                                         101                                                                         101
                                                           SA-LM2015                                                                   SA-LM2015                                                             SA-LM2015
                        10−1                               SA-sLM2015                               10−1                               SA-sLM2015                               10−1                         SA-sLM2015
Relative Residual, d

                                                                            Relative Residual, d

                                                                                                                                                        Relative Residual, d
                        10−3                                                                        10−3                                                                        10−3
                        10−5                                                                        10−5                                                                        10−5
                        10−7                                                                        10−7                                                                        10−7
                        10−9                                                                        10−9                                                                        10−9
                   10−11                                                                       10−11                                                                       10−11
                   10−13                                                                       10−13                                                                       10−13
                   10−15       0.0      0.5   1.0   1.5   2.0   2.5   3.0                      10−15       0.0      0.5   1.0   1.5   2.0   2.5   3.0                      10−15 0.0      0.5   1.0   1.5   2.0   2.5   3.0
                                     Equivalent Residual Evaluations 1e5                                         Equivalent Residual Evaluations 1e5                                   Equivalent Residual Evaluations 1e5
                                          (g) aoa = 7◦                                                                (h) aoa = 8◦                                                          (i) aoa = 9◦

Figure 4.4: S809 residual convergence histories for the SA-LM2015 and SA-sLM2015 transition models on
the ‘fine’ grid level. The airfoil is simulated at M = 0.1, Re = 2.0 × 106 , T u = 0.07% and over a range of
angles of attack.

4.1.3                                 NASA NLF2-0415 Infinite Swept Wing

The NASA NLF2-0415 Infinite Swept Wing consists of an NLF2-0415 airfoil extruded with
a 45◦ sweep angle. Experimental measurements were obtained at an angle of attack of −4◦ ,
over a range of Mach and Reynolds numbers [38], and including a variety of surface rough-
ness levels [148]. Transition locations were determined using the naphthalene visualization
technique. Transition locations for a roughness level, h, of 3.3µm, representing a painted
surface, were obtained using the SA-LM2015 and SA-sLM2015 models and are compared
with experimental values [148] in Figure 4.5. The Reynolds number based on freestream
velocity magnitude and chord length was varied from 1.8 × 106 to 3.5 × 106 . The freestream
turbulence intensity for the wind tunnel was measured as 0.20%. Transition locations are
presented only for the upper surface of the wing, as experimental results are not available
for the lower surface.

<!-- ===== PDF page 79 ===== -->

4.1. SUBSONIC TRANSITION TEST CASES                                                                                                       59

                                                                       SA-LM2015                                             SA-LM2015
                                                                       SA-sLM2015   0.006                                    SA-sLM2015
                                                                       Experiment
                                    0.8

        transition location (x/c)   0.6                                             0.004

                                                                                    Cf
                                    0.4
                                                                                    0.002

                                          2E+06             3E+06                           0              0.5                    1
                                              Reynolds number                                              x/c
                                            (a) Transition locations                            (b) Skin friction profiles

Figure 4.5: Experimental and computational results for the NLF2-0415 infinite swept wing [148] simulated
using the SA-LM2015 and SA-sLM2015 transition models at M = 0.15, α = −4◦ , T u = 0.20%, and over a
range of Reynolds numbers from 1.8–3.5 × 106 .

   The computational grid used to obtain these results consists of 270,336 nodes with 90,112
nodes per cross-section, and an off-wall spacing of 2.5 × 10−6 chords, resulting in average
and maximum y+ values of 0.21 and 0.46, respectively. The analytical full Jacobian is used
in the inexact-Newton phase for this case and the TU Braunschweig Sickle Wing. Periodic
boundary conditions are used to simulate the infinite swept wing. The results in Figure 4.5
demonstrate that both transition models predict transition accurately, while plots of the
residual convergence histories in Figure 4.6 demonstrate that the smooth variant provides
significantly improved convergence to steady state.
   The results in Figure 4.5 demonstrate that the SA-LM2015 model is able to closely predict
the location of transition for all Reynolds numbers except for 2×106 , where there is an error of
approximately 14% between the experimental and predicted values. This error is also present
in the work of Jung and Baeder [76], where they coupled a semi-local correlation-based
transition model formulation with the SA model using the crossflow correlations developed
by Müller and Herbst [125], and in the original model developed by Langtry et al. [93]. In
the work produced by Radeztsky et al. [148], the results demonstrate that there is significant
experimental error associated with this case, which is reflected in the error bars in Figure 4.5
and may help to explain this discrepancy.
   It is important to note that the choice of numerical dissipation strongly affects the ac-
curacy of the crossflow transition model. Both the scalar dissipation model developed by
Jameson et al. [75] and later refined by Pulliam [146], and the matrix-based dissipation model
developed by Swanson and Turkel [184] were investigated for this case. On the grid used,
the scalar dissipation model drastically under-predicts the crossflow strength and helicity
values. This was determined to be due to the increased dissipation of the scalar dissipation
model which affected the vorticity profile in the laminar boundary layer.

<!-- ===== PDF page 80 ===== -->

60                                                                                                                                                                              CHAPTER 4. RESULTS: ANALYSIS

                        101                                                                         101                                                                         101
                                                           SA-LM2015                                                                   SA-LM2015                                                             SA-LM2015
                        10−1                               SA-sLM2015                               10−1                               SA-sLM2015                               10−1                         SA-sLM2015
Relative Residual, d

                                                                            Relative Residual, d

                                                                                                                                                        Relative Residual, d
                        10−3                                                                        10−3                                                                        10−3
                        10−5                                                                        10−5                                                                        10−5
                        10−7                                                                        10−7                                                                        10−7
                        10−9                                                                        10−9                                                                        10−9
                   10−11                                                                       10−11                                                                       10−11
                   10−13                                                                       10−13                                                                       10−13
                   10−15       0.0      0.5   1.0   1.5   2.0   2.5   3.0                      10−15       0.0      0.5   1.0   1.5   2.0   2.5   3.0                      10−15 0.0      0.5   1.0   1.5   2.0   2.5   3.0
                                     Equivalent Residual Evaluations 1e5                                         Equivalent Residual Evaluations 1e5                                   Equivalent Residual Evaluations 1e5
                                      (a) Re = 1.8 × 106                                                          (b) Re = 2.0 × 106                                                    (c) Re = 2.2 × 106
                        101                                                                         101                                                                         101
                                                           SA-LM2015                                                                   SA-LM2015                                                             SA-LM2015
                        10−1                               SA-sLM2015                               10−1                               SA-sLM2015                               10−1                         SA-sLM2015
Relative Residual, d

                                                                            Relative Residual, d

                                                                                                                                                        Relative Residual, d
                        10−3                                                                        10−3                                                                        10−3
                        10−5                                                                        10−5                                                                        10−5
                        10−7                                                                        10−7                                                                        10−7
                        10−9                                                                        10−9                                                                        10−9
                   10−11                                                                       10−11                                                                       10−11
                   10−13                                                                       10−13                                                                       10−13
                   10−15 0.0            0.5   1.0   1.5   2.0   2.5   3.0                      10−15 0.0            0.5   1.0   1.5   2.0   2.5   3.0                      10−15 0.0      0.5   1.0   1.5   2.0   2.5   3.0
                                     Equivalent Residual Evaluations 1e5                                         Equivalent Residual Evaluations 1e5                                   Equivalent Residual Evaluations 1e5
                                      (d) Re = 2.5 × 106                                                          (e) Re = 3.0 × 106                                                    (f) Re = 3.5 × 106

Figure 4.6: NLF2-0415 residual convergence histories for the SA-LM2015 and SA-sLM2015 transition models.
The infinite swept wing is simulated at M = 0.15, α = −4◦ , T u = 0.20%, and over a range of Reynolds
numbers.

4.1.4                                 TU Braunschweig Sickle Wing

The TU Braunschweig Sickle Wing [138] was developed as a test case for crossflow transition
on a complex geometry. The model consists of five main sections designed to create three-
dimensional boundary layers with increasing crossflow towards the wing tip. The first section
raises the model away from the turbulent boundary layer of the wind tunnel wall, with the
three main wing sections swept at 30, 45, and 55◦ with distinct kinks in the leading edge
sweep. As a result, large spanwise gradients are present and the assumptions of linear local
stability theory are challenged [138]. The wing is capped with a rounded wing tip. At
the experimental flow conditions, the upper surface of the wing has extended regions of
strong favourable pressure gradients with crossflow instabilities developing at the section
kinks. This results in mixed mode transition scenarios with TS and crossflow instabilities
triggering transition, except on the upper surface near the root of the wing, where a laminar
separation bubble triggers transition at approximately 80% chord.
   The Sickle Wing is simulated at a Reynolds number of 2.75 × 106 , Mach number of 0.16,
and −2.6◦ angle of attack. The RMS surface roughness of the experimental model was
measured to be 1.47µm with a five-point average of the maximum peak to peak height of
9.78µm [138]. In order to be consistent with results obtained using the LM2015 model [93],
the surface roughness for the simulations is set to the peak to peak height value of 9.78µm
with a freestream turbulence intensity of 0.17%, which is consistent with the turbulence
intensity measured in the wind tunnel [138].

<!-- ===== PDF page 81 ===== -->

4.2. TRANSONIC TRANSITION TEST CASES                                                            61

                         Table 4.3: Sickle Wing structured O-O grid dimensions.
       grid level   total node count, N   N-span   N-chord   ∆s (chord)    average/maximum y+
        coarse           2.32 × 106        171       111     7.81 × 10−6         0.70/2.64
        medium           5.19 × 106        217       139     6.02 × 10−6         0.54/1.55
          fine          10.19 × 106        273       175     4.29 × 10−6         0.39/1.13

    Three levels of O-O topology structured multi-block grids are simulated, with the grid
characteristics presented in Table 4.3, including the number of spanwise and streamwise
nodes on each of the upper and lower surfaces. Figure 4.7 illustrates the skin friction profiles
on the upper and lower surfaces of the Sickle Wing obtained using both the SA-LM2015 and
SA-sLM2015 transition models for all three grid levels overlaid with experimental transition
locations [138]. Both models predict the transition fronts reasonably well.
    As the grid is refined, the transition front on the lower surface, which is dominated by
natural transition, moves slightly downstream and agrees well with the experiment. The
transition front on the upper surface, where crossflow-induced transition is the primary
mechanism, appears more sensitive to the grid refinement. The transition locations at the
kink regions move upstream as the streamwise vortices are better resolved. The lift and drag
grid convergence for both models is presented in Figure 4.8b.
    The pressure coefficient profiles are extracted at the midspan of each sweep section for the
‘fine’ grid, corresponding to spanwise sections at 25%, 55%, and 85% span, and are compared
with the experimental profiles presented by Kruse et al. [90] in Figure 4.8a. The pressure coef-
ficient profiles from the simulations lie on top of one another and the results demonstrate the
‘wiggles’ observed in the experiment at approximately x/c = 0.35 corresponding to natural-
transition on the lower surface of the Sickle Wing are accurately reproduced [90]. There is
a slight deviation between the simulations and the experiment at 25% span on the upper
surface of the Sickle Wing where both transition models predict crossflow-induced transition
upstream of the transition location in the experiment, which is triggered at x/c = 0.80 by a
laminar separation bubble [138, 90]. This discrepancy is also illustrated in Figures 4.7a–4.7b
and is consistent with results obtained using the LM2015 transition model [188].
    Residual and drag coefficient convergence histories are presented in Figure 4.9. The flow
complexity in this case slows convergence for both models; however the smooth variant dis-
plays improved residual convergence behaviour in Figures 4.9a–4.9c. The drag convergence,
the magnitude of the difference between the current drag coefficient and the converged drag
coefficient, plotted in Figures 4.9d–4.9f shows that the smooth variant produces substantially
faster convergence of the drag coefficient.

4.2    Transonic Transition Test Cases

The smooth transition model with and without the compressibility corrections, SA-
sLM2015cc and SA-sLM2015, respectively, is applied to two- and three-dimensional tran-

<!-- ===== PDF page 82 ===== -->

62                                                                             CHAPTER 4. RESULTS: ANALYSIS

               coarse                               medium                                  fine

                                         (a) SA-LM2015 : upper surface

                                         (b) SA-sLM2015 : upper surface

                                         (c) SA-LM2015 : lower surface

                                         (d) SA-sLM2015 : lower surface

Figure 4.7: Skin friction profiles on the upper and lower surfaces of the Sickle Wing for all grid levels overlaid
with experimental transition locations [138]. Results are obtained using the SA-LM2015 and SA-sLM2015
transition models at M = 0.16, α = −2.6◦ , T u = 0.17%, h = 9.78µm, and a Reynolds number of 2.75 × 106 .

<!-- ===== PDF page 83 ===== -->

4.2. TRANSONIC TRANSITION TEST CASES                                                                                                                                                                                            63

                                                                                                                                                                                      SA-LM2015cl
                                                              lower                                                                                                                   SA-LM2015cd
                                                                                                                                                                                      SA-sLM2015cl
                                    -0.4                                                                                             .0064                                            SA-sLM2015cd
                                                                                                               upper
                                                                                                                                                                                                                -0.064
                                    -0.2
                                                                                                                                     .0062
                                        0                       SA-LM201525%

                                                                                                                                  CD
                                  Cp

                                                                                                                                                                                                                  CL
                                                                SA-LM201555%
                                       0.2                      SA-LM201585%
                                                                SA-sLM201525%                                                        .0060
                                                                SA-sLM201555%
                                       0.4                      SA-sLM201585%                                                                                                                                   -0.068
                                                                Experiment25%
                                       0.6                      Experiment55%                                                        .0058
                                                                Experiment85%
                                             0          0.2     0.4         0.6                               0.8       1                       2E-05                                4E-05           6E-05
                                                                      x/c                                                                                                            N-2/3
                                  (a) Pressure coefficient profiles obtained on the                                                      (b) Grid convergence for lift and drag
                                  ‘fine’ grid at three spanwise sections compared                                                                        filler
                                  with the experiment [90]                                                                                               filler

Figure 4.8: Experimental and ‘fine’ grid pressure coefficient profiles, and grid convergence for the SA-LM2015
and SA-sLM2015 transition models. The Sickle Wing is simulated at M = 0.16, α = −2.6◦ , T u = 0.17%,
and a Reynolds number of 2.75 × 106 [90].

                                                 coarse                                                                     medium                                                                     fine
                            101                                                                               101                                                                     101
                                                               SA-LM2015                                                                     SA-LM2015                                                             SA-LM2015
                           10−1                                SA-sLM2015                                    10−1                                                                    10−1
                                                                                                                                             SA-sLM2015                                                            SA-sLM2015
Relative Residual, d

                                                                                  Relative Residual, d

                                                                                                                                                          Relative Residual, d
                           10−3                                                                              10−3                                                                    10−3
                           10−5                                                                              10−5                                                                    10−5
                           10−7                                                                              10−7                                                                    10−7
                           10−9                                                                              10−9                                                                    10−9
                       10−11                                                                             10−11                                                                   10−11
                       10−13                                                                             10−13                                                                   10−13
                       10−150.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0                                          10−150.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0                                10−150.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0
                                  Equivalent Residual Evaluations 1e5                                               Equivalent Residual Evaluations 1e5                                      Equivalent Residual Evaluations 1e5
                                                  (a)                                                                          (b)                                                                      (c)

                           104                                 SA-LM2015                                     104                             SA-LM2015                               104                           SA-LM2015
Drag Convergence(counts)

                                                                                  Drag Convergence(counts)

                                                                                                                                                          Drag Convergence(counts)

                                                               SA-sLM2015                                                                    SA-sLM2015                                                            SA-sLM2015
                           102                                                                               102                                                                     102
                           100                                                                               100                                                                     100
                           10−2                                                                              10−2                                                                    10−2
                           10−4                                                                              10−4                                                                    10−4
                           10−6                                                                              10−6                                                                    10−6
                           10−8                                                                              10−8                                                                    10−8
                      10−100.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0                                          10−100.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0                                10−100.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0
                                  Equivalent Residual Evaluations 1e5                                               Equivalent Residual Evaluations 1e5                                      Equivalent Residual Evaluations 1e5
                                                  (d)                                                                          (e)                                                                      (f)

Figure 4.9: Relative residual and drag coefficient convergence histories for the SA-LM2015 and SA-sLM2015
transition models on all grid levels. The Sickle Wing is simulated at M = 0.16, α = −2.6◦ , T u = 0.17%,
h = 9.78µm, and a Reynolds number of 2.75 × 106 .

sonic transition test cases. These cases consist of the CAST10-2 and VA-2 airfoil test cases
and the transonic NASA CRM-NLF wing-body geometry. The Reynolds number increases
between each test case. As the Reynolds number increases and transition occurs in regions
of more favourable pressure gradient, the stabilizing effect of flow compressibility is expected

<!-- ===== PDF page 84 ===== -->

64                                                                   CHAPTER 4. RESULTS: ANALYSIS

to become more significant (see Section 2.3). For the CRM-NLF simulations, the crossflow
correlations are activated after an initial five order of magnitude residual drop, with the
switch to the inexact-Newton phase delayed until the residual meets the drop tolerance with
the crossflow correlations active.
   The CAST10-2 test case, which is evaluated at a Reynolds number of 2 × 106 , demon-
strates the decreased stabilizing effect of flow compressibility when transition occurs in re-
gions of decelerating flow and therefore the need to accurately model the relationship between
pressure gradient and flow compressibility. The VA-2 airfoil is investigated at a Reynolds
number of 10 × 106 and features transition in regions of both adverse and favourable pres-
sure gradients. The results demonstrate that the TS instability compressibility correction
successfully reproduces the relationship between pressure gradient and the stabilizing effect
of flow compressibility in both flow environments. In addition to transition locations, a
comparison between experimental and simulated skin friction coefficient profiles is provided.
The cases culminate with the NASA CRM-NLF, which is evaluated at a Reynolds number of
approximately 15 × 106 . At this higher Reynolds number, the TS instability compressibility
correction provides a significant improvement to the predicted upper surface transition front
relative to experimental data. Moreover, the stationary crossflow instability compressibil-
ity correction prevents premature transition on the upper surface of the wing and delays
transition on the lower surface where crossflow instabilities are dominant.

4.2.1   CAST10-2 Airfoil

The CAST10-2 airfoil, designed for transonic commercial transport aircraft, features a small
pitching moment and thick trailing edge in order to incorporate a flap system [175]. The
airfoil was investigated by Hebler et al. [66] in the Transonic Wind Tunnel Göttingen (DNW-
TWG) at a Mach number of 0.75, Reynolds number of 2 × 106 , and angles of attack ranging
from −0.79 to 1.41 degrees. The CAST10-2 airfoil has also been investigated numerically
by several researchers [8, 52, 51, 181]. The turbulence intensity of the DNW-TWG wind
tunnel was not reported in the experiment. However, Fehrs et al. [52] suggest a range of
turbulence intensities for the wind tunnel of between 0.25–0.40%, while Ströer et al. [181]
assumed a value of 0.22%. Hebler et al. simulated the CAST10-2 airfoil with the 2D coupled
Euler-Boundary-layer solver MSES [44], which is capable of predicting transition using a
database eN method, and found that an N-factor of 6, which using Mack’s relation [107]
corresponds to a turbulence intensity of approximately 0.24%, produced a transition front
that agreed well with the experiment [66].
   The blunt trailing-edge airfoil coordinates are provided by Dress et al. [45]. As discussed
by Fehrs et al. [52], the trailing edge of the airfoil is located below the reference line, with the
airfoil coordinates exhibiting a 0.88◦ angle of attack. All angles of attack are provided with
respect to this reference system. A Mach number shift of −0.01 and angle-of-attack shift
of −0.30 degrees were identified by Hebler et al. [66] in order to match the experimental

<!-- ===== PDF page 85 ===== -->

4.2. TRANSONIC TRANSITION TEST CASES                                                                                              65

                                  Table 4.4: CAST10-2 structured O-grid dimensions.
                grid level      chord × off-wall nodes            avg/max ∆s × 10−6 (chord)             avg/max y+
                   L0                 541 × 121                          5.06 / 5.52                     0.27 / 0.84
                   L1                 761 × 171                          3.58 / 3.82                     0.19 / 0.59
                   L2                1081 × 241                          2.52 / 2.63                     0.13 / 0.42
                   L3                1521 × 341                          1.79 / 1.86                     0.11 / 0.32

        alpha = -0.39                                                  alpha = 0.82
                                         SA-sLM2015: Cl                                             SA-sLM2015: Cl
                                         SA-sLM2015: Cd                                             SA-sLM2015: Cd
        .0060                            SA-sLM2015cc: Cl
                                                               .4920   .0094                        SA-sLM2015cc: Cl
                                                                                                                          .7120
                                         SA-sLM2015cc: Cd                                           SA-sLM2015cc: Cd

        .0058                                                  .4900   .0092                                              .7100

        .0056                                                  .4880   .0090                                              .7080
        Cd

                                                                       Cd
                                                                 Cl

                                                                                                                            Cl
        .0054                                                  .4860   .0088                                              .7060

        .0052                                                  .4840   .0086                                              .7040

        .0050                                                  .4820   .0084                                              .7020
                0       5E-06      1E-05       1.5E-05      2E-05              0      5E-06   1E-05       1.5E-05      2E-05
                                        -1                                                         -1
                                    N                                                          N
                          (a) −0.39 degrees                                              (b) 0.82 degrees

Figure 4.10: CAST10-2 grid-convergence results produced at −0.39 and 0.82 degrees angle of attack and
M = 0.74, Re = 2 × 106 , and T u = 0.25%.

data with free-air simulations. These shifts are applied in the current work, along with
an assumed turbulence intensity of 0.25% based on the lower bound estimated by Fehrs et
al. [52]. However, the results are plotted and referenced using the uncorrected angles of
attack in order to be consistent with the results presented in [66].
   Four grid levels are investigated using structured multi-block grids, with the grid char-
acteristics presented in Table 4.4. A grid-refinement study is performed using the smooth
transition model both with and without the compressibility corrections at the −0.39 and 0.82
degree angles of attack, with the lift and drag grid convergence illustrated in Figure 4.10.
The largest difference in the Cl and Cd between the second-finest and finest grid levels is less
than 0.30% for both transition model variants at each angle of attack.
   An angle-of-attack sweep is performed on the three finest grid levels both fully turbulent
using the SA turbulence model and with free transition using the smooth transition model
with and without the compressibility corrections, SA-sLM2015cc and SA-sLM2015, respec-
tively. The pressure and upper-surface skin friction coefficient profiles are overlaid with the
pressure profiles from the experiment in Figure 4.11. At low angles of attack, the primary
transition mechanism is a laminar separation bubble. As the angle of attack increases, an
adverse pressure gradient at mid-chord grows, pushing transition significantly upstream.
Above an angle of attack of 0.64 degrees the upper surface pressure coefficient plateaus with
a shock forming at approximately 70% chord. At these higher angles of attack, the weak

<!-- ===== PDF page 86 ===== -->

66                                                                                                    CHAPTER 4. RESULTS: ANALYSIS

                                                                                                                            0.015
                                                                             -1

                                  SAL1
                                  SAL2                                  -0.5
                                                                                                                            0.01
                                  SAL3
                                  SA-sLM2015L1
                                  SA-sLM2015L2                               0

                                                                                                                                Cf,x
                                                                       Cp
                                  SA-sLM2015L3
                                  SA-sLM2015ccL1                                                                            0.005
                                                                            0.5
                                  SA-sLM2015ccL2
                                  SA-sLM2015ccL3
                                  Hebler et al (2013)

                                                                            1.5
                                                                                  0   0.2       0.4         0.6   0.8   1
                                                                                                      x/c

                                                                                            (a) −0.39 degrees

                                                            0.015                                                           0.015
              -1                                                             -1

          -0.5                                                          -0.5
                                                            0.01                                                            0.01

              0                                                              0
                                                                Cf,x

                                                                                                                                Cf,x
        Cp

                                                                       Cp

                                                            0.005                                                           0.005
             0.5                                                            0.5

              1                                                              1
                                                            0                                                               0

             1.5                                                            1.5
                   0   0.2      0.4         0.6   0.8   1                         0   0.2       0.4         0.6   0.8   1
                                      x/c                                                             x/c

                             (b) 0.05 degrees                                               (c) 0.82 degrees

Figure 4.11: CAST10-2 pressure and upper-surface skin friction coefficient profiles produced at M = 0.74,
Re = 2 × 106 , T u = 0.25%, and at three angles of attack overlaid with the pressure profiles from the
experiment [66].

pressure gradients delay transition until the laminar boundary layer separates, which triggers
transition.
   In general, the free-transition pressure coefficient profiles obtained using both transition
model variants demonstrate good agreement with the experiment, especially for the angles
of attack in the linear lift region where the transition locations are relatively insensitive to
changes in the angle of attack. However, for the 0.05 degree angle of attack case, in the
nonlinear lift region, the fully turbulent pressure profiles provide better agreement with the
experiment. The results demonstrate the significant impact of the extent of the laminar
boundary layer on the pressure profiles at these flight conditions.
   The upper-surface skin friction coefficient profiles demonstrate that the two transition
model variants predict the onset of transition at similar streamwise locations, with the transi-
tion location moving slightly downstream with increased grid refinement. The TS instability
compressibility correction does not have a significant effect on the transition onset locations,

<!-- ===== PDF page 87 ===== -->

4.2. TRANSONIC TRANSITION TEST CASES                                                                  67

                                                         SA
                                                         SA-sLM2015
                                   1.4                   SA-sLM2015cc
                                                         Hebler et al (2013)   0.8

                                   1.2

                                    1                                          0.6

                                                                                 x/c
                              Cl
                                   0.8
                                                                               0.4

                                   0.6

                                   0.4                                         0.2

                                   0.2
                                         -1      0       1               2
                                                 alpha
Figure 4.12: CAST10-2 upper-surface transition locations (dashed) and lift curve (solid) produced on the
L3 grid at M = 0.74, Re = 2 × 106 , T u = 0.25%, and over a range of angles of attack compared with the
results from the experiment [66].

as transition primarily occurs either due to a laminar separation bubble or a strong adverse
pressure gradient. The results from Section 2.3 demonstrate that flow compressibility is not
predicted to have a strong stabilizing effect in these flow environments. This behaviour is
also demonstrated in the work by Ströer et al. [181], who simulated the CAST10-2 airfoil
at similar flight conditions using a transport-equation-based formulation of the compressible
AHD criterion.
   The predicted transition locations on the upper surface of the CAST10-2 airfoil and lift
curves, both obtained on the L3 grid, are presented in Figure 4.12 and are compared with
the lift curve and transition locations from Hebler et al. [66]. The results demonstrate
that the smooth transition model both with and without the compressibility corrections
does a reasonable job of predicting the nonlinear transition and lift curves produced by
the experiment. However, at angles of attack between −0.17 and 0.64 degrees, the free-
transition simulations under-predict the drop in the lift curve by over-predicting the extent
of the laminar boundary layer relative to experimental data. This could be a result of the
uncertainty in the turbulence intensity for the wind tunnel at these conditions.
    As illustrated in Figure 4.12, the CAST10-2 airfoil produces highly nonlinear transition
and lift curves, with transition moving forward significantly on the upper surface of the air-
foil over a small range of angles of attack. This feature makes the airfoil an ideal candidate
for investigating the aeroelastic behaviour of NLF designs [66]. However, for the 0.05 degree
angle of attack case, steady-state simulations using the SA-sLM2015cc smooth transition
model produced an oscillatory force history, with the upper surface transition front moving
up and downstream in a periodic manner. The red areas in Figure 4.11 represent the re-
gions between the pressure and skin friction coefficient profiles for the most upstream and
downstream transition locations on each grid level, while the error bars in Figure 4.12 repre-
sent the range in the transition location and lift coefficient over one cycle of this oscillatory

<!-- ===== PDF page 88 ===== -->

68                                                                 CHAPTER 4. RESULTS: ANALYSIS

behaviour.
   The results in Figures 4.11 and 4.12 demonstrate the importance of including the ef-
fects of pressure gradient in the TS instability compressibility correction. Due to the low
Reynolds number, transition occurs in regions of strong adverse pressure gradient, where, as
demonstrated in Section 2.3, stability analysis predicts that flow compressibility produces
a weaker stabilizing effect. In such regions where the flow is decelerating, a compressibility
correction that does not take into account the effects of pressure gradient will over-predict
the stabilizing effect of flow compressibility. This results in the transition model predicting
transition downstream of the second suction peak at near-zero angles of attack, located at
approximately 60% chord, resulting in predicted transition and lift curves that do not cap-
ture the nonlinear trends produced by Hebler et al. [66]. As the Reynolds number increases,
the stabilizing effects of flow compressibility become more significant. This is demonstrated
in the following test cases.

4.2.2   VA-2 Supercritical Airfoil

The VA-2 supercritical airfoil was recently investigated in the DNW-TWG wind tunnel by
Costantini et al. [34]. The airfoil was examined at a Mach number of 0.72, Reynolds number
of 10 × 106 , and angles of attack from −0.4 to 2.0 degrees. The upper-surface skin friction
coefficient distributions and transition locations were determined using a global luminescent
oil-film skin friction field estimation method (GLOFSFE), which measures the development
of the thickness of an oil film, the distribution of which can be used to calculate the skin
friction, based on its luminescent intensity. The pressure distributions used for comparison
in the current work were obtained before applying the oil film.
    Two transition locations were identified in the experiment: the transition onset location,
identified where the laminar skin friction coefficient increased beyond a value of 4 × 10−4 ,
and the end of the transition region, defined as the location where the skin friction coefficient
reached a value of 3 × 10−3 . While the authors state that the end of the transition region
can only be considered as qualitative due to uncertainties in the turbulent skin friction
estimation originating from thick oil-film distributions in these regions, quantitative results
were obtained for transition onset locations with lower uncertainty. Transition locations were
also determined using an automatic temperature-sensitive paint (TSP) transition detection
method, where the recorded transition locations are expected to be between the transition
onset and end locations, with the transition locations agreeing well with those obtained using
the GLOFSFE method [34, 33].
    The Mach number and angle of attack shifts identified by Hebler et al. for the DNW-
TWG wind tunnel [66], −0.01 Mach and −0.30 degrees, respectively, are applied for the
VA-2 simulations. The results are presented and referenced using the uncorrected angles of
attack in order to be consistent with the results presented in [34]. The turbulence intensity
for the VA-2 experimental investigations in the DNW-TWG was not provided. Therefore, the

<!-- ===== PDF page 89 ===== -->

4.2. TRANSONIC TRANSITION TEST CASES                                                                                               69

                                      Table 4.5: VA-2 structured O-grid dimensions.
                 grid level      chord × off-wall nodes            avg/max ∆s × 10−6 (chord)             avg/max y+
                    L0                 541 × 121                          1.01 / 1.07                     0.25 / 0.49
                    L1                 761 × 171                          0.71 / 0.75                     0.18 / 0.35
                    L2                1081 × 241                          0.50 / 0.52                     0.13 / 0.25
                    L3                1521 × 341                          0.36 / 0.37                     0.09 / 0.19

         alpha = -0.40                                                  alpha = 1.80
                                          SA-sLM2015: Cl                                             SA-sLM2015: Cl
                                          SA-sLM2015: Cd                                             SA-sLM2015: Cd
                                          SA-sLM2015cc: Cl
                                                                .3600                                SA-sLM2015cc: Cl
                                                                                                                           .8100
                                          SA-sLM2015cc: Cd                                           SA-sLM2015cc: Cd

         .0070                                                          .0085
                                                                .3500                                                      .8000

         .0065                                                  .3400   .0080                                              .7900
         Cd

                                                                        Cd
                                                                  Cl

                                                                                                                             Cl
                                                                .3300                                                      .7800
         .0060                                                          .0075

                                                                .3200                                                      .7700
         .0055                                                          .0070
                                                                .3100                                                      .7600
                 0       5E-06      1E-05       1.5E-05      2E-05              0      5E-06   1E-05       1.5E-05      2E-05
                                         -1                                                         -1
                                     N                                                          N
                           (a) −0.40 degrees                                              (b) 1.80 degrees

Figure 4.13: VA-2 grid-convergence results produced at −0.40 and 1.80 degrees angle of attack and M = 0.71,
Re = 10 × 106 , and T u = 0.25%.

turbulence intensity is assumed to be consistent with the CAST10-2 DNW-TWG conditions,
with a value of 0.25% used in the current work.
    Four grid levels are investigated using structured multi-block grids with their characteris-
tics presented in Table 4.5. Grid-refinement studies are performed for the smooth transition
model with and without the compressibility corrections at angles of attack of −0.40 and
1.80 degrees, with the results presented in Figure 4.13. The results demonstrate that both
models produce similar grid convergence. The difference in the drag coefficient between the
finest grid level and the grid-converged value was evaluated using Richardson extrapolation
for each transition model variant at both angles of attack. For the −0.40 degree case, the
differences are 0.18% and 0.80% for the smooth transition model with and without the com-
pressibility corrections, respectively, while for the 1.80 degree case the differences are 0.24%
and 0.37%. A more thorough investigation of the grid convergence of the SA-sLM2015cc
smooth transition model for this test case is presented in Appendix B.
   An angle-of-attack sweep is performed on the three finest grid levels using the smooth
transition model both with and without the compressibility corrections. The pressure and
upper-surface skin friction coefficient profiles are overlaid with the profiles from the exper-
iment [34] in Figure 4.14. For the lower angles of attack, −0.40 degrees to 1.20 degrees,
there is good agreement between the simulated pressure coefficient profiles and the profiles
from the experiment. As the angle of attack is increased past 1.20 degrees, the pressure

<!-- ===== PDF page 90 ===== -->

70                                                                 CHAPTER 4. RESULTS: ANALYSIS

profiles deviate from the experimental data. Costantini et al. observed that the algorithm
used to control the DNW-TWG adaptive wind tunnel walls, which are used to reduce the
effect of the walls on the pressure distributions, failed to converge at angles of attack above
1.20 degrees [34]. In their work, they state that the pressure distributions above this angle
of attack are not representative of free-air conditions, which helps to explain the difference
in the simulated and experimental pressure profiles observed at angles of attack of 1.50 and
2.00 degrees. The increased laminar extent of the boundary layer pushes the shock location
aft for the angles of attack of 1.20 and above. Both the GLOFSFE data [34] and TSP
data obtained in a separate experimental campaign by Costantini et al. [32] demonstrate
a separation bubble on the upper surface at approximately 20% chord for the 1.50 degree
angle of attack case. This feature is illustrated by the negative experimental skin friction
coefficient values in Figure 4.14f but is not present in the simulated results, which may be
due to differences between the free-air simulations and the wind tunnel environment.
    The upper-surface skin friction coefficient profiles produced by the smooth transition
model both with and without the compressibility corrections are plotted and compared
with the profiles obtained using the GLOFSFE method [34]. Although Costantini et al.
demonstrate that repeatability for the skin friction coefficient profiles was limited, they
identified skin friction profiles with lower uncertainty [34], which are the profiles plotted in
the current work. The experimental skin friction coefficient profiles are plotted both with
the raw data, which is clipped at a value of 6 × 10−3 , and with the data scaled by a factor of
17. Costantini et al. note that the oil-film thickness in the turbulent boundary-layer regions
was larger than the height of the viscous sublayer, which produced a hydraulically rough
surface that increased the skin friction relative to the clean configuration [34]. This can
explain the difference in turbulent skin friction between the simulations and the experiment,
as well as the increased uncertainty of the experimental locations for the end of the transition
region. The clipped experimental data provides a quantitative comparison of the laminar
skin friction coefficient profiles, which were not as affected by the oil film due to the thinner
oil-film thickness in these regions, and the transition onset locations, while the scaled data
allows for a qualitative comparison of the skin friction coefficient profiles in the transition
and turbulent regions. However, upstream portions of the experimental laminar skin friction
profiles are also clipped for angles of attack above 1.20 degrees, where overestimations were
observed in the experiment due to a build-up of the oil film [34]. The clipped experimental
skin friction coefficient profiles agree well with the simulated profiles in the laminar boundary
layer.
   At the lower angles of attack, below 0.80 degrees, there is a weak favourable pressure
gradient on the upper surface of the airfoil that develops into a weak adverse pressure
gradient as the angle of attack increases. For the more favourable pressure gradients, the
compressibility correction successfully pushes transition aft to better agree with the scaled
GLOFSFE-estimated skin friction coefficient profiles [34]. As the angle of attack increases

<!-- ===== PDF page 91 ===== -->

4.2. TRANSONIC TRANSITION TEST CASES                                                                                                          71

                                                                               -1.5                                                0.015

                                                                                    -1
                               SA-sLM2015L1
                               SA-sLM2015L2                                                                                        0.01
                                                                               -0.5
                               SA-sLM2015L3
                               SA-sLM2015ccL1

                                                                                                                                       Cf,x
                                                                              Cp
                               SA-sLM2015ccL2                                       0
                               SA-sLM2015ccL3
                               Costantini et al (2021)                                                                             0.005
                               Costantini et al (2021)scaled                       0.5
                               Costantini et al (2021)

                                                                                   1.5
                                                                                         0   0.2       0.4         0.6   0.8   1
                                                                                                             x/c

                                                                                                   (a) −0.40 degrees

          -1.5                                                     0.015       -1.5                                                0.015

              -1                                                                    -1

          -0.5                                                     0.01        -0.5                                                0.01
                                                                       Cf,x

                                                                                                                                       Cf,x
        Cp

                                                                              Cp

              0                                                                     0

                                                                   0.005                                                           0.005
             0.5                                                                   0.5

              1                                                                     1

                                                                   0                                                               0
             1.5                                                                   1.5
                   0   0.2      0.4         0.6     0.8        1                         0   0.2       0.4         0.6   0.8   1
                                      x/c                                                                    x/c

                             (b) 0.00 degrees                                                      (c) 0.40 degrees

          -1.5                                                     0.015       -1.5                                                0.015

              -1                                                                    -1

          -0.5                                                     0.01        -0.5                                                0.01
                                                                       Cf,x

                                                                                                                                       Cf,x
        Cp

                                                                              Cp

              0                                                                     0

                                                                   0.005                                                           0.005
             0.5                                                                   0.5

              1                                                                     1

                                                                   0                                                               0
             1.5                                                                   1.5
                   0   0.2      0.4         0.6     0.8        1                         0   0.2       0.4         0.6   0.8   1
                                      x/c                                                                    x/c

                             (d) 0.80 degrees                                                      (e) 1.20 degrees

Figure 4.14: VA-2 pressure and upper-surface skin friction coefficient profiles produced at M = 0.71, Re =
10 × 106 , T u = 0.25%, and over a range of angles of attack overlaid with the results from the experiment [34].

<!-- ===== PDF page 92 ===== -->

72                                                                                                                                CHAPTER 4. RESULTS: ANALYSIS

          -1.5                                                            0.015            -1.5                                                         0.015

              -1                                                                            -1

          -0.5                                                            0.01             -0.5                                                         0.01

                                                                               Cf,x

                                                                                                                                                            Cf,x
        Cp

                                                                                      Cp
              0                                                                                 0

                                                                          0.005                                                                         0.005
             0.5                                                                           0.5

              1                                                                                 1

                                                                          0                                                                             0
             1.5                                                                           1.5
                   0   0.2      0.4         0.6           0.8         1                             0         0.2           0.4         0.6   0.8   1
                                      x/c                                                                                         x/c

                             (f) 1.50 degrees                                                                       (g) 1.80 degrees

                                             -1.5                                                                       0.015

                                              -1

                                             -0.5                                                                       0.01

                                                                                                                             Cf,x
                                        Cp

                                                                                                                        0.005
                                             0.5

                                             1.5
                                                      0         0.2           0.4         0.6           0.8         1
                                                                                    x/c

                                                                      (h) 2.00 degrees

Figure 4.14 (cont.): VA-2 pressure and upper-surface skin friction coefficient profiles produced at M = 0.71,
Re = 10 × 106 , T u = 0.25%, and over a range of angles of attack overlaid with the results from the
experiment [34].

and the flow decelerates on the upper surface, the stabilizing effect of flow compressibility
is reduced, and the transition front moves upstream to match both the smooth transition
model without the compressibility corrections and experimental data. Above an angle of
attack of 1.20 degrees, the upper surface pressure coefficient profiles plateau and the strong
adverse pressure gradient moves aft in the form of a shock wave. The stabilizing effect of
compressibility increases as the adverse pressure gradient is reduced with increasing angle of
attack, pushing transition downstream. Although wind-tunnel effects at the 1.50 and 2.00
degree angles of attack appear to be significant, good agreement is achieved for both the
pressure and skin friction coefficient profiles at 1.80 degrees.
    The predicted transition locations on the L3 grid are compared with values from the
experiments using the GLOFSFE and TSP methods [34] in Figure 4.15. The results demon-
strate that the smooth transition model without the compressibility corrections produces a
flatter transition curve that is upstream of the experimental GLOFSFE-estimated transition

<!-- ===== PDF page 93 ===== -->

4.2. TRANSONIC TRANSITION TEST CASES                                                                                 73

                                           1.5

                                                                   SA-sLM2015
                                            1                      SA-sLM2015cc

                                   alpha
                                                                   Costantini et al (2021)GLOFSFE,L
                                                                   Costantini et al (2021)GLOFSFE,T
                                                                   Costantini et al (2021)TSP
                                           0.5

                                       -0.5

                                                 0.1   0.2   0.3      0.4           0.5
                                                             x/c
Figure 4.15: VA-2 upper-surface transition locations produced on the L3 grid at M = 0.71, Re = 10 × 106 ,
T u = 0.25%, and over a range of angles of attack compared with results from the experiments [34].

onset locations for all angles of attack except −0.40 and 1.20 degrees. However, at −0.40
degrees the transition location produced by the SA-sLM2015 transition model is upstream of
the location produced by the TSP method. The TSP transition locations should lie between
the GLOFSFE-estimated transition onset and end locations, such as for the −0.40, 0.80, and
1.80 degrees angle of attack, demonstrating that there is increased uncertainty at the 0.00
and 0.40 degree angles of attack. This is likely due to the near-zero pressure gradients on
the upper surface at these conditions. The TS instability compressibility correction in the
SA-sLM2015cc transition model delays transition in regions of favourable pressure gradient.
The predicted transition locations agree best with the GLOFSFE-estimated end of transition
region locations, except for 1.80 and 2.00 degrees angle of attack, where they are closer to
the GLOFSE-estimated transition onset locations.

4.2.3     NASA CRM-NLF Wing-Body Geometry

The NASA CRM-NLF configuration was recently investigated as part of the First AIAA
Transition Modelling and Prediction Workshop2 . Transition visualizations, pressure coef-
ficient profiles, and integrated forces were provided from wind-tunnel tests at the NASA
Langley National Transonic Facility (NTF) over a range of angles of attack [104, 105]. The
test conditions for the data provided are presented in Table 4.6. The TS N-factor in the
NTF tunnel was estimated to vary from 4 to 8 [35]. A critical N-factor of 6 was assumed
for the workshop. Using Mack’s relation [107] this corresponds to a turbulence intensity of
approximately 0.24%, which is the value used in the current work. Surface roughness has a
significant effect on crossflow instabilities [148, 38]. To prevent roughness elements destabi-
lizing the boundary layer, the surface of the CRM-NLF wind tunnel model was frequently
sanded and polished, with the average size of the roughness elements after testing measured
   2 https://transitionmodeling.larc.nasa.gov/wp-content/uploads/sites/109/2020/02/TransitionMPW CaseDescriptions.pdf,

accessed June 2021

<!-- ===== PDF page 94 ===== -->

74                                                                                CHAPTER 4. RESULTS: ANALYSIS

                              Table 4.6: CRM-NLF wind tunnel test conditions2 .
              test case    angle of attack (◦ )   Mach number       Reynolds number ×106 (MAC)
                2523            1.44848            0.856489                   14.97197
                2524            1.98031            0.856491                   14.94591
                2525            2.46141            0.856051                   14.90909
                2526            2.93787            0.855801                   14.85308

                             Table 4.7: CRM-NLF structured grid characteristics3 .
          grid level   # of nodes     average/maximum ∆s × 10−6 (chord)            average/maximum y+
             L0        8,893,456                  1.00 / 2.24                           0.36 / 1.50
             L1        16,691,200                 0.78 / 1.75                           0.28 / 1.06
             L2        32,787,200                 0.60 / 1.35                           0.22 / 0.78
             L3        64,330,000                 0.47 / 1.06                           0.18 / 0.64

to vary from 0.83 to 1.10µin [104]. A value of 1.00µin is assumed in the current work.
   Four grid levels are investigated using structured multi-block grids following the gridding
guidelines provided by the workshop committee3 , with the grid characteristics given in Ta-
ble 4.7. A grid-refinement study is performed at the 2524 test conditions (see Table 4.6) using
three modelling strategies: fully turbulent using the SA turbulence model [174, 3] with the
QCR2000 correction [173] in order to accurately capture flow separation at the wing-body
junction, and free transition using the three-equation smooth local correlation-based transi-
tion model both with, SA-QCR2000-sLM2015cc, and without, SA-QCR2000-sLM2015, the
TS and stationary crossflow instability compressibility corrections. The original Dscf cross-
flow instability source term (Equation 2.15) is used for the cases without the compressibility
corrections, with the new Fonset,scf source term presented in Section 2.3 used for the cases with
the compressibility corrections. For both cases, the activation of the crossflow correlations
is delayed until the total residual converges five orders of magnitude.
   The grid-refinement study results at the 2524 test conditions (see Table 4.6) obtained fully
turbulent and with free transition using the smooth transition model both with and without
the compressibility corrections are presented in Figure 4.16. All three models produce similar
grid convergence for lift and drag. The difference in the drag coefficient on the finest grid level
relative to the grid-converged value was calculated for the fully turbulent simulations and
the free-transition simulations with and without the compressibility corrections to be 0.78%,
0.65%, and 0.66%, respectively. A more thorough investigation of the grid convergence of
the SA-QCR2000-sLM2015cc model for this test case is presented in Appendix B.
    Residual convergence histories for the fully turbulent and free-transition simulations are
presented in Figure 4.17. The residuals are plotted against equivalent residual evaluations,
and therefore the additional computational cost associated with increased refinement is not
explicitly illustrated; however, the increased difficulty associated with decreased grid spacing
is evident, especially for the smooth transition model without the compressibility corrections.
  3 https://transitionmodeling.larc.nasa.gov/wp-content/uploads/sites/109/2020/02/CRM-NLF GridGuidelines.pdf, accessed

June 2021

<!-- ===== PDF page 95 ===== -->

4.2. TRANSONIC TRANSITION TEST CASES                                                                                                                                                                                               75

                                        SA-QCR2000                                    .0200
                                        SA-QCR2000-sLM2015
               0.43                     SA-QCR2000-sLM2015cc                                                                                                        -0.08

               0.42                                                                                                                                                 -0.09
                                                                                      .0190
               0.41                                                                                                                                                 -0.10

                                                                                                                                                                    CM
                                                                                        CD
            CL

               0.40                                                                                                                                                 -0.11
                                                                                      .0180
               0.39                                                                                                                                                 -0.12

               0.38                                                                                                                                                 -0.13
                                                                                      .0170
                  0.0E+00          1.0E-05           2.0E-05     3.0E-05                0.0E+00                1.0E-05           2.0E-05     3.0E-05                      0.0E+00          1.0E-05           2.0E-05     3.0E-05
                                             N-2/3                                                                       N-2/3                                                                       N-2/3

                                        (a) Lift                                                                    (b) Drag                                                           (c) Pitching moment

                                Figure 4.16: CRM-NLF grid-convergence results at the 2524 test conditions (α ≈ 2.0◦ ).

                        101                                                                         101                                                                         101
                                                                     L0                                                                          L0                                                                          L0
                        10−1                                         L1                             10−1                                         L1                             10−1                                         L1
Relative Residual, d

                                                                            Relative Residual, d

                                                                                                                                                        Relative Residual, d
                        10−3                                         L2                             10−3                                         L2                             10−3                                         L2
                                                                     L3                                                                          L3                                                                          L3
                        10−5                                                                        10−5                                                                        10−5
                        10−7                                                                        10−7                                                                        10−7
                        10−9                                                                        10−9                                                                        10−9
                   10−11                                                                       10−11                                                                       10−11
                   10−13                                                                       10−13                                                                       10−13
                   10−150.0       0.5        1.0       1.5     2.0    2.5                      10−150.0       0.5        1.0       1.5     2.0    2.5                      10−150.0       0.5        1.0       1.5     2.0    2.5
                               Equivalent Residual Evaluations 1e5                                         Equivalent Residual Evaluations 1e5                                         Equivalent Residual Evaluations 1e5

                                 (a) SA-QCR2000                                                        (b) SA-QCR2000-sLM2015                                                     (c) SA-QCR2000-sLM2015cc

Figure 4.17: CRM-NLF grid-refinement study residual convergence histories at the 2524 test conditions
(α ≈ 2.0◦ ).

Furthermore, the cost for evaluating the free-transition residual is approximately twice that
of the fully turbulent residual. For these simulations, where we also use the original cross-
flow source term, Dscf , the residual spikes where we activate the crossflow correlations are
significant. Smaller increases in the residual are produced by the smooth transition model
at the crossflow activation switch tolerance with the compressibility corrections and the new
Fonset,scf source term. As demonstrated in Figure 4.17, the total residual is converged to
machine-precision for each case. The results demonstrate that the compressibility correc-
tions and new Fonset,scf source term do not significantly affect the numerical behaviour of the
model.
    The pressure and skin friction coefficient profiles for the three finest grids are extracted
at nine spanwise stations across the wing and compared with the experimental profiles in
Figure 4.18. The pressure and skin friction coefficient profiles for the fully turbulent solution
are similar across the grid levels, with the shocks becoming better resolved and the skin
friction increasing as the grid is refined. The upper surface shocks help to differentiate
the upper and lower surface skin friction coefficient profiles, as the shock causes a local
reduction in the upper-surface skin friction coefficient due to the rapid flow deceleration.
The transition locations produced by both transition model variants on the upper surface of

<!-- ===== PDF page 96 ===== -->

76                                                                CHAPTER 4. RESULTS: ANALYSIS

the wing move slightly aft with increased refinement. At the outboard stations, where the
chord Reynolds number is reduced, the smooth transition model with the compressibility
corrections produces a laminar boundary layer up to the shock. As the shock is better
resolved with increased refinement, the transition location moves further aft. For both
transition model variants, the transition locations on the lower surface of the wing, which
are dominated by crossflow instabilities, move upstream as the streamwise vortices are better
resolved, except at the η = 0.910 station where the transition front produced by the smooth
transition model with the compressibility corrections moves downstream with refinement.
    The pressure coefficient profiles produced by the fully turbulent simulations provide bet-
ter agreement with the experiment at the inboard and midspan sections (η < 0.640) relative
to the simulations with free transition. This counter-intuitive behaviour can be explained by
the TSP images in [104, 105], which demonstrate a significant amount of bypass transition
produced by surface imperfections, including the leading-edge pressure ports. The transi-
tion front from the experiment is estimated behind the bypass transition-induced turbulent
wedges. As the simulated transition front moves downstream to better match the estimated
experimental transition front, worse agreement is seen with the experimental pressure coef-
ficient profiles. This behaviour is also demonstrated by Paredes et al. [132], Krimmelbein
and Krumbein [85], and Helm et al. [67]. At the further outboard stations, where there is
less bypass transition observed in the experiment due to the reduced chord Reynolds num-
ber, the increased laminar extent of the boundary layer produced by the smooth transition
model with the compressibility corrections moves the shock locations aft to better agree with
experimental data. The skin friction coefficient profiles demonstrate that the TS instabil-
ity compressibility correction pushes the upper surface transition front, which is dominated
by TS instabilities, downstream in regions of favourable pressure gradient, where flow com-
pressibility is predicted to have a more significant stabilizing effect (see Section 2.3). The
stationary crossflow instability compressibility correction moves the lower surface transition
front, which is dominated by crossflow instabilities, downstream and prevents an inaccurate,
upstream transition front on the upper surface.
   To better compare the effects of the compressibility corrections on the predicted transition
fronts, the skin friction coefficient distributions on the upper surface of the wing produced
by the smooth transition model both with and without the compressibility corrections are
overlaid with the estimated transition fronts provided by [104] in Figure 4.19. The results
again demonstrate that the transition fronts are largely insensitive to increased grid refine-
ment, with the transition fronts produced by both transition model variants moving slightly
downstream on the upper surface of the wing with increased grid refinement. However, the
transition fronts produced by the smooth transition model without the compressibility cor-
rections remain significantly upstream of the experiment. The TS instability compressibility
correction successfully delays transition in regions of favourable pressure gradient to better
agree with the experiment, particularly in the outboard regions. The wing lower surface and

<!-- ===== PDF page 97 ===== -->

4.2. TRANSONIC TRANSITION TEST CASES                                                                                                                                              77

                                                                                  SA-QCR2000L1
                                                                                  SA-QCR2000L2
                                                                                  SA-QCR2000L3
                                                                                  SA-QCR2000-sLM2015L1
                                                                                  SA-QCR2000-sLM2015L2
                                                                                  SA-QCR2000-sLM2015L3
                                                                                  SA-QCR2000-sLM2015ccL1
                                                                                  SA-QCR2000-sLM2015ccL2
                                                                                  SA-QCR2000-sLM2015ccL3
                                                                                  Lynde et al (2019)

                                                   0.009                                                        0.009                                                        0.009

  -0.5                                                        -0.5                                                         -0.5

                                                   0.006                                                        0.006                                                        0.006
      0                                                            0                                                            0
                                                       Cf,x

                                                                                                                    Cf,x

                                                                                                                                                                                 Cf,x
 Cp

                                                              Cp

                                                                                                                           Cp
  0.5                                                          0.5                                                          0.5
                                                   0.003                                                        0.003                                                        0.003

      1                                                            1                                                            1

                                                   0                                                            0                                                            0
          0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1
                             x/c                                                          x/c                                                          x/c

                    (a) η = 0.163                                                (b) η = 0.252                                                (c) η = 0.370

                                                   0.009                                                        0.009                                                        0.009

  -0.5                                                        -0.5                                                         -0.5

                                                   0.006                                                        0.006                                                        0.006
      0                                                            0                                                            0
                                                       Cf,x

                                                                                                                    Cf,x

                                                                                                                                                                                 Cf,x
 Cp

                                                              Cp

  0.5                                                          0.5                                                         Cp
                                                                                                                            0.5
                                                   0.003                                                        0.003                                                        0.003

      1                                                            1                                                            1

                                                   0                                                            0                                                            0
          0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1
                             x/c                                                          x/c                                                          x/c

                    (d) η = 0.460                                                (e) η = 0.550                                                (f) η = 0.640

                                                   0.009                                                        0.009                                                        0.009

  -0.5                                                        -0.5                                                         -0.5

                                                   0.006                                                        0.006                                                        0.006
      0                                                            0                                                            0
                                                       Cf,x

                                                                                                                    Cf,x

                                                                                                                                                                                 Cf,x
 Cp

                                                              Cp

                                                                                                                           Cp

  0.5                                                          0.5                                                          0.5
                                                   0.003                                                        0.003                                                        0.003

      1                                                            1                                                            1

                                                   0                                                            0                                                            0
          0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1
                             x/c                                                          x/c                                                          x/c

                    (g) η = 0.730                                                (h) η = 0.820                                                (i) η = 0.910

Figure 4.18: CRM-NLF grid-refinement study pressure and skin friction coefficient profiles produced at the
2524 test conditions (α ≈ 2.0◦ ) compared with the pressure profiles from the experiment [104] at varying
spanwise stations, η.

<!-- ===== PDF page 98 ===== -->

78                                                                                CHAPTER 4. RESULTS: ANALYSIS

fuselage nose boundary layers were tripped in the experiment. Therefore, transition fronts
on the fuselage and lower surface of the wing are not available for comparison.
   An angle-of-attack sweep is performed across the test conditions presented in Table 4.6 on
the L1 grid fully turbulent and with free transition both with and without the compressibility
corrections. The skin friction coefficient profiles produced by the free-transition simulations
are presented in Figure 4.20. Again, the transition fronts produced by the smooth transition
model without the compressibility corrections are consistently upstream of the experiment,
especially in regions of favourable pressure gradient. The compressibility corrections delay
transition, producing transition fronts that better agree with experimental data. Regions
of the transition fronts produced by the smooth transition model with the compressibility
corrections are still upstream of the experiment. This is particularly evident for the 2523
test case, where transition is upstream of experimental data inboard near the Yehudi break.
However, this behaviour is also present to a lesser extent in the work of Paredes et al. [132]
and Krimmelbein and Krumbein [85], who both investigate the CRM-NLF configuration
using a variety of stability analysis methods. These differences could be due to uncertainties
in the turbulence intensity or flow conditions. In general, there is better agreement between
the smooth transition model with the compressibility corrections and experimental data
outboard.
   The angle-of-attack sweep curves for lift, drag, and pitching moment coefficient produced
by the fully turbulent simulations and the free-transition simulations both with and without
the compressibility corrections on the L1 grid are compared with experimental data [104] in
Figure 4.21. All three modelling strategies produce increased lift relative to the experiment,
with higher lift, drag, and pitching moment coefficient slopes. The free-transition simula-
tions produce increased lift and decreased drag and pitching moment relative to both the
fully turbulent simulations and the experimental data as the laminar extent of the boundary
layer increases. These trends are consistent with the majority of results presented at the
First AIAA Transition Prediction and Modeling Workshop4 [67]. As previously discussed,
the transition fronts from the experiment demonstrate significant regions of bypass transi-
tion [104]. Furthermore, the wing lower surface and fuselage nose were explicitly tripped.
Therefore, there is more laminar flow in the free-transition simulations than the experiment,
which, combined with possible wind tunnel blockage effects, could help to explain the differ-
ence in the simulated and experimental results. Although explicitly tripping the flow could
improve the agreement between the simulated and experimental integrated loads, simula-
tions are performed in the current work with free transition globally in order to investigate
the impact of the compressibility corrections across the geometry, such as on the stationary
crossflow instabilities on the wing lower surface.
   The angle-of-attack sweep residual convergence histories for the fully turbulent simula-
tions and the free-transition simulations are presented in Figure 4.22. The total residual is
     4 https://transitionmodeling.larc.nasa.gov/workshop i/, accessed July 2021

<!-- ===== PDF page 99 ===== -->

4.2. TRANSONIC TRANSITION TEST CASES                                                                      79

               (a) L1 : SA-QCR2000-sLM2015                   (b) L1 : SA-QCR2000-sLM2015cc

               (c) L2 : SA-QCR2000-sLM2015                   (d) L2 : SA-QCR2000-sLM2015cc

               (e) L3 : SA-QCR2000-sLM2015                   (f) L3 : SA-QCR2000-sLM2015cc

Figure 4.19: CRM-NLF grid-refinement study upper-surface skin friction coefficient profiles at the 2524 test
conditions (α ≈ 2.0◦ ) overlaid with the estimated transition front from the experiment [104].

<!-- ===== PDF page 100 ===== -->

80                                                                        CHAPTER 4. RESULTS: ANALYSIS

         (a) 2523 (α ≈ 1.5◦ ): SA-QCR2000-sLM2015      (b) 2523 (α ≈ 1.5◦ ): SA-QCR2000-sLM2015cc

         (c) 2525 (α ≈ 2.5◦ ): SA-QCR2000-sLM2015      (d) 2525 (α ≈ 2.5◦ ): SA-QCR2000-sLM2015cc

         (e) 2526 (α ≈ 3.0◦ ): SA-QCR2000-sLM2015      (f) 2526 (α ≈ 3.0◦ ): SA-QCR2000-sLM2015cc

Figure 4.20: CRM-NLF angle-of-attack sweep upper-surface skin friction coefficient profiles obtained on the
L1 grid overlaid with the estimated transition front from the experiment [104].

<!-- ===== PDF page 101 ===== -->

4.2. TRANSONIC TRANSITION TEST CASES                                                                                                                                                                                                          81

                                       SA-QCR2000                                           0.028                                                                             -0.07
               0.65                    SA-QCR2000-sLM2015
                                       SA-QCR2000-sLM2015cc
                                       Lynde et al (2019)                                   0.026                                                                             -0.08
                    0.6

               0.55                                                                         0.024                                                                             -0.09

                    0.5                                                                     0.022                                                                                  -0.1

                                                                                                                                                                        CM
                                                                                        CD
      CL

               0.45                                                                              0.02                                                                         -0.11

                    0.4                                                                     0.018                                                                             -0.12

               0.35
                                                                                            0.016                                                                             -0.13
                    0.3
                                                                                            0.014                                                                             -0.14
                         1.2 1.4 1.6 1.8   2   2.2 2.4 2.6 2.8    3   3.2                       1.2 1.4 1.6 1.8            2   2.2 2.4 2.6 2.8    3   3.2                         1.2 1.4 1.6 1.8          2   2.2 2.4 2.6 2.8    3   3.2
                                      angle of attack                                                                 angle of attack                                                                 angle of attack

                                           (a) Lift                                                                     (b) Drag                                                                 (c) Pitching moment

Figure 4.21: CRM-NLF angle-of-attack sweep force curves for lift, drag, and pitching moment coefficient
produced on the L1 grid compared with results from the experiment [104].

                        101                                                                               101                                                                             101
                                                                      2523                                                                            2523                                                                            2523
                        10−1                                          2525                                10−1                                        2525                                10−1                                        2525
Relative Residual, d

                                                                                  Relative Residual, d

                                                                                                                                                                  Relative Residual, d
                        10−3                                          2526                                10−3                                        2526                                10−3                                        2526
                        10−5                                                                              10−5                                                                            10−5
                        10−7                                                                              10−7                                                                            10−7
                        10−9                                                                              10−9                                                                            10−9
                   10−11                                                                             10−11                                                                           10−11
                   10−13                                                                             10−13                                                                           10−13
                   10−150.0          0.5         1.0        1.5             2.0                      10−150.0        0.5         1.0        1.5             2.0                      10−150.0        0.5         1.0        1.5             2.0
                               Equivalent Residual Evaluations 1e5                                               Equivalent Residual Evaluations 1e5                                             Equivalent Residual Evaluations 1e5

                                  (a) SA-QCR2000                                                             (b) SA-QCR2000-sLM2015                                                         (c) SA-QCR2000-sLM2015cc

                  Figure 4.22: CRM-NLF angle-of-attack sweep residual convergence histories produced on the L1 grid.

converged to machine-zero for all three modelling strategies at each angle of attack. Again,
the smooth transition model without the compressibility corrections produces larger residual
spikes at the crossflow activation switch relative to the smooth model with the corrections
and the new Fonset,scf source term. In general, the smooth transition model with and without
the compressibility corrections produces similar numerical behaviour.
   The pressure and skin friction coefficient profiles produced by the fully turbulent simu-
lations and the free-transition simulations both with and without the compressibility cor-
rections at the 2523, 2525, and 2526 test conditions (see Table 4.6) are extracted at nine
spanwise locations on the L1 grid and illustrated in Figures 4.23, 4.24, and 4.25, respec-
tively. Similar to the 2524 results presented in Figure 4.18, there is a consistent trend with
the compressibility corrections delaying crossflow transition on the lower surface, and signifi-
cantly delaying TS transition on the upper surface in regions of favourable pressure gradient.
Unfortunately, there is again significant bypass transition observed in the experiment at the
2523, 2525, and 2526 test conditions, especially at the leading-edge pressure ports [104, 105],
which makes a detailed comparison of the pressure coefficient profiles difficult.

<!-- ===== PDF page 102 ===== -->

82                                                                                                                            CHAPTER 4. RESULTS: ANALYSIS

                                                                                    SA-QCR2000L1
                                                                                    SA-QCR2000-sLM2015L1
                                                                                    SA-QCR2000-sLM2015ccL1
                                                                                    Lynde et al (2019)

                                                     0.009                                                        0.009                                                        0.009

     -0.5                                                       -0.5                                                         -0.5

                                                     0.006                                                        0.006                                                        0.006
       0                                                             0                                                            0
                                                         Cf,x

                                                                                                                      Cf,x

                                                                                                                                                                                   Cf,x
 Cp

                                                                Cp

                                                                                                                             Cp
     0.5                                                         0.5                                                          0.5
                                                     0.003                                                        0.003                                                        0.003

       1                                                             1                                                            1

                                                     0                                                            0                                                            0
            0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1
                               x/c                                                          x/c                                                          x/c

                      (a) η = 0.163                                                (b) η = 0.252                                                (c) η = 0.370

                                                     0.009                                                        0.009                                                        0.009

     -0.5                                                       -0.5                                                         -0.5

                                                     0.006                                                        0.006                                                        0.006
       0                                                             0                                                            0
                                                         Cf,x

                                                                                                                      Cf,x

                                                                                                                                                                                   Cf,x
 Cp

                                                                Cp

                                                                                                                             Cp

     0.5                                                         0.5                                                          0.5
                                                     0.003                                                        0.003                                                        0.003

       1                                                             1                                                            1

                                                     0                                                            0                                                            0
            0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1
                               x/c                                                          x/c                                                          x/c

                      (d) η = 0.460                                                (e) η = 0.550                                                (f) η = 0.640

                                                     0.009                                                        0.009                                                        0.009

     -0.5                                                       -0.5                                                         -0.5

                                                     0.006                                                        0.006                                                        0.006
       0                                                             0                                                            0
                                                         Cf,x

                                                                                                                      Cf,x

                                                                                                                                                                                   Cf,x
 Cp

                                                                Cp

                                                                                                                             Cp

     0.5                                                         0.5                                                          0.5
                                                     0.003                                                        0.003                                                        0.003

       1                                                             1                                                            1

                                                     0                                                            0                                                            0
            0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1
                               x/c                                                          x/c                                                          x/c

                      (g) η = 0.730                                                (h) η = 0.820                                                (i) η = 0.910

Figure 4.23: CRM-NLF pressure and skin friction coefficient profiles produced at the 2523 test conditions
(α ≈ 1.5◦ ) on the L1 grid compared with the pressure profiles from the experiment [104] at varying spanwise
stations, η.

<!-- ===== PDF page 103 ===== -->

4.2. TRANSONIC TRANSITION TEST CASES                                                                                                                                              83

                                                                                  SA-QCR2000L1
                                                                                  SA-QCR2000-sLM2015L1
                                                                                  SA-QCR2000-sLM2015ccL1
                                                                                  Lynde et al (2019)

                                                   0.009                                                        0.009                                                        0.009

  -0.5                                                        -0.5                                                         -0.5

                                                   0.006                                                        0.006                                                        0.006
      0                                                            0                                                            0

                                                       Cf,x

                                                                                                                    Cf,x

                                                                                                                                                                                 Cf,x
 Cp

                                                              Cp

                                                                                                                           Cp
   0.5                                                         0.5                                                          0.5
                                                   0.003                                                        0.003                                                        0.003

      1                                                            1                                                            1

                                                   0                                                            0                                                            0
          0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1
                             x/c                                                          x/c                                                          x/c

                    (a) η = 0.163                                                (b) η = 0.252                                                (c) η = 0.370

                                                   0.009                                                        0.009                                                        0.009

  -0.5                                                        -0.5                                                         -0.5

                                                   0.006                                                        0.006                                                        0.006
      0                                                            0                                                            0
                                                       Cf,x

                                                                                                                    Cf,x

                                                                                                                                                                                 Cf,x
 Cp

                                                              Cp

                                                                                                                           Cp
   0.5                                                         0.5                                                          0.5
                                                   0.003                                                        0.003                                                        0.003

      1                                                            1                                                            1

                                                   0                                                            0                                                            0
          0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1
                             x/c                                                          x/c                                                          x/c

                    (d) η = 0.460                                                (e) η = 0.550                                                (f) η = 0.640

                                                   0.009                                                        0.009                                                        0.009

  -0.5                                                        -0.5                                                         -0.5

                                                   0.006                                                        0.006                                                        0.006
      0                                                            0                                                            0
                                                       Cf,x

                                                                                                                    Cf,x

                                                                                                                                                                                 Cf,x
 Cp

                                                              Cp

                                                                                                                           Cp

   0.5                                                         0.5                                                          0.5
                                                   0.003                                                        0.003                                                        0.003

      1                                                            1                                                            1

                                                   0                                                            0                                                            0
          0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1
                             x/c                                                          x/c                                                          x/c

                    (g) η = 0.730                                                (h) η = 0.820                                                (i) η = 0.910

Figure 4.24: CRM-NLF pressure and skin friction coefficient profiles produced at the 2525 test conditions
(α ≈ 2.5◦ ) on the L1 grid compared with the pressure profiles from the experiment [104] at varying spanwise
stations, η.

<!-- ===== PDF page 104 ===== -->

84                                                                                                                            CHAPTER 4. RESULTS: ANALYSIS

                                                                                    SA-QCR2000L1
                                                                                    SA-QCR2000-sLM2015L1
                                                                                    SA-QCR2000-sLM2015ccL1
                                                                                    Lynde et al (2019)

      -1                                             0.009        -1                                              0.009        -1                                              0.009

     -0.5                                                       -0.5                                                         -0.5

                                                     0.006                                                        0.006                                                        0.006
       0                                                             0                                                            0
                                                         Cf,x

                                                                                                                      Cf,x

                                                                                                                                                                                   Cf,x
 Cp

                                                                Cp

                                                                                                                             Cp
     0.5                                                         0.5                                                          0.5
                                                     0.003                                                        0.003                                                        0.003

       1                                                             1                                                            1

                                                     0                                                            0                                                            0
            0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1
                               x/c                                                          x/c                                                          x/c

                      (a) η = 0.163                                                (b) η = 0.252                                                (c) η = 0.370

      -1                                             0.009        -1                                              0.009        -1                                              0.009

     -0.5                                                       -0.5                                                         -0.5

                                                     0.006                                                        0.006                                                        0.006
       0                                                             0                                                            0
                                                         Cf,x

                                                                                                                      Cf,x

                                                                                                                                                                                   Cf,x
 Cp

                                                                Cp

                                                                                                                             Cp

     0.5                                                         0.5                                                          0.5
                                                     0.003                                                        0.003                                                        0.003

       1                                                             1                                                            1

                                                     0                                                            0                                                            0
            0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1
                               x/c                                                          x/c                                                          x/c

                      (d) η = 0.460                                                (e) η = 0.550                                                (f) η = 0.640

      -1                                             0.009        -1                                              0.009        -1                                              0.009

     -0.5                                                       -0.5                                                         -0.5

                                                     0.006                                                        0.006                                                        0.006
       0                                                             0                                                            0
                                                         Cf,x

                                                                                                                      Cf,x

                                                                                                                                                                                   Cf,x
 Cp

                                                                Cp

                                                                                                                             Cp

     0.5                                                         0.5                                                          0.5
                                                     0.003                                                        0.003                                                        0.003

       1                                                             1                                                            1

                                                     0                                                            0                                                            0
            0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1                       0   0.2      0.4         0.6   0.8   1
                               x/c                                                          x/c                                                          x/c

                      (g) η = 0.730                                                (h) η = 0.820                                                (i) η = 0.910

Figure 4.25: CRM-NLF pressure and skin friction coefficient profiles produced at the 2526 test conditions
(α ≈ 3.0◦ ) on the L1 grid compared with the pressure profiles from the experiment [104] at varying spanwise
stations, η.

<!-- ===== PDF page 105 ===== -->

Chapter 5

Results: Optimization

In this chapter, the free-transition optimization framework presented in Chapter 3 is applied
to single-point lift-constrained drag minimizations of airfoils and an infinite swept wing
geometry. Specifically, the smooth transition model with compressibility corrections, SA-
sLM2015cc, is used for the results presented. Airfoil optimizations are performed using an
initial geometry based on a blunt trailing-edge RAE2822 airfoil. For the infinite swept wing
optimization, the airfoil is extruded one chord length with 11 nodes in the spanwise direction
and with periodic boundary conditions applied at the wing root and tip. The airfoil and
wing surfaces are parameterized using B-spline surface patches, which are controlled using
FFD design variables. For the results presented, 6 streamwise FFD design variables are
applied to each of the upper and lower surfaces, with the infinite swept wing root and tip
design variables constrained to be equal and the leading- and trailing-edge design variables
constrained to move symmetrically, for a total of 10 effective geometric design variables plus
angle of attack.

5.1    Airfoil Optimization

For the two-dimensional airfoil optimizations, three design conditions are considered that are
representative of aircraft classes ranging from light aircraft up to single-aisle aircraft based
on the work of Rashad and Zingg [150], as presented in Table 5.1. The Cessna 172R and
De Havilland Dash8-Q400 aircraft have nominally zero sweep; however, the Boeing 737-800
design Mach number is corrected from a cruise Mach number of 0.785 to an effective Mach
number of 0.71 based on a wing sweep of 25 degrees. Although the two-dimensional Boeing
737-800 design condition does not include the effects of crossflow instabilities, and therefore
is not representative of the disturbance environment for a single-aisle aircraft wing, drag
minimizations at these design conditions are valuable for investigating the streamwise grid
resolution requirements for NLF optimizations at large Mach and Reynolds numbers, as well
as for comparing the optimized designs with designs previously developed by Rashad and
Zingg [150]. A drag minimization of an infinite swept wing that includes the effects of both

<!-- ===== PDF page 106 ===== -->

86                                                                        CHAPTER 5. RESULTS: OPTIMIZATION

Table 5.1: Design conditions for the two-dimensional lift-constrained drag minimizations [150]. For each case
turbulence intensity, T u, is specified as 0.07%.
              design condition         lift coefficient   Mach number      Reynolds number (MAC)
               Cessna 172R                   0.30              0.19                5.6 × 106
          De Havilland Dash8-Q400            0.42              0.60               15.7 × 106
              Boeing 737-800                 0.50          0.71 (corr.)           20.3 × 106

Table 5.2: Structured multi-block O-grid dimensions for the two-dimensional optimizations at the Cessna
172R design conditions (see Table 5.1).
             grid level   chord x off-wall nodes    avg/max ∆s × 10−6 (chord)      avg/max y+
                L0               281x121                   1.09 / 1.19              0.17 / 0.54
                L1               561x121                        "                        "
                L2              1121x121                        "                        "

TS and stationary crossflow instabilities is investigated in Section 5.2.

5.1.1    Cessna 172R Skyhawk

Preliminary optimizations revealed that the capability of the optimization algorithm to delay
the onset of boundary-layer transition when applied to free-transition design problems was
sensitive to the streamwise resolution of the grid. To investigate this further, three structured
multi-block O-type topology grids with varying levels of streamwise grid resolution were con-
sidered for the Cessna 172R drag minimizations, with the grid details presented in Table 5.2.
Lift-constrained drag minimizations are performed at each grid level with the cross-sectional
area of the airfoil constrained so that it cannot decrease and with loose thickness-to-chord
ratio constraints enforced at each FFD design variable pair. The optimization problem,
which will be referred to as the baseline optimization problem, is given by,

                                        min Cd (X ),                                                   (5.1)
                                         X

                                        s.t. Cl = Cl∗ ,                                                (5.2)
                                                A ≥ Ainit ,                                            (5.3)
                                                t/c ≥ 0.15(t/c)init ,                                  (5.4)

where X represents the vector of design variables, Cd and Cl are the drag and lift coefficients,
respectively, Cl∗ represents the design lift coefficient presented in Table 5.1, A represents the
cross-sectional area, t/c represents the thickness-to-chord ratio evaluated at each FFD design
variable pair, and the subscript ‘init’ denotes a quantity associated with the initial design.
   The results from the baseline drag minimizations are presented in Figure 5.1. The re-
sults suggest that the large gradients in the transition region can be under-resolved due to
poor streamwise grid resolution. Instead of delaying transition, the results in Figure 5.1a
demonstrate that for the L0 grid level the optimization algorithm generates significant ad-
verse pressure gradients downstream of the transition locations by reducing the trailing-edge

<!-- ===== PDF page 107 ===== -->

5.1. AIRFOIL OPTIMIZATION                                                                                                                                                                                                                                                                                                                           87

                                                   Optimality                                                                                    2
                                                                                                                                                                    Optimality                                                                                       2
                                                                                                                                                                                                                                                                                         Optimality
                           10                      Feasibility                                            .0060                             10                      Feasibility                                               .0060                             10                       Feasibility                                             .0060

                                                                                     Merit function

                                                                                                                                                                                                         Merit function

                                                                                                                                                                                                                                                                                                                            Merit function
                                                   Merit function                                                                                                   Merit function                                                                                                       Merit function
                                                                                                          .0055                                                                                                               .0055                                                                                                              .0055

 Optimality, Feasibility

                                                                                                                  Optimality, Feasibility

                                                                                                                                                                                                                                      Optimality, Feasibility
                           100                                                                                                              100                                                                                                                 100
                                                                                                          .0050                                                                                                               .0050                                                                                                              .0050

                                -2                                                                        .0045                             10
                                                                                                                                                 -2                                                                           .0045                             10
                                                                                                                                                                                                                                                                     -2                                                                          .0045

                                                                                                          .0040                                                                                                               .0040                                                                                                              .0040
                           10-4                                                                           .0035                             10-4                                                                              .0035                             10-4                                                                             .0035

                                -6
                                                                                                          .0030                                  -6
                                                                                                                                                                                                                              .0030                                  -6
                                                                                                                                                                                                                                                                                                                                                 .0030
                           10                                                                                                               10                                                                                                                  10
                                                                                                          .0025                                                                                                               .0025                                                                                                              .0025

                           10-8                                                                           .0020                             10-8                                                                              .0020                             10-8                                                                             .0020
                                     0         5                10              15                                                                    0   20       40         60          80      100                                                                     0   20       40          60     80      100       120
                                               Design iteration                                                                                                   Design iteration                                                                                                     Design iteration

                                                                         Optimized                                                                                                             Optimized                                                                                                          Optimized
                                                                         Initial                          .015                                                                                 Initial                        .015                                                                                Initial                        .015

                       -0.6                                                                                                             -0.6                                                                                                                -0.6
                                                                                                          .012                                                                                                                .012                                                                                                               .012

                           0.0                                                                                                              0.0                                                                                                                 0.0
                                                                                                          .009                                                                                                                .009                                                                                                               .009
 Cp

                                                                                                                  Cp

                                                                                                                                                                                                                                      Cp
                           0.6                                                                                                              0.6                                                                                                                 0.6
                                                                                                            Cf

                                                                                                                                                                                                                                Cf

                                                                                                                                                                                                                                                                                                                                                   Cf
                                                                                                          .006                                                                                                                .006                                                                                                               .006
                           1.2                                                                                                              1.2                                                                                                                 1.2
                                          upper                                                                                                                 upper                                                                                                                upper
                                                                                                          .003                                                                                                                .003                                                                                                               .003
                           1.8                       upper                                                                                  1.8                           lower                                                                                 1.8                            lower
                                                      lower                                                                                                                                      lower                                                                                                              lower
                                                                 lower                                                                                                               upper                                                                                                                upper
                           2.4                                                                            .000                              2.4                                                                               .000                              2.4                                                                              .000
                                     0   0.2        0.4         0.6       0.8                         1                                               0   0.2           0.4         0.6         0.8                       1                                               0    0.2           0.4         0.6       0.8                       1
                                                          x/c                                                                                                                 x/c                                                                                                                  x/c

                                                     (a) L0                                                                                                             (b) L1                                                                                                               (c) L2

Figure 5.1: Optimization convergence histories and cross-sectional profiles of the initial and optimized designs
produced by drag minimizations with the baseline optimization problem (Equations 5.1–5.4) at the Cessna
172R conditions with varying streamwise grid resolution.

thickness, where the flow deceleration reduces the turbulent skin friction coefficient and skin
friction drag. This behaviour is also produced in the Dash8-Q400 optimizations presented in
the following section (Figure 5.3a). The optimization algorithm produces more efficient de-
signs with larger extents of laminar flow as the streamwise resolution of the grid is increased.
Figure 5.1b illustrates that at the L1 grid level the optimization algorithm is able to delay
the onset of boundary-layer transition on both surfaces significantly. Diminishing returns are
encountered by doubling the streamwise resolution going from the L1 to the L2 grid level,
which provides a modest reduction in drag but significantly increases the computational cost
of the optimizations.
   The optimized designs in Figures 5.1b and 5.1c feature concave pressure recovery re-
gions on the airfoil upper surface and extended regions of favourable pressure gradient on
the lower surface. Instead of a flat pressure plateau on the upper surface, the optimiza-
tion algorithm has produced a flow deceleration section followed by flow acceleration. The
optimization algorithm has achieved these designs by increasing the angle of attack, decreas-
ing the trailing-edge thickness, and moving the aerodynamic loading upstream. However,
from a manufacturability standpoint the thin trailing edge is undesirable. Furthermore, the
design features a reflexed trailing edge that reduces lift in order to satisfy the target lift
coefficient, and steep adverse pressure gradients that could reduce the performance of the

<!-- ===== PDF page 108 ===== -->

88                                                                   CHAPTER 5. RESULTS: OPTIMIZATION

Table 5.3: Aerodynamic performance of the initial and optimized designs at the Cessna 172R conditions on
the L1 grid level.
                                                     Cd (cnts.)    Cl      Cm       L/D     aoa (deg.)
                         Initial                       52.04      0.30   −0.060    57.647     0.971
             Baseline Optimization Problem             29.04      0.30    0.034   103.294     2.542
                 Additional Constraints                28.88      0.30   −0.060   103.872     0.551
     Re-optimized w/ Baseline Optimization Problem     27.75      0.30   −0.058   108.127     0.600

airfoil at higher angles of attack. Rashad and Zingg [150] produced a more optimal design
at the same design conditions using a discrete-adjoint optimization algorithm coupled with
a stability-analysis framework that features more aft-loading with increased trailing-edge
thickness. To investigate whether the optimization algorithm can recover a geometry similar
to that of Rashad and Zingg [150] and to address the undesirable characteristics mentioned
above, the optimization problem was modified to include a pitching-moment constraint and
more conservative thickness-to-chord ratio constraints. The optimization problem with the
additional constraints is presented below,

                                      min Cd (X ),                                                 (5.5)
                                       X

                                       s.t. Cl = Cl∗ ,                                             (5.6)
                                             Cm = Cm,init ,                                        (5.7)
                                             A ≥ Ainit ,                                           (5.8)
                                             t/c ≥ 0.85(t/c)init ,                                 (5.9)

where Cm is the pitching-moment coefficient.
   Based on the results from the streamwise grid-resolution study presented in Figure 5.1, the
drag minimization with the optimization problem with additional constraints is performed
on the L1 grid level. The results from the drag minimizations with the baseline and more
constrained optimization problems are compared in Figure 5.2. The optimization with the
more constrained optimization problem (Equations 5.5–5.9) produces a design with extended
regions of near-zero pressure gradient on both the upper and lower surfaces of the airfoil that
delay boundary-layer transition and increase the laminar extent of the boundary layer. The
optimization algorithm produces this design by maintaining more aft-loading relative to the
baseline optimized design, with the new design also featuring increased trailing-edge thick-
ness. The resulting design closely resembles the design produced by the lift-constrained drag
minimization performed by Rashad and Zingg [150]. Moreover, the optimization history
demonstrates that the design produced by the optimization with additional constraints pro-
duces lower drag than the design produced by the baseline optimization problem, despite
lying in the design space of the latter. The results suggest that the Cessna 172R design
space is multi-modal. Details of the aerodynamic performance of the initial blunt trailing-
edge RAE2822 airfoil and the two optimized designs are presented in Table 5.3.

<!-- ===== PDF page 109 ===== -->

5.1. AIRFOIL OPTIMIZATION                                                                                                                                                                                                                                                                                                                             89

                                                   Optimality                                                                                       2
                                                                                                                                                                           Optimality                                                                                       2
                                                                                                                                                                                                                                                                                            Optimality
                           10                      Feasibility                                               .0060                             10                          Feasibility                                               .0060                             10                   Feasibility                                            .0060

                                                                                        Merit function

                                                                                                                                                                                                                Merit function

                                                                                                                                                                                                                                                                                                                              Merit function
                                                   Merit function                                                                                                          Merit function                                                                                                   Merit function
                                                                                                             .0055                                                                                                                   .0055                                                                                                         .0055

 Optimality, Feasibility

                                                                                                                     Optimality, Feasibility

                                                                                                                                                                                                                                             Optimality, Feasibility
                           100                                                                                                                 100                                                                                                                     100
                                                                                                             .0050                                                                                                                   .0050                                                                                                         .0050

                                -2                                                                           .0045                             10
                                                                                                                                                    -2                                                                               .0045                             10
                                                                                                                                                                                                                                                                            -2                                                                     .0045

                                                                                                             .0040                                                                                                                   .0040                                                                                                         .0040
                           10-4                                                                              .0035                             10-4                                                                                  .0035                             10-4                                                                        .0035

                                -6
                                                                                                             .0030                                  -6
                                                                                                                                                                                                                                     .0030                                  -6
                                                                                                                                                                                                                                                                                                                                                   .0030
                           10                                                                                                                  10                                                                                                                      10
                                                                                                             .0025                                                                                                                   .0025                                                                                                         .0025

                           10-8                                                                              .0020                             10-8                                                                                  .0020                             10-8                                                                        .0020
                                     0   20       40         60          80      100                                                                     0            20             40                 60                                                                       0           10                     20
                                                 Design iteration                                                                                                       Design iteration                                                                                                   Design iteration

                                                                              Optimized                                                                                                             Optimized                                                                                                      Optimized
                                                                              Initial                        .015                                                                                   Initial                          .015                                                                          Initial                         .015

                       -0.6                                                                                                                -0.6                                                                                                                    -0.6
                                                                                                             .012                                                                                                                    .012                                                                                                          .012

                           0.0                                                                                                                 0.0                                                                                                                     0.0
                                                                                                             .009                                                                                                                    .009                                                                                                          .009
 Cp

                                                                                                                     Cp

                                                                                                                                                                                                                                             Cp
                           0.6                                                                                                                 0.6                                                                                                                     0.6
                                                                                                               Cf

                                                                                                                                                                                                                                       Cf

                                                                                                                                                                                                                                                                                                                                                     Cf
                                                                                                             .006                                                                                                                    .006                                                                                                          .006
                           1.2                                                                                                                 1.2                                                                                                                     1.2
                                               upper                                                                                                                  upper
                                                                                                             .003                                                                                                                    .003                                                                                                          .003
                           1.8                           lower                                                                                 1.8                              lower                                                                                  1.8                                      upper
                                                                                lower                                                                                                           upper                                                                                                           upper lower
                                                                    upper                                                                                                                               lower                                                                                                         lower
                           2.4                                                                               .000                              2.4                                                                                   .000                              2.4                                                                         .000
                                     0   0.2           0.4         0.6         0.8                       1                                               0      0.2           0.4         0.6           0.8                      1                                               0   0.2      0.4         0.6         0.8                      1
                                                             x/c                                                                                                                    x/c                                                                                                             x/c

(a) Baseline Optimization Problem                                                                                                                            (b) Additional Constraints                                                      (c) Re-optimized w/ Baseline Opti-
filler                                                                                                                                                                  filler                                                               mization Problem

Figure 5.2: Optimization convergence histories and cross-sectional profiles of the initial and optimized designs
produced by drag minimizations with the baseline (Equations 5.1–5.4) and more constrained optimization
problem (Equations 5.5–5.9) at the Cessna 172R conditions on the L1 grid level.

   To confirm that multiple local minima exist in the baseline design space, the more con-
strained optimized design was used as an initial geometry with the relaxed constraints of the
baseline optimization problem. The results in Figure 5.2c and Table 5.3 demonstrate that
the optimization algorithm is able to reduce drag by an additional count by decreasing the
trailing-edge thickness to create stronger adverse pressure gradients near the trailing edge,
which reduce the skin friction in the turbulent boundary layer. However, the new design
maintains a similar profile to the initial design produced with the additional constraints. The
results in Figure 5.2c and Table 5.3 confirm that the Cessna 172R design space has at least
two distinct local minima, one with a higher design angle of attack and extended regions of
laminar flow on the lower surface of the airfoil and one with more aft-loading and more bal-
anced extents of laminar flow on the upper and lower surfaces, with the latter outperforming
the former. The results demonstrate the significant risk multi-modality presents at these
design conditions, as the baseline optimization problem produces an inferior design, and
highlights the importance of efficient global optimization techniques, such as gradient-based
multi-start methods [23, 179].

<!-- ===== PDF page 110 ===== -->

90                                                                                                                                                                                               CHAPTER 5. RESULTS: OPTIMIZATION

Table 5.4: Structured multi-block O-grid dimensions for the two-dimensional optimizations at the Dash8-
Q400 design conditions (see Table 5.1).
              grid level                                     chord x off-wall nodes                                          avg/max ∆s × 10−6 (chord)                                                                  avg/max y+
                 L1                                                 561x121                                                         0.65 / 0.74                                                                          0.26 / 0.53
                 L2                                                1121x121                                                              "                                                                                    "

                                                                         Optimality                                                                                2
                                                                                                                                                                                           Optimality
                                            10                           Feasibility                                        .0060                             10                           Feasibility                                         .0060

                                                                                                       Merit function

                                                                                                                                                                                                                          Merit function
                                                                         Merit function                                                                                                    Merit function
                                                                                                                            .0055                                                                                                              .0055
                  Optimality, Feasibility

                                                                                                                                    Optimality, Feasibility
                                            100                                                                                                               100

                                                                                                                            .0050                                                                                                              .0050
                                                 -2                                                                                                                -2
                                            10                                                                                                                10
                                                                                                                            .0045                                                                                                              .0045
                                            10-4                                                                                                              10-4
                                                                                                                            .0040                                                                                                              .0040

                                                 -6                                                                                                                -6
                                            10                                                                                                                10
                                                                                                                            .0035                                                                                                              .0035

                                            10-8                                                                            .0030                             10-8                                                                             .0030
                                                      0           5       10         15         20     25                                                               0           5       10       15            20     25
                                                                        Design iteration                                                                                                  Design iteration

                                                                                                Optimized                                                                                                          Optimized
                                                                                                Initial                     .015                                                                                   Initial                     .015
                                        -1.2                                                                                                              -1.2

                                        -0.6                                                                                .012                          -0.6                                                                                 .012

                                            0.0                                                                                                               0.0
                                                                                                                            .009                                                                                                               .009
                                            0.6                                                                                                               0.6
                  Cp

                                                                                                                                    Cp
                                                                                                                              Cf

                                                                                                                                                                                                                                                 Cf
                                            1.2                                                                             .006                              1.2                                                                              .006

                                            1.8                                                                                                               1.8
                                                          upper                                                                                                             upper
                                                                        upper
                                                                                                                            .003                                                                                                               .003
                                                                                                                                                                                           upper
                                            2.4                           lower                                                                               2.4                                  lower
                                                                                  lower                                                                                                                    lower
                                            3.0                                                                             .000                              3.0                                                                              .000
                                                      0           0.2       0.4           0.6    0.8                    1                                               0           0.2      0.4         0.6        0.8                    1
                                                                                  x/c                                                                                                              x/c

                                                                             (a) L1                                                                                                          (b) L2

Figure 5.3: Optimization convergence histories and cross-sectional profiles of the initial and optimized designs
produced by drag minimizations with the baseline optimization problem (Equations 5.1–5.4) at the De
Havilland Dash8-Q400 conditions with varying streamwise grid resolution.

5.1.2    De Havilland Dash8-Q400

Based on the results from the Cessna 172R drag minimizations, grids with streamwise reso-
lutions representative of the L1 and L2 grid levels are investigated for the drag minimizations
at the higher Mach number (0.60), Reynolds number (15.7 × 106 ), and lift coefficient (0.42)
of the Dash8-Q400 design conditions. The Dash8-Q400 design conditions are presented in
Table 5.1, with the grid details for the L1 and L2 grid levels presented in Table 5.4. Similar
to the Cessna 172R case, a streamwise grid-resolution study is performed using the L1 and
L2 grids and the baseline optimization problem (Equations 5.1–5.4) without the pitching-
moment and tighter thickness constraints in order to evaluate the capabilities of the opti-
mization algorithm to explore the Dash8-Q400 design space. The optimization histories and
cross-sectional profiles for the initial and optimized designs are presented in Figure 5.3.
  The results demonstrate that there is a significant difference in both airfoil profile shape
and performance of the designs produced by the optimization algorithm on the L1 and L2

<!-- ===== PDF page 111 ===== -->

5.1. AIRFOIL OPTIMIZATION                                                                    91

grid levels. Where the L1 grid level provides adequate resolution for the lower Reynolds
number drag minimizations at the Cessna 172R design conditions, the optimizations at
the higher Mach and Reynolds numbers of the Dash8-Q400 design conditions require finer
streamwise grid resolutions. There is a decrease in the transition length predicted by the
transition model at the Dash8-Q400 design conditions, illustrated in Figure 5.3, relative
to the transition length produced at the Cessna 172R conditions, which are presented in
Figure 5.1. This decrease in transition length due to the increase in Reynolds number can
help to explain the increased sensitivity of the optimizations to streamwise grid resolution,
as finer streamwise grid spacing is required in order to maintain a similar resolution in the
transition region to the Cessna 172R optimizations. Instead of delaying transition on the
L1 grid, the optimization algorithm prioritizes developing strong adverse pressure gradients
downstream of the transition location, which reduce the turbulent skin friction coefficient.
This is similar to the behaviour of the drag minimization at the Cessna 172R design condi-
tions on the L0 grid level (Figure 5.1a). The results suggest that the streamwise resolution
of the grid must be increased as the Reynolds number increases in order to accurately re-
solve the large gradients in the transition region. However, this can lead to a significant
increase in computational cost, especially for the higher Reynolds numbers associated with
single-aisle transport aircraft. An alternative approach is investigated in Appendix C, where
the transition region length is increased in order to increase the effective grid resolution in
the transition region. The results demonstrate that increasing the transition region length
greatly reduces the streamwise grid requirements for the free-transition optimizations.

   Similar to the designs in Figure 5.1, the designs in Figure 5.3 optimized using the baseline
optimization problem produce excessively thin trailing edges. In order to improve manufac-
turability and to further explore the design space, drag minimizations are performed on
the L2 grid level using the optimization problem with the additional constraints (pitching-
moment and tighter thickness-to-chord ratio constraints) introduced in Section 5.1.1 (Equa-
tions 5.5–5.9). The optimized designs produced using the two optimization problems on the
L2 grid level are presented in Figure 5.4. The results demonstrate that the optimization
with the additional constraints produces a more aft-loaded design with increased trailing-
edge thickness and a more gradual concave pressure recovery. However, in contrast to the
behaviour observed at the Cessna 172R design conditions (Figure 5.2), the more constrained
design produces increased drag relative to the design optimized with the baseline optimiza-
tion problem. Although the two optimized designs both delay transition to approximately
45% chord on the upper surface of the airfoil, the design optimized with the baseline optimiza-
tion problem maintains a larger region of laminar flow on the lower surface. In addition, the
steep pressure recovery produced by the baseline optimization problem decelerates the flow,
producing decreased skin friction in the turbulent boundary layer, resulting in a decrease
in skin friction drag relative to the more constrained design. Details of the aerodynamic
performance of the initial blunt trailing-edge RAE2822 airfoil and two optimized designs on

<!-- ===== PDF page 112 ===== -->

92                                                                                                                                                                                                  CHAPTER 5. RESULTS: OPTIMIZATION

                                                                         Optimality                                                                                 2
                                                                                                                                                                                                Optimality
                                            10                           Feasibility                                         .0060                             10                               Feasibility                                         .0060

                                                                                                        Merit function

                                                                                                                                                                                                                               Merit function
                                                                         Merit function                                                                                                         Merit function
                                                                                                                             .0055                                                                                                                  .0055

                  Optimality, Feasibility

                                                                                                                                     Optimality, Feasibility
                                            100                                                                                                                100

                                                                                                                             .0050                                                                                                                  .0050
                                                 -2                                                                                                                 -2
                                            10                                                                                                                 10
                                                                                                                             .0045                                                                                                                  .0045
                                            10-4                                                                                                               10-4
                                                                                                                             .0040                                                                                                                  .0040

                                                 -6                                                                                                                 -6
                                            10                                                                                                                 10
                                                                                                                             .0035                                                                                                                  .0035

                                                 -8                                                                                                                 -8
                                            10                                                                               .0030                             10                                                                                   .0030
                                                      0           5       10       15            20     25                                                               0                 10                  20         30
                                                                        Design iteration                                                                                                   Design iteration

                                                                                                 Optimized                                                                                                            Optimized
                                                                                                 Initial                     .015                                                                                     Initial                       .015
                                        -1.2                                                                                                               -1.2

                                        -0.6                                                                                 .012                          -0.6                                                                                     .012

                                            0.0                                                                                                                0.0
                                                                                                                             .009                                                                                                                   .009
                                            0.6                                                                                                                0.6
                  Cp

                                                                                                                                     Cp
                                                                                                                               Cf

                                                                                                                                                                                                                                                      Cf
                                            1.2                                                                              .006                              1.2                                                                                  .006

                                            1.8                                                                                                                1.8
                                                          upper                                                                                                              upper
                                                                                                                             .003                                                                                                                   .003
                                                                         upper                                                                                                                   lower
                                            2.4                                  lower                                                                         2.4                                       upper
                                                                                         lower                                                                                                           lower

                                            3.0                                                                              .000                              3.0                                                                                  .000
                                                      0           0.2      0.4         0.6        0.8                    1                                               0           0.2          0.4          0.6     0.8                      1
                                                                                 x/c                                                                                                                     x/c

                                (a) Baseline Optimization Problem                                                                                                            (b) Additional Constraints

Figure 5.4: Optimization convergence histories and cross-sectional profiles of the initial and optimized designs
produced by drag minimizations with the baseline (Equations 5.1–5.4) and more constrained optimization
problem (Equations 5.5–5.9) at the De Havilland Dash8-Q400 conditions on the L2 grid level.

Table 5.5: Aerodynamic performance of the initial and optimized designs at the Dash8-Q400 conditions on
the L2 grid level.
                                                                                                        Cd (cnts.)                                                        Cl                Cm                        L/D                       aoa (deg.)
                        Initial                                                                           54.85                                                          0.42              −0.074                    76.438                       1.329
            Baseline Optimization Problem                                                                 35.85                                                          0.42              −0.017                   117.124                       1.998
                Additional Constraints                                                                    40.18                                                          0.42              −0.074                   104.521                       0.722

the L2 grid level are presented in Table 5.5. The results demonstrate that the optimization
algorithm was able to reduce drag by approximately 19 and 15 drag counts using the base-
line and more constrained optimization problems, respectively. Therefore, adding a pitching
moment constraint and more conservative thickness-to-chord ratio constraints incurs a drag
penalty of 4 drag counts. Further investigation is required to determine if the design pro-
duced using the optimization problem with additional constraints is a local minimum of the
baseline optimization problem design space.
   It is important to note that although the results in Figure 5.4 demonstrate that the merit
function, here representing the drag coefficient, has greatly reduced, optimality has not
converged. Instead, the optimization algorithm exits after being unable to further improve
the design. These results suggest the possibility of noise in the gradient at this higher
Reynolds number design condition. Preliminary results indicate that this is due to the
Cartesian velocity gradient-based pressure gradient parameter (Equations 2.10–2.14) and

<!-- ===== PDF page 113 ===== -->

5.1. AIRFOIL OPTIMIZATION                                                                            93

Table 5.6: Structured multi-block O-grid dimensions for the two-dimensional optimizations at the sweep-
corrected Boeing 737-800 design conditions (see Table 5.1).
             grid level   chord x off-wall nodes   avg/max ∆s × 10−6 (chord)   avg/max y+
                L1               561x121                  0.54 / 0.60           0.25 / 0.50
                L2              1121x121                       "                     "

will be investigated in more detail in future work.

5.1.3    Boeing 737-800

To explore the performance of the free-transition optimization framework when applied to air-
foil optimizations at design conditions similar to that of a single-aisle aircraft, lift-constrained
drag minimizations were performed at the sweep-corrected Boeing 737-800 design conditions
presented in Table 5.1. It is important to note again that because these optimizations are
two-dimensional the effects of crossflow instabilities are not included. However, the purpose
of this design problem is to evaluate the capabilities of the optimization algorithm to re-
duce TS instabilities at the high Mach number (0.71 corr.), Reynolds number (20.3 × 106 ),
and lift coefficient (0.50) typical of a single-aisle transport aircraft. Drag minimizations
are performed on grids with streamwise grid resolutions representative of the L1 and L2
grid levels previously investigated, with the grid details presented in Table 5.6. Optimiza-
tions are performed using the baseline optimization problem (Equations 5.1–5.4) without
the pitching-moment and tighter thickness-to-chord ratio constraints.
   The results, which are presented in Figure 5.5, demonstrate that the grid-resolution re-
quirements identified in the Cessna 172R and De Havilland Dash8-Q400 drag minimizations
are amplified at the higher Mach and Reynolds numbers of the Boeing 737-800 conditions.
Specifically, as the Reynolds number increases the transition region decreases further. There-
fore, the drag minimization on the L1 grid (Figure 5.5a) fails to significantly move the
transition fronts aft as the effective grid resolution in the transition region decreases. The
optimization algorithm appears to encounter difficulty trading a decrease in viscous drag
associated with delaying boundary-layer transition with the increase in wave drag associated
with the development of a favourable pressure gradient to delay transition. This is also
demonstrated by the optimization history in Figure 5.5a, which illustrates that convergence
of the merit function has stalled. The performance of the optimization algorithm improves
significantly going from the L1 to the L2 grid level. The optimization algorithm successfully
produces a favourable pressure gradient on the upper surface of the airfoil to delay transi-
tion and in the process produces a shock wave at approximately 50% chord. The results in
Figure 5.5b demonstrate that provided enough grid resolution in the transition region the
optimization algorithm is able to successfully trade a decrease in viscous drag produced by
delaying boundary-layer transition with an increase in wave drag.
   The aerodynamic performance metrics for the initial and optimized design on the L2 grid

<!-- ===== PDF page 114 ===== -->

94                                                                                                                                                                                                 CHAPTER 5. RESULTS: OPTIMIZATION

                                                                    Optimality                                                                                        2
                                                                                                                                                                                               Optimality
                                            10                      Feasibility                                                .0060                             10                            Feasibility                                             .0060

                                                                                                          Merit function

                                                                                                                                                                                                                                  Merit function
                                                                    Merit function                                                                                                             Merit function
                                                                                                                               .0055                                                                                                                   .0055

                  Optimality, Feasibility

                                                                                                                                       Optimality, Feasibility
                                            100                                                                                                                  100

                                                                                                                               .0050                                                                                                                   .0050
                                                 -2                                                                                                                   -2
                                            10                                                                                                                   10
                                                                                                                               .0045                                                                                                                   .0045
                                            10-4                                                                                                                 10-4
                                                                                                                               .0040                                                                                                                   .0040

                                                 -6                                                                                                                   -6
                                            10                                                                                                                   10
                                                                                                                               .0035                                                                                                                   .0035

                                                 -8                                                                                                                   -8
                                            10                                                                                 .0030                             10                                                                                  .0030
                                                      0     5       10        15         20          25                                                                    0   2           4        6         8       10     12                    14
                                                                  Design iteration                                                                                                         Design iteration

                                                                                              Optimized                                                                                                                Optimized
                                                                                              Initial                          .015                                                                                    Initial                         .015
                                        -1.2                                                                                                                 -1.2

                                        -0.6                                                                                   .012                          -0.6                                                                                      .012

                                            0.0                                                                                                                  0.0
                                                                                                                               .009                                                                                                                    .009
                                            0.6                                                                                                                  0.6
                  Cp

                                                                                                                                       Cp
                                                                                                                                 Cf

                                                                                                                                                                                                                                                         Cf
                                            1.2                                                                                .006                              1.2                                                                                   .006

                                            1.8                                                                                                                  1.8
                                                          upper                                                                                                                    upper
                                                                  upper                                                        .003                                                                                                                    .003
                                                                                                                                                                                                lower
                                            2.4                            lower                                                                                 2.4                           lower
                                                                           lower                                                                                                                              upper

                                            3.0                                                                                .000                              3.0                                                                                   .000
                                                      0     0.2       0.4          0.6         0.8                         1                                               0       0.2           0.4          0.6          0.8                     1
                                                                            x/c                                                                                                                         x/c

                                                                          (a) L1                                                                                                                  (b) L2

Figure 5.5: Optimization convergence histories and cross-sectional profiles of the initial and optimized designs
produced by drag minimizations with the baseline optimization problem (Equations 5.1–5.4) at the sweep-
corrected Boeing 737-800 conditions with varying streamwise grid resolution.

Table 5.7: Aerodynamic performance of the initial and optimized designs at the sweep-corrected Boeing
737-800 conditions on the L2 grid level.
                                                                Cd (cnts.)               Cd,p                              Cd,f                          Cl                     Cm                        L/D                     aoa (deg.)
                Initial                                           54.33                  18.62                             35.71                        0.50                   −0.086                    92.018                     1.291
               Optimized                                          42.49                  13.30                             29.19                        0.50                   −0.088                    117.519                    1.133

level are presented in Table 5.7, including drag breakdowns for the initial and optimized
designs. The results demonstrate that the optimized design with a larger extent of laminar
flow on the upper surface of the airfoil produces decreased skin friction drag as well as a
net decrease in pressure drag, as the reduction in the pressure component of viscous drag is
larger than the increase in wave drag.
   Similar to the results at the Dash8-Q400 design conditions presented in Figure 5.4, the
optimization histories presented in Figure 5.5 demonstrate that the optimization algorithm
at this high Reynolds number design condition struggles to reduce optimality, and instead
exits without converging to a local minimum. Although the results presented in Appendix C
demonstrate that increasing the transition region length reduces grid requirements and com-
putational cost for higher Reynolds numbers, it does not appear to improve the convergence
of optimality. Again, preliminary investigations suggest that this is due to the Cartesian ve-
locity gradient-based pressure gradient parameter (Equations 2.10–2.14). Future work will
focus on alternative formulations of the pressure gradient parameter in order to improve the

<!-- ===== PDF page 115 ===== -->

5.2. INFINITE SWEPT WING OPTIMIZATION                                                         95

convergence of optimality and to ensure that the optimization algorithm can produce a local
minimum at these higher Reynolds number conditions.

5.2    Infinite Swept Wing Optimization

An infinite swept wing optimization was performed at design conditions similar to a transonic
strut-braced wing aircraft [22]. The goal of the drag minimization at this design condition
is to evaluate the ability of the optimization algorithm to delay both TS and stationary
crossflow instabilities. Specifically, a lift-constrained drag minimization was conducted using
an initial geometry consisting of an RAE2822 airfoil extruded with a 30-degree sweep. The
infinite swept wing O-grid consists of 857, 103, and 11 nodes in the streamwise, off-wall,
and spanwise directions, respectively, with an average y+ value of 0.50. Design conditions
were specified as a Mach number of 0.80, Reynolds number of 12.12 × 106 , turbulence in-
tensity of 0.10%, surface roughness of 1µinch, with a target cruise lift coefficient of 0.56.
An optimization problem similar to the baseline optimization problem (Equations 5.1–5.4)
was used for the infinite swept wing optimization; however, to expand the design space the
cross-sectional area constraint was removed and more conservative thickness bounds were
placed on the FFD design variable pairs. The infinite swept wing optimization problem is
defined as follows:

                                   min CD (X ),                                          (5.10)
                                    X

                                    s.t. CL = 0.56,                                      (5.11)
                                         t/c ≥ 0.85(t/c)init .                           (5.12)

Similar to the airfoil optimizations, 6 streamwise FFD design variables are applied to each
of the upper and lower surfaces, with the infinite swept wing root and tip design variables
constrained to be equal and the leading- and trailing-edge design variables constrained to
move symmetrically, for a total of 10 effective geometric design variables plus angle of attack.
   The results from the drag minimization are presented in Figure 5.6. The optimization
history illustrated in Figure 5.6a demonstrates that the optimization algorithm successfully
reduces drag by approximately 12 counts. The optimized cross-sectional profiles in Fig-
ure 5.6b demonstrate that the optimization algorithm produces near-zero pressure gradients
on the upper and lower surfaces of the wing, which delay transition to 60% and 40% chord,
respectively. While the initial design already featured a near-zero upper surface pressure
gradient, it is important to note that the optimization algorithm has removed the favourable
pressure gradient on the lower surface of the wing to prevent the growth of crossflow in-
stabilities and therefore delay transition on the lower surface. The optimization algorithm
produced these near-zero pressure gradients on the upper and lower surfaces by reducing the
wing thickness and moving the aerodynamic loading aft.

<!-- ===== PDF page 116 ===== -->

96                                                                                                                                                       CHAPTER 5. RESULTS: OPTIMIZATION

                                                                    Optimality                                                                                               Optimized
                                            10                      Feasibility                               .0070                                                          Initial         .015

                                                                                             Merit function
                                                                                                                       -1.2
                                                                    Merit function
                                                                                                              .0065    -0.6

                  Optimality, Feasibility
                                            100                                                                                                                                              .012

                                                                                                              .0060        0.0
                                                 -2
                                            10                                                                                                                                               .009
                                                                                                                           0.6

                                                                                                                      Cp

                                                                                                                                                                                               Cf
                                                                                                              .0055
                                            10-4                                                                           1.2                                                               .006
                                                                                                              .0050
                                                                                                                           1.8                                       upper   upper
                                                 -6
                                            10                                                                                                 lower                                         .003
                                                                                                              .0045                                    lower
                                                                                                                           2.4

                                            10-8                                                              .0040        3.0                                                               .000
                                                      0        50            100       150                                       0       0.2           0.4           0.6        0.8      1
                                                               Design iteration                                                                                x/c

                                                          (a) Optimization history                                                   (b) Cross-sectional profiles

Figure 5.6: Infinite swept wing drag minimization at the transonic strut-braced wing aircraft design condi-
tions.

Table 5.8: Aerodynamic performance of the initial and optimized infinite swept wing designs at the transonic
strut-braced wing aircraft design conditions.
                                                              CD (cnts.)             CD,p                     CD,f          CL             CM                     L/D                 aoa (deg.)
                Initial                                         58.10                25.97                    32.13         0.56          −0.23                   96.38                  2.10
               Optimized                                        46.12                18.99                    27.13         0.56          −0.28                  121.42                  0.94

   The integrated forces for the initial and optimized infinite swept wing designs are pre-
sented in Table 5.8. The results demonstrate that the pressure and skin friction drag compo-
nents account for approximately 60% and 40% of the total drag reduction, respectively. The
relatively larger reduction in pressure drag can be explained by Figure 5.6b, which illustrates
that the optimized design simultaneously reduces the strength of the shock while increasing
the laminar extent of the boundary layer, which the optimizer achieved by reducing the
maximum thickness-to-chord ratio and the cross-sectional area. Future work will investigate
this optimization problem with a more conservative cross-sectional area constraint, similar
to that used in the airfoil optimizations.

<!-- ===== PDF page 117 ===== -->

Chapter 6

Contributions, Conclusions, and
Recommendations

6.1    Contributions and Conclusions

The objective of this work is to develop and investigate a high-fidelity aerodynamic shape
optimization framework based on the RANS equations that includes the effects of boundary-
layer transition which can be applied to the design of subsonic and transonic NLF wings.
Here we revisit the intermediate goals presented in Chapter 1 which were completed in order
to satisfy this primary objective, and discuss the contributions presented in the current work.

  1. Extend a three-dimensional RANS solver with a local transition predic-
     tion model. The γ-Re    ˜ θt LM2015 transition model with helicity-based crossflow cor-
     relations was modified and coupled to the one-equation SA turbulence model, while
     maintaining a fully local formulation, reducing the computational cost of the original
     four-equation system (SA-LM2015) [139]. Smooth approximations to non-differentiable
     and stiff source-term functions were introduced to eliminate discontinuities and improve
     the numerical behaviour of the model (SA-sLM2015) [139]. Compressibility corrections
     were developed and applied to include the stabilizing effects of flow compressibility on
     TS and stationary crossflow instabilities, extending the LM2009 and LM2015 empirical
     correlations, respectively, to transonic flight regimes typical of a commercial transport
     aircraft (SA-sLM2015cc) [142].

  2. Develop a framework for achieving deep and efficient iterative convergence.
     Modifications to a Newton-Krylov-Schur flow solver were introduced to achieve robust,
     efficient convergence of the turbulence and transition model equations with a fully cou-
     pled, fully implicit solution strategy [139, 140]. A source-term time step restriction,
     equation and variable scaling measures, and a solution update damping algorithm, in-
     cluding a physics-based restriction and an unsteady residual backtracking line search,
     were developed and implemented to prevent unstable solution updates [139, 140]. Con-

<!-- ===== PDF page 118 ===== -->

98                                CHAPTER 6. CONTRIBUTIONS, CONCLUSIONS, AND RECOMMENDATIONS

       vergence studies presented in Appendix B with varying levels of coupling for the tur-
       bulence and transition model equations demonstrate that the fully coupled solution
       strategy with the source-term time step restriction provides robust nonlinear conver-
       gence [140]. Furthermore, the source-term time step restriction has been found to
       improve the robustness of fully turbulent simulations using the SA turbulence model.

     3. Validate the model for a range of subsonic and transonic flow regimes. Two-
        and three-dimensional subsonic transition test cases demonstrate that the SA-LM2015
        and SA-sLM2015 transition models accurately predict boundary-layer transition due to
        a variety of transition mechanisms, with the SA-sLM2015 model demonstrating signifi-
        cantly improved numerical behaviour [139]. The compressibility corrections and the un-
        derlying smooth transition model were applied to two- and three-dimensional transonic
        test cases for which transition onset has been measured, with the results demonstrat-
        ing that the compressibility corrections produce substantially improved agreement with
        the experimental transition locations, particularly for higher Reynolds number applica-
        tions. Grid-refinement studies presented in Appendix B demonstrate that the smooth
        transition model with and without the compressibility corrections, SA-sLM2015cc and
        SA-sLM2015, respectively, produces similar grid convergence to fully turbulent simula-
        tions performed using the SA turbulence model [140].

     4. Integrate the model in a discrete-adjoint gradient-based optimization algo-
        rithm. The flow and mesh adjoint problems contain partial derivatives of the total
        residual, namely, the flow Jacobian and the metric linearization. To obtain these par-
        tial derivatives, the modified turbulence model and the transition model equations are
        mostly differentiated analytically, with the linearization of the simultaneous approxima-
        tion terms, which are difficult to linearize analytically, formed using the complex-step
        method. The partial derivatives and the gradient are verified by comparing with ap-
        proximations formed using the complex-step method and with a finite-difference approx-
        imation, respectively, enabling free-transition aerodynamic shape optimization [141].

     5. Apply the framework to drag-minimization studies over a range of design
        conditions. The aerodynamic shape optimization framework incorporating transi-
        tion effects was applied to airfoil and infinite swept wing drag minimizations at design
        conditions ranging from light to single-aisle aircraft [141]. The results highlight the
        risks posed by multi-modality at design conditions representative of a light aircraft
        and demonstrate that the ability of the gradient-based optimization algorithm to delay
        boundary-layer transition is sensitive to the streamwise grid resolution in the transition
        regions. The streamwise grid requirements increase as the transition length decreases
        with an increasing Reynolds number. Drag minimizations at transonic flight conditions
        demonstrate that the optimization framework can successfully trade a decrease in vis-
        cous drag with an increase in wave drag, and an optimization of an infinite swept wing

<!-- ===== PDF page 119 ===== -->

6.2. RECOMMENDATIONS                                                                       99

      demonstrates that the framework can delay both TS and stationary crossflow insta-
      bilities. A new intermittency source-term formulation is investigated in Appendix C
      to reduce the streamwise grid requirements for NLF optimization by increasing the
      transition region length.

6.2     Recommendations

Recommendations for future work are presented in order to improve the predictive capability
and performance of the free-transition optimization framework, and to further evaluate its
capabilities for exploring the design of transonic NLF wings.

  • Improve turnaround times for free-transition simulations. The introduction of
    smooth source terms significantly improved the numerical behaviour of the correlation-
    based transition model and combined with the modifications to the Newton-Krylov
    solution strategy, including the source-term time step restriction, facilitated robust
    iterative convergence over a range of flight conditions. However, the smooth transition
    model consists of two additional transport equations, which significantly increase the
    memory requirements and computational cost of free-transition simulations relative to
    those performed assuming fully turbulent flow. Furthermore, free-transition simulations
    still often require significantly more linear and nonlinear iterations to converge due to
    the complex physical and numerical aspects of boundary-layer transition. To improve
    turnaround times and facilitate rapid NLF trade space exploration, one option is to
    improve and automate the domain decomposition strategy in Diablo, which is currently
    performed manually during grid generation, while another is to simplify the transition
    model. Recently, one-equation and algebraic transition models have been introduced
    in the literature. Simplifying the smooth transition model using similar approaches
    in order to reduce the number of transport equations while maintaining its predictive
    capability could reduce the wallclock times and computational resources required for
    these free-transition simulations.

  • Validate the smooth transition model for higher Reynolds number applica-
    tions. The smooth transition model has been validated for Mach and Reynolds numbers
    ranging from 0.10 ≤ M ≤ 0.856 and 1.80 × 106 ≤ Re ≤ 15 × 106 , respectively. However,
    as the experimental campaign by Lynde et al. [104, 105] demonstrates, free-transition
    experimental data from wind tunnels at flight-scale Reynolds numbers is difficult to ob-
    tain due to bypass transition resulting from high freestream turbulence and the relative
    significance of surface imperfections due to thin boundary layers. It is important to
    validate the smooth transition model over a larger range of flow conditions at higher
    Reynolds numbers as this data becomes available, as the ability of the optimization
    algorithm to produce NLF designs is dependent on its predictive capability.

<!-- ===== PDF page 120 ===== -->

100                               CHAPTER 6. CONTRIBUTIONS, CONCLUSIONS, AND RECOMMENDATIONS

      • Extend the framework to include additional transition mechanisms. Although
        the current framework can accurately predict natural, bypass, separation-induced, and
        crossflow-induced transition, the effects of Görtler and attachment-line instabilities are
        not explicitly modelled. These transition mechanisms could be accounted for by either
        extending the empirical correlations in the smooth transition model to predict these
        instabilities or by developing geometric constraints for the optimization framework to
        prevent their growth.

      • Reduce streamwise grid requirements for free-transition optimizations. The
        results in Chapter 5 demonstrate that the ability of the optimization algorithm to delay
        boundary-layer transition is dependent on the streamwise grid resolution in the tran-
        sition region, and that these grid requirements increase as the transition region length
        decreases with increasing Reynolds number. Appendix C presents promising results
        for a modification to the intermittency transport equation source terms to maintain a
        consistent transition length across a range of flight conditions, which can reduce the
        computational cost of the free-transition optimizations. However, future work must be
        completed to extend this modification to include the stationary crossflow source terms.

      • Improve convergence of the optimization algorithm for high Reynolds num-
        ber applications. The drag minimizations in Chapter 5 demonstrate that the opti-
        mization algorithm struggles to reduce optimality at higher Reynolds numbers. While
        the results in Appendix C demonstrate that increasing the transition region length re-
        duces streamwise grid requirements, it does not appear to improve the convergence of
        optimality. Future work must be completed to improve the convergence of the opti-
        mization algorithm at higher Reynolds numbers. This will require identifying possible
        sources of noise in the gradient, including investigating different formulations of the
        local pressure gradient parameter.

      • Evaluate the transonic NLF wing design space. The application of the free-
        transition optimization framework to the design of transonic NLF wings could help
        expand our understanding of the design space. Specifically, understanding the potential
        drag reduction benefit of NLF designs, the tradeoffs between viscous drag and wave
        drag, and the balance of wing thickness and sweep can help to accelerate the design of
        transonic NLF wings.

      • Investigate the design of robust transonic NLF wings. There are several chal-
        lenges that must be overcome for the adoption of NLF wings on commercial transport
        aircraft, such as surface imperfections resulting from dirt and insect strikes, access pan-
        els, manufacturing tolerances, leading-edge high-lift devices, and anti-icing coatings.
        However, the design of NLF configurations sufficiently robust to changes in the distur-
        bance environment and flow conditions using multi-point optimizations could help to
        reduce the barrier to entry of this technology.

<!-- ===== PDF page 121 ===== -->

References

 [1] B. J. Abu-Ghannam and R. Shaw. Natural transition of boundary layers, the effects
     of turbulence, pressure gradient, and flow history. Journal of Mechanical Engineering
     Science, 22(5):213–228, 1980.
 [2] S. Allmaras. Multigrid for the 2-D compressible Navier-Stokes equations. 14th Com-
     putational Fluid Dynamics Conference, AIAA 1999-3336, 1999.
 [3] S. R. Allmaras and F. T. Johnson. Modifications and clarifications for the implemen-
     tation of the Spalart-Allmaras turbulence model. Seventh International Conference on
     Computational Fluid Dynamics, pages 1–11, ICCFD7 Paper 1902, 2012.
 [4] O. Amoignon, J. Pralits, A. Hanifi, M. Berggren, and D. Henningson. Shape opti-
     mization for delay of laminar-turbulent transition. AIAA Journal, 44(5):1009–1024,
     2006.
 [5] D. Arnal. Laminar-turbulent transition problems in supersonic and hypersonic flows.
     Special Course on Aerothermodynamics of Hypersonic Vehicles, AGARD, AGARD-R-
     761, 1989.
 [6] D. Arnal. Boundary layer transition: predictions based on linear theory. Special Course
     on Progress in Transition Modelling, AGARD, AGARD-R-793, 1993.
 [7] D. Arnal. Transition prediction in transonic flows. Symposium Transsonicum III. In-
     ternational Union of Theoretical and Applied Mechanics, Springer, Berlin, Heidelberg,
     1988.
 [8] D. Arnal, G. Casalis, and R. Houdeville. Practical transition prediction methods:
     Subsonic and transonic flows. VKI Lecture Series: Advances in Laminar-Turbulent
     Transition Modeling, pages 1–34, 2009.
 [9] D. Arnal, M. Habiballah, and E. Coutols. Laminar instability theory and transition
     criteria in two and three-dimensional flow. La Recherche Aérospatiale (English Edi-
     tion)(ISSN 0379-380X), pages 45–63, No. 2, 1984.
[10] D. Arnal, R. Houdeville, A. Seraudie, and O. Vermeersch.                  Overview of
     laminar-turbulent transition investigations at ONERA Toulouse. 41st AIAA Fluid
     Dynamics Conference and Exhibit, 2011.
[11] B. Aupoix, D. Arnal, H. Bezard, B. Chaouat, and F. Chedevergne. Transition and
     turbulence modeling. AerospaceLab, pages 1–13, 2011.
[12] B. Baldwin and H. Lomax. Thin-layer approximation and algebraic model for separated

<!-- ===== PDF page 122 ===== -->

102                                                                               REFERENCES

      turbulent flows. 16th Aerospace Sciences Meeting, 1978.
 [13] D. P. Bertsekas. Nondifferentiable optimization via approximation. Mathematical Pro-
      gramming Study 3, pages 1–25, 1975.
 [14] M. Blanco and D. W. Zingg. Newton-Krylov algorithm with a loosely coupled turbu-
      lence model for aerodynamic flows. AIAA Journal, 45(5):980–987, 2007.
 [15] D. M. Bushnell. Overview of aircraft drag reduction technology. AGARD Report 786
      (Special Course on Skin Friction Drag Reduction), 1992.
 [16] R. L. Campbell and M. N. Lynde. Building a practical natural laminar flow design
      capability. 35th AIAA Applied Aerodynamics Conference, AIAA 2017-3059, 2017.
 [17] R. L. Campbell and M. N. Lynde. Natural laminar flow design for wings with mod-
      erate sweep. 34th AIAA Applied Aerodynamics Conference, AIAA Paper 2016-4326,
      Washington D.C., 2016.
 [18] J. A. Carnes and J. G. Coder.                  Effect of crossflow transition on the
      pressure-sensitive-paint rotor in hover. Journal of Aircraft, 59(1):29–46, 2022.
 [19] T. Cebeci and J. Cousteix. Modelling and Computation of Boundary-Layer Flows.
      Springer, 2005.
 [20] I. B. Celik, U. Ghia, P. J. Roache, C. J. Freitas, H. Coleman, and P. E. Raad. Procedure
      for estimation and reporting of uncertainty due to discretization in CFD applications.
      Journal of Fluids Engineering, 130(7), July 2008.
 [21] A. Chakraborty, A. Sinha Roy, and B. Dasgupta. Non-parametric smoothing for gra-
      dient methods in non-differentiable optimization problems. 2016 IEEE International
      Conference on Systems, Man, and Cybernetics, 2016.
 [22] T. Chau and D. W. Zingg. Aerodynamic design optimization of a transonic
      strut-braced-wing regional aircraft. Journal of Aircraft, 59(1):253–271, 2022.
 [23] O. Chernukhin and D. W. Zingg. Multimodality and global optimization in aerody-
      namic design. AIAA Journal, 51(6):1342–1354, 2013.
 [24] T. T. Chisholm and D. W. Zingg. A Jacobian-free Newton-Krylov algorithm for com-
      pressible turbulent fluid flows. Journal of Computational Physics, 228(9):3490–3507,
      2009.
 [25] J. H. Choi and O. J. Kwon. Enhancement of a correlation-based transition turbulence
      model for simulating crossflow instability. AIAA Journal, 53(10):3063–3072, 2015.
 [26] J. G. Coder. Enhancement of the amplification factor transport transition modeling
      framework. 55th AIAA Aerospace Sciences Meeting, AIAA Paper 2017-1709, 2017.
 [27] J. G. Coder. Standard test cases for transition model verification and validation in
      computational fluid dynamics. 2018 AIAA Aerospace Sciences Meeting, AIAA Paper
      2018-0029, 2018.
 [28] J. G. Coder. Further development of the amplification factor transport transition
      model for aerodynamic flows. 2019 AIAA Aerospace Sciences Meeting, AIAA Paper
      2019-0039, 2019.

<!-- ===== PDF page 123 ===== -->

REFERENCES                                                                                103

 [29] J. G. Coder and M. D. Maughmer. Computational fluid dynamics compatible transition
      modeling using an amplification factor transport equation. AIAA Journal, 52(11):2506–
      2512, 2014.
 [30] J. G. Coder and M. D. Maughmer. A CFD-compatible transition model using an am-
      plification factor transport equation. 51st AIAA Aerospace Sciences Meeting including
      the New Horizons Forum and Aerospace Exposition, AIAA Paper 2013-0253, 2013.
 [31] J. G. Coder, T. H. Pulliam, and J. C. Jensen. Contributions to HiLiftPW-3 using
      structured, overset grid methods. 2018 AIAA Aerospace Sciences Meeting, AIAA Paper
      2018-1039, 2018.
 [32] M. Costantini, U. Henne, C. Klein, and M. Miozzi. Skin-friction-based identification of
      the critical lines in a transonic, high Reynolds number flow via temperature-sensitive
      paint. Sensors, 21(15), 2021.
 [33] M. Costantini, U. Henne, S. Risius, and C. Klein. A robust method for reliable transi-
      tion detection in temperature-sensitive paint data. Aerospace Science and Technology,
      113:106702, 2021.
 [34] M. Costantini, T. Lee, T. Nonomura, K. Asai, and C. Klein. Feasibility of skin-friction
      field measurements in a transonic wind tunnel using a global luminescent oil film.
      Experiments in Fluids, 62(1):21, 2021.
 [35] J. Crouch, M. Sutanto, D. Witkowski, A. Watkins, M. Rivers, and R. Campbell. As-
      sessment of the national transonic facility for natural laminar flow testing. 48th AIAA
      Aerospace Sciences Meeting including the New Horizons Forum and Aerospace Expo-
      sition, 2010.
 [36] J. D. Crouch. Transition prediction and control for airplane applications. 28th Fluid
      Dynamics Conference, AIAA Paper 1997-1907, 1997.
 [37] J. D. Crouch. Boundary-layer transition prediction for laminar flow control. 45th AIAA
      Fluid Dynamics Conference, AIAA Paper 2015-2472, Dallas, TX, 2015.
 [38] J. Dagenhart and W. Saric. Crossflow stability and transition experiments in
      swept-wing flow. NASA Langley Technical Report Server, NASA/TP-1999-209344,
      1999.
 [39] J. R. Dagenhart. Amplified crossflow disturbances in the laminar boundary layer on
      swept wings with suction. NASA Technical Paper (1902), 1981.
 [40] M. Denison and T. H. Pulliam. Implementation and assessment of the amplification
      factor transport laminar-turbulent transition model. 2018 AIAA Aerospace Sciences
      Meeting, AIAA Paper 2018-3382, 2018.
 [41] B. Diskin, J. L. Thomas, C. L. Rumsey, and A. Schwöppe. Grid convergence for
      turbulent flows (invited). 53rd AIAA Aerospace Sciences Meeting, AIAA 2015-1746,
      2015.
 [42] S. S. Dodbele. Design optimization of natural laminar flow bodies in compressible flow.
      Journal of aircraft, 29(3):343–347, 1992.

<!-- ===== PDF page 124 ===== -->

104                                                                                REFERENCES

 [43] M. Drela. Design and optimization method for multi-element airfoils. Aerospace Design
      Conference, 1993.
 [44] M. Drela and M. B. Giles. Viscous-inviscid analysis of transonic and low Reynolds
      number airfoils. AIAA Journal, 25(10):1347–1355, 1987.
 [45] D. A. Dress, C. B. Johnson, P. D. McGuire, E. Stanewsky, and E. J. Ray. High
      Reynolds number tests of the CAST 10-2/DOA 2 airfoil in the Langley 0.3-meter
      transonic cryogenic tunnel, phase I. NASA, TM-84620, Hampton, 1983.
 [46] J. Driver and D. W. Zingg. Numerical aerodynamic optimization incorporating
      laminar-turbulent transition prediction. AIAA Journal, 45(8):1810–1818, 2007.
 [47] H. L. Dryden. Recent advances in the mechanics of boundary layer flow. Advances in
      Appl. Mech. I, pages 1–40, ed. R. von Mises and Th. von Karman, New York, 1948.
 [48] L. Eça and M. Hoekstra. Evaluation of numerical error estimation based on grid re-
      finement studies with the method of the manufactured solutions. Computers & Fluids,
      38(8):1580–1591, 2009.
 [49] Z. Fan, X. Yu, N. Qin, and M. Zhu. Benefit assessment of low-sweep transonic natural
      laminar flow wing for commercial aircraft. Journal of Aircraft, 58(6):1294–1301, 2021.
 [50] M. Fehrs. One-equation transition model for airfoil and wing aerodynamics. In: Dill-
      mann A. et al. (eds) New Results in Numerical and Experimental Fluid Mechanics
      XI. Notes on Numerical Fluid Mechanics and Multidisciplinary Design, 136:199–208,
      Springer, Cham. 2018.
 [51] M. Fehrs, S. Helm, and C. Kaiser. Numerical investigation of unsteady transitional
      boundary layer flows. IFASD 2019 -International Forum on Aeroelasticity and Struc-
      tural Dynamics, 2019.
 [52] M. Fehrs, A. C. van Rooij, and J. Nitzsche. Influence of boundary layer transition on
      the flutter behavior of a supercritical airfoil. CEAS Aeronautical Journal, 6(2):291–303,
      2015.
 [53] A. E. Forest. Engineering predictions of transitional boundary-layers. AGARD CP
      224, 1977.
 [54] M. Fujino, Y. Yoshizaki, and Y. Kawamura. Natural-laminar-flow airfoil development
      for a lightweight business jet. Journal of Aircraft, 40(4):609–615, 2003.
 [55] H. Gagnon and D. W. Zingg. Two-level free-form and axial deformation for exploratory
      aerodynamic shape optimization. AIAA Journal, 53(7):2015–2026, 2015.
 [56] M. Germano, U. Piomelli, P. Moin, and W. H. Cabot. A dynamic subgridscale eddy
      viscosity model. Physics of Fluids A: Fluid Dynamics, 3(7), 1991.
 [57] P. E. Gill, W. Murray, and M. A. Saunders. SNOPT: an SQP algorithm for large-scale
      constrained optimization. SIAM Journal on Optimization, 12(4):979–1006, 2002.
 [58] C. Grabe and A. Krumbein.                   Correlation-based transition modeling for
      three-dimensional aerodynamic configurations. Journal of Aircraft, 50(5):1533–1539,
      2013.

<!-- ===== PDF page 125 ===== -->

REFERENCES                                                                                 105

 [59] C. Grabe, N. Shengyang, and A. Krumbein. Transition transport modeling for the
      prediction of crossflow transition. AIAA Journal, 56(12):3167–3178, 2018.
 [60] C. Grabe, N. Shengyang, and A. Krumbein. Transition transport modeling for the
      prediction of crossflow transition. 34th AIAA Applied Aerodynamics Conference, AIAA
      Paper 2016-3572, Washington D.C., June, 2016.
 [61] P. S. Granville. The calculation of the viscous drag of bodies of revolution. Technical
      report, David Taylor Model Basin Rept. 849, 1953.
 [62] B. E. Green, J. L. Whitesides, R. L. Campbell, and R. E. Mineck. Method for the
      constrained design of natural laminar flow airfoils. Journal of Aircraft, 34(6):706–712,
      1997.
 [63] J. E. Green. Civil aviation and the environment: the next frontier for the aerodynam-
      icist. Aeronautical Journal, 110(1110):469–486, 2006.
 [64] G. L. O. Halila, J. R. R. A. Martins, and K. J. Fidkowski. Adjoint-based aerodynamic
      shape optimization including transition to turbulence effects. Aerospace Science and
      Technology, 107:1–15, 2020.
 [65] Z.-H. Han, J. Chen, K.-S. Zhang, Z.-M. Xu, Z. Zhu, and W.-P. Song. Aerody-
      namic shape optimization of natural-laminar-flow wing using surrogate-based ap-
      proach. AIAA Journal, 56(7):2579–2593, 2018.
 [66] A. Hebler, L. Schojda, and H. Mai. Experimental investigation of the aeroelastic
      behavior of a laminar airfoil in transonic flow. Proc. IFASD 2013, IFASD, Bristol,
      2013.
 [67] S. Helm, M. Fehrs, N. Krimmelbein, and A. Krumbein. Transition prediction and
      analysis of the CRM-NLF wing with the DLR TAU code. In: Dillmann A., Heller
      G., Krämer E., Wagner C. (eds) New Results in Numerical and Experimental Fluid
      Mechanics XIII. STAB/DGLR Symposium 2020. Notes on Numerical Fluid Mechanics
      and Multidisciplinary Design, 151, Springer, Cham., 2021.
 [68] J. E. Hicken. Efficient algorithms for future aircraft design: contributions to aerody-
      namic shape optimization. PhD thesis, Graduate Department of Aerospace Science
      and Engineering, University of Toronto, 2009.
 [69] J. E. Hicken and D. W. Zingg. Parallel Newton-Krylov solver for the Euler equations
      discretized using simultaneous approximation terms. AIAA Journal, 46(11):2773–2786,
      2008.
 [70] J. E. Hicken and D. W. Zingg. Aerodynamic optimization algorithm with integrated
      geometry parameterization and mesh movement. AIAA Journal, 48(2):400–413, 2010.
 [71] J. E. Hicken and D. W. Zingg. A simplified and flexible variant of GCROT for solving
      nonsymmetric linear systems. SIAM Journal on Scientific Computing, 32(3):1672–
      1694, 2010.
 [72] E. H. Hirschel, J. Cousteix, and W. Kordulla. Three-Dimensional Attached Viscous
      Flow. Springer, 2014.

<!-- ===== PDF page 126 ===== -->

106                                                                                 REFERENCES

 [73] J. Hollom and N. Qin. Uncertainty analysis and robust shape optimisation for laminar
      flow aerofoils. The Aeronautical Journal, 125(1284):365–388, 2021.
 [74] F. Ilinca and D. Pelletier. Positivity preservation and adaptive solution of two-equation
      models of turbulence. International Journal of Thermal Sciences, 38(7):560–571, 1999.
 [75] A. Jameson, W. Schmidt, and E. Turkel. Numerical solution of the Euler equations
      by finite volume methods using Runge-Kutta time-stepping schemes. 14th Fluid and
      Plasma Dynamics Conference, AIAA Paper 81-1259, 1981.
 [76] Y. S. Jung and J. D. Baeder. γ-Re        ˜ θt -SA with crossflow transition model using
      Hamiltonian-Strand approach. 2018 AIAA Aerospace Sciences Meeting, AIAA Paper
      2018-1040, 2018.
 [77] H. Kaya and I. H. Tuncer. Discrete adjoint-based aerodynamic shape optimization
      framework for natural laminar flows. AIAA Journal, 2021.
 [78] G. J. Kennedy and J. E. Hicken. Improved constraint-aggregation methods. Computer
      Methods in Applied Mechanics and Engineering, 289:332–354, 2015.
 [79] G. K. W. Kenway and J. R. R. A. Martins. Multipoint high-fidelity aerostructural
      optimization of a transport aircraft configuration. Journal of Aircraft, 51(1):144–160,
      2014.
 [80] P. Khayatzadeh and S. K. Nadarajah. Laminar-turbulent flow simulation for wind
                                    ˜ θt transition model. Wind Energy, 17(6):901–918, 2014.
      turbine profiles using the γ-Re
 [81] P. Khayatzadeh and S. K. Nadarajah. Aerodynamic shape optimization via discrete
      viscous adjoint equations for the k-ω SST turbulence and γ-Re  ˜ θt transition models. 49th
      AIAA Aerospace Sciences Meeting including the New Horizons Forum and Aerospace
      Exposition, AIAA Paper 2011-1247, Orlando, Florida, January, 2011.
 [82] P. Khayatzadeh and S. K. Nadarajah. Aerodynamic shape optimization of natural
      laminar flow (NLF) airfoils. 50th AIAA Aerospace Sciences Meeting including the
      New Horizons Forum and Aerospace Exposition, AIAA Paper 2012-0061, Nashville,
      Tennessee, January, 2012.
 [83] D. Knoll and D. Keyes. Jacobian-free Newton-Krylov methods: a survey of approaches
      and applications. Journal of Computational Physics, 193(2):357–397, 2004.
 [84] G. Kreisselmeier and R. Steinhauser. Systematic control design by optimizing a vec-
      tor performance index. In M. CUENOD, editor, Computer Aided Design of Control
      Systems, pages 113–117. Pergamon, 1980.
 [85] N. Krimmelbein and A. Krumbein. Determination of critical N-factors for the
      CRM-NLF wing. In: Dillmann A., Heller G., Krämer E., Wagner C. (eds) New
      Results in Numerical and Experimental Fluid Mechanics XIII. STAB/DGLR Sympo-
      sium 2020. Notes on Numerical Fluid Mechanics and Multidisciplinary Design, 151,
      Springer, Cham., 2021.
 [86] N. Krimmelbein and R. Radespiel. Transition prediction for three-dimensional flows
      using parallel computation. Computers and Fluids Journal, 38(1):121–136, 2009.

<!-- ===== PDF page 127 ===== -->

REFERENCES                                                                                   107

 [87] I. Kroo and P. Sturdza. Design-oriented aerodynamic analysis for supersonic laminar
      flow wings. 41st Aerospace Sciences Meeting and Exhibit, 2003.
 [88] A. Krumbein, N. Krimmelbein, and G. Schrauf. Automatic transition prediction in hy-
      brid flow solver, part 1: Methodology and sensitivities. Journal of Aircraft, 46(4):1176–
      1190, 2009.
 [89] A. Krumbein, N. Krimmelbein, and G. Schrauf. Automatic transition prediction in
      hybrid flow solver, part 2: Practical application. Journal of Aircraft, 46(4):1191–1199,
      2009.
 [90] M. Kruse, F. Munoz, and R. Radespiel. Transition prediction results for sickle wing
      and NLF(1)-0416 test cases. 2018 AIAA Aerospace Sciences Meeting, AIAA Paper
      2018-0537, 2018.
 [91] R. B. Langtry. A correlation-based transition model using local variables for unstruc-
      tured parallelized CFD codes. PhD thesis, DLR, 2006.
 [92] R. B. Langtry and F. R. Menter. Correlation-based transition modeling for unstruc-
      tured parallelized computational fluid dynamics codes. AIAA Journal, 47(12):2894–
      2906, 2009.
 [93] R. B. Langtry, K. Sengupta, D. T. Yeh, and A. J. Dorgan. Extending the γ-Re           ˜ θt
      correlation based transition model for crossflow effects. 45th AIAA Fluid Dynamics
      Conference, AIAA Paper 2015-2474, 2015.
 [94] A. Lee, T. Reist, and D. W. Zingg. Further exploration of regional-class hybrid
      wing-body aircraft through multifidelity optimization. AIAA Scitech 2021 Forum,
      AIAA 2021-0014, 2021.
 [95] C. Lee, D. Koo, and D. W. Zingg. Comparison of B-spline surface and free-form
      deformation geometry control methods for aerodynamic shape optimization. AIAA
      Journal, 55(1):228–240, 2017.
 [96] J. Lee and A. Jameson. Natural-laminar-flow airfoil and wing design by adjoint method
      and automatic transition prediction. 47th AIAA Aerospace Sciences Meeting and Ex-
      hibit, AIAA Paper 2009-897, Orlando, Florida, January, 2009.
 [97] L. Lees and C. C. Lin. Investigation of the stability of the laminar boundary layer in
      a compressible fluid. National Advisory Committee for Aeronautics, NACA TN No.
      1115, 1946.
 [98] X. Li. An entropy-based aggregate method for minimax optimization. Engineering
      Optimization, 18(4):277–285, 1997.
 [99] C. Lian, G. Xia, and C. L. Merkle. Solution-limited time stepping to enhance reliability
      in CFD applications. Journal of Computational Physics, 228(13):4836–4857, 2009.
[100] C. Lian, G. Xia, and C. L. Merkle. Impact of source terms on reliability of CFD
      algorithms. Computers and Fluids Journal, 39:1909–1922, 2010.
[101] R. H. Liebeck. A class of airfoils designed for high lift in incompressible flow. Journal
      of Aircraft, 10(10):610–617, 1973.

<!-- ===== PDF page 128 ===== -->

108                                                                               REFERENCES

[102] R. Lopes, L. Eça, and G. Vaz. On the numerical behavior of RANS-based transition
      models. Journal of Fluids Engineering, 142(5), 2020.
[103] R. Lopes, L. Eça, G. Vaz, and M. Kerkvliet. Assessing numerical aspects of transitional
      flow simulations using the RANS equations. International Journal of Computational
      Fluid Dynamics, Published Online, January 14, 2021.
[104] M. N. Lynde, R. L. Campbell, M. B. Rivers, S. A. Viken, D. T. Chan, N. A. Watkins,
      and S. L. Goodliff. Preliminary results from an experimental assessment of a natural
      laminar flow design method. 2019 AIAA Aerospace Sciences Meeting, AIAA Paper
      2019-2298, 2019.
[105] M. N. Lynde, R. L. Campbell, and S. A. Viken. Additional findings from the common
      research model natural laminar flow wind tunnel test. AIAA Aviation 2019 Forum,
      AIAA Paper 2019-3292, 2019.
[106] L. M. Mack. Linear stability theory and the problem of supersonic boundary-layer
      transition. AIAA Journal, 13(3):278–289, 1975.
[107] L. M. Mack. Transition and laminar instability. Jet Propulsion Laboratory Publication,
      NASA-CP-153203, 1977.
[108] M. R. Malik, P. Balakumar, and C. Chang. Linear stability of hypersonic boundary
      layers. 10th National Aero-Space Plane Symposium, (189), 1991.
[109] M. R. Malik, J. D. Crouch, W. S. Saric, J. C. Lin, and E. A. Whalen. Application
      of drag reduction techniques to transport aircraft. In Encyclopedia of Aerospace Engi-
      neering (eds R. Blockley and W. Shyy), 2015.
[110] J. R. R. A. Martins and N. M. K. Poon. On structural optimization using constraint
      aggregation. In VI World Congress on Structural and Multidisciplinary Optimization
      WCSMO6, Rio de Janeiro, Brasil. Citeseer, 2005.
[111] D. J. Mavriplis, Z. Yang, and E. M. Anderson. Adjoint based optimization of a slotted
      natural laminar flow wing for ultra efficient flight. AIAA Scitech 2020 Forum, 2020.
[112] E. Mayda. Boundary Layer Transition Prediction for Reynolds-averaged Navier-Stokes
      Methods. PhD thesis, University of California, Davis, 2007.
[113] R. E. Mayle. The role of laminar-turbulent transition in gas turbine engines. Journal
      of Turbomachinery, 113(4):509–536, 1991.
[114] R. E. Mayle, K. Dullenkopf, and A. Schulz. The turbulence that matters. ASME 1997
      International Gas Turbine and Aeroengine Congress and Exhibition, 1, 1997.
[115] S. Medida. Correlation-based Transition Modeling for External Aerodynamic Flows.
      PhD thesis, University of Maryland, College Park, 2014.
[116] S. Medida and J. D. Baeder. Application of the correlation-based γ-Re    ˜ θt transition
      model to the Spalart-Allmaras turbulence model. 20th AIAA Computational Fluid
      Dynamics Conference, AIAA Paper 2011-3979, Honolulu, HI, 2011.
[117] F. R. Menter. Two-equation eddy-viscosity turbulence models for engineering applica-
      tions. AIAA Journal, 32(8):1598–1605, 1994.

<!-- ===== PDF page 129 ===== -->

REFERENCES                                                                                109

[118] F. R. Menter, P. E. Smirnov, T. Liu, and R. Avancha. A one-equation local
      correlation-based transition model. Flow, Turbulence and Combustion, 95:583–619,
      2015.
[119] V. Michelassi, J. G. Wissink, J. Frohlich, and W. Rodi. Large-eddy simulation of flow
      around low-pressure turbine blade with incoming wakes. AIAA Journal, 41(11):2143–
      2156, 2003.
[120] J. M. Modisette. An automated reliable method for two-dimensional Reynolds-averaged
      Navier-Stokes simulations. PhD thesis, Massachusetts Institute of Technology, 2011.
[121] F. Moens, J. Perraud, A. Krumbein, T. Toulorge, P. Iannelli, and A. Hanifi. Transition
      prediction and impact on a three-dimensional high-lift-wing configuration. Journal of
      Aircraft, 45(5):1751–1766, 2008.
[122] Y. Moryossef and Y. Levy. Unconditionally positive implicit procedure for two-equation
      turbulence models: application to k-ω turbulence models. Journal of Computational
      Physics, 220:88–108, 2006.
[123] A. Mosahebi and E. Laurendeau. Convergence characteristics of fully and loosely
      coupled numerical approaches for transition models. AIAA Journal, 53(5):1399–1404,
      2015.
[124] A. Mosahebi and E. Laurendeau. Introduction of a modified segregated numerical
                                              ˜ θt transition model. International Journal of
      approach for efficient simulation of γ-Re
      Computational Fluid Dynamics, 29:357–375, 2015.
[125] C. Müller and F. Herbst. Modelling of crossflow-induced transition based on local
      variables. 6th European Conference on Computational Fluid Dynamics (ECFD), 2014.
[126] S. Nie, N. Krimmelbein, A. Krumbein, and C. Grabe. Extension of a Reynolds-
      stress-based transition transport model for crossflow transition. Journal of Aircraft,
      55(4):1641–1654, 2018.
[127] L. Osusky. A Numerical Methodology for Aerodynamic Shape Optimization in Turbu-
      lent Flow Enabling Large Geometric Variation. PhD thesis, Graduate Department of
      Aerospace Science and Engineering, University of Toronto, 2014.
[128] L. Osusky, H. P. Buckley, T. A. Reist, and D. W. Zingg. Drag minimization based on the
      Navier-Stokes equations using a Newton-Krylov approach. AIAA Journal, 53(6):1555–
      1577, 2015.
[129] M. Osusky. A Parallel Newton-Krylov-Schur Algorithm for the Reynolds-averaged
      Navier-Stokes Equations. PhD thesis, Graduate Department of Aerospace Science and
      Engineering, University of Toronto, 2013.
[130] M. Osusky, P. D. Boom, and D. W. Zingg. Results from the Fifth AIAA Drag Prediction
      Workshop obtained with a parallel Newton-Krylov-Schur flow solver discretized using
      summation-by-parts operators. 31st AIAA Applied Aerodynamics Conference, AIAA
      2013-2511, 2013.
[131] M. Osusky and D. W. Zingg. Parallel Newton-Krylov-Schur flow solver for the

<!-- ===== PDF page 130 ===== -->

110                                                                             REFERENCES

      Navier-Stokes equations. AIAA Journal, 51(12):2833–2851, 2013.
[132] P. Paredes, B. S. Venkatachari, M. Choudhari, F. Li, N. Hildebrand, and C. L. Chang.
      Transition analysis for the CRM-NLF wind tunnel configuration. 2021 AIAA Aerospace
      Sciences Meeting, AIAA Paper 2021-1431, 2021.
[133] L. Pascal, G. Delattre, H. Deniau, and J. Cliquet. Stability-based transition model
      using transport equations. AIAA Journal, 58(7):2933–2942, 2020.
[134] D. Pasquale, A. Rona, and S. J. Garrett. A selective review of CFD transition models.
      39th AIAA Fluid Dynamics Conference, AIAA Paper 2009-3812, San Antonio, Texas,
      2009.
[135] E. Y. Pee. On solving large-scale finite minimax problems using exponential smoothing.
      Journal of Optimization Theory and Applications, 148(2):390–421, 2011.
[136] J. Perraud, H. Deniau, and G. Casalis. Overview of transition prediction tools in the
      elsA software. ECCOMAS 2014, 2014.
[137] J. Perraud and A. Durant. Stability-based Mach zero to four longitudinal transition
      prediction criterion. Journal of Spacecraft and Rockets, 53(4):730–742, 2016.
[138] R. Petzold and R. Radespiel. Transition on a wing with spanwise varying crossflow
      and linear stability analysis. AIAA Journal, 53(2):321–335, 2015.
[139] M. G. H. Piotrowski and D. W. Zingg. Smooth local correlation-based transition model
      for the Spalart-Allmaras turbulence model. AIAA Journal, 59(2):474–492, 2021.
[140] M. G. H. Piotrowski and D. W. Zingg. Numerical behaviour of a smooth local
      correlation-based transition model in a Newton-Krylov flow solver. AIAA Scitech 2022
      Forum, AIAA 2022-0909, 2022.
[141] M. G. H. Piotrowski and D. W. Zingg. Investigation of a smooth local correlation-based
      transition model in a discrete-adjoint aerodynamic shape optimization algorithm.
      AIAA Scitech 2022 Forum, AIAA 2022-1865, 2022.
[142] M. G. H. Piotrowski and D. W. Zingg. Compressibility corrections to extend a smooth
      local correlation-based transition model to transonic flows. In Preparation.
[143] E. Polak, J. O. Royset, and R. S. Womersley. Algorithms with adaptive smooth-
      ing for finite minimax problems. Journal of Optimization Theory and Applications,
      119(3):459–484, 2003.
[144] D. Poll. Transition in the infinite swept attachment line boundary layer. Aeronautical
      Quarterly, 30(4):607–629, 1979.
[145] J. O. Pralits. Optimal design of natural and hybrid laminar flow control on wings. PhD
      thesis, Department of Mechanics, Royal Institute of Technology, Stockholm, Sweden,
      2003.
[146] T. H. Pulliam. Efficient solution methods for the Navier-Stokes equations. Tech. rep.,
      Lecture Notes for the von Karman Inst. for Fluid Dynamics Lecture Series: Numerical
      Techniques for Viscous Flow Computation in Turbomachinery Bladings, 1986.
[147] T. H. Pulliam and D. W. Zingg. Fundamental algorithms in computational fluid dy-

<!-- ===== PDF page 131 ===== -->

REFERENCES                                                                                    111

      namics. Springer International Publishing, 2014.
[148] R. H. Radeztsky, M. S. Reibert, and W. S. Saric. Effect of micron-sized roughness on
      transition in swept-wing flows. 31st Aerospace Sciences Meeting and Exhibit, AIAA
      Paper 93-0076, 1993.
[149] R. Rashad. High-Fidelity Aerodynamic Shape Optimization for Natural Laminar Flow.
      PhD thesis, University of Toronto, Toronto, Ontario, Canada, 2016.
[150] R. Rashad and D. W. Zingg. Aerodynamic shape optimization for natural laminar flow
      using a discrete-adjoint approach. AIAA Journal, 54(11):3321–3337, 2016.
[151] H. L. Reed and W. S. Saric. Transition mechanisms for transport aircraft. 38th Fluid
      Dynamics Conference and Exhibit, AIAA Paper 2008-3743, Seattle, Washington, 2008.
[152] T. A. Reist, D. Koo, D. W. Zingg, P. Bochud, P. Castonguay, and D. Leblond. Cross
      validation of aerodynamic shape optimization methodologies for aircraft wing-body
      optimization. AIAA Journal, 58(6):2581–2595, 2020.
[153] T. A. Reist and D. W. Zingg. High-fidelity aerodynamic shape optimization of a
      lifting-fuselage concept for regional aircraft. Journal of Aircraft, 54(3):1085–1097, 2016.
[154] T. A. Reist, D. W. Zingg, M. Rakowitz, G. Potter, and S. Banerjee. Multifidelity opti-
      mization of hybrid wing-body aircraft with stability and control requirements. Journal
      of Aircraft, 56(2):442–456, 2018.
[155] F. J. Richards. A flexible growth function for empirical use. Journal of Experimental
      Botany, 10(2):290–301, 1959.
[156] S. Risuis, M. Costantini, S. Koch, S. Hein, and C. Klein. Unit Reynolds num-
      ber, Mach number and pressure gradient effects on laminar-turbulent transition in
      two-dimensional boundary layers. Experiments in Fluids, 59(5), 2018.
[157] R. L. Rivest. Game tree searching by min/max approximation. Artificial Intelligence,
      34(1):77–96, 1988.
[158] P. J. Roache. Perspective: A method for uniform reporting of grid refinement studies.
      Journal of Fluids Engineering, 116(3):405–413, 1994.
[159] M. Robitaille, A. Mosahebi, and E. Laurendeau. Design of adaptive transonic laminar
                           ˜ θt transition model. Aerospace Science and Technology, 46:60–71,
      airfoils using the γ-Re
      2015.
[160] Y. Saad and M. H. Schultz. GMRES: A generalized minimal residual algorithm for
      solving nonsymmetric linear systems. SIAM Journal on Scientific and Statistical Com-
      puting, 7(3):856–869, 1986.
[161] Y. Saad and M. Sosonkina. Distributed Schur complement techniques for general sparse
      linear systems. SIAM Journal of Scientific Computing, 21:1337–1357, 1999.
[162] C. Sabater, P. Bekemeyer, and S. Görtz. Robust design of transonic natural laminar
      flow wings under environmental and operational uncertainties. AIAA Scitech 2021
      Forum, 2021.
[163] G. B. Schubauer and H. K. Skramstad. Laminar-boundary-layer oscillations and tran-

<!-- ===== PDF page 132 ===== -->

112                                                                                REFERENCES

      sition on a flat plate. NACA TN 909, 1948.
[164] J. Schucker.          Development of a three-equation γ-Re         ˜ θt -Spalart-Allmaras
      turbulence-transition model. Master’s thesis, DLR, 2012.
[165] C. Seyfert and A. Krumbein. Evaluation of a correlation-based transition model and
      comparison with the eN method. Journal of Aircraft, 49(6):1765–1773, 2012.
[166] Y. Shi, R. Gross, C. A. Mader, and J. Martins. Transition prediction based on linear
      stability theory with the RANS solver for three-dimensional configurations. 2018 AIAA
      Aerospace Sciences Meeting, AIAA Paper 2018-0819, 2018.
[167] Y. Shi, C. A. Mader, S. He, G. L. O. Halila, and J. R. R. A. Martins. Natural
      laminar-flow airfoil optimization design using a discrete adjoint approach. AIAA Jour-
      nal, 58(11):4702–4722, 2020.
[168] Y. Shi, C. A. Mader, and J. R. R. A. Martins. Natural laminar flow wing optimiza-
      tion using a discrete adjoint approach. Structural and Multidisciplinary Optimization,
      64(2):541–562, 2021.
[169] J. Smagorinsky. General circulation experiments with the primitive equations: I. the
      basic experiment. General Circulation Research Laboratory, U.S. Weather Bureau,
      Washington, D.C, 1963.
[170] A. M. O. Smith and N. Gamberoni. Transition, pressure gradient, and stability theory.
      Douglas Aircraft Report ES-26388, 1956.
[171] D. M. Somers. Design and experimental results for a natural-laminar-flow airfoil for
      general aviation applications. NASA TP-1861, 1981.
[172] D. M. Somers. Design and experimental results for the S809 airfoil. National Renewable
      Energy Laboratory, Golden, Colorado, NRELSR-440-6918, 1989.
[173] P. R. Spalart. Strategies for turbulence modelling and simulation. International Journal
      of Heat and Fluid Flow, 21(3):252–263, 2000.
[174] P. R. Spalart and S. R. Allmaras. A one-equation turbulence model for aerodynamic
      flows. 30th AIAA Aerospace Sciences Meeting and Exhibit, AIAA Paper 092-0439,
      Reno, Nevada, United States 1992.
[175] E. Stanewsky and H. Zimmer. Development and wind tunnel investigation of three
      supercritical airfoil profiles for transport aircraft. NASA, TM-75840, Washington,
      D.C., 1980.
[176] D. Stefanski, R. Glasby, J. T. Erwin, and J. G. Coder. Development of a predictive
      capability for laminar-turbulent transition in HPCMP CREATETM-AV kestrel com-
      ponent COFFE using the amplification factor transport model. 2018 AIAA Aerospace
      Sciences Meeting, AIAA Paper 2018-1041, 2018.
[177] H. W. Stock. eN transition prediction in three-dimensional boundary layers on inclined
      prolate spheroids. AIAA Journal, 44(1):108–118, 2006.
[178] T. Streit, H. Horstmann, G. Shraut, S. Hein, U. Fey, Y. Egmai, J. Perraud, I. Salah
      El Din, U. Cella, and J. Quest. Complementary numerical and experimental data

<!-- ===== PDF page 133 ===== -->

REFERENCES                                                                                  113

      analysis of the ETW telfona pathfinder wing transition tests. 49th Aerospace Sciences
      Meeting and Exhibit, AIAA Paper 2011-881, Orlando, Florida, January, 2011.
[179] G. M. Streuber and D. W. Zingg. Evaluating the risk of local optima in aerodynamic
      shape optimization. AIAA Journal, 59(1):75–87, 2021.
[180] P. Ströer, N. Krimmelbein, A. Krumbein, and C. Grabe. Stability-based transition
      transport modeling for unstructured computational fluid dynamics including convec-
      tion effects. AIAA Journal, 58(4):1506–1517, 2020.
[181] P. Ströer, N. Krimmelbein, A. Krumbein, and C. Grabe. Stability-based transition
      transport modeling for unstructured computational fluid dynamics at transonic flow
      conditions. AIAA Journal, 59(9):3585–3597, 2021.
[182] K. Suluksna, P. Dechaumphai, and E. Juntasaro. Correlations for modeling transitional
      boundary layers under influences of freestream turbulence and pressure gradient. In-
      ternational Journal of Heat and Fluid Flow, 30(1):66–75, 2009.
[183] A. Suzen, Y. Bora. Predictions of separated and transitional boundary layers un-
      der low-pressure turbine airfoil conditions using an intermittency transport equation.
      Journal of Turbomachinery, 125(3):455–464, 2003.
[184] R. C. Swanson and E. Turkel. On central-difference and upwind schemes. Journal of
      Computational Physics, 101:292–306, 1992.
[185] Z. Tang, Y. Chen, and L. Zhang. Natural laminar flow shape optimization in transonic
      regime with competitive Nash game strategy. Applied Mathematical Modelling, 48:534–
      547, 2017.
[186] E. R. Van Driest. Calculation of the stability of the laminar boundary layer in a com-
      pressible fluid on a flat plate with heat transfer. Journal of the Aeronautical Sciences,
      19(12):801–812, 1952.
[187] J. L. Van Ingen. A suggested semi-empirical method for the calculation of the bound-
      ary layer transition region. Technische Hogeschool Delft, Vliegtuigbouwkunde, Rapport
      VTH-74, 1956.
[188] B. S. Venkatachari, P. Paredes, J. M. Derlaga, P. Buning, M. Choudhari, F. Li, and
      C. L. Chang. Assessment of transition modeling capability in overflow with emphasis
      on swept-wing configurations. 2020 AIAA Aerospace Sciences Meeting, AIAA Paper
      2020-1034, 2020.
[189] B. S. Venkatachari, P. Paredes, J. M. Derlaga, P. Buning, M. Choudhari, F. Li, and
      C. L. Chang. Assessment of RANS-based transition models based on experimental
      data of the common research model with natural laminar flow. 2021 AIAA Aerospace
      Sciences Meeting, AIAA Paper 2021-1430, 2021.
[190] H. K. Versteeg and W. Malalasekera. An introduction to computational fluid dynamics:
      the finite volume method. Pearson Education, 2007.
[191] F. M. White. Viscous Fluid Flow. McGraw-Hill, Inc., third edition, 2005.
[192] X. Wu and P. A. Durbin. DNS of fully turbulent flow in a LPT passage. International

<!-- ===== PDF page 134 ===== -->

114                                                                                 REFERENCES

      Journal of Heat and Fluid Flow, 24(4):636–644, 2003.
[193] J. Xu, X. Han, L. Qiao, J. Bai, and Y. Zhang. Fully local amplification factor transport
      equation for stationary crossflow instabilities. AIAA Journal, 57(7):2682–2693, 2019.
[194] J. Xu, L. Qiao, and J. Bai. Improved local amplification factor transport equation
      for stationary crossflow instability in subsonic and transonic flows. Chinese Journal of
      Aeronautics, 33(12):3073–3081, 2020.
[195] J. K. Xu, J. Q. Bai, L. Qiao, and Y. Zhang. Correlation-based transition transport
      modeling for simulating crossflow instabilities. Journal of Applied Fluid Mechanics,
      9(5):2435–2442, 2016.
[196] Z. Yang and D. J. Mavriplis. Discrete adjoint formulation for turbulent flow problems
      with transition modelling on unstructured meshes. AIAA Scitech 2019 Forum, 2019.
[197] A. Yildirim, G. K. W. Kenway, C. A. Mader, and J. R. R. A. Martins. A Jacobian-free
      approximate Newton-Krylov startup strategy for RANS simulations. Journal of Com-
      putational Physics, 397, 2019.
[198] Y. Zhang, X. Fang, H. Chen, S. Fu, Z. Duan, and Y. Zhang. Supercritical natural
      laminar flow airfoil optimization for regional aircraft wing design. Aerospace Science
      and Technology, 43:152–164, 2015.
[199] K. Zhao, Z. Gao, and J. Huang. Robust design of natural laminar flow supercriti-
      cal airfoil by multi-objective evolution method. Applied Mathematics and Mechanics,
      35(2):191–202, 2014.
[200] X. Zheng, X. Liu, F. Liu, and C. Yang. Turbulent transition simulation using the k-ω
      model. International Journal for Numerical Methods in Engineering, 42:902–976, 1998.
[201] M. Zhu and N. Qin. Balancing laminar extension and wave drag for transonic swept
      wings. AIAA Journal, 59(5):1660–1672, 2021.
[202] D. W. Zingg. An approach to the design of airfoils with high lift to drag ratios. Master’s
      thesis, University of Toronto, Toronto, Ontario, Canada, 1981.
[203] D. W. Zingg. Grid studies for thin-layer Navier-Stokes computations of airfoil flowfields.
      AIAA Journal, 30(10):2561–2564, 1992.
[204] D. W. Zingg. Viscous airfoil computations using Richardson extrapolation. 10th Com-
      putational Fluid Dynamics Conference, AIAA 1991-1559, 1991.
[205] D. W. Zingg and P. Godin. A perspective on turbulence models for aerodynamic flows.
      International Journal of Computational Fluid Dynamics, 23(4):327–335, 2009.

<!-- ===== PDF page 135 ===== -->

Appendix A

Fonset,scf Source-Term Validation

In this chapter, the TU Braunschweig Sickle Wing [138] transition test case is revisited to
validate the new Fonset,scf source term (Equation 2.81). A thorough study of the Sickle Wing
using the SA-LM2015 and SA-sLM2015 transition models with the original Dscf source term
(Equation 2.15) is presented in Section 4.1.4 [139].

A.1     TU Braunschweig Sickle Wing

The TU Braunschweig Sickle Wing [138] was developed as a test case for crossflow transition
on a complex geometry. The three structured multi-block grids presented previously in
Table 4.3 are simulated using the SA-sLM2015 transition model with the two crossflow
source-term strategies, the original Dscf source term and the new Fonset,scf source term. The
skin friction coefficient profiles on both the upper and lower surfaces of the Sickle Wing are
overlaid with transition locations from the experiment [138] in Figure A.1. The transition
front on the lower surface of the wing, where transition is dominated by TS instabilities,
is consistent between the two crossflow source-term approaches. On the upper surface,
where stationary crossflow instabilities are the primary transition mechanism, the two source-
term approaches produce similar transition fronts. Both strategies do a reasonable job of
reproducing the transition front from the experiment.

<!-- ===== PDF page 136 ===== -->

116                                                   APPENDIX A. FONSET,SCF SOURCE-TERM VALIDATION

              coarse                            medium                              fine

                                         (a) Dscf : upper surface

                                       (b) Fonset,scf : upper surface

                                         (c) Dscf : lower surface

                                       (d) Fonset,scf : lower surface

Figure A.1: Sickle Wing grid-refinement study upper and lower surface skin friction coefficient profiles
produced at M = 0.16, α = −2.6◦ , T u = 0.17%, h = 9.78µm, and Re = 2.75 × 106 overlaid with transition
locations from the experiment [138].

<!-- ===== PDF page 137 ===== -->

Appendix B

Further Evaluation of the Iterative
and Grid Convergence of the Smooth
Transition Model

The development of best practices for solving transport-equation-based transition models
in implicit solvers facilitates the deep iterative convergence necessary for investigating their
grid convergence [158]. Eça et al. [48] demonstrated that iterative error must be reduced
two to three orders of magnitude below the discretization error to prevent contamination,
while the elimination of numerical error, including both iterative and discretization error, is
necessary in order to properly evaluate the modelling error of these RANS-based transition
models [205].

   The goal of this Appendix is to investigate in more detail the iterative and grid convergence
of the SA-sLM2015 smooth local correlation-based transition model presented in Section 2.2
with the compressibility corrections presented in Section 2.3 in the implicit Newton-Krylov
flow solver presented in Section 3.1. To achieve this goal, the numerical behaviour of free-
transition simulations of the subsonic NLF0416 airfoil (Section 4.1.1), transonic VA-2 airfoil
(Section 4.2.2), and NASA CRM-NLF geometry (Section 4.2.3) are evaluated relative to fully
turbulent simulations performed using the SA turbulence model. Other than the NLF0416
test case, the numerical behaviour is investigated using the grids previously presented for
each test case. For the free-transition simulations, fully, loosely coupled, and decoupled
linearization strategies for the turbulence and transition model equations are investigated
along with the source-term time step (STTS) restriction presented in Section 3.1.4 in order
to determine best practices for solving transport equations with large and highly non-linear
source terms in a strong implicit solver. Grid-refinement studies are performed to investigate
the grid convergence of the smooth transition model relative to fully turbulent simulations.

<!-- ===== PDF page 138 ===== -->

118                          APPENDIX B. FURTHER EVALUATION OF ITERATIVE AND GRID CONVERGENCE

      (a) fully coupled        (b) turb-trans coupled       (c) trans coupled            (d) decoupled

Figure B.1: Block-diagonal Jacobians produced by the four linearization strategies for the turbulence and
transition model equations with varying levels of coupling, with the grey elements representing filled Jacobian
entries.

B.1      Methodology

B.1.1     Linearization Strategies

Four linearization strategies for the turbulence and transition model equations with varying
levels of coupling are considered: a fully coupled approach, a loosely coupled approach with
the turbulence and transition model equations coupled, a loosely coupled approach with
the transition model equations coupled with an independent linearization for the turbulence
model, and a decoupled approach with an independent linearization for each of the turbu-
lence and transition model equations. The four linearization strategies are investigated each
with and without the source-term time step restriction presented in Section 3.1.4, with the
coupling of the source-term Jacobian used to evaluate the source-term eigenvalues mirroring
the solver linearization strategy. Other than the source-term time step restriction, no special
considerations are applied to the turbulence and transition model production and destruc-
tion terms. Loosely coupled and decoupled linearization strategies for the turbulence and
transition model equations are implemented by neglecting off-diagonal components of the
fully coupled 8x8 Jacobian. The block-diagonal Jacobians for the four linearization strategies
for a three-dimensional simulation are presented in Figure B.1.
   The source terms in the turbulence and transition models are most destabilizing during
the early stages of convergence. To recover the nonlinear convergence properties of Newton’s
method, the loosely coupled and decoupled linearization strategies are only adopted in the
approximate-Newton phase of the solver, with the mean-flow, turbulence, and transition
model equations fully coupled in the inexact-Newton phase. The emphasis of this study
is on the impact of decoupling the turbulence and transition model equations on nonlinear
convergence relative to a fully coupled linearization strategy, in particular the impact on
robustness. Although using a sequential segregated algorithm can reduce wall-clock times
for these loosely coupled and decoupled approaches in the approximate-Newton phase, the
current implementation isolates the effects of varying coupling strategies and avoids differ-
ences introduced by using a different solution algorithm, such as the introduction of different
preconditioners.

<!-- ===== PDF page 139 ===== -->

B.1. METHODOLOGY                                                                                119

B.1.2     Estimating Grid Convergence

The grid convergence of the fully turbulent and free-transition simulations is evaluated us-
ing Richardson extrapolation [158]. There are several important requirements for using
Richardson extrapolation to report grid-converged values; for example, the grids must be
sufficiently fine to be in the asymptotic range of convergence in order to have full confidence
in the results [158]. Previous work has demonstrated the difficulty in achieving asymptotic
convergence for simple test cases with fully turbulent RANS solutions [204, 203, 41]. The
presence of discontinuities, such as shocks, and singularities can complicate grid-convergence
studies [158]. In the current work, Richardson extrapolation is used to evaluate the relative
grid convergence of the fully turbulent and free-transition simulations. The grid-refinement
results are evaluated using the apparent order of convergence, p, the actual fractional error
(AFE), and the grid-convergence index (GCI), each calculated for the finest grid level. These
grid-convergence metrics are calculated according to ASME standards [20] with some minor
modifications5 as follows:
                                   ln |(f3 − f2 )/(f2 − f1 )| + q(p)
                              p=                                     ,                       (B.1)
                                                 ln(r21 )
                                   p                                   
                                    r21 − s                       f3 − f2
                        q(p) = ln p           , s = 1 · sign                ,                (B.2)
                                    r32 − s                       f2 − f1
                        fexact − f1                r p f1 − f2               f1 − f2
                  AFE =              , fexact = 21p            , EFE =               ,       (B.3)
                           fexact                    r21 − 1                    f1

where f1 , f3 , and fexact represent the values on the finest and coarsest of the three grid levels
being investigated, and the grid-converged value estimated using Richardson extrapolation,
respectively, r is the refinement factor, and EFE is the estimated fractional error.
   For the results presented, the total number of grid points is approximately doubled be-
tween refinement levels, with the refinement factors calculated using the total number of
grid nodes at each grid level: r21 = (N1 /N2 )1/2 and r21 = (N1 /N2 )1/3 for two- and three-
dimensional cases, respectively. Using a consistent family of grids is important for evaluating
grid convergence. This can be done by removing every second node in each coordinate di-
rection repeatedly from a fine mesh. However, it can be preferable to use a refinement factor
less than two. Grid sequences for the NLF0416 airfoil and CRM-NLF geometry are produced
using an automatic grid-refinement algorithm that uses an analytical representation of the
grid based on B-spline volumes [128]. The VA-2 airfoil grid family is obtained by generating
a sequence of meshes with suitably varying grid spacings and growth ratios. Eça and Hoek-
stra introduced modifications to the GCI calculation to prevent under- and over-estimating
discretization error due to scatter in the data and the possibility of the data being outside
the asymptotic range [48]. The following modification is used in the current work according
to the Turbulence Modelling Resources website guidelines5 :
  5 https://turbmodels.larc.nasa.gov/uncertainty summary.pdf, accessed November 2021

<!-- ===== PDF page 140 ===== -->

120                       APPENDIX B. FURTHER EVALUATION OF ITERATIVE AND GRID CONVERGENCE

                   Table B.1: NLF0416 structured multi-block C-grid dimensions.
                   grid level   chord × off-wall nodes   ∆s (chord)    avg/max y+
                      L0              443 × 79           5.00 × 10−6    0.61 / 1.65
                      L1              613 × 109          3.54 × 10−6    0.42 / 1.03
                      L2              885 × 157          2.50 × 10−6    0.30 / 0.63
                      L3             1225 × 217          1.77 × 10−6    0.21 / 0.42
                      L4             1769 × 313          1.25 × 10−6    0.15 / 0.28
                      L5             2449 × 433          0.89 × 10−6    0.13 / 0.21

                                                    
                                max 1.25·EFE  1.25∆M
                                              , |f1 |   p > 3.05,
                                
                                       r3 −1
                     GCI =                   21
                                                                                        (B.4)
                                 1.25·EFE
                                
                                     p                        0.95 ≤ p ≤ 3.05,
                                   r21 −1

where,

                            ∆M = max(|f2 − f1 |, |f3 − f2 |, |f3 − f1 |).               (B.5)

B.2      Results

B.2.1    NLF0416 General Aviation Airfoil

As discussed in Section 4.1.1, the grids provided by the TMPCS committee [27] have high
aspect-ratio cells emanating from the sharp trailing-edge that extend deep into the far field,
which increases the stiffness of the linear system and delays convergence. This problem
is amplified with increased grid refinement. To facilitate deep and efficient iterative con-
vergence, new structured C-topology multi-block grids with characteristics similar to the
grids provided by the TMPCS committee are used in the following section, but with in-
creased splaying and therefore reduced cell aspect ratios in the far field, with leading- and
trailing-edge aspect ratios of approximately 200 used for the initial grid-refinement study.
The dimensions of the structured C-grids are presented in Table B.1.
   Grid-refinement studies are presented at zero angle of attack using the grids in Table B.1
for fully turbulent flow (SA) and with free transition (SA-sLM2015). Low-Mach-number
preconditioning is not used for the results presented; however, simulations at a higher Mach
number of 0.30 were performed with the results demonstrating similar iterative and grid
convergence. The residual convergence histories for the fully turbulent and free-transition
grid-refinement studies are presented in Figure B.2. The total residual converges approx-
imately 12–13 orders of magnitude for each case, with the solver exiting after hitting the
absolute residual tolerance of 10−10 . Although the simulations performed on the L5 grid
level satisfy the same convergence tolerance, the results are not included in Figure B.2. The
additional stiffness from the small grid spacings significantly affects the linear solver per-
formance, requiring larger Krylov subspace sizes and significantly more equivalent residual
evaluations to converge. The free-transition simulations converge with approximately 1.2

<!-- ===== PDF page 141 ===== -->

B.2. RESULTS                                                                                                                                     121

                               103                                                                  103
                                                                     L0                                                                   L0
                               101                                   L1                             101                                   L1

       Relative Residual, d

                                                                            Relative Residual, d
                               10−1                                  L2                             10−1                                  L2
                                                                     L3                                                                   L3
                               10−3                                  L4                             10−3                                  L4
                               10−5                                                                 10−5
                               10−7                                                                 10−7
                               10−9                                                                 10−9
                          10−11                                                                10−11
                          10−130.0        0.5       1.0     1.5      2.0                       10−130.0        0.5      1.0      1.5      2.0
                                      Equivalent Residual Evaluations 1e5                                  Equivalent Residual Evaluations 1e5
                                                (a) SA                                                         (b) SA-sLM2015

Figure B.2: Residual convergence histories for the NLF0416 airfoil simulated fully turbulent and with free
transition at M a = 0.10, Re = 4.0 × 106 , T u = 0.15%, and α = 0◦ using the grids presented in Table B.1.

to 1.6 times the number of equivalent residual evaluations required for the fully turbulent
simulations.
   A linearization study for the free-transition simulations on the L4 grid level is presented,
with the residual convergence histories illustrated in Figure B.3. With the source-term
time step restriction, the fully coupled solution strategy requires more equivalent residual
evaluations to converge relative to the loosely coupled and decoupled approaches. Although
the increased coupling can provide a more accurate solution update, the linearization of a
large source term can require a smaller local time step to prevent unstable solution updates,
as demonstrated in Figure 3.1. This behaviour is illustrated by comparing Figures B.3a and
B.3b, where the fully coupled solution without the source-term time step restriction fails to
converge. Counter to the behaviour of the fully coupled simulations, the loosely coupled and
decoupled strategies require more equivalent residual evaluations to converge with the source-
term time step restriction. This behaviour can likely be attributed to the relative simplicity
of the problem. As the flow complexity increases, the more accurate update provided by the
fully coupled approach improves convergence, as demonstrated in Sections B.2.2 and B.2.3.
   Diskin et al. [41] demonstrate that the grid aspect ratio at geometric singularities, such as
a sharp trailing edge, is the primary factor affecting accuracy and grid convergence of fully
turbulent simulations. To investigate the grid convergence of free-transition simulations,
grid-refinement studies were performed with trailing-edge grid aspect ratios of 200, 20, and 2.
The grids in the current work are generated with streamwise cluster locations only specified
at the leading and trailing edges. Therefore, as the trailing-edge aspect ratio decreases,
the streamwise spacings on the upper and lower surfaces of the airfoil increase away from
the trailing edge. While this does not have a noticeable impact on fully turbulent results,
the transition onset locations move downstream as the streamwise grid spacing decreases in
the vicinity of the transition points with increasing trailing-edge aspect ratio. These small
changes in the transition onset location can significantly affect both lift and drag. A trailing-

<!-- ===== PDF page 142 ===== -->

122                                                       APPENDIX B. FURTHER EVALUATION OF ITERATIVE AND GRID CONVERGENCE

                               103                                                                                103
                                                                 fully coupled                                                                    fully coupled
                               101                               turb-trans coupled                               101                             turb-trans coupled
       Relative Residual, d

                                                                                          Relative Residual, d
                               10−1                              trans coupled                                    10−1                            trans coupled
                                                                 decoupled                                                                        decoupled
                               10−3                                                                               10−3
                               10−5                                                                               10−5
                               10−7                                                                               10−7
                               10−9                                                                               10−9
                          10−11                                                                              10−11
                          10−130.0              0.5         1.0           1.5       2.0                      10−130.0               0.5         1.0      1.5             2.0
                                         Equivalent Residual Evaluations 1e5                                                Equivalent Residual Evaluations 1e5
                                          (a) w/o STTS restriction                                                            (b) w/ STTS restriction

Figure B.3: Linearization study residual convergence histories for the NLF0416 airfoil simulated on the L4
grid (Table B.1) using the SA-sLM2015 transition model at M a = 0.10, Re = 4.0 × 106 , T u = 0.15%, and
α = 0◦ .

                           .0097                                            SA                                                                           SA
                                          pL3-L5 = 1.07                                                                  pL3-L5 = 1.01
                                                                            SA-sLM2015        0.488                                                      SA-sLM2015

                           .0096
                                                                                              0.487
                                        AFE = 0.37%                                                                    AFE = 0.09%
                                        GCI = 0.47%                                                                    GCI = 0.11%
                           .0095
                         Cd

                                   //                                                                             //
                                                                                             Cl

                           .0053                                                              0.452
                                        AFE = 0.20%
                                        GCI = 0.25%                                                                      pL3-L5 = 1.01

                           .0052                                                              0.451

                                          pL3-L5 = 1.20                                                                AFE = 0.03%
                                                                                                                       GCI = 0.04%
                           .0051                                                              0.450
                                0                  5E-06   N-1    1E-05         1.5E-05            0                              5E-06   N-1    1E-05         1.5E-05

Figure B.4: Grid-convergence results for the NLF0416 airfoil simulated fully turbulent and with free transi-
tion at M a = 0.10, Re = 4.0 × 106 , T u = 0.15%, and α = 0◦ .

edge aspect ratio of approximately 20 is used for the grid-convergence results presented, as
it provides sufficient streamwise grid resolution at the transition locations and at the trailing
edge. The leading-edge aspect ratio is maintained at approximately 200. The fully turbulent
and free-transition grid-convergence results produced at zero angle of attack are presented
in Figure B.4.
   The fully turbulent and free-transition lift and drag coefficient grid-convergence metrics
for the finest grid level are calculated using the methods described in Section B.1.2. The
free-transition grid convergence for lift is similar to the fully turbulent results, but with
higher AFE and GCI values on the finest grid level. The drag coefficient decreases with
increasing grid refinement for the free-transition simulations and increases with refinement for
the simulations with fully turbulent flow, with the free-transition results producing smaller
AFE and GCI values. The NLF0416 sharp trailing-edge transition test case highlights the
impact of streamwise grid spacing at the transition locations in addition to the trailing-edge

<!-- ===== PDF page 143 ===== -->

B.2. RESULTS                                                                                                                                          123

                                103                                                                   103
                                                                        L0                                                                      L0
                                101                                     L1                            101                                       L1

        Relative Residual, d

                                                                              Relative Residual, d
                                10−1                                    L2                            10−1                                      L2
                                                                        L3                                                                      L3
                                10−3                                                                  10−3
                                10−5                                                                  10−5
                                10−7                                                                  10−7
                                10−9                                                                  10−9
                           10−11                                                                 10−11
                           10−13
                              0.00                  0.25         0.50                            10−13
                                                                                                    0.00                    0.25         0.50
                                       Equivalent Residual Evaluations 1e5                                     Equivalent Residual Evaluations 1e5
                                          (a) α = −0.40◦ : SA                                                (b) α = −0.40◦ : SA-sLM2015cc
                                103                                                                   103
                                                                        L0                                                                      L0
                                101                                     L1                            101                                       L1
        Relative Residual, d

                                                                              Relative Residual, d
                                10−1                                    L2                            10−1                                      L2
                                                                        L3                                                                      L3
                                10−3                                                                  10−3
                                10−5                                                                  10−5
                                10−7                                                                  10−7
                                10−9                                                                  10−9
                           10−11                                                                 10−11
                           10−130.0           0.5          1.0          1.5                      10−130.0             0.5          1.0          1.5
                                       Equivalent Residual Evaluations 1e5                                     Equivalent Residual Evaluations 1e5
                                           (c) α = 1.80◦ : SA                                                (d) α = 1.80◦ : SA-sLM2015cc

Figure B.5: Residual convergence histories for the VA-2 airfoil simulated fully turbulent and with free
transition at M = 0.71, Re = 10 × 106 , T u = 0.25%, and α = −0.40◦ , 1.80◦ .

geometric singularity on the grid convergence of free-transition simulations. The results
also demonstrate the potential benefits of using adaptive mesh refinement and redistribution
methods for free-transition simulations. The blunt trailing-edge VA-2 and CRM-NLF cases
presented in the following sections are simulated using O-topology grids with trailing-edge
aspect ratios on the order of 100. The need for a small grid aspect ratio at the trailing edge
is not as acute for these cases.

B.2.2                 VA-2 Supercritical Airfoil

Grid-refinement studies are presented for fully turbulent flow (SA) and with free transition
(SA-sLM2015cc) at two angles of attack, −0.40 and 1.80 degrees, using the grids presented
in Section 4.2.2. The residual convergence histories for the total residual for the fully turbu-
lent and free-transition grid-refinement studies are illustrated in Figure B.5. Similar to the
NLF0416 results presented in Section B.2.1, all cases exit after hitting the absolute residual
convergence tolerance of 10−10 . The free-transition results demonstrate a significant increase
in the computational cost for the 1.80 degree angle of attack case, which features shock-
induced transition on the upper surface of the airfoil at all grid levels, relative to the −0.40

<!-- ===== PDF page 144 ===== -->

124                                              APPENDIX B. FURTHER EVALUATION OF ITERATIVE AND GRID CONVERGENCE

                               103                                                                     103
                                                         fully coupled                                                           fully coupled
                               101                       turb-trans coupled                            101                       turb-trans coupled
       Relative Residual, d

                                                                               Relative Residual, d
                               10−1                      trans coupled                                 10−1                      trans coupled
                                                         decoupled                                                               decoupled
                               10−3                                                                    10−3
                               10−5                                                                    10−5
                               10−7                                                                    10−7
                               10−9                                                                    10−9
                          10−11                                                                   10−11
                          10−130.0                 0.5                  1.0                       10−130.0                 0.5                 1.0
                                      Equivalent Residual Evaluations 1e5                                     Equivalent Residual Evaluations 1e5
                                (a) α = −0.40◦ : w/o STTS restriction                                   (b) α = −0.40◦ : w/ STTS restriction
                               103                                                                     103
                                                         fully coupled                                                           fully coupled
                               101                       turb-trans coupled                            101                       turb-trans coupled
       Relative Residual, d

                                                                               Relative Residual, d
                               10−1                      trans coupled                                 10−1                      trans coupled
                                                         decoupled                                                               decoupled
                               10−3                                                                    10−3
                               10−5                                                                    10−5
                               10−7                                                                    10−7
                               10−9                                                                    10−9
                          10−11                                                                   10−11
                          10−130.0         0.5      1.0        1.5       2.0                      10−130.0        0.5      1.0         1.5       2.0
                                      Equivalent Residual Evaluations 1e5                                     Equivalent Residual Evaluations 1e5
                                 (c) α = 1.80◦ : w/o STTS restriction                                    (d) α = 1.80◦ : w/ STTS restriction

Figure B.6: Linearization study residual convergence histories for the VA-2 airfoil simulated on the L3 grid
(Table 4.5) using the SA-sLM2015cc transition model at M = 0.71, Re = 10 × 106 , T u = 0.25%, and
α = −0.40◦ , 1.80◦ .

degrees angle of attack case that has a near-zero pressure gradient on the upper surface.
    A linearization study for the free-transition simulations at each angle of attack, −0.40
and 1.80 degrees, on the L3 grid is presented, with the total residual convergence histories
illustrated in Figure B.6. For the −0.40 degree case, with the source-term time step restric-
tion all four linearization strategies produce similar convergence behaviour. The number of
equivalent residual evaluations required to converge increases for each linearization strategy
without the source-term time step restriction, demonstrating that the effect of the source
terms is more significant for this case relative to the NLF0416 simulations. At 1.80 degrees
angle of attack, with the source-term time step restriction the fully coupled linearization
strategy provides improved nonlinear convergence relative to the loosely coupled and decou-
pled approaches. The fully coupled linearization strategy provides a more complete Jacobian
and produces a more accurate solution update, which significantly improves convergence for
this more challenging case. All solutions fail to converge within the allocated compute time
at 1.80 degrees angle of attack without the source-term time step restriction.
   The grid-convergence results for the fully turbulent and free-transition simulations are
presented in Figure B.7. For both angles of attack, the fully turbulent and free-transition

<!-- ===== PDF page 145 ===== -->

B.2. RESULTS                                                                                                              125

                       AFE = 0.05%                 SA                                                      SA
                       GCI = 0.07%                 SA-sLM2015cc   0.360                                    SA-sLM2015cc
          .0089                                                                   pL1-L3 = 1.35

                                                                  0.350
          .0088
                          pL1-L3 = 2.26
                                                                  0.340        AFE = 1.17%
          .0087                                                                GCI = 1.48%
         Cd
                  //                                                      //

                                                                  Cl
                       AFE = 0.18%
          .0060        GCI = 2.07%                                0.320           pL1-L3 = 1.02

          .0059                                                   0.310

          .0058           pL1-L3 = 3.44                           0.300        AFE = 1.50%
                                                                               GCI = 1.91%
                             5E-06        N-1   1.5E-05                              5E-06        N-1   1.5E-05
                                                                          ◦
                                                          (a) α = −0.40

          .0104                                                   0.810
                       AFE = 0.12%                 SA                                                      SA
                       GCI = 0.14%                 SA-sLM2015cc                   pL1-L3 = 1.12            SA-sLM2015cc
          .0103                                                   0.800

          .0102                                                   0.790
                          pL1-L3 = 1.91                                        AFE = 0.88%
                                                                               GCI = 1.10%
          .0101                                                   0.780
         Cd

                  //                                                      //
                                                                  Cl

                                                                  0.750
                       AFE = 0.24%                                                pL1-L3 = 1.10
          .0076        GCI = 0.30%
                                                                  0.740
          .0075
                                                                  0.730
                          pL1-L3 = 2.25                                        AFE = 0.71%
          .0074                                                                GCI = 0.90%
                                                                  0.720
                             5E-06        N-1   1.5E-05                              5E-06        N-1   1.5E-05

                                                          (b) α = 1.80◦

Figure B.7: Grid-convergence results for the VA-2 airfoil simulated fully turbulent and with free transition
at M = 0.71, Re = 10 × 106 , T u = 0.25%, and α = −0.40◦ , 1.80◦ .

simulations produce similar lift coefficient grid convergence. At −0.40 degrees angle of
attack, the free-transition results produce smaller AFE and GCI values and a higher order
of convergence relative to the fully turbulent results; however, at 1.80 degrees the results
are reversed, with the free-transition results producing larger AFE and GCI values. At both
angles of attack, the free-transition drag coefficient grid-convergence results demonstrate an
increase in the discretization error produced on the finest grid level relative to the fully
turbulent results. The results demonstrate that additional mesh resolution is required for
free-transition simulations in order to achieve a similar level of accuracy in the drag coefficient
to the fully turbulent simulations.

B.2.3     NASA CRM-NLF Wing-Body Geometry

The NASA CRM-NLF grid-refinement study iterative convergence at the 2524 run condi-
tions with both fully turbulent flow and with free transition was presented previously in Sec-

<!-- ===== PDF page 146 ===== -->

126                                               APPENDIX B. FURTHER EVALUATION OF ITERATIVE AND GRID CONVERGENCE

                              103                                                                  103
                                                      fully coupled                                101
                                                                                                                          fully coupled
                                                      turb-trans coupled                                                  turb-trans coupled
      Relative Residual, d

                                                                           Relative Residual, d
                              10−1                    trans coupled                                10−1                   trans coupled
                              10−3                    decoupled                                    10−3                   decoupled
                              10−5                                                                 10−5
                              10−7                                                                 10−7
                              10−9                                                                 10−9
                         10−11                                                                10−11
                         10−13                                                                10−13
                         10−150.0           0.5          1.0         1.5                      10−150.0           0.5         1.0         1.5
                                     Equivalent Residual Evaluations 1e5                                  Equivalent Residual Evaluations 1e5
                                     (a) w/o STTS restriction                                              (b) w/ STTS restriction

Figure B.8: Linearization study residual convergence histories for the NASA CRM-NLF simulated on the
L2 grid (Table 4.7) using the SA-QCR2000-sLM2015cc transition model at approximately M = 0.86, Re =
15 × 106 , T u = 0.24%, and α = 2.0◦ .

tion 4.2.3. In this section, a linearization study for the free-transition simulations performed
using the SA-sLM2015cc transition model on the L2 grid is presented, with the residual
convergence histories illustrated in Figure B.8. As demonstrated in Figure B.8a, without
the source-term time step restriction all linearization strategies fail to converge and instead
stall after a 4 order of magnitude drop in the total residual. The fully coupled lineariza-
tion strategy without the source-term time step restriction exits after approximately 5 × 104
equivalent residual evaluations after the solver is unable to prevent an unphysical solution
update. However, with the source-term time step restriction (Figure B.8b) the fully coupled
linearization strategy provides the most robust nonlinear convergence. The turb-trans cou-
pled and trans coupled solutions stall after the total residual drops 5 orders of magnitude.
The decoupled solution strategy appears to be converging but exits due to a wall-clock con-
straint. The results in Figure B.8 demonstrate that the fully coupled linearization strategy
provides a more accurate solution update that improves nonlinear convergence; however, the
source-term time step restriction is required in order to prevent unstable solution updates.
   The fully turbulent and free-transition grid-convergence results are presented in Fig-
ure B.9. While the fully turbulent and free-transition simulations produce similar drag
coefficient grid convergence, the free-transition lift coefficient grid convergence produces a
larger GCI value. The high free-transition lift coefficient order of convergence demonstrates
that the results are not in the asymptotic range of convergence [48], and indicates that
additional grid-refinement levels are required in order to draw further conclusions.

<!-- ===== PDF page 147 ===== -->

                                         SA-QCR2000                 0.430                           SA-QCR2000
                                         SA-QCR2000-sLM2015cc                                       SA-QCR2000-sLM2015cc

                                                                                    pL1-L3 = 4.15
                         pL1-L3 = 1.71                              0.420
         .0187

                 AFE = 0.78%                                        0.410   AFE = 0.18%
                 GCI = 0.99%                                                GCI = 1.21%
          CD

                                                                    CL

                         pL1-L3 = 1.64                              0.400           pL1-L3 = 2.81
         .0182

                                                                    0.390
                 AFE = 0.65%                                                AFE = 0.39%
                 GCI = 0.82%                                                GCI = 0.48%
         .0177                                                      0.380
           0.0E+00       1.0E-05              2.0E-05   3.0E-05       0.0E+00       1.0E-05           2.0E-05     3.0E-05
                                       -2/3
                                   N                                                          N-2/3
Figure B.9: Grid-convergence results for the NASA CRM-NLF simulated fully turbulent and with free
transition at approximately M = 0.86, Re = 15 × 106 , T u = 0.24%, and α = 2.0◦ .

<!-- ===== PDF page 149 ===== -->

Appendix C

Transition Length Modification

The airfoil optimizations presented in Section 5.1 demonstrate that the capability of the
optimization algorithm to delay boundary-layer transition is sensitive to the streamwise grid
resolution in the transition regions, and that streamwise grid requirements increase as the
transition length decreases with increasing Reynolds number. However, instead of increasing
the number of streamwise nodes, an alternative approach is to increase the transition region
length. Rashad demonstrated that slowing the ramp-up of the eddy viscosity, which slows
the growth of the turbulent skin friction coefficient and therefore increases the transition
length, improved the smoothness of the design space, and found that a transition length of
10% chord was sufficient for NLF optimization [149]. However, Rashad used an algebraic
intermittency function which provides greater control over the transition length relative to
transport-equation-based intermittency models. The following chapter presents an initial
attempt at increasing the transition region length produced by the intermittency transport
equation to the value of 10% chord used by Rashad [149].

C.1     Modified Intermittency Source Terms

As a first step, the intermittency production and destruction terms were combined into a
single source term, Sγ , given by,
                                  √       
                                 M∞ M∞ Re∞
                   Sγ = φ−300 Ω,             (0.98Fonset − γ + 0.02).                    (C.1)

The source term was simplified by removing the Flength and Fturb functions. The first term in
Equation C.1 is the vorticity limiting procedure described previously in Section 2.2.2. Similar
to the original source terms in the SA-LM2015 and SA-sLM2015 transition models, the new
source term is controlled by Fonset and implicitly enforces the stable bounds for intermittency
of 0.02 and 1. While the transition length in the original model is controlled by the Flength
function, which scales the intermittency production term, Pγ , decreasing Flength in an attempt
to enforce a transition length of approximately 10% chord affected the predicted transition

<!-- ===== PDF page 150 ===== -->

130                                                     APPENDIX C. TRANSITION LENGTH MODIFICATION

           1                                              1

          0.8                                            0.8

          0.6                                            0.6

          0.4                                            0.4

          0.2                                            0.2

           0                                              0
                0    1       2        3        4    5          0     10         20            30   40

                    (a) Transition criterion                       (b) Eddy viscosity ratio

Figure C.1: Sensitivity of the different Fonset formulations to the transition criterion (Equation 2.35) and
the eddy viscosity ratio.

onset location, which moved downstream with decreasing Flength . Instead of controlling the
transition length using Flength , the Fonset function was modified as follows to remain active
over a larger portion of the transition process:

                     Fonset = φ−300 (φ300 (Fonset 2 , 0), 1),                                           (C.2)
                                                       
                    Fonset 2 = 0.1575 Fonset 1 − 1 ,                                                    (C.3)
                                           2                r      32  12
                                       ReS                    3     µ
                    Fonset 1 = 10                 + 3 × 10 RT                     .                     (C.4)
                                     2.6Reθc                       ρUe

   The Fonset function for the SA-LM2015 (Equations 2.32–2.35) and SA-sLM2015 (Equa-
tions 2.49–2.50) transition models, along with the modified function defined above, are il-
lustrated in Figure C.1. The sensitivity of the Fonset function to the transition criterion
(Equation 2.35) was reduced to increase the transition length in the transition onset region,
while the sensitivity to the eddy viscosity ratio was greatly reduced. As noted by Langtry and
Menter [92], vorticity decreases once the transition process begins; therefore, downstream
of transition onset the growth of the Fonset function is controlled by the growth of the eddy
viscosity. By reducing the sensitivity to the eddy viscosity ratio, the new Fonset formulation
significantly increases the transition length.
   The growth of the eddy viscosity ratio, and therefore the growth of the Fonset function
downstream of the transition onset location, is sensitive to the freestream Mach and Reynolds
numbers. To ensure a constant transition length over a range of flight conditions, an attempt
was made to normalize the eddy viscosity ratio in the new Fonset function. The inverse of the
square root of the velocity magnitude at the edge of the boundary layer is included to nor-
malize the eddy viscosity ratio to Mach number, while the square root of the kinematic eddy
viscosity is included to normalize eddy viscosity ratio to Reynolds number. To verify that

<!-- ===== PDF page 151 ===== -->

C.1. MODIFIED INTERMITTENCY SOURCE TERMS                                                                                                       131

                                                         Re = 5e6                                                           Ma = 0.66
         0.015                                                -1.5
                                                         Re = 10e6          0.015                                           Ma = -1.5
                                                                                                                                  0.71
                                                         Re = 15e6                                                          Ma = 0.76
                                                                  -1                                                                 -1

          0.01                                                    -0.5       0.01                                                    -0.5

         Cf,x

                                                                            Cf,x
                                                                       Cp

                                                                                                                                          Cp
                                                                  0                                                                  0

         0.005                                                              0.005
                                                                  0.5                                                                0.5

                                                                  1                                                                  1

                0                                                                  0
                                                                  1.5                                                                1.5
                    0   0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9   1                        0   0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9   1
                                       x/c                                                                x/c

           (a) Varying Reynolds number; M a = 0.71                            (b) Varying Mach number; Re = 10 × 106

Figure C.2: Pressure and upper-surface skin friction coefficient profiles produced by simulations of the VA-2
airfoil on the L0 grid level (Table 4.5) with varying Mach and Reynolds number; α = 0.00, T u = 0.25%.

the transition length produced by the new Fonset function (Equations C.2–C.4) and intermit-
tency source-term formulation (Equation C.1) is insensitive to freestream flow quantities,
simulations are performed of the VA-2 airfoil with varying Mach and Reynolds numbers
using the L0 grid presented in Table 4.5 and are compared with results obtained using the
smooth Pγ and Eγ source terms presented in Section 2.2.2, including the original, smooth
Fonset function (Equations 2.49–2.50), with the results presented in Figure C.2. The upper
surface skin friction coefficient profiles demonstrate that the new source-term formulation
with the modified Fonset function produces an increased transition length of approximately
10% chord, and that the transition length is insensitive to changes in the freestream Mach
and Reynolds number, while producing similar transition onset locations.
   Further work is required to develop an equivalent Fonset,scf source term, as previous at-
tempts using a similar approach to Equations C.2–C.4 were unsuccessful. However, to verify
that the new intermittency source-term formulation produces a consistent transition length
for both two- and three-dimensional flows, simulations are presented of an RAE2822 airfoil
extruded one chord length with and without a 25-degree sweep at the Boeing 737-800 design
conditions, with sweep corrections applied to the unswept design. The grids are representa-
tive of the L1 grid level presented in Table 5.6, while the design conditions for both configu-
rations are presented in Table C.1. It is important to note that these sweep corrections are
more thorough than the sweep-corrected Mach number applied to the airfoil optimizations
at the Boeing 737-800 design conditions (Table 5.1) in order to obtain a one-to-one match of
the pressure distributions between the unswept and swept designs, with corrections applied
to the lift coefficient, Mach number, Reynolds number, and maximum thickness-to-chord
ratio based on the leading-edge sweep angle. The results in Figure C.3 demonstrate that the
new intermittency source-term formulation produces very similar transition onset locations
and transition lengths for both the swept and sweep-corrected configurations.

<!-- ===== PDF page 152 ===== -->

132                                                                  APPENDIX C. TRANSITION LENGTH MODIFICATION

                                                               B737-800
                                                                                         0.015
                                                               B737-800: corrected
                                        -1

                                      -0.5
                                                                                         0.01

                                                                                             Cf,x
                                      Cp
                                       0.5                                               0.005

                                       1.5
                                               0   0.2   0.4         0.6   0.8       1
                                                               x/c

Figure C.3: Pressure and skin friction coefficient profiles produced by simulations of a swept and sweep-
corrected RAE2822 airfoil at the Boeing 737-800 design conditions (Table C.1) with the new intermittency
source-term formulation.

Table C.1: Design conditions for the extruded RAE2822 airfoil simulated at the swept and sweep-corrected
Boeing 737-800 design conditions. For each configuration turbulence intensity, T u, is specified as 0.07%.
          sweep    lift coefficient        Mach number         Reynolds number (MAC)                relative (t/c)max
           25◦          0.5000               0.78500                20.3000 × 106                            –
            0◦          0.6087               0.71145                16.6743 × 106                         1.1034

C.2      Airfoil Optimization

The airfoil drag minimizations presented in Section 5.1 are revisited in order to investigate
the impact of the new intermittency source-term formulation on the streamwise grid require-
ments for NLF optimization. The infinite swept wing optimization is not investigated as the
transition length modification has not yet been extended to include the stationary crossflow
instability source term. It is important to note that the parameters for the optimizations
are not changed from those presented in Section 5.1 in order to isolate any differences due
to the new source terms.

C.2.1    Cessna 172R Skyhawk

Drag minimizations are performed at the Cessna light aircraft design conditions presented
in Table 5.1 using the L1 grid level (Table 5.2). The optimizations are performed using
the baseline optimization problem (Equations 5.1–5.4) and the optimization problem with
additional constraints (Equations 5.5–5.9), with the results presented in Figure C.4. The
results demonstrate that the new source-term formulation produces similar designs to those
presented in Figure 5.2, with the new formulation producing similar local minima. The
design produced with additional constraints is re-optimized using the baseline optimization
problem with the results also presented in Figure C.4. Similar to the results presented
in Section 5.1.1, the light aircraft design space appears to be multi-modal with the new

<!-- ===== PDF page 153 ===== -->

C.2. AIRFOIL OPTIMIZATION                                                                                                                                                                                                                                                                                                                                         133

                                                     Optimality                                                                                        2
                                                                                                                                                                              Optimality                                                                                          2
                                                                                                                                                                                                                                                                                                      Optimality
                           10                        Feasibility                                                .0060                             10                          Feasibility                                                  .0060                             10                       Feasibility                                               .0060

                                                                                           Merit function

                                                                                                                                                                                                                      Merit function

                                                                                                                                                                                                                                                                                                                                           Merit function
                                                     Merit function                                                                                                           Merit function                                                                                                          Merit function
                                                                                                                .0055                                                                                                                      .0055                                                                                                                .0055

 Optimality, Feasibility

                                                                                                                        Optimality, Feasibility

                                                                                                                                                                                                                                                   Optimality, Feasibility
                           100                                                                                                                    100                                                                                                                        100
                                                                                                                .0050                                                                                                                      .0050                                                                                                                .0050

                                -2                                                                              .0045                             10
                                                                                                                                                       -2                                                                                  .0045                             10
                                                                                                                                                                                                                                                                                  -2                                                                            .0045

                                                                                                                .0040                                                                                                                      .0040                                                                                                                .0040
                           10-4                                                                                 .0035                             10-4                                                                                     .0035                             10-4                                                                               .0035

                                -6
                                                                                                                .0030                                  -6
                                                                                                                                                                                                                                           .0030                                  -6
                                                                                                                                                                                                                                                                                                                                                                .0030
                           10                                                                                                                     10                                                                                                                         10
                                                                                                                .0025                                                                                                                      .0025                                                                                                                .0025

                           10-8                                                                                 .0020                             10-8                                                                                     .0020                             10-8                                                                             .0020
                                     0   5    10     15      20      25    30      35      40                                                               0            10                  20                  30                                                                    0         5            10             15                             20
                                                   Design iteration                                                                                                      Design iteration                                                                                                            Design iteration

                                                                                Optimized                                                                                                             Optimized                                                                                                              Optimized
                                                                                Initial                         .015                                                                                  Initial                              .015                                                                              Initial                            .015

                       -0.6                                                                                                                   -0.6                                                                                                                       -0.6
                                                                                                                .012                                                                                                                       .012                                                                                                                 .012

                           0.0                                                                                                                    0.0                                                                                                                        0.0
                                                                                                                .009                                                                                                                       .009                                                                                                                 .009
 Cp

                                                                                                                        Cp

                                                                                                                                                                                                                                                   Cp
                           0.6                                                                                                                    0.6                                                                                                                        0.6
                                                                                                                  Cf

                                                                                                                                                                                                                                             Cf

                                                                                                                                                                                                                                                                                                                                                                  Cf
                                                                                                                .006                                                                                                                       .006                                                                                                                 .006
                           1.2                                                                                                                    1.2                                                                                                                        1.2
                                                   upper                                                                                                                 upper
                                                                                                                .003                                                                                                                       .003                                                                                                                 .003
                           1.8                               lower                                                                                1.8                                lower                                                                                   1.8                                          upper
                                                                           upper                                                                                                                   upper
                                                                                   lower                                                                                                                                                                                                                                  upper
                                                                                                                                                                                                           lower                                                                                                                   lower
                                                                                                                                                                                                                                                                                                                                  lower
                           2.4                                                                                  .000                              2.4                                                                                      .000                              2.4                                                                                .000
                                     0       0.2       0.4           0.6         0.8                        1                                               0      0.2         0.4           0.6           0.8                         1                                               0   0.2          0.4         0.6           0.8                       1
                                                             x/c                                                                                                                     x/c                                                                                                                      x/c

(a) Baseline Optimization Problem                                                                                                                               (b) Additional Constraints                                                         (c) Re-optimized w/ Baseline Opti-
filler                                                                                                                                                                     filler                                                                  mization Problem

Figure C.4: Optimization convergence histories and cross-sectional profiles of the initial and optimized
designs produced by drag minimizations with the baseline (Equations 5.1–5.4) and more constrained op-
timization problem (Equations 5.5–5.9) at the Cessna 172R conditions (Table 5.1) on the L1 grid level
(Table 5.2).

intermittency source-term formulation. Although optimizations were not performed with
grids with coarser streamwise grid resolutions, these results confirm that the new formulation
produces similar designs to the original source-term formulation using the same grid with
the same design conditions.

C.2.2                                        De Havilland Dash8-Q400

Drag minimizations at the De Havilland turbo-prop aircraft design conditions presented in
Section 5.1.2 demonstrate that the L1 grid level (Table 5.4) provides insufficient stream-
wise grid resolution in the transition region produced by the original SA-sLM2015 transition
model intermittency source terms, and that 1121 streamwise nodes are required to adequately
resolve the streamwise gradients. Simulations using the L1 grid level with 561 streamwise
nodes are performed using the new source-term formulation with both the baseline optimiza-
tion problem (Equations 5.1–5.4) and the optimization problem with additional constraints
(Equations 5.5–5.9), with the results presented in Figure C.5. Comparing with the results
in Figure 5.3, the results in Figure C.5 demonstrate that the new source-term formulation
produces designs with significantly more laminar flow on the L1 grid level. The increased
transition length reduces the streamwise grid requirements in order to maintain adequate

<!-- ===== PDF page 154 ===== -->

134                                                                                                                                                   APPENDIX C. TRANSITION LENGTH MODIFICATION

                                                                     Optimality                                                                                     2
                                                                                                                                                                                            Optimality
                                           10                        Feasibility                                             .0060                             10                           Feasibility                                           .0060

                                                                                                        Merit function

                                                                                                                                                                                                                             Merit function
                                                                     Merit function                                                                                                         Merit function
                                                                                                                             .0055                                                                                                                .0055

                 Optimality, Feasibility

                                                                                                                                     Optimality, Feasibility
                                           100                                                                                                                 100

                                                                                                                             .0050                                                                                                                .0050
                                                -2                                                                                                                  -2
                                           10                                                                                                                  10
                                                                                                                             .0045                                                                                                                .0045
                                           10-4                                                                                                                10-4
                                                                                                                             .0040                                                                                                                .0040

                                                -6                                                                                                                  -6
                                           10                                                                                                                  10
                                                                                                                             .0035                                                                                                                .0035

                                                -8                                                                                                                  -8
                                           10                                                                                .0030                             10                                                                                 .0030
                                                     0           5               10                15                                                                    0             5            10            15           20
                                                                 Design iteration                                                                                                          Design iteration

                                                                                                Optimized                                                                                                             Optimized
                                                                                                Initial                      .015                                                                                     Initial                     .015
                                       -1.2                                                                                                                -1.2

                                       -0.6                                                                                  .012                          -0.6                                                                                   .012

                                           0.0                                                                                                                 0.0
                                                                                                                             .009                                                                                                                 .009
                                           0.6                                                                                                                 0.6
                 Cp

                                                                                                                                     Cp
                                                                                                                               Cf

                                                                                                                                                                                                                                                    Cf
                                           1.2                                                                               .006                              1.2                                                                                .006

                                           1.8                                                                                                                 1.8
                                                         upper                                                                                                                 upper
                                                                                                                             .003                                                                                                                 .003
                                           2.4                        lower                                                                                    2.4                            lower
                                                                              upper                                                                                                                   upper
                                                                                                                                                                                                              lower
                                                                                        lower
                                           3.0                                                                               .000                              3.0                                                                                .000
                                                     0     0.2        0.4             0.6        0.8                     1                                               0       0.2          0.4           0.6        0.8                    1
                                                                              x/c                                                                                                                     x/c

                               (a) Baseline Optimization Problem                                                                                                             (b) Additional Constraints

Figure C.5: Optimization convergence histories and cross-sectional profiles of the initial and optimized
designs produced by drag minimizations with the baseline (Equations 5.1–5.4) and more constrained opti-
mization problem (Equations 5.5–5.9) at the De Havilland Dash8-Q400 conditions (Table 5.1) on the L1 grid
level (Table 5.4).

streamwise resolution in the transition region. In addition, the new designs produce more
laminar flow than the designs in Figure 5.4 produced using the L2 grid level with the orig-
inal source terms. Although the merit function appears to converge for both cases, similar
to Chapter 5 the results demonstrate that additional work is required to facilitate deep
convergence of optimality at these design conditions.

C.2.3    Boeing 737-800

Finally, an airfoil drag minimization is performed using the sweep-corrected Boeing single-
aisle aircraft design conditions presented in Table 5.1. Similar to the previous cases, the
optimization is performed on the L1 grid level (Table 5.6) with 561 streamwise nodes, with
the results presented in Figure C.6. While the results presented in Figure 5.5 demonstrate
that this grid level significantly under-resolves the transition region using the original source
terms, the results in Figure C.6 demonstrate that the new source-term formulation with an
increased transition length produces a similar design to that produced on the L2 grid level
with the original source terms, with comparable extents of laminar flow on both the upper
and lower surfaces.

<!-- ===== PDF page 155 ===== -->

C.3. SUMMARY                                                                                                                                                                          135

                                                                Optimality                                                                                     Optimized
                                           10                   Feasibility                            .0060                                                   Initial         .015

                                                                                      Merit function
                                                                                                                -1.2
                                                                Merit function
                                                                                                       .0055    -0.6

                 Optimality, Feasibility
                                           100                                                                                                                                 .012

                                                                                                       .0050        0.0
                                                -2
                                           10                                                                                                                                  .009
                                                                                                                    0.6

                                                                                                               Cp

                                                                                                                                                                                 Cf
                                                                                                       .0045
                                           10-4                                                                     1.2                                                        .006
                                                                                                       .0040
                                                                                                                    1.8
                                                -6                                                                                upper
                                           10                                                                                                                                  .003
                                                                                                       .0035                               lower
                                                                                                                    2.4                   lower
                                                                                                                                                     upper

                                           10-8                                                        .0030        3.0                                                        .000
                                                     0     2     4       6       8   10                                   0       0.2      0.4           0.6    0.8        1
                                                               Design iteration                                                                    x/c

                                                         (a) Optimization history                                             (b) Cross-sectional profiles

Figure C.6: Drag minimization with the baseline optimization problem (Equations 5.1–5.4) at the sweep-
corrected Boeing 737-800 conditions (Table 5.1) on the L1 grid level (Table 5.6).

C.3     Summary

The results demonstrate that increasing the transition length is an effective method for
reducing the streamwise grid requirements and computational cost of free-transition opti-
mizations. However, the modification does not address the poor convergence of optimality
at higher Reynolds numbers, which was also encountered in Chapter 5. Future work is re-
quired to address this problem. In addition, future work will involve developing an equivalent
transition length modification compatible with the crossflow onset function, Fonset,scf .
