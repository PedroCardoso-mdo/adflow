# Piotrowski & Zingg (2023, Aeronautical Journal) — Compressibility corrections to extend a smooth local correlation-based transition model to transonic flows

> **Full transcript** — verbatim text extracted from the source PDF with
> `pdftotext -layout`, page-delimited, nothing removed except bare page-number
> lines. Layout is preserved, so tables and multi-column material read as they
> appear in the PDF. **Equations are raw PDF extractions**: superscripts,
> subscripts and math symbols may be flattened or split across lines — for any
> equation you intend to implement, confirm against the PDF itself (path below)
> or the digest companion file.
>
> Source PDF: `/home/mdo/Desktop/ARTIGOS SMOTH/PiotrowskiandZingg2023-AeroJ (2).pdf`

---

<!-- ===== PDF page 1 ===== -->

            The Aeronautical Journal (2023), 1–30
            doi:10.1017/aer.2022.105

            R E G U L A R PA P E R

            Compressibility corrections to extend a smooth local
            correlation-based transition model to transonic flows
            M.G.H. Piotrowski1, 2,∗              and D.W. Zingg2
            1 Computational Aerosciences Branch, NASA Ames Research Center, Mountain View, CA, USA and 2 University of Toronto,

            Institute for Aerospace Studies, Toronto, ON, Canada
            ∗ Corresponding author. Email: m.piotrowski@mail.utoronto.ca

            Received: 17 May 2022; Revised: 25 November 2022; Accepted: 29 November 2022
            Keywords: boundary-layer transition prediction; computational ﬂuid dynamics; transonic ﬂow

            Abstract
            This paper presents progress towards a transition modelling capability for use in the numerical solution of the
            Reynolds-averaged Navier-Stokes equations that provides accurate predictions for transonic ﬂows and is thus suit-
            able for use in the design of wings for aircraft ﬂying at transonic speeds. To this end, compressibility corrections
            are developed and investigated to extend commonly used empirical correlations to transonic ﬂight conditions
            while retaining their accuracy at low speeds. A compressibility correction for Tollmien-Schlichting instabili-
            ties is developed and applied to a smooth local correlation-based transition model and a stationary crossﬂow
            instability compressibility correction is included by adding a new crossﬂow source term function. Two- and
            three-dimensional transonic transition test cases demonstrate that the Tollmien-Schlichting compressibility cor-
            rection produces substantially improved agreement with the experimental transition locations, particularly for
            higher Reynolds number applications where the eﬀects of ﬂow compressibility are expected to be more signif-
            icant, such as the NASA CRM-NLF wing-body conﬁguration, while the crossﬂow compressibility correction
            prevents an inaccurate, upstream transition front. The compressibility corrections and modiﬁcations do not sig-
            niﬁcantly aﬀect the numerical behaviour of the model, which provides an eﬃcient alternative to non-local and
            higher-ﬁdelity approaches, and can be applied to other transport-equation-based transition models with low-speed
            empirical correlations without aﬀecting their predictive capability in the incompressible regime.

            Nomenclature
             a                sound speed
             d                distance to nearest solid wall
             l                characteristic length
             M                Mach number
             N                number of grid nodes
             p                static pressure
             Re               Reynolds number based on the characteristic length, ρUl/μ
             R̃eθ t           transported quantity of the transition onset momentum-thickness Reynolds number
             Reθ c            critical momentum-thickness Reynolds number, ρUθc /μ
             Reθ t            transition onset momentum-thickness Reynolds number, ρUθt /μ
             ReS              strain-rate magnitude Reynolds number, ρd2 S/μ
             RT               eddy viscosity ratio, μt /μ
             s                streamwise coordinate
             S                strain-rate magnitude, (2Sij Sij )1/2
             Sij              strain-rate tensor, 0.5(∂ui /∂xj + ∂uj /∂xi )
             S̃               modiﬁed strain-rate magnitude
             Tu               turbulence intensity
             uτ               friction velocity, (τw /ρ)1/2

            
            C The Author(s), 2023. Published by Cambridge University Press on behalf of Royal Aeronautical Society. This is an Open Access article, distributed

            under the terms of the Creative Commons Attribution licence (http://creativecommons.org/licenses/by/4.0/), which permits unrestricted re-use,
            distribution and reproduction, provided the original article is properly cited.
https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 2 ===== -->

           2        Piotrowski and Zingg

            ui               Cartesian velocity component
            U                local velocity magnitude, (u2 + v2 + w2 )1/2
            y+               non-dimensional wall distance, ρuτ d/μ

           Greek symbol
                            vorticity magnitude, (2ij ij )1/2
            ij              vorticity tensor, 0.5(∂ui /∂xj − ∂uj /∂xi )
            γ                transported quantity of intermittency
            ρ                density
            μ                molecular viscosity
            μt               eddy viscosity
            ν̃               modiﬁed eddy viscosity
            λθ               pressure gradient parameter, (ρθ 2 /μ)(dU/ds)
            θ                momentum thickness
            θc               critical momentum thickness
            θt               transition onset momentum thickness
            φ                smooth maximum/minimum function
            ψ                compressibility correction for Tollmien-Schlichting instabilities
            ψscf             compressibility correction for stationary crossﬂow instabilities
            κ                heat capacity ratio
            τw               wall shear stress, μ(∂U/∂d)d=0

           Subscripts
            e                boundary-layer edge
            comp             compressible
            ∞                far ﬁeld
            scf              stationary crossﬂow

           1.0 Introduction
           The location of boundary-layer transition from laminar to turbulent ﬂow can have a signiﬁcant impact on
           the aerodynamic performance of a transonic wing. An upstream laminar boundary layer directly inﬂu-
           ences shock location and strength, while laminar separation bubbles can lead to severe drag penalties,
           stability and control issues, and adverse stall characteristics. Recent studies suggest the application of
           laminar ﬂow control to large commercial aircraft can reduce aerodynamic drag by approximately 10%
           [1]. In order to design wings optimised for low drag for commercial transport aircraft, which typically
           ﬂy at transonic Mach numbers and high Reynolds numbers, it is necessary to be able to accurately and
           eﬃciently predict boundary-layer transition in transonic ﬂow regimes.
               The primary transition mechanisms for a swept transport-aircraft wing result from Tollmien-
           Schlichting wave growth leading to natural transition, crossﬂow instabilities resulting from highly swept
           wings, concave curvature producing Görtler instabilities, and attachment-line instabilities as a result of
           large leading-edge radius and sweep at the root of the wing [2]. While the last two mechanisms can
           be prevented by appropriate proﬁle design, balancing natural and crossﬂow-induced transition can be
           diﬃcult, as favourable pressure gradients used to stabilise streamwise instabilities destabilise crossﬂow
           instabilities.
               The eﬀects of boundary-layer transition can be included in RANS-based simulations using vari-
           ous transition prediction and modelling techniques, consisting of either local or non-local operations
           [3–5]. Local correlation-based transition models, initially developed by Langtry and Menter [6], provide
           a framework for predicting boundary-layer transition that can be easily integrated in modern highly par-
           allel solution algorithms. However, the empirical correlations of Langtry and Menter (LM2009) [6] and
           Langtry et al. (LM2015) [7] for two-dimensional mechanisms (natural, bypass and separation-induced
           transition) and stationary crossﬂow instabilities, respectively, were developed based on results from

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 3 ===== -->

                                                                                      The Aeronautical Journal       3

            subsonic experiments and stability analysis. Recent studies have demonstrated that models based on
            these correlations severely under-predict the extent of laminar ﬂow when applied to transonic test cases
            [8, 9], such as the NASA CRM-NLF [10–12]. This behaviour was foreseen by Arnal who stated, “sim-
            ple transition criteria developed for low-speed ﬂows cannot be used with conﬁdence in conﬁgurations
            where compressibility eﬀects become signiﬁcant” [13].
                The stabilising eﬀect of Mach number was recently observed experimentally by Risuis et al. [14]
            for high subsonic Mach numbers. However, this phenomenon can be more eﬃciently studied using sta-
            bility theory due to the diﬃculties of obtaining detailed transition measurements at transonic and high
            Reynolds number ﬂow conditions [15]. One of the ﬁrst investigations into the eﬀects of ﬂow compress-
            ibility on the stability of laminar boundary layers was performed by Lees and Lin [16] who developed
            results for the generalised inﬂection point and for the eﬀects of wall heating and cooling in subsonic
            and transonic ﬂow. More recently, the stabilising eﬀect of ﬂow compressibility on laminar boundary
            layers was investigated ﬁrst by Mack [17] and then by Arnal [13, 18]. By applying the eN method with
            linear stability theory to adiabatic ﬂat plates, Mack [17] and Arnal [18] observed that compressibility
            has a strong stabilising eﬀect in transonic ﬂow, before becoming destabilising from Mach 2 to 3.5, and
            stabilising again for hypersonic ﬂow. However, applying stability analysis to three-dimensional conﬁg-
            urations, such as in the design of commercial transport aircraft, can be challenging: the framework for
            applying stability analysis is complex and can require signiﬁcant modiﬁcations to the ﬂow solver, while
            simulations often require intervention from experienced and knowledgeable users.
                The non-local AHD transition criterion [19] was developed as a simpler alternative to the eN method
            [20, 21], based on stability analysis of Falkner-Skan attached self-similar boundary-layer velocity
            proﬁles. The criterion was extended to include the stabilising eﬀect of compressibility on Tollmien-
            Schlichting instabilities based on stability analysis of compressible local similarity solutions ﬁrst up to
            Mach 1.6 [22] and, more recently, up to Mach 4 [23]. However, the AHD transition criterion requires
            non-local ﬂow quantities, which can complicate its integration in highly parallel and unstructured solu-
            tion algorithms. Two transport-equation-based transition models have recently been developed that
            create a local framework for evaluating the AHD criterion [24–26]. However, both models consist of four
            transport equations coupled to the two-equation Shear Stress Tensor (SST) turbulence model [27], result-
            ing in six-equation coupled systems. Two one-equation AHD-based transition models, each coupled to
            the SST turbulence model resulting in three-equation coupled systems, have recently been developed;
            however, they have not yet been extended to include crossﬂow instabilities [28, 29].
                Crossﬂow instabilities are dominated by the properties of the inﬂection point, which are not aﬀected
            by ﬂow compressibility to the same degree as viscous instabilities, such as Tollmien-Schlichting waves
            [30, 31]. Malik et al. [32] developed a compressibility correction for a crossﬂow Reynolds number crite-
            rion. This correction was used by Kroo and Sturdza [33] to develop a compressible crossﬂow Reynolds
            number criterion for the design of laminar supersonic swept wings, which was subsequently used by Lee
            and Jameson [34] to perform aerodynamic shape optimisation of transonic natural-laminar-ﬂow (NLF)
            wings. More recently, the C1 criterion [19], a crossﬂow criterion developed based on stability analysis
            of the solutions of the Falkner-Skan and Cooke equations, was extended to take into account compress-
            ibility eﬀects using this correction [22]. A similar correction was also used by Xu et al. [35] to extend a
            local model for stationary crossﬂow instabilities to transonic ﬂows.
                The SA-sLM2015 smooth transition model [36] is a local two-equation model which couples the
            LM2009 and LM2015 empirical correlations with the Spalart-Allmaras (SA) one-equation turbulence
            model [37], resulting in a three-equation coupled system. Smooth approximations were introduced
            to replace discontinuous and stiﬀ functions in the transition model source terms, which improves
            the iterative convergence of the model and provides continuous gradients that facilitate integration
            in gradient-based optimisation algorithms [36, 38]. The model provides a less expensive and simpler
            alternative to non-local and higher-ﬁdelity approaches described above.
                The goal of this work is to develop a transition model suitable for use in the design of commercial
            aircraft. This is achieved by extending the LM2009 [6] and LM2015 [7] empirical correlations in the
            SA-sLM2015 transition model to transonic ﬂight regimes typical of a commercial transport aircraft.

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 4 ===== -->

           4        Piotrowski and Zingg

           A compressibility correction for Tollmien-Schlichting instabilities is developed based on results from
           the compressible extension to the AHD criterion [22, 39] and the stationary crossﬂow instability com-
           pressibility correction developed by Malik et al. [32] is implemented through a new crossﬂow source
           term function. The transition model is investigated with and without the compressibility corrections
           through comparison with two- and three-dimensional transonic test cases for which transition onset has
           been measured.
              The work is presented as follows: the development of the Tollmien-Schlichting compressibility cor-
           rection and the modiﬁcations made to the transition model to integrate the stationary crossﬂow instability
           compressibility correction are presented in Section 2. The comparisons of the model with and without
           the compressibility corrections to experimental data for two- and three-dimensional transonic transition
           test cases are presented in Section 3, along with a discussion of the results. Conclusions are presented
           in Section 4.

           2.0 Compressibility corrections
           This section presents the compressibility correction developed for the Tollmien-Schlichting instability
           empirical correlation, ψ, and the framework created for applying the crossﬂow compressibility correc-
           tion developed by Malik et al. [32], ψscf . The SA-sLM2015 transition model equations are presented in
           Appendix A, while a full description of the model, including validation with two- and three-dimensional
           subsonic transition test cases, can be found in Piotrowski and Zingg [36] and Piotrowski [40]. All
           equations are presented in non-dimensional form to be consistent with previous work [36]. The non-
           dimensionalisation procedure follows that of the Navier-Stokes equations, which is described by Pulliam
           and Zingg [41]. The corrections extend the domain of applicability for these commonly used empirical
           correlations and do not aﬀect the predictive capabilities of the model in the incompressible ﬂow regime,
           as both compressibility corrections, ψ and ψscf , approach unity for low Mach numbers.

           2.1 Tollmien-Schlichting instabilities
           In an eﬀort to compare the behaviour of the LM2009 empirical correlation [6] to stability analysis,
           Perraud et al. [39] developed a simpliﬁed model based on the database of stability analysis results that
           were used to develop the compressible extension to the AHD criterion [22]. The model was developed
           for Mach numbers up to 1.1 and low turbulence intensities (Tu < 1%). Their model, which approximates
           the sensitivity of the linear stability equations to varying turbulence intensities, Tu, pressure gradients,
           λθ,e , and ﬂow compressibility, Me , is deﬁned as,
                                                                                      ρθ 2 dUe
                                          Reθt = f (Tu, λθ , Me ),        λθ,e =               Re∞ ,   Tu = Tu∞ ,   (1)
                                                                                       μ ds
                                                                                             
                           Reθt = − 177Me2 − 22Me + 210 ln((7Me + 4.8)Tu/100) exp (5Me + 27)λθ,e ,                  (2)

           where Re∞ is introduced as part of the non-dimensionalisation procedure.
              A comparison of the sensitivities of the LM2009 empirical correlation and the stability-based model
           to pressure gradient and Mach number is illustrated in Fig. 1. The LM2009 empirical correlation [6]
           is based on the empirical correlations by Abu-Ghannam et al. [42], which were developed based on
           experimental data with a focus on higher turbulence intensity conditions. In addition, the AHD criterion
           is calculated using an averaged pressure gradient parameter that is computed by integrating from the
           critical point to the location being evaluated, while the LM2009 empirical correlation is calculated using
           a fully local formulation. These factors can help to explain the diﬀerent sensitivities to pressure gradient
           in Fig. 1(a). However, as Fig. 1(b) demonstrates, the LM2009 empirical correlation under-predicts the
           transition onset momentum-thickness Reynolds number at higher Mach numbers by not including the
           stabilising eﬀect of ﬂow compressibility.

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 5 ===== -->

                                                                                            The Aeronautical Journal    5

                      (a)                                                             (b)

                                              Pressure gradient                             Mach number

            Figure 1. Sensitivity of the LM2009 empirical correlation [6], the stability-based model (Equations (1)
            and (2)), and stability analysis [39] to pressure gradient and Mach number. A higher transition onset
            momentum-thickness Reynolds number delays boundary-layer transition.

                     (a)                                                              (b)

            Figure 2. Sensitivity of the stability-based model (Equations (1) and (2)) and stability analysis [39] to
            Mach number with varying turbulence intensity and pressure gradient. The results are normalised by
            the values at a Mach number of zero to isolate the eﬀects of ﬂow compressibility.

                To better understand the stabilising eﬀect of ﬂow compressibility as predicted by linear stability
            theory, the stability-based model is investigated over a range of boundary-layer edge Mach numbers,
            Me , at various turbulence intensities, Fig. 2(a), and pressure gradients, Fig. 2(b). To isolate the eﬀects of
            ﬂow compressibility, the values of transition onset momentum-thickness Reynolds number predicted by
            both the stability-based model and stability analysis [39] are normalised by the value at a boundary-layer
            edge Mach number of zero.
                The results demonstrate that an increasing Mach number produces a more signiﬁcant stabilising eﬀect
            at low turbulence intensities and in regions with less adverse and more favourable pressure gradients
            (here corresponding to a more positive λθ,e ). This behaviour is consistent with the analysis by Ströer et al.
            [26], who investigated the compressible AHD criterion on a ﬂat plate with varying pressure gradients.
            As the Reynolds number increases, transition often occurs in regions of favourable pressure gradient.
            Therefore, the stabilising eﬀect of ﬂow compressibility is most relevant at the cruise conditions typical

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 6 ===== -->

           6        Piotrowski and Zingg

                     (a)                                                          (b)

           Figure 3. Sensitivity of the stability-based model (blue) (Equations (1) and (2)), stability analysis [39],
           and the initial compressibility correction, ψinit (grey) (Equations (3) and (4)), to Mach number with
           varying turbulence intensity and pressure gradient. The stability analysis-based results are normalised
           by the values at a Mach number of zero to isolate the eﬀects of ﬂow compressibility.

           of commercial transport aircraft, which is characterised by low turbulence intensity and high Reynolds
           number ﬂow.
              To include the stabilising eﬀect of ﬂow compressibility on the LM2009 empirical correlation for
           Tollmien-Schlichting instabilities, an initial compressibility correction, ψinit , is introduced. The correc-
           tion was designed to reproduce the trends produced by the normalised stability analysis-based results
           (Fig. 2). Speciﬁcally, it was developed to approximate the stabilising eﬀect of ﬂow compressibility on
           Tollmien-Schlichting instabilities, and its relationship with the ﬂow disturbance environment, including
           the eﬀects of turbulence intensity (Fig. 2(a)) and pressure gradient (Fig. 2(b)). The correction is evalu-
           ated using the Mach number at the edge of the boundary layer, Me , the freestream turbulence intensity,
           Tu∞ , and a local approximation to the boundary-layer edge pressure gradient, λθ,e , as follows:
                                                                                 
                            ψinit = a1 Me2 + a2 Me + a3 exp (b1 λθ,e )Me exp c1 + c2 Tu∞ Me ,                        (3)

                                a1 = 0.44, a2 = −0.38, a3 = 1.00;                 b1 = 5.00;   c1 = 0.41, c2 = −0.27.   (4)

              It is compared with the normalised stability analysis-based results in Fig. 3, with the results presented
           in grey. The results demonstrate that the correction does a reasonable job of reproducing the trends
           predicted by stability analysis.
              The correction was initially applied to the local transition onset momentum-thickness Reynolds
           number, Reθt , which is convected and diﬀused into the boundary layer by the R̃eθt transport equation.
           However, simulations using the test cases presented in Section 3 demonstrated that this approach sig-
           niﬁcantly over-predicted the laminar extent of the boundary layer for the CAST10-2 and VA-2 aerofoils.
           The rapid acceleration at the upper surface leading edge produced large upstream values of the com-
           pressibility correction, and therefore R̃eθt . These large values convect downstream in the boundary layer
           and delay transition. This behaviour is not reproduced by methods based on the AHD criterion where
           the pressure gradient parameter is integrated from the critical point.
              Applying the correction to the critical momentum-thickness Reynolds number, Reθc (Equation
           (A.33)), produced transition locations that agree better with the CAST10-2 and VA-2 experimental
           results by avoiding this behaviour, while still agreeing well with the CRM-NLF test case. This approach
           is similar to scaling the local value of the transported momentum-thickness Reynolds number, R̃eθt ,

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 7 ===== -->

                                                                                                               The Aeronautical Journal    7

                      (a)                                                             (b)

            Figure 4. Sensitivity of the initial, ψinit (grey) (Equations (3) and (4)), and modiﬁed, ψ (red)
            (Equations 5–7), compressibility corrections to Mach number with varying turbulence intensity and
            pressure gradient.

            as Reθc is a near-linear function of R̃eθt . A consequence of this approach is that the Flength empirical cor-
            relation (Equation (A.34)) is not scaled by the compressibility correction. However, the stability-based
            model and the AHD criterion do not predict the eﬀects of compressibility on the transition region length.
            For example, the AHD-based transition models developed by Pascal et al. [24] and Ströer et al. [25, 26]
            rely on user-speciﬁed inputs for the length of the transition region. Further work is required to develop
            a correlation for the transition length in compressible ﬂows.
               Using this approach, simulations revealed that the smooth transition model with the initial com-
            pressibility correction, ψinit , still over-predicted the extent of the laminar boundary layer for simulations
            of the VA-2 aerofoil at lower angles of attack, and appeared to over-estimate the eﬀect of favourable
            and adverse pressure gradients on the stabilising eﬀect of compressibility. This could again be related
            to diﬀerences between the local model and the AHD criterion, such as the use of a fully local pres-
            sure gradient parameter versus an integrated value. Therefore, the correction was modiﬁed to improve
            agreement with the experimental results. Speciﬁcally, the a1 and b1 constants were reduced relative to
            the initial correction in order to reduce the sensitivity to Mach number and pressure gradient, respec-
            tively. Care was taken to ensure that these modiﬁcations did not have a signiﬁcant adverse eﬀect on the
            simulated CAST10-2 and CRM-NLF transition fronts. The modiﬁed compressibility correction, ψ, pro-
            duces transition locations that agree well with experimental results for the test cases and ﬂow conditions
            investigated (0.71 ≤ M ≤ 0.856, 2 × 106 ≤ Re ≤ 15 × 106 ), and is given by,
                                                                      Reθc,comp = ψReθc ,                                                 (5)
                                                                             
                                   ψ = a1 Me2 + a2 Me + a3 exp b1 λθ,e Me exp c1 + c2 Tu∞ Me ,                                            (6)

                                 a1 = 0.34, a2 = −0.38, a3 = 1.00;                    b1 = 3.00;   c1 = 0.41, c2 = −0.27.                 (7)

               The correction is compared with the initial correction, ψinit (Equations (3) and (4)), in Fig. 4.
               The compressibility corrections are evaluated using the Mach number at the edge of the boundary
            layer, Me , which can be approximated locally using isentropic relations and the compressible Bernoulli
            equation as follows [43]:
                                                                                                   1− κ1
                                                              Ue                        κp∞ p
                                                          Me = ,           ae =                            ,                              (8)
                                                              ae                        ρ∞ p∞

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 8 ===== -->

           8        Piotrowski and Zingg

                                                                                                          
                                                         
                                                               2κ   p       p
                                                                                                    1− κ1

                                                    Ue = U∞
                                                           2 +        ∞
                                                                         1−                                     ,    (9)
                                                               κ − 1 ρ∞     p∞

           where κ is the heat capacity ratio.
              The pressure gradient parameter, λθ (Equations (A.13)–(A.17)), used in the LM2009 empirical cor-
           relation is based on the local streamwise velocity gradient formed using Cartesian velocity gradients,
           which is not valid in the boundary layer. To avoid problems with this formulation, the LM2009 transition
           model introduces the Fθt function (Equation (A.4)) in the Pθt source term (Equation (A.2)) to disable
           the Tollmien-Schlichting empirical correlation inside the boundary layer. However, the compressibility
           corrections are most active in the middle of the boundary layer where the transition onset process begins.
           To address this, the pressure gradient parameter used to evaluate the compressibility correction is cal-
           culated using a local approximation to the boundary-layer edge pressure gradient, λθ,e , adopted from
           Grabe et al. [44, 45], which is approximated using Cartesian pressure gradients and isentropic relations
           as follows:
                                                             ρθ 2 dUe
                                                      λθ,e =          Re∞ ,                                      (10)
                                                              μ ds
                                                                                    
                                        dUe       u dUe  v dUe  w dUe
                                             =              +            +              ,                        (11)
                                         ds       U dx          U dy         U dz
                                                                                      − κ1
                                                              dUe      1    p                dp
                                                                  =−                            ,                   (12)
                                                               dx    ρ∞ Ue p∞                dx
                                                                                      − κ1
                                                              dUe      1    p                dp
                                                                  =−                            ,                   (13)
                                                               dy    ρ∞ Ue p∞                dy
                                                                                      − κ1
                                                              dUe      1    p                dp
                                                                  =−                            .                   (14)
                                                               dz    ρ∞ Ue p∞                dz

               Similar to the original pressure gradient parameter, λθ (Equations (A.13)–(A.17)), λθ,e is not Galilean
           invariant.
               It is important to note that the boundary-layer edge pressure gradient parameter, λθ,e , deﬁned above
           is not used to calculate the LM2009 empirical correlation (Equations (A.8)–(A.12)). Although the two
           pressure gradient formulations produce similar results for favourable and adverse pressure gradients,
           the original, velocity-gradient-based pressure gradient parameter, λθ , is used to be consistent with the
           Langtry-Menter model [6], as it provides better agreement with experimental data for zero-pressure-
           gradient ﬂat plate cases. This is likely due to the calibration of the model. For numerical robustness, λθ,e
           is limited using the bounds introduced by Langtry and Menter [6] for λθ (−0.1 ≤ λθ,e ≤ 0.1).

           2.2 Stationary crossflow instabilities
           The stationary crossﬂow transition onset momentum-thickness Reynolds number correlation, Reθt,scf ,
           produced by the LM2015 empirical correlation (Equations (A.18)–(A.23)) was originally implemented
           using an additional source term in the R̃eθt transport equation, Dscf [7]. The LM2009 empirical corre-
           lation inﬂuences R̃eθt in the freestream, which is then convected and diﬀused into the boundary layer
           where the stationary crossﬂow source term is active. The combined eﬀect of these correlations is realised
           through the Fonset (Equation (A.31)) function in the intermittency source terms through the critical
           momentum-thickness Reynolds number, Reθc .
               In order to include the stabilising eﬀects of compressibility on each set of correlations, and sub-
           sequent Reynolds numbers, Reθt and Reθt,scf , separately, Fonset is separated into two functions, one for
           two-dimensional mechanisms and one for stationary crossﬂow instabilities. A similar approach to that
https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 9 ===== -->

                                                                                                        The Aeronautical Journal     9

            developed by Grabe et al. [45] and Carnes and Coder [46] is adopted. A new onset function, Fonset,scf , is
            introduced for the intermittency transport equation, replacing the Dscf source term in the R̃eθt transport
            equation. The Fonset,scf function was developed using the same approach as the original Fonset function,
            and therefore the smooth approximation introduced in previous work [36] (Equation (A.31)), leading to
            the following:
                                                           tanh (6(Fonset,scf,1 − 1.35)) + 1
                                              Fonset,scf =                                   ,                    (15)
                                                           ReS
                                        Fonset,scf,1 =               + (RT )2 , Reθc,scf = 0.623Reθt,scf ,               (16)
                                                      2.60Reθc,scf
            where ReS and RT represent the strain-rate magnitude Reynolds number and eddy viscosity ratio, respec-
            tively. The transition process is triggered when Fonset,scf,1 exceeds unity and is sustained as the eddy
            viscosity ratio, RT , increases. The stationary crossﬂow critical momentum-thickness Reynolds number,
            Reθc,scf , correlation in the new Fonset,scf function was calibrated using the NLF2-0415 Inﬁnite Swept Wing
            [47, 48] transition test case. The value of 0.623 applied to Reθt,scf in order to initiate the transition process
            upstream of the predicted transition location closely resembles a linear ﬁt of the smooth Reθc function
            (Equation (A.33)) and is consistent with the value of 0.62 used by Medida [49].
               The Fonset,scf function is incorporated in the intermittency source terms (Equations (A.29) and (A.30))
            by taking a smooth maximum of the two onset functions:
                                                                Fonset ← φ300 (Fonset , Fonset,scf ),                              (17)
            where φ±300 represent smooth approximations of minimum/maximum operators introduced in previ-
            ous work [36]. A positive value of the smoothing parameter reproduces a maximum, while a negative
            value approximates a minimum, and a larger magnitude value more closely reproduces the mini-
            mum/maximum operator. A value of ±300 is used for the results presented. The vorticity limiting
            procedure in the intermittency source terms introduced in previous work to improve the iterative con-
            vergence of the model [36] is re-written in Appendix A to clarify the implementation in dimensional
            solvers. This formulation is equivalent to that presented by Piotrowski and Zingg [36] which was
            non-dimensionalised following the procedure for the Navier-Stokes equations(41), where velocity and
            distance are normalised by the sound speed and reference length, l, respectively. The new crossﬂow
            source term formulation was validated using the TU Braunschweig Sickle Wing transition test case
            [50], with the results presented by Piotrowski [40].
                Simulations of transonic swept wings revealed that both crossﬂow source-term approaches, Dscf and
            Fonset,scf , produce a transition front that is dependent on the initialisation of the ﬂow ﬁeld. A solution ini-
            tialised with far-ﬁeld conditions converges to a transition front upstream of the converged transition front
            produced when initialised with a converged solution obtained with the crossﬂow correlation inactive. It
            is important to note that this behaviour does not appear in subsonic transition test cases with crossﬂow,
            such as the NASA NLF2-0415 Inﬁnite Swept Wing [47, 48] and the TU Braunschweig Sickle Wing [50].
            To prevent these inaccurate upstream transition fronts, for transonic cases the crossﬂow source term is
            activated after the total residual drops several orders of magnitude without the crossﬂow source term,
            Fonset,scf , active. A relative residual drop tolerance of ﬁve orders of magnitude for the total residual was
            found to be suﬃcient and is used for the results presented in the current work. This strategy produces the
            same results as initialising with a converged simulation performed without the crossﬂow correlations
            active.
                The crossﬂow Reynolds number compressibility correction developed by Malik et al. [32], ψscf ,
            is applied by scaling the incompressible stationary crossﬂow transition onset momentum-thickness
            Reynolds number, Reθt,scf , in the Fonset,scf function as follows:
                                                                   Reθt,scf,comp = ψscf Reθt,scf ,                                 (18)
                                                                                      κ −1 2
                                                                    ψscf = 1 +            Me .                                     (19)

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 10 ===== -->

           10        Piotrowski and Zingg

           Figure 5. Eﬀects of Mach number on the Tollmien-Schlichting and stationary crossﬂow instability
           compressibility corrections, ψ and ψscf , respectively.

              The compressibility corrections for Tollmien-Schlichting and stationary crossﬂow instabilities, ψ
           and ψscf , respectively, are presented in Fig. 5 for Tu∞ = 0.05% and λθ,e = 0.04. As expected, the com-
           pressibility correction for Tollmien-Schlichting instabilities produces a stronger stabilising eﬀect than
           the crossﬂow instability compressibility correction, especially in a favourable ﬂow environment.

           3.0 Results
           The transition model with and without the compressibility corrections, SA-sLM2015cc and SA-
           sLM2015, respectively, is applied to two- and three-dimensional transonic transition test cases. These
           cases consist of the CAST10-2 and VA-2 aerofoil test cases and the NASA CRM-NLF wing-body
           geometry, with the Reynolds number increasing for each successive test case. As the Reynolds num-
           ber increases and transition occurs in regions with less adverse and more favourable pressure gradients,
           the stabilising eﬀect of ﬂow compressibility on Tollmien-Schlichting instabilities is expected to become
           more signiﬁcant (see Fig. 3(b)).
               A second-order discretisation is used for the mean-ﬂow equations using a matrix-based dissipation
           model [51], with a ﬁrst-order upwind scheme applied to the turbulence and transition model convec-
           tive terms. Care is taken for each test case to ensure that all results are suﬃciently iteratively and grid
           converged to produce the low levels of numerical error required to investigate the modelling error and
           validity of the transition model variants. Machine-zero residual convergence for the total residual is
           achieved for the cases presented unless stated otherwise. More details on the solution strategy and the
           iterative and grid convergence of the smooth transition model with the compressibility corrections can
           be found in previous work [40, 52].

           3.1 CAST10-2 aerofoil
           The CAST10-2 aerofoil [53] was investigated by Hebler et al. [54] in the Transonic Wind Tunnel
           Göttingen (DNW-TWG) at a Mach number of 0.75, Reynolds number of 2 × 106 , and angles of attack
           ranging from −0.79 to 1.41 degrees. The aerofoil has also been investigated numerically by several
           researchers [3, 26, 55, 56]. The turbulence intensity of the DNW-TWG wind tunnel was not reported
           in the experiment. However, Fehrs et al. [55] suggest a range of turbulence intensities for the wind tun-
           nel of between 0.25–0.40%, while Ströer et al. [26] assumed a value of 0.22%. Hebler et al. simulated
           the aerofoil with the 2D coupled Euler-boundary-layer solver MSES [57], which is capable of predict-
           ing transition using a database eN method. The results demonstrate that an N-factor of 6 produced a

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 11 ===== -->

                                                                                                    The Aeronautical Journal    11

                                                 Table 1. CAST10-2 structured O-grid dimensions
                 Grid level             Chord × oﬀ-wall nodes                     Avg/max     s × 10−6 (chord)     Avg/max y+
                 L0                           541 × 121                                     5.06/5.52               0.27/0.84
                 L1                           761 × 171                                     3.58/3.82               0.19/0.59
                 L2                         1,081 × 241                                     2.52/2.63               0.13/0.42
                 L3                         1,521 × 341                                     1.79/1.86               0.11/0.32

                     (a)                                                              (b)

            Figure 6. Grid-convergence results for the CAST10-2 aerofoil simulations at −0.39 and 0.82 degrees
            angle-of-attack and M = 0.74, Re = 2 × 106 , and Tu = 0.25%.

            transition front that agrees well with the experiment, which, using Mack’s relation [58], corresponds to
            a turbulence intensity of approximately 0.24% [54].
                The blunt trailing-edge aerofoil coordinates are provided by Dress et al. [59]. As discussed by Fehrs
            et al. [55], the trailing edge of the aerofoil is located below the reference line, with the aerofoil coordi-
            nates exhibiting a 0.88◦ angle-of-attack. All angles of attack are provided with respect to this reference
            system. A Mach number shift of −0.01 and angle-of-attack shift of −0.30 degrees were identiﬁed
            by Hebler et al. [54] in order to match the experimental data with free-air simulations. These shifts
            are applied in the current work, along with an assumed turbulence intensity of 0.25% based on the
            lower bound estimated by Fehrs et al. [55]. However, the results are plotted and referenced using the
            uncorrected angles of attack in order to be consistent with the results presented by Hebler et al. [54].
                Four grid levels are investigated using structured multi-block grids with the characteristics presented
            in Table 1. A grid-reﬁnement study is presented using the transition model both with and without the
            compressibility corrections at the −0.39 and 0.82 degree angles of attack, with the lift and drag grid
            convergence illustrated in Fig. 6. The largest diﬀerence in the Cl and Cd between the second-ﬁnest and
            ﬁnest grid levels is less than 0.30% for both transition model variants at each angle-of-attack.
                An angle-of-attack sweep is presented for the three ﬁnest grid levels both fully turbulent using the SA
            turbulence model and with free transition using the transition model with and without the compressibility
            corrections. The pressure and upper-surface skin friction coeﬃcient proﬁles at angles of attack of −0.39,
            0.05, and 0.82 degrees are overlaid with the pressure proﬁles from the experiment in Fig. 7.
                The results demonstrate that the pressure and skin friction proﬁles are suﬃciently grid converged. The
            free-transition pressure coeﬃcient proﬁles demonstrate good agreement with the experiment at −0.39
            and 0.82 degrees angle-of-attack. However, for the 0.05 degree angle-of-attack case, in the nonlinear
            lift region (see Fig. 8), the fully turbulent pressure proﬁles appear to provide better agreement with the

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 12 ===== -->

           12        Piotrowski and Zingg

                                                                                (a)

                     (b)                                                              (c)

           Figure 7. Pressure and upper-surface skin friction coeﬃcient proﬁles for the CAST10-2 aerofoil pro-
           duced at M = 0.74, Re = 2 × 106 , and Tu = 0.25% at three angles of attack overlaid with the pressure
           proﬁles from the experiment [54].

           experiment. In general, the wind-tunnel corrections suggested by Hebler et al. [54] appear to reproduce
           the wind-tunnel environment well.
               As expected, the Tollmien-Schlichting instability compressibility correction does not have a signiﬁ-
           cant eﬀect on the transition onset locations for this low Reynolds number case, as demonstrated by the
           upper-surface skin friction coeﬃcient proﬁles, as transition primarily occurs either due to a laminar sep-
           aration bubble or a strong adverse pressure gradient. Figure 3(b) demonstrates that ﬂow compressibility
           is not predicted to have a strong stabilising eﬀect in these ﬂow environments. This behaviour is also
           demonstrated in the work by Ströer et al. [26], who simulated the CAST10-2 aerofoil at similar ﬂight
           conditions using a transport-equation-based formulation of the compressible AHD criterion.
               The predicted transition locations on the upper surface of the CAST10-2 aerofoil and lift curves,
           both obtained on the L3 grid, are presented in Fig. 8 and are compared with the experimental results
           from Hebler et al. [54]. The results demonstrate that the free-transition simulations do a reasonable job
           of predicting the nonlinear transition and lift curves produced by the experiment. However, at angles
           of attack between −0.17 and 0.64 degrees, the free-transition simulations under-predict the drop in the
           lift curve by over-predicting the extent of the laminar boundary layer relative to experimental data. This
           could be a result of the uncertainty in the turbulence intensity for the wind tunnel at these conditions.
               As illustrated in Fig. 8, the CAST10-2 aerofoil produces highly nonlinear transition and lift curves,
           with transition moving forward signiﬁcantly on the upper surface of the aerofoil over a small range of
           angles of attack. For the 0.05 degree angle-of-attack case, steady-state simulations using the transition
           model with the compressibility corrections produced an oscillatory force history, with the upper-surface

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 13 ===== -->

                                                                                      The Aeronautical Journal        13

            Figure 8. Upper-surface transition locations (dashed) and lift curves (solid) for the CAST10-2 aerofoil
            obtained on the L3 grid at M = 0.74, Re = 2 × 106 , and Tu = 0.25% over a range of angles of attack
            compared with the results from the experiment [54].

            transition front moving upstream and downstream in a periodic manner. The red areas in Fig. 7 repre-
            sent the regions between the pressure and skin friction coeﬃcient proﬁles for the most upstream and
            downstream transition locations on each grid level, while the error bars in Fig. 8 represent the range in
            the transition location and lift coeﬃcient over one cycle of this oscillatory behaviour.
                The results in Figs 7 and 8 demonstrate the importance of including the relationship between pres-
            sure gradient and ﬂow compressibility in the Tollmien-Schlichting instability compressibility correction.
            Due to the low Reynolds number, transition occurs in regions of adverse pressure gradient where sta-
            bility analysis does not predict the eﬀects of ﬂow compressibility to be signiﬁcant. A compressibility
            correction that does not take into account the eﬀects of pressure gradient, such as that investigated by
            Venkatachari et al. [60], will over-predict the stabilising eﬀect of ﬂow compressibility at these ﬂow con-
            ditions, resulting in predicted transition and lift curves that do not capture the nonlinear trends produced
            by Hebler et al. [54]. As the Reynolds number increases, the stabilising eﬀect of ﬂow compressibility
            becomes more signiﬁcant. This is demonstrated in the following test cases.

            3.2 VA-2 supercritical aerofoil
            The VA-2 supercritical aerofoil was recently investigated in the DNW-TWG wind tunnel by Costantini
            et al. [61]. The aerofoil was examined at a Mach number of 0.72, Reynolds number of 10 × 106 , and
            angles of attack from −0.4 to 2.0 degrees. The upper-surface skin friction coeﬃcient distributions and
            transition locations were determined using a global luminescent oil-ﬁlm skin friction ﬁeld estimation
            method (GLOFSFE), which measures the development of the thickness of an oil ﬁlm, the distribution
            of which can be used to calculate the skin friction, based on its luminescent intensity. The pressure
            distributions used for comparison in the current work were obtained before applying the oil ﬁlm.
                Two transition locations were identiﬁed in the experiment: the transition onset location, identiﬁed
            where the skin friction coeﬃcient increased beyond a value of 4 × 10−4 , and the end of the transition
            region, deﬁned as the location where the skin friction coeﬃcient reached a value of 3 × 10−3 . While the
            authors state that the end of the transition region can only be considered as qualitative due to uncertain-
            ties in the turbulent skin friction estimation originating from thick oil-ﬁlm distributions in these regions,

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 14 ===== -->

           14        Piotrowski and Zingg

                                                   Table 2. VA-2 structured O-grid dimensions
              Grid level               Chord × oﬀ-wall nodes                      Avg/max     s × 10−6 (chord)   Avg/max y+
              L0                             541 × 121                                      1.01/1.07             0.25/0.49
              L1                             761 × 171                                      0.71/0.75             0.18/0.35
              L2                           1,081 × 241                                      0.50/0.52             0.13/0.25
              L3                           1,521 × 341                                      0.36/0.37             0.09/0.19

                     (a)                                                              (b)

           Figure 9. Grid-convergence results for the VA-2 aerofoil simulations at −0.40 and 1.80 degrees
           angle-of-attack and M = 0.71, Re = 10 × 106 , and Tu = 0.25%.

           quantitative results were obtained for transition onset locations with lower uncertainty. Transition loca-
           tions were also determined using an automatic temperature-sensitive paint (TSP) transition detection
           method, where the recorded transition locations are expected to be between the transition onset and end
           locations [61, 62].
              The Mach number and angle-of-attack shifts identiﬁed by Hebler et al. for the DNW-TWG wind tun-
           nel [54], −0.01 Mach and −0.30 degrees, respectively, are applied for the VA-2 simulations. Similar
           to the previous case, the results are presented and referenced using the uncorrected angles of attack.
           The turbulence intensity for the VA-2 experimental investigations in the DNW-TWG was not pro-
           vided. Therefore, the turbulence intensity is assumed to be consistent with the CAST10-2 DNW-TWG
           conditions, with a value of 0.25% used in the current work.
              Four grid levels are investigated using structured multi-block grids with the characteristics presented
           in Table 2. Grid-reﬁnement studies are presented both with and without the compressibility corrections
           at angles of attack of −0.40 and 1.80 degrees in Fig. 9. The diﬀerence in the drag coeﬃcient between
           the ﬁnest grid level and the grid-converged value was evaluated using Richardson extrapolation for each
           transition model variant at both angles of attack. For the −0.40 degree case, the diﬀerences are 0.18%
           and 0.80% for the transition model with and without the compressibility corrections, respectively, while
           for the 1.80 degree case the diﬀerences are 0.24% and 0.37%. The iterative and grid convergence for
           this test case is presented in more detail in previous work [52].
              An angle-of-attack sweep is presented for the three ﬁnest grid levels both with and without the com-
           pressibility corrections. The pressure and upper-surface skin friction coeﬃcient proﬁles are overlaid
           with the proﬁles from the experiment [61] in Fig. 10. For the lower angles of attack, −0.40 degrees
           to 1.20 degrees, there is good agreement between the simulated pressure coeﬃcient proﬁles and the
           proﬁles from the experiment. However, Costantini et al. observed that the algorithm used to control the
           DNW-TWG adaptive wind tunnel walls, which are used to reduce the eﬀect of the walls on the pressure

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 15 ===== -->

                                                                                      The Aeronautical Journal    15

            Figure 10. Pressure and upper-surface skin friction coeﬃcient proﬁles for the VA-2 aerofoil produced
            at M = 0.71, Re = 10 × 106 , and Tu = 0.25% over a range of angles of attack overlaid with the results
            from the experiment [61].

            distributions, failed to converge at angles of attack above 1.20 degrees [61], which helps to explain the
            diﬀerence in the simulated and experimental pressure proﬁles observed at angles of attack of 1.50 and
            2.00 degrees. The results conﬁrm that the Mach number and angle-of-attack corrections applied in the
            simulations accurately reproduce the wind-tunnel environment.
               The simulated upper-surface skin friction coeﬃcient proﬁles are plotted and compared with the pro-
            ﬁles obtained using the GLOFSFE method [61]. The experimental skin friction coeﬃcient proﬁles are

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 16 ===== -->

           16        Piotrowski and Zingg

                                                                 Figure 10. Continued.

           plotted both with the raw data, which is clipped at a value of 6 × 10−3 , and with the data scaled by a
           factor of 17. Costantini et al. note that the oil-ﬁlm thickness in the turbulent boundary-layer regions
           was larger than the height of the viscous sublayer, which produced a hydraulically rough surface that
           increased the skin friction relative to the clean conﬁguration [61]. The clipped experimental data pro-
           vides a quantitative comparison of the laminar skin friction coeﬃcient proﬁles and the transition onset
           locations, which were not as aﬀected by the oil ﬁlm due to the thinner oil-ﬁlm thickness in these regions,
           while the scaled data allows for a qualitative comparison of the skin friction coeﬃcient proﬁles in the
           transition and turbulent regions. The clipped experimental skin friction coeﬃcient proﬁles agree well
           with the simulated proﬁles in the laminar boundary layer.
              At the lower angles of attack, below 0.80 degrees, there is a weak favourable pressure gradient on
           the upper surface of the aerofoil that develops into a weak adverse pressure gradient as the angle-of-
           attack increases. In these regions, where Fig. 3(b) demonstrates that ﬂow compressibility is expected
           to be stabilising, the compressibility correction pushes transition aft to better agree with the scaled
           GLOFSFE-estimated skin friction coeﬃcient proﬁles [61]. As the angle-of-attack increases and the
           ﬂow decelerates on the upper surface, the stabilising eﬀect of ﬂow compressibility is reduced, and the
           transition front moves upstream to match both the transition model without the compressibility correc-
           tions and the experimental data. Above an angle-of-attack of 1.20 degrees, the upper-surface pressure
           coeﬃcient proﬁles plateau and the strong adverse pressure gradient moves aft in the form of a shock
           wave. The stabilising eﬀect of compressibility increases with increasing angle-of-attack as the adverse
           pressure gradient weakens, pushing transition downstream. Although wind-tunnel eﬀects at the 1.50 and

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 17 ===== -->

                                                                                               The Aeronautical Journal            17

            Figure 11. Upper-surface transition locations for the VA-2 aerofoil obtained on the L3 grid at M = 0.71,
            Re = 10 × 106 , and Tu = 0.25% over a range of angles of attack compared with results from the
            experiments [61].

            2.00 degree angles of attack appear to be signiﬁcant, good agreement is achieved for both the pressure
            and skin friction coeﬃcient proﬁles at 1.80 degrees.
                The predicted transition locations on the L3 grid are compared with values from the experiments
            performed using the GLOFSFE and TSP methods [61] in Fig. 11. The results demonstrate that the tran-
            sition model without the compressibility corrections produces a ﬂatter transition curve that is upstream
            of the experimental GLOFSFE-estimated transition onset locations for all angles of attack except −0.40
            and 1.20 degrees. However, at −0.40 degrees the transition location produced without the compressibil-
            ity corrections is upstream of the location produced by the TSP method. The TSP transition locations
            should lie between the GLOFSFE-estimated transition onset and end locations [61], such as for the
            −0.40, 0.80, and 1.80 degrees angle-of-attack. There is increased uncertainty where this is not the case,
            such as at the 0.00 and 0.40 degree angles of attack. In addition, at 1.20 and 1.50 degrees the experimen-
            tal upper-surface shocks are upstream of the simulated results, which explains the upstream transition
            locations. The Tollmien-Schlichting instability compressibility correction delays the predicted transition
            onset locations to agree best with the GLOFSFE-estimated end of transition region locations, except for
            1.80 and 2.00 degrees angle-of-attack, where they are closer to the GLOFSE-estimated transition onset
            locations.

            3.3 NASA CRM-NLF wing-body geometry
            The NASA CRM-NLF conﬁguration was recently investigated as part of the First AIAA Transition
            Modelling and Prediction Workshop∗ . Transition visualisations, pressure coeﬃcient proﬁles, and inte-
            grated forces were provided from the experiment over a range of angles of attack [11, 12]. The test
            conditions for the data provided are presented in Table 3. The turbulence intensity in the NASA
            Langley National Transonic Facility (NTF) wind tunnel was determined through comparison with results
            obtained using stability analysis [63]. The Tollmien-Schlichting N-factor was found to vary from 4 to 8
            [63], with a critical N-factor of 6 assumed for the workshop. Using Mack’s relation [58], this corresponds
            to a turbulence intensity of approximately 0.24%, which is the value used in the current work. Surface
              ∗ https://transitionmodeling.larc.nasa.gov/wp-content/uploads/sites/109/2020/02/TransitionMPWCaseDescriptions.pdf,

            accessed June 2021

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 18 ===== -->

           18        Piotrowski and Zingg

                                                Table 3. CRM-NLF wind tunnel test conditions∗
                  Test case            Angle-of-attack (◦ )               Mach number   Reynolds number ×106 (MAC)
                  2523                     1.44848                         0.856489               14.97197
                  2524                     1.98031                         0.856491               14.94591
                  2525                     2.46141                         0.856051               14.90909
                  2526                     2.93787                         0.855801               14.85308

                                              Table 4. CRM-NLF structured grid characteristics†
                 Grid level         # of nodes           Average/maximum s × 10−6 (chord)          Average/maximum y+
                 L0                  8,893,456                       1.00/2.24                           0.36/1.50
                 L1                 16,691,200                       0.78/1.75                           0.28/1.06
                 L2                 32,787,200                       0.60/1.35                           0.22/0.78
                 L3                 64,330,000                       0.47/1.06                           0.18/0.64

                     (a)                                     (b)                        (c)

           Figure 12. Grid-convergence results for the CRM-NLF simulations at the 2524 test conditions
           (α ≈ 2.0◦ ).

           roughness has a signiﬁcant eﬀect on crossﬂow instabilities [47, 48]. To prevent roughness elements from
           destabilising the boundary layer, the surface of the CRM-NLF wind tunnel model was frequently sanded
           and polished, with the average size of the roughness elements after testing measured to vary from 0.83
           to 1.10μin [11]. A value of 1.00μin is assumed in the current work.
              Four grid levels are investigated using structured multi-block grids following the gridding guidelines
           provided by the workshop committee† , with the grid characteristics listed in Table 4. A grid-reﬁnement
           study is presented at the 2524 (≈ 2.0◦ ) test conditions (see Table 3) using three modelling strategies: fully
           turbulent using the SA turbulence model [37, 64] with the QCR2000 correction [65], and free transition
           using the smooth transition model both with, SA-QCR2000-sLM2015cc, and without, SA-QCR2000-
           sLM2015, the Tollmien-Schlichting and stationary crossﬂow instability compressibility corrections. The
           original Dscf crossﬂow instability source term (Equation (A.3)) is used for the cases without the com-
           pressibility corrections, with the new Fonset,scf source term presented in Section 2.2 (Equations (15) and
           (16)) used for the cases with the compressibility corrections.
              The grid-reﬁnement study results at the 2524 (≈ 2.0◦ ) test conditions obtained with free transition
           using both transition model variants are presented in Fig. 12. The diﬀerence in the drag coeﬃcient on the
           ﬁnest grid level relative to the grid-converged value was calculated for the simulations with and without
           the compressibility corrections to be 0.65% and 0.66%, respectively. The iterative and grid convergence
           for this test case is also presented in more detail in previous work [52].
              The residual convergence histories, illustrated by the norm of the total residual normalised by the
           value at the ﬁrst iteration, for the free-transition simulations are presented in Fig. 13. The residuals are
             †
             https://transitionmodeling.larc.nasa.gov/wp-content/uploads/sites/109/2020/02/CRM-NLF_GridGuidelines.pdf,   accessed
           June 2021.

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 19 ===== -->

                                                                                      The Aeronautical Journal      19

            Figure 13. Grid-reﬁnement study residual convergence histories for the CRM-NLF simulations at the
            2524 test conditions (α ≈ 2.0◦ ).

            plotted against equivalent residual evaluations, which is the total wallclock time normalised by the aver-
            age cost for computing a residual evaluation. It is important to note that the same solver parameters are
            used for both transition model variants on all four grids. The results demonstrate that the compressibility
            corrections and new Fonset,scf source term do not signiﬁcantly aﬀect the iterative and grid convergence of
            the model.
                The pressure and skin friction coeﬃcient proﬁles for the three ﬁnest grids are extracted at nine span-
            wise stations across the wing and compared with the experimental proﬁles in Fig. 14. The upper-surface
            shocks help to diﬀerentiate the upper and lower surface skin friction coeﬃcient proﬁles, as the shock
            causes a local reduction in the upper-surface skin friction coeﬃcient due to the rapid ﬂow deceleration.
            The compressibility corrections delay transition on the lower surface, which is dominated by crossﬂow
            instabilities, and on the upper surface in regions with less adverse and more favourable pressure gradi-
            ents. The upper-surface transition locations appear suﬃciently grid converged. However, at the η = 0.910
            station the lower-surface transition front produced by the transition model with the compressibility
            corrections moves downstream from the L1 to the L2 grid level.
                The pressure coeﬃcient proﬁles produced by the fully turbulent simulations provide better agreement
            with the experiment at the inboard and midspan sections (η < 0.640) relative to the simulations with
            free transition. This counter-intuitive behaviour can be explained by the TSP images by Lynde et al.
            [11, 12], which demonstrate a signiﬁcant amount of bypass transition produced by surface imperfections,
            including the leading-edge pressure ports. Lynde et al. [11] attempted to estimate the natural transition
            front behind these bypass transition-induced turbulent wedges. As the simulated transition front moves
            downstream to better match this estimated natural transition front, a worse agreement is seen with the
            experimental pressure coeﬃcient proﬁles. This behaviour was investigated in detail by Helm et al. [66].
            At the further outboard stations, where there is less bypass transition observed in the experiment due
            to the reduced chord Reynolds number, the increased laminar extent of the boundary layer produced by
            the transition model with the compressibility corrections moves the shock locations aft to better agree
            with experimental data.
                The free-transition skin friction coeﬃcient proﬁles on the upper surface of the wing obtained on the
            three ﬁnest grid levels are overlaid with the estimated natural transition front from Lynde et al. [11] in
            Fig. 15. The compressibility corrections successfully push the transition front aft to better agree with
            the experimental results, particularly in the outboard regions. The results again demonstrate that the
            upper-surface transition front is suﬃciently grid converged. The wing lower surface and fuselage nose
            boundary layers were tripped in the experiment. Therefore, transition fronts on the fuselage and lower
            surface of the wing are not available for comparison.

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 20 ===== -->

           20        Piotrowski and Zingg

                    (a)                                       (b)                     (c)

                    (d)                                       (e)                     (f)

                    (g)                                       (h)                     (i)

           Figure 14. Pressure and skin friction coeﬃcient proﬁles for the CRM-NLF grid-reﬁnement study at the
           2524 test conditions (α ≈ 2.0◦ ) compared with the pressure proﬁles from the experiment [11] at varying
           spanwise stations η.

               Skin friction coeﬃcient proﬁles produced by free-transition simulations at the 2523 (≈ 1.5◦ ), 2525
           (≈ 2.5◦ ), and 2526 (≈ 3.0◦ ) test conditions obtained on the L1 grid are illustrated in Fig. 16. Again, the
           compressibility corrections delay transition, producing transition fronts that better agree with exper-
           imental data. Regions of the transition fronts produced with the compressibility corrections remain
           upstream of the experiment. This is particularly evident for the 2523 (≈ 1.5◦ ) test case, where transition
           is upstream of experimental data inboard near the Yehudi break, and at approximately mid-span at the
           2525 (≈ 2.5◦ ) and 2526 (≈ 3.0◦ ) conditions. However, this behaviour is also present to a lesser extent
           in the work of Paredes et al. [67], Venkatachari et al. [60], and Krimmelbein and Krumbein [68], who
           each investigate the CRM-NLF conﬁguration using a variety of stability analysis methods. These diﬀer-
           ences could be due to uncertainties in the turbulence intensity, ﬂow conditions, or aeroelastic deﬂection
           of the wind-tunnel model. In general, there is better agreement between the transition model with the
           compressibility corrections and the experimental data over the outboard portion of the wing.
               The residual convergence histories for the free-transition simulations are presented in Fig. 17. The
           same solver parameters are used for both transition model variants at all three angles of attack, which

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 21 ===== -->

                                                                                            The Aeronautical Journal   21

                     (a)                                                              (b)

                     (c)                                                              (d)

                     (e)                                                              (f)

            Figure 15. Upper-surface skin friction coeﬃcient proﬁles for the CRM-NLF grid-reﬁnement study
            at the 2524 test conditions (α ≈ 2.0◦ ) overlaid with the estimated natural transition front from the
            experiment [11].

            also match the solver settings used in the grid-reﬁnement study at the 2524 test conditions shown in
            Fig. 13. The results again demonstrate similar iterative convergence between the two transition model
            variants.
               The pressure and skin friction coeﬃcient proﬁles produced by the fully turbulent simulations and the
            free-transition simulations both with and without the compressibility corrections at the 2523 (≈ 1.5◦ ),

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 22 ===== -->

           22        Piotrowski and Zingg

                     (a)                                                              (b)

                     (c)                                                              (d)

                     (e)                                                              (f)

           Figure 16. Upper-surface skin friction coeﬃcient proﬁles for the CRM-NLF obtained on the L1 grid
           overlaid with the estimated natural transition fronts from the experiment [11].

           2525 (≈ 2.5◦ ), and 2526 (≈ 3.0◦ ) test conditions are extracted at nine spanwise locations on the L1 grid
           and illustrated in Figs 18, 19, and 20, respectively. Similar to the 2524 results presented in Fig. 14, there
           is a consistent trend with the compressibility corrections delaying crossﬂow transition on the lower
           surface, and signiﬁcantly delaying Tollmien-Schlichting transition on the upper surface in regions of
           favourable pressure gradient. There is again signiﬁcant bypass transition observed in the experiment at
           the 2523, 2525, and 2526 test conditions, especially at the leading-edge pressure ports [11, 12], which
           makes a detailed comparison of the pressure coeﬃcient proﬁles diﬃcult.

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 23 ===== -->

                                                                                        The Aeronautical Journal   23

                     Figure 17. Residual convergence histories for the CRM-NLF simulations on the L1 grid.

                      (a)                                      (b)                    (c)

                      (d)                                      (e)                    (f)

                      (g)                                      (h)                    (i)

            Figure 18. Pressure and skin friction coeﬃcient proﬁles for the CRM-NLF at the 2523 test conditions
            (α ≈ 1.5◦ ) obtained on the L1 grid compared with the pressure proﬁles from the experiment [11] at
            varying spanwise stations η.

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 24 ===== -->

           24        Piotrowski and Zingg

                    (a)                                       (b)                     (c)

                    (d)                                       (e)                     (f)

                    (g)                                       (h)                     (i)

           Figure 19. Pressure and skin friction coeﬃcient proﬁles for the CRM-NLF at the 2525 test conditions
           (α ≈ 2.5◦ ) obtained on the L1 grid compared with the pressure proﬁles from the experiment [11] at
           varying spanwise stations η.

           4.0 Conclusions
           Compressibility corrections have been developed and applied to extend the SA-sLM2015 smooth local
           correlation-based transition model to transonic ﬂow regimes. A compressibility correction for Tollmien-
           Schlichting instabilities has been developed to reproduce trends predicted by stability analysis and a
           crossﬂow source term function has been developed to apply a compressibility correction for stationary
           crossﬂow instabilities. These corrections and modiﬁcations do not impact the predictive capability of
           the model in the incompressible ﬂow regime and do not have a signiﬁcant impact on its iterative and
           grid convergence behaviour.
              The smooth transition model both with and without the compressibility corrections is applied to a
           range of transonic external aerodynamic transition test cases. The results demonstrate that at higher
           Reynolds numbers the compressibility corrections successfully delay the predicted transition onset
           locations to better agree with experimental data. Speciﬁcally, the Tollmien-Schlichting compressibility
           correction delays transition in regions of favourable pressure gradient where stability analysis predicts a
           strong stabilising eﬀect. The correction is less stabilising in adverse pressure gradients, enabling accurate
           prediction of transition locations at lower Reynolds number conditions. The crossﬂow compressibility

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 25 ===== -->

                                                                                                     The Aeronautical Journal               25

                     (a)                                       (b)                                (c)

                     (d)                                       (e)                                (f)

                     (g)                                       (h)                                (i)

            Figure 20. Pressure and skin friction coeﬃcient proﬁles for the CRM-NLF at the 2526 test conditions
            (α ≈ 3.0◦ ) obtained on the L1 grid compared with the pressure proﬁles from the experiment [11] at
            varying spanwise stations η.

            correction prevents an inaccurate, upstream transition front from forming on the upper surface of the
            CRM-NLF wing where crossﬂow instabilities are not expected to be dominant.
               The results demonstrate that the compressibility corrections successfully extend the empirical cor-
            relations in the smooth local correlation-based transition model to transonic ﬂow regimes. Speciﬁcally,
            the corrections and the underlying transition model are investigated over a range of Mach numbers from
            0.71 to 0.856, and for a range of Reynolds numbers from 2 × 106 to 15 × 106 . However, further work is
            required to validate RANS-based transition models over a wider range of ﬂow conditions, speciﬁcally
            at the higher Reynolds numbers typical of commercial transport aircraft — the increased availabil-
            ity of high-quality data from experiments at these ﬂight conditions is a necessary prerequisite for this
            development.
            Acknowledgments. This work was partially funded by the NASA Transformational Tools and Technologies (TTT) and Advanced
            Air Transport Technology (AATT) projects, the Natural Sciences and Engineering Research Council (NSERC), and the University
            of Toronto. Computations were performed on the Niagara supercomputer at the SciNet HPC Consortium, a part of Compute
            Canada. The authors gratefully acknowledge the input and support provided by Dr. Thomas Reist, Dr. Holger Mai, and Dr. Cetin
            Kiris. In addition, the authors sincerely appreciate the generous support provided by Dr. Marco Costantini in providing high-quality
            data for the VA-2 aerofoil experiment.

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 26 ===== -->

           26        Piotrowski and Zingg

            References
            [1] Malik, M.R., Crouch, J.D., Saric, W.S., Lin, J.C. and Whalen, E.A. Application of drag reduction techniques to transport
                aircraft, In R. Blockley and W. Shyy (eds) Encyclopedia of Aerospace Engineering, 2015.
            [2] Reed, H.L. and Saric, W.S. Transition mechanisms for transport aircraft, 38th Fluid Dynamics Conference and Exhibit,
                AIAA Paper 2008-3743, Seattle, Washington, 2008.
            [3] Arnal, D., Casalis, G. and Houdeville, R. Practical Transition Prediction Methods: Subsonic and Transonic Flows, VKI
                Lecture Series: Advances in Laminar-Turbulent Transition Modeling, 2009, pp. 1–34.
            [4] Aupoix, B., Arnal, D., Bezard, H., Chaouat, B. and Chedevergne, F. Transition and Turbulence Modeling, AerospaceLab,
                2011, pp 1–13.
            [5] Pasquale, D., Rona, A. and Garrett, S.J. A selective review of CFD transition models, 39th AIAA Fluid Dynamics Conference,
                AIAA Paper 2009-3812, San Antonio, Texas, 2009.
            [6] Langtry, R.B. and Menter, F.R. Correlation-based transition modeling for unstructured parallelized computational ﬂuid
                dynamics codes, AIAA J, 2009, 47, (12), pp. 2894–2906.
            [7] Langtry, R.B., Sengupta, K., Yeh, D.T. and Dorgan, A.J. Extending the γ R̃eθt Correlation based Transition Model for
                Crossﬂow Eﬀects, 45th AIAA Fluid Dynamics Conference, AIAA Paper 2015-2474, 2015.
            [8] Venkatachari, B.S., Paredes, P., Derlaga, J.M., Buning, P., Choudhari, M., Li, F. and Chang, C.L. Assessment of RANS-
                based transition models based on experimental data of the common research model with natural laminar ﬂow, 2021 AIAA
                Aerospace Sciences Meeting, AIAA Paper 2021-1430, 2021.
            [9] Fehrs, M. One-equation transition model for airfoil and wing aerodynamics, In Dillmann A. et al. (eds) New Results in
                Numerical and Experimental Fluid Mechanics XI. Notes on Numerical Fluid Mechanics and Multidisciplinary Design, Vol.
                136, Springer, Cham, 2018, pp. 199–208.
           [10] Campbell, R.L. and Lynde, M.N. Natural Laminar ﬂow design for wings with moderate sweep, 34th AIAA Applied
                Aerodynamics Conference, AIAA Paper 2016-4326, Washington D.C., 2016.
           [11] Lynde, M.N., Campbell, R.L., Rivers, M.B., Viken, S.A., Chan, D.T., Watkins, N.A. and Goodliﬀ, S.L. Preliminary results
                from an experimental assessment of a natural laminar ﬂow design method, 2019 AIAA Aerospace Sciences Meeting, AIAA
                Paper 2019-2298, 2019.
           [12] Lynde, M.N., Campbell, R.L. and Viken, S.A. Additional ﬁndings from the common research model natural laminar ﬂow
                wind tunnel test, AIAA Aviation 2019 Forum, AIAA Paper 2019-3292, 2019.
           [13] Arnal, D. Transition Prediction in Transonic Flows, Symposium Transsonicum III. International Union of Theoretical and
                Applied Mechanics, Springer, Berlin, Heidelberg, 1988.
           [14] Risuis, S., Costantini, M., Koch, S., Hein, S. and Klein, C. Unit Reynolds number, Mach number and pressure gradient
                eﬀects on laminar-turbulent transition in two-dimensional boundary layers, Exp Fluids, 59, (5), 2018.
           [15] Van Driest, E.R. Calculation of the stability of the laminar boundary layer in a compressible ﬂuid on a ﬂat plate with heat
                transfer, J Aeronaut Sci, 1952, 19, (12), pp. 801–812.
           [16] Lees, L. and Lin, C.C. Investigation of the stability of the laminar boundary layer in a compressible ﬂuid, National Advisory
                Committee for Aeronautics, NACA TN No. 1115, 1946.
           [17] Mack, L.M. Linear stability theory and the problem of supersonic boundary-layer transition, AIAA J, 1975, 13, (3), pp.
                278–289.
           [18] Arnal, D. Laminar-turbulent transition problems in supersonic and hypersonic ﬂows, Special Course on
                Aerothermodynamics of Hypersonic Vehicles, AGARD, AGARD-R-761, 1989.
           [19] Arnal, D., Habiballah, M. and Coutols, E. Laminar instability theory and transition criteria in two and three-dimensional
                ﬂow, La Recherche Aérospatiale (English Edition)(ISSN 0379-380X), 1984, 2, pp. 45–63.
           [20] Smith, A.M.O. and Gamberoni, N. Transition, Pressure Gradientand Stability Theory, Douglas Aircraft Report ES-26388,
                1956.
           [21] Van Ingen, J.L. A Suggested Semi-Empirical Method for the Calculation of the Boundary Layer Transition Region,
                Technische Hogeschool Delft, Vliegtuigbouwkunde, Rapport VTH-74, 1956.
           [22] Arnal, D., Houdeville, R., Seraudie, A. and Vermeersch, O. Overview of Laminar-Turbulent transition investigations at
                ONERA toulouse, 41st AIAA Fluid Dynamics Conference and Exhibit, 2011.
           [23] Perraud, J. and Durant, A. Stability-based Mach zero to four longitudinal transition prediction criterion, J Spacecraft Rockets,
                2016, 53, (4), pp. 730–742.
           [24] Pascal, L., Delattre, G., Deniau, H. and Cliquet, J., Stability-Based transition model using transport equations, AIAA J, 58,
                (7), 2020, pp. 2933–2942.
           [25] Ströer, P., Krimmelbein, N., Krumbein, A. and Grabe, C. Stability-based transition transport modeling for unstructured
                computational ﬂuid dynamics including convection eﬀects, AIAA J, 2020, 58, (4), pp. 1506–1517.
           [26] Ströer, P., Krimmelbein, N., Krumbein, A. and Grabe, C. Stability-based transition transport modeling for unstructured
                computational ﬂuid dynamics at transonic ﬂow conditions, AIAA J, 2021, 59, (9), pp. 3585–3597.
           [27] Menter, F.R. Two-Equation eddy-viscosity turbulence models for engineering applications, AIAA J, 1994, 32, (8), pp. 1598–
                1605.
           [28] François, D.G., Krumbein, A., Krimmelbein, N. and Grabe, C. Simpliﬁed stability-based transition transport modeling for
                unstructured computational ﬂuid dynamics, AIAA SciTech 2022 Forum, AIAA 2022-1543, 2022.
           [29] Ströer, P., Krimmelbein, N., Krumbein, A. and Grabe, C. Galilean-invariant stability-based transition transport modeling
                framework, AIAA J, Article in Advance, 2022.
           [30] Dagenhart, J.R. Ampliﬁed crossﬂow disturbances in the laminar boundary layer on swept wings with suction, NASA
                Technical Paper (1902), 1981.

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 27 ===== -->

                                                                                                     The Aeronautical Journal               27

            [31] Arnal, D. Boundary layer transition: Predictions based on linear theory, Special Course on Progress in Transition Modelling,
                 AGARD, AGARD-R-793, 1993.
            [32] Malik, M.R., Balakumar, P. and Chang, C. Linear stability of hypersonic boundary layers, 10th National Aero-Space Plane
                 Symposium, No. 189, 1991.
            [33] Kroo, I. and Sturdza, P. Design-oriented aerodynamic analysis for supersonic laminar ﬂow wings, 41st Aerospace Sciences
                 Meeting and Exhibit, 2003.
            [34] Lee, J. and Jameson, A. Natural-laminar-ﬂow airfoil and wing design by adjoint method and automatic transition prediction,
                 47th AIAA Aerospace Sciences Meeting and Exhibit, AIAA Paper 2009-897, Orlando, Florida, January 2009.
            [35] Xu, J., Qiao, L. and Bai, J. Improved local ampliﬁcation factor transport equation for stationary crossﬂow instability in
                 subsonic and transonic ﬂows, Chin J Aeronaut, 2020, 33, (12), pp. 3073–3081.
            [36] Piotrowski, M.G.H. and Zingg, D.W. Smooth local correlation-based transition model for the Spalart-Allmaras turbulence
                 model, AIAA J, 2021, 59, (2), pp. 474–492.
            [37] Spalart, P.R. and Allmaras, S.R. A one-equation turbulence model for aerodynamic ﬂows, 30th AIAA Aerospace Sciences
                 Meeting and Exhibit, AIAA Paper 092-0439, Reno, Nevada, United States 1992.
            [38] Piotrowski, M.G.H. and Zingg, D.W. Investigation of a smooth local correlation-based transition model in a discrete-adjoint
                 aerodynamic shape optimization algorithm, AIAA Scitech 2022 Forum, AIAA 2022-1865, 2022.
            [39] Perraud, J., Deniau, H. and Casalis, G., Overview of Transition Prediction Tools in the elsA Software, ECCOMAS 2014,
                 2014.
            [40] Piotrowski, M.G.H. Development of a transition prediction methodology suitable for aerodynamic shape optimization, PhD
                 thesis, Graduate Department of Aerospace Science and Engineering, University of Toronto, 2022.
            [41] Pulliam, T.H. and Zingg, D.W. Fundamental Algorithms in Computational Fluid Dynamics, Springer International
                 Publishing, 2014.
            [42] Abu-Ghannam, B.J. and Shaw, R. Natural transition of boundary layers, the eﬀects of turbulence, pressure gradient, and
                 ﬂow history, J Mech Eng Sci, 1980, 22, (5), pp. 213–228.
            [43] Coder, J.G. and Maughmer, M.D. Computational ﬂuid dynamics compatible transition modeling using an ampliﬁcation
                 factor transport equation, AIAA J, 2014, 52, (11), pp. 2506–2512.
            [44] Grabe, C., Shengyang, N. and Krumbein, A. Transition transport modeling for the prediction of crossﬂow transition, 34th
                 AIAA Applied Aerodynamics Conference, AIAA Paper 2016-3572, Washington D.C., June 2016.
            [45] Grabe, C., Shengyang, N. and Krumbein, A. Transition transport modeling for the prediction of crossﬂow transition, AIAA
                 J, 2018, 56, (12), pp. 3167–3178.
            [46] Carnes, J.A. and Coder, J.G. Eﬀect of crossﬂow transition on the pressure-sensitive-paint rotor in Hover, J Aircr, 2022, 59,
                 (1), pp. 29–46.
            [47] Dagenhart, J. and Saric, W. Crossﬂow stability and transition experiments in swept-wing ﬂow, NASA Langley Technical
                 Report Server, NASA/TP-1999-209344, 1999.
            [48] Radeztsky, R.H., Reibert, M.S. and Saric, W.S. Eﬀect of micron-sized roughness on transition in swept-wing ﬂows, 31st
                 Aerospace Sciences Meeting and Exhibit, AIAA Paper 93-0076, 1993.
            [49] Medida, S. Correlation-based transition modeling for external aerodynamic ﬂows, PhD thesis, University of Maryland,
                 College Park, 2014.
            [50] Petzold, R. and Radespiel, R. Transition on a Wing with Spanwise varying crossﬂow and linear stability analysis, AIAA J,
                 2015, 53, (2), pp. 321–335.
            [51] Swanson, R.C. and Turkel, E. On central-diﬀerence and upwind schemes, J Comput Phys, 1992, 101, pp. 292–306.
            [52] Piotrowski, M.G.H. and Zingg, D.W. Numerical Behaviour of a smooth local correlation-based transition model in a Newton-
                 Krylov Flow Solver, AIAA Scitech 2022 Forum, AIAA 2022-0909, 2022.
            [53] Stanewsky, E. and Zimmer, H. Development and wind tunnel investigation of three supercritical airfoil proﬁles for transport
                 aircraft, NASA, TM-75840, Washington, D.C., 1980.
            [54] Hebler, A., Schojda, L. and Mai, H. Experimental investigation of the aeroelastic behavior of a Laminar airfoil in transonic
                 ﬂow, Proc. IFASD 2013, IFASD, Bristol, 2013.
            [55] Fehrs, M., van Rooij, A.C. and Nitzsche, J. Inﬂuence of boundary layer transition on the ﬂutter behavior of a supercritical
                 airfoil, CEAS Aeronaut J, 2015, 6, (2), pp. 291–303.
            [56] Fehrs, M., Helm, S. and Kaiser, C. Numerical Investigation of Unsteady Transitional Boundary Layer Flows, IFASD 2019
                 -International Forum on Aeroelasticity and Structural Dynamics, 2019.
            [57] Drela, M. and Giles, M.B. Viscous-inviscid analysis of transonic and low Reynolds number airfoils, AIAA J, 1987, 25, (10),
                 pp. 1347–1355.
            [58] Mack, L.M. Transition and laminar instability, Jet Propulsion Laboratory Publication, NASA-CP-153203, 1977.
            [59] Dress, D.A., Johnson, C.B., McGuire, P.D., Stanewsky, E. and Ray, E.J. High Reynolds Number Tests of the CAST 10-
                 2/DOA 2 Airfoil in the Langley 0.3-Meter Transonic Cryogenic Tunnel, Phase I, NASA, TM-84620, Hampton, 1983.
            [60] Venkatachari, B.S., Paredes, P., Choudhari, M.M., Li, F. and Chang, C.-L., Transition analysis for the CRM-NLF wind
                 tunnel conﬁguration using transport equation models and linear stability correlations, AIAA SciTech 2022 Forum, AIAA
                 2022-1542, 2022.
            [61] Costantini, M., Lee, T., Nonomura, T., Asai, K. and Klein, C., Feasibility of skin-friction ﬁeld measurements in a transonic
                 wind tunnel using a global luminescent oil ﬁlm, Exp Fluids, 2021, 62, (1), p 21.
            [62] Costantini, M., Henne, U., Risius, S. and Klein, C., A robust method for reliable transition detection in temperature-sensitive
                 paint data, Aerosp Sci Technol, 2021, 113, p 106702.

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 28 ===== -->

           28        Piotrowski and Zingg

           [63] Crouch, J., Sutanto, M., Witkowski, D., Watkins, A., Rivers, M. and Campbell, R. Assessment of the National Transonic
                Facility for Natural Laminar Flow Testing, 48th AIAA Aerospace Sciences Meeting including the New Horizons Forum and
                Aerospace Exposition, 2010.
           [64] Allmaras, S.R. and Johnson, F.T. Modiﬁcations and clariﬁcations for the implementation of the Spalart-Allmaras turbulence
                model, Seventh International Conference on Computational Fluid Dynamics, ICCFD7 Paper 1902, 2012, pp. 1–11.
           [65] Spalart, P.R. Strategies for turbulence modelling and simulation, Int J Heat Fluid Flow, 2000, 21, (3), pp. 252–263.
           [66] Helm, S., Fehrs, M., Krimmelbein, N. and Krumbein, A. Transition prediction and analysis of the CRM-NLF wing with the
                DLR TAU code, In Dillmann A., Heller G., Krämer E. and Wagner C. (eds) New Results in Numerical and Experimental
                Fluid Mechanics XIII. STAB/DGLR Symposium 2020. Notes on Numerical Fluid Mechanics and Multidisciplinary Design,
                Vol. 151, Springer, Cham., 2021.
           [67] Paredes, P., Venkatachari, B.S., Choudhari, M., Li, F., Hildebrand, N. and Chang, C.L. Transition Analysis for the CRM-NLF
                Wind Tunnel Conﬁguration, 2021 AIAA Aerospace Sciences Meeting, AIAA Paper 2021-1431, 2021.
           [68] Krimmelbein, N. and Krumbein, A. Determination of Critical N-Factors for the CRM-NLF Wing, In Dillmann A., Heller G.,
                Krämer E. and Wagner C. (eds) New Results in Numerical and Experimental Fluid Mechanics XIII. STAB/DGLR Symposium
                2020. Notes on Numerical Fluid Mechanics and Multidisciplinary Design, Vol. 151, Springer, Cham, 2021.

           Appendix A: SA-sLM2015 transition model equations
           A.1 Transition onset momentum-thickness Reynolds number transport equation
                                                                                                 
                               ∂ R̃eθt      ∂ R̃eθt                 1 ∂                   ∂ R̃eθt
                                       + uj         = Pθt + Dscf +          σθt (ν + νt )                                          (A.1)
                                 ∂t          ∂xj                   Re∞ ∂xj                 ∂xj
                                                                         cθt
                                                             Pθt =           (Reθt − R̃eθt )(1 − Fθt )                             (A.2)
                                                                          t
                                                              cθt
                                                    Dscf =        ccrossflow φ300 (Reθt,scf − R̃eθt , 0)(Fθt )                     (A.3)
                                                               t
                                                                                                             ρd2 S
                                          Fθt = Fwake e−( δ ) ,
                                                              d   4                          −6
                                                                          Fwake = e−ReS ×10 ,        ReS =         Re∞             (A.4)
                                                                                                              μ
                                                     50d                         15                  R̃eθt μ 1
                                               δ=         δBL ,           δBL =      θBL ,    θBL =                                (A.5)
                                                       U                           2                    ρU Re∞

                                                       cθt = 0.03,          σθt = 2.0,       ccrossflow = 0.6                      (A.6)

                                                                                500μ 1
                                                                           t=                                                      (A.7)
                                                                                ρU 2 Re∞
                                               ⎧                                
                                               ⎪                        0.2196
                                     ρUθ       ⎨ 1173.51 − 589.428Tu∞ +             F(λθ ) Tu∞ ≤ 1.3
                             Reθ t =     Re∞ =                            Tu2∞                                                     (A.8)
                                      μ        ⎪
                                               ⎩
                                                    331.50[Tu∞ − 0.5658]−0.671 F(λθ )      Tu∞ > 1.3
                                                                                         Tu∞
                                                      F(λθ )1 = 1 + 0.275 1 − e[−35λθ ] e−[ 0.5 ]                                  (A.9)

                                                                      F(λθ )2 = φ300 (F(λθ )1 , 1)                               (A.10)

                                                                                         Tu∞ 1.5
                                       F(λθ )3 = 1 − −12.986λθ − 123.66λ2θ − 405.689λ3θ e−[ 1.5 ]                                (A.11)

                                                             F(λθ ) = φ−300 (F(λθ )2 , F(λθ )3 )                                 (A.12)

                                                            ρθ 2 dU                                             1
                                                     λθ =           Re∞ ,           U = (u2 + v2 + w2 ) 2                        (A.13)
                                                             μ ds

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 29 ===== -->

                                                                                                             The Aeronautical Journal      29

                                                                                                                
                                                        dU      u          dU  v            dU  w          dU
                                                           =                  +                 +                                       (A.14)
                                                        ds      U          dx   U            dy   U          dz
                                                                                      
                                                      dU             2 − 12   du dv dw
                                                         = (u + v + w ) · u + v + w
                                                             2   2
                                                                                                                                        (A.15)
                                                      dx                      dx dx dx
                                                                                      
                                                      dU             2 − 12   du dv dw
                                                         = (u + v + w ) · u + v + w
                                                             2   2
                                                                                                                                        (A.16)
                                                      dy                      dy dy dy
                                                                                            
                                                      dU                    1    du  dv   dw
                                                         = (u2 + v2 + w2 )− 2 · u + v + w                                               (A.17)
                                                      dz                         dz  dz   dz
                                 U 
                               ρ 0.82 θt                  h                                             +
                                                                                                                        −
                                                                                                                                    
                    Reθt,scf =           Re∞ = −35.088 ln                         + 319.51 + f          Hcrossflow   − f Hcrossflow     (A.18)
                                  μ                       θt

                                                                                                                      μt
                                                     Hcrossflow = Hcrossflow (1.0 + φ−300 (RT , 0.4)),       RT =                       (A.19)
                                                                                                                      μ

                                                              +
                                                             Hcrossflow = φ300 (0.1066 − Hcrossflow , 0)                                (A.20)

                                                       +
                                                                          +
                                                                                               +
                                                                                                           2
                                             f         Hcrossflow   = 6200 Hcrossflow   + 50000 Hcrossflow                              (A.21)

                                                           −
                                                          Hcrossflow = φ300 (− (0.1066 − Hcrossflow ) , 0)                              (A.22)

                                                                    −
                                                                                                  −
                                                                                                Hcrossflow
                                                            f       Hcrossflow   = 75tanh                                               (A.23)
                                                                                                0.0125

                                                                    u                    v           w
                                            U= √                             ,√                ,√                                       (A.24)
                                                           u2 + v2 + w2           u2 + v2 + w2   u2 + v2 + w2

                                                                        ∂w ∂v ∂u ∂w ∂v ∂u
                                                            =             − ,    −      −                                              (A.25)
                                                                        ∂y  ∂z ∂z   ∂x ∂x ∂y

                                                                          streamwise = |U · |                                         (A.26)

                                                                                         dstreamwise
                                                                          Hcrossflow =                                                  (A.27)
                                                                                             U

            A.2 Intermittency transport equation
                                                                                             
                                            ∂γ          ∂γ                  1 ∂       νt ∂γ
                                                 + uj        = Pγ − Eγ +           ν+                                                   (A.28)
                                             ∂t         ∂xj              Re∞ ∂xj      σf ∂xj
                                                                            √       
                                                                         U∞ M∞ Re∞     √
                                       Pγ = ca1 Flength Fonset φ−300 ,                  γ (1 − ce1 γ )                                 (A.29)
                                                                          l     20

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press

<!-- ===== PDF page 30 ===== -->

           30        Piotrowski and Zingg

                                                                      √      
                                                                    U∞ M∞ Re∞
                                            Eγ = ca2 Fturb φ−300 ,             γ (ce2 γ − 1)                                       (A.30)
                                                                     l   20
                                                                                      
                                                                  tanh 6 Fonset,1 − 1.35 + 1
                                                         Fonset =                                                                   (A.31)

                                                                                ReS
                                                           Fonset,1 =                          + (RT )2                             (A.32)
                                                                             2.60Reθc
                                                                                       
                                                                             R̃eθt
                                                   Reθc = 0.67R̃eθt + 24 sin       + 0.5 + 14                                       (A.33)

                                                                   44 − (0.50 − 3 · 10−4 (R̃eθt − 596))
                                                Flength = 44 −                                     1
                                                                                                                                    (A.34)
                                                                                  (1 + Flength,1 ) 6
                                                                                               
                                                        Flength,1 = exp −3 · 10−2 (R̃eθt − 460)                                     (A.35)

                                                              Fturb = (1 − Fonset ) exp(−RT )                                       (A.36)

                                         ce1 = 1.0,       ca1 = 2.0,       σf = 1.0,       ce2 = 50,          ca2 = 0.06            (A.37)

           A.3 Spalart-Allmaras turbulence model (SA-neg-noft2)
                                                                                                                               
                           ∂ ν̃      ∂ ν̃                 1    ∂           ∂ ν̃                                     ∂ ν̃ ∂ ν̃
                                + uj      = Pν̃ − Dν̃ +           (ν + ν̃)                                    + cb2                 (A.38)
                            ∂t       ∂xj                σ Re∞ ∂xj          ∂xj                                      ∂xi ∂xi

                                                                                          cw1 fw ν̃
                                                                    cb1
                                                         Pν̃ = γ        S̃ν̃,     Dν̃ =                                             (A.39)
                                                                   Re∞                    Re∞ d

           Cite this article: Piotrowski M.G.H. and Zingg D.W. Compressibility corrections to extend a smooth local correlation-based
           transition model to transonic ﬂows. The Aeronautical Journal, https://doi.org/10.1017/aer.2022.105

https://doi.org/10.1017/aer.2022.105 Published online by Cambridge University Press
