<!--
NOTA DE EXTRAÇÃO (provenance / limitações) — ler antes de usar
Fonte: Piotrowski, M. G. H., e Zingg, D. W., "Smooth Local Correlation-Based
Transition Model for the Spalart–Allmaras Turbulence Model", AIAA Journal,
DOI: 10.2514/1.J059784 (versão publicada / version of record).
O preprint PiotrowskiandZingg2020 é o MESMO artigo (conteúdo ~95% idêntico);
diferenças menores em agrupamento de citações e formatação.

O que é FIEL aqui:
- Prosa: reproduzida do texto do PDF (ordem de leitura e espaços corrigidos).
- Tabelas: reconstruídas dos valores do PDF.
- Referências: as 91, reproduzidas do PDF.
- Equações: transcritas para LaTeX e cruzadas com a extração de layout.
  -> Devem ser VERIFICADAS contra o PDF. Onde um símbolo é ambíguo na fonte,
     está marcado com ⚠ em vez de adivinhado.

O que NÃO pode ser reproduzido (limite real, não inventei nada):
- As FIGURAS (imagens). Só as legendas são reproduzidas, verbatim, mais uma
  descrição do que a figura contém segundo os seus próprios rótulos. Os dados
  traçados (curvas de convergência, etc.) NÃO existem em texto e não foram
  fabricados.
-->

# Smooth Local Correlation-Based Transition Model for the Spalart–Allmaras Turbulence Model

**Michael G. H. Piotrowski**\*
NASA Ames Research Center, Moffett Field, California 94035
and
**David W. Zingg**†
University of Toronto, Toronto, Ontario M3H 5T6, Canada

DOI: https://doi.org/10.2514/1.J059784 — *AIAA Journal* (Article in Advance)

> \* NASA Pathways Ph.D. Candidate, University of Toronto, NASA Computational Aeroscience Branch; currently Institute for Aerospace Studies, University of Toronto, ON M3H 5T6; m.piotrowski@mail.utoronto.ca. Student Member AIAA.
> † Distinguished Professor of Computational Aerodynamics and Sustainable Aviation, Institute for Aerospace Studies; dwz@oddjob.utias.utoronto.ca. Associate Fellow AIAA.
>
> Received 12 May 2020; revision received 12 August 2020; accepted for publication 31 August 2020; published online 22 October 2020. This material is declared a work of the U.S. Government and is not subject to copyright protection in the United States.

---

## Abstract

A smooth version of the $\gamma\text{-}\tilde{Re}_{\theta t}$ local correlation-based transition model (LM2015) for the Spalart–Allmaras (SA) turbulence model is presented. The LM2015 model with helicity-based crossflow correlations is modified and coupled to the SA turbulence model, designated SA-LM2015. The LM2015 and SA-LM2015 transition models include source terms that contain stiff and nonsmooth functions. Approximations to these functions are introduced to eliminate discontinuities and improve the numerical behavior of the model, with the smooth model designated SA-sLM2015. Deep convergence is achieved using a fully coupled, implicit Newton–Krylov algorithm globalized using an efficient pseudotransient continuation strategy. Modifications to the Newton–Krylov algorithm are introduced, including a source-term time stepping strategy, to address the large sources introduced by the turbulence and transition model equations. Two- and three-dimensional transition test cases demonstrate that both models, SA-LM2015 and SA-sLM2015, are able to predict transition due to a variety of mechanisms accurately, with the smooth variant displaying significantly improved numerical behavior. Several of the strategies presented, including the smoothing techniques and source-term time stepping, could be useful in the context of other transition models as well.

---

## Nomenclature

| Symbol | Definition |
|---|---|
| $d$ | minimum off-wall distance |
| $H_{crossflow}$ | crossflow strength parameter |
| $h$ | surface roughness height (rms) |
| $M$ | Mach number |
| $p$ | smooth maximum/minimum function parameter |
| $p_{switch}$ | proximity switch cutoff value, $10^{-15}$ |
| $R_T$ | eddy-viscosity ratio |
| $Re$ | Reynolds number |
| $Re_S$ | strain-rate magnitude Reynolds number |
| $Re_{\theta t}$ | transition onset momentum-thickness Reynolds number |
| $Re_\Omega$ | vorticity Reynolds number |
| $\tilde{Re}_{\theta t}$ | transported quantity of local transition onset momentum-thickness Reynolds number |
| $S$ | strain rate magnitude |
| $Tu$ | turbulence intensity |
| $U$ | local velocity magnitude |
| $\mathbf{U}$ | unit velocity vector |
| $u_i$ | Cartesian velocity component |
| $\gamma$ | intermittency |
| $\gamma_{source}$ | largest positive source-term Jacobian eigenvalue |
| $\delta$ | boundary-layer thickness |
| $\theta_t$ | transition onset momentum thickness |
| $\mu$ | molecular viscosity |
| $\mu_t$ | eddy viscosity |
| $\lambda_\theta$ | pressure gradient parameter |
| $\rho$ | density |
| $\phi_p$ | smooth maximum/minimum function |
| $\phi_{switch}$ | proximity switch for smooth maximum/minimum function |
| $\Omega$ | vorticity magnitude |
| $\Omega_{streamwise}$ | streamwise vorticity, helicity |
| $\boldsymbol{\Omega}$ | vorticity vector |

---

## I. Introduction

Laminar-flow wing designs produce an increased laminar extent of the boundary layer, resulting in a decrease in viscous drag [1]. These drag savings have the potential to significantly reduce the fuel burn of transport aircraft, where viscous drag constitutes approximately 50% of the total drag [2]. The potential of laminar-flow designs has been demonstrated by various studies which suggest the application of laminar flow control to large commercial aircraft can reduce aerodynamic drag by approximately 10% [3].

In the near term, natural laminar flow (NLF) is being investigated for winglet, tail, and nacelle designs, with the Boeing 787-8 and 777X representing the first commercial applications of NLF nacelles to large transport aircraft and the 737 Max aircraft exploiting laminar flow for both the nacelle and winglet designs [4]. In the future, the design of wings which exploit significant regions of laminar flow can result in more significant drag reduction. The Airbus BLADE project is investigating this, while experimental studies of the NASA common research model NLF demonstrate that it is possible to design transport wings with high sweep and Reynolds numbers with significant regions of laminar flow [5–7]. However, there are several operational challenges that need to be addressed. For example, the Boeing 757 EcoDemonstrator is investigating the effect of anti-insect devices and surface coatings for the application of laminar flow on wing leading edge surfaces.

Efficient high-fidelity numerical optimization algorithms provide the designers of next generation aircraft with powerful tools for the development of more fuel-efficient designs [8–10]. The integration of a flow solver with transition prediction capability into an aerodynamic shape optimization algorithm will allow for the study of various tradeoffs in the design of novel NLF wing configurations. The design of aircraft configurations using design tools capable of incorporating and exploiting NLF, which requires the ability to predict laminar–turbulent transition efficiently, could play a key role in reducing the environmental impact of aviation.

There are several mechanisms for boundary-layer transition. The method of transition depends on the application being studied, such as flow past wings, airframes, rotor blades, and wind turbine blades, as well as the characteristics of the flow, such as Reynolds number, angle of attack, pressure gradient, and turbulence intensity. For a swept transport-aircraft wing, the dominant mechanisms for transition result from Tollmien–Schlichting wave growth leading to natural transition, crossflow instabilities resulting from highly swept wings, concave curvature producing Görtler instabilities, and attachment-line instabilities as a result of large leading-edge radius and sweep at the root of the wing [11]. While the last two mechanisms can be prevented by appropriate profile design, balancing natural and crossflow-induced transition can be difficult as favorable pressure gradients used to stabilize streamwise instabilities destabilize crossflow instabilities [11]. The effects of laminar–turbulent transition can be included in numerical simulations using various transition prediction and modeling techniques, consisting of either local or nonlocal operations [12–14].

Flow solvers based on the Reynolds-averaged Navier–Stokes (RANS) equations provide an excellent balance between accuracy, robustness, and efficiency, making them a suitable choice for practical engineering applications that require the solution of large and complex geometries. However, the turbulence models typically used in RANS solvers do not have the stand-alone capability to predict boundary-layer transition. To predict transition, one must apply a transition criterion. Reviews by Arnal et al. [12], Aupoix et al. [13], and Pasquale et al. [14] describe the advantages and disadvantages of several approaches for coupling transition prediction criteria into RANS solvers. The most widely implemented transition prediction strategies for practical engineering applications are based on either the $e^N$ criterion or on transition onset functions.

The $e^N$ method for predicting transition is a semi-empirical approach, with the $N$-factor depending on freestream flow conditions, based on linear stability theory. To apply the $e^N$ criterion, one must first approximate the $N$-factor curves which represent the amplification ratios of the unstable frequencies of the disturbances in the boundary layer. This can be achieved using stability analysis methods such as the linearized or parabolized stability equations (PSEs) [15,16], or through the use of simplified $e^N$ envelope methods, developed by Drela and Giles [17], and used by Rashad and Zingg [18] and Mayda [19], depending on the level of fidelity required in a simulation.

A boundary-layer code can be used to increase the accuracy in determining the boundary-layer edge, which is required to satisfy the requirements of the stability analysis [20,21]. However, with sufficient grid resolution, it is possible to define the boundary-layer edge and extract the boundary-layer properties based on the RANS solution accurately, without using a boundary-layer solver [18,21]. The influence of the grid resolution on the stability analysis was studied in detail by Krimmelbein and Radespiel [21].

Transition models based on the $e^N$ method are able to accurately predict natural and separation-induced transition. In addition, correlations for crossflow instabilities (such as the C1 criterion [22]) have been successfully combined with these approaches, with experimental validation demonstrating accurate transition prediction on transonic swept wings in three dimensions [23–27]. However, $e^N$ methods require nonlocal boundary-layer quantities that are not directly accessible in RANS-based computational fluid dynamics (CFD) codes [28].

In the context of the current Paper, the requirement that the transition model be formulated locally is significant. Modern CFD codes use domain decomposition and parallel computation to increase the efficiency of simulations. Nonlocal operations, such as integrating quantities over the boundary layer, introduce difficulties to this process. Of the many transition models developed for RANS-based turbulence models, two which satisfy the requirement of a fully local formulation are the empirical correlation-based $\gamma\text{-}\tilde{Re}_{\theta t}$ local correlation-based transition model developed by Langtry and Menter (LM2009) [29] and the Amplification Factor Transport (AFT) transition model developed by Coder and Maughmer [30,31] and later revised by Coder [32], Stefanski et al. [33], Coder et al. [34], and Coder [35].

The AFT transition model is based on the approximate envelope method for Tollmien–Schlichting wave growth, developed by Drela and Giles [17], and is able to predict natural and separation-induced transition in subsonic and transonic flow regimes. The model consists of a single transport equation for the envelope amplification factor $\tilde{n}$ with a locally derived approximation for the boundary-layer shape factor. The transition model was originally coupled directly to the one-equation Spalart–Allmaras (SA) eddy-viscosity model [36], resulting in a two-equation coupled system [30–33]. The most recent model [34] is coupled to an additional, intermittency equation, inspired by the work of Langtry and Menter [29], resulting in a three-equation coupled system. Recent work by Xu et al. [37] extended the AFT model to include crossflow effects, using the Falkner–Skan and Cooke equations to develop correlations for local approximations of nonlocal boundary-layer quantities, with promising results. This work involved coupling the model to the Shear Stress Transport (SST) turbulence model, resulting in a four-equation system.

The transition model formulated by Langtry and Menter [29] is composed of two transport equations. The first is an equation for the intermittency $\gamma$, which is used to trigger the turbulence transition process by controlling the level of turbulence in the boundary layer. The second is an equation for the transition momentum-thickness Reynolds number $\tilde{Re}_{\theta t}$, which in the freestream is equal to the value produced by the empirical correlations $Re_{\theta t}$ and is diffused into the boundary layer using a standard diffusion term. Of the various crossflow transition models developed for the $\gamma\text{-}\tilde{Re}_{\theta t}$ transition model, two maintain a fully local formulation, methods based on the local C1 criterion [38–41] and local helicity approaches [42,43]. The local C1 approach uses the solutions of the Falkner–Skan and Cooke equations and is therefore only expected to work for winglike geometries with high aspect ratios [40], unlike methods based on the local helicity.

Crossflow instability correlations based on the local helicity were initially developed by Muller and Herbst [42], with more recent correlations, which include a framework for introducing roughness effects, developed by Langtry et al. (LM2015) [43]. Grabe et al. [40] conducted simulations using correlations based on both the C1 and helicity approaches. The C1 approach yielded more accurate results in terms of the comparison of predicted and measured transition locations on infinite swept wing flows, while the local helicity approach produced more accurate transition locations for nonwinglike geometries [40].

The local helicity-based correlations developed by Muller and Herbst [42] and Langtry et al. [43] both consist of an additional source term in the $\tilde{Re}_{\theta t}$ equation that acts as a sink to reduce the value of $\tilde{Re}_{\theta t}$ in regions of significant crossflow, triggering transition farther upstream. While the crossflow correlations of Muller and Herbst and Langtry et al. are similar, the latter represents a more thoroughly validated model [44]. Furthermore, the ability of the Langtry model to include roughness effects expands the applicability of the model.

The original $\gamma\text{-}\tilde{Re}_{\theta t}$ LM2009 model developed by Langtry and Menter [29] was coupled with the two-equation SST turbulence model [45]. Subsequent extensions of this model to include crossflow instabilities with correlations based on the local helicity by Grabe et al. [40], Muller and Herbst [42], and Langtry et al. [43] also included this coupling, while Nie et al. [46] investigated the performance of the Langtry–Menter model with C1 and helicity-based crossflow correlations coupled with a Reynolds stress model.

Medida and Baeder [47] modified the $\gamma\text{-}\tilde{Re}_{\theta t}$ model to be compatible with the SA turbulence model and developed a new set of correlations validated using several test cases, designating this model the $\gamma\text{-}\tilde{Re}_{\theta t}$-SA transition model. Recent work demonstrated the $\gamma\text{-}\tilde{Re}_{\theta t}$-SA model coupled with the crossflow correlations developed by Muller and Herbst [42] was able to predict transition for several three-dimensional test cases, including a swept wing, requiring only a slight modification [48]. However, their current formulation contains a function $G_{onset}$, which requires nonlocal flow information. More recently, Carnes and Coder [49] successfully coupled the algebraic helicity-based crossflow correlations developed by Langtry et al. [43] to the AFT model, with the results comparing well to experiment.

The goal of this Paper is to develop a RANS-based Newton–Krylov flow solver with laminar–turbulent boundary-layer transition prediction capabilities, to be used in a high-fidelity discrete-adjoint gradient-based aerodynamic shape optimization algorithm. To achieve this goal, the $\gamma\text{-}\tilde{Re}_{\theta t}$ LM2015 transition model with helicity-based crossflow correlations [43] is modified and coupled to the one-equation SA turbulence model, while maintaining a fully local formulation, reducing the computational cost of the original four-equation system (Sec. II). This model will be referred to as the SA-LM2015 transition model. Smooth approximations of nondifferentiable functions are developed to eliminate discontinuities in the model. Additionally, modifications to numerically stiff functions are produced in order to improve numerical convergence. This new, smooth model is designated the SA-sLM2015 transition model (Sec. III). Numerical methods for addressing the strong source terms present in the SA-LM2015 and SA-sLM2015 models, including a source-term time stepping strategy, are developed, and their implementation in the Newton–Krylov framework is discussed (Sec. IV). A range of two- and three-dimensional transition test cases is simulated, with the results demonstrating that the SA-LM2015 and SA-sLM2015 transition models are able to accurately predict transition due to a variety of transition mechanisms, with the SA-sLM2015 model demonstrating significantly improved convergence characteristics (Sec. V).

---

## II. SA-LM2015 Transition Model

The $\gamma\text{-}\tilde{Re}_{\theta t}$ LM2015 model with helicity-based crossflow correlations [43] is modified and coupled to the SA turbulence model [36]. The LM2015 model is decoupled from the SST turbulence model using modifications to the $R_T$ and $F_{wake}$ functions. Schücker [50] suggests replacing the specific turbulence dissipation rate $\omega$ in the original $F_{wake}$ function with velocity gradients. This approach is adopted in the current Paper using the strain-rate magnitude. A simplification of the $F_{\theta t}$ function is introduced, which reduces the amount of nonsmooth minimum/maximum operators present in the model. An adjustment to the $F_{onset1}$ term was applied to recalibrate the model, and the strain-rate magnitude in the original intermittency production term is replaced by vorticity magnitude to make the model more stable near laminar separation bubbles. The separation-induced transition modification developed by Langtry [51] was neglected, as the delay in the growth of turbulence after laminar separation bubbles produced by the SST turbulence model is not reproduced by the SA turbulence model, which is demonstrated by the results in Sec. V.

The remaining terms and coefficients remain unchanged from the original model. The modified LM2015 model with helicity-based crossflow correlations coupled to the SA turbulence model, designated SA-LM2015, is presented as follows.

### A. Momentum-Thickness Reynolds Number Transport Equation

The transport equation for the momentum-thickness Reynolds number, nondimensionalized by freestream quantities and written in nonconservative form to be consistent with the SA turbulence model, is defined as

$$
\frac{\partial \tilde{Re}_{\theta t}}{\partial t} + u_j \frac{\partial \tilde{Re}_{\theta t}}{\partial x_j}
= P_{\theta t} + D_{scf}
+ \frac{1}{Re}\frac{\partial}{\partial x_j}\!\left[\sigma_{\theta t}(\nu + \nu_t)\frac{\partial \tilde{Re}_{\theta t}}{\partial x_j}\right]
\tag{1}
$$

$$
P_{\theta t} = \frac{c_{\theta t}}{t}\,(Re_{\theta t} - \tilde{Re}_{\theta t})(1 - F_{\theta t})
\tag{2}
$$

$$
F_{\theta t} = F_{wake}\, e^{-(d/\delta)^4}
\tag{3}
$$

$$
\delta = \frac{50\,d\,\Omega}{U}\,\delta_{BL}; \qquad
\delta_{BL} = \frac{15}{2}\theta_{BL}; \qquad
\theta_{BL} = \frac{\tilde{Re}_{\theta t}\,\mu}{\rho U}\frac{1}{Re}
\tag{4}
$$

$$
F_{wake} = e^{-(Re_S/1\mathrm{E}{+}6)}; \qquad
Re_S = \frac{\rho d^2 S}{\mu}Re
\tag{5}
$$

$$
c_{\theta t} = 0.03; \qquad \sigma_{\theta t} = 2.0
\tag{6}
$$

$$
t = \frac{500\,\mu}{\rho U^2}\frac{1}{Re}
\tag{7}
$$

where $Re = (\rho_\infty a_\infty l)/\mu_\infty$, $a$ is the sound speed, $l$ is the characteristic (reference) length, and $\infty$ denotes a freestream quantity. The nondimensionalization procedure for the transition model equations is consistent with the procedure for the Navier–Stokes equations, which is described in [52] and is used, for example, in the well-known OVERFLOW flow solver.‡ For the cases presented in the current Paper, the reference length is specified as the root chord. The natural and bypass transition empirical correlations are defined as

$$
Re_{\theta t} = \frac{\rho U \theta}{\mu}Re =
\begin{cases}
\left[1173.51 - 589.428\,Tu + \dfrac{0.2196}{Tu^2}\right]F(\lambda_\theta) & Tu \le 1.3 \\[2ex]
331.50\,[Tu - 0.5658]^{-0.671}\,F(\lambda_\theta) & Tu > 1.3
\end{cases}
\tag{8}
$$

$$
F(\lambda_\theta) =
\begin{cases}
1 - \left[-12.986\lambda_\theta - 123.66\lambda_\theta^2 - 405.689\lambda_\theta^3\right]e^{-[Tu/1.5]^{1.5}} & \lambda_\theta \le 0 \\[1.5ex]
1 + 0.275\left[1 - e^{-35\lambda_\theta}\right]e^{-[Tu/0.5]} & \lambda_\theta > 0
\end{cases}
\tag{9}
$$

where the local pressure gradient parameter is given by

$$
\lambda_\theta = \frac{\rho \theta^2}{\mu}\frac{dU}{ds}Re; \qquad U = (u^2 + v^2 + w^2)^{1/2}
\tag{10}
$$

$$
\frac{dU}{ds} = \left(\frac{u}{U}\right)\frac{dU}{dx} + \left(\frac{v}{U}\right)\frac{dU}{dy} + \left(\frac{w}{U}\right)\frac{dU}{dz}
\tag{11}
$$

$$
\frac{dU}{dx} = (u^2+v^2+w^2)^{-1/2}\cdot\left[u\frac{du}{dx} + v\frac{dv}{dx} + w\frac{dw}{dx}\right]
\tag{12}
$$

$$
\frac{dU}{dy} = (u^2+v^2+w^2)^{-1/2}\cdot\left[u\frac{du}{dy} + v\frac{dv}{dy} + w\frac{dw}{dy}\right]
\tag{13}
$$

$$
\frac{dU}{dz} = (u^2+v^2+w^2)^{-1/2}\cdot\left[u\frac{du}{dz} + v\frac{dv}{dz} + w\frac{dw}{dz}\right]
\tag{14}
$$

The local turbulence intensity is not directly available when using the SA turbulence model. However, Suluksana et al. [53] demonstrated it can be redundant to use the local value of the turbulence intensity $Tu$ as well as the pressure gradient parameter $\lambda_\theta$, because the effects of pressure gradient are reduced by the decay of the local turbulence intensity. A favorable pressure gradient accelerates the flow leading to a decrease in the turbulence intensity, and therefore momentum-thickness Reynolds number, while an adverse pressure gradient leads to a reduction in turbulence intensity. The freestream turbulence intensity value $Tu_\infty$ is used throughout the domain in the current Paper, which is consistent with the approach taken by Medida and Baeder [47].

The helicity-based stationary crossflow correlations developed by Langtry et al. [43] are defined as follows:

$$
D_{scf} = \frac{c_{\theta t}}{t}\,c_{crossflow}\,\min\!\left(Re_{scf} - \tilde{Re}_{\theta t},\, 0\right)(F_{\theta t})
\tag{15}
$$

$$
c_{crossflow} = 0.6
\tag{16}
$$

$$
Re_{scf} = \frac{\rho (U/0.82)\theta_t}{\mu}Re = -35.088\,\ln\!\left(\frac{h}{\theta_t}\right) + 319.51 + f(\Delta H^{+}_{crossflow}) - f(\Delta H^{-}_{crossflow})
\tag{17}
$$

$$
\Delta H_{crossflow} = H_{crossflow}\left(1.0 + \min[R_T, 0.4]\right); \qquad R_T = \frac{\mu_t}{\mu}
\tag{18}
$$

$$
\Delta H^{+}_{crossflow} = \max\!\left(0.1066 - \Delta H_{crossflow},\, 0\right)
\tag{19}
$$

$$
f(\Delta H^{+}_{crossflow}) = 6200(\Delta H^{+}_{crossflow}) + 50000(\Delta H^{+}_{crossflow})^2
\tag{20}
$$

$$
\Delta H^{-}_{crossflow} = \max\!\left(-(0.1066 - \Delta H_{crossflow}),\, 0\right)
\tag{21}
$$

$$
f(\Delta H^{-}_{crossflow}) = 75\,\tanh\!\left(\frac{\Delta H^{-}_{crossflow}}{0.0125}\right)
\tag{22}
$$

Helicity, alternatively known as the streamwise vorticity $\Omega_{streamwise}$, is the magnitude of the scalar product of the local values of the unit velocity vector $\mathbf{U}$ and vorticity vector $\boldsymbol{\Omega}$,

$$
\mathbf{U} = \left(\frac{u}{\sqrt{u^2+v^2+w^2}},\ \frac{v}{\sqrt{u^2+v^2+w^2}},\ \frac{w}{\sqrt{u^2+v^2+w^2}}\right)
\tag{23}
$$

$$
\boldsymbol{\Omega} = \left(\frac{\partial w}{\partial y} - \frac{\partial v}{\partial z},\ \frac{\partial u}{\partial z} - \frac{\partial w}{\partial x},\ \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}\right)
\tag{24}
$$

and can be represented using the nondimensional crossflow strength $H_{crossflow}$ as

$$
\Omega_{streamwise} = |\mathbf{U}\cdot\boldsymbol{\Omega}|
\tag{25}
$$

$$
H_{crossflow} = \frac{d\,\Omega_{streamwise}}{U}
\tag{26}
$$

### B. Intermittency Transport Equation

The transport equation for intermittency, also nondimensionalized and written in nonconservative form, is

$$
\frac{\partial \gamma}{\partial t} + u_j \frac{\partial \gamma}{\partial x_j}
= P_\gamma - E_\gamma + \frac{1}{Re}\frac{\partial}{\partial x_j}\!\left[\left(\nu + \frac{\nu_t}{\sigma_f}\right)\frac{\partial \gamma}{\partial x_j}\right]
\tag{27}
$$

$$
P_\gamma = c_{a1}\,F_{length}\,F_{onset}\,\Omega\,\sqrt{\gamma}\,(1 - c_{e1}\gamma)
\tag{28}
$$

$$
F_{onset} = \max\!\left(\sqrt{F_{onset2} - F_{onset3}},\, 0\right)
\tag{29}
$$

> ⚠ **Verificar Eq. (29):** na fonte, o radical (√) e o seu alcance estão parcialmente ambíguos; a leitura mais provável do layout é $\sqrt{F_{onset2} - F_{onset3}}$. Confirmar contra o PDF se o radical cobre a diferença toda ou só $F_{onset2}$.

$$
F_{onset3} = \max\!\left(1 - \left(\frac{R_T}{2.5}\right)^3,\, 0\right)
\tag{30}
$$

$$
F_{onset2} = \min\!\left(\max\!\left(F_{onset1},\, F_{onset1}^4\right),\, 2\right)
\tag{31}
$$

$$
F_{onset1} = \frac{Re_S}{2.6\,Re_{\theta c}}
\tag{32}
$$

$$
c_{e1} = 1.0; \qquad c_{a1} = 2.0; \qquad \sigma_f = 1.0
\tag{33}
$$

$$
E_\gamma = c_{a2}\,F_{turb}\,\Omega\,\gamma\,(c_{e2}\gamma - 1)
\tag{34}
$$

$$
F_{turb} = e^{-(R_T/4)^4}
\tag{35}
$$

$$
c_{e2} = 50; \qquad c_{a2} = 0.06
\tag{36}
$$

The $F_{length}$ and $Re_{\theta c}$ correlations are as follows:

$$
F_{length} =
\begin{cases}
398.189\cdot10^{-1} + (-119.270\cdot10^{-4})\tilde{Re}_{\theta t} + (-132.567\cdot10^{-6})\tilde{Re}_{\theta t}^2 & \tilde{Re}_{\theta t} < 400 \\[1ex]
263.404 + (-123.939\cdot10^{-2})\tilde{Re}_{\theta t} + (194.548\cdot10^{-5})\tilde{Re}_{\theta t}^2 + (-101.695\cdot10^{-8})\tilde{Re}_{\theta t}^3 & 400 \le \tilde{Re}_{\theta t} < 596 \\[1ex]
0.5 - (\tilde{Re}_{\theta t} - 596.0)\cdot3\cdot10^{-4} & 596 \le \tilde{Re}_{\theta t} < 1200 \\[1ex]
0.3188 & 1200 \le \tilde{Re}_{\theta t}
\end{cases}
\tag{37}
$$

$$
Re_{\theta c} =
\begin{cases}
\tilde{Re}_{\theta t} - \big(396.035\cdot10^{-2} + (-120.656\cdot10^{-4})\tilde{Re}_{\theta t} + (868.230\cdot10^{-6})\tilde{Re}_{\theta t}^2 \\
\qquad\quad + (-696.506\cdot10^{-9})\tilde{Re}_{\theta t}^3 + (174.105\cdot10^{-12})\tilde{Re}_{\theta t}^4\big) & \tilde{Re}_{\theta t} \le 1870 \\[1ex]
\tilde{Re}_{\theta t} - \big(593.11 + 0.482(\tilde{Re}_{\theta t} - 1870.0)\big) & \tilde{Re}_{\theta t} > 1870
\end{cases}
\tag{38}
$$

### C. Coupling to SA Turbulence Model

The SA-neg-noft2 variant of the SA turbulence model [36,54] was used in the current Paper and is presented in nonconservative, nondimensional form:

$$
\frac{\partial \tilde{\nu}}{\partial t} + u_j \frac{\partial \tilde{\nu}}{\partial x_j}
= P_{\tilde{\nu}} - D_{\tilde{\nu}}
+ \frac{1}{\sigma Re}\left[\frac{\partial}{\partial x_j}\!\left((\nu + \tilde{\nu})\frac{\partial \tilde{\nu}}{\partial x_j}\right) + c_{b2}\frac{\partial \tilde{\nu}}{\partial x_i}\frac{\partial \tilde{\nu}}{\partial x_i}\right]
\tag{39}
$$

$$
P_{\tilde{\nu}} = \frac{c_{b1}}{Re}\tilde{S}\tilde{\nu}; \qquad
D_{\tilde{\nu}} = \frac{c_{w1}f_w}{Re}\left(\frac{\tilde{\nu}}{d}\right)^2
\tag{40}
$$

Osusky [55] presents the nondimensionalization procedure for the SA turbulence model, which is consistent with that used for the Navier–Stokes equations. The transition model is coupled to the production term $P_{\tilde{\nu}}$ through the following:

$$
\tilde{P}_{\tilde{\nu}} = \gamma P_{\tilde{\nu}}
\tag{41}
$$

---

## III. SA-sLM2015 Transition Model

The LM2015 and SA-LM2015 models contain several nonsmooth functions in the form of minimum/maximum (min/max) operators and conditional statements. These nonsmooth functions, referred to as simple kinks [56,57], produce discontinuities in the Jacobian. This presents a problem for both the Newton–Krylov solution algorithm and the discrete-adjoint gradient-based optimization algorithm. In the current Paper, this has been addressed through the application of smooth approximations to min/max operators and conditional statements. In addition, modifications are introduced, such as new $F_{onset}$ and $F_{turb}$ functions, in order to reduce numerical stiffness.

The modifications below describe the simplification and smoothing of the SA-LM2015 transition model, presented in Sec. II. For the remainder of this Paper, it will be designated SA-sLM2015.

### A. Minimum/Maximum Operators

Several examples of smooth approximations to minimum and maximum operators can be found in the literature. Rivest [58] suggested using a generalized mean-value operator, such as the $p$-mean:

$$
\phi_p(x) \triangleq \frac{\left(\sum_i^n (g_i(x))^p\right)^{1/p}}{n}
\tag{42}
$$

As $p$ approaches positive and negative infinity, we recover the maximum and minimum operators, respectively. However, this approximation only holds if $g = (g_1,\ldots,g_n)$ is a vector of positive real numbers, which limits its applicability. Another approach is to use an exponential penalty function [59–61], defined as

$$
\phi_p(x) \triangleq \frac{\log\left(\sum_i^n \exp(p\,g_i(x))\right)}{p}
\tag{43}
$$

Again, as $p$ approaches positive infinity, the function produces the maximum of $g_i$, while $p$ approaching negative infinity produces the minimum.

The advantage of this approach is that $g = (g_1,\ldots,g_n)$ can be a vector of any real numbers. A drawback is the presence of log and exponential functions, which can increase the computational cost. This is especially true as $g_i$ approaches values close to zero. This can introduce denormalized floating point numbers which significantly slow the algorithm. To prevent this, a proximity switch has been developed [Eq. (43)] to activate the exponential penalty function only in the region near the discontinuity. In addition, when considering a two-variable array, the exponential penalty function can be rewritten as

$$
\phi_p(x) \triangleq g_1(x) + \frac{\log(1 + \exp(p(g_2(x) - g_1(x))))}{p}
\tag{44}
$$

In combination with the proximity switch, this form of the exponential penalty function prevents both underflow and overflow. The proximity switch is given by

$$
\phi_{switch} \triangleq -\frac{\log(|p|\cdot p_{switch})}{|p|}
\tag{45}
$$

where $p_{switch}$ represents a cutoff value for the proximity check, which for the current Paper is set to a value of $10^{-15}$. If the magnitude of the difference between the two variables, $|g_1(x) - g_2(x)|$, is greater than $\phi_{switch}$, demonstrating we are far from the discontinuity, a simple min/max operator is used. As we approach the kink, if the difference between the two variables is less than $\phi_{switch}$, the exponential penalty function is used.

Algorithm 1 presents the complete method for approximating a minimum/maximum operator, based on the value of the smoothing parameter $p$ for two variables $g_1(x)$ and $g_2(x)$ using the exponential penalty functions including the proximity switch. Figure 1 illustrates the effect of this smoothing by approximating the maximum of a variable and a constant using various smoothing values $p$. For large values of $p$, the region of smoothing is small, and the sharp region is closely modeled. As $p$ is reduced, the domain where the exponential penalty functions are active increases, increasing computational cost and producing a smoother approximation of the kink region. Although adaptively updating the value of $p$ has been shown to help prevent ill-conditioned systems [60], all min/max operators present in the SA-sLM2015 model were approximated using this algorithm with fixed magnitudes of $p$ of 300, $\phi_{\pm300}$, which was found through numerical experimentation to provide a good balance between convergence properties and solution accuracy.

**Algorithm 1: Smooth Maximum/Minimum Function of Two Variables**

```
1:  function φ_p(g1(x), g2(x))
2:      p_switch ← 10^-15                                  ▷ cutoff value for proximity check
3:      a ← max(g1(x), g2(x))                              ▷ save larger value
4:      b ← min(g1(x), g2(x))                              ▷ save smaller value
5:      if p > 0 then                                      ▷ determine maximum
6:          if |a - b| > -log(|p|·p_switch)/|p| then       ▷ proximity switch
7:              φ_p(x) ← a                                 ▷ return maximum
8:              return
9:          else
10:             φ_p(x) ← { a + log(1 + exp(p(b - a)))/p }  ▷ exponential penalty function
11:     else                                               ▷ determine minimum
12:         if |a - b| > -log(|p|·p_switch)/|p| then       ▷ proximity switch
13:             φ_p(x) ← b                                 ▷ return minimum
14:             return
15:         else
16:             φ_p(x) ← { b + log(1 + exp(p(a - b)))/p }  ▷ exponential penalty function
```

> **Fig. 1** Approximation of $\max(0, x-1)$ using $\phi_p$ with various values for $p$.
> *[Figure not reproducible — line plot of $\max(0,x-1)$ vs $x$ compared against $\phi_{60}$, $\phi_{100}$, $\phi_{300}$, with an inset zoom near $x=1$. Plotted data not available in text.]*

### B. Intermittency Source Terms

#### 1. $F_{onset}$

The $F_{onset}$ function [Eqs. (29–32)] couples the SA, intermittency, and momentum-thickness Reynolds number equations. It was developed using several layers of minimum and maximum operators to activate the production of intermittency in the boundary layer, triggering the transition process. Solving the mean flow, turbulence, and transition equations in a fully coupled manner results in the $F_{onset}$ function being particularly problematic, and often the bottleneck in flow solver convergence. Rashad [62] investigated the design space of several transition region models and noted that the smoothness of the ramp up of the eddy viscosity was found to be particularly important to the smoothness of the design space. He also demonstrated the importance of achieving deep numerical convergence of the transition model residual. Therefore, special care was taken in ensuring that the $F_{onset}$ function be as smooth as possible, while maintaining the same characteristics of the original function.

The transition process is triggered when $F_{onset1}$ [Eq. (32)] exceeds unity and is sustained as the eddy-viscosity ratio increases. To smooth and simplify this process, the $F_{onset}$ function was reformulated using a hyperbolic tangent function as follows:

$$
F_{onset1} = \sqrt{\left(\frac{Re_S}{2.6\,Re_{\theta c}}\right)^2 + (R_T)^2}
\tag{46}
$$

$$
F_{onset} = \frac{\tanh(6(F_{onset1} - 1.35)) + 1}{2}
\tag{47}
$$

Similar to the original, the new $F_{onset}$ function is activated when $F_{onset1}$ exceeds unity; however, once a turbulent boundary layer has formed, $F_{onset1}$ remains active. This increases the stability of the model around laminar separation bubbles, where the solution can oscillate as the bubble grows and shrinks, but prevents the model from being able to predict relaminarization. For flight conditions typical of transport aircraft, we do not expect relaminarization to occur; however, the reader should be aware of this limitation.

#### 2. $F_{turb}$

In the SA-LM2015 model, the $F_{turb}$ function [Eq. (35)] remains active until a threshold eddy-viscosity ratio is reached, approximately a value of 4. Until this threshold value is reached, the destruction term remains active and interferes with the production term, slowing convergence and often producing oscillations. To remedy this, a new function was developed which deactivates the destruction term when $F_{onset}$ becomes active and rapidly deactivates $F_{turb}$ as turbulent eddy viscosity $R_T$ is produced. It is given as follows:

$$
F_{turb} = (1 - F_{onset})\exp(-R_T)
\tag{48}
$$

#### 3. $F_{length}$

The original $F_{length}$ function [Eq. (37)] consists of conditional statements and is nondifferentiable at the boundaries. In the current Paper, this is avoided by approximating $F_{length}$ using a modified, flexible version of a sigmoid function known as a generalized logistic function [63]. The new $F_{length}$ function is defined in the following and plotted in Fig. 2:

$$
F_{length1} = \exp(-3\cdot10^{-2}(\tilde{Re}_{\theta t} - 460))
\tag{49}
$$

$$
F_{length} = 44 - \frac{44 - (0.50 - 3\cdot10^{-4}(\tilde{Re}_{\theta t} - 596))}{(1 + F_{length1})^{1/6}}
\tag{50}
$$

> **Fig. 2** Original and smooth $F_{length}$ functions.
> *[Figure not reproducible — $F_{length}$ vs $\tilde{Re}_{\theta t}$, LM2015 (solid) vs sLM2015 (dotted), with inset zoom over $540 \le \tilde{Re}_{\theta t} \le 640$. Plotted data not available in text.]*

#### 4. $Re_{\theta c}$

In previous work, approximations of the critical momentum thickness Reynolds number $Re_{\theta c}$ as linear functions of $\tilde{Re}_{\theta t}$ are introduced [47,48,50]. In the current work, a sinusoidal approximation is developed in order to better represent the lower $\tilde{Re}_{\theta t}$ behavior of the original correlation [Eq. (38)]. The critical momentum-thickness Reynolds number is simplified from the original version as follows:

$$
Re_{\theta c} = 0.67\,\tilde{Re}_{\theta t} + 24\sin\!\left(\frac{\tilde{Re}_{\theta t}}{240} + 0.5\right) + 14
\tag{51}
$$

The original and smooth $Re_{\theta c}$ functions are plotted in Fig. 3. For large values of $\tilde{Re}_{\theta t}$, the new correlation produces smaller $Re_{\theta c}$ values relative to the original model; however, this difference has been accounted for in the model calibration.

> **Fig. 3** Original and smooth $Re_{\theta c}$ functions.
> *[Figure not reproducible — $Re_{\theta c}$ vs $\tilde{Re}_{\theta t}$, LM2015 (solid) vs sLM2015 (dotted), with inset zoom over $480 \le \tilde{Re}_{\theta t} \le 520$. Plotted data not available in text.]*

#### 5. $P_\gamma$ and $E_\gamma$

The maximum value of vorticity magnitude in the intermittency production and destruction terms $P_\gamma$ [Eq. (28)] and $E_\gamma$ [Eq. (34)] was limited to improve stability. Vorticity magnitude is included to obtain the correct dimensions for the source terms; however, it was found that large values of vorticity, especially near laminar separation bubbles, can lead to the solution oscillating. Limiting the maximum vorticity value in the intermittency source terms restricts the relative magnitude of the sources, improving the numerical behavior of the model. The limiting is a function of the freestream Mach number in order to account for the nondimensionalization of velocity by the speed of sound $a_\infty$ and the freestream Reynolds number to account for the scaling of vorticity in laminar boundary layers with respect to the square root of the Reynolds number. The Reynolds number is multiplied by the Mach number to recover the Reynolds number based on the freestream velocity. It is important to note that, while the presence of the Mach number is the result of the particular nondimensionalization used, the square root of the freestream Reynolds number based on the freestream velocity is a physical scaling that is independent of the nondimensionalization. Furthermore, the vorticity limiting is introduced only to improve the numerical behavior of the model and does not aim to improve its predictive capability. The new source terms are given as

$$
P_\gamma = c_{a1}\,F_{length}\,F_{onset}\left[\phi_{-300}\!\left(\Omega,\ \frac{M\sqrt{MRe}}{20}\right)\right]\sqrt{\gamma}\,(1 - c_{e1}\gamma)
\tag{52}
$$

$$
E_\gamma = c_{a2}\,F_{turb}\left[\phi_{-300}\!\left(\Omega,\ \frac{M\sqrt{MRe}}{20}\right)\right]\gamma\,(c_{e2}\gamma - 1)
\tag{53}
$$

and are illustrated in Figs. 4a–4d alongside the original SA-LM2015 source terms. The total source terms, $S_\gamma = P_\gamma - E_\gamma$, are plotted with fixed $F_{length}$ and vorticity magnitude $\Omega$ values of unity versus the transition criterion [Eq. (32)], intermittency $\gamma$, and the eddy-viscosity ratio $R_T$. The reformulated intermittency source terms provide a smoother design space, while the results in Sec. V demonstrate that the predicted transition locations with these modifications match well with the SA-LM2015 model and with experiment.

> **Fig. 4** Original and smooth intermittency source terms, $S_\gamma = P_\gamma - E_\gamma$, plotted vs the transition criterion [Eq. (32)], intermittency $\gamma$, and the eddy-viscosity ratio $R_T$, with $F_{length} = \Omega = 1$.
> *[Figure not reproducible — four 3D surface plots: (a) LM2015 at $R_T=0.1$; (b) sLM2015 at $R_T=0.1$; (c) LM2015 at $\gamma=0.5$; (d) sLM2015 at $\gamma=0.5$. Surface data not available in text.]*

### C. Momentum-Thickness Reynolds Number Empirical Correlations

Similar to the $F_{length}$ function, the $F(\lambda_\theta)$ correlations [Eq. (9)] are conditional, depending on the sign of the pressure gradient parameter $\lambda_\theta$. Here, $F(\lambda_\theta)$ is approximated using exponential penalty functions:

$$
F(\lambda_\theta)_1 = 1 + 0.275\left[1 - e^{-35\lambda_\theta}\right]e^{-[Tu/0.5]}
\tag{54}
$$

$$
F(\lambda_\theta)_2 = \phi_{300}(F(\lambda_\theta)_1, 1)
\tag{55}
$$

$$
F(\lambda_\theta)_3 = 1 - \left[-12.986\lambda_\theta - 123.66\lambda_\theta^2 - 405.689\lambda_\theta^3\right]e^{-[Tu/1.5]^{1.5}}
\tag{56}
$$

$$
F(\lambda_\theta) = \phi_{-300}\!\left(F(\lambda_\theta)_2, F(\lambda_\theta)_3\right)
\tag{57}
$$

The original and smooth $F(\lambda_\theta)$ functions are plotted in Fig. 5.

> **Fig. 5** Original and smooth $F(\lambda_\theta)$ functions for a turbulence intensity of 0.05%.
> *[Figure not reproducible — $F(\lambda_\theta)$ vs $\lambda_\theta$ over $-0.3 \le \lambda_\theta \le 0.2$, LM2015 (solid) vs sLM2015 (dotted), with inset zoom near $\lambda_\theta = 0$. Plotted data not available in text.]*

---

## IV. Solution Algorithm

### A. Parallel Newton–Krylov–Schur Flow Solver

An efficient, robust, and accurate flow solver is crucial to the performance of an optimization algorithm, as flow analysis is conducted many times over the course of the optimization process, with a variety of geometries expected. The flow solver used in the current Paper is a three-dimensional multiblock structured finite-difference solver developed by Hicken and Zingg [64] for the solution of the Euler equations and extended to the RANS equations coupled to the one-equation SA turbulence model by Osusky and Zingg [65]. The governing equations are spatially discretized using summation-by-parts operators, with simultaneous approximation terms applied to enforce boundary conditions and interblock coupling. To decrease the computational time required to complete a flow solution, the computational domain is decomposed into multiple blocks, resulting in multiblock structured grids, which allow for efficient parallel computations.

A fully coupled, implicit Newton–Krylov–Schur solution algorithm making use of a pseudotransient continuation strategy is applied to the set of nonlinear algebraic equations produced by the spatial discretization of the governing equations to drive the solution from an initial guess to a converged steady-state solution. The solution strategy consists of two phases, an approximate-Newton phase followed by an inexact-Newton phase. Globalization is achieved using the approximate-Newton phase, where an implicit Euler time marching strategy with local time linearization is applied using an approximate, analytically derived Jacobian. After the residual drops several orders of magnitude, Newton's method using a full second-order Jacobian, formed either analytically or using a forward difference, is used to converge the system to a residual norm of machine zero. The large system of linear equations generated at each iteration in both phases is solved using the GMRES linear solver [66], preconditioned using an approximate-Schur parallel preconditioner [67].

### B. Considerations for SA-LM2015 and SA-sLM2015 Transition Models

Unless suitable steps are taken to modify the basic algorithm, the implementation of the SA-LM2015 and SA-sLM2015 transition model equations in a Newton–Krylov solver can lead to numerical instability. The source terms in the $\gamma\text{-}\tilde{Re}_{\theta t}$ equations are large and highly nonlinear. This leads to numerical stiffness through the presence of restrictive source-term time scales [68,69]. Several authors [68,70,71] have suggested that fully implicit approaches can lead to numerical instabilities due to the formation of Jacobians with unfavorable characteristics. More specifically, it has been demonstrated that the implicit treatment of sources can lead to unbounded solution updates and that a modified explicit treatment of sources can improve solver performance [69]. Additionally, recent work has suggested that fully coupled approaches may lead to reduced performance due to increased numerical stiffness [72–74]. It has been shown that relaxing the coupling between the turbulence and mean flow, or transition and turbulence equations, can reduce this stiffness, leading to reduced computational cost [74] and improved convergence [72,73].

However, the rate of convergence of Newton's method is closely linked to the accuracy of the Jacobian. Neglecting terms or otherwise manipulating the Jacobian can prevent quadratic convergence [75]. Yildirim et al. [76] demonstrated that segregating the turbulence model in an approximate-Newton phase to globalize Newton's method can be effective for fully turbulent simulations. Similar approaches have been investigated for the SA-LM2015 and SA-sLM2015 models; however, this can result in the solution stalling in the approximate-Newton phase. The transition model equations require a tight coupling in order to achieve deep convergence, often requiring stricter linear solver tolerances relative to fully turbulent simulations. Moreover, the transition model source terms remain active throughout the inexact-Newton phase. Rather than modifying the Jacobian matrix, Lian et al. [69] demonstrated that restricting the local time step based on the nondimensional source-term number, defined as the product of the largest positive eigenvalue of the source-term Jacobian matrix and the local time step, can be an effective way of stabilizing large sources in a fully coupled, implicit algorithm. A source-term time step restriction based on the ideas developed by Lian et al. [69] has been implemented in the current Paper.

Additional modifications to the Newton–Krylov algorithm were implemented to improve robustness. These include equation and variable scaling measures, a solution update damping routine, and an unsteady line-search backtracking algorithm.

#### 1. Equation and Variable Scaling

The turbulence and transition model variables can vary significantly in magnitude from the mean-flow variables and do not contain inherent geometric scaling [75]. In addition, large off-diagonal entries are introduced through the linearization of the turbulence and transition model equations with respect to the mean-flow variables. If not accounted for, these factors can result in improper scaling of the linear system and poor performance of the linear solver. To address this, row and column scaling measures, including residual autoscaling, are used based on the work by Osusky and Zingg [55,65] and Chisholm and Zingg [75]. Instead of solving the linear system of equations directly, a scaled version of the equations is solved, given by

$$
S_a S_r\left(\frac{\mathcal{I}}{\Delta t} + \mathcal{A}^{(n)}\right)S_c S_c^{-1}\Delta \mathcal{Q}^{(n)} = -S_a S_r \mathcal{R}^{(n)}
\tag{58}
$$

where $S_a$ represents the residual autoscaling matrix and

$$
S_{ri} =
\begin{bmatrix}
J_i^{2/3} & & & & \\
& \ddots & & & \\
& & J_i^{2/3} & & \\
& & & \tilde{\nu}_{max}^{-1}J_i^{-1/3} & \\
& & & & \gamma_{max}^{-1}J_i^{-1/3} \\
& & & & & \tilde{Re}_{\theta t\,max}^{-1}J_i^{-1/3}
\end{bmatrix},
\qquad
S_{ci} =
\begin{bmatrix}
1 & & & & \\
& \ddots & & & \\
& & 1 & & \\
& & & \tilde{\nu}_{max} & \\
& & & & \gamma_{max} \\
& & & & & \tilde{Re}_{\theta t\,max}
\end{bmatrix}
$$

are the $i$th blocks of the diagonal row and column scaling matrices, respectively. Through a thorough algorithm optimization procedure, Osusky and Zingg determined an optimal normalization value for $\tilde{\nu}_{max}$ of $10^3$, which accounts for the maximum turbulence value that is likely to be encountered in the flow solve. Following a similar approach, $\gamma_{max}$ and $\tilde{Re}_{\theta t\,max}$ values of 10 and $10^4$, respectively, were selected. The partially scaled residual $S_r\mathcal{R}^{(n)}$ is used when monitoring convergence. Full details on the autoscaling, row scaling, and column scaling measures can be found in the work by Osusky and Zingg.

#### 2. Solution Update

The transition model equations can become unstable if the intermittency and transition onset momentum-thickness Reynolds number are not bounded. Solution-limited time stepping methods, such as those developed by Chisholm and Zingg [75] and Lian et al. [69,77], were investigated. However, these methods can slow convergence for cases with stable updates and lead to restrictive time steps for variables approaching zero. To ensure stable physical transition model updates, a solution update damping procedure was implemented, as described in Algorithm 2. Clipping the transition model variables at these bounds was also investigated; however, placing a hard limit on variable bounds can potentially lead to stalling of the nonlinear solver. It is important to note that with the implementation of the source-term time stepping strategy, which prevents unstable solution updates, in combination with a more dissipative first-order upwind discretization strategy for the SA, SA-LM2015, and SA-sLM2015 equations, the upper and lower bounds specified in Algorithm 2 are rarely reached. Furthermore, steady-state solutions for intermittency respect the bounds of 0.02 and 1 implicitly defined in the source terms. In addition to the solution update damping, a backtracking line-search method was implemented based on the work of Modisette [78].

**Algorithm 2: Transition Model Solution Update Damping**

```
1:  θ_fac ← 0.99                                          ▷ damping factor
2:  for j = 1, n_nodes do                                 ▷ loop over local nodes
3:      for l = 1, n_var do                               ▷ loop over variables
4:          Q^(n+1)_{j,l} = Q^(n)_{j,l} + ΔQ^(n)_{j,l}    ▷ update solution
5:          m = 0                                          ▷ initialize count variable
6:          while Q^(n+1)_{j,7} > 2 or Q^(n+1)_{j,7} < 10^-10 do   ▷ limit intermittency
7:              m = m + 1                                  ▷ count number of damping iterations
8:              Q^(n+1)_{j,7} = Q^(n)_{j,7} + (θ_fac)^m ΔQ^(n)_{j,7}   ▷ damp update
9:          m = 0                                          ▷ initialize count variable
10:         while Q^(n+1)_{j,8} < 20 do                    ▷ limit R̃e_θt
11:             m = m + 1                                  ▷ count number of damping iterations
12:             Q^(n+1)_{j,8} = Q^(n)_{j,8} + (θ_fac)^m ΔQ^(n)_{j,8}   ▷ damp update
```

#### 3. Source-Term Time Stepping

Lian et al. [69] demonstrated that an implicit scheme becomes unstable with source-term time step values greater than unity. Maintaining a source-term time step less than unity can prevent unbounded solution updates, where the source-term time step is determined by the product of the largest positive eigenvalue of the source-term Jacobian matrix $\mathcal{A}_{source}$ and the local time step $\lambda_{source}\Delta t^{(n)}_{j,k,m}$. For the current Paper, the eigenvalues of the block diagonal source-term Jacobian are determined using a QR algorithm. The local time step is restricted to ensure the following source-term time step condition is met,

$$
\mathcal{A}_{source} =
\begin{bmatrix}
\dfrac{\partial S_{\tilde{\nu}}}{\partial \tilde{\nu}} & \dfrac{\partial S_{\tilde{\nu}}}{\partial \gamma} & \dfrac{\partial S_{\tilde{\nu}}}{\partial \tilde{Re}_{\theta t}} \\[2ex]
\dfrac{\partial S_{\gamma}}{\partial \tilde{\nu}} & \dfrac{\partial S_{\gamma}}{\partial \gamma} & \dfrac{\partial S_{\gamma}}{\partial \tilde{Re}_{\theta t}} \\[2ex]
\dfrac{\partial S_{\tilde{Re}_{\theta t}}}{\partial \tilde{\nu}} & \dfrac{\partial S_{\tilde{Re}_{\theta t}}}{\partial \gamma} & \dfrac{\partial S_{\tilde{Re}_{\theta t}}}{\partial \tilde{Re}_{\theta t}}
\end{bmatrix},
\qquad
\lambda_{source}\Delta t^{(n)}_{j,k,m} \le 0.9
\tag{59}
$$

where $\mathcal{S}$ represents the local source-term vector, in this case consisting of the SA-LM2015, SA-sLM2015, and SA turbulence model source terms. Several source-term time step values were investigated; however, 0.9 was determined to be an effective balance between speed and robustness. The spatially varying time step and the time step ramping procedure are described in [65].

As demonstrated by Lian et al. [69], source-term eigenvalues can remain large throughout convergence. This can cause the local time step to remain excessively small throughout the solution process, delaying convergence. However, for most cases, the large sources are only destabilizing in the early stages of convergence. To accelerate convergence, a switch was developed for deactivating the source-term time stepping in the inexact-Newton phase. The source-term time stepping is active throughout the approximate-Newton phase and is deactivated after five successive inexact-Newton iterations for which backtracking is not active. Source-term time stepping is reactivated if backtracking is required, or if the relative residual drop [Eq. (60)] increases beyond the residual drop tolerance for the switch from the approximate- to the inexact-Newton phase, which is typically a value of $10^{-5}$.

---

## V. Results

Both transition models, SA-LM2015 and SA-sLM2015, are validated in the following section using a range of two- and three-dimensional transition test cases. These test cases include natural and separation-induced transition as well as transition due to crossflow instabilities, with the results compared to experimental values. Convergence histories are illustrated by the L2-norm of the fully coupled residual, including the row (equation) scaling measure presented in Sec. IV, normalized by the initial residual,

$$
\mathcal{R}^{(n)}_d = \frac{\|\mathcal{R}^{(n)}\|_2}{\|\mathcal{R}^{(0)}\|_2}
\tag{60}
$$

versus the computational time in seconds. All simulations were conducted on Intel Skylake processors with a core clock speed of 2.4 GHz and one Message Passing Interface thread per computational block. All results were obtained using a second-order discretization for the mean-flow equations with matrix-based dissipation using $V_l$ and $V_n$ values of 0 and a fourth-difference dissipation coefficient of 0.04 [79]. First-order upwinding was used for the SA, SA-LM2015, and SA-sLM2015 models. The same solution algorithm parameters were used for the two- and three-dimensional cases, respectively, demonstrating the robustness of the solver. These settings are similar to the values recommended by Osusky and Zingg [55], with the settings for the source-term time stepping strategy and scaling measures described previously in Sec. IV.

### A. NLF0416 General Aviation Airfoil

The NLF0416 airfoil is a natural-laminar-flow general aviation airfoil designed by Somers [80]. Experimental results of the NLF0416 airfoil were obtained in the low turbulence pressure tunnel (LTPT) at NASA Langley Research Center [80]. Transition locations were measured for several angles of attack at a Mach number of 0.1 and a Reynolds number of $4.0\times10^6$. The transition mechanisms for these conditions are natural and separation-induced transition. The turbulence intensity in the wind tunnel used for the experimental study was not specified; however, a value of $Tu = 0.15\%$ was assumed in order to be consistent with results obtained using the Transition Modelling and Predictive Capabilities Seminar (TCMPS) guidelines [81].

The coarse, medium, and fine NLF0416 family of grids provided by the TMPCS were used in the current work [81]. The dimensions of the structured C grids are presented in Table 1. A grid-refinement study was conducted for the 0 deg angle of attack case using both transition models and is presented in Fig. 6. The refinement study demonstrates a clear second-order convergence trend, as the inverse of the number of grid points $N^{-1}$ is proportional to the square of the mesh spacing, although they converge to slightly different lift and drag values due to a difference in predicted transition location. Figure 6 also illustrates the computational results obtained using the SA-LM2015 and SA-sLM2015 transition models on the fine grid level compared to experimental values. The results demonstrate that both the SA-LM2015 model and the SA-sLM2015 model do a reasonable job of matching the forces and transition locations from experiment. Figure 6e illustrates that both models underpredict drag at low angles of attack; however, this is consistent with results obtained from previous studies using variations of the $\gamma\text{-}\tilde{Re}_{\theta t}$ model [48,82,83] and does not appear to be caused by incorrect transition point locations.

**Table 1 — TMPCS NLF0416 structured C-grid dimensions [81]**

| Grid level | Dimensions | $\Delta s$ (chord) | Average/maximum $y^+$ |
|---|---|---|---|
| Coarse | 529 × 73 | $4.47\times10^{-6}$ | 0.57 / 1.19 |
| Medium | 705 × 97 | $3.25\times10^{-6}$ | 0.51 / 0.89 |
| Fine | 1057 × 145 | $2.10\times10^{-6}$ | 0.26 / 0.45 |

> **Fig. 6** NLF0416 grid convergence at 0 deg angle of attack and angle of attack sweep simulated using the SA-LM2015 and SA-sLM2015 transition models at $M = 0.1$, $Re = 4.0\times10^6$, and $Tu = 0.15\%$. The fine grid level was used for the sweep, with the experimental results obtained at the LTPT wind tunnel [80].
> *[Figure not reproducible — 6 subplots: (a) skin friction grid convergence ($C_f$ vs $x/c$, upper/lower), (b) grid convergence for lift and drag ($C_d$, $C_l$ vs $N^{-1}$), (c) drag polar, (d) lift curve, (e) drag curve, (f) transition locations vs $C_l$. Plotted data not available in text.]*

The residual convergence histories for these cases are presented in Fig. 7. For these simulations, the grids were split into eight computational blocks. The convergence histories demonstrate the improved numerical behavior of the SA-sLM2015 model. At angles of attack of −6 and −4 deg, a large separation bubble forms on the lower surface of the airfoil, increasing the stiffness of the linear system. This is evident in the decreased slope of the residual curve for the SA-sLM2015 model, as each outer inexact-Newton iteration required more inner GMRES linear solver iterations. It should also be noted that small trailing edge spacings on the C grids provided by the TMPCS propagated to the far field of the domain with little splaying, further increasing the stiffness of the linear system and increasing the computational cost of the simulations.

> **Fig. 7** NLF0416 residual convergence histories for the SA-LM2015 and SA-sLM2015 transition models on the fine grid level split into eight computational blocks. The airfoil was simulated at $M = 0.1$, $Re = 4.0\times10^6$, and $Tu = 0.15\%$ and over a range of angles of attack.
> *[Figure not reproducible — 7 subplots (relative residual $\mathcal{R}_d$ vs time) at angles of attack of −6, −4, −2, 0, 2, 4, 6 deg. Plotted data not available in text.]*

### B. S809 Wind Turbine Airfoil

The S809 airfoil is a natural-laminar-flow wind turbine airfoil designed by Somers [84] and studied in the low-turbulence wind tunnel of the Delft University of Technology Low Speed Laboratory. Experimental results for the S809 airfoil were obtained at a Reynolds number of $2.0\times10^6$ at angles of attack ranging from 0 to 9 deg. The freestream turbulence intensity in the wind tunnel varied between $Tu = 0.02\%$ and $Tu = 0.04\%$. For a Reynolds number of $2.0\times10^6$ and an angle of attack of 0 deg, Somers noted that the primary transition mechanism on both the upper and lower surfaces of the airfoil was separation-induced transition due to the presence of large laminar separation bubbles. As the angle of attack increased, the separation bubble on the upper surface decreased in length, disappearing completely at an angle of attack of 5.13 deg, with the transition location moving upstream rapidly on the upper surface and the lower surface separation bubble fixing transition near midchord over the range of angles of attack. Numerical results were obtained using a Mach number of 0.10, Reynolds number of $2.0\times10^6$, and a freestream turbulence intensity of 0.07%, as recommended by the TCMPS guidelines [81].

Computational results were obtained using the S809 grid families supplied through the TCMPS, which were generated using the same conformal mapping method as for the NLF0416 [81], with the coarse, medium, and fine grid dimensions presented in Table 2. A grid-refinement study was conducted for the S809 airfoil at 6 deg angle of attack and is plotted in Fig. 8. This condition was chosen as it represents the point beyond the upper corner of the drag bucket where the transition location is rapidly moving, which is expected to amplify the effects of grid refinement on the predictions [81]. Figure 8 also illustrates the forces and transition locations obtained using the SA-LM2015 and SA-sLM2015 transition models on the fine grid level. Although the transition locations match experiment well, both models underpredict drag before and after the upper corner of the drag bucket. Again, this trend has been identified by previous studies using variations of the $\gamma\text{-}\tilde{Re}_{\theta t}$ model [48,82,83], as well as studies using the AFT transition model [35,82], and does not appear to be caused by inaccurate transition prediction.

**Table 2 — TMPCS S809 structured C-grid dimensions [81]**

| Grid level | Dimensions | $\Delta s$ (chord) | Average/maximum $y^+$ |
|---|---|---|---|
| Coarse | 529 × 73 | $8.90\times10^{-6}$ | 0.60 / 1.73 |
| Medium | 705 × 97 | $6.47\times10^{-6}$ | 0.52 / 0.81 |
| Fine | 1057 × 145 | $4.19\times10^{-6}$ | 0.28 / 0.54 |

> **Fig. 8** S809 grid convergence at 6 deg angle of attack and angle of attack sweep simulated using the SA-LM2015 and SA-sLM2015 transition models at $M = 0.1$, $Re = 2.0\times10^6$, and $Tu = 0.07\%$. The fine grid level was used for the sweep, with the experimental results obtained at the Delft University of Technology wind tunnel [84].
> *[Figure not reproducible — 6 subplots analogous to Fig. 6 (skin friction, lift/drag grid convergence, drag polar, lift curve, drag curve, transition locations), compared with TU Delft experiment. Plotted data not available in text.]*

The residual convergence histories for these cases are presented in Fig. 9 for computations using eight processors. Similar to the NLF0416 airfoil, the large separation bubbles at low angles of attack led to stiff linear systems, which slowed convergence. For all angles of attack, the SA-sLM2015 model displays substantially improved convergence. In particular, for the 9 deg angle of attack case, the SA-LM2015 model was unable to converge due to the unsteady interaction of the transition model with the laminar separation bubble on the upper surface of the airfoil. This behavior was also identified by Denison and Pulliam [82] using the LM2015 model. Limiting vorticity in the transition model source terms prevented this oscillation without affecting the accuracy of the solution.

> **Fig. 9** S809 residual convergence histories for the SA-LM2015 and SA-sLM2015 transition models on the fine grid level split into eight computational blocks. The airfoil was simulated at $M = 0.1$, $Re = 2.0\times10^6$, and $Tu = 0.07\%$ and over a range of angles of attack.
> *[Figure not reproducible — 9 subplots (relative residual $\mathcal{R}_d$ vs time) at angles of attack of 1 through 9 deg. Plotted data not available in text.]*

### C. NASA NLF2-0415 Infinite Swept Wing

The NASA NLF2-0415 Infinite Swept Wing consists of an NLF2-0415 airfoil extruded with a 45 deg sweep angle. Experimental measurements were obtained at an angle of attack of −4 deg, over a range of Mach and Reynolds numbers [85], and including a variety of surface roughness levels [86]. Transition locations were determined using the naphthalene visualization technique. Transition locations for a roughness level $h$ of 3.3 μm, representing a painted surface, were obtained using the SA-LM2015 and SA-sLM2015 models and are compared with experimental values [86] in Fig. 10. The Reynolds number based on freestream velocity magnitude and chord length was varied from $1.8\times10^6$ to $3.5\times10^6$. The freestream turbulence intensity for the wind tunnel was measured as 0.20%. Transition locations are presented only for the upper surface of the wing, as experimental results are not available for the lower surface.

> **Fig. 10** Experimental and computational results for the NLF2-0415 infinite swept wing [86] simulated using the SA-LM2015 and SA-sLM2015 transition models at $M = 0.15$, $\alpha = -4°$, and $Tu = 0.20\%$ and over a range of Reynolds numbers from $1.8\text{–}3.5\times10^6$.
> *[Figure not reproducible — (a) transition location $x/c$ vs Reynolds number, LM2015 / sLM2015 / experiment; (b) skin friction profiles $C_f$ vs $x/c$. Plotted data not available in text.]*

The computational grid used to obtain these results consists of 270,336 nodes split into 128 blocks, with 90,112 nodes per cross-section and an off-wall spacing of $2.5\times10^{-6}$ chords, resulting in average and maximum $y^+$ values of 0.21 and 0.46, respectively. The analytical full Jacobian was used in the inexact-Newton phase for this case and the Technical University of Braunschweig Sickle Wing. Periodic boundary conditions were used to simulate the infinite swept wing. The results in Fig. 10 demonstrate that both SA-LM2015 and SA-sLM2015 transition models predict transition accurately, while plots of the residual convergence histories in Fig. 11 demonstrate that the SA-sLM2015 model provides significantly improved convergence to steady state.

The results demonstrate that the SA-LM2015 model is able to closely predict the location of transition for all Reynolds numbers except for $2\times10^6$, where there is an error of approximately 14% between the experimental and predicted values. This error is also present in the work of Jung and Baeder [48], who coupled a semilocal correlation-based transition model formulation with the SA model using the crossflow correlations developed by Muller and Herbst [42] and in the original model developed by Langtry et al. [43]. In the work produced by Radeztsky et al. [86], the results demonstrate that there is significant experimental error associated with this case, which may help to explain this discrepancy.

> **Fig. 11** NLF2-0415 residual convergence histories for the SA-LM2015 and SA-sLM2015 transition models, with the grid split into 128 computational blocks. The infinite swept wing was simulated at $M = 0.15$, $\alpha = -4°$, and $Tu = 0.20\%$ and over a range of Reynolds numbers.
> *[Figure not reproducible — 6 subplots (relative residual $\mathcal{R}_d$ vs time) at $Re = 1.8, 2.0, 2.2, 2.5, 3.0, 3.5\times10^6$. Plotted data not available in text.]*

It is important to note that the choice of numerical dissipation strongly affects the accuracy of the crossflow transition model. Both the scalar dissipation model developed by Jameson et al. [87] and later refined by Pulliam [88] and the matrix-based dissipation model developed by Swanson and Turkel [79] were investigated for this case. On the mesh used, the scalar dissipation model drastically underpredicts the crossflow strength and helicity values. This was determined to be due to the increased dissipation of the scalar dissipation model in and above the laminar boundary layer, which affected the vorticity profile within the boundary layer.

### D. Technical University of Braunschweig Sickle Wing

The Technical University of Braunschweig Sickle Wing [89] was developed as a test case for crossflow transition on a complex geometry. The model consists of five main sections designed to create three-dimensional boundary layers with increasing crossflow toward the wing tip. The first section raises the model away from the turbulent boundary layer of the wind tunnel wall, with the three main wing sections swept at 30, 45, and 55 deg with distinct kinks in the leading-edge sweep. As a result, large spanwise gradients are present, and the assumptions of linear local stability theory are challenged [89]. The wing is capped with a rounded wing tip. At the experimental flow conditions, the upper surface of the wing has extended regions of strong favorable pressure gradients with crossflow instabilities developing at the section kinks. This results in mixed mode transition scenarios with Tollmien–Schlichting and crossflow instabilities triggering transition, except on the upper surface near the root of the wing, where a laminar separation bubble triggers transition at approximately 80% chord.

The Sickle Wing was simulated at a Reynolds number of $2.75\times10^6$, Mach number of 0.16, and −2.6 deg angle of attack. The rms surface roughness of the experimental model was measured to be 1.47 μm with a five-point average of the maximum peak-to-peak height of 9.78 μm [89]. To be consistent with results obtained using the LM2015 model [43,44], the surface roughness for the simulations was set to the peak-to-peak height value of 9.78 μm with a freestream turbulence intensity of 0.17%, which is consistent with the turbulence intensity measured in the wind tunnel [89].

Three levels of O-O topology structured multiblock grids were simulated with the grid details presented in Table 3, including the number of spanwise and chordwise nodes on each of the upper and lower surfaces. Each grid level was split into 592 computational blocks. Figure 12 illustrates the skin friction profiles on the upper and lower surfaces of the Sickle Wing using both the SA-LM2015 and SA-sLM2015 transition models for all three grid levels overlaid with experimental transition locations [89]. Both models predict the transition fronts reasonably well.

**Table 3 — Sickle Wing structured O–O grid dimensions**

| Grid level | Total node count, $N$ | $N$ span | $N$ chord | $\Delta s$ (chord) | Average/maximum $y^+$ |
|---|---|---|---|---|---|
| Coarse | $2.32\times10^6$ | 171 | 111 | $7.81\times10^{-6}$ | 0.70 / 2.64 |
| Medium | $5.19\times10^6$ | 217 | 139 | $6.02\times10^{-6}$ | 0.54 / 1.55 |
| Fine | $10.19\times10^6$ | 273 | 175 | $4.29\times10^{-6}$ | 0.39 / 1.13 |

> **Fig. 12** Skin friction profiles on the upper and lower surfaces of the Sickle Wing for all grid levels overlaid with experimental transition locations [89]. Results were obtained using the SA-LM2015 and SA-sLM2015 transition models at $M = 0.16$, $\alpha = -2.6°$, $Tu = 0.17\%$, $h = 9.78$ μm, and a Reynolds number of $2.75\times10^6$.
> *[Figure not reproducible — 12-panel $C_f$ contour montage (coarse/medium/fine × SA-LM2015/SA-sLM2015 × upper/lower surface), $C_f$ colour scale 0.001–0.006. Contour data not available in text.]*

As the grid is refined, the transition front on the lower surface, which is dominated by natural transition, moves slightly downstream and agrees well with experiment. The transition front on the upper surface, where crossflow-induced transition is the primary mechanism, appears more sensitive to the grid refinement. The transition locations at the kink regions move upstream as the stationary crossflow vortices are better resolved. The lift and drag grid convergence for both models is presented in Fig. 13b. The pressure coefficient profiles are extracted at the midspan of each sweep section for the fine grid, corresponding to spanwise sections at 25%, 55%, and 85% span, and are compared with the experimental profiles presented by Kruse et al. [90] in Fig. 13a. The pressure coefficient profiles from the simulations lie on top of one another, and the results demonstrate the wiggles observed in experiment at approximately $x/c = 0.35$, corresponding to natural transition on the lower surface of the Sickle Wing, are accurately reproduced [90]. There is a slight deviation between the simulations and experiment at 25% span on the upper surface of the Sickle Wing where both transition models, SA-LM2015 and SA-sLM2015, predict crossflow-induced transition upstream of the transition location in the experiment, which is triggered at $x/c = 0.80$ by a laminar separation bubble [89,90]. This deviation is also illustrated in Figs. 12a–12f and is consistent with results obtained using the LM2015 transition model [44].

> **Fig. 13** Experimental and fine grid pressure coefficient profiles and grid convergence for the SA-LM2015 and SA-sLM2015 transition models. The Sickle Wing was simulated at $M = 0.16$, $\alpha = -2.6°$, $Tu = 0.17\%$, and a Reynolds number of $2.75\times10^6$ [90].
> *[Figure not reproducible — (a) $C_p$ vs $x/c$ at 25%/55%/85% span, simulation vs experiment; (b) grid convergence for lift and drag ($C_D$, $C_L$ vs $N^{-2/3}$). Plotted data not available in text.]*

The relative residual and drag coefficient convergence histories are presented in Fig. 14. The flow complexity in this case slows convergence for both models; however, the smooth-variant displays improved residual convergence behavior in Figs. 14a–14c. The drag convergence, the magnitude of the difference between the current drag coefficient, and the converged drag coefficient, plotted in Figs. 14d–14f show that the SA-sLM2015 model leads to faster convergence of the drag coefficient.

> **Fig. 14** Relative residual and drag coefficient convergence histories for the SA-LM2015 and SA-sLM2015 transition models on all grid levels, with the grids split into 592 computational blocks. The Sickle Wing was simulated at $M = 0.16$, $\alpha = -2.6°$, $Tu = 0.17\%$, $h = 9.78$ μm, and a Reynolds number of $2.75\times10^6$.
> *[Figure not reproducible — 6 subplots: (a–c) relative residual $\mathcal{R}_d$ vs time for coarse/medium/fine; (d–f) drag convergence (counts) vs time for coarse/medium/fine. Plotted data not available in text.]*

---

## VI. Conclusions

A novel coupling of the $\gamma\text{-}\tilde{Re}_{\theta t}$ LM2015 transition model with helicity-based crossflow correlations [43] and the one-equation SA turbulence model has been developed and designated SA-LM2015. The modifications maintain the fully local formulation pioneered by Langtry and Menter [29]. Smooth approximations to nondifferentiable and stiff functions in the SA-LM2015 model are introduced in order to eliminate discontinuities and improve numerical behavior, with the smooth model designated SA-sLM2015. Deep convergence is achieved through the use of a fully coupled, fully implicit Newton–Krylov algorithm with modifications to the solution algorithm, including a source-term time stepping strategy, introduced to address the large source terms in the turbulence and transition model equations. Several of the strategies presented, including the smoothing techniques and source-term time stepping, could be useful in the context of other transition models, specifically transition models which use similar intermittency transport equations, such as the following [35,39,40,48,91].

Simulations of the S809 and NLF0416 natural-laminar-flow airfoil test cases demonstrate that both models, SA-LM2015 and SA-sLM2015, are able to predict natural and separation-induced transition accurately, with the SA-sLM2015 model displaying improved numerical behavior. Specifically, the new smooth formulations improve convergence, and combined with the new source-term formulation, including limiting vorticity magnitude, oscillations near laminar separation bubbles produced by the LM2015 and SA-LM2015 models are eliminated. Three-dimensional simulations of the NLF2-0415 infinite swept wing and Sickle Wing test cases demonstrate that the model is able to capture the effects of crossflow accurately.

The SA-LM2015 and SA-sLM2015 transition models are validated in the current Paper using popular transition test cases under subsonic conditions. Future work will involve applying both transition models, SA-LM2015 and SA-sLM2015, to transonic transition test cases in order to confirm that the trends produced in subsonic flow extend to higher Mach numbers.

---

## Acknowledgments

This Paper was partially funded by the NASA Transformational Tools and Technologies and Advanced Air Transport Technology projects, the Natural Sciences and Engineering Research Council, and the University of Toronto. Computations were performed on the Niagara supercomputer at the SciNet High performance computing Consortium, a part of Compute Canada. The authors gratefully acknowledge the input and support provided by Masayuki Yano, Thomas Reist, Gaetan Kenway, Thomas Pulliam, James Coder, and Cetin Kiris.

---

## Footnotes

‡ Data available online at https://overflow.larc.nasa.gov/home/users-manual-for-overflow-2-2/ [retrieved 10 August 2020].

---

## References

[1] Green, J. E., "Civil Aviation and the Environment: The Next Frontier for the Aerodynamicist," *Aeronautical Journal*, Vol. 110, No. 1110, 2006, pp. 469–486. https://doi.org/10.1017/S0001924000001378

[2] Bushnell, D. M., "Overview of Aircraft Drag Reduction Technology," *Special Course on Skin Friction Drag Reduction*, AGARD Rept. 786, 1992.

[3] Malik, M. R., Crouch, J. D., Saric, W. S., Lin, J. C., and Whalen, E. A., "Application of Drag Reduction Techniques to Transport Aircraft," *Encyclopedia of Aerospace Engineering*, Wiley, Hoboken, NJ, 2015, pp. 1–10. https://doi.org/10.1002/9780470686652.eae1013

[4] Crouch, J. D., "Boundary-Layer Transition Prediction for Laminar Flow Control," *45th AIAA Fluid Dynamics Conference*, AIAA Paper 2015-2472, 2015. https://doi.org/10.2514/6.2015-2472

[5] Campbell, R. L., and Lynde, M. N., "Natural Laminar Flow Design for Wings with Moderate Sweep," *34th AIAA Applied Aerodynamics Conference*, AIAA Paper 2016-4326, 2016. https://doi.org/10.2514/6.2016-4326

[6] Lynde, M. N., Campbell, R. L., Rivers, M. B., Viken, S. A., Chan, D. T., Watkins, N. A., and Goodliff, S. L., "Preliminary Results from an Experimental Assessment of a Natural Laminar Flow Design Method," *2019 AIAA Aerospace Sciences Meeting*, AIAA Paper 2019-2298, 2019. https://doi.org/10.2514/6.2019-2298

[7] Lynde, M. N., Campbell, R. L., and Viken, S. A., "Additional Findings from the Common Research Model Natural Laminar Flow Wind Tunnel Test," *AIAA Aviation 2019 Forum*, AIAA Paper 2019-3292, 2019. https://doi.org/10.2514/6.2019-3292

[8] Kenway, G. K. W., and Martins, J. R. R. A., "Multipoint High-Fidelity Aerostructural Optimization of a Transport Aircraft Configuration," *Journal of Aircraft*, Vol. 51, No. 1, 2014, pp. 144–160. https://doi.org/10.2514/1.C032150

[9] Reist, T. A., and Zingg, D. W., "High-Fidelity Aerodynamic Shape Optimization of a Lifting-Fuselage Concept for Regional Aircraft," *Journal of Aircraft*, Vol. 54, No. 3, 2017, pp. 1085–1097. https://doi.org/10.2514/1.C033798

[10] Reist, T. A., Zingg, D. W., Rakowitz, M., Potter, G., and Banerjee, S., "Multifidelity Optimization of Hybrid Wing-Body Aircraft with Stability and Control Requirements," *Journal of Aircraft*, Vol. 56, No. 2, 2019, pp. 442–456. https://doi.org/10.2514/1.C034703

[11] Reed, H. L., and Saric, W. S., "Transition Mechanisms for Transport Aircraft," *38th Fluid Dynamics Conference and Exhibit*, AIAA Paper 2008-3743, 2008. https://doi.org/10.2514/6.2008-3743

[12] Arnal, D., Casalis, G., and Houdeville, R., "Practical Transition Prediction Methods: Subsonic and Transonic Flows," *VKI Lecture Series: Advances in Laminar-Turbulent Transition Modeling*, von Karman Inst. for Fluid Dynamics, Belgium, 2009, pp. 1–34.

[13] Aupoix, B., Arnal, D., Bezard, H., Chaouat, B., and Chedevergne, F., "Transition and Turbulence Modeling," *AerospaceLab*, No. 2, 2011, pp. 1–13.

[14] Pasquale, D., Rona, A., and Garrett, S. J., "A Selective Review of CFD Transition Models," *39th AIAA Fluid Dynamics Conference*, AIAA Paper 2009-3812, 2009. https://doi.org/10.2514/6.2009-3812

[15] Savill, A. M., "By-Pass Transition Using Conventional Closures," *Closure Strategies for Turbulent and Transitional Flows*, Vol. 17, Cambridge Univ. Press, Cambridge, England, U.K., 2002, pp. 464–492. https://doi.org/10.1017/CBO9780511755385.019

[16] Halila, G. L. O., Chen, G., Shi, Y., Fidkowksi, K. J., Martins, J. R. R. A., and de Mendonca, M. T., "High-Reynolds Number Transitional Flow Simulation via Parabolized Stability Equations with an Adaptive RANS Solver," *Aerospace Science and Technology*, Vol. 91, Aug. 2019, pp. 321–336. https://doi.org/10.1016/j.ast.2019.05.018

[17] Drela, M., and Giles, M. B., "Viscous-Inviscid Analysis of Transonic and Low Reynolds Number Airfoils," *AIAA Journal*, Vol. 25, No. 10, 1987, pp. 1347–1355. https://doi.org/10.2514/3.9789

[18] Rashad, R., and Zingg, D. W., "Aerodynamic Shape Optimization for Natural Laminar Flow Using a Discrete-Adjoint Approach," *AIAA Journal*, Vol. 54, No. 11, 2016, pp. 3321–3337. https://doi.org/10.2514/1.J054940

[19] Mayda, E., "Boundary Layer Transition Prediction for Reynolds-Averaged Navier-Stokes Methods," Ph.D. Thesis, Univ. of California, Davis, Davis, CA, 2007.

[20] Stock, H. W., "$e^N$ Transition Prediction in Three-Dimensional Boundary Layers on Inclined Prolate Spheroids," *AIAA Journal*, Vol. 44, No. 1, 2006, pp. 108–118. https://doi.org/10.2514/1.16026

[21] Krimmelbein, N., and Radespiel, R., "Transition Prediction for Three-Dimensional Flows Using Parallel Computation," *Computers and Fluids Journal*, Vol. 38, No. 1, 2009, pp. 121–136. https://doi.org/10.1016/j.compfluid.2008.01.004

[22] Arnal, D., Habiballah, M., and Coutols, E., "Laminar Instability Theory and Transition Criteria in Two and Three-Dimensional Flow," *La Recherche Aérospatiale*, No. 2, 1984, pp. 125–143.

[23] Krumbein, A., Krimmelbein, N., and Schrauf, G., "Automatic Transition Prediction in Hybrid Flow Solver, Part 1: Methodology and Sensitivities," *Journal of Aircraft*, Vol. 46, No. 4, 2009, pp. 1176–1190. https://doi.org/10.2514/1.39736

[24] Krumbein, A., Krimmelbein, N., and Schrauf, G., "Automatic Transition Prediction in Hybrid Flow Solver, Part 2: Practical Application," *Journal of Aircraft*, Vol. 46, No. 4, 2009, pp. 1191–1199. https://doi.org/10.2514/1.39738

[25] Moens, F., Perraud, J., Krumbein, A., Toulorge, T., Iannelli, P., and Hanifi, A., "Transition Prediction and Impact on a Three-Dimensional High-Lift-Wing Configuration," *Journal of Aircraft*, Vol. 45, No. 5, 2008, pp. 1751–1766. https://doi.org/10.2514/1.36238

[26] Streit, T., Horstmann, H., Shraut, G., Hein, S., Fey, U., Egmai, Y., Perraud, J., Salah El Din, I., Cella, U., and Quest, J., "Complementary Numerical and Experimental Data Analysis of the ETW Telfona Pathfinder Wing Transition Tests," *49th Aerospace Sciences Meeting and Exhibit*, AIAA Paper 2011-881, Jan. 2011. https://doi.org/10.2514/6.2011-881

[27] Shi, Y., Gross, R., Mader, C. A., and Martins, J., "Transition Prediction Based on Linear Stability Theory with the RANS Solver for Three-Dimensional Configurations," *2018 AIAA Aerospace Sciences Meeting*, AIAA Paper 2018-0819, 2018. https://doi.org/10.2514/6.2018-0819

[28] Seyfert, C., and Krumbein, A., "Evaluation of a Correlation-Based Transition Model and Comparison with the $e^N$ Method," *Journal of Aircraft*, Vol. 49, No. 6, 2012, pp. 1765–1773. https://doi.org/10.2514/1.C031448

[29] Langtry, R. B., and Menter, F. R., "Correlation-Based Transition Modeling for Unstructured Parallelized Computational Fluid Dynamics Codes," *AIAA Journal*, Vol. 47, No. 12, 2009, pp. 2894–2906. https://doi.org/10.2514/1.42362

[30] Coder, J. G., and Maughmer, M. D., "A CFD-Compatible Transition Model Using an Amplification Factor Transport Equation," *51st AIAA Aerospace Sciences Meeting Including the New Horizons Forum and Aerospace Exposition*, AIAA Paper 2013-0253, 2013. https://doi.org/10.2514/6.2013-253

[31] Coder, J. G., and Maughmer, M. D., "Computational Fluid Dynamics Compatible Transition Modeling Using an Amplification Factor Transport Equation," *AIAA Journal*, Vol. 52, No. 11, 2014, pp. 2506–2512. https://doi.org/10.2514/1.J052905

[32] Coder, J. G., "Enhancement of the Amplification Factor Transport Transition Modeling Framework," *55th AIAA Aerospace Sciences Meeting*, AIAA Paper 2017-1709, 2017. https://doi.org/10.2514/6.2017-1709

[33] Stefanski, D., Glasby, R., Erwin, J. T., and Coder, J. G., "Development of a Predictive Capability for Laminar-Turbulent Transition in HPCMP CREATE™-AV Kestrel Component COFFE Using the Amplification Factor Transport Model," *2018 AIAA Aerospace Sciences Meeting*, AIAA Paper 2018-1041, 2018. https://doi.org/10.2514/6.2018-1041

[34] Coder, J. G., Pulliam, T. H., and Jensen, J. C., "Contributions to HiLiftPW-3 Using Structured, Overset Grid Methods," *2018 AIAA Aerospace Sciences Meeting*, AIAA Paper 2018-1039, 2018. https://doi.org/10.2514/6.2018-1039

[35] Coder, J. G., "Further Development of the Amplification Factor Transport Transition Model for Aerodynamic Flows," *2019 AIAA Aerospace Sciences Meeting*, AIAA Paper 2019-0039, 2019. https://doi.org/10.2514/6.2019-0039

[36] Spalart, P. R., and Allmaras, S. R., "A One Equation Turbulence Model for Aerodynamic Flows," *30th AIAA Aerospace Sciences Meeting and Exhibit*, AIAA Paper 1992-0439, 1992. https://doi.org/10.2514/6.1992-439

[37] Xu, J., Han, X., Qiao, L., Bai, J., and Zhang, Y., "Fully Local Amplification Factor Transport Equation for Stationary Crossflow Instabilities," *AIAA Journal*, Vol. 57, No. 7, 2019, pp. 2682–2693. https://doi.org/10.2514/1.J057502

[38] Choi, J. H., and Kwon, O. J., "Enhancement of a Correlation-Based Transition Turbulence Model for Simulating Crossflow Instability," *AIAA Journal*, Vol. 53, No. 10, 2015, pp. 3063–3072. https://doi.org/10.2514/1.J053887

[39] Grabe, C., and Krumbein, A., "Correlation-Based Transition Modeling for Three-Dimensional Aerodynamic Configurations," *Journal of Aircraft*, Vol. 50, No. 5, 2013, pp. 1533–1539. https://doi.org/10.2514/1.C032063

[40] Grabe, C., Shengyang, N., and Krumbein, A., "Transition Transport Modeling for the Prediction of Crossflow Transition," *34th AIAA Applied Aerodynamics Conference*, AIAA Paper 2016-3572, 2016. https://doi.org/10.2514/6.2016-3572

[41] Xu, J. K., Bai, J. Q., Qiao, L., and Zhang, Y., "Correlation-Based Transition Transport Modeling for Simulating Crossflow Instabilities," *Journal of Applied Fluid Mechanics*, Vol. 9, No. 5, 2016, pp. 2435–2442. https://doi.org/10.18869/acadpub.jafm.68.236.25356

[42] Muller, C., and Herbst, F., "Modelling of Crossflow-Induced Transition Based on Local Variables," *6th European Conference on Computational Fluid Dynamics (ECFD)*, Leibniz Univ. Hannover, Inst. of Turbomachinery and Fluid Dynamics, Germany, July 2014, pp. 6358–6369.

[43] Langtry, R. B., Sengupta, K., Yeh, D. T., and Dorgan, A. J., "Extending the $\gamma\text{-}\tilde{Re}_{\theta t}$ Correlation Based Transition Model for Crossflow Effects," *45th AIAA Fluid Dynamics Conference*, AIAA Paper 2015-2474, 2015. https://doi.org/10.2514/6.2015-2474

[44] Venkatachari, B. S., Paredes, P., Derlaga, J. M., Buning, P., Choudhari, M., Li, F., and Chang, C. L., "Assessment of Transition Modeling Capability in OVERFLOW with Emphasis on Swept-Wing Configurations," *2020 AIAA Aerospace Sciences Meeting*, AIAA Paper 2020-1034, 2020. https://doi.org/10.2514/6.2020-1034

[45] Menter, F. R., "Two-Equation Eddy-Viscosity Turbulence Models for Engineering Applications," *AIAA Journal*, Vol. 32, No. 8, 1994, pp. 1598–1605. https://doi.org/10.2514/3.12149

[46] Nie, S., Krimmelbein, N., Krumbein, A., and Grabe, C., "Extension of a Reynolds-Stress-Based Transition Transport Model for Crossflow Transition," *Journal of Aircraft*, Vol. 55, No. 4, 2018, pp. 1641–1654. https://doi.org/10.2514/1.C034586

[47] Medida, S., and Baeder, J. D., "Application of the Correlation-Based $\gamma\text{-}\tilde{Re}_{\theta t}$ Transition Model to the Spalart-Allmaras Turbulence Model," *20th AIAA Computational Fluid Dynamics Conference*, AIAA Paper 2011-3979, 2011. https://doi.org/10.2514/6.2011-3979

[48] Jung, Y. S., and Baeder, J. D., "$\gamma\text{-}\tilde{Re}_{\theta t}$-SA with Crossflow Transition Model Using Hamiltonian-Strand Approach," *2018 AIAA Aerospace Sciences Meeting*, AIAA Paper 2018-1040, 2018. https://doi.org/10.2514/6.2018-1040

[49] Carnes, J. A., and Coder, J. G., "Effect of Crossflow on the S-76 and PSP Rotors in Hover," *2020 AIAA Aerospace Sciences Meeting*, AIAA Paper 2020-0773, 2020. https://doi.org/10.2514/6.2020-0773

[50] Schucker, J., "Development of a Three-Equation $\gamma\text{-}\tilde{Re}_{\theta t}$-Spalart-Allmaras Turbulence-Transition Model," M.S. Thesis, DLR, German Aerospace Center, Göttingen, Germany, 2012.

[51] Langtry, R. B., "A Correlation-Based Transition Model Using Local Variables for Unstructured Parallelized CFD Codes," Ph.D. Thesis, DLR, German Aerospace Center, Göttingen, Germany, 2006.

[52] Pulliam, T. H., and Zingg, D. W., *Fundamental Algorithms in Computational Fluid Dynamics*, Springer International, Switzerland, 2014, Chap. 3.

[53] Suluksna, K., Dechaumphai, P., and Juntasaro, E., "Correlations for Modeling Transitional Boundary Layers Under Influences of Freestream Turbulence and Pressure Gradient," *International Journal of Heat and Fluid Flow*, Vol. 30, No. 1, 2009, pp. 66–75. https://doi.org/10.1016/j.ijheatfluidflow.2008.09.004

[54] Allmaras, S. R., and Johnson, F. T., "Modifications and Clarifications for the Implementation of the Spalart-Allmaras Turbulence Model," *Proceedings of the 7th International Conference on Computational Fluid Dynamics*, Paper ICCFD7-1902, Big Island, HI, July 2012, pp. 1–11.

[55] Osusky, M., "A Parallel Newton-Krylov-Schur Algorithm for the Reynolds-Averaged Navier-Stokes Equations," Ph.D. Thesis, Graduate Dept. of Aerospace Science and Engineering, Univ. of Toronto, Toronto, Ontario, Canada, 2013.

[56] Bertsekas, D. P., "Nondifferentiable Optimization via Approximation," *Mathematical Programming Study*, Vol. 3, Springer, Berlin, 1975, pp. 1–25. https://doi.org/10.1007/BFb0120696

[57] Chakraborty, A., Sinha Roy, A., and Dasgupta, B., "Non-Parametric Smoothing for Gradient Methods in Non-Differentiable Optimization Problems," *2016 IEEE International Conference on Systems, Man, and Cybernetics*, IEEE Publ., Piscataway, NJ, 2016, pp. 3759–3764. https://doi.org/10.1109/SMC.2016.7844819

[58] Rivest, R. L., "Game Tree Searching by Min/Max Approximation," *Artificial Intelligence*, Vol. 34, No. 1, 1987, pp. 77–96. https://doi.org/10.1016/0004-3702(87)90004-X

[59] Pee, E. Y., "On Solving Large-Scale Finite Minimax Problems Using Exponential Smoothing," *Journal of Optimization Theory and Applications*, Vol. 148, No. 2, 2011, pp. 390–421. https://doi.org/10.1007/s10957-010-9759-1

[60] Polak, E., Royset, J. O., and Womersley, R. S., "Algorithms with Adaptive Smoothing for Finite Minimax Problems," *Journal of Optimization Theory and Applications*, Vol. 119, No. 3, 2003, pp. 459–484. https://doi.org/10.1023/B:JOTA.0000006685.60019.3e

[61] Xingsi, L., "An Entropy-Based Aggregate Method for Minimax Optimization," *Engineering Optimization*, Vol. 18, No. 4, 1997, pp. 277–285. https://doi.org/10.1080/03052159208941026

[62] Rashad, R., "High-Fidelity Aerodynamic Shape Optimization for Natural Laminar Flow," Ph.D. Thesis, Univ. of Toronto, Toronto, Ontario, Canada, 2016.

[63] Richards, F. J., "A Flexible Growth Function for Empirical Use," *Journal of Experimental Botany*, Vol. 10, No. 2, 1959, pp. 290–301. https://doi.org/10.1093/jxb/10.2.290

[64] Hicken, J. E., and Zingg, D. W., "Parallel Newton-Krylov Solver for the Euler Equations Discretized Using Simultaneous Approximation Terms," *AIAA Journal*, Vol. 46, No. 11, 2008, pp. 2773–2786. https://doi.org/10.2514/1.34810

[65] Osusky, M., and Zingg, D. W., "Parallel Newton-Krylov-Schur Flow Solver for the Navier-Stokes Equations," *AIAA Journal*, Vol. 51, No. 12, 2013, pp. 2833–2851. https://doi.org/10.2514/1.J052487

[66] Saad, Y., and Schultz, M. H., "GMRES: A Generalized Minimal Residual Algorithm for Solving Nonsymmetric Linear Systems," *SIAM Journal on Scientific and Statistical Computing*, Vol. 7, No. 3, 1986, pp. 856–869. https://doi.org/10.1137/0907058

[67] Saad, Y., and Sosonkina, M., "Distributed Schur Complement Techniques for General Sparse Linear Systems," *SIAM Journal of Scientific Computing*, Vol. 21, No. 4, 1999, pp. 1337–1356. https://doi.org/10.1137/S1064827597328996

[68] Moryossef, Y., and Levy, Y., "Unconditionally Positive Implicit Procedure for Two-Equation Turbulence Models: Application to k-ω Turbulence Models," *Journal of Computational Physics*, Vol. 220, No. 1, 2006, pp. 88–108. https://doi.org/10.1016/j.jcp.2006.05.001

[69] Lian, C., Xia, G., and Merkle, C. L., "Impact of Source Terms on Reliability of CFD Algorithms," *Computers & Fluids*, Vol. 39, No. 10, 2010, pp. 1909–1922. https://doi.org/10.1016/j.compfluid.2010.06.021

[70] Mor-Yossef, Y., and Levy, Y., "The Unconditionally Positive-Convergent Implicit Time Integration Scheme for Two-Equation Turbulence Models: Revisited," *Computers and Fluids Journal*, Vol. 38, No. 10, 2009, pp. 1984–1994. https://doi.org/10.1016/j.compfluid.2009.06.005

[71] Mor-Yossef, Y., "Unconditionally Stable Time Marching Scheme for Reynolds Stress Models," *Journal of Computational Physics*, Vol. 276, Nov. 2014, pp. 635–664. https://doi.org/10.1016/j.jcp.2014.07.047

[72] Mosahebi, A., and Laurendeau, E., "Convergence Characteristics of Fully and Loosely Coupled Numerical Approaches for Transition Models," *AIAA Journal*, Vol. 53, No. 5, 2015, pp. 1399–1404. https://doi.org/10.2514/1.J053722

[73] Mosahebi, A., and Laurendeau, E., "Introduction of a Modified Segregated Numerical Approach for Efficient Simulation of $\gamma\text{-}\tilde{Re}_{\theta t}$ Transition Model," *International Journal of Computational Fluid Dynamics*, Vol. 29, Nos. 6–8, 2015, pp. 357–375. https://doi.org/10.1080/10618562.2015.1093624

[74] Candler, G. V., Subbareddy, P. K., and Nompelis, I., "Decoupled Implicit Method for Aerothermodynamics and Reacting Flows," *AIAA Journal*, Vol. 51, No. 5, 2013, pp. 1245–1254. https://doi.org/10.2514/1.J052070

[75] Chisholm, T. T., and Zingg, D. W., "A Jacobian-Free Newton-Krylov Algorithm for Compressible Turbulent Fluid Flows," *Journal of Computational Physics*, Vol. 228, No. 9, 2009, pp. 3490–3507. https://doi.org/10.1016/j.jcp.2009.02.004

[76] Yildirim, A., Kenway, G. K. W., Mader, C. A., and Martins, J. R. R. A., "A Jacobian-Free Approximate Newton-Krylov Startup Strategy for RANS Simulations," *Journal of Computational Physics*, Vol. 397, Nov. 2019, Paper 108741. https://doi.org/10.1016/j.jcp.2019.06.018

[77] Lian, C., Xia, G., and Merkle, C. L., "Solution-Limited Time Stepping to Enhance Reliability in CFD Applications," *Journal of Computational Physics*, Vol. 228, No. 13, 2009, pp. 4836–4857. https://doi.org/10.1016/j.jcp.2009.03.040

[78] Modisette, J. M., "An Automated Reliable Method for Two-Dimensional Reynolds-Averaged Navier-Stokes Simulations," Ph.D. Thesis, Massachusetts Inst. of Technology, Cambridge, MA, 2011.

[79] Swanson, R. C., and Turkel, E., "On Central-Difference and Upwind Schemes," *Journal of Computational Physics*, Vol. 101, No. 2, 1992, pp. 292–306. https://doi.org/10.1016/0021-9991(92)90007-L

[80] Somers, D. M., "Design and Experimental Results for a Natural-Laminar-Flow Airfoil for General Aviation Applications," NASA TP-1861, 1981.

[81] Coder, J. G., "Standard Test Cases for Transition Model Verification and Validation in Computational Fluid Dynamics," *2018 AIAA Aerospace Sciences Meeting*, AIAA Paper 2018-0029, 2018. https://doi.org/10.2514/6.2018-0029

[82] Denison, M., and Pulliam, T. H., "Implementation and Assessment of the Amplification Factor Transport Laminar-Turbulent Transition Model," *2018 AIAA Aerospace Sciences Meeting*, AIAA Paper 2018-3382, 2018. https://doi.org/10.2514/6.2018-3382

[83] Khayatzadeh, P., and Nadarajah, S. K., "Laminar-Turbulent Flow Simulation for Wind Turbine Profiles Using the $\gamma\text{-}\tilde{Re}_{\theta t}$ Transition Model," *Wind Energy*, Vol. 17, No. 6, 2014, pp. 901–918. https://doi.org/10.1002/we.1606

[84] Somers, D. M., "Design and Experimental Results for the S809 Airfoil," National Renewable Energy Lab. Rept. NRELSR-440-6918, Golden, CO, 1989. https://doi.org/10.2172/437668

[85] Dagenhart, J., and Saric, W., "Crossflow Stability and Transition Experiments in Swept-Wing Flow," NASA Langley Technical Report Server, NASA TP-1999-209344, 1999.

[86] Radeztsky, R. H., Reibert, M. S., and Saric, W. S., "Effect of Micron-Sized Roughness on Transition in Swept-Wing Flows," *31st Aerospace Sciences Meeting and Exhibit*, AIAA Paper 1993-0076, 1993. https://doi.org/10.2514/6.1993-76

[87] Jameson, A., Schmidt, W., and Turkel, E., "Numerical Solution of the Euler Equations by Finite Volume Methods Using Runge-Kutta Time-Stepping Schemes," *14th Fluid and Plasma Dynamics Conference*, AIAA Paper 1981-1259, 1981. https://doi.org/10.2514/6.1981-1259

[88] Pulliam, T. H., *Efficient Solution Methods for the Navier-Stokes Equations*, Lecture Series: Numerical Techniques for Viscous Flow Computation in Turbomachinery Bladings, von Karman Inst. for Fluid Dynamics, Belgium, 1986.

[89] Petzold, R., and Radespiel, R., "Transition on a Wing with Spanwise Varying Crossflow and Linear Stability Analysis," *AIAA Journal*, Vol. 53, No. 2, 2015, pp. 321–335. https://doi.org/10.2514/1.J053127

[90] Kruse, M., Munoz, F., and Radespiel, R., "Transition Prediction Results for Sickle Wing and NLF(1)-0416 Test Cases," *2018 AIAA Aerospace Sciences Meeting*, AIAA Paper 2018-0537, 2018. https://doi.org/10.2514/6.2018-0537

[91] Menter, F. R., Smirnov, P. E., Liu, T., and Avancha, R., "A One-Equation Local Correlation-Based Transition Model," *Flow, Turbulence and Combustion*, Vol. 95, No. 4, 2015, pp. 583–619. https://doi.org/10.1007/s10494-015-9622-4

*S. Fu, Associate Editor*
