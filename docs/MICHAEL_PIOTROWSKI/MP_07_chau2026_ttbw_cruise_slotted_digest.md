# MP_07 digest — Chau, Piotrowski & Duensing (2026, J. Aircraft), Cruise-slotted TTBW optimization

Companion to [`MP_07_chau2026_ttbw_cruise_slotted_full.md`](MP_07_chau2026_ttbw_cruise_slotted_full.md).

## Relevance to this branch: low. Read this section and stop unless you need ASO methodology.

**This paper is fully turbulent.** Piotrowski is a co-author but the transition
model is not used — transition appears only in the closing paragraph as future
work. There is no SA-sLM2015 content, no transition solver detail, no crossflow.
Nothing here bears on our convergence problem or on the transition model.

What it actually is: aerodynamic shape optimization of a Mach 0.800 single-aisle
**transonic truss-braced wing** with a **cruise-slotted** wing, versus a
non-slotted TTBW reference. Gradient-based drag minimization (Jetstream), with
post-optimization analysis in NASA's **LAVA Curvilinear** solver on structured
overset grids.

## The few things worth knowing

**Result.** The optimizer produces efficient designs for both configurations,
but the cruise-slotted TTBW carries **3.0 % more drag** at M 0.80 than the
non-slotted one (hence 3.0 % lower L/D at equal lift), driven by **higher
skin-friction drag across the span, "likely arising from the slotted region"**.
The slotted variant does show a **more gradual drag rise** with Mach number, so
the authors argue its case would be stronger under multipoint optimization,
where wave-drag trades return. The single-point optimizations were "unexpectedly"
successful at eliminating nearly all wave drag, which left the slotted wing's
advantage with nothing to buy.

**Why the skin-friction result is interesting to us anyway:** a configuration
whose penalty is almost entirely viscous is exactly the case where natural
laminar flow would change the answer. The authors say so — NLF variants of the
cruise-slotted wing are named as the follow-on, citing Somers, Coder & Somers,
and Hiller et al., and pointing at transition-enabled gradient-based ASO
(Husain et al. for transonic swept wings; Saadeh et al. for transonic
strut-braced wings). This is the paper that motivates the *next* one, not one we
can act on.

**Two-solver cross-check as practice.** Optimization in Jetstream on structured
multiblock grids; performance verification in LAVA on structured overset grids,
with grid-convergence studies and Richardson extrapolation to `L∞`. Comparing
optimized pressure distributions between two independent solvers before quoting
performance numbers is a good habit and cheap to imitate if we ever publish AR5
numbers.

**Geometry-control limitation, stated plainly.** Their FFD/axial scheme treats
the cruise-slotted wing as an integrated system and "is limited in its ability
to directly control the shape of the slot and aft element"; they suggest more
explicit schemes (Hiller et al.) that manipulate the aft element's shape,
location, and orientation directly. Generic ASO lesson: the parameterization
can silently cap the achievable benefit.

## Bottom line

Filed for completeness of the PDF set. If you are working on the transition
model or its solver, skip this one entirely — MP_03, MP_05 and MP_06 are where
the content is.
