# Debugging Derivatives

Part of the high cost of implementing derivatives is verification and debugging of derivatives. In fact, as a student in the MDO Lab, it may feel like your job boils down to this very task. We will look at debugging the derivatives of an implicit set of equations since this is often the more difficult case and explicit equations can always be recast as a set of implicit ones (see UDE derivations).

In the end what we need are accurate total derivatives. Unfortunately there are many pieces that are used to compute the total derivative. To debug issues with total derivatives it is necessary to verify that each set of derivatives used to compute the total is correct. We can do this by climbing the derivative checking ladder.

```
+-----------------------------------------------------------------+
|                 check serial and parallel match                 |
+-----------------------------------------------------------------+
|      CS check of totals derivatives (adjoint and/or direct)     |  <--- very expensive
+-----------------------------------------------------------------+
|    dot product check of partials derivatives in reverse mode    |
+-----------------------------------------------------------------+
|    FD & CS check of partial derivatives in forward mode         |
+-----------------------------------------------------------------+
```

## FD & CS check of partial derivatives in forward mode

The first step is to verify the forward mode of the partial derivatives. This should be done with both complex step and finite difference (ideally centered finite difference).

If the finite difference matches the complex step to the expected value of precision then the derivatives are smooth in that region and the program was correctly complexified. If the complex step and the FWD routines match then the program was complexified correctly and FWD routines are correct.

If the finite difference and the derivative routines match, but they don't match CS then the code may be incorrectly complexified. You can then use the FWD routines to debug the CS. Check that intrinsic functions, such as `min`, `max`, `norm`, are behaving as you would expect.

If nothing matches, you may have one or more errors. Try a different step size for the complex step. See if the FD approximation is sensitive to the form (central, backwards, forward, etc).

## Dot product check of partials derivatives in reverse mode

To check that the reverse mode routines (REV) are consistent with the FWD mode, one can use the dot product test.

> **Note:** The dot product test only checks that the FWD and REV routines are consistent, which is insufficient to prove that either is correct. You need to check that at least one is correct (see the first step).

The mathematical definition of the dot product rule is given below. Let $\dot{x}$ and $\bar{y}$ be the input seeds to the FWD and REV routines respectively, and let $\dot{y}$ and $\bar{x}$ be the output seeds of the FWD and REV routines respectively:

$$\dot{y}^T \bar{y} = \bar{x}^T \dot{x}$$

*(Note: the exact symbols for the input/output seed vectors were not legible in the source document — the relation above is the standard dot-product/adjoint consistency identity: the dot product of the FWD-routine output with an arbitrary REV-routine input seed equals the dot product of the REV-routine output with the corresponding FWD-routine input seed.)*

The dot product test results from the combination of two ideas. The first is that derivative propagation is linear and thus can be reformulated as a matrix-vector product. The second is that vector-matrix-vector products are associative (you can do the multiplication in any order).

If the FWD and the REV routines are consistent then the dot product of the input of the FWD routine and the output of the REV routine should match the dot product of the output of the REV routine and the input of the FWD routine for a random set of inputs.

If the dot product test is failing, bisect the program and perform the dot product check with intermediate variables to find the place where the derivatives become inconsistent in your program.

> **Note:** When bisecting the program, the outputs of the FWD routine should match the variables used as inputs of the REV routine and vice-versa. Thus when inspecting intermediate variables you will need to look in opposite places (e.g. variables at the beginning of the FWD routine will have corresponding variables at the end of the REV routine).

> **Note:** Check that both the FWD and REV modes are still accurate after repeated calls, if you plan on solving for the total derivatives using a matrix-free approach.

## CS check of totals derivatives (adjoint and/or direct)

The total derivatives are formed by combining the partial derivatives. There are a number of things that can go wrong in the process of "combining".

The first thing to do is to verify that the linear system can be solved accurately for the adjoint and direct vector. If you cannot solve the linear system:

- verify that the matrix-vector products are the same inside and outside of the solver.
- check the assembly of the `pRpu` matrix.

If the linear system can be solved in one of the two modes (adjoint or direct) then the `pRpu` matrix is likely fine in the case of a matrix-based approach.

- check that the RHS is set correctly (be aware that a solver may be overwriting the value of the RHS during the solve and a new RHS copy may be required on each solve).

If the linear system is solved but the total derivatives are wrong:

- check that the sign of the adjoint/direct vector is correct when combining with the other derivatives.
- check that the adjoint/direct vector is being stored and passed to the FWD or REV routine properly.

## Check serial and parallel match

Solving for the derivatives can be just as expensive as solving for the primal flow. Thus it is important that they are also scalable.

Check that when run in parallel the total derivatives are still the same. If they are not, then you can use the serial routines to debug the parallel routines.

## Coupled Derivatives

At its heart, checking coupled derivatives is no different than checking derivatives for single disciplines. However before checking the coupled derivatives one must check the derivatives of each discipline independently first. Afterwards the "ladder" is much like that of a single discipline.

```
+-----------------------------------------------------------------+
|                 check serial and parallel match                 |
+-----------------------------------------------------------------+
|              CS check of totals derivatives                     |  <--- very very expensive
+-----------------------------------------------------------------+
|   dot product test of REV partials derivative coupling terms    |
+-----------------------------------------------------------------+
|     FD & CS of FWD partial derivatives of coupling terms        |
+-----------------------------------------------------------------+
|          All derivative checks for each component               |
+-----------------------------------------------------------------+
```

## Things that get in the way of good debugging

### Tools aren't set up

It is common to rely too heavily on finite difference because complexified versions of the code have not been compiled. Although it can be a frustrating detour to complexify a code when you are in a rush to find a bug, it is worth it. Not only will complexified code make it easier in the future to debug, it will also give you more confidence in your derivatives.

Similarly, setting up and learning a debugger will make all future debugging easier, but has an initial cost. It is advised that you set up and use these tools before they are needed.

### Don't know the discipline

If you are unfamiliar with the discipline, you may miss some obvious signs pointing towards the error. Be sure to read the relevant literature and don't be afraid to ask for help from other members with more experience.

### No tests to work with

Tests help you narrow down what is going wrong by testing the range of capabilities of the program. Furthermore, the tests provide assurance that a feature was indeed working originally. If there is no existing test that covers the bug, you should add a test that does so you can use it in debugging and ensure that the bug will be caught immediately if it reappears again.
