module ADjointPETSc

    !      This module contains the objects used by PETSc for the
    !      solution of the discrete adjoint equations.
    !
    use constants
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    Mat dRdWT, dRdWPreT

    ! These are empty vectors
    Vec w_like1, w_like2, psi_like1, psi_like2, psi_like3, x_like
    ! This logical is used to indicate whether the vectors have been created
    logical :: adjointPETScPreProcVarsAllocated
    PetscErrorCode PETScIerr
    PetscFortranAddr matfreectx(1)

    !adjointKSP   Linear solver (Krylov subspace method) context
    KSP adjointKSP

    ! Initial, start and final adjoint residuals
    real(kind=alwaysRealType) :: adjResInit
    real(kind=alwaysRealType) :: adjResStart
    real(kind=alwaysRealType) :: adjResFinal
    logical :: adjointPETScVarsAllocated

    ! Buffer of intermediate adjoint solution estimates, sampled every
    ! psiHistoryStep KSP iterations during solveAdjoint, used to report
    ! how the total derivative evolves as the adjoint converges. Only
    ! psiHistory(:, 1:psiHistoryCount) / psiHistoryIters(1:psiHistoryCount)
    ! hold valid data. Populated in MyKSPMonitor, (re)allocated/reset in
    ! solveAdjoint.
    real(kind=realType), dimension(:, :), allocatable :: psiHistory
    integer(kind=intType), dimension(:), allocatable :: psiHistoryIters
    real(kind=alwaysRealType), dimension(:), allocatable :: psiHistoryResid
    integer(kind=intType) :: psiHistoryCount

end module ADjointPETSc
