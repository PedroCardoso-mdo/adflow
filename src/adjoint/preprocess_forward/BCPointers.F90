











module BCPointers

! Thiss module contains data structures used to apply BCs.

    use constants, only: intType, realType
    implicit none
    save

    real(kind=realType), dimension(:, :, :), pointer :: ww0, ww1, ww2, ww3
    real(kind=realType), dimension(:, :), pointer :: pp0, pp1, pp2, pp3
    real(kind=realType), dimension(:, :), pointer :: rlv0, rlv1, rlv2, rlv3
    real(kind=realType), dimension(:, :), pointer :: rev0, rev1, rev2, rev3
    real(kind=realType), dimension(:, :), pointer :: gamma0, gamma1, gamma2, gamma3
    real(kind=realType), dimension(:, :, :), pointer :: ssi, ssj, ssk
    real(kind=realType), dimension(:, :, :), pointer :: ss, xx
    real(kind=realType), dimension(:, :), pointer :: dd2wall, sFace
    integer(kind=intType), dimension(:, :), pointer :: gcp

    integer(kind=intType) :: iStart, iEnd, iSize
    integer(kind=intType) :: jStart, jEnd, jSize


end module BCPointers

