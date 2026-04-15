










module actuatorRegion

    use constants
    use communication, only: commType, internalCommType
    use actuatorRegionData
    implicit none

contains

    subroutine computeActuatorRegionVolume(nn, iRegion)
        use blockPointers, only: nDom, vol
        implicit none

        ! Inputs
        integer(kind=intType), intent(in) :: nn, iRegion

        ! Working
        integer(kind=intType) :: iii
        integer(kind=intType) :: i, j, k

        ! Loop over the region for this block
        do iii = actuatorRegions(iRegion)%blkPtr(nn - 1) + 1, actuatorRegions(iRegion)%blkPtr(nn)
            i = actuatorRegions(iRegion)%cellIDs(1, iii)
            j = actuatorRegions(iRegion)%cellIDs(2, iii)
            k = actuatorRegions(iRegion)%cellIDs(3, iii)

            ! Sum the volume of each cell within the region on this proc
            actuatorRegions(iRegion)%volLocal = actuatorRegions(iRegion)%volLocal + vol(i, j, k)
        end do

    end subroutine computeActuatorRegionVolume

    ! ----------------------------------------------------------------------
    !                                                                      |
    !                    No Tapenade Routine below this line               |
    !                                                                      |
    ! ----------------------------------------------------------------------

end module actuatorRegion
