










module surfaceFamilies

    use constants

    ! Special BC array's that are sometime required for reducitons.
    real(kind=realType), dimension(:, :), allocatable, target :: zeroCellVal
    real(kind=realType), dimension(:, :), allocatable, target :: oneCellVal
    real(kind=realType), dimension(:, :), allocatable, target :: zeroNodeVal

end module surfaceFamilies
