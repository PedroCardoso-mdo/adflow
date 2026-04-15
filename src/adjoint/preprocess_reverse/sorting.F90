










module sorting

    use utils, only: terminate
contains

    function famInList(famID, famList)
        use constants
        implicit none
        integer(kind=intType), intent(in) :: famID, famList(:)
        logical :: famInList
        famInLIst = .False.
        if (bsearchIntegers(famID, famList) > 0) then
            famInList = .True.
        end if
    end function famInList

    function bsearchIntegers(key, base)
        !
        !       bsearchIntegers returns the index in base where key is stored.
        !       A binary search algorithm is used here, so it is assumed that
        !       base is sorted in increasing order. In case key appears more
        !       than once in base, the result is arbitrary. If key is not
        !       found, a zero is returned.
        !
        use precision
        implicit none
        !
        !      Function type
        !
        integer(kind=intType) :: bsearchIntegers
        !
        !      Function arguments.
        !
        integer(kind=intType), intent(in) :: key
        integer(kind=intType), dimension(:), intent(in) :: base
        integer(kind=intType) :: nn
        !
        !      Local variables.
        !
        integer(kind=intType) :: ii, pos, start
        logical :: entryFound

        ! Initialize some values.

        start = 1
        ii = size(base)
        entryFound = .false.

        ! Binary search to find key.

        do
            ! Condition for breaking the loop

            if (ii == 0) exit

            ! Determine the position in the array to compare.

            pos = start + ii / 2

            ! In case this is the entry, break the search loop.

            if (base(pos) == key) then
                entryFound = .true.
                exit
            end if

            ! In case the search key is larger than the current position,
            ! only parts to the right must be searched. Remember that base
            ! is sorted in increasing order. Nothing needs to be done if the
            ! key is smaller than the current element.

            if (key > base(pos)) then
                start = pos + 1
                ii = ii - 1
            end if

            ! Modify ii for the next branch to search.

            ii = ii / 2
        end do

        ! Set bsearchIntegers. This depends whether the key was found.

        if (entryFound) then
            bsearchIntegers = pos
        else
            bsearchIntegers = 0
        end if

    end function bsearchIntegers

    ! ----------------------------------------------------------------------
    !                                                                      |
    !                    No Tapenade Routine below this line               |
    !                                                                      |
    ! ----------------------------------------------------------------------

end module sorting
