










module communication
!
!       Contains the variable definition of the processor number,
!       myID and the number of processors, nProc, which belong to the
!       group defined by the communicator ADflow_comm_world. The range
!       of processor numbers is <0..Nproc-1>, i.e. the numbering
!       starts at 0. This is done for compatibility with MPI.
!       Furthermore this module contains the communication pattern for
!       all the multigrid levels.
!
    use constants, only: intType, realType
    implicit none
    save
!
!       The definition of the derived data type commListType, which
!       stores the i,j and k indices as well as the block id of the
!       data to be communicated. Send lists may contain interpolants
!       since the indices may refer to a stencil, while the receive
!       list does not. All interpolations should be done on the send
!       side to keep message sizes to a minimum.
!
end module communication
