      module mod_sps
      save
      integer :: nsps  ! number of sp orbit
      type sps
        integer :: qnw  ! quantum number w
        integer :: qnn  ! quantum number n
        integer :: qnl  ! quantum number l
        integer :: qnj  ! quantum number 2j
        integer :: qnm  ! quantum number 2m
        integer :: qntau  ! quantum number 2tau
      end type
      type(sps),pointer :: sporbit(:)
      integer*8,pointer :: sporbitb(:)


      contains
      subroutine readsps_pairing()
      implicit none
      integer :: i1
      open(unit=101,file='sps.dat')
      read(101,*) nsps
      allocate(sporbit(1:nsps),sporbitb(1:nsps))
      sporbit(:)%qnw=-999999
      sporbit(:)%qnn=-999999
      sporbit(:)%qnl=-999999
      sporbit(:)%qnj=-999999
      sporbit(:)%qnm=-999999
      sporbit(:)%qntau=-999999
      sporbitb=0
      do i1=1,nsps
        read(101,*) sporbit(i1)%qnn,sporbit(i1)%qnm
        sporbitb(i1)=ibset(sporbitb(i1),i1)
      enddo
      close(101)
      end subroutine


      end module mod_sps
