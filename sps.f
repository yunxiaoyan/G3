      module mod_sps
      save
      integer :: nsps  ! number of sp orbit
      type sps
        integer :: qnw  ! quantum number w
        integer :: qnn  ! quantum number n
        integer :: qnl  ! quantum number l
        integer :: qnj  ! quantum number 2j
        integer :: qnm  ! quantum number m or 2m
      end type
      type(sps),pointer :: sporbit(:)

      contains
      subroutine readsps_pairing()
      implicit none
      integer :: i1
      open(unit=101,file='sps.dat')
      read(101,*) nsps
      allocate(sporbit(1:nsps))
      sporbit(:)%qnw=-999999
      sporbit(:)%qnn=-999999
      sporbit(:)%qnl=-999999
      sporbit(:)%qnj=-999999
      sporbit(:)%qnm=-999999
      do i1=1,nsps
        read(101,*) sporbit(i1)%qnn,sporbit(i1)%qnm
      enddo
      close(101)
      end subroutine

      end module mod_sps



 
!      program one
!      implicit none
!      integer :: i1
!      integer :: a
!      integer*8 :: b,c
!
!      b=0
!      b=ibset(b,1)
!      b=ibset(b,2)
!      b=ibset(b,3)
!      b=ibset(b,64)
!      b=ishftc(b,1)
!      do i1=1,64
!        write(*,*) btest(b,i1)
!      enddo
!
!
!      end
