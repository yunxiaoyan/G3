      include 'sps.f'

      module mod_CCDsimple
      save
      integer :: nparticle
      integer :: nspshl,nspshu,nspspl,nspshu
      real*8,pointer :: tijab(:,:,:,:)


      contains
      subroutine calCCDsimple
      use mod_sps
      implicit none
      integer :: i0,i1,i2
      real*8 :: a1,a2
      character(len=1) :: c1
      
      open(unit=101,file='parameter.dat')
      read(101,*) nparticle
      read(101,*) c1
      close(101)
      nspshl=1
      nspshu=nparticle
      nspspl=nparticle+1
      nspspu=nsps
      allocate( tijab(nspspl:nspspu,nspspl:nspspu,1:nspshu,1:nspshu) )
      tijab=0d0


      end subroutine


      end module mod_CCDsimple



      module mod_CCD
      save
      integer :: nparticle
      integer :: nspshl,nspshu,nspspl,nspshu
      type tCC
        integer :: i,j,a,b
      end type
      type(tCC),pointer :: ltijab(:)
      real*8,pointer :: tijab(:)
      integer :: ntijab

      end module mod_CCD
