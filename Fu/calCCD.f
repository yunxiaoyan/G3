      include 'sps.f'

      module mod_CCDsimple
      save
      integer :: nparticle
      integer :: nspshl,nspshu,nspspl,nspspu
      real*8,pointer :: tijab(:,:,:,:)
      real*8 :: v(1:nsps,1;nsps,1:nsps,1:nsps),f(1:nsps,1:nsps)


      contains
      subroutine calCCDsimple
      use mod_sps
      implicit none
      integer :: i0,i1,i2
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2
      character(len=1) :: c1
      integer :: ntbme
      real*8 :: aspe(1:nsps)
      real*8 :: term1,term2,term3,term4,term5,term6,term7,term8,term9
     $,term10
      
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

      v=0d0
      f=0d0
      open(unit=101,file='TBME.b')
      read(101,*)
      read(101,*)
      read(101,*)
      read(101,*) ntbme,aspe(:)
      read(101,*)
      do i1=1,ntbme
        read(101,*) i2,i3,i4,i5,v(i2,i3,i4,i5)
      enddo
      close(101)
      do i1=1,nsps
        do i2=1,nsps
          if (i1==i2) f(i2,i1)=f(i2,i1)+aspe(i2)
          do i3=1,nspshu
            f(i2,i1)=f(i2,i1)+v(i1,i3,i2,i3)
          enddo
        enddo
      enddo

      
      end subroutine


      function CCDterm1(i,j,a,b)
      implicit none
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2,a3
      real*8 :: CCDterm1

      CCDterm1=v(a,b,i,j)

      endfunction


      function CCDterm2(i,j,a,b)
      implicit none
      integer :: i1
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2,a3
      real*8 :: CCDterm2
      real*8 :: vA(nspspl:nspspu)
      real*8 :: vB(nspspl:nspspu)

      i1=nsps-nparticle

      vA(nspspl:nspspu)=f(nspspl:nspspu,b)
      vB(nspspl:nspspu)=tijab(i,j,a,nspspl:nspspu)
      a1=dsdot(i1,vA,1,vB,1)
      CCDterm2=a1

      vA(nspspl:nspspu)=f(nspspl:nspspu,a)
      vB(nspspl:nspspu)=tijab(i,j,b,nspspl:nspspu)
      a1=dsdot(i1,vA,1,vB,1)
      CCDterm2=CCDterm2-a1

      endfunction


      function CCDterm3(i,j,a,b)
      implicit none
      integer :: i1
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2,a3
      real*8 :: CCDterm3
      real*8 :: vA(nspshl:nspshu)
      real*8 :: vB(nspshl:nspshu)

      i1=nparticle

      vA(nspshl:nspshu)=f(j,nspshl:nspshu)
      vB(nspshl:nspshu)=tijab(i,nspshl:nspshu,a,b)
      a1=dsdot(i1,vA,1,vB,1)
      CCDterm3=-a1

      vA(nspshl:nspshu)=f(i,nspshl:nspshu)
      vB(nspshl:nspshu)=tijab(j,nspshl:nspshu,a,b)
      a1=dsdot(i1,vA,1,vB,1)
      CCDterm3=CCDterm3+a1

      endfunction


      function CCDterm4(i,j,a,b)
      implicit none
      integer :: i1,i2,i3,i4
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2,a3
      real*8 :: CCDterm4
      real*8 :: vA(1:((nsps-nparticle)*(nsps-nparticle)))
      real*8 :: vB(1:((nsps-nparticle)*(nsps-nparticle)))

      i1=nsps-nparticle
      i1=i1*i1

      i4=0
      do i2=nspspl,nspspu
        do i3=nspspl,nspspu
          i4=i4+1
          vA(i4)=v(a,b,i2,i3)
          vB(i4)=tijab(i,j,i2,i3)
        enddo
      enddo
      a1=dsdot(i1,vA,1,vB,1)
      CCDterm4=a1/2d0

      endfunction


      function CCDterm5(i,j,a,b)
      implicit none
      integer :: i1,i2,i3,i4
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2,a3
      real*8 :: CCDterm5
      real*8 :: vA(1:nparticle*nparticle)
      real*8 :: vB(1:nparticle*nparticle)

      i1=nparticle*nparticle

      i4=0
      do i2=1,nspshu
        do i3=1,nspshu
          i4=i4+1
          vA(i4)=v(i2,i3,i,j)
          vB(i4)=tijab(i2,i3,a,b)
        enddo
      enddo
      a1=dsdot(i1,vA,1,vB,1)
      CCDterm5=a1/2d0

      endfunction


!      function CCDterm6(i,j,a,b)  ! Not useful for pairing problems


      function CCDterm7(i,j,a,b)
      implicit none
      integer :: i1,i2,i3,i4
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2,a3
      real*8 :: CCDterm7
      real*8 :: ch1(nspspl:nspspu,1:nspshu,1:nspshu,nspspl:nspspu)
      real*8 :: vA(1:nparticle*(nsps-nparticle))
      real*8 :: vB(1:nparticle*(nsps-nparticle))

      i1=nparticle*(nsps-nparticle)

      i4=0
      do i2=1,nspshu
        do i3=nspspl,nspspu
          i4=i4+1
          vA(i4)=v(
          vB(i4)=tijab(
        enddo
      enddo
      a1=dsdot(i1,vA,1,vB,1)

      endfunction


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
