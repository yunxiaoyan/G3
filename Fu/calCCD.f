      include 'sps.f'

      module mod_CCDsimple
      save
      integer :: nparticle
      integer :: nspshl,nspshu,nspspl,nspspu
      real*8,pointer :: tijab(:,:,:,:)
      real*8 :: v(1:nsps,1;nsps,1:nsps,1:nsps),f(1:nsps,1:nsps)
      real*8,pointer :: chi7(:,:,:,:)
      real*8,pointer :: chi8(:,:)
      real*8,pointer :: chi9(:,:)
      real*8,pointer :: chi10(:,:,:,:)


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

      end function


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

      end function


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

      end function


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

      end function


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

      end function


!      function CCDterm6(i,j,a,b)  ! Not useful for pairing problems


      function CCDterm7(i,j,a,b)
      implicit none
      integer :: i1,i2,i3,i4,i5,i6,i7,i8
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2,a3
      real*8 :: CCDterm7
      real*8 :: vA(1:nparticle*(nsps-nparticle))
      real*8 :: vB(1:nparticle*(nsps-nparticle))

      i1=nparticle*(nsps-nparticle)

      i8=0   !                            +1/2
      do i2=1,nspshu  ! l
        do i3=nspspl,nspspu  ! d
          i8=i8+1
          vA(i8)=chi7(i3,i2,a,i)
          vB(i8)=tijab(i2,j,i3,b)
        enddo
      enddo
      a1=dsdot(i1,vA,1,vB,1)
      CCDterm7=a1/2d0

      i8=0  !                            -1/2 Pab
      do i2=1,nspshu  ! l
        do i3=nspspl,nspspu  ! d
          i8=i8+1
          vA(i8)=chi7(i3,i2,b,i)
          vB(i8)=tijab(i2,j,i3,a)
        enddo
      enddo
      a1=dsdot(i1,vA,1,vB,1)
      CCDterm7=CCDterm7-a1/2d0

      i8=0  !                            -1/2 Pij
      do i2=1,nspshu  ! l
        do i3=nspspl,nspspu  ! d
          i8=i8+1
          vA(i8)=chi7(i3,i2,a,j)
          vB(i8)=tijab(i2,i,i3,b)
        enddo
      enddo
      a1=dsdot(i1,vA,1,vB,1)
      CCDterm7=CCDterm7-a1/2d0

      i8=0  !                             +1/2 Pab Pij
      do i2=1,nspshu  ! l
        do i3=nspspl,nspspu  ! d
          i8=i8+1
          vA(i8)=chi7(i3,i2,b,j)
          vB(i8)=tijab(i2,i,i3,a)
        enddo
      enddo
      a1=dsdot(i1,vA,1,vB,1)
      CCDterm7=CCDterm7+a1/2d0

      endfunction


      subroutine calchi7()
      implicit none
      integer :: i1,i2,i3,i4,i5,i6,i7,i8
      real*8 :: a1,a2,a3
      real*8 :: vA(1:nparticle*(nsps-nparticle))
      real*8 :: vB(1:nparticle*(nsps-nparticle))

      allocate( chi7(nspspl:nspspu,1:nspshu,nspspl:nspspu,1:nspshu) )
      i1=nparticle*(nsps-nparticle)

      do i2=nspspl,nspspu  ! d        
        do i3=1,nspshu  ! l
          do i4=nspspl,nspspu  ! a
            do i5=1,nspshu  ! i
              
              i8=0
              do i7=nspspl,nspspu  ! c
                do i6=1,nspshu  ! k
                  i8=i8+1
                  vA(i8)=v(i6,i3,i7,i2)
                  vB(i8)=tijab(i5,i6,i4,i7)
                enddo
              enddo
              a1=dsdot(i1,vA,1,vB,1)
              chi7(i2,i3,i4,i5)=a1

            enddo
          enddo
        enddo
      enddo
    
      end subroutine


      function CCDterm8(i,j,a,b)
      implicit none
      integer :: i1,i2,i3,i4,i5,i6,i7,i8
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2,a3
      real*8 :: CCDterm8
      real*8 :: vA(1:nparticle)
      real*8 :: vB(1:nparticle)

      i1=nparticle

      i8=0   !                            +1/2
      do i2=1,nspshu  ! l
        i8=i8+1
        vA(i8)=chi8(i2,i)
        vB(i8)=tijab(i2,j,a,b)
      enddo
      a1=dsdot(i1,vA,1,vB,1)
      CCDterm8=a1/2d0

      i8=0   !                            -1/2 Pij
      do i2=1,nspshu  ! l
        i8=i8+1
        vA(i8)=chi8(i2,j)
        vB(i8)=tijab(i2,i,a,b)
      enddo
      a1=dsdot(i1,vA,1,vB,1)
      CCDterm8=CCDterm8-a1/2d0
    
      endfunction


      subroutine calchi8()
      implicit none
      integer :: i1,i2,i3,i4,i5,i6,i7,i8
      real*8 :: a1,a2,a3
      real*8 :: vA(1:nparticle*(nsps-nparticle)*(nsps-nparticle))
      real*8 :: vB(1:nparticle*(nsps-nparticle)*(nsps-nparticle))

      allocate( chi8(1:nparticle,1:nparticle) )
      i1=nparticle*(nsps-nparticle)*(nsps-nparticle)

      do i3=1,nspshu  ! l
        do i5=1,nspshu  ! i
              
          i8=0
          do i2=nspspl,nspspu  ! d        
            do i7=nspspl,nspspu  ! c
              do i6=1,nspshu  ! k
                i8=i8+1
                vA(i8)=v(i6,i3,i7,i2)
                vB(i8)=tijab(i5,i6,i7,i2)
              enddo
            enddo
          enddo
          a1=dsdot(i1,vA,1,vB,1)
          chi8(i3,i5)=a1

        enddo
      enddo
    
      end subroutine


      function CCDterm9(i,j,a,b)
      implicit none
      integer :: i1,i2,i3,i4,i5,i6,i7,i8
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2,a3
      real*8 :: CCDterm9
      real*8 :: vA(1:(nsps-nparticle))
      real*8 :: vB(1:(nsps-nparticle))

      i1=nsps-nparticle

      i8=0   !                            +1/2
      do i2=nspspl,nspspu  ! d        
        i8=i8+1
        vA(i8)=chi9(i2,a)
        vB(i8)=tijab(i,j,i2,b)
      enddo
      a1=dsdot(i1,vA,1,vB,1)
      CCDterm9=a1/2d0

      i8=0   !                            -1/2 Pab
      do i2=nspspl,nspspu  ! d        
        i8=i8+1
        vA(i8)=chi9(i2,b)
        vB(i8)=tijab(i,j,i2,a)
      enddo
      a1=dsdot(i1,vA,1,vB,1)
      CCDterm9=CCDterm9-a1/2d0
    
      endfunction


      subroutine calchi9()
      implicit none
      integer :: i1,i2,i3,i4,i5,i6,i7,i8
      real*8 :: a1,a2,a3
      real*8 :: vA(1:nparticle*nparticle*(nsps-nparticle))
      real*8 :: vB(1:nparticle*nparticle*(nsps-nparticle))

      allocate( chi9(nspspl:nspspu,nspspl:nspspu) )
      i1=nparticle*nparticle*(nsps-nparticle)

      do i2=nspspl,nspspu  ! d        
        do i4=nspspl,nspspu  ! a
              
          i8=0
          do i7=nspspl,nspspu  ! c
            do i3=1,nspshu  ! l
              do i6=1,nspshu  ! k
                i8=i8+1
                vA(i8)=v(i6,i3,i7,i2)
                vB(i8)=tijab(i6,i3,i4,i7)
              enddo
            enddo
          enddo
          a1=dsdot(i1,vA,1,vB,1)
          chi9(i2,i4)=a1

        enddo
      enddo
    
      end subroutine


      function CCDterm10(i,j,a,b)
      implicit none
      integer :: i1,i2,i3,i4,i5,i6,i7,i8
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2,a3
      real*8 :: CCDterm10
      real*8 :: vA(1:nparticle*nparticle)
      real*8 :: vB(1:nparticle*nparticle)

      i1=nparticle*nparticle

      i8=0   !                            +1/4
      do i2=1,nspshu  ! l
        do i6=1,nspshu  ! k
          i8=i8+1
          vA(i8)=chi10(i6,i2,i,j)
          vB(i8)=tijab(i6,i2,a,b)
        enddo
      enddo
      a1=dsdot(i1,vA,1,vB,1)
      CCDterm10=a1/4d0

      endfunction


      subroutine calchi10()
      implicit none
      integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9
      real*8 :: a1,a2,a3
      real*8 :: vA(1:((nsps-nparticle)*(nsps-nparticle)))
      real*8 :: vB(1:((nsps-nparticle)*(nsps-nparticle)))

      allocate( chi10(1:nspshu,1:nspshu,1:nspshu,1:nspshu) )
      i1=nsps-nparticle
      i1=i1*i1

      do i6=1,nspshu  ! k
        do i3=1,nspshu  ! l
          do i5=1,nspshu  ! i
            do i9=1,nspshu  ! j
              
              i8=0
              do i2=nspspl,nspspu  ! d        
                do i7=nspspl,nspspu  ! c
                  i8=i8+1
                  vA(i8)=v(i6,i3,i7,i2)
                  vB(i8)=tijab(i5,i9,i7,i2)
                enddo
              enddo
              a1=dsdot(i1,vA,1,vB,1)
              chi10(i6,i3,i5,i9)=a1

            enddo
          enddo
        enddo
      enddo
    
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
