      include 'sps.f'

      module mod_CCDsimplepairing
      save
      real*8 :: deltaE,weight,Delta
      integer :: nparticle
      integer :: nspshl,nspshu,nspspl,nspspu
      integer :: niter
      real*8,pointer :: tijab(:,:,:,:)
      real*8,pointer :: v(:,:,:,:),f(:,:)
      real*8,pointer :: chi7(:,:,:,:)
      real*8,pointer :: chi8(:,:)
      real*8,pointer :: chi9(:,:)
      real*8,pointer :: chi10(:,:,:,:)
      real*8 :: Ecorr,Embpt2,Ehf


      contains
      function CCDEc()
      use mod_sps
      implicit none
      integer :: i,j,k,l,a,b,c,d
      integer :: i1,i2
      real*8 :: a1,a2,a3
      real*8 :: CCDEc,ddot
      real*8 :: vA(1:nparticle*nparticle*(nsps-nparticle)
     $*(nsps-nparticle))
      real*8 :: vB(1:nparticle*nparticle*(nsps-nparticle)
     $*(nsps-nparticle))

      i1=nparticle*(nsps-nparticle)
      i1=i1*i1

      i2=0
      do b=nspspl,nspspu  
      do a=nspspl,nspspu  
      do j=1,nspshu   
      do i=1,nspshu   
        i2=i2+1
        vA(i2)=v(i,j,a,b)
        vB(i2)=tijab(i,j,a,b)
      enddo
      enddo
      enddo
      enddo
      a1=ddot(i1,vA,1,vB,1)

      CCDEc=a1/4d0

      end function


      function CCDterm1(i,j,a,b)
      implicit none
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2,a3
      real*8 :: CCDterm1

      CCDterm1=v(a,b,i,j)

      end function


      function CCDterm2(i,j,a,b)
      use mod_sps
      implicit none
      integer :: i1
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2,a3
      real*8 :: CCDterm2,ddot
      real*8 :: vA(nspspl:nspspu)
      real*8 :: vB(nspspl:nspspu)

      i1=nsps-nparticle

      vA(nspspl:nspspu)=f(nspspl:nspspu,b)
      vB(nspspl:nspspu)=tijab(i,j,a,nspspl:nspspu)
      a1=ddot(i1,vA,1,vB,1)
      CCDterm2=a1

      vA(nspspl:nspspu)=f(nspspl:nspspu,a)
      vB(nspspl:nspspu)=tijab(i,j,b,nspspl:nspspu)
      a1=ddot(i1,vA,1,vB,1)
      CCDterm2=CCDterm2-a1

      end function


      function CCDterm3(i,j,a,b)
      implicit none
      integer :: i1
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2,a3
      real*8 :: CCDterm3,ddot
      real*8 :: vA(nspshl:nspshu)
      real*8 :: vB(nspshl:nspshu)

      i1=nparticle

      vA(nspshl:nspshu)=f(j,nspshl:nspshu)
      vB(nspshl:nspshu)=tijab(i,nspshl:nspshu,a,b)
      a1=ddot(i1,vA,1,vB,1)
      CCDterm3=-a1

      vA(nspshl:nspshu)=f(i,nspshl:nspshu)
      vB(nspshl:nspshu)=tijab(j,nspshl:nspshu,a,b)
      a1=ddot(i1,vA,1,vB,1)
      CCDterm3=CCDterm3+a1

      end function


      function CCDterm4(i,j,a,b)
      use mod_sps
      implicit none
      integer :: i1,i2,i3,i4
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2,a3
      real*8 :: CCDterm4,ddot
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
      a1=ddot(i1,vA,1,vB,1)
      CCDterm4=a1/2d0

      end function


      function CCDterm5(i,j,a,b)
      implicit none
      integer :: i1,i2,i3,i4
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2,a3
      real*8 :: CCDterm5,ddot
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
      a1=ddot(i1,vA,1,vB,1)
      CCDterm5=a1/2d0

      end function


!      function CCDterm6(i,j,a,b)  ! Not useful for pairing problems


      function CCDterm7(i,j,a,b)
      use mod_sps
      implicit none
      integer :: i1,i2,i3,i4,i5,i6,i7,i8
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2,a3
      real*8 :: CCDterm7,ddot
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
      a1=ddot(i1,vA,1,vB,1)
      CCDterm7=a1/2d0
!      write(*,*) chi7(5,1,5,1)
!      write(*,*) i,j,a,b
!      write(*,*) CCDterm7

      i8=0  !                            -1/2 Pab
      do i2=1,nspshu  ! l
        do i3=nspspl,nspspu  ! d
          i8=i8+1
          vA(i8)=chi7(i3,i2,b,i)
          vB(i8)=tijab(i2,j,i3,a)
        enddo
      enddo
      a1=ddot(i1,vA,1,vB,1)
      CCDterm7=CCDterm7-a1/2d0
!      write(*,*) CCDterm7

      i8=0  !                            -1/2 Pij
      do i2=1,nspshu  ! l
        do i3=nspspl,nspspu  ! d
          i8=i8+1
          vA(i8)=chi7(i3,i2,a,j)
          vB(i8)=tijab(i2,i,i3,b)
        enddo
      enddo
      a1=ddot(i1,vA,1,vB,1)
      CCDterm7=CCDterm7-a1/2d0
!      write(*,*) CCDterm7

      i8=0  !                             +1/2 Pab Pij
      do i2=1,nspshu  ! l
        do i3=nspspl,nspspu  ! d
          i8=i8+1
          vA(i8)=chi7(i3,i2,b,j)
          vB(i8)=tijab(i2,i,i3,a)
        enddo
      enddo
      a1=ddot(i1,vA,1,vB,1)
      CCDterm7=CCDterm7+a1/2d0
!      write(*,*) CCDterm7
!      stop


      end function


      subroutine calchi7()
      use mod_sps
      implicit none
      integer :: i1,i2,i3,i4,i5,i6,i7,i8
      real*8 :: a1,a2,a3
      real*8 :: vA(1:nparticle*(nsps-nparticle))
      real*8 :: vB(1:nparticle*(nsps-nparticle))
      real*8 :: ddot

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
              a1=ddot(i1,vA,1,vB,1)
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
      real*8 :: CCDterm8,ddot
      real*8 :: vA(1:nparticle)
      real*8 :: vB(1:nparticle)

      i1=nparticle

      i8=0   !                            +1/2
      do i2=1,nspshu  ! l
        i8=i8+1
        vA(i8)=chi8(i2,i)
        vB(i8)=tijab(i2,j,a,b)
      enddo
      a1=ddot(i1,vA,1,vB,1)
      CCDterm8=a1/2d0

      i8=0   !                            -1/2 Pij
      do i2=1,nspshu  ! l
        i8=i8+1
        vA(i8)=chi8(i2,j)
        vB(i8)=tijab(i2,i,a,b)
      enddo
      a1=ddot(i1,vA,1,vB,1)
      CCDterm8=CCDterm8-a1/2d0
    
      end function


      subroutine calchi8()
      use mod_sps
      implicit none
      integer :: i1,i2,i3,i4,i5,i6,i7,i8
      real*8 :: a1,a2,a3
      real*8 :: vA(1:nparticle*(nsps-nparticle)*(nsps-nparticle))
      real*8 :: vB(1:nparticle*(nsps-nparticle)*(nsps-nparticle))
      real*8 :: ddot

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
          a1=ddot(i1,vA,1,vB,1)
          chi8(i3,i5)=a1

        enddo
      enddo
    
      end subroutine


      function CCDterm9(i,j,a,b)
      use mod_sps
      implicit none
      integer :: i1,i2,i3,i4,i5,i6,i7,i8
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2,a3
      real*8 :: CCDterm9,ddot
      real*8 :: vA(1:(nsps-nparticle))
      real*8 :: vB(1:(nsps-nparticle))

      i1=nsps-nparticle

      i8=0   !                            +1/2
      do i2=nspspl,nspspu  ! d        
        i8=i8+1
        vA(i8)=chi9(i2,a)
        vB(i8)=tijab(i,j,i2,b)
      enddo
      a1=ddot(i1,vA,1,vB,1)
      CCDterm9=a1/2d0

      i8=0   !                            -1/2 Pab
      do i2=nspspl,nspspu  ! d        
        i8=i8+1
        vA(i8)=chi9(i2,b)
        vB(i8)=tijab(i,j,i2,a)
      enddo
      a1=ddot(i1,vA,1,vB,1)
      CCDterm9=CCDterm9-a1/2d0
    
      end function


      subroutine calchi9()
      use mod_sps
      implicit none
      integer :: i1,i2,i3,i4,i5,i6,i7,i8
      real*8 :: a1,a2,a3
      real*8 :: vA(1:nparticle*nparticle*(nsps-nparticle))
      real*8 :: vB(1:nparticle*nparticle*(nsps-nparticle))
      real*8 :: ddot

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
          a1=ddot(i1,vA,1,vB,1)
          chi9(i2,i4)=a1

        enddo
      enddo
    
      end subroutine


      function CCDterm10(i,j,a,b)
      implicit none
      integer :: i1,i2,i3,i4,i5,i6,i7,i8
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2,a3
      real*8 :: CCDterm10,ddot
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
      a1=ddot(i1,vA,1,vB,1)
      CCDterm10=a1/4d0

      end function


      subroutine calchi10()
      use mod_sps
      implicit none
      integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9
      real*8 :: a1,a2,a3
      real*8 :: vA(1:((nsps-nparticle)*(nsps-nparticle)))
      real*8 :: vB(1:((nsps-nparticle)*(nsps-nparticle)))
      real*8 :: ddot

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
              a1=ddot(i1,vA,1,vB,1)
              chi10(i6,i3,i5,i9)=a1

            enddo
          enddo
        enddo
      enddo
    
      end subroutine


      subroutine calCCDsimplepairing()
      use mod_sps
      implicit none
      integer :: i0,i1,i2,i3,i4,i5
      integer :: i,j,k,l,a,b,c,d
      real*8 :: a1,a2
      character(len=1) :: c1
      integer :: ntbme
      real*8 :: aspe(1:nsps)
      real*8,pointer :: tijabnew(:,:,:,:)
      real*8 :: Ecorrnew
      
      write(*,*) 'This CCDsimple version works only for pairing 
     $problems. '

      open(unit=101,file='parameter.dat')
      read(101,*) nparticle
      read(101,*) c1
      close(101)
      nspshl=1
      nspshu=nparticle
      nspspl=nparticle+1
      nspspu=nsps
      allocate( tijab(1:nspshu,1:nspshu,nspspl:nspspu,nspspl:nspspu) )
      tijab=0d0
      i1=nparticle*(nsps-nparticle)
      if (i1>46340) then
        write(*,*) 'Particle and orbit number too large!!! STOP! '
        stop
      endif

      allocate( v(1:nsps,1:nsps,1:nsps,1:nsps), f(1:nsps,1:nsps) )
      v=0d0
      f=0d0
      open(unit=101,file='TBME.b')
      read(101,*)
      read(101,*)
      read(101,*)
      read(101,*) ntbme,aspe(:)
      read(101,*)
      do i1=1,ntbme
        read(101,*) i2,i3,i4,i5,a1
        v(i2,i3,i4,i5)=a1
        v(i3,i2,i4,i5)=-a1
        v(i2,i3,i5,i4)=-a1
        v(i3,i2,i5,i4)=a1
        v(i4,i5,i2,i3)=a1
        v(i5,i4,i2,i3)=-a1
        v(i4,i5,i3,i2)=-a1
        v(i5,i4,i3,i2)=a1
      enddo
      close(101)
      do i1=1,nsps
        do i2=1,nsps
          if (i1==i2) f(i2,i1)=f(i2,i1)+aspe(i2)
          do i3=1,nspshu
            f(i2,i1)=f(i2,i1)+v(i1,i3,i2,i3)
          enddo
        enddo
!        write(*,*) f(i1,i1)
      enddo

      a1=0d0
      do i1=1,nspshu
        do i2=1,nspshu
          a1=a1+v(i2,i1,i2,i1)
        enddo
      enddo
      a1=a1/2d0
      do i1=1,nspshu
        a1=a1+aspe(i1)
      enddo
      Ehf=a1
      write(*,'(a7,f20.10)')  '  Ehf =',Ehf

      !--------------- Initializing and calculating tijab --------------
      write(*,*) '  Calculating tijab... '
      do b=nspspl,nspspu  
      do a=nspspl,b-1  
      do j=1,nspshu   
      do i=1,j-1   
        a1=v(a,b,i,j)/(f(i,i)+f(j,j)-f(a,a)-f(b,b))
        tijab(i,j,a,b)=a1
        tijab(j,i,a,b)=-a1
        tijab(i,j,b,a)=-a1
        tijab(j,i,b,a)=a1
      enddo
      enddo
      enddo
      enddo
      Ecorr=CCDEc()
      Embpt2=Ecorr
      write(*,'(a14,f20.10)')  '  Init Ecorr =',Embpt2
      niter=0

      allocate(tijabnew(1:nspshu,1:nspshu,nspspl:nspspu,nspspl:nspspu))
      tijabnew=0d0

201   call calchi7()
      call calchi8()
      call calchi9()
      call calchi10()
      do b=nspspl,nspspu  
      do a=nspspl,b-1  
      do j=1,nspshu   
      do i=1,j-1   
        a1=CCDterm1(i,j,a,b)
     $+CCDterm2(i,j,a,b)
     $+CCDterm3(i,j,a,b)
     $+CCDterm4(i,j,a,b)
     $+CCDterm5(i,j,a,b)
     $+CCDterm7(i,j,a,b)  ! No term 6 for pairing problems
     $+CCDterm8(i,j,a,b)
     $+CCDterm9(i,j,a,b)
     $+CCDterm10(i,j,a,b)
        a1=tijab(i,j,a,b)+a1/(f(i,i)+f(j,j)-f(a,a)-f(b,b)-Delta)
        tijabnew(i,j,a,b)=a1
        tijabnew(j,i,a,b)=-a1
        tijabnew(i,j,b,a)=-a1
        tijabnew(j,i,b,a)=a1
      enddo
      enddo
      enddo
      enddo
      tijab=tijabnew*weight+(1d0-weight)*tijab
      Ecorrnew=CCDEc()
      niter=niter+1
      if (mod(niter,10)==0) then
        write(*,'(a12,i8)')  '  Iteration:',niter
        write(*,'(a11,f20.10)')  '    Ecorr =',Ecorrnew
      endif

      if (abs(Ecorrnew-Ecorr)>deltaE) then
        Ecorr=Ecorrnew
        go to 201
      else
        Ecorr=Ecorrnew
      endif

      deallocate( chi7,chi8,chi9,chi10,tijabnew )
      
      end subroutine


      end module mod_CCDsimplepairing



!      module mod_CCD
!      save
!      integer :: nparticle
!      integer :: nspshl,nspshu,nspspl,nspspu
!      type tCC
!        integer :: i,j,a,b
!      end type
!      type(tCC),pointer :: ltijab(:)
!      real*8,pointer :: tijab(:)
!      integer :: ntijab
!
!      end module mod_CCD
