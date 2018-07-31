      include 'sps.f'

      module mod_CCDcommon
      save
      real*8 :: deltaE,weight,Delta
      real*8 :: Ecorr,Embpt2,Ehf
      integer :: niter
      end module mod_CCDcommon


      module mod_CCDsimplepairing
      use mod_CCDcommon
      save
      integer,private :: nparticle
      integer,private :: nspshl,nspshu,nspspl,nspspu
      real*8,pointer,private :: tijab(:,:,:,:)
      real*8,pointer,private :: v(:,:,:,:),f(:,:)
      real*8,pointer,private :: chi7(:,:,:,:)
      real*8,pointer,private :: chi8(:,:)
      real*8,pointer,private :: chi9(:,:)
      real*8,pointer,private :: chi10(:,:,:,:)


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
      a1=ddot(i1,va,1,vb,1)
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
      integer :: ntbme
      real*8 :: aspe(1:nsps)
      real*8,pointer :: tijabnew(:,:,:,:)
      real*8 :: Ecorrnew
      
      write(*,*) 'This CCDsimple version works only for pairing 
     $problems. '

      open(unit=101,file='parameter.dat')
      read(101,*) nparticle
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
      write(*,'(a7,f20.15)')  '  Ehf =',Ehf

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
      write(*,'(a14,f20.15)')  '  Init Ecorr =',Embpt2
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
!        a1=tijab(i,j,a,b)-a1/(f(i,i)+f(j,j)-f(a,a)-f(b,b)-Delta)**2
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
        write(*,'(a11,f25.15)')  '    Ecorr =',Ecorrnew
      endif

      if (abs(Ecorrnew-Ecorr)>deltaE) then
        Ecorr=Ecorrnew
        go to 201
      else
        Ecorr=Ecorrnew
      endif
!      write(*,*) weight,Delta,deltaE

      deallocate( chi7,chi8,chi9,chi10,tijabnew )
      
      end subroutine


      end module mod_CCDsimplepairing


      module mod_CCDnuclearmatter
      use mod_CCDcommon
      save
      integer,private :: nparticle
      integer,private :: nspshl,nspshu,nspspl,nspspu
      integer,private :: nhhpp
      integer,pointer,private :: rhhpp(:,:,:)
      integer,pointer,private :: btablehhpp(:,:)
      type tv
        integer,private :: nnhh
        integer,private :: nnpp
        integer,pointer,private :: rhh(:,:,:,:,:,:,:)
        integer,pointer,private :: rpp(:,:,:,:,:,:,:)
        integer,pointer,private :: tablehh(:,:)
        integer,pointer,private :: tablepp(:,:)
        real*8,pointer,private :: tpphh(:,:)
        real*8,pointer,private :: tpphhnew(:,:)
        real*8,pointer,private :: vpphh(:,:)
        real*8,pointer,private :: vpppp(:,:)
        real*8,pointer,private :: vhhhh(:,:)
      end type
      type(tv),pointer :: stablehhpp(:)
      real*8,pointer,private :: f(:,:)

!      real*8,pointer,private :: chi7(:,:,:,:)
!      real*8,pointer,private :: chi8(:,:)
!      real*8,pointer,private :: chi9(:,:)
!      real*8,pointer,private :: chi10(:,:,:,:)


      contains
      subroutine calCCDnuclearmatter()
      use mod_sps
      implicit none
      integer :: i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10
      integer :: i11,i12,i13,i14,i15,i16,i17,i18,i19
      integer,pointer :: itmp1(:,:,:),itmp2(:,:,:)
      integer :: itmp11(1:3),itmp12(1:7),itmp13(1:7)
      integer :: itmp21(1:5),itmp22(1:5),itmp23(1:5),itmp24(1:5)
      real*8 :: a1,a2
      integer :: i,j,k,l,a,b,c,d
      integer :: Nmax,ntbme
      integer :: bPx,bPy,bPz
      integer :: spx,spy,spz,s1,s2,t1,t2
      real*8 :: aspe(1:nsps)
      real*8 :: Ecorrnew
      
      write(*,*) 'This CCDsimple version works only for nuclear matter.'

      open(unit=101,file='parameter.dat')
      read(101,*)
      read(101,*)
      read(101,*)
      read(101,*) nparticle,i1,Nmax
      close(101)
      nspshl=1
      nspshu=nparticle
      nspspl=nparticle+1
      nspspu=nsps

      allocate(                                      !          building big table  
     $itmp1(-2*Nmax:2*Nmax,-2*Nmax:2*Nmax,-2*Nmax:2*Nmax) 
     $,itmp2(-2*Nmax:2*Nmax,-2*Nmax:2*Nmax,-2*Nmax:2*Nmax) 
     $,rhhpp(-2*Nmax:2*Nmax,-2*Nmax:2*Nmax,-2*Nmax:2*Nmax) 
     $)
      itmp1=0
      itmp2=0
      rhhpp=0
      do i1=1,nparticle
        do i2=1,nparticle
          bPx=sporbit(i2)%qnnm(1)+sporbit(i1)%qnnm(1)
          bPy=sporbit(i2)%qnnm(2)+sporbit(i1)%qnnm(2)
          bPz=sporbit(i2)%qnnm(3)+sporbit(i1)%qnnm(3)
          itmp1(bPx,bPy,bPz)=itmp1(bPx,bPy,bPz)+1
        enddo
      enddo
      do i1=nspspl,nspspu
        do i2=nspspl,nspspu
          bPx=sporbit(i2)%qnnm(1)+sporbit(i1)%qnnm(1)
          bPy=sporbit(i2)%qnnm(2)+sporbit(i1)%qnnm(2)
          bPz=sporbit(i2)%qnnm(3)+sporbit(i1)%qnnm(3)
          if (itmp1(bPx,bPy,bPz)/=0) then
            itmp2(bPx,bPy,bPz)=itmp2(bPx,bPy,bPz)+1
          endif
        enddo
      enddo

      i1=0
      do bPz=-2*Nmax,2*Nmax  ! Pz
        do bPy=-2*Nmax,2*Nmax  ! Py
          do bPx=-2*Nmax,2*Nmax  ! Px
            if ( (itmp1(bPx,bPy,bPz)/=0)
     $.and.(itmp2(bPx,bPy,bPz)/=0) ) then
              i1=i1+1
              rhhpp(bPx,bPy,bPz)=i1
            endif
          enddo
        enddo
      enddo
      nhhpp=i1

      allocate( btablehhpp(1:3,nhhpp),stablehhpp(1:nhhpp) )  ! building small table
      i1=0
      do bPz=-2*Nmax,2*Nmax  ! Pz
        do bPy=-2*Nmax,2*Nmax  ! Py
          do bPx=-2*Nmax,2*Nmax  ! Px
            if ( (itmp1(bPx,bPy,bPz)/=0)
     $.and.(itmp2(bPx,bPy,bPz)/=0) ) then
              i1=i1+1
              btablehhpp(1:3,i1)=[bPx,bPy,bPz]
              i2=itmp1(bPx,bPy,bPz)  ! nnhh
              i3=itmp2(bPx,bPy,bPz)  ! nnpp
              stablehhpp(i1)%nnhh=i2
              stablehhpp(i1)%nnpp=i3
              allocate( 
     $stablehhpp(i1)%rhh(-2*Nmax:2*Nmax,-2*Nmax:2*Nmax,-2*Nmax:2*Nmax
     $,-1:1,-1:1,-1:1,-1:1)
     $,stablehhpp(i1)%rpp(-2*Nmax:2*Nmax,-2*Nmax:2*Nmax,-2*Nmax:2*Nmax
     $,-1:1,-1:1,-1:1,-1:1)
     $,stablehhpp(i1)%tablehh(1:7,1:i2)
     $,stablehhpp(i1)%tablepp(1:7,1:i3)
     $,stablehhpp(i1)%tpphh(1:i3,1:i2)
     $,stablehhpp(i1)%tpphhnew(1:i3,1:i2)
     $,stablehhpp(i1)%vpphh(1:i3,1:i2)
     $,stablehhpp(i1)%vpppp(1:i3,1:i3)
     $,stablehhpp(i1)%vhhhh(1:i2,1:i2)
     $)
              stablehhpp(i1)%rhh=0
              stablehhpp(i1)%rpp=0
              stablehhpp(i1)%tpphh=0d0
              stablehhpp(i1)%tpphhnew=0d0
              stablehhpp(i1)%vpphh=0d0
              stablehhpp(i1)%vpppp=0d0
              stablehhpp(i1)%vhhhh=0d0
              i6=0
              do i4=1,nparticle
                do i5=1,nparticle
                  spx=sporbit(i4)%qnnm(1)+sporbit(i5)%qnnm(1)
                  spy=sporbit(i4)%qnnm(2)+sporbit(i5)%qnnm(2)
                  spz=sporbit(i4)%qnnm(3)+sporbit(i5)%qnnm(3)
                 if ( (spx/=bPx).or.(spy/=bPy).or.(spz/=bPz) ) go to 202
                  spx=2*sporbit(i4)%qnnm(1)-bPx
                  spy=2*sporbit(i4)%qnnm(2)-bPy
                  spz=2*sporbit(i4)%qnnm(3)-bPz
                  s1=sporbit(i4)%qnnm(4)
                  s2=sporbit(i5)%qnnm(4)
                  t1=sporbit(i4)%qnnm(5)
                  t2=sporbit(i5)%qnnm(5)
                  i6=i6+1
                  stablehhpp(i1)%tablehh(1:7,i6)
     $=[spx,spy,spz,s1,s2,t1,t2]
                  stablehhpp(i1)%rhh(spx,spy,spz,s1,s2,t1,t2)=i6
202             enddo
              enddo
              if (i6/=i2) then
                write(*,*) '  Building table error 1!!! STOP! '
                stop
              endif
              i6=0
              do i4=nspspl,nspspu
                do i5=nspspl,nspspu
                  spx=sporbit(i4)%qnnm(1)+sporbit(i5)%qnnm(1)
                  spy=sporbit(i4)%qnnm(2)+sporbit(i5)%qnnm(2)
                  spz=sporbit(i4)%qnnm(3)+sporbit(i5)%qnnm(3)
                 if ( (spx/=bPx).or.(spy/=bPy).or.(spz/=bPz) ) go to 203
                  spx=2*sporbit(i4)%qnnm(1)-bPx
                  spy=2*sporbit(i4)%qnnm(2)-bPy
                  spz=2*sporbit(i4)%qnnm(3)-bPz
                  s1=sporbit(i4)%qnnm(4)
                  s2=sporbit(i5)%qnnm(4)
                  t1=sporbit(i4)%qnnm(5)
                  t2=sporbit(i5)%qnnm(5)
                  i6=i6+1
                  stablehhpp(i1)%tablepp(1:7,i6)
     $=[spx,spy,spz,s1,s2,t1,t2]
                  stablehhpp(i1)%rpp(spx,spy,spz,s1,s2,t1,t2)=i6
203             enddo
              enddo
              if (i6/=i3) then
                write(*,*) '  Building table error 2!!! STOP! '
                stop
              endif
            endif
          enddo
        enddo
      enddo
      deallocate( itmp1,itmp2 )

      allocate( f(1:nsps,1:nsps) )
      f=0d0
      open(unit=101,file='fock.b')
      read(101,*) ntbme
      do i1=1,ntbme
        read(101,*) i2,i3,a1
        f(i2,i3)=a1
      enddo
      close(101)
      open(unit=101,file='TBME.b')
      read(101,*)
      read(101,*)
      read(101,*)
      read(101,*) ntbme
      read(101,*)
      do i0=1,ntbme
        read(101,*) i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,a1
        do i18=1,nhhpp
          i19=stablehhpp(i18)%rpp(i1,i2,i3,i4,i5,i6,i7)
          i15=stablehhpp(i18)%rpp(i8,i9,i10,i11,i12,i13,i14)
          i16=stablehhpp(i18)%rhh(i1,i2,i3,i4,i5,i6,i7)
          i17=stablehhpp(i18)%rhh(i8,i9,i10,i11,i12,i13,i14)
          if ( (i16/=0).and.(i17/=0) ) then
            stablehhpp(i18)%vhhhh(i16,i17)=a1
            stablehhpp(i18)%vhhhh(i17,i16)=a1
          endif
          if ( (i19/=0).and.(i17/=0) ) then
            stablehhpp(i18)%vpphh(i19,i17)=a1
          endif
          if ( (i15/=0).and.(i16/=0) ) then
            stablehhpp(i18)%vpphh(i15,i16)=a1
          endif
          if ( (i19/=0).and.(i15/=0) ) then
            stablehhpp(i18)%vpppp(i19,i15)=a1
            stablehhpp(i18)%vpppp(i15,i19)=a1
          endif
        enddo
      enddo
      close(101)
   
      !--------------- Initializing and calculating tijab --------------
      write(*,*) '  Calculating tijab... '
      do i1=1,nhhpp
        i2=stablehhpp(i1)%nnhh
        i3=stablehhpp(i1)%nnpp
        itmp11(1:3)=btablehhpp(1:3,i1)                ! bigP
        do i4=1,i2  ! ij
          itmp12(1:7)=stablehhpp(i1)%tablehh(1:7,i4)  ! smallp
          itmp21(1:3)=( itmp11(1:3)+itmp12(1:3) )/2
          itmp22(1:3)=( itmp11(1:3)-itmp12(1:3) )/2
          itmp21(4)=itmp12(4)
          itmp21(5)=itmp12(6)
          itmp22(4)=itmp12(5)
          itmp22(5)=itmp12(7)
          i6=rspsnm( itmp21(1),itmp21(2),itmp21(3),itmp21(4),itmp21(5) )
          i7=rspsnm( itmp22(1),itmp22(2),itmp22(3),itmp22(4),itmp22(5) )
          do i5=1,i3  ! ab
            itmp13(1:7)=stablehhpp(i1)%tablepp(1:7,i4)  ! smallp
            itmp23(1:3)=( itmp11(1:3)+itmp13(1:3) )/2
            itmp24(1:3)=( itmp11(1:3)-itmp13(1:3) )/2
            itmp23(4)=itmp13(4)
            itmp23(5)=itmp13(6)
            itmp24(4)=itmp13(5)
            itmp24(5)=itmp13(7)
            i8=rspsnm(itmp23(1),itmp23(2),itmp23(3),itmp23(4),itmp23(5))
            i9=rspsnm(itmp24(1),itmp24(2),itmp24(3),itmp24(4),itmp24(5))

            stablehhpp(i1)%tpphh(i5,i4)
     $=stablehhpp(i1)%vpphh(i5,i4)/(f(i6,i6)+f(i7,i7)-f(i8,i8)-f(i9,i9))

          enddo
        enddo
      enddo


!      do i1=1,nhhpp
!        write(*,*)  i1,btablehhpp(1:3,i1)
!      enddo
!      write(*,*)  stablehhpp(13)%tpphh
!      stop

      Ecorr=CCDEcnm()
      Embpt2=Ecorr
      write(*,'(a14,f20.10)')  '  Init Ecorr =',Embpt2
      niter=0

!      call calchi7()
!      call calchi8()
!      call calchi9()
!      call calchi10()
204   continue
      call CCDterm145()
      call CCDterm23()
      do i1=1,nhhpp
        i2=stablehhpp(i1)%nnhh
        i3=stablehhpp(i1)%nnpp
        itmp11(1:3)=btablehhpp(1:3,i1)                ! bigP
        do i4=1,i2  ! ij
          itmp12(1:7)=stablehhpp(i1)%tablehh(1:7,i4)  ! smallp
          itmp21(1:3)=( itmp11(1:3)+itmp12(1:3) )/2
          itmp22(1:3)=( itmp11(1:3)-itmp12(1:3) )/2
          itmp21(4)=itmp12(4)
          itmp21(5)=itmp12(6)
          itmp22(4)=itmp12(5)
          itmp22(5)=itmp12(7)
          i6=rspsnm( itmp21(1),itmp21(2),itmp21(3),itmp21(4),itmp21(5) )
          i7=rspsnm( itmp22(1),itmp22(2),itmp22(3),itmp22(4),itmp22(5) )
          do i5=1,i3  ! ab
            itmp13(1:7)=stablehhpp(i1)%tablepp(1:7,i4)  ! smallp
            itmp23(1:3)=( itmp11(1:3)+itmp13(1:3) )/2
            itmp24(1:3)=( itmp11(1:3)-itmp13(1:3) )/2
            itmp23(4)=itmp13(4)
            itmp23(5)=itmp13(6)
            itmp24(4)=itmp13(5)
            itmp24(5)=itmp13(7)
            i8=rspsnm(itmp23(1),itmp23(2),itmp23(3),itmp23(4),itmp23(5))
            i9=rspsnm(itmp24(1),itmp24(2),itmp24(3),itmp24(4),itmp24(5))

            stablehhpp(i1)%tpphh(i5,i4)=stablehhpp(i1)%tpphh(i5,i4)
     $+weight*stablehhpp(i1)%tpphhnew(i5,i4)
     $/(f(i6,i6)+f(i7,i7)-f(i8,i8)-f(i9,i9))

          enddo
        enddo
      enddo

      Ecorrnew=CCDEcnm()
      niter=niter+1
!      if (mod(niter,10)==0) then
        write(*,'(a12,i8)')  '  Iteration:',niter
        write(*,'(a11,f20.10)')  '    Ecorr =',Ecorrnew
!      endif

      if (abs(Ecorrnew-Ecorr)>deltaE) then
        Ecorr=Ecorrnew
        go to 204
      else
        Ecorr=Ecorrnew
      endif

!      write(*,*) ' HAHAHAHAHAHAHAHAHAHAHAHAHHAHAHAHA'


      end subroutine


      function CCDEcnm()
      use mod_sps
      implicit none
      integer :: i,j,k,l,a,b,c,d
      integer :: i1,i2,i3,i4,i5,i6
      real*8 :: a1
      real*8 :: CCDEcnm,ddot
      real*8,pointer :: vA(:),vB(:)

      CCDEcnm=0d0
      do i1=1,nhhpp
        i2=stablehhpp(i1)%nnhh
        i3=stablehhpp(i1)%nnpp
        allocate( vA(1:i2*i3),vB(1:i2*i3) )
        i6=0
        do i4=1,i2  ! ij
          do i5=1,i3  ! ab
            i6=i6+1
            vA(i6)=stablehhpp(i1)%vpphh(i5,i4)
            vB(i6)=stablehhpp(i1)%tpphh(i5,i4)
          enddo
        enddo
        a1=ddot(i6,vA,1,vB,1)
        deallocate( vA,vB )
        CCDEcnm=CCDEcnm+a1
      enddo
      CCDEcnm=CCDEcnm/4d0

      end function


      subroutine CCDterm145()
      implicit none
      integer :: i1,i2,i3,i4
      real*8,pointer :: A(:,:),B(:,:),C(:,:)

      do i1=1,nhhpp
        stablehhpp(i1)%tpphhnew=0d0
        i2=stablehhpp(i1)%nnhh
        i3=stablehhpp(i1)%nnpp
        allocate( C(i3,i2) )
        C=stablehhpp(i1)%vpphh

!        if (i1==2) then
!        do i4=1,16
!          write(*,'(8f7.3)') C(i4,:)
!        enddo
!        write(*,*) 
!        write(*,'(8i3)') btablehhpp(1:3,i1)
!        do i4=1,16
!          write(*,'(8i3)') stablehhpp(i1)%tablepp(:,i4)
!        enddo
!        write(*,*) 
!        do i4=1,8
!          write(*,'(8i3)') stablehhpp(i1)%tablehh(:,i4)
!        enddo
!        write(*,*) 
!        stop
!        endif

        allocate( A(i3,i3),B(i3,i2) )
        A=stablehhpp(i1)%vpppp
        B=stablehhpp(i1)%tpphh

!        if (i1==2) then
!        do i4=1,16
!          write(*,'(16f7.3)') A(i4,:)
!        enddo
!        write(*,*) 
!        do i4=1,16
!          write(*,'(8f9.5)') B(i4,:)
!        enddo
!        write(*,*) 
!        stop
!        endif


        call dgemm('N','N',i3,i2,i3,0.5d0,A,i3,B,i3,1d0,C,i3)

!        if (i1==2) then
!        do i4=1,16
!          write(*,'(8f7.3)') C(i4,:)
!        enddo
!        stop
!        write(*,*) 
!        endif

        deallocate( A )
        allocate( A(i2,i2) )
        A=stablehhpp(i1)%vhhhh
        B=stablehhpp(i1)%tpphh
 
!        if (i1==2) then
!        do i4=1,8
!          write(*,'(8f7.3)') A(i4,:)
!        enddo
!        write(*,*) 
!        do i4=1,16
!          write(*,'(8f7.3)') B(i4,:)
!        enddo
!        write(*,*) 
!        stop
!        endif

        call dgemm('N','N',i3,i2,i2,0.5d0,B,i3,A,i2,1d0,C,i3)
 
!        if (i1==2) then
!        do i4=1,16
!          write(*,'(8f7.3)') C(i4,:)
!        enddo
!        stop
!        write(*,*) 
!        endif

        stablehhpp(i1)%tpphhnew=stablehhpp(i1)%tpphhnew+C
        deallocate( A,B,C )
      enddo

      end subroutine


      subroutine CCDterm23()
      use mod_sps
      implicit none
      integer :: i1,i2,i3,i4,i5,i6,i7,i8
      integer :: itmp1(1:3),itmp2(1:7),itmp3(1:3),itmp4(1:7)

      do i1=1,nhhpp                                  ! term2
        itmp1(1:3)=btablehhpp(1:3,i1)                ! bigP
        i2=stablehhpp(i1)%nnhh
        i3=stablehhpp(i1)%nnpp
        do i4=1,i3  ! ab

          itmp2(1:7)=stablehhpp(i1)%tablepp(1:7,i4)  ! smallp
          itmp3(1:3)=( itmp1(1:3)-itmp2(1:3) )/2     ! pb
          i8=rspsnm(itmp3(1),itmp3(2),itmp3(3),itmp2(5),itmp2(7))
          do i5=1,i2  ! ij
            do i6=nspspl,nspspu  ! c
              if ( (itmp3(1)==sporbit(i6)%qnnm(1))
     $.and.(itmp3(2)==sporbit(i6)%qnnm(2))
     $.and.(itmp3(3)==sporbit(i6)%qnnm(3)) ) then
                itmp4(1:3)=itmp2(1:3)
                itmp4(4)=itmp2(4)
                itmp4(5)=sporbit(i6)%qnnm(4)
                itmp4(6)=itmp2(6)
                itmp4(7)=sporbit(i6)%qnnm(5)
                i7=stablehhpp(i1)%rpp( itmp4(1),itmp4(2),itmp4(3)
     $,itmp4(4),itmp4(5),itmp4(6),itmp4(7) )
                stablehhpp(i1)%tpphhnew(i4,i5)
     $=stablehhpp(i1)%tpphhnew(i4,i5)
     $+f(i8,i6)*stablehhpp(i1)%tpphh(i7,i5)
              endif
            enddo
          enddo

          itmp3(1:3)=( itmp1(1:3)+itmp2(1:3) )/2     ! pa
          i8=rspsnm(itmp3(1),itmp3(2),itmp3(3),itmp2(4),itmp2(6))
          do i5=1,i2  ! ij
            do i6=nspspl,nspspu  ! c
              if ( (itmp3(1)==sporbit(i6)%qnnm(1))
     $.and.(itmp3(2)==sporbit(i6)%qnnm(2))
     $.and.(itmp3(3)==sporbit(i6)%qnnm(3)) ) then
                itmp4(1:3)=-itmp2(1:3)
                itmp4(4)=itmp2(5)
                itmp4(5)=sporbit(i6)%qnnm(4)
                itmp4(6)=itmp2(7)
                itmp4(7)=sporbit(i6)%qnnm(5)
                i7=stablehhpp(i1)%rpp( itmp4(1),itmp4(2),itmp4(3)
     $,itmp4(4),itmp4(5),itmp4(6),itmp4(7) )
                stablehhpp(i1)%tpphhnew(i4,i5)
     $=stablehhpp(i1)%tpphhnew(i4,i5)
     $-f(i8,i6)*stablehhpp(i1)%tpphh(i7,i5)
              endif
            enddo
          enddo

        enddo
      enddo

      do i1=1,nhhpp                                  ! term3
        itmp1(1:3)=btablehhpp(1:3,i1)                ! bigP
        i2=stablehhpp(i1)%nnhh
        i3=stablehhpp(i1)%nnpp
        do i4=1,i2  ! ij

          itmp2(1:7)=stablehhpp(i1)%tablehh(1:7,i4)  ! smallp
          itmp3(1:3)=( itmp1(1:3)-itmp2(1:3) )/2     ! pj
          i8=rspsnm(itmp3(1),itmp3(2),itmp3(3),itmp2(5),itmp2(7))
          do i5=1,i3  ! ab
            do i6=1,nspshu  ! k
              if ( (itmp3(1)==sporbit(i6)%qnnm(1))
     $.and.(itmp3(2)==sporbit(i6)%qnnm(2))
     $.and.(itmp3(3)==sporbit(i6)%qnnm(3)) ) then
                itmp4(1:3)=itmp2(1:3)
                itmp4(4)=itmp2(4)
                itmp4(5)=sporbit(i6)%qnnm(4)
                itmp4(6)=itmp2(6)
                itmp4(7)=sporbit(i6)%qnnm(5)
                i7=stablehhpp(i1)%rhh( itmp4(1),itmp4(2),itmp4(3)
     $,itmp4(4),itmp4(5),itmp4(6),itmp4(7) )
                stablehhpp(i1)%tpphhnew(i5,i4)
     $=stablehhpp(i1)%tpphhnew(i5,i4)
     $-f(i6,i8)*stablehhpp(i1)%tpphh(i5,i7)
              endif
            enddo
          enddo

          itmp2(1:7)=stablehhpp(i1)%tablehh(1:7,i4)  ! smallp
          itmp3(1:3)=( itmp1(1:3)+itmp2(1:3) )/2     ! pi
          i8=rspsnm(itmp3(1),itmp3(2),itmp3(3),itmp2(4),itmp2(6))
          do i5=1,i3  ! ab
            do i6=1,nspshu  ! k
              if ( (itmp3(1)==sporbit(i6)%qnnm(1))
     $.and.(itmp3(2)==sporbit(i6)%qnnm(2))
     $.and.(itmp3(3)==sporbit(i6)%qnnm(3)) ) then
                itmp4(1:3)=-itmp2(1:3)
                itmp4(4)=itmp2(5)
                itmp4(5)=sporbit(i6)%qnnm(4)
                itmp4(6)=itmp2(7)
                itmp4(7)=sporbit(i6)%qnnm(5)
                i7=stablehhpp(i1)%rhh( itmp4(1),itmp4(2),itmp4(3)
     $,itmp4(4),itmp4(5),itmp4(6),itmp4(7) )
                stablehhpp(i1)%tpphhnew(i5,i4)
     $=stablehhpp(i1)%tpphhnew(i5,i4)
     $+f(i6,i8)*stablehhpp(i1)%tpphh(i5,i7)
              endif
            enddo
          enddo

        enddo
      enddo

      end subroutine


!      function CCDterm6(i,j,a,b)  ! Not useful for pairing problems


      end module mod_CCDnuclearmatter
