      module mod_mergesort
      save
      integer,pointer :: spl(:,:)
      real*8,pointer :: aspe(:)


      contains
      subroutine Merger(left,mid,right)
      implicit none
      integer :: left,mid,right
      integer :: length,i,j,k,ind
      integer,pointer :: temp(:,:)
      real*8,pointer :: atemp(:)

      length=right-left+1
      allocate( atemp(1:length),temp(1:5,1:length) )
      ind=1
      i=left
      j=mid+1

      do while ( (i<=mid).and.(j<=right) )
        if (aspe(i)<=aspe(j)) then
          atemp(ind)=aspe(i)
          temp(1:5,ind)=spl(1:5,i)
          ind=ind+1
          i=i+1
        else
          atemp(ind)=aspe(j)
          temp(1:5,ind)=spl(1:5,j)
          ind=ind+1
          j=j+1
        endif
      enddo

      do while (i<=mid)
        atemp(ind)=aspe(i)
        temp(1:5,ind)=spl(1:5,i)
        ind=ind+1
        i=i+1
      enddo

      do while (j<=right)
        atemp(ind)=aspe(j)
        temp(1:5,ind)=spl(1:5,j)
        ind=ind+1
        j=j+1
      enddo

      do k=1,length
        aspe(left)=atemp(k)
        spl(1:5,left)=temp(1:5,k)
        left=left+1
      enddo

      end subroutine


      subroutine MergeSorter(length)
      implicit none
      integer :: left,mid,right
      integer :: length,i,j,k

      i=1
      do while (i<length)
        left=1
        do while (left+i<=length)
          mid=left+i-1
          if (mid+i<=length) then
            right=mid+i
          else
            right=length
          endif
          call Merger(left,mid,right)
          left=right+1
        enddo
        i=i*2
      enddo
      
      end subroutine


      end module mod_mergesort


      program generate sps and tbme for nuclear matter
      use mod_mergesort
      implicit none
      real*8,parameter :: Pi=3.14159265358979323846d0
      real*8,parameter :: mp=938.2720813d0
!      real*8,parameter :: mn=939.5654133d0
!      real*8,parameter :: hbar=197.326978795182d0
      real*8,parameter :: mn=938.92d0    ! for benchmarks
      real*8,parameter :: hbar=197.33d0  ! for benchmarks
      integer :: i0,i1,i2,i3,i4,i5,i6,i7,i8
      real*8 :: a1,a2,a3,a4
      integer :: itmpl(1:4),itmpr(1:4)
      real*8 :: atmpl(1:3),atmpr(1:3)
      !------------------ end implicit -------------------------------
      integer :: nsps
      real*8 :: rho,kF,L,thetax,thetay,thetaz
      real*8 :: Minnesota
      integer :: A,Nmax
      character(len=1) :: c1
      integer :: ntb,ntbme,nfock
      real*8,pointer :: smallp(:,:)
      integer,pointer :: st(:,:),smallpi(:,:)
      integer,pointer :: tbmel(:,:),fockl(:,:)
      real*8,pointer :: tbmev(:),fockv(:)

      open(unit=101,file='parameter.dat')
      read(101,*)
      read(101,*)
      read(101,*)
      read(101,*) A,rho,Nmax
      read(101,*) c1  ! nuclear matter type
      read(101,*) thetax,thetay,thetaz
      close(101)
      thetax=thetax*Pi
      thetay=thetay*Pi
      thetaz=thetaz*Pi
      L=(dfloat(A)/rho)**(0.33333333333333333333333333d0)
      kF=(3d0*rho*Pi*Pi)**(0.33333333333333333333333333d0)
      a1=Nmax*2d0*Pi/L
      write(*,'(a21,f15.7)') '  Fermi momentum kF =',kF
      write(*,'(a16,f15.7)') '  Box length L =',L
      write(*,'(a9,f15.7)') '  knmax =',a1

      if ( (c1=='n').or.(c1=='N') ) then
        nsps=(2*Nmax+1)
        nsps=nsps*nsps*nsps*2
        allocate( spl(1:5,1:nsps) )
        i0=0
        do i1=0,3*Nmax*Nmax
          do i2=-Nmax,Nmax  ! nz
          do i3=-Nmax,Nmax  ! ny
          do i4=-Nmax,Nmax  ! nx
            i5=i2*i2+i3*i3+i4*i4
            if (i5==i1) then
              i0=i0+1
              spl(:,i0)=[i4,i3,i2,1,1]
              i0=i0+1
              spl(:,i0)=[i4,i3,i2,-1,1]
            endif
          enddo
          enddo
          enddo
        enddo
        if (i0/=nsps) then
          write(*,*) 'Single-particle states wrong!!! STOP! '
          stop
        endif
      elseif ( (c1=='s').or.(c1=='S') ) then
        nsps=(2*Nmax+1)
        nsps=nsps*nsps*nsps*4
        allocate( spl(1:5,1:nsps) )
        i0=0
        do i1=0,3*Nmax*Nmax
          do i2=-Nmax,Nmax  ! nz
          do i3=-Nmax,Nmax  ! ny
          do i4=-Nmax,Nmax  ! nx
            i5=i2*i2+i3*i3+i4*i4
            if (i5==i1) then
              i0=i0+1
              spl(:,i0)=[i4,i3,i2,1,1]
              i0=i0+1
              spl(:,i0)=[i4,i3,i2,-1,1]
              i0=i0+1
              spl(:,i0)=[i4,i3,i2,1,-1]
              i0=i0+1
              spl(:,i0)=[i4,i3,i2,-1,-1]
            endif
          enddo
          enddo
          enddo
        enddo
        if (i0/=nsps) then
          write(*,*) 'Single-particle states wrong!!! STOP! '
          stop
        endif
      endif

      allocate( aspe(1:nsps) )
      a1=2d0*mn*L*L
      a1=hbar*hbar/a1
      a2=2d0*mp*L*L
      a2=hbar*hbar/a2
      do i1=1,nsps
        a3=0d0
        a4=2d0*Pi*dfloat(spl(1,i1))+thetax
        a3=a3+a4*a4
        a4=2d0*Pi*dfloat(spl(2,i1))+thetay
        a3=a3+a4*a4
        a4=2d0*Pi*dfloat(spl(3,i1))+thetaz
        a3=a3+a4*a4
        if (spl(5,i1)==1) then
          aspe(i1)=a3*a1
        elseif (spl(5,i1)==-1) then
          aspe(i1)=a3*a2
        endif
      enddo

      call MergeSorter(nsps)
      open(unit=101,file='sps.b')
      write(101,'(i6,a21)') nsps,' nx ny nz sz tz spe'
      do i1=1,nsps
        write(101,'(5i5,f25.15)') spl(1:5,i1),aspe(i1)
      enddo
      close(101)

      if ( (c1=='n').or.(c1=='N') ) then
        ntb=(4*Nmax+1)*(4*Nmax+1)*(4*Nmax+1)*4
      elseif ( (c1=='s').or.(c1=='S') ) then
        ntb=(4*Nmax+1)*(4*Nmax+1)*(4*Nmax+1)*16
      endif
      allocate( smallp(1:3,1:ntb),smallpi(1:3,1:ntb),st(1:4,1:ntb) )
      i1=0
      do i2=-2*Nmax,2*Nmax  ! nz
        a3=2d0*Pi*dfloat(i2)/L
      do i3=-2*Nmax,2*Nmax  ! ny
        a2=2d0*Pi*dfloat(i3)/L
      do i4=-2*Nmax,2*Nmax  ! nx
        a1=2d0*Pi*dfloat(i4)/L
        if ( (c1=='n').or.(c1=='N') ) then
          do i5=1,-1,-2
          do i6=1,-1,-2
            i1=i1+1
            smallp(1,i1)=a1
            smallp(2,i1)=a2
            smallp(3,i1)=a3
            smallpi(1:3,i1)=[i4,i3,i2]
            st(1:4,i1)=[i5,i6,1,1]
          enddo
          enddo
        elseif ( (c1=='s').or.(c1=='S') ) then
          do i5=1,-1,-2
          do i6=1,-1,-2
          do i7=1,-1,-2
          do i8=1,-1,-2
            i1=i1+1
            smallp(1,i1)=a1
            smallp(2,i1)=a2
            smallp(3,i1)=a3
            smallpi(1:3,i1)=[i4,i3,i2]
            st(1:4,i1)=[i5,i6,i7,i8]
          enddo
          enddo
          enddo
          enddo
        endif
      enddo
      enddo
      enddo
      if (i1/=ntb) then
        write(*,*) ' Two body wrong!!! STOP! '
        stop
      endif
      write(*,*) ' S.p. and t.b. states ready...' 

      ntbme=0
      do i1=1,ntb
        do i2=i1,ntb
          if ( (mod(abs(smallpi(1,i1)-smallpi(1,i2)),2)/=0)
     $.or.(mod(abs(smallpi(2,i1)-smallpi(2,i2)),2)/=0)
     $.or.(mod(abs(smallpi(3,i1)-smallpi(3,i2)),2)/=0)
     $.or.((st(1,i1)+st(2,i1))/=(st(1,i2)+st(2,i2)))
     $.or.((st(3,i1)+st(4,i1))/=(st(3,i2)+st(4,i2))) ) go to 202
          ntbme=ntbme+1
202     enddo
      enddo
      write(*,'(a33,i10)') ' Number of possible nonzero tbme:',ntbme
      allocate( tbmel(1:14,1:ntbme),tbmev(1:ntbme) )

      ntbme=0
      do i1=1,ntb
        do i2=i1,ntb
          if ( (mod(abs(smallpi(1,i1)-smallpi(1,i2)),2)/=0)
     $.or.(mod(abs(smallpi(2,i1)-smallpi(2,i2)),2)/=0)
     $.or.(mod(abs(smallpi(3,i1)-smallpi(3,i2)),2)/=0)
     $.or.((st(1,i1)+st(2,i1))/=(st(1,i2)+st(2,i2)))
     $.or.((st(3,i1)+st(4,i1))/=(st(3,i2)+st(4,i2))) ) go to 201
          atmpl(1:3)=smallp(1:3,i1)
          atmpr(1:3)=smallp(1:3,i2)
          itmpl(1:4)=st(1:4,i1)
          itmpr(1:4)=st(1:4,i2)
          a1=Minnesota(atmpl,itmpl,atmpr,itmpr,L)
          atmpr(1:3)=-smallp(1:3,i2)
          itmpr(1:4)=[st(2,i2),st(1,i2),st(4,i2),st(3,i2)]
          a1=a1-Minnesota(atmpl,itmpl,atmpr,itmpr,L)
          if (abs(a1)>1d-13)  then
            ntbme=ntbme+1
            tbmel(1:14,ntbme)=[smallpi(1,i1),smallpi(2,i1),smallpi(3,i1)
     $,st(1,i1),st(2,i1),st(3,i1),st(4,i1)
     $,smallpi(1,i2),smallpi(2,i2),smallpi(3,i2)
     $,st(1,i2),st(2,i2),st(3,i2),st(4,i2)]
            tbmev(ntbme)=a1
          endif
201     enddo
      enddo
      write(*,'(a33,i10)') ' Number of nonzero tbme:',ntbme

      open(unit=101,file='TBME.b')
      write(101,*) '###'
      write(101,*) '###'
      write(101,*) '###'
      write(101,'(i10)',advance='no') ntbme
!      do i1=1,nsps
!        write(101,'(f25.15)',advance='no') aspe(i1)
!      enddo
      write(101,*) 
      write(101,'(a48,a25)') 
     $' 2px 2py 2pz s1 s2 t1 t2 2px 2py 2pz s3 s4 t3 t4','v'
      do i1=1,ntbme
        write(101,'(3i4,4i3,3i4,4i3,f25.15)') 
     $tbmel(1:14,i1),tbmev(i1)
      enddo
      close(101)

      write(*,'(a39)') ' Calculating fock matrix and Ehf...'
      a2=2d0*Pi/L
      nfock=0
      allocate( fockl(1:2,1:nsps*nsps),fockv(1:nsps*nsps) )
      do i1=1,nsps
        do i2=1,nsps
          a1=0d0
          if ( (spl(1,i1)==spl(1,i2)).and.(spl(2,i1)==spl(2,i2))
     $.and.(spl(3,i1)==spl(3,i2)) ) then
            do i3=1,A
              atmpl(1)=dfloat(spl(1,i1)-spl(1,i3))*a2
              atmpl(2)=dfloat(spl(2,i1)-spl(2,i3))*a2
              atmpl(3)=dfloat(spl(3,i1)-spl(3,i3))*a2
              itmpl(1)=spl(4,i1)
              itmpl(2)=spl(4,i3)
              itmpl(3)=spl(5,i1)
              itmpl(4)=spl(5,i3)
              atmpr(1)=dfloat(spl(1,i2)-spl(1,i3))*a2
              atmpr(2)=dfloat(spl(2,i2)-spl(2,i3))*a2
              atmpr(3)=dfloat(spl(3,i2)-spl(3,i3))*a2
              itmpr(1)=spl(4,i2)
              itmpr(2)=spl(4,i3)
              itmpr(3)=spl(5,i2)
              itmpr(4)=spl(5,i3)
              a1=a1+Minnesota(atmpl,itmpl,atmpr,itmpr,L)
              atmpr(1)=-dfloat(spl(1,i2)-spl(1,i3))*a2
              atmpr(2)=-dfloat(spl(2,i2)-spl(2,i3))*a2
              atmpr(3)=-dfloat(spl(3,i2)-spl(3,i3))*a2
              itmpr(1)=spl(4,i3)
              itmpr(2)=spl(4,i2)
              itmpr(3)=spl(5,i3)
              itmpr(4)=spl(5,i2)
              a1=a1-Minnesota(atmpl,itmpl,atmpr,itmpr,L)
            enddo
          endif
          if (i1==i2) a1=a1+aspe(i1)
          if (abs(a1)>1d-13) then
            nfock=nfock+1
            fockl(1:2,nfock)=[i1,i2]
            fockv(nfock)=a1
          endif
        enddo
      enddo

      a2=2d0*Pi/L
      a3=0d0
      do i1=1,A
        do i2=1,A
          atmpl(1)=dfloat(spl(1,i1)-spl(1,i2))*a2
          atmpl(2)=dfloat(spl(2,i1)-spl(2,i2))*a2
          atmpl(3)=dfloat(spl(3,i1)-spl(3,i2))*a2
          itmpl(1)=spl(4,i1)
          itmpl(2)=spl(4,i2)
          itmpl(3)=spl(5,i1)
          itmpl(4)=spl(5,i2)
          atmpr=atmpl
          itmpr=itmpl
          a3=a3+Minnesota(atmpl,itmpl,atmpr,itmpr,L)/2d0
          atmpr(1)=-atmpr(1)
          atmpr(2)=-atmpr(2)
          atmpr(3)=-atmpr(3)
          itmpr(1)=spl(4,i2)
          itmpr(2)=spl(4,i1)
          itmpr(3)=spl(5,i2)
          itmpr(4)=spl(5,i1)
          a3=a3-Minnesota(atmpl,itmpl,atmpr,itmpr,L)/2d0
        enddo
        a3=a3+aspe(i1)
      enddo

      open(unit=101,file='fock.b')
      write(101,'(i10,a8,f25.15)') nfock,'Ehf = ',a3
      do i1=1,nfock
        write(101,'(2i5,f25.15)') fockl(1:2,i1),fockv(i1)
      enddo
      close(101)

      deallocate( aspe,spl,smallp,st,smallpi )

      stop
      end


      function Minnesota(spl,stl,spr,str,L)
      implicit none
      real*8,parameter :: Pi=3.14159265358979323846d0
      real*8,parameter :: VR=200d0
      real*8,parameter :: VS=-91.85d0
      real*8,parameter :: VT=-178d0
!      real*8,parameter :: VR=0d0
!      real*8,parameter :: VS=0d0
!      real*8,parameter :: VT=0d0
      real*8,parameter :: kapaR=1.487d0
      real*8,parameter :: kapaS=0.465d0
      real*8,parameter :: kapaT=0.639d0
      integer :: i1,i2,i3
      real*8 :: a1,a2
      real*8 :: atmp(1:3),q2,L3
      real*8 :: spl(1:3),spr(1:3)
      integer :: stl(1:4),str(1:4)   ! szp,szq,tzp,tzq
      real*8 :: L,Minnesota

      Minnesota=0d0

      i1=0
      if ( (stl(1)==str(1)).and.(stl(2)==str(2)).and.(stl(3)==str(3))
     $.and.(stl(4)==str(4)) ) i1=i1+1
      if ( (stl(1)==str(2)).and.(stl(2)==str(1)).and.(stl(3)==str(4))
     $.and.(stl(4)==str(3)) ) i1=i1-1

      i2=0
      if ( (stl(1)==str(1)).and.(stl(2)==str(2)).and.(stl(3)==str(3))
     $.and.(stl(4)==str(4)) ) i2=i2+1
      if ( (stl(1)==str(2)).and.(stl(2)==str(1)).and.(stl(3)==str(3))
     $.and.(stl(4)==str(4)) ) i2=i2+1
      if ( (stl(1)==str(1)).and.(stl(2)==str(2)).and.(stl(3)==str(4))
     $.and.(stl(4)==str(3)) ) i2=i2-1
      if ( (stl(1)==str(2)).and.(stl(2)==str(1)).and.(stl(3)==str(4))
     $.and.(stl(4)==str(3)) ) i2=i2-1

      i3=0
      if ( (stl(1)==str(1)).and.(stl(2)==str(2)).and.(stl(3)==str(3))
     $.and.(stl(4)==str(4)) ) i3=i3+1
      if ( (stl(1)==str(2)).and.(stl(2)==str(1)).and.(stl(3)==str(3))
     $.and.(stl(4)==str(4)) ) i3=i3-1
      if ( (stl(1)==str(1)).and.(stl(2)==str(2)).and.(stl(3)==str(4))
     $.and.(stl(4)==str(3)) ) i3=i3+1
      if ( (stl(1)==str(2)).and.(stl(2)==str(1)).and.(stl(3)==str(4))
     $.and.(stl(4)==str(3)) ) i3=i3-1

      if ( (i1/=0).or.(i2/=0).or.(i3/=0) ) then
        atmp(1)=(spl(1)-spr(1))/2d0
        atmp(2)=(spl(2)-spr(2))/2d0
        atmp(3)=(spl(3)-spr(3))/2d0
        q2=atmp(1)*atmp(1)+atmp(2)*atmp(2)+atmp(3)*atmp(3)
        L3=1d0/(L*L*L)
      endif
      if (i1/=0) then
        Minnesota=Minnesota+(VR*L3)*((Pi/kapaR)**(1.5d0))
     $*dexp(-q2/kapaR/4d0)*dfloat(i1)/2d0
      endif
      if (i2/=0) then
        Minnesota=Minnesota+(VT*L3)*((Pi/kapaT)**(1.5d0))
     $*dexp(-q2/kapaT/4d0)*dfloat(i2)/4d0
      endif
      if (i3/=0) then
        Minnesota=Minnesota+(VS*L3)*((Pi/kapaS)**(1.5d0))
     $*dexp(-q2/kapaS/4d0)*dfloat(i3)/4d0
      endif

      end function


