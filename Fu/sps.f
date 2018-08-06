      module mod_sps
      save
      integer :: nsps  ! number of sp orbit
      logical :: lspspairing
      type sps
        integer :: qnn  ! quantum number n
        integer :: qnm  ! quantum number 2m
        integer :: qnnm(1:5) ! quantum number for nuclear matter nx,ny,nz,sz,tz
!        integer :: qnl  ! quantum number l
!        integer :: qnj  ! quantum number 2j
!        integer :: qntau  ! quantum number 2tau
!        integer :: qnw  ! quantum number w
      end type
      type(sps),pointer :: sporbit(:)
      integer*8,pointer :: sporbitb(:)
      integer,pointer :: rspsnm(:,:,:,:,:)


      contains
      subroutine readsps_pairing()
      implicit none
      integer :: i1
      open(unit=101,file='sps.b')
      read(101,*) nsps
      allocate(sporbit(1:nsps),sporbitb(1:nsps))
      sporbit(:)%qnn=-999999
      sporbit(:)%qnm=-999999
      sporbitb=0
      do i1=1,nsps
        read(101,*) sporbit(i1)%qnn,sporbit(i1)%qnm
        sporbitb(i1)=ibset(sporbitb(i1),i1)
      enddo
      close(101)

      lspspairing=.true.
      if (mod(nsps,2)/=0) then
        lspspairing=.false.
        go to 201
      endif
      do i1=2,nsps,2
        if ( (sporbit(i1-1)%qnn/=sporbit(i1)%qnn)
     $.or.(sporbit(i1-1)%qnm/=-sporbit(i1)%qnm) ) then
          lspspairing=.false.
          go to 201
        endif
      enddo

201   continue
      end subroutine


      subroutine readsps_nuclearmatter()
      implicit none
      integer :: i1,i2,i3,i4,i5,i6,Nmax

      open(unit=101,file='parameter.dat')
      read(101,*)
      read(101,*)
      read(101,*)
      read(101,*) i1,i1,Nmax
      close(101)
      allocate( rspsnm(-Nmax:Nmax,-Nmax:Nmax,-Nmax:Nmax,-1:1,-1:1) )
      rspsnm=0

      open(unit=101,file='sps.b')
      read(101,*) nsps
      allocate( sporbit(1:nsps) )
      do i1=1,nsps
        read(101,*) i2,i3,i4,i5,i6
        sporbit(i1)%qnnm(1:5)=[i2,i3,i4,i5,i6]
        rspsnm(i2,i3,i4,i5,i6)=i1
      enddo
      close(101)

      end subroutine


      end module mod_sps
