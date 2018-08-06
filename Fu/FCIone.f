      include 'calH.f'

      program FCIone
      use mod_sps
      use mod_basis
      use mod_Hmatrix
      implicit none
      integer :: i1,i2
      integer :: ntmp
      real*8 :: a1,a2
      real*8,pointer :: eigenv(:),work(:)
      character(len=1) :: c1

      write(*,*) '#####################################################'
      write(*,*) '#               This is FCIone program              #'
      write(*,*) '#####################################################'

      open(unit=101,file='parameter.dat')
      read(101,*)
      read(101,*) c1
      close(101)

      call readsps_pairing()
      if (  (c1=='e').or.(c1=='E') ) then
        call setbasis_exact()
      elseif (  (c1=='p').or.(c1=='P') ) then
        call setbasis_pairing()
      endif
      write(*,*) 'Basis vectors generated... '
      write(*,'(a13,i11)') '  dimension =',nbasis

      call calHmatrix_SM()
      write(*,*) 'H matrix elements calculated... '

      allocate( eigenv(1:nbasis),work(1:3*nbasis) )
      a1=H(1,1)
!      write(*,*) a1,H(1,1:10)
      call dsyev('N','U',nbasis,H,nbasis,eigenv,work,3*nbasis,ntmp)
      a2=eigenv(1)

      open(unit=101,file='level_FCIone.dat')
      write(101,'(a11,i11)') 'dimension =',nbasis
      write(*,'(5x,3a15)') 'Eigenv','Ex','DeltaE'
      write(101,'(5x,3a15)') 'Eigenv','Ex','DeltaE'
      do i1=1,min(100,nbasis)
        write(*,'(i5)',advance='no') i1
        write(*,'(f15.7)',advance='no') eigenv(i1)
        write(*,'(f15.7)',advance='no') eigenv(i1)-a2
        write(*,'(f15.7)',advance='no') eigenv(i1)-a1
        write(*,*)
        write(101,'(i5)',advance='no') i1
        write(101,'(f15.7)',advance='no') eigenv(i1)
        write(101,'(f15.7)',advance='no') eigenv(i1)-a2
        write(101,'(f15.7)',advance='no') eigenv(i1)-a1
        write(101,*)
      enddo
      close(101)

!      do i1=1,nbasis
!        write(*,'(i8)',advance='no') i1
!        do i2=1,nsps
!          write(*,'(l2)',advance='no') btest(vectorb(i1),i2) 
!        enddo
!        write(*,*)
!      enddo


      deallocate( eigenv,work)
      stop
      end
