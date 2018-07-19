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

      call readsps_pairing()
      call setbasis_exact()
      call calHmatrix_SM()


      allocate( eigenv(1:nbasis),work(1:3*nbasis) )
      a1=H(1,1)
      call dsyev('N','U',nbasis,H,nbasis,eigenv,work,3*nbasis,ntmp)
      write(*,*) eigenv(1),eigenv(1)-a1

      do i1=1,nbasis
        write(*,'(i8)',advance='no') i1
        do i2=1,nsps
          write(*,'(l2)',advance='no') btest(vectorb(i1),i2) 
        enddo
        write(*,*)
      enddo


      deallocate( eigenv,work)
      stop
      end
