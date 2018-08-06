      include 'calCCD.f'

      program CCone
      use mod_sps
      use mod_CCDcommon
      use mod_CCDsimplepairing
      use mod_CCDnuclearmatter
      implicit none
      integer :: i1,i2
      real*8 :: a1,a2
      character(len=1) :: c1

      write(*,*) '#####################################################'
      write(*,*) '#                This is CCone program              #'
      write(*,*) '#####################################################'

      open(unit=101,file='parameter.dat')
      read(101,*) i1,deltaE,weight,Delta
      read(101,*) c1
      close(101)

      if ( (c1=='p').or.(c1=='P') ) then
        call readsps_pairing()
        call calCCDsimplepairing()
        open(unit=101,file='level_CCone.dat')
        write(*,'(a12,i8)')  '  Iteration:',niter
        write(*,'(a11,f25.15)')  '    Ecorr =',Ecorr
        write(101,'(a7,f25.15)') '  Ehf =',Ehf 
        write(101,'(a14,f25.15)') '  Init Ecorr =',Embpt2
        write(101,'(a12,i8)')  '  Iteration:',niter
        write(101,'(a11,f25.15)')  '    Ecorr =',Ecorr
        close(101)
      elseif ( (c1=='n').or.(c1=='N') ) then
        call readsps_nuclearmatter()
        call calCCDnuclearmatter()
        open(unit=101,file='level_CCone.dat')
        write(*,'(a12,i8)')  '  Iteration:',niter
        write(*,'(a11,f25.15)')  '    Ecorr =',Ecorr
        write(101,'(a14,f25.15)') '  Init Ecorr =',Embpt2
        write(101,'(a12,i8)')  '  Iteration:',niter
        write(101,'(a11,f25.15)')  '    Ecorr =',Ecorr
        close(101)
      else
        write(*,*) 'No CC function for other uses! STOP! '
        stop
      endif

      stop
      end
