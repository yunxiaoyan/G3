      include 'calCCD.f'

      program CCone
      use mod_sps
      use mod_CCDsimple
      implicit none
      integer :: i1,i2
      real*8 :: a1,a2

      write(*,*) '#####################################################'
      write(*,*) '#              This is CCsimple program             #'
      write(*,*) '#####################################################'

      call calCCDsimple()

      open(unit=101,file='level_CCsimple.dat')
      write(*,'(a12,i8)')  '  Iteration:',niter
      write(*,'(a11,f20.10)')  '    Ecorr =',Ecorr
      write(101,'(a12,i8)')  '  Iteration:',niter
      write(101,'(a11,f20.10)')  '    Ecorr =',Ecorr
      close(101)

      end
