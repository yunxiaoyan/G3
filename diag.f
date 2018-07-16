      program using diagonalization
      implicit none
      !----------------- begin implicit -------------------------------
      integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9
      real*8 :: a1,a2,a3,a4,a5,a6,a7,a8,a9
      integer :: ntmp
      !------------------ end implicit -------------------------------
      real*8 :: delta,g
      real*8 :: H(1:2,1:2),eigenv(1:2),work(1:6)

      write(*,*) 'Input delta & g:'
      read(*,*) delta,g

      H=0d0
      H(1,1)=-g
      H(1,2)=-g
      H(2,1)=-g
      H(2,2)=2d0*delta-g
      
      call dsyev('V','U',2,H,2,eigenv,work,3*2,ntmp)
      write(*,*) eigenv

      end
