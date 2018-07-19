      program test
      implicit none
      integer :: i1
      integer*8 :: a,b,c,d
      integer :: e

      a=0
      b=0
      a=ibset(a,3)
      a=ibset(a,4)
      a=ibset(a,5)
      b=ibset(b,8)
      c=a+b
      d=ibits(a,3,3)
      e=1
      

      do i1=0,64
        write(*,*) i1,btest(a,i1),btest(b,i1),btest(c,i1),btest(d,i1)
     $,btest(e,i1)
      enddo

      end
