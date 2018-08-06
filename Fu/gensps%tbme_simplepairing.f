      program generate sps and tbme for simple pairing problem
      implicit none
      integer :: i1,i2
      real*8 :: a1
      integer :: nsps,ntbme
      integer :: a(1:100000),b(1:100000),c(1:100000),d(1:100000)
      real*8 :: aspe(1:100000),atbme(1:100000)
      real*8 :: delta,g
      integer :: result,systemqq

      write(*,*) 'Please input delta and g: '
      read(*,*) delta,g

      open(unit=101,file='sps.b')
      read(101,*) nsps
      close(101)

      open(unit=101,file='sps.b')
      write(101,'(i3)') nsps
      do i1=1,nsps
        write(101,'(2i5)') (i1+1)/2,2*mod(i1,2)-1
      enddo
      close(101)

      do i1=1,nsps
        aspe(i1)=((i1-1)/2)*delta
      enddo
      ntbme=0
      do i1=1,nsps/2
        do i2=i1,nsps/2
          ntbme=ntbme+1
          a(ntbme)=i1*2-1
          b(ntbme)=i1*2
          c(ntbme)=i2*2-1
          d(ntbme)=i2*2
          atbme(ntbme)=-g
        enddo
      enddo
      
      open(unit=101,file='TBME.b')
      write(101,*) '###'
      write(101,*) '###'
      write(101,*) '###'
      write(101,'(i3)',advance='no') ntbme
      do i1=1,nsps
        write(101,'(f10.5)',advance='no') aspe(i1)
      enddo
      write(101,*) 
      write(101,'(a31)') '   a   b   c   d                        v'
      do i1=1,ntbme
        write(101,'(4i4,f25.15)') a(i1),b(i1),c(i1),d(i1),atbme(i1) 
      enddo
      close(101)

      stop
      end
