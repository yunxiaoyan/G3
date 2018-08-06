      module mod_mergesort
      save
      real*8,pointer :: A(:)


      contains
      subroutine Merger(left,mid,right)
      implicit none
      integer :: left,mid,right
      integer :: length,i,j,k,ind
      real*8,pointer :: temp(:)

      length=right-left+1
      allocate( temp(1:length) )
      ind=1
      i=left
      j=mid+1

      do while ( (i<=mid).and.(j<=right) )
        if (A(i)<=A(j)) then
          temp(ind)=A(i)
          ind=ind+1
          i=i+1
        else
          temp(ind)=A(j)
          ind=ind+1
          j=j+1
        endif
      enddo

      do while (i<=mid)
        temp(ind)=A(i)
        ind=ind+1
        i=i+1
      enddo

      do while (j<=right)
        temp(ind)=A(j)
        ind=ind+1
        j=j+1
      enddo

      do k=1,length
        A(left)=temp(k)
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


      program test
      use mod_mergesort
      implicit none
      integer :: i1
      integer :: length
      real*8 :: a1
      real*8,pointer :: B(:)

      write(*,*) 'Input length of the vector:'
      read(*,*) length

      allocate( A(1:length),B(1:length) )
      call random_seed()
      do i1=1,length
        call random_number(a1)
        a1=a1-0.5d0
        A(i1)=a1
        B(i1)=a1
      enddo
      call MergeSorter(length)

      do i1=1,length
        write(*,*) A(i1),B(i1)
      enddo

      end



