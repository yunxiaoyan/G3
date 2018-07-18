      include 'sps.f'

      module mod_basis
      save
      integer :: nparticle  ! number of particles
      integer :: nbasis  ! number of states
      integer*8,pointer :: vectorb(:)
      integer,pointer :: vectori(:,:)

      
      contains
      subroutine setbasis_exact(totalm)
      use mod_sps
      implicit none 
      integer :: i1,i2,i3
      integer :: totalm
      integer,pointer :: vectortmp(:)

      open(unit=101,file='parameter.dat')
      read(101,*) nparticle
      close(101)
      nbasis=1
      do i1=nparticle+1,nsps
        nbasis=nbasis*i1
      enddo
      do i1=1,nparticle
        nbasis=nbasis/i1
      enddo
      
      allocate( vectorb(1:nbasis),vectori(1:nparticle,1:nbasis) )
      allocate( vectortmp(1:nparticle) )
      vectorb=0
      vectori=0
      vectortmp=0
      i1=1
      do i2=1,nparticle
        vectortmp(i2)=i2
      enddo
      vectori(:,1)=vectortmp

201   i2=nparticle
      i1=i1+1
      if (vectortmp(i2)==nsps) then
        go to 202
      else
        go to 203
      endif

202   if (i2-1<1) then
        go to 204
      else
        go to 205
      endif

203   vectortmp(i2)=vectortmp(i2)+1
      do i3=1,nparticle-i2
        vectortmp(i2+i3)=vectortmp(i2)+i3
      enddo
      vectori(:,i1)=vectortmp
      go to 201

205   i2=i2-1
      if (vectortmp(i2)+1==vectortmp(i2+1)) then
        go to 202
      else
        go to 203
      endif
      
204   continue
      if (i1/=nbasis+1) write(*,*) 'Generating basis wrong!'

      do i1=1,nbasis
        do i2=1,nparticle
          vectorb(i1)=ibset(vectorb(i1),vectori(i2,i1))
        enddo
      enddo


      deallocate(vectortmp)
      end subroutine


      end module mod_basis
