      include 'sps.f'

      module mod_basis
      save
      integer,parameter :: nbasistmp=16000
      integer :: nparticle  ! number of particles
      integer :: nbasis  ! number of states
      integer*8,pointer :: vectorb(:)
      integer,pointer :: vectori(:,:)

      
      contains
      subroutine setbasis_exact()  ! Not using totalm yet!!!
      use mod_sps
      implicit none 
      integer :: i1,i2,i3
!      integer :: totalm
      integer,pointer :: vectortmp(:)

      open(unit=101,file='parameter.dat')
      read(101,*) nparticle
!      read(101,*) totalm
      close(101)
      nbasis=1
      do i1=nsps-nparticle+1,nsps
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

201   vectori(:,i1)=vectortmp
      i2=nparticle
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


      deallocate(vectortmp,vectori)
      end subroutine


      subroutine setbasis_pairing()  ! pairing in quantum number m
      use mod_sps
      implicit none 
      integer :: i1,i2,i3,i4,i5
!      integer :: totalm
      integer,pointer :: vectortmp(:)

      write(*,*) 'Please make sure your sps is orgnized in pairs 
     $12,34,56,... '
      if (.not.lspspairing) then
        write(*,*) 'Your sps is not orgnized in pairs!!! STOP! '
        stop
      endif

      open(unit=101,file='parameter.dat')
      read(101,*) nparticle
!      read(101,*) totalm
      close(101)
      nbasis=1
      do i1=nsps-nparticle+1,nsps
        nbasis=nbasis*i1
      enddo
      do i1=1,nparticle
        nbasis=nbasis/i1
      enddo
      
      allocate( vectortmp(1:nparticle) )
      allocate( vectori(1:nparticle,1:nbasistmp) )
      vectori=0
      vectortmp=0
      i5=mod(nparticle,2)
      i1=1
      do i2=1,nparticle
        vectortmp(i2)=i2
      enddo

      !----- Check pairing when building basis ---------- 
201   i3=1
      i4=0
      do while (i3<=nparticle)
        if (i3==nparticle) then
          i4=i4+1
          i3=i3-1
        elseif ( ((vectortmp(i3)+1)/2)/=((vectortmp(i3+1)+1)/2) ) then
          i4=i4+1
          i3=i3-1
        endif
        i3=i3+2
      enddo
      if (i4/=i5) go to 206
      !--- End check pairing when building basis ---------- 
      vectori(:,i1)=vectortmp
      i1=i1+1
206   i2=nparticle
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
      go to 201

205   i2=i2-1
      if (vectortmp(i2)+1==vectortmp(i2+1)) then
        go to 202
      else
        go to 203
      endif
      
204   continue
      nbasis=i1-1

      allocate( vectorb(1:nbasis) )
      vectorb=0
      do i1=1,nbasis
        do i2=1,nparticle
          vectorb(i1)=ibset(vectorb(i1),vectori(i2,i1))
        enddo
      enddo


      deallocate(vectortmp,vectori)
      end subroutine


      end module mod_basis
