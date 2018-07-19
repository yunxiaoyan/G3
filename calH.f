      include 'basis.f'

      module mod_Hmatrix
      save
      real*8,pointer :: H(:,:)


      contains
      ! You have to check difference between left and right vectors
      ! before using this function
      ! do i1=1,nbasis
      !   vtmp1=vectorb(i1)  ! b1
      ! do i2=1,nbasis 
      !   vtmp2=vectorb(i2)  ! b2
      !   vtmp3=iand(vectorb(i1),vectorb(i2))  ! And(b1,b2)
      !   i7=popcnt(vtmp3)
      !   if (ltmp<nparticle-2) go to 201
      function tbdm(b1,b2,p,q,s,r)
      use mod_sps
      implicit none
      integer :: i1
      integer :: tbdm,p,q,r,s
      integer*8 :: b1,b2

      tbdm=0
      if (.not.btest(b2,r)) go to 201
      i1=poppar(ibits(b2,1,r-1))
      b2=b2-iand(b2,sporbitb(r))
      if (.not.btest(b2,s)) go to 201
      i1=ieor( i1,poppar(ibits(b2,1,s-1)) )
      b2=b2-iand(b2,sporbitb(s))

      if (btest(b2,q)) go to 201
      i1=ieor( i1,poppar(ibits(b2,1,q-1)) )
      b2=b2+sporbitb(q)
      if (btest(b2,p)) go to 201
      i1=ieor( i1,poppar(ibits(b2,1,p-1)) )
      b2=b2+sporbitb(p)

      if (b1/=b2) go to 201
      if (btest(i1,0)) then
        tbdm=-1
      else
        tbdm=1
      endif

201   continue
      end function


      function obdm(b1,b2,p,q)
      use mod_sps
      implicit none
      integer :: i1
      integer :: obdm,p,q
      integer*8 :: b1,b2

      obdm=0
      if (.not.btest(b2,q)) go to 202
      i1=poppar(ibits(b2,1,q-1))
      b2=b2-iand(b2,sporbitb(q))

      if (btest(b2,p)) go to 202
      i1=ieor( i1,poppar(ibits(b2,1,p-1)) )
      b2=b2+sporbitb(p)

      if (b1/=b2) go to 202
      if (btest(i1,0)) then
        obdm=-1
      else
        obdm=1
      endif

202   continue
      end function


      subroutine calHmatrix_SM()
      use mod_sps
      use mod_basis
      implicit none
      integer :: i0,i1,i2,i3,i4,i5,i6,i7
      integer :: i11,i12,i13,i14
      real*8 :: a1,a2
      integer :: ntbme
      integer :: bitdiff(1:nbasis,1:nbasis)
      integer*8 :: vtmp1,vtmp2,vtmp3
      real*8 :: aspe(1:nsps)
      integer,pointer :: a(:),b(:),c(:),d(:)
      real*8,pointer :: atbme(:)

      bitdiff=999999
      do i1=1,nbasis
        do i2=i1,nbasis 
          vtmp1=iand(vectorb(i1),vectorb(i2))  ! And(b1,b2)
          i7=popcnt(vtmp1)
          bitdiff(i1,i2)=nparticle-i7
        enddo
      enddo

      allocate(H(1:nbasis,1:nbasis))
      H=0d0
      open(unit=101,file='TBME')
      read(101,*)
      read(101,*)
      read(101,*)
      read(101,*) ntbme,aspe(:)
      read(101,*)
      do i1=1,nbasis
        vtmp1=vectorb(i1)
        do i0=1,nsps
          H(i1,i1)=H(i1,i1)+obdm(vtmp1,vtmp1,i0,i0)*aspe(i0)
        enddo
      enddo

      allocate( a(1:ntbme),b(1:ntbme),c(1:ntbme),d(1:ntbme)
     $,atbme(1:ntbme) )
      a=0
      b=0
      c=0
      d=0
      atbme=0d0
      do i0=1,ntbme
        read(101,*) a(i0),b(i0),c(i0),d(i0),atbme(i0)
      enddo
      do i1=1,nbasis
        vtmp1=vectorb(i1)  ! b1
        do i2=i1,nbasis 
          vtmp2=vectorb(i2)  ! b2
          if (bitdiff(i1,i2)>2) go to 203
          do i0=1,ntbme
            H(i1,i2)=H(i1,i2)
     $+tbdm(vtmp1,vtmp2,a(i0),b(i0),d(i0),c(i0))*atbme(i0)
     $+tbdm(vtmp1,vtmp2,c(i0),d(i0),b(i0),a(i0))*atbme(i0)
          enddo
203     enddo
      enddo

      close(101)
      deallocate(a,b,c,d,atbme)
      end subroutine


      end module mod_Hmatrix
