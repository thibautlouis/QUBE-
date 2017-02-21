!FFLAGS="-fopenmp -fPIC -Ofast -ffree-line-length-none" f2py -c -m mcm_code mcm_code.f90 wigner3j_sub.f -lgomp


subroutine calc_mcm(ml,bl, mcm, mcm_p, mcm_pp, mcm_mm)
  implicit none
  real(8), intent(in)    :: ml(:),bl(:)
  real(8), intent(inout) :: mcm(:,:),mcm_p(:,:),mcm_pp(:,:),mcm_mm(:,:)
  real(8), allocatable   :: thrcof0(:),thrcof1(:)
  real(8), parameter     :: pi = 3.14159265358979323846264d0
  integer :: l1, l2, l3, info, nlmax, lmin, lmax, count
  real(8) :: l1f(2), fac
  nlmax = size(mcm,1)-1
  allocate(thrcof0(2*nlmax+2))
  allocate(thrcof1(2*nlmax+2))

  do l3 = 2, nlmax
    write(*,*) l3
    fac=(2*l3+1)/(4*pi)*bl(l3)**2
    do l2 = 2, nlmax
      ! Get all the 3j symbols for the given l2,l3,m2=0,m3=0.
      call drc3jj(dble(l2),dble(l3),0d0,0d0,l1f(1),l1f(2),thrcof0, size(thrcof0),info)
      call drc3jj(dble(l2),dble(l3),2d0,-2d0,l1f(1),l1f(2),thrcof1, size(thrcof1),info)

      lmin=INT(l1f(1))
      lmax=MIN(nlmax,INT(l1f(2)))
!      mcm(l2-1,l3-1) = (2*l3+1)/(4*pi)*sum(ml(lmin:lmax)*thrcof0(1:lmax-lmin+1)**2d0)

      count=1
      do l1=lmin,lmax
         mcm(l2-1,l3-1) =mcm(l2-1,l3-1)+ fac*(ml(l1)*thrcof0(count)**2d0)
         mcm_p(l2-1,l3-1) =mcm_p(l2-1,l3-1)+ fac*(ml(l1)*thrcof0(count)*thrcof1(count))
         mcm_pp(l2-1,l3-1) =mcm_pp(l2-1,l3-1)+ fac*(ml(l1)*thrcof1(count)**2*(1+(-1)**(l1+l2+l3))/2)
         mcm_mm(l2-1,l3-1) =mcm_mm(l2-1,l3-1)+ fac*(ml(l1)*thrcof1(count)**2*(1-(-1)**(l1+l2+l3))/2)

         count=count+1
      end do

    end do
  end do
end subroutine

subroutine bin_mcm(mcm, bsize, mbb)
  ! Bin the given mode coupling matrix mcm(0:lmax,0:lmax) into
  ! mbb(nbin,nbin) using bins of the given binsize
  implicit none
  real(8), intent(in)    :: mcm(:,:)
  integer, intent(in)    :: bsize
  real(8), intent(inout) :: mbb(:,:)
  integer :: b1, b2, l3, l2, lmax
  lmax = size(mcm,1)-1
  mbb  = 0

  do b2=1,size(mbb,2)
    do b1=1,size(mbb,1)
      ! bin b goes from l=(b-1)*bsize to l=b*bsize-1 inclusive
       do l3=2+(b2-1)*bsize,2+b2*bsize-1
          do l2=2+(b1-1)*bsize,2+b1*bsize-1
             mbb(b1,b2)=mbb(b1,b2) + l2*(l2+1d0)/(l3*(l3+1d0))*mcm(l2-1,l3-1)
          end do
       end do
       mbb(b1,b2) = mbb(b1,b2) / bsize

    end do
  end do
end subroutine
