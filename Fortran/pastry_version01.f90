!##########################################
!############ Fortran code version 01 for PASRY and LS by EUNJI KIM 
!############ Read summary statistics(beta,se) and apply PASTRY and LS approaches
!############ This version works only for the same number of cases and controls between the studies
!############ Main code starts from the line containing "program assoc_pastry"
!##########################################



MODULE equations
contains
   ! Returns the inverse of a matrix calculated by finding the LU
   ! decomposition.  Depends on LAPACK.
   function inv(A,nstudy) result(Ainv)
      integer,intent(in) :: nstudy
      DOUBLE PRECISION, dimension(1:nstudy,1:nstudy), intent(in) :: A
      DOUBLE PRECISION, dimension(1:nstudy,1:nstudy) :: Ainv
      DOUBLE PRECISION, dimension(1:nstudy) :: work  ! work array for LAPACK
      integer, dimension(1:nstudy) :: ipiv   ! pivot indices
      integer :: n, info, i,j
      ! External procedures defined in LAPACK
      external DGETRF
      external DGETRI
      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      Ainv = A
      n = size(A,1)
      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call DGETRF(nstudy, nstudy, Ainv, nstudy, ipiv, info)
      !call DGETRF(n, n, Ainv, n, ipiv, info)
      if (info /= 0) then
         stop 'Matrix is numerically singular!'
      end if
      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call DGETRI(nstudy, Ainv, nstudy, ipiv, work, nstudy, info)
      !call DGETRI(n, Ainv, n, ipiv, work, n, info)
      if (info /= 0) then
          stop 'Matrix inversion failed!'
      end if
   end function inv
   
   ! Correlation
   double precision function pearsonr(x, y) result(r)
   ! given two arrays x and y, this function computes their Pearson correlation coefficient r
      implicit none
      double precision, dimension(:) :: x, y
      double precision, dimension(size(x)) :: xt, yt
      double precision :: ax, ay, df, sxx, sxy, syy
      integer :: n
     
      if (size(x) /= size(y)) STOP 'Dimension mismatch in pearsonr'
      n = size(x)
      ! find the means and subtract them
      ax = sum(x)/n;ay = sum(y)/n
      xt = x - ax;yt = y - ay
      sxx = dot_product(xt,xt);syy = dot_product(yt,yt)
      sxy = dot_product(xt,yt)
      r = sxy / sqrt(sxx*syy)
      ! abs(r) cannot be > 1, except for artifacts of floating point arithmetic
      r = max(min(r, 1.d0), -1.d0) 
   end function pearsonr

END MODULE equations

!read file
module strmod
contains
   subroutine parse ( line, words, nw )
      implicit none
      character(*), intent(in)  :: line
      character(*), intent(out) :: words(:)
      integer,      intent(out) :: nw
      character(len(words)) :: buf( size(words) )
      integer :: k, ios
      nw = 0 ; words(:) = ""
      do k = 1, len(line)
          read( line, *, iostat=ios ) buf( 1 : k )
          if ( ios /= 0 ) exit
          nw = k
          words( 1: nw ) = buf( 1 : nw )
      enddo
   endsubroutine
endmodule


!!!!! main program 
program assoc_pastry
    use equations
    use strmod
   
    implicit none
    character*10000 :: line,firstline
    character*10000 :: casemaf_file,contmaf_file
    character*10000, allocatable :: data_file(:)
    integer :: endline,i,j,k,l,m,n,stat
    integer :: nshared,ncase,ncont
    integer :: novalue,nstudy
    double precision :: beta_LS,zstat_LS,meta_beta_LS,meta_se_LS
    double precision :: zstat_pastry,meta_beta_pastry,meta_se_pastry
    double precision :: total,controlpart,denomat,summat,nummat
    double precision :: casemaf,contmaf
    double precision  :: exp_alpha_i,exp_alpha_j,C,x,x0i,x0j,x1i,x1j
    double precision  :: mat21,mat22,pastrycor
    double precision  :: x_0(1:2,1:2),x_1(1:2,1:2),info_i(1:2,1:2),info_j(1:2,1:2)
    double precision  :: cov00(1:2,1:2),cov11(1:2,1:2),covmatrix(1:2,1:2)
    double precision  :: invi(1:2,1:2),invj(1:2,1:2),covfinal(1:2,1:2),midfinal(1:2,1:2)
    double precision, allocatable :: R(:,:),invR(:,:),pastryR(:,:),pastryinvR(:,:)
    double precision, allocatable :: mat(:,:),midmat(:,:),rowsum(:),betamat(:)
    double precision, allocatable :: comatrix_lin(:,:),invcomatrix_lin(:,:)
    double precision, allocatable :: beta(:),se(:),stder_LS(:),stder_pastry(:)
    
    call getarg(1,line)
    read(line,'(i10)') nstudy
    call getarg(2,line)
    read(line,'(i10)') endline
    call getarg(3,line)
    read(line,'(i10)') ncase
    call getarg(4,line)
    read(line,'(i10)') ncont
    call getarg(5,line)
    read(line,'(i10)') nshared
    call getarg(6,line)
    read(line,'(A1000)') casemaf_file
    call getarg(7,line)
    read(line,'(A1000)') contmaf_file

    allocate(R(nstudy,nstudy),invR(nstudy,nstudy))
    allocate(rowsum(nstudy),mat(nstudy,nstudy))
    allocate(betamat(nstudy),stder_LS(nstudy))
    allocate(comatrix_lin(nstudy,nstudy))
    allocate(invcomatrix_lin(nstudy,nstudy))
    allocate(data_file(nstudy))
    do i=1,nstudy
      call getarg(i+7,line)
      read(line,'(A5000)') data_file(i)
    end do
    allocate(beta(nstudy),se(nstudy))
    allocate(stder_pastry(nstudy))
    allocate(pastryR(nstudy,nstudy),pastryinvR(nstudy,nstudy))
    stat=0
    !initialization
    do i=1,nstudy
       do j=1,nstudy
          R(i,j)=1.0
          invR(i,j)=1.0
       enddo
    enddo
    !correlation from LS 
    do i=1,(nstudy-1)
       do j=(i+1),nstudy
             total=sqrt((real(ncase+ncont+nshared)*real(ncase+ncont+nshared)))
             controlpart=real(nshared)*sqrt(real(ncase*ncase)/(real(ncont+nshared)*real(ncont+nshared)))
             R(i,j)=real(controlpart/total)
             R(j,i)=real(controlpart/total)
       enddo
    enddo
    invR=inv(R,nstudy)
    !To skip first line
    do  i=1,nstudy
       open(i,FILE= data_file(i)) 
       read(i,*,iostat=stat) firstline 
    enddo
    open(1000,FILE= casemaf_file) 
    read(1000,*,iostat=stat) firstline 
    open(1001,FILE= contmaf_file) 
    read(1001,*,iostat=stat) firstline 

    ! output file
    open(500,FILE="Zstat_LS.txt") 
    open(700,FILE="Zstat_PASTRY.txt") 
    open(1000,FILE= casemaf_file) 
    open(1001,FILE= contmaf_file) 
    write(500,*) "zstat"
    write(700,*) "zstat"
    
    do k=1,(endline-1)
       novalue=0
       read(1000,*,iostat=stat) casemaf
       read(1001,*,iostat=stat) contmaf
       do i=1,nstudy
          read(i,*,iostat=stat) beta(i),se(i)
          !pass the line with casemaf or contmaf=0 or beta or se value with NA (from summary statistics)
          if((casemaf==0).OR.(contmaf==0).OR.(se(i)==0.00000) .OR. (beta(i)==0.00000))then
             novalue=1
          endif
       enddo 
    
       if(novalue==0)then
          mat=0;rowsum=0
          nummat=0;betamat=0
          do i=1,nstudy
             summat=0
             do j=1,nstudy
                mat(i,j)=(1/se(i))*invR(i,j)*(1/se(j))
                summat=summat+mat(i,j)
             enddo
             rowsum(i)=summat
             betamat(i)=summat*beta(i)
          enddo
         
          do i=1,nstudy
             stder_LS(i)=sqrt(1/rowsum(i))
          enddo 
          beta_LS=(sum(betamat))/(sum(rowsum))
          meta_beta_LS=beta_LS
          meta_se_LS=sqrt(1.00/(sum(rowsum)))
          !PASTRY
          info_i=0;info_j=0;invi=0;invj=0
          x_0=0;x_0(1,1)=1;x_1=1
          x=real((real(ncase)*casemaf+real(nshared+ncont)*contmaf)/(real(ncase)+nshared+ncont))
          exp_alpha_i=exp(log(real(ncase)/real(nshared+ncont))-beta_LS*x)
          x=real((real(ncase)*casemaf+real(nshared+ncont)*contmaf)/(real(ncase)+nshared+ncont))
          exp_alpha_j=exp(log(real(ncase)/real(nshared+ncont))-beta_LS*x)
          x0i=exp_alpha_i;x0j=exp_alpha_j
          x1i=exp_alpha_i*exp(real(beta_LS));x1j=exp_alpha_j*exp(real(beta_LS))
    
          do l=1,2
             do m=1,2
               C=real((ncont+nshared))*((1-contmaf)*(x0i/(1+x0i)**2)*x_0(l,m) + contmaf*(x1i/(1+x1i)**2)*x_1(l,m)) 
               info_i(l,m)=real(ncase)*((1-casemaf)*(x0i/(1+x0i)**2)*x_0(l,m)+casemaf*(x1i/(1+x1i)**2)*x_1(l,m))+C
               C=real((ncont+nshared))*((1-contmaf)*(x0j/(1+x0j)**2)*x_0(l,m) + contmaf*(x1j/(1+x1j)**2)*x_1(l,m))
               info_j(l,m)=real(ncase)*((1-casemaf)*(x0j/(1+x0j)**2)*x_0(l,m)+casemaf*(x1j/(1+x1j)**2)*x_1(l,m))+C
               cov00(l,m)=(1-contmaf)*(real(x0i*x0j)/real((1+x0i)*(1+x0j)))*x_0(l,m)
               cov11(l,m)=contmaf*(real(x1i*x1j)/real((1+x1i)*(1+x1j)))*x_1(l,m)
               covmatrix(l,m)=real(nshared)*real(cov00(l,m)+cov11(l,m))
             enddo
          enddo
          invi=inv(info_i,2);invj=inv(info_j,2)
          mat21=invi(2,1)*covmatrix(1,1)+invi(2,2)*covmatrix(2,1)
          mat22=invi(2,1)*covmatrix(1,2)+invi(2,2)*covmatrix(2,2)
          covfinal(2,2)=mat21*invj(1,2)+mat22*invj(2,2)
          pastrycor=covfinal(2,2)/sqrt(real(invi(2,2)*invj(2,2)))
          
          do i=1,nstudy
             do j=1,nstudy
                pastryR(i,j)=1.0
                pastryinvR(j,i)=1.0
             enddo
          enddo
          
          do i=1,(nstudy-1) 
             do j=(i+1),nstudy
               pastryR(i,j)=pastrycor
               pastryR(j,i)=pastrycor
             enddo
          enddo
          print *, "pastrycor"
          print *, pastrycor
          pastryinvR=inv(pastryR,nstudy)
          mat=0;rowsum=0;betamat=0
          do i=1,nstudy
             summat=0
             do j=1,nstudy
                mat(i,j)=(1/se(i))*pastryinvR(i,j)*(1/se(j))
                summat=summat+mat(i,j)
             enddo
             rowsum(i)=summat
             betamat(i)=summat*beta(i)
          enddo

         do i=1,nstudy
            stder_pastry(i)=sqrt(1/rowsum(i))
         enddo
         meta_beta_pastry=(sum(betamat))/(sum(rowsum))
         meta_se_pastry=sqrt(1.00/(sum(rowsum)))
    endif
  
    if(novalue==1)then
       meta_beta_LS=0;meta_se_LS=0;zstat_LS=0
       meta_beta_pastry=0;meta_se_pastry=0;zstat_pastry=0
    else
       zstat_LS=meta_beta_LS/meta_se_LS
       zstat_pastry=meta_beta_pastry/meta_se_pastry
    endif
  
       write(500,*)  zstat_LS
       write(700,*)  zstat_pastry
  end do
    do i=1,nstudy
       close (i)
    enddo
 
  close(500)
  close(700)
end program assoc_pastry
