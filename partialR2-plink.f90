program partialR2

!
!This program is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!any later version.

!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.

!You should have received a copy of the GNU General Public License
!along with this program.  If not, see <http://www.gnu.org/licenses/>.

!Pierre FAUX (pierrefaux@gmail.com), 2020
!Version 0.2, last update on 2020-12-10
!
!Additional credits:
!       Subroutine ginv2 was borrowed from Ignacy Misztal's
!       blupf90 collection of programs.
!       See http://nce.ads.uga.edu/html/projects/programs/
!
implicit none
integer:: i, j, k, l, m, n, u, v, nc, m0, nsnps
integer, allocatable:: w1(:), w2(:,:), pos(:), chr(:), nmis(:), &
                kept(:), miss(:)
integer*1, allocatable:: gi1(:,:)
real:: r
real, allocatable:: x(:,:), y(:), xy(:), xx(:,:), &
        ixx(:,:,:), ixy(:,:), xxi(:,:), xi(:), &
        yr(:), yr0(:), z(:), zr(:), b(:), xz(:), &
        pr2(:), xx1(:,:), xy1(:)
character*5:: ncc
character*20:: fmt1
character*100:: file1, file2, file3, junk1, junk2
character*1000:: junkrs
character*1000, allocatable:: rs(:)

! ASSUMPTIONS:
!       - No missing phenotypes
!       - No missing covariates
!       - Phenotypes and covariates files to follow sample's order of FAM

! TODO:
!	- Check how good are approximations considering max 5% missing, 
!       without updating xxi and yr
!	- OMP parallelization
!	- Integrate missing phenotypes/covariates


call getarg(1,file1)    ! prefix of BED-BIM-FAM
call getarg(2,file2)    ! phenotype file (no header; 3 columns: FID, IID, pheno)
call getarg(3,file3)    ! covariates file (no header; n columns)
call getarg(4,ncc)      ! number of covariates (excluding sex)

! (1) Load files and compute basic matrices
!       * load X (x) and y (y) for further use
!       * Compute X'.X (xx) and X'.y (xy)
!       * Compute individual contributions to X'.X (ixx) and X'.y (ixy)
!       * Invert X'.X (xxi)
open(newunit=u, file=file1(:len_trim(file1))//'.fam')
m=count_lines(u)
!nc=7
read(ncc,'(i5)') nc
n=nc+2
print *, '   #samples =', m
print *, '#covariates =', n
allocate(x(m,n), y(m), xy(n), xx(n,n), ixx(n,n,m), &
        ixy(n,m), xxi(n,n), xi(n), xx1(n,n), xy1(n))
x=0.d0
do i=1, m
        read(u,*) junk1, junk2, j, k, l
        x(i,l)=1
enddo
close(u)

open(newunit=u, file=file2)
open(newunit=v, file=file3)
xx=0.d0
xy=0.d0
do i=1, m
        xi(1:2)=x(i,1:2)
        read(u,*) junk1, junk2, y(i)
        read(v,*) junk1, junk2, xi(3:n)
        x(i,:)=xi
        do j=1, n
                ixx(j,:,i)=xi(j)
        enddo
        do j=1, n
                ixx(:,j,i)=ixx(:,j,i)*xi(j)
        enddo
        ixy(:,i)=xi*y(i)
        xx=xx+ixx(:,:,i)
        xy=xy+ixy(:,i)
enddo
close(u)
close(v)
xxi=xx
call ginv2(xxi,n)

open(newunit=u, file=file1(:len_trim(file1))//'.bim')
nsnps=count_lines(u)
m0=ceiling(real(m)*.25)*4
print *, ' dims=', m, ' by ', nsnps
print *, ' actual #rows = ', m0
allocate(gi1(m0,nsnps),chr(nsnps),pos(nsnps), rs(nsnps))
do i=1, nsnps
        read(u,*) chr(i), rs(i), r, pos(i)
enddo
close(u)
call readBed2(file1,m0,m,nsnps,gi1)
call translateBed(gi1,m0,nsnps)

! (2) Compute linear regression of X on y and get yr, where
!       * X is the m-by-n matrix of covariates (from file3)
!       * y is the m-by-1 vector of phenotypes (from file2)
!       * yr is the m-by-1 vector of residuals: yr = y - X.b
!       * b is the n-by-1 vector of effect estimates, obtained  by OLS
allocate(yr(m),yr0(m),b(n))
yr=y
do i=1, n
        b(i)=dot_product(xxi(i,:),xy)
        yr=yr-x(:,i)*b(i)
enddo
yr0=yr

! (3) For each SNP in the BED, compute linear regression of
!       X on z and get zr, where
!       * X is the m-by-n matrix of covariates (from file3)
!       * z is the m-by-1 vector of genotypes (from BED)
!       * zr is the m-by-1 vector of residuals: zr = z - X.b
!       * b is the n-by-1 vector of effect estimates, obtained  by OLS
!       ... starting by SNPs with no missing values and then others ...

allocate(xz(n), pr2(nsnps), z(m), nmis(nsnps), zr(m))
pr2=0.d0
do i=1, nsnps
        nmis(i)=count(gi1(1:m,i)==9)
        if (nmis(i)>0) cycle
        z=real(gi1(1:m,i))
        xz=0
        do j=1, n
            xz(j)=dot_product(x(:,j),z)
        enddo
        zr=z    
        do j=1, n
                b(j)=dot_product(xxi(j,:),xz)
                zr=zr-x(:,j)*b(j)
        enddo
        pr2(i)=sqcorr(yr,zr,m)
enddo
deallocate(z,zr,yr)
!       ... for SNPs with NAs, the computation sequence is:
!               * Get missing/kept samples addresses
!               * Subtract their ixx(:,:,i) from original xx (xx1)
!               * Invert the update xx1 --> xxi
!               * Subtract their ixy(:,i) from original xy (xy1)
!               * Get b=xxi*xy1 and yr=y-Xb
!               * Get xz=Xz, b=xxi*xz and zr=z-Xb
!               * Get pr2(i)
allocate(w1(m))
w1=(/(i,i=1,m)/)
do i=1, nsnps
        if (nmis(i)==0) cycle
        if (nmis(i)==m) then
                pr2(i)=0
                cycle
        endif
        l=m-nmis(i)
        allocate(kept(l), miss(nmis(i)), yr(l), zr(l), z(l))
        kept=pack(w1,gi1(:,i)<9)
        miss=pack(w1,gi1(:,i)==9)
        xx1=xx
        xy1=xy
        do j=1, nmis(i)
                xx1=xx1-ixx(:,:,miss(j))
                xy1=xy1-ixy(:,miss(j))
        enddo
        xxi=xx1
        call ginv2(xxi,n)
        yr=y(kept)
        do j=1, n
                b(j)=dot_product(xxi(j,:),xy1)
                yr=yr-x(kept,j)*b(j)
        enddo
        z=real(gi1(kept,i))
        xz=0
        do j=1, n
                xz(j)=dot_product(x(kept,j),z)
        enddo
        zr=z
        do j=1, n
                b(j)=dot_product(xxi(j,:),xz)
                zr=zr-x(kept,j)*b(j)
        enddo
        pr2(i)=sqcorr(yr,zr,l)
        deallocate(kept, miss, yr, zr, z)
enddo

! (4) Output results to file partial_r2.txt
!       ($1: chr, $2: pos, $3: rs, $4: #samples, $5: partial R2)
print *, 'Results are written to file partial_r2.txt:'
print *, '     $1:      CHR'
print *, '     $2:      POS'
print *, '     $3:      RSID'
print *, '     $4:      #SAMPLES'
print *, '     $5:      PARTIAL_R2'
open(newunit=u, file='partial_r2.txt')
do i=1, nsnps
        junkrs=rs(i)
        write(u,'(i2,i10,a1,a,i7,e12.5)') chr(i), pos(i), ' ', &
                junkrs(:len_trim(junkrs)), m-nmis(i), pr2(i)
enddo
close(u)

contains

function count_lines(un_file) result(n)
implicit none
integer:: un_file, n, k
character*1:: junk
rewind(un_file); n=0
do
        read(un_file,'(a1)',iostat=k) junk
        if(k/=0)exit
        n=n+1
enddo
rewind(un_file)
end function count_lines


subroutine ginv2(a,n)
! returns generalized inverse of x(n,n). tol is working zero 
! and irank returns the rank of the matrix. rework of rohan fernando's 
! f77 subroutine by i. misztal 05/05/87-05/23/00
implicit none
real :: a(n,n),w(n),t1,t2, tol
integer  :: ii, i, n, j, rank
tol=1d-12
rank=n
do i=1,n
        do  j=1,i-1
                a(i:n,i)=a(i:n,i)-a(i,j)*a(i:n,j)
        enddo
        if (abs(a(i,i)).lt.tol) then
                a(i:n,i)=0.0
                rank=rank-1
        elseif( a(i,i) < 0.d0) then
                print*,'ginv2: Negative pivot:',i,a(i,i)
                stop
        else
                a(i,i)=sqrt(a(i,i))
                a(i+1:n,i)=a(i+1:n,i)/a(i,i)
        endif
enddo
do i=1,n
        if (a(i,i).eq.0.) then
                a(i+1:n,i)=0
        else
                a(i,i)=1.0/ a(i,i)
                w(i+1:n)=0
                do  ii=i+1,n
                        w(ii:n)=w(ii:n)-a(ii:n,ii-1)*a(ii-1,i)
                        if (a(ii,ii).eq.0.) then
                                a(ii,i)=0.
                        else
                                a(ii,i)=w(ii)/a(ii,ii)
                        endif
                enddo
        endif
enddo
do j=1,n
        do i=j,n
                a(i,j)=dot_product(a(i:n,j),a(i:n,i))
        enddo
enddo
do i=1,n
        a(i,i+1:n)=a(i+1:n,i)
enddo
end subroutine ginv2

subroutine readBed2(bf, n0, n, m, g)
implicit none
integer:: i, j, k, l, m, n, u, v, n0
integer*1:: junk(3), g(n0, m)
integer*1, allocatable:: bed(:,:)
character*100:: bf

l=ceiling(real(n)*.25)
!print *, 'nb integers1=', l
allocate(bed(l,m))
open(newunit=u, file=bf(:len_trim(bf))//'.bed',form='unformatted',access='stream')
read(u) junk
!print *, 'starting values =', junk
do i=1, m
        read(u) bed(:,i)
enddo
close(u)
k=0
do i=1, l
        do j=0, 6, 2
                k=k+1
                g(k,:)=ibits(bed(i,:),j,2)
        enddo
enddo
deallocate(bed)
end subroutine readBed2

subroutine translateBed(g, n0, m)
implicit none
integer:: i, j, k, n0, m
integer*1:: g(n0,m)
do i=1, n0
        forall(j=1:m, g(i,j)==1) g(i,j)=9
        forall(j=1:m, g(i,j)==2) g(i,j)=1
        forall(j=1:m, g(i,j)==0) g(i,j)=2
        forall(j=1:m, g(i,j)==3) g(i,j)=0
enddo
end subroutine translateBed

function sqcorr(a,b,n) result(r2)
implicit none
integer:: i, j, k, n
real:: a(n), b(n), r2, ma, mb
ma=sum(a)/real(n)
mb=sum(b)/real(n)
a=a-ma
b=b-mb
r2=(dot_product(a,b))**2/(dot_product(a,a)*dot_product(b,b))
end function sqcorr

end program
