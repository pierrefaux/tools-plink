program cutHits

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
!Version 2.0, last update on 2020-04-13
!

! This version includes global FDR computation according to
! Benjamini & Hochberg (1995)
! THIS VERSION IS STILL IN DEVELOPMENT ! 
! TODO:
!       - Extend output to blocks below FDR threshold
!	- set alpha as parameter
!	- formatted output/input


implicit none
type mydata
        integer:: pm, nb
        integer*2, allocatable:: mes(:)
        real, allocatable:: pv(:)
end type
type mydata2
        integer:: n
        integer, allocatable:: chr(:), pos(:)
        real, allocatable:: pv(:)
end type
integer:: i, j, k, l, m, n, p, u, v, w, tbp, &
        rstart, rend, nhits, rsl(50000,6), nbm, &
        nbdum, bk
integer, allocatable:: pos(:)
real:: r, rslpv(50000), r2, alpha, fdr5, fdr1
real, allocatable:: pv(:)
character*3:: nbmc
character*9:: tbpc
character*100:: fname
type(mydata):: chr(23)
type(mydata2), allocatable:: measure(:)

!! alpha by default
alpha=0.05

call getarg(1,fname)
call getarg(2,tbpc)
call getarg(3,nbmc)

! input file FNAME contains 4 columns: measurement - chr - posBP - pval
! it can therefore includes measurement from only one or from various
! measurements

! parameter TBPC is the size of the block increment in bp (e.g. 50000)
! parameter NBMC is the maximal number of trait numbers (e.g. 59)

read(tbpc,'(i9)') tbp
open(newunit=u, file=fname)
n=count_lines(u)
print *, '    data in file : ', trim(adjustl(fname))
print *, '       # of rows : ', n
read(nbmc,'(i3)') nbm
print *, 'number of traits : ', nbm
nbdum=nbm+1

! (1) Parsing input file:
!       First the file is read to determine the maximum size of each chromosome,
!       i.e. the maximum position of any SNPs in that chromosome, from chr 1 to
!       23 (24 and other not taken into account so far). Besides, the number of
!       SNPs/measurement is also determined.
!       Then, 23 real vectors and 23 integer vectors are allocated with the
!       maximum size of the chromosome; the real ones are filled in with the
!       value 2.0, corresponding to "no SNP at that position". The input file is
!       sequentially read and the vectors are filled in with p-value and
!       measurement# if the p-value at that SNPs is lower than a previously
!       recorded p-value. Besides, 3 vectors are allocated per measurement (in
!       the variable "measure"); they will store the chr, position and p-values
!       data per measurement.
!       Eventually a few statistics are made per chromosome: minimum p-value,
!       most significative measurement, etc.
allocate(measure(nbdum))
chr(:)%pm=0
measure(:)%n=0
do i=1, n
        read(u,*) j, k, l, r
        if (chr(k)%pm<l) chr(k)%pm=l
        measure(j)%n=measure(j)%n+1
enddo
print *, '    -- Maximum positions per chromosome :'
do i=1, 23
        print *, '        ', i , '   ', chr(i)%pm
        allocate(chr(i)%pv(chr(i)%pm), chr(i)%mes(chr(i)%pm))
        chr(i)%pv=2.0
enddo
do i=1, nbm
        allocate(measure(i)%chr(measure(i)%n), measure(i)%pos(measure(i)%n), &
                measure(i)%pv(measure(i)%n))
enddo
rewind(u)
measure(:)%n=0
do i=1, n
        read(u,*) j, k, l, r
        if (chr(k)%pv(l)>r) then
                chr(k)%pv(l)=r
                chr(k)%mes(l)=j
        endif
        measure(j)%n=measure(j)%n+1
        measure(j)%chr(measure(j)%n)=k
        measure(j)%pos(measure(j)%n)=l
        measure(j)%pv(measure(j)%n)=r
enddo
j=measure(1)%n
allocate(measure(nbdum)%chr(j), measure(nbdum)%pos(j), measure(nbdum)%pv(j))
measure(nbdum)%chr=measure(1)%chr
measure(nbdum)%pos=measure(1)%pos
measure(nbdum)%pv=measure(1)%pv

! (1b) FDR thresholds computation
!       The FDR threshold overall measures is computed first by gathering all
!       p-values into one single vector, then per measurement
!
allocate(pv(n))
k=0
do i=1, nbm
        pv(k+1:k+measure(i)%n)=measure(i)%pv
        k=k+measure(i)%n
enddo
print *, ' FDR threshold over all measures (alpha 5%, 1%) = ', fdr_threshold(pv,alpha,n), &
        fdr_threshold(pv,0.01,n)
deallocate(pv)
do i=1, nbm
        fdr5=fdr_threshold(measure(i)%pv,alpha,measure(i)%n)
        fdr1=fdr_threshold(measure(i)%pv,0.01,measure(i)%n)
        print *, i, fdr5, fdr1
enddo


stop

print *, '    -- Minimum p-value, measurement, number of SNPs below significance '
print *, '       threshold'
do i=1, 23
        print *, '        ', i , '   ', minval(chr(i)%pv), '   ', &
                count(chr(i)%pv<=1E-5)
enddo 

! (2) Defining the blocks
!       For each chromosome, the real vector of p-values is searched for its
!       lowest value: if that value is smaller than the threshold of 1E-5, a
!       region is defined around that position by incrementing down- and upwards
!       blocks of the specified length in TBPC if the block contain at least one
!       SNP with value lower than 1E-5.
!       Eventually, the incrementation stops if (1) we get to chromosome edges
!       or (2) if the block to increment does not include any SNP with p-v<1E-5
nhits=0
do i=1, 23
        l=nhits
        allocate(pv(chr(i)%pm))
        pv=chr(i)%pv
        r=minval(pv)
        if (r>1E-5) then
                print '(a,i3,a)', '    chr ', i, ' (no SNP with p-value lower than 1E-5)'
                deallocate(pv)
                cycle
        else
                do
                        r=minval(pv)
                        if (r>1E-5) exit
                        nhits=nhits+1
                        rslpv(nhits)=r
                        rstart=minloc(pv,1)
                        rsl(nhits,1)=chr(i)%mes(rstart)
                        rsl(nhits,2)=rstart
                        rsl(nhits,6)=i
                        rend=rstart+tbp
                        if (rend>chr(i)%pm) rend=chr(i)%pm
                        rstart=rstart-tbp
                        if (rstart<1) rstart=1
                        ! upwards
                        do
                                r=minval(pv(rend-tbp:rend))
                                if (r>1E-5) then
                                        rend=rend-tbp
                                        exit
                                endif
                                rend=rend+tbp
                                if (rend>chr(i)%pm) then
                                        rend=chr(i)%pm
                                        exit
                                endif
                        enddo
                        !downwards
                        do
                                r=minval(pv(rstart:rstart+tbp))
                                if (r>1E-5) then
                                        rstart=rstart+tbp
                                        exit
                                endif
                                rstart=rstart-tbp
                                if (rstart<1) then
                                        rstart=1
                                        exit
                                endif
                        enddo
                        pv(rstart:rend)=2.0
                        rsl(nhits,3:4)=(/rstart,rend/)
                        rsl(nhits,5)=rend-rstart+1
                enddo
                deallocate(pv)
                print '(a,i3,a,i5,a)', '    chr ', i, ' found ', nhits-l, ' blocks :'
                do j=l+1, nhits
                        print '(a,i3,4i10,e12.4)', '             ', rsl(j,1:5), rslpv(j)
                enddo
        endif
enddo
print *, '   TOTAL : ', nhits , ' blocks found, average size = ', &
        1E-6*real(sum(rsl(1:nhits,5)))/real(nhits), ' Mb'

! (3) Output block details
!       For each block found, the collected information (leading SNP chromosome
!       and position, block extent, block length, measurement and p-value
!       of the leading SNP) is outputted along with a criterion that is the
!       maximal distance between the leading SNP and any other significative SNP
!       falling in the block extent. The measurement and p-value of that other
!       SNP are also outputted.
open(newunit=v, file='blocks_details')
do i=1, nhits
        l=rsl(i,2)
        rstart=rsl(i,3)
        rend=rsl(i,4)
        r=rslpv(i)
        allocate(pv(rend))
        pv=2.0
        pv(rstart:rend)=chr(rsl(i,6))%pv(rstart:rend)
        m=0
        do k=rstart, rend
                if (pv(k)>5E-8) cycle
                if (abs(k-l)>m) then
                        r2=pv(k)
                        p=chr(rsl(i,6))%mes(k)
                        m=abs(k-l)
                        bk=k
                endif
        enddo
        if (m==0) then
                write(v,'(i3,i4,4i10,e12.4,i10,i4,e12.4,i10)') rsl(i,6),rsl(i,1:5),rslpv(i),0,0,0,0
        else
                write(v,'(i3,i4,4i10,e12.4,i10,i4,e12.4,i10)') rsl(i,6),rsl(i,1:5),rslpv(i),m,p,r2,bk
        endif
        deallocate(pv)
enddo
close(v)

! (4) Output measurement assignation and p-values along the genome
!       First reset the chr(.)%mes vectors to match with the measurements
!       assigned in blocks and set measurement NBDUM (= arbitrary, = 1 so far) to
!       any other positions.
!       Then, run along the genome to find the p-values corresponding to the
!       measurement corresponding to a given block.
!       Eventually, output a file with this data: chr - position - measurement -
!       p-value

do i=1, 23
        chr(i)%mes=nbdum
        chr(i)%pv=-1.0
        j=count(measure(nbdum)%chr==i)
        allocate(pv(j), pos(j))
        pv=pack(measure(nbdum)%pv,measure(nbdum)%chr==i)
        pos=pack(measure(nbdum)%pos,measure(nbdum)%chr==i)
        chr(i)%pv(pos)=pv
        deallocate(pv, pos)
enddo
do i=1, nhits
        j=rsl(i,1) ! measurement
        k=rsl(i,6) ! chromosome
        rstart=rsl(i,3)
        rend=rsl(i,4)
        chr(k)%mes(rstart:rend)=j
        l=count(measure(j)%chr==k.and.measure(j)%pos>=rstart.and.measure(j)%pos<=rend)
        allocate(pv(l), pos(l))
        pv=pack(measure(j)%pv,measure(j)%chr==k.and.measure(j)%pos>=rstart.and.measure(j)%pos<=rend)
        pos=pack(measure(j)%pos,measure(j)%chr==k.and.measure(j)%pos>=rstart.and.measure(j)%pos<=rend)
        chr(k)%pv(pos)=pv
        deallocate(pv,pos)
enddo
open(newunit=v, file='pvalues_accross_measurements')
do i=1, 23
        do j=1, chr(i)%pm
                if (chr(i)%pv(j)<0.0) cycle
                write(v,'(i3,i10,i4,e12.4)') i, j, chr(i)%mes(j), chr(i)%pv(j)
        enddo
enddo
close(v)


contains

function fdr_threshold(pvs0,alpha,ntest) result(t)
! This function implements Benjamini-Hochberg procedure
implicit none
integer:: i, j, k, ntest
real:: pvs0(ntest), alpha, t, p, mpvs
real, allocatable:: pvs(:)

! reducing the size of the pv vector in order to improve fetching speed
k=count(pvs0<.05)
allocate(pvs(k))
pvs=pack(pvs0,pvs0<.05)
! TODO: a condition should here make sure that some elements remain in the pvs
! vector

! print *, '    (will compute FDR threshold on the ', k, ' lowest p-values)'

! computing the FDR threshold
t=0.0
p=alpha/real(ntest)
mpvs=minval(pvs)
if (mpvs<p) then
!        print *, '    ', 1, p, mpvs
        ! keep iterating from 2nd-lowest p-value
        j=minloc(pvs,1)
        pvs(j)=1.0
        do i=2, k
                t=mpvs
                mpvs=minval(pvs)
                p=real(i)*alpha/real(ntest)
!                print *, '    ', i, p, mpvs
                if (mpvs>p) then
                        t=p
                        exit
                endif
                j=minloc(pvs,1)
                pvs(j)=1.0
        enddo
endif
end function fdr_threshold

function count_lines(un_file) result(n)
implicit none
integer:: un_file, n, k
character*1:: junk
rewind(un_file); n=0
do
        read(un_file,'(a1)',iostat=k) junk
        if (k/=0) exit
        n=n+1
enddo
rewind(un_file)
end function count_lines

end program
