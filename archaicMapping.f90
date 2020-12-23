program archaicMapping

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
!Version 0, last update on 2020-05-15
!

! THIS VERSION IS IN PRE-DEVELOPMENT ! 
! Next parts to include: 
!	- association testing
!	- output archaic-coded BED files
!	- OMP paralellization 

implicit none
integer:: i, j, k, l, m, n, u, v, npos, nind, from, to, chr
integer, allocatable:: w1(:), w2(:,:), w3(:), w5(:), pos(:), gc(:,:), &
        inds(:)
real:: a, r
character*1:: tab, junk(1000)
character*1, allocatable:: gcped(:)
character*2:: chrc
character*20:: fmt1
character*100:: file1, file2, file3
character*1000:: junk1000

tab=char(9)

! admixtureHMM input file: npos+1 lines
call getarg(1,file1)
file1=trim(adjustl(file1))
! admixtureHMM output file: nind lines
call getarg(2,file2)
file2=trim(adjustl(file2))
call getarg(3,chrc)
read(chrc,'(i2)') chr

open(newunit=u, file=file1)
open(newunit=v, file=file2)
npos=count_lines(u)-1
nind=count_lines(v)/2
print *, '    # of markers = ', npos
print *, '    # of samples = ', nind

! read hmm input to get marker positions
read(u,'(a)') junk1000 ! read junk line in hmm input (header line)
allocate(pos(npos), gc(nind,npos), gcped(npos*2))
do i=1, npos
        read(u,*) pos(i)
enddo
from=minval(pos)
to=maxval(pos)
pos=pos-from+1
l=to-from+1
close(u)

! read hmm output (variable line length), 2 haplotypes at a time, and store the
! genotype as 0 (modern homozygote), 1 (heterozygote) and 2 (archaic homozygote)
! and write a file (1 row/sample; 2 columns/marker) with the gene content as A
! (archaic) or C (modern)
open(newunit=u, file='archaic.gcped')
write(fmt1,'(a,i6,a)') '(a7,', npos*2, 'a2)'
gc=0
allocate(w3(l), inds(npos))
inds=(/(i,i=1,npos)/)
do i=1,nind
        gcped='C'
        do j=1,2
                read(v,'(a)') junk1000
                read(junk1000,'(1000a1)') junk
                m=count(iachar(junk)==9)
                if (m==1) then
                        read(junk1000,*) file3
                        cycle
                endif
                allocate(w1(m-1))
                read(junk1000,*) file3, w1
                w3=0
                w1=w1-from+1
                do k=1,m-1,2
                        w3(w1(k):w1(k+1))=1
                enddo
                n=count(w3(pos)==1) ! as m>1, n cannot be equal to 0
                allocate(w5(n))
                w5=pack(inds,w3(pos)==1)
                gcped((w5-1)*2+j)='A'
                gc(i,w5)=gc(i,w5)+1
                deallocate(w1, w5)
        enddo
        write(u,fmt1) file3(:7), gcped
enddo
deallocate(w3)
close(v)
close(u)


! output frequencies and map file
open(newunit=u, file='archaic.frq')
open(newunit=v, file='archaic.map')
a=0.0
do i=1, npos
        write(v,'(i2,a,i3,i10)') chr, ' marker', 0, pos(i)+from-1
        j=sum(gc(:,i))
        r=.5*real(j)/real(nind)
        if (r>a) a=r
        write(u,'(i10,f9.4)') pos(i), r
enddo
close(u)
close(v)
print *, ' max archaic frequency = ', a


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

end program
