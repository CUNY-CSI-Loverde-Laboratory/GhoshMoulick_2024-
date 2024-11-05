!Program for calculating nucleosomal basepair displacement with respect to COM DNA 
!heavy atoms in the crystal structure

!The code is written by Prabir Khatua

!Date: 26th April, 2022
!Please cite Khatua et al. JPCB, 2024, 128 (13), 3090-3101
!Usage: f95 BpairDisplWrtCrystal.F90 or gfortran BpairDisplWrtCrystal.F90

!       ./a.out  (for interactive run and then enter the inputs interactivly 
!                 as the program will ask)
!       
!       ./a.out<input.inp >output.log &  
!                            (input.inp is the input file containing all the 
!                            required inputs that has to be prepared by the 
!                            user. For this, one needs to open the code or 
!                            better compile the code once in interactive mode 
!                            to know what set of inputs will be required and 
!                            prepare the input file. output.log will print 
!                            general information or error) 

!The code is written to analyse the dcd trajectories
!------------------------------------------------------------------------------
program OrderParameter
real,dimension(:),allocatable:: x,y,z,xref,yref,zref,mass
real(kind=8):: RefX_Tail1,RefY_Tail1,RefZ_Tail1
real(kind=8):: RefX_Tail2,RefY_Tail2,RefZ_Tail2
real(kind=8):: SimX_Tail1,SimY_Tail1,SimZ_Tail1
real(kind=8):: SimX_Tail2,SimY_Tail2,SimZ_Tail2
real(kind=8),dimension(30):: ComAX,ComAY,ComAZ,ComBX,ComBY,ComBZ
real(kind=8):: SenseBaseComX,SenseBaseComY,SenseBaseComZ
real(kind=8):: AntiSenseBaseComX,AntiSenseBaseComY,AntiSenseBaseComZ
real(kind=8),dimension(5):: BaseComX,BaseComY,BaseComZ
integer,dimension(:),allocatable:: SenseDnaResNatom,AntiSenseDnaResNatom
integer,dimension(:),allocatable:: ResID,DnaAtomId
integer,dimension(:,:),allocatable:: SenseDnaAtomId,AntiSenseDnaAtomId
integer,dimension(5):: SelectedBaseID
integer:: dummyi,punit,SHL0,ResNo
character(len=200),dimension(:),allocatable:: dcdfile
character(len=200):: pdbfile,Prefix
character(len=300):: RTail1,RTail2,ATail1,ATail2,InnerGyre
character(len=4),dimension(:),allocatable:: AtomName,ResName
character(len=1):: ElementName
logical:: there,ANS
real(kind=8):: dx,dy,dz,bx,by,bz
real(kind=8):: DistTail1,DistTail2,RefDistTail1,RefDistTail2
real(kind=8):: AngleXY_Tail1,AngleYZ_Tail1,AngleZX_Tail1,Angle_Tail1
real(kind=8):: AngleXY_Tail2,AngleYZ_Tail2,AngleZX_Tail2,Angle_Tail2
real(kind=8):: ComX,ComY,ComZ,cx,cy,cz
real(kind=8):: time,dt,RefGyre,SimGyre,Dist
data inunit,punit /10,11/
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the prefix for all output files'
read(*,'(a)')Prefix
write(*,*)Prefix

RTail1 = trim(Prefix)//"_Tail1_SHL0_Breathing_distance_time.out"
RTail2 = trim(Prefix)//"_Tail2_SHL0_Breathing_distance_time.out"
ATail1 = trim(Prefix)//"_Tail1_Breathing_Angle_Displ_wrt_crystal_time.out" 
ATail2 = trim(Prefix)//"_Tail2_Breathing_Angle_Displ_wrt_crystal_time.out"
InnerGyre = trim(Prefix)//"_Gyre_Breathing_wrt_crystal_dist_time.out"

open(12,file=RTail1,status='unknown')
open(13,file=RTail2,status='unknown')
open(14,file=ATail1,status='unknown')
open(15,file=ATail2,status='unknown')
open(17,file=InnerGyre,status='unknown')
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the no of total atoms, natom'
read(*,*)natom
write(*,*)'natom=====>',natom
write(*,'(1x,a)')'Enter the first and last atom no for DNA, dna_first, dna_last'
read(*,*)ifirst_dna,ilast_dna
write(*,*)'ifirst_dna,ilast_dna=====>',ifirst_dna,ilast_dna
write(*,'(1x,a)')'Enter the total no of bases present in DNA'
read(*,*)NumDnaBase
write(*,*)'NumDnaBase====>',NumDnaBase
write(*,'(1x,a)')'Choose from the following options'
write(*,'(1x,a)')'Type 0 if the trajectory is saved as wrapped coordinates'
write(*,'(1x,a)')'Type 1 if the trajectory is not saved as wrapped coordinates'
read(*,*)iwrap
write(*,*)'You have selected the option',iwrap
write(*,'(1x,a)')'Enter the timestep in ps, dt'
read(*,*)dt
write(*,*)'dt=====>',dt
write(*,'(1x,a)')'Enter the skip value, nskip'
read(*,*)nskip
write(*,*)'nskip====>',nskip
write(*,'(1x,a)')'Enter the SHL0 base pair no index'
read(*,*)SHL0
write(*,*)'SHL0====>',SHL0

NumDnaBasePair = NumDnaBase/2
write(*,*)'# of Nucleosomal DNA base pair===>',NumDnaBasePair

nDNA = ilast_DNA - ifirst_DNA + 1

write(*,*)'# of DNA atoms========>',nDNA

SelectedBaseID(1) = SHL0
SelectedBaseID(2) = 1
SelectedBaseID(3) = 10
SelectedBaseID(4) = NumDnaBasePair
SelectedBaseID(5) = NumDnaBasePair - 9

allocate(x(natom),y(natom),z(natom),mass(natom))
allocate(xref(natom),yref(natom),zref(natom))
allocate(SenseDnaResNatom(NumDnaBasePair),AntiSenseDnaResNatom(NumDnaBasePair))
allocate(DnaAtomId(nDNA))
allocate(SenseDnaAtomId(NumDnaBasePair,100),AntiSenseDnaAtomId(NumDnaBasePair,100))
allocate(AtomName(natom),ResName(natom),ResID(natom))
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the name of pdbfile'
read(*,'(a)')pdbfile
write(*,*)pdbfile

inquire(file=pdbfile,exist=there)
if(.not.there)then
write(*,'(a)')'ERROR: pdbfile NOT FOUND'
stop
endif

open(punit,file=pdbfile,status='old',form='formatted',iostat=ios)

read(punit,111)AtomName(1),ResName(1),ResNo,xref(1),yref(1),zref(1)

close(punit)

open(punit,file=pdbfile,status='old',form='formatted',iostat=ios)

Npair = 0

SenseDnaResNatom(:) = 0
AntiSenseDnaResNatom(:) = 0

do i=1,natom
read(punit,111)AtomName(i),ResName(i),Num,xref(i),yref(i),zref(i)
111   format(12x,a4,1x,a4,1x,i4,4x,3f8.3)

ElementName=adjustl(AtomName(i))

ResID(i) = Num - ResNo + 1
!write(*,*) ResID(i)
if(ElementName == 'H')mass(i) = 1.0080
if(ElementName == 'C')mass(i) = 12.0110
if(ElementName == 'N')mass(i) = 14.0070
if(ElementName == 'O')mass(i) = 15.9990
if(ElementName == 'S')mass(i) = 32.0600
if(ElementName == 'P')mass(i) = 30.9740

if(ElementName == 'H')cycle

if(i >= ifirst_dna.and.i <= ilast_dna)then

        Npair = Npair + 1
        DnaAtomID(Npair) = i
	ResID(i) = ResID(i) - 956
	!write(*,*) ResID(i)
        if(ResID(i) <= NumDnaBasePair)then
                SenseDnaResNatom(ResID(i)) = SenseDnaResNatom(ResID(i)) + 1
                SenseDnaAtomId(ResID(i),SenseDnaResNatom(ResID(i))) = i
        else
                IDBasePair = NumDnaBase + 1 - ResID(i)
                AntiSenseDnaResNatom(IDBasePair) = AntiSenseDnaResNatom(IDBasePair) + 1
                AntiSenseDnaAtomId(IDBasePair,AntiSenseDnaResNatom(IDBasePair)) = i
        endif
endif

enddo

write(*,*)'# of selected DNA heavy atom====>',Npair

call CENTER(natom,Npair,DnaAtomID,xref,yref,zref,mass,ComX,ComY,ComZ)

do ii=1,5
i = SelectedBaseID(ii)
call DNA_BASE_COM(i,NumDnaBasePair,SenseDnaResNatom(i),SenseDnaAtomId,mass,xref,yref,zref,&
        &SenseBaseComX,SenseBaseComY,SenseBaseComZ)
call DNA_BASE_COM(i,NumDnaBasePair,AntiSenseDnaResNatom(i),AntiSenseDnaAtomId,mass,&
        &xref,yref,zref,AntiSenseBaseComX,AntiSenseBaseComY,AntiSenseBaseComZ)

BaseComX(ii) = 0.5*(SenseBaseComX + AntiSenseBaseComX)
BaseComY(ii) = 0.5*(SenseBaseComY + AntiSenseBaseComY)
BaseComZ(ii) = 0.5*(SenseBaseComZ + AntiSenseBaseComZ)

enddo

RefX_Tail1 = BaseComX(2) - BaseComX(3)
RefY_Tail1 = BaseComY(2) - BaseComY(3)
RefZ_Tail1 = BaseComZ(2) - BaseComZ(3)

RefX_Tail2 = BaseComX(4) - BaseComX(5)
RefY_Tail2 = BaseComY(4) - BaseComY(5)
RefZ_Tail2 = BaseComY(4) - BaseComY(5)

call DISTANCE(BaseComX(1),BaseComY(1),BaseComZ(1),BaseComX(2)&
        &,BaseComY(2),BaseComZ(2),0,10.d0,10.d0,10.d0,RefDistTail1)
call DISTANCE(BaseComX(1),BaseComY(1),BaseComZ(1),BaseComX(4)&
        &,BaseComY(4),BaseComZ(4),0,10.d0,10.d0,10.d0,RefDistTail2)
!-----------------------------------------------------------------------------
iGyreFirst = SHL0 - 56
iGyreLast = SHL0 - 27

jGyreFirst = SHL0 + 21
jGyreLast = SHL0 + 50

nCounter = 0
do i=iGyreFirst,iGyreLast
nCounter = nCounter + 1
call DNA_BASE_COM(i,NumDnaBasePair,SenseDnaResNatom(i),SenseDnaAtomId,mass,xref,yref,zref,&
        &SenseBaseComX,SenseBaseComY,SenseBaseComZ)
call DNA_BASE_COM(i,NumDnaBasePair,AntiSenseDnaResNatom(i),AntiSenseDnaAtomId,mass,&
        &xref,yref,zref,AntiSenseBaseComX,AntiSenseBaseComY,AntiSenseBaseComZ)

ComAX(nCounter) = 0.5*(SenseBaseComX + AntiSenseBaseComX)
ComAY(nCounter) = 0.5*(SenseBaseComY + AntiSenseBaseComY)
ComAZ(nCounter) = 0.5*(SenseBaseComZ + AntiSenseBaseComZ)

enddo

nCounter = 0

do i=jGyreFirst,jGyreLast
nCounter = nCounter + 1
call DNA_BASE_COM(i,NumDnaBasePair,SenseDnaResNatom(i),SenseDnaAtomId,mass,xref,yref,zref,&
        &SenseBaseComX,SenseBaseComY,SenseBaseComZ)
call DNA_BASE_COM(i,NumDnaBasePair,AntiSenseDnaResNatom(i),AntiSenseDnaAtomId,mass,&
        &xref,yref,zref,AntiSenseBaseComX,AntiSenseBaseComY,AntiSenseBaseComZ)

ComBX(nCounter) = 0.5*(SenseBaseComX + AntiSenseBaseComX)
ComBY(nCounter) = 0.5*(SenseBaseComY + AntiSenseBaseComY)
ComBZ(nCounter) = 0.5*(SenseBaseComZ + AntiSenseBaseComZ)

enddo

RefGyre = 0.d0

do i = 1, nCounter
call DISTANCE(ComAX(i),ComAY(i),ComAZ(i),ComBX(i),ComBY(i),ComBZ(i),0,10.d0,10.d0,10.d0,Dist)
RefGyre = RefGyre + Dist
enddo

RefGyre = RefGyre/dble(nCounter)

write(*,*)'Average Inter-Gyre COM distance in the crystal structure==>',RefGyre
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the number of trajectory files you want to analyze.'
read(*,*) ninput
write(*,*)'ninput==========>',ninput

allocate(dcdfile(ninput))

do inp=1,ninput
write(*,*)'Enter the name of trajectory file',inp
read(*,'(a)')dcdfile(inp)
write(*,*)dcdfile(inp)
inquire(file=dcdfile(inp),exist=there)
if(.not.there)then
write(*,'(a)')'ERROR: dcdfile NOT FOUND'
stop
endif
enddo
!------------------------------------------------------------------------------
nread = 0
nused = 0

time = 0.d0
dt = dt*nskip

do inp=1,ninput
write(*,'(1x,a,i5,a)')'Enter name of trajectory file',inp
write(*,*)dcdfile(inp)
open(inunit,file=dcdfile(inp),status='old',form='unformatted')
read(inunit)dummyc,nframes,(dummyi,i=1,8),dummyr,(dummyi,i=1,9)
read(inunit)dummyi,dummyr
write(*,*)"dummyi,dummyr====>",dummyi,dummyr
read(inunit)natom
write(*,*)"natom====>",natom
write(*,*)"nframes===>",nframes
do iframe=1,nframes
read(inunit)dx,bx,dy,by,bz,dz
read(inunit)(x(j),j=1,natom)
read(inunit)(y(j),j=1,natom)
read(inunit)(z(j),j=1,natom)
nread=nread+1
if(mod(nread,nskip) /= 0)cycle
nused=nused+1
!------------------------------------------------------------------------------
time = time + dt

call CENTER(natom,Npair,DnaAtomID,x,y,z,mass,cx,cy,cz)

call ALIGN(natom,Npair,DnaAtomID,x,y,z,xref,yref,zref)

do ii=1,5
i = SelectedBaseID(ii)
call DNA_BASE_COM(i,NumDnaBasePair,SenseDnaResNatom(i),SenseDnaAtomId,mass,x,y,z,&
        &SenseBaseComX,SenseBaseComY,SenseBaseComZ)
call DNA_BASE_COM(i,NumDnaBasePair,AntiSenseDnaResNatom(i),AntiSenseDnaAtomId,mass,&
        &x,y,z,AntiSenseBaseComX,AntiSenseBaseComY,AntiSenseBaseComZ)

BaseComX(ii) = 0.5*(SenseBaseComX + AntiSenseBaseComX)
BaseComY(ii) = 0.5*(SenseBaseComY + AntiSenseBaseComY)
BaseComZ(ii) = 0.5*(SenseBaseComZ + AntiSenseBaseComZ)
enddo

SimX_Tail1 = BaseComX(2) - BaseComX(3)
SimY_Tail1 = BaseComY(2) - BaseComY(3)
SimZ_Tail1 = BaseComZ(2) - BaseComZ(3)

SimX_Tail2 = BaseComX(4) - BaseComX(5)
SimY_Tail2 = BaseComY(4) - BaseComY(5)
SimZ_Tail2 = BaseComY(4) - BaseComY(5)

call DISTANCE(BaseComX(1),BaseComY(1),BaseComZ(1),BaseComX(2)&
        &,BaseComY(2),BaseComZ(2),iwrap,dx,dy,dz,DistTail1)

call DISTANCE(BaseComX(1),BaseComY(1),BaseComZ(1),BaseComX(4)&
        &,BaseComY(4),BaseComZ(4),iwrap,dx,dy,dz,DistTail2)


call COMP_ANGLE(RefX_Tail1,RefY_Tail1,RefZ_Tail1,SimX_Tail1,SimY_Tail1,SimZ_Tail1,&
        &AngleXY_Tail1,AngleYZ_Tail1,AngleZX_Tail1,Angle_Tail1)

call COMP_ANGLE(RefX_Tail2,RefY_Tail2,RefZ_Tail2,SimX_Tail2,SimY_Tail2,SimZ_Tail2,&
        &AngleXY_Tail2,AngleYZ_Tail2,AngleZX_Tail2,Angle_Tail2)


write(14,'(f12.6,1x,4f8.3)')time/1000000.0,Angle_Tail1,AngleXY_Tail1,AngleYZ_Tail1,AngleZX_Tail1
write(15,'(f12.6,1x,4f8.3)')time/1000000.0,Angle_Tail2,AngleXY_Tail2,AngleYZ_Tail2,AngleZX_Tail2

write(12,'(f12.6,1x,2f8.3)')time/1000000.0,DistTail1 - RefDistTail1,RefDistTail1
write(13,'(f12.6,1x,2f8.3)')time/1000000.0,DistTail2 - RefDistTail2,RefDistTail2

!Once again this part could have been written in subroutine

nCounter = 0

do i=iGyreFirst,iGyreLast
nCounter = nCounter + 1
call DNA_BASE_COM(i,NumDnaBasePair,SenseDnaResNatom(i),SenseDnaAtomId,mass,x,y,z,&
        &SenseBaseComX,SenseBaseComY,SenseBaseComZ)
call DNA_BASE_COM(i,NumDnaBasePair,AntiSenseDnaResNatom(i),AntiSenseDnaAtomId,mass,&
        &x,y,z,AntiSenseBaseComX,AntiSenseBaseComY,AntiSenseBaseComZ)

ComAX(nCounter) = 0.5*(SenseBaseComX + AntiSenseBaseComX)
ComAY(nCounter) = 0.5*(SenseBaseComY + AntiSenseBaseComY)
ComAZ(nCounter) = 0.5*(SenseBaseComZ + AntiSenseBaseComZ)

enddo

nCounter = 0

do i=jGyreFirst,jGyreLast
nCounter = nCounter + 1
call DNA_BASE_COM(i,NumDnaBasePair,SenseDnaResNatom(i),SenseDnaAtomId,mass,x,y,z,&
        &SenseBaseComX,SenseBaseComY,SenseBaseComZ)
call DNA_BASE_COM(i,NumDnaBasePair,AntiSenseDnaResNatom(i),AntiSenseDnaAtomId,mass,&
        &x,y,z,AntiSenseBaseComX,AntiSenseBaseComY,AntiSenseBaseComZ)

ComBX(nCounter) = 0.5*(SenseBaseComX + AntiSenseBaseComX)
ComBY(nCounter) = 0.5*(SenseBaseComY + AntiSenseBaseComY)
ComBZ(nCounter) = 0.5*(SenseBaseComZ + AntiSenseBaseComZ)

enddo

SimGyre = 0.d0
do i=1,nCounter
call DISTANCE(ComAX(i),ComAY(i),ComAZ(i),ComBX(i),ComBY(i),ComBZ(i),iwrap,dx,dy,dz,Dist)
SimGyre = SimGyre + Dist
enddo

SimGyre = SimGyre/dble(nCounter)

write(17,'(f12.6,1x,2f8.3)')time/1000000.0,SimGyre - RefGyre,RefGyre

enddo
enddo

end program OrderParameter
!------------------------------------------------------------------------------
subroutine COMP_ANGLE(aX,aY,aZ,bX,bY,bZ,AngleXY,AngleYZ,AngleZX,Angle)
real(kind=8),intent(in):: aX,aY,aZ,bX,bY,bZ
real(kind=8),intent(out):: AngleXY,AngleYZ,AngleZX,Angle


pi = 4.0 * datan(1.d0)

AngleXY = 180.0*acos((aX*bX + aY*bY)/(dsqrt(aX*aX + aY*aY) * dsqrt(bX*bX + bY*bY)))/pi
AngleYZ = 180.0*acos((aZ*bZ + aY*bY)/(dsqrt(aZ*aZ + aY*aY) * dsqrt(bZ*bZ + bY*bY)))/pi
AngleZX = 180.0*acos((aX*bX + aZ*bZ)/(dsqrt(aX*aX + aZ*aZ) * dsqrt(bX*bX + bZ*bZ)))/pi
Angle = 180.0*acos((aX*bX + aY*bY + aZ*bZ)/(dsqrt(aX*aX + aY*aY + aZ*aZ) * dsqrt(bX*bX + bY*bY + bZ*bZ)))/pi

end subroutine COMP_ANGLE
!------------------------------------------------------------------------------
subroutine DISTANCE(aX,aY,aZ,bX,bY,bZ,iwrap,BoxX,BoxY,BoxZ,Modab)  ! AB = B - A
integer,intent(in):: iwrap
real(kind=8),intent(in):: aX,aY,aZ,bX,bY,bZ,BoxX,BoxY,BoxZ
real(kind=8):: abX,abY,abZ
real(kind=8),intent(out):: Modab

abX = bX - aX
abY = bY - aY
abZ = bZ - aZ

if(iwrap == 1)then
abX = abX - anint(abX/BoxX)*BoxX
abY = abY - anint(abY/BoxY)*BoxY
abZ = abZ - anint(abZ/BoxZ)*BoxZ
endif

Modab = sqrt(abX*abX + abY*abY + abZ*abZ)

end subroutine DISTANCE
!------------------------------------------------------------------------------
subroutine DNA_BASE_COM(i,Nbase,ResNatom,AtomId,mass,x,y,z,ComX,ComY,ComZ)
integer,intent(in):: i,Nbase,ResNatom
integer,dimension(Nbase,*),intent(in):: AtomId
real,dimension(*),intent(in):: x,y,z
real,dimension(*),intent(in):: mass
real(kind=8),intent(out):: ComX,ComY,ComZ
real(kind=8):: TotalMass

TotalMass = 0.d0
ComX = 0.d0
ComY = 0.d0 
ComZ = 0.d0

do j=1,ResNatom

ii = AtomId(i,j)

ComX = ComX + x(ii)*mass(ii)
ComY = ComY + y(ii)*mass(ii)
ComZ = ComZ + z(ii)*mass(ii)

TotalMass = TotalMass + mass(ii)

enddo

ComX = ComX/TotalMass
ComY = ComY/TotalMass
ComZ = ComZ/TotalMass

end subroutine DNA_BASE_COM
!------------------------------------------------------------------------------
subroutine COM(Npair,AtomId,mass,x,y,z,ComX,ComY,ComZ)
real,dimension(*),intent(in):: x,y,z
real,dimension(*),intent(in):: mass
integer,dimension(*),intent(in):: AtomId
integer,intent(in):: Npair
real(kind=8),intent(out):: ComX,ComY,ComZ
real(kind=8):: TotalMass

TotalMass = 0.d0
ComX = 0.d0
ComY = 0.d0
ComZ = 0.d0

do i=1,Npair

ii = AtomId(i)

ComX = ComX + x(ii)*mass(ii)
ComY = ComY + y(ii)*mass(ii)
ComZ = ComZ + z(ii)*mass(ii)

TotalMass = TotalMass + mass(ii)

enddo

ComX = ComX/TotalMass
ComY = ComY/TotalMass
ComZ = ComZ/TotalMass

end subroutine COM
!------------------------------------------------------------------------------
subroutine ALIGN(npr,np,atp,xa,ya,za,xr,yr,zr)
real,dimension(*),intent(inout):: xa,ya,za
real,dimension(*),intent(in):: xr,yr,zr
integer,dimension(*),intent(in):: atp
integer,intent(in):: np,npr
real,dimension(4):: q
real,dimension(4,4):: m

m(1,1)=0.0
m(1,2)=0.0
m(1,3)=0.0
m(1,4)=0.0
m(2,2)=0.0
m(2,3)=0.0
m(2,4)=0.0
m(3,3)=0.0
m(3,4)=0.0
m(4,4)=0.0
do j=1,np
i=atp(j)
m(1,1)=m(1,1)+xa(i)*xr(i)+ya(i)*yr(i)+za(i)*zr(i)
m(2,2)=m(2,2)+xa(i)*xr(i)-ya(i)*yr(i)-za(i)*zr(i)
m(3,3)=m(3,3)+ya(i)*yr(i)-za(i)*zr(i)-xa(i)*xr(i)
m(4,4)=m(4,4)+za(i)*zr(i)-xa(i)*xr(i)-ya(i)*yr(i)
m(1,2)=m(1,2)+ya(i)*zr(i)-za(i)*yr(i)
m(1,3)=m(1,3)+za(i)*xr(i)-xa(i)*zr(i)
m(1,4)=m(1,4)+xa(i)*yr(i)-ya(i)*xr(i)
m(2,3)=m(2,3)+xa(i)*yr(i)+ya(i)*xr(i)
m(2,4)=m(2,4)+za(i)*xr(i)+xa(i)*zr(i)
m(3,4)=m(3,4)+ya(i)*zr(i)+za(i)*yr(i)
end do
m(2,1)=m(1,2)
m(3,1)=m(1,3)
m(4,1)=m(1,4)
m(3,2)=m(2,3)
m(4,2)=m(2,4)
m(4,3)=m(3,4)

call JACOBI(q,m,4)
call ROTN(xa,ya,za,q,npr)

end subroutine ALIGN
!------------------------------------------------------------------------------
subroutine ROTN(xa,ya,za,q,n)
integer,intent(in):: n
real,dimension(*),intent(inout):: xa,ya,za
real,dimension(*),intent(in):: q
real,dimension(n):: a,b,c
real:: q0,q2,qd,qv,q1

do i=1,n
a(i)=xa(i)
b(i)=ya(i)
c(i)=za(i)
end do
q0=q(1)
q1=2*q0
q0=q0**2
q2=q(2)**2+q(3)**2+q(4)**2
qd=q0-q2
do i=1,n
qv=q(2)*a(i)+q(3)*b(i)+q(4)*c(i)
qv=2*qv
xa(i)=qd*a(i)+qv*q(2)+q1*(q(3)*c(i)-q(4)*b(i))
ya(i)=qd*b(i)+qv*q(3)+q1*(q(4)*a(i)-q(2)*c(i))
za(i)=qd*c(i)+qv*q(4)+q1*(q(2)*b(i)-q(3)*a(i))
end do

end subroutine ROTN
!------------------------------------------------------------------------------
SUBROUTINE JACOBI(Q,A,N)
INTEGER::I,J,IP,IQ,M,N
INTEGER::NROT,SWEEP
REAL,DIMENSION(1:4,1:4),INTENT(INOUT)::A
REAL,DIMENSION(1:4,1:4)::V
REAL,DIMENSION(1:4),INTENT(OUT)::Q
REAL,DIMENSION(1:4)::D(4),B(4),Z(4)
REAL::THETA,TAU,TRESH,C,S,T,G,H,SM,DMX
!
!     ----Jacobi Method---
!     ----Ref: Numerical Recipe in Fortran; page 460 
!     ----Initialization of eigenvector martix as unit matrix---
DO IP=1,N
DO IQ=1,N
V(IP,IQ)=0.0
END DO
V(IP,IP)=1.0
END DO
!     ----Initialization of eigenvector as the diagonal elements---
DO IP=1,N
B(IP)=A(IP,IP)
D(IP)=B(IP)
Z(IP)=0.0
END DO
!     ----Iteration---
NROT=0
SWEEP=0
112  SWEEP=SWEEP+1
IF(SWEEP>150)THEN
PRINT*,"PROGRAM STOP: TOO MANY ITERATION IN THE CALCULATION."
STOP
END IF
SM=0.0                                     !Sum of absolute value
DO IP=1,N-1                                !of super off-diagonal
DO IQ=IP+1,N                             !elements. Program will
SM=SM+ABS(A(IP,IQ))                      !stop if the sum become
END DO                                   !zero.
END DO
IF(SM==0.0)GOTO 132
!     -----------------------------------------------------------------
IF(SWEEP<4)THEN                         !Read the paragraph in
TRESH=0.2*SM/N**2                        !Numerical Recipes in
ELSE                                     !Fortan, on the page 
TRESH=0.0                                !No. 459
END IF
!     -----------------------------------------------------------------
DO IP=1,N-1
DO IQ=IP+1,N
!        ----------calculation in a sweep-----
G=100.0*ABS(A(IP,IQ))
!After 4 step skip the rotation if of diagonal element is small
IF((SWEEP>4).AND.(ABS(D(IP))+G==ABS(D(IP))).AND.(ABS(D(IQ))+G==ABS(D(IQ))))THEN
A(IP,IQ)=0.0
ELSE IF(ABS(A(IP,IQ))>TRESH)THEN    
H=D(IQ)-D(IP)                        
IF(ABS(H)+G==ABS(H))THEN           
T=A(IP,IQ)/H
ELSE                                     
THETA=0.5*H/A(IP,IQ)               !Rotation parameters calculation
T=1.0/(ABS(THETA)+SQRT(1.0+ THETA**2))  
IF(THETA<0.0)T=-T               
END IF                               
C=1.0/SQRT(1.0+T**2)                 
S=T*C                                
TAU=S/(1.0+C)
!     ----operation corresponding rotation-----------------------------
H=T*A(IP,IQ)
Z(IP)=Z(IP)-H
Z(IQ)=Z(IQ)+H
D(IP)=D(IP)-H
D(IQ)=D(IQ)+H
A(IP,IQ)=0.0
!          -----------------------
DO J=1,IP-1                          !Changes in ip coulumn
G=A(J,IP)                          !and iq coulumn upto
H=A(J,IQ)                          !ip row. Stop before
A(J,IP)=G-S*(H+G*TAU)              !ip row. 
A(J,IQ)=H+S*(G-H*TAU)              !
END DO                               !
!           -----------------------
DO J=IP+1,IQ-1                       !Change in ip row from
G=A(IP,J)                          !ip coulumn to iq coulu
H=A(J,IQ)                          !mn and change in iq co
A(IP,J)=G-S*(H+G*TAU)              !ulumn from ip row to 
A(J,IQ)=H+S*(G-H*TAU)              !iq row.
END DO
!           -----------------------
DO J=IQ+1,N                          !change in ip row from
G=A(IP,J)                          !iq coulumn to end coul
H=A(IQ,J)                          !umn and change in iq 
A(IP,J)=G-S*(H+G*TAU)              !row from iq coulumn 
A(IQ,J)=H+S*(G-H*TAU)              !to end coulumn.
END DO                               !
!           -----------------------
DO J=1,N                             !
G=V(J,IP)                          !Update of the eigen-
H=V(J,IQ)                          !vector matrix.
V(J,IP)=G-S*(H+G*TAU)              !
V(J,IQ)=H+S*(G-H*TAU)              !
END DO                               !
!           -----------------------
NROT=NROT+1
END IF
END DO
END DO
!     ----------------
DO IP=1,N
B(IP)=B(IP)+Z(IP)                          !Update of d with sum 
D(IP)=B(IP)                                !of ta(p,q) and reinit
Z(IP)=0.0                                  !ialization of z
end do
!     ----------------
SWEEP=SWEEP+1
GOTO 112
132   DMX=D(1)
M=1
DO IP=1,N
IF(D(IP)>DMX)THEN
DMX=D(IP)
M=IP
END IF
END DO
DO IP=1,N
Q(IP)=V(IP,M)
END DO
END SUBROUTINE JACOBI
!------------------------------------------------------------------------------
subroutine CENTER(npr,np,atp,x,y,z,xm,cx,cy,cz)
real,dimension(*),intent(inout):: x,y,z
real,dimension(*),intent(in):: xm
integer,dimension(*),intent(in):: atp
integer,intent(in):: np,npr
real(kind=8),intent(out):: cx,cy,cz

cx=0.0
cy=0.0
cz=0.0
tm=0.0
do i=1,np
ii=atp(i)
cx=cx+xm(ii)*x(ii)
cy=cy+xm(ii)*y(ii)
cz=cz+xm(ii)*z(ii)
tm=tm+xm(ii)
enddo

cx=cx/tm
cy=cy/tm
cz=cz/tm

do i=1,npr
x(i)=x(i)-cx
y(i)=y(i)-cy
z(i)=z(i)-cz
enddo

end subroutine CENTER
!------------------------------------------------------------------------------
