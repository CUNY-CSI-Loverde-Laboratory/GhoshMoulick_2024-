!Program for calculating Rotational Order Parameter of Nucleosome

!Reference: PNAS, 2017, 114, E9197-E9205
!           PRL, 2014, 113, 168101
!           JCP, 2019, 150, 215102
!			JPCB, 2024, 128 (13), 3090-3101

!The code is written by Prabir Khatua 

!Date: 1st November, 2021
!Finalized, Generalized, and Correctd on: 10th November, 2021
!Fitting version writen on: 3rd May, 2022

!Usage: f95 NucleosomeRepositionOrderHistoneFitting.F90 or gfortran NucleosomeRepositionOrderHistoneFitting.F90

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
real,dimension(:),allocatable:: x,y,z,xref,yref,zref
real(kind=8),dimension(:),allocatable:: BxBaseSim,ByBaseSim,BzBaseSim
real(kind=8),dimension(:),allocatable:: PxBaseSim,PyBaseSim,PzBaseSim,mass
real(kind=8),dimension(:),allocatable:: PxBaseRef,PyBaseRef,PzBaseRef
real(kind=8):: SenseBaseComX,SenseBaseComY,SenseBaseComZ
real(kind=8):: AntiSenseBaseComX,AntiSenseBaseComY,AntiSenseBaseComZ
real(kind=8):: BaseComX,BaseComY,BaseComZ,ModP
integer,dimension(:),allocatable:: SenseDnaResNatom,AntiSenseDnaResNatom
integer,dimension(:),allocatable:: ResID
integer,dimension(:),allocatable:: H31_AtomId,H32_AtomId,HistoneAtomId,SelectedBaseID
integer,dimension(:,:),allocatable:: SenseDnaAtomId,AntiSenseDnaAtomId
integer:: dummyi,ounit,punit,NumDnaBase,NumDnaBasePair,CalculationType
character(len=100),dimension(:),allocatable:: dcdfile
character(len=100):: pdbfile,outfile
character(len=4),dimension(:),allocatable:: AtomName,ResName
character(len=1):: ElementName
logical:: there,ANS
real(kind=8):: ComX5Prime,ComY5Prime,ComZ5Prime
real(kind=8):: ComX3Prime,ComY3Prime,ComZ3Prime
real(kind=8):: DxStrand,DyStrand,DzStrand
real(kind=8):: UX,UY,UZ
real(kind=8):: dx,dy,dz,bx,by,bz,Srot,SrotX,SrotY,Strans
real(kind=8):: HistoneComX,HistoneComY,HistoneComZ
real(kind=8):: H31_ComX,H31_ComY,H31_ComZ
real(kind=8):: H32_ComX,H32_ComY,H32_ComZ
real(kind=8):: WX,WY,WZ,FX,FY,FZ
real(kind=8):: time,dt,Lambda,pi
real(kind=8):: CrosX,CrosY,CrosZ
data inunit,punit,ounit /10,11,12/
!------------------------------------------------------------------------------
Lambda = 0.08
pi = 4.d0*datan(1.d0)
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the name of output file for contact map'
read(*,'(a)')outfile
write(*,*)outfile
open(ounit,file=outfile,status='unknown',form='formatted',iostat=ios)
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the no of total atoms, natom'
read(*,*)natom
write(*,*)'natom=====>',natom
write(*,'(1x,a)')'Enter the first and last atom no for DNA, dna_first, dna_last'
read(*,*)ifirst_dna,ilast_dna
write(*,*)'ifirst_dna,ilast_dna=====>',ifirst_dna,ilast_dna
write(*,'(1x,a)')'Enter the first and last atom nos for core regions of each of eight histone subunits'

write(*,'(1x,a)')'For H31:'
read(*,*)ifirst_H31,ilast_H31
write(*,*)'ifirst_H31,ilast_H31======>',ifirst_H31,ilast_H31

write(*,'(1x,a)')'For H41:'
read(*,*)ifirst_H41,ilast_H41
write(*,*)'ifirst_H41,ilast_H41======>',ifirst_H41,ilast_H41

write(*,'(1x,a)')'For H2A1:'
read(*,*)ifirst_H2A1,ilast_H2A1
write(*,*)'ifirst_H2A1,ilast_H2A1======>',ifirst_H2A1,ilast_H2A1

write(*,'(1x,a)')'For H2B1:'
read(*,*)ifirst_H2B1,ilast_H2B1
write(*,*)'ifirst_H2B1,ilast_H2B1======>',ifirst_H2B1,ilast_H2B1

write(*,'(1x,a)')'For H32:'
read(*,*)ifirst_H32,ilast_H32
write(*,*)'ifirst_H32,ilast_H32======>',ifirst_H32,ilast_H32

write(*,'(1x,a)')'For H42:'
read(*,*)ifirst_H42,ilast_H42
write(*,*)'ifirst_H42,ilast_H42======>',ifirst_H42,ilast_H42

write(*,'(1x,a)')'For H2A2:'
read(*,*)ifirst_H2A2,ilast_H2A2
write(*,*)'ifirst_H2A2,ilast_H2A2======>',ifirst_H2A2,ilast_H2A2

write(*,'(1x,a)')'For H2B2:'
read(*,*)ifirst_H2B2,ilast_H2B2
write(*,*)'ifirst_H2B2,ilast_H2B2======>',ifirst_H2B2,ilast_H2B2

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
write(*,'(1x,a)')'Enter how many bases will be considered for definig a DNA turn'
read(*,*)NbaseStep
write(*,*)'NbaseStep=====>',NbaseStep
write(*,'(1x,a)')'Choose from the following options'
write(*,'(1x,a)')'Type 0 if the calculation to be done as average over all DNA base pairs'
write(*,'(1x,a)')'Type 1 if the calculation to be done as average over specific DNA base pairs'
read(*,*)CalculationType
write(*,*)'You have chosen the option',CalculationType

NumDnaBasePair = NumDnaBase/2
write(*,*)'# of Nucleosomal DNA base pair===>',NumDnaBasePair

nHISTONE = ilast_H2B2 - ifirst_H31 + 1
nDNA = ilast_DNA - ifirst_DNA + 1

write(*,*)'# of DNA atoms========>',nDNA
write(*,*)'# of Histone atoms====>',nHISTONE

if(CalculationType == 1)then
        write(*,'(1x,a)')'Enter how many specific base pairs to be considered'
        read(*,*)nSelectedBase
        write(*,*)'nSelectedBase=====>',nSelectedBase

        allocate(SelectedBaseID(nSelectedBase))
        write(*,'(1x,a)')'Enter the basepair index for each of these bases in asci format'
        read(*,*)(SelectedBaseID(i),i=1,nSelectedBase)

        do i=1,nSelectedBase
        write(*,*)'Base Pair No==>',i,'Base Pair Index===>',SelectedBaseID(i)
        enddo

else
        nSelectedBase = NumDnaBasePair
        allocate(SelectedBaseID(NumDnaBasePair))
        do i=1,NumDnaBasePair
        SelectedBaseID(i) = i
        write(*,*)'Base Pair No==>',i,'Base Pair Index===>',SelectedBaseID(i)
        enddo
endif


allocate(x(natom),y(natom),z(natom),mass(natom))
allocate(xref(natom),yref(natom),zref(natom))
allocate(BxBaseSim(NumDnaBasePair),ByBaseSim(NumDnaBasePair))
allocate(BzBaseSim(NumDnaBasePair),PxBaseRef(NumDnaBasePair))
allocate(PyBaseRef(NumDnaBasePair),PzBaseRef(NumDnaBasePair))
allocate(PxBaseSim(NumDnaBasePair),PyBaseSim(NumDnaBasePair))
allocate(PzBaseSim(NumDnaBasePair),SenseDnaResNatom(NumDnaBasePair))
allocate(AntiSenseDnaResNatom(NumDnaBasePair))
allocate(HistoneAtomId(nHISTONE),H31_AtomId(nHISTONE),H32_AtomId(nHISTONE))
allocate(SenseDnaAtomId(NumDnaBasePair,100))
allocate(AntiSenseDnaAtomId(NumDnaBasePair,100))
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

NpairHistone = 0
NpairH31 = 0
NpairH32 = 0

SenseDnaResNatom(:) = 0
AntiSenseDnaResNatom(:) = 0


do i=1,natom
read(punit,111)AtomName(i),ResName(i),ResID(i),xref(i),yref(i),zref(i)
111   format(12x,a4,1x,a4,1x,i4,4x,3f8.3)

ElementName=adjustl(AtomName(i))


if(ElementName == 'H')mass(i) = 1.0080
if(ElementName == 'C')mass(i) = 12.0110
if(ElementName == 'N')mass(i) = 14.0070
if(ElementName == 'O')mass(i) = 15.9990
if(ElementName == 'S')mass(i) = 32.0600
if(ElementName == 'P')mass(i) = 30.9740

if(ElementName == 'H')cycle

if(i >= ifirst_dna.and.i <= ilast_dna)then
	ResID(i) = ResID(i) - 974
        if(ResID(i) <= NumDnaBasePair)then
                SenseDnaResNatom(ResID(i)) = SenseDnaResNatom(ResID(i)) + 1
                SenseDnaAtomId(ResID(i),SenseDnaResNatom(ResID(i))) = i
        else
                IDBasePair = NumDnaBase + 1 - ResID(i)
                AntiSenseDnaResNatom(IDBasePair) = AntiSenseDnaResNatom(IDBasePair) + 1
                AntiSenseDnaAtomId(IDBasePair,AntiSenseDnaResNatom(IDBasePair)) = i
        endif
endif

if((i >= ifirst_H31.and.i <= ilast_H31).or.(i >= ifirst_H32.and.i <= ilast_H32).or. &      ! H3
  &(i >= ifirst_H41.and.i <= ilast_H41).or.(i >= ifirst_H42.and.i <= ilast_H42).or. &      ! H4
  &(i >= ifirst_H2A1.and.i <= ilast_H2A1).or.(i >= ifirst_H2A2.and.i <= ilast_H2A2).or. &  ! H2A
  &(i >= ifirst_H2B1.and.i <= ilast_H2B1).or.(i >= ifirst_H2B2.and.i <= ilast_H2B2))then   ! H2B 
        NpairHistone = NpairHistone + 1
        HistoneAtomId(NpairHistone) = i

endif

if(i >= ifirst_H31.and.i <= ilast_H31)then    ! H31
        NpairH31 = NpairH31 + 1
        H31_AtomId(NpairH31) = i
endif

if(i >= ifirst_H32.and.i <= ilast_H32)then    ! H32
        NpairH32 = NpairH32 + 1
        H32_AtomId(NpairH32) = i
endif

enddo

write(*,*)'# of selected histone core heavy atom====>',NpairHistone
write(*,*)'# of selected H31 core heavy atom=======>',NpairH31
write(*,*)'# of selected H32 core heavy atom=======>',NpairH32
write(*,*)


call CENTER(natom,NpairHistone,HistoneAtomId,xref,yref,zref,mass)

call COM(NpairHistone,HistoneAtomId,mass,xref,yref,zref,HistoneComX,HistoneComY,HistoneComZ)

!---------------Calculating P_0 vector for the reference structure-----

do ii=1,nSelectedBase
i = SelectedBaseID(ii)

call DNA_BASE_COM(i,NumDnaBasePair,SenseDnaResNatom(i),SenseDnaAtomId,mass,xref,yref,zref,&
        &SenseBaseComX,SenseBaseComY,SenseBaseComZ)
call DNA_BASE_COM(i,NumDnaBasePair,AntiSenseDnaResNatom(i),AntiSenseDnaAtomId,mass,&
        &xref,yref,zref,AntiSenseBaseComX,AntiSenseBaseComY,AntiSenseBaseComZ)

BaseComX = 0.5*(SenseBaseComX + AntiSenseBaseComX)
BaseComY = 0.5*(SenseBaseComY + AntiSenseBaseComY)
BaseComZ = 0.5*(SenseBaseComZ + AntiSenseBaseComZ)

call UNIT_VECTOR(BaseComX,BaseComY,BaseComZ,HistoneComX,HistoneComY&
        &,HistoneComZ,PxBaseRef(i),PyBaseRef(i),PzBaseRef(i))
enddo
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
nread=0
nused=0

time=0.d0
dt=dt*nskip
write(*,*)'dt====>',dt

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
call CENTER(natom,NpairHistone,HistoneAtomId,x,y,z,mass)
call ALIGN(natom,NpairHistone,HistoneAtomId,x,y,z,xref,yref,zref)

call COM(NpairHistone,HistoneAtomId,mass,x,y,z,HistoneComX,HistoneComY,HistoneComZ)

call COM(NpairH31,H31_AtomId,mass,x,y,z,H31_ComX,H31_ComY,H31_ComZ)
call COM(NpairH32,H32_AtomId,mass,x,y,z,H32_ComX,H32_ComY,H32_ComZ)

WX = 0.5*(H31_ComX + H32_ComX)    ! W(WX,WY,WZ) = Middle points of clusters
WY = 0.5*(H31_ComY + H32_ComY)    !              formed by H31 and H32
WZ = 0.5*(H31_ComZ + H32_ComZ)

call UNIT_VECTOR(HistoneComX,HistoneComY,HistoneComZ,WX,WY,WZ&
        &,UX,UY,UZ) ! Dyad Vector U(UX,UY,UZ) = W - Histone COM

call UNIT_VECTOR(H32_ComX,H32_ComY,H32_ComZ,H31_ComX,H31_ComY,H31_ComZ,&
        &WX,WY,WZ)    ! W(WX,WY,WZ) = COM of H32 - COM of H31

call CROSS_PROD(UX,UY,UZ,WX,WY,WZ,FX,FY,FZ) ! F = U X W 

SrotX = 0.d0      ! SrotX = X-compoent of rotational order parameter
SrotY = 0.d0      ! SrotY = Y-compoent of rotational order parameter
Strans = 0.d0     ! Strans = translocation order parameter

do ii=1,nSelectedBase

i = SelectedBaseID(ii)

call DNA_BASE_COM(i,NumDnaBasePair,SenseDnaResNatom(i),SenseDnaAtomId,mass,x,y,z,&
        &SenseBaseComX,SenseBaseComY,SenseBaseComZ)
call DNA_BASE_COM(i,NumDnaBasePair,AntiSenseDnaResNatom(i),AntiSenseDnaAtomId,mass,&
        &x,y,z,AntiSenseBaseComX,AntiSenseBaseComY,AntiSenseBaseComZ)

BaseComX = 0.5*(SenseBaseComX + AntiSenseBaseComX)
BaseComY = 0.5*(SenseBaseComY + AntiSenseBaseComY)
BaseComZ = 0.5*(SenseBaseComZ + AntiSenseBaseComZ)

!---------------Calculating P vector for the simulated structure-----

call UNIT_VECTOR(BaseComX,BaseComY,BaseComZ,HistoneComX,HistoneComY&
        &,HistoneComZ,PxBaseSim(i),PyBaseSim(i),PzBaseSim(i))

call DOT_PROD(PxBaseSim(i),PyBaseSim(i),PzBaseSim(i)&
        &,PxBaseRef(i),PyBaseRef(i),PzBaseRef(i),Angle)

call CROSS_PROD(PxBaseSim(i),PyBaseSim(i),PzBaseSim(i)&
        &,PxBaseRef(i),PyBaseRef(i),PzBaseRef(i),CrosX,CrosY,CrosZ)

call DOT_PROD(CrosX,CrosY,CrosZ,FX,FY,FZ,SignDet)

if(SignDet <= 0.0)then
        Angle = acos(Angle)
else
        Angle = -acos(Angle)
endif


!This is where the sign of the angle has to be decided 

Strans = Strans + Angle

call UNIT_VECTOR(SenseBaseComX,SenseBaseComY,SenseBaseComZ,AntiSenseBaseComX&
        &,AntiSenseBaseComY,AntiSenseBaseComZ,&
        &BxBaseSim(i),ByBaseSim(i),BzBaseSim(i))

call DOT_PROD(PxBaseSim(i),PyBaseSim(i),PzBaseSim(i)&
        &,BxBaseSim(i),ByBaseSim(i),BzBaseSim(i),Angle)

call CROSS_PROD(PxBaseSim(i),PyBaseSim(i),PzBaseSim(i)&
        &,BxBaseSim(i),ByBaseSim(i),BzBaseSim(i),CrosX,CrosY,CrosZ)


!Defining the vector D in 5' to 3' Direction; this D will be used for 
!determining the sign of the rotational location of the selected basepair

if(i <= NumDnaBasePair - NbaseStep + 1)then

ComX5Prime = SenseBaseComX ! Here the base on the sense strand for which 
ComY5Prime = SenseBaseComY ! rotational location to be determined is considered
ComZ5Prime = SenseBaseComZ ! as 5' base for defining the D vector

j = i + NbaseStep-1  ! 3' base is considered as the base on the strand strand 
                     ! NbaseStep - 1 ahead from the selected base for which the 
                     ! rotational location is being calculated

call DNA_BASE_COM(j,NumDnaBasePair,SenseDnaResNatom(j),SenseDnaAtomId,mass,&
        &x,y,z,ComX3Prime,ComY3Prime,ComZ3Prime)

call UNIT_VECTOR(ComX5Prime,ComY5Prime,ComZ5Prime,ComX3Prime,ComY3Prime&
        &,ComZ3Prime,DxStrand,DyStrand,DzStrand)
endif

call DOT_PROD(CrosX,CrosY,CrosZ,DxStrand,DyStrand,DzStrand,SignDet)

if(SignDet <= 0.0)then           ! Positive angle means that minor groove 
        Angle = acos(Angle)      ! of the selected base is oriented away from
else                             ! Histone; whereas the negative angle indicates
        Angle = -acos(Angle)     ! that the minor groove is facing toward 
endif                            ! histone core

SrotX = SrotX + cos(Angle)
SrotY = SrotY + sin(Angle)
enddo

Strans = Strans/(dble(nSelectedBase)*Lambda)
SrotX = SrotX/dble(nSelectedBase)
SrotY = SrotY/dble(nSelectedBase)

Srot = datan(SrotY/SrotX)/pi

time=time+dt
write(ounit,'(3f12.6)')time/1000000.0,Strans,Srot

enddo
enddo

end program OrderParameter
!------------------------------------------------------------------------------
subroutine CROSS_PROD(ax,ay,az,bx,by,bz,cx,cy,cz)   ! C = A X B
real(kind=8),intent(in):: ax,ay,az,bx,by,bz
real(kind=8),intent(out):: cx,cy,cz

cx = ay*bz - by*az                        
cy = bx*az - ax*bz
cz = ax*by - bx*ay

end subroutine CROSS_PROD
!------------------------------------------------------------------------------
subroutine DOT_PROD(aX,ay,aZ,bX,bY,bZ,Angle)   ! A.B = |A|*|B| cos(Angle)
real(kind=8),intent(in):: aX,ay,aZ,bX,bY,bZ    !     = aX*bX + aY*bY + aZ*BZ
real,intent(out):: Angle
real(kind=8):: DotX,DotY,DotZ

DotX = aX*bX
DotY = aY*bY
DotZ = aZ*bZ

Angle = DotX + DotY + DotZ

end subroutine DOT_PROD
!------------------------------------------------------------------------------
subroutine UNIT_VECTOR(aX,aY,aZ,bX,bY,bZ,abX,abY,abZ)  ! AB = B - A
real(kind=8),intent(in):: aX,aY,aZ,bX,bY,bZ
real(kind=8),intent(out):: abX,abY,abZ
real(kind=8):: Modab

abX = bX - aX
abY = bY - aY
abZ = bZ - aZ

Modab = sqrt(abX*abX + abY*abY + abZ*abZ)

abX=abX/Modab
abY=abY/Modab
abZ=abZ/Modab

end subroutine UNIT_VECTOR
!------------------------------------------------------------------------------
subroutine DNA_BASE_COM(i,Nbase,ResNatom,AtomId,mass,x,y,z,ComX,ComY,ComZ)
integer,intent(in):: i,Nbase,ResNatom
integer,dimension(Nbase,*),intent(in):: AtomId
real,dimension(*),intent(in):: x,y,z
real(kind=8),dimension(*),intent(in):: mass
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
real(kind=8),dimension(*),intent(in):: mass
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
subroutine CENTER(npr,np,atp,x,y,z,xm)
real,dimension(*),intent(inout):: x,y,z
real(kind=8),dimension(*),intent(in):: xm
integer,dimension(*),intent(in):: atp
integer,intent(in):: np,npr

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
