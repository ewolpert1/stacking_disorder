program graphite_Ising
use iso_fortran_env, only: dp=>real64, int64
implicit none
character(len=100) :: line
character(len=1024) :: restart_path
real(dp) :: start_T=1000.0_dp, floor_T=0.0_dp
integer :: seed,E3,E3flip,AcceptedQ,AcceptedT,Rejected,Accepted0,AcceptedE,E4,E4flip
integer a,n,nn,atest,axtest,nn1,minus,iter2,iostatus,pos
integer(int64) :: iter
real(dp) J1,E1,E2,E1flip,E2flip,Energychange,i,j,k,x,y,z,H,Etot,m,rn1,rn2,Q,T,x1,x2,y1,y2,Etotm,J5,a_lat_pat,c_lat_pat
integer,allocatable:: coords(:,:),coordsnew(:,:)
character (len=90) :: filename
integer :: nmoves,nstages,accepted_stage
real(dp) :: Hcheck,drift,cooling_ratio
nmoves=1000000;nstages=135;cooling_ratio=.95_dp
call get_command_argument(1,line);read(line,*)seed
call get_command_argument(2,line)
if(len_trim(line)>0)read(line,*)nmoves
call get_command_argument(3,line)
if(len_trim(line)>0)read(line,*)nstages
call get_command_argument(4,line)
if(len_trim(line)>0)read(line,*)cooling_ratio
if(cooling_ratio<=0.or.cooling_ratio>=1)error stop 'Invalid cooling ratio'
call get_command_argument(5,restart_path)
call get_command_argument(6,line)
if(len_trim(line)>0)read(line,*)start_T
call get_command_argument(7,line)
if(len_trim(line)>0)read(line,*)floor_T
if(start_T<=0.or.floor_T<0.or.floor_T>start_T)error stop 'Invalid temperatures'
call Randomseed(seed)

! Read in the parameters for the simulations from the file 'parameters_Ising.txt'
open(unit=10, file='parameters_Ising.txt', status='old', action='read', iostat=iostatus)
    if (iostatus /= 0) then
        print *, 'Error opening file.'
        stop
    endif

    do a = 1, 4
        read(10, '(A)') line
        pos = scan(line, "=")  ! Find the position of the equals sign
        select case(a)
            case(1)
                read(line(pos+1:), *) n
            case(2)
                read(line(pos+1:), *) J1
            case(3)
                read(line(pos+1:), *) a_lat_pat
            case(4)
                read(line(pos+1:), *) c_lat_pat
        end select
    end do

    close(10)

    print *, 'Read n = ', n
    print *, 'Read J1 = ', J1
    print *, 'Read a_lat_pat = ', a_lat_pat
    print *, 'Read c_lat_pat = ', c_lat_pat

m=n

allocate(coords(n,1),coordsnew(n,1))


! Initalise random system
do a=1,n
	call random_number(rn1)
	if (rn1<0.5_dp) then
		coords(a,1)=-1
	else
		coords(a,1)=1
	endif
enddo


if(len_trim(restart_path)>0)then
 open(23,file=trim(restart_path),access='stream',form='unformatted',status='old',action='read')
 read(23)coords
 close(23)
endif
open(23,file='initial_state.bin',access='stream',form='unformatted',status='replace')
write(23)coords
close(23)
open(21,file='stats.csv',status='replace')
write(21,'(A)')'stage,T,J1,J2,energy_K,energy_check_K,energy_drift_K,acceptance'
open(22,file='states.bin',access='stream',form='unformatted',status='replace')

T=start_T


! Calculate energy

H=0
do a=1,n
	nn=a+1
	if (nn>n) then
		nn=nn-n
	endif
	H=J1*coords(a,1)*coords(nn,1)+H
enddo
Etot=H
open(23,file='initial_energy_K.txt',status='replace')
write(23,'(ES25.16)')H
close(23)

iter=0
iter2=0
AcceptedQ=0
Accepted0=0
AcceptedE=0
do
	iter=iter+1
	! Pick random slip direction to change
	call random_number(rn1)
		atest=rn1*n+1
		axtest=1
	coordsnew(atest,:)=coords(atest,:)
	coordsnew(atest,axtest)=-coords(atest,axtest)
	nn=atest+1
	if (nn>n) then
		nn=nn-n
	endif
	nn1=atest-1
	if (nn1<1) then
		nn1=nn1+n
	endif

	E1=0
	E2=0
	E3=0
	E4=0
	E4flip=0
	E2flip=0
	E1flip=0

	! Calculate the current energy
	E1=J1*coords(atest,1)*coords(nn,1)+J1*coords(atest,1)*coords(nn1,1)

	! Calculate the energy of the new state
	E1flip=J1*coordsnew(atest,1)*coords(nn,1)+J1*coordsnew(atest,1)*coords(nn1,1)

	! Calculate energy change
	Energychange=E1flip-E1

	! Metropolis criterion
	Q=exp(-max(Energychange,0.0_dp)/T)
	call random_number(rn1)
	if (Energychange<0) then
		Etot=Etot+Energychange
		coords(atest,axtest)=coordsnew(atest,axtest)
		AcceptedE=AcceptedE+1
	elseif (Energychange==0) then
		coords(atest,axtest)=coordsnew(atest,axtest)
		Accepted0=Accepted0+1
	elseif (Q>rn1) then
		Etot=Etot+Energychange
		coords(atest,axtest)=coordsnew(atest,axtest)
		AcceptedQ=AcceptedQ+1
	endif
	if (mod(iter,int(nmoves,int64))==0) then
		iter2=iter2+1

        call check_energy(coords,n,J1,0.0_dp,Hcheck)
        drift=Etot-Hcheck
        accepted_stage=AcceptedQ+Accepted0+AcceptedE
        write(21,'(I0,7(",",ES25.16))')iter2,T,J1,0.0_dp,Etot,Hcheck,drift,real(accepted_stage,dp)/nmoves
        flush(21)
        write(22)coords
        flush(22)
        if(abs(drift)>1.0e-5_dp*max(1.0_dp,abs(Hcheck)))error stop 'Energy drift exceeds tolerance'

		print*,'Number of energetically unfavourable moves accepted = ',AcceptedQ
		print*,'Number of moves accepted with no energy change = ',Accepted0
		print*,'Number of energetically favourable moves accepted = ',AcceptedE
		print*,'Total Energy = ',Etot
		print*,'Temperature = ',T
		AcceptedQ=0
		Accepted0=0
		AcceptedE=0
		! Write out file for the configuration
		write(filename,fmt = '(A12,F12.6,A4)')'phase_Ising_',T,'.cif'
		open(10, file = filename, form='formatted')
		write(10,'(A25,F12.6)') ' data_ Ising_temperature_',T
		write(10,*) '_audit_creation_method         "stacking_disorder Fortran; deterministic seeded revision"'
		write(10,*) '_cell_length_a                   ',a_lat_pat
		write(10,*) '_cell_length_b                   ',a_lat_pat
		write(10,*) '_cell_length_c                   ',c_lat_pat*n
		write(10,*) '_cell_angle_alpha               90.0000(0)'
		write(10,*) '_cell_angle_beta                90.0000(0)'
		write(10,*) '_cell_angle_gamma               120.0000(0)'

		write(10,*) '_symmetry_space_group_name_H-M     "P 1"'
		write(10,*) '_symmetry_Int_Tables_number         1'
		write(10,*) '_symmetry_cell_setting             triclinic'
		write(10,*) 'loop_'
		write(10,*) '_symmetry_equiv_pos_as_xyz'
		write(10,*) '+x,+y,+z'

		write(10,*) 'loop_'
		write(10,*) '_atom_site_label'
		write(10,*) '_atom_site_type_symbol'
		write(10,*) '_atom_site_occupancy'
		write(10,*) '_atom_site_fract_x'
		write(10,*) '_atom_site_fract_y'
		write(10,*) '_atom_site_fract_z'


		m=n
		x=(1.0_dp/3.0_dp)
		y=(2.0_dp/3.0_dp)
		do a=1,n
			if (coords(a,1)==1) then
				x=x+(1.0_dp/3.0_dp)
				y=y+(2.0_dp/3.0_dp)
				write(10,*)'C	C	1.0		',x,y,a/m
				write(10,*)'C	C	1.0		',x-(1.0_dp/3.0_dp),y+(1.0_dp/3.0_dp),a/m
			else
				x=x-(1.0_dp/3.0_dp)
				y=y-(2.0_dp/3.0_dp)
				write(10,*)'C	C	1.0		',x,y,a/m
				write(10,*)'C	C	1.0		',x-(1.0_dp/3.0_dp),y+(1.0_dp/3.0_dp),a/m
			endif
		enddo
		close(10)
		T=max(floor_T,T*cooling_ratio)
	endif
	if (iter2==nstages) then
		print*,'Final total Energy from MC = ',Etot
	exit
	endif
enddo


! This is a check to make sure the Monte Carlo is calculating the energy properly
H=0
do a=1,n
	nn=a+1
	if (nn>n) then
		nn=nn-n
	endif
	H=J1*coords(a,1)*coords(nn,1)+H
enddo
print*,'Final total Energy calculated = ',H
close(21);close(22)
write(*,'(A)')'COMPLETED'
end program

! This ensures that everytime you run a simulation you get new random numbers. For trouble shooting it is best to comment this section out

subroutine RandomSeed(x)
use iso_fortran_env, only: int64
implicit none
integer,intent(in) :: x
integer :: i,n
integer,allocatable :: seed(:)
call random_seed(size=n)
allocate(seed(n))
do i=1,n
 seed(i)=int(modulo(int(x,int64)+104729_int64*int(i,int64),2147483646_int64))
enddo
call random_seed(put=seed)
end subroutine RandomSeed

subroutine check_energy(coords,n,j1,j2,h)
use iso_fortran_env,only:dp=>real64
implicit none
integer,intent(in)::n
integer,intent(in)::coords(n,1)
real(dp),intent(in)::j1,j2
real(dp),intent(out)::h
real(dp)::dx,dy,pi
integer::a,b,c
pi=acos(-1.0_dp);h=0.0_dp
do a=1,n
 b=modulo(a,n)+1
 h=h+j1*coords(a,1)*coords(b,1)
enddo
end subroutine check_energy
