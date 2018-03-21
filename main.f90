program main 

	use mpi
	use utils

	implicit none
	integer :: ierr, rank, num_cores, k, i, j, num2receive, tag
	real(kind=kind(0.0d0)), parameter :: pi = 4*atan(1.0d0)
	real(kind=kind(0.0d0)) :: M, T, S, PM, PT, PS, HM, HT, HS, stime, etime
	real(kind=kind(0.0d0)), dimension(:), allocatable :: mypts, allpts
	integer, dimension(mpi_status_size) :: mystatus
	integer, dimension(:), allocatable :: num2send, displs

	call mpi_init(ierr)
	call mpi_comm_rank(mpi_comm_world,rank,ierr)
	call mpi_comm_size(mpi_comm_world,num_cores,ierr)
	stime = mpi_wtime()

	open(unit=13,file="m.dat",action="write",status="replace")
	open(unit=14,file="t.dat",action="write",status="replace")
	open(unit=15,file="s.dat",action="write",status="replace")
	
	if (rank == 0) then
		write(*,*) "Using", num_cores, "cores."
		write(*,"(a5, a30, a30, a30)") "K", "Midpoint_H(K)", "Trapezoidal_H(K)", "Simpson's_H(K)"
	endif
	do k = 100, 10000
		allocate(allpts(0:k))	
		allocate(mypts(0:k/num_cores))
		allocate(num2send(0:num_cores-1))
		allocate(displs(0:num_cores-1))

		allpts = (/(j*pi/k,j=0,k)/)
		mypts = (/(-1.0,j=0,k/num_cores)/)
		displs = (/(0,j=0,num_cores-1)/)
		do i=0,num_cores-1
			num2send(i) = k/num_cores
			if (i<=modulo(k,num_cores)) then
				num2send(i) = num2send(i)+1
			endif
			if (i /= 0) then
				displs(i) = displs(i-1) + num2send(i)
			endif
		enddo
		num2receive = num2send(rank)
	
		call mpi_scatterv(allpts,num2send,displs,mpi_double_precision,&
						  mypts,num2receive,mpi_double_precision,0,&
						  mpi_comm_world,ierr)

		PM = 0; PT = 0; PS = 0
		do i=0,size(mypts)-1
			M = 0; T = 0; S =0
			if (mypts(i) > -1) then
				call mts(k,(/(mypts(i)+j*pi/100/k,j=0,100)/),101,M,T,S)
				PM = PM + M
				PT = PT + T
				PS = PS + S
			endif
		enddo
		deallocate(allpts,mypts,num2send,displs)
		
		call mpi_barrier(mpi_comm_world, ierr)

		tag = 0			
		if (rank == 0) then
			HM = PM
			HT = PT
			HS = PS
			do i = 1, num_cores-1
				call mpi_recv(PM,1,mpi_double_precision,i,tag,mpi_comm_world,mystatus,ierr)
				call mpi_recv(PT,1,mpi_double_precision,i,tag,mpi_comm_world,mystatus,ierr)
				call mpi_recv(PS,1,mpi_double_precision,i,tag,mpi_comm_world,mystatus,ierr)
				HM = HM + PM
				HT = HT + PT
				HS = HS + PS
			enddo
		else 
			call mpi_send(PM,1,mpi_double_precision,0,tag,mpi_comm_world,ierr)
			call mpi_send(PT,1,mpi_double_precision,0,tag,mpi_comm_world,ierr)
			call mpi_send(PS,1,mpi_double_precision,0,tag,mpi_comm_world,ierr)
		endif
		
		if (rank == 0) then
			if (modulo(k, 1000) == 0) then
				write(*,"(i5, es30.15e3, es30.15e3, es30.15e3)") k, HM, HT, HS
			endif
			write(13,*) k, HM
			write(14,*) k, HT
			write(15,*) k, HS
		endif
	enddo

	close(13); close(14); close(15)
	
	etime = mpi_wtime()

	if (rank == 0) then
		write(*,*) "Time taken: ", etime-stime
	endif

	call mpi_finalize(ierr)
	
endprogram main
