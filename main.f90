program testDFLgen
	implicit none

	integer :: n, m
	integer	:: nf, nf_max, iprint, istop
	real*8  :: fob, alfatol
	real	:: tbegin, tend
	logical :: varscaling

	integer, allocatable  :: index_int(:)
	real*8,  allocatable  :: x(:),bl(:),bu(:),step(:)      

!----------------------------------------------
!   set problem dimensions 
!   n  : number of variables
!   m  : number of inequality constraints
!----------------------------------------------
	call setdim(n,m)

!----------------------------------------------
!   allocate necessary vectors
!----------------------------------------------
	allocate(index_int(n))
	allocate(x(n),bl(n),bu(n),step(n))

!----------------------------------------------
!   set variable bounds
!   bl  : array of lower bounds
!   bu  : array of upper bounds
!----------------------------------------------
	call setbounds(n,bl,bu)
!----------------------------------------------
!   set the initial (starting) point
!----------------------------------------------
	call startp(n,x)
!----------------------------------------------
!   define the integer variables and their step
!----------------------------------------------
	call which_integer(n,bl,bu,index_int,step,x)

!-----------------------------------------------------------------------
!  varscaling : variables scaling
!  		ON   if varscaling = .TRUE.
!  		OFF  if varscaling = .FALSE.
!
!  alfatol    : convergence tolerance
!  nf_max     : maximum allowed number of f.evaluations
!  iprint     : printing verbosity
!-----------------------------------------------------------------------
	varscaling  = .false.
	alfatol     = 1.d-6
	nf_max      = 5000
	iprint      = -1

	call cpu_time(tbegin)

	write(*,*)
	write(*,*) ' call the optimizer ...'

	call mixed_optimizer(n,m,index_int,step,varscaling,x,fob,bl,bu,alfatol,nf_max,nf,iprint,istop)

!--------------------------------------------------------
! Upon normal termination istop is equal to 0
! an istop value > 0 signals an abnormal termination
!--------------------------------------------------------
	write(*,*) '... optimizer ended with istop =',istop
	write(*,*)

	call cpu_time(tend)
	
	write(*,*) ' total time:',tend-tbegin

end program testDFLgen
