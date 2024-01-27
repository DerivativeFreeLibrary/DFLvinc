!===========================================================
! Problem n. 31 in Hock and Schittkowski collection
!===========================================================

subroutine setdim(n,m)
	implicit none
	integer	:: n ,m

	!----------------------------------------------
	! set the number of variables
	!----------------------------------------------
	n = 3

	!----------------------------------------------
	! set the number of inequality constraints
	! Be aware that an equality constraint h(x) = 0
	! must be written as h(x) <= 0 and -h(x) <= 0
	!----------------------------------------------
	m = 1

	return

end subroutine setdim

subroutine funob(n, x, f)
	implicit none
	integer		:: n
	real*8		:: x(n), f

	!----------------------------------------------------
	! this routine computes the objective function value
	!----------------------------------------------------

	f = 9.d0*x(1)**2.d0 + x(2)**2.d0 + 9.d0*x(3)**2.d0

	return
end subroutine funob

subroutine fconstr(n, m, x, constr)
	implicit none
	integer		:: n, m
	real*8		:: x(n), constr(m)

	integer		:: i

	!------------------------------------------------------
	! this routine computes the constraint function values
	!------------------------------------------------------

	constr(1) = 1.d0 - x(1)*x(2)

	return
end subroutine fconstr

subroutine startp(n,x)
	implicit none
	integer		n
	real*8		x(n)

	!----------------------------------------------
	! this routine sets the starting point
	!----------------------------------------------
	
	X(1) = 1.D0
	X(2) = 1.D0
	X(3) = 1.D0

	return
end subroutine startp

subroutine setbounds(n,lb,ub)
	implicit none
	integer		n
	real*8		lb(n), ub(n)

	!----------------------------------------------
	! this routine sets the bounds on the variables
	!----------------------------------------------

	lb(1) = -10.0d0  
	ub(1) = 10.0d0
	lb(2) = 1.0d0  
	ub(2) = 10.0d0
	lb(3) = -10.0d0
	ub(3) = 1.0d0

	return
end subroutine setbounds

subroutine which_integer(n,A,B,is_integer,step,x)
!----------------------------------------------------
! this routine defines which variables are discrete,
! sets the smallest step between two consecutive 
! discrete values for all discrete variables and, 
! possibly, changes the starting point to make it 
! satisfy the integrality constraint
!----------------------------------------------------
	implicit none
	integer, intent(IN)	:: n
	real*8,  intent(INOUT)	:: A(n), B(n)
	integer, intent(OUT)	:: is_integer(n)
	real*8,  intent(OUT)	:: step(n)
	real*8,  intent(INOUT)	:: x(n)
	integer			:: i
	
	is_integer = 0
	step       = 0.d0
    
	is_integer(2) = 1
	step(2)       = ( B(2) - A(2) ) / 20.0d0
	x(2)          = ( A(2) + B(2) ) /  2.0d0
	
	return

end subroutine which_integer
