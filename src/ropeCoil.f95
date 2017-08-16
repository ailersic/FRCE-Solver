program ropeCoil
	use arcLength
	use zEquation
	use zEquationSOR
	use xyEquations
	implicit none
	
	integer													:: i, n
	real(SELECTED_REAL_KIND(15))							:: beta, g, sf = 2, w = 1, omega = 6.28
	real(SELECTED_REAL_KIND(15)), allocatable, dimension(:)	:: s, x, y, z
	
	! Request user input
	!print *, "Number of nodes:"
	!read *, n
	!print *, "Beta value:"
	!read *, beta
	!print *, "Unitless gravity constant:"
	!read *, g
	
	n = 10
	beta = 5
	g = 0
	
	! Allocate coordinate arrays with n elements
	allocate(s(n))
	allocate(x(n))
	allocate(y(n))
	allocate(z(n))
	
	! Calculate arc length of fluid rope
	call solveArcLength(sf, beta, g)
	
	! Populate arc length array with linearly spaced nodes
	do i = 1,n
		s(i) = (sf*(i - 1))/(n - 1)
		z(i) = (i-1.0)/(n-1.0)
	end do
	
	! Use analytical solution for z
	!call solveZ(z, s, beta, g)
	
	!do i = 1,n
	!	z(i) = z(i) + 4*1e-2*(i - 1.0)*(n - i)/((n - 1.0)*(n - 1.0))
	!end do
	
	!print *, "z:"
	print *, z
	
	! Use numerical solution for z
	call solveZSOR(z, s, beta)
	
	!print *, "z:"
	!print *, z(n/2+1)
	
	! Solve x and y by successive over-relaxation
	!call solveXY(x, y, s, beta, omega)
	
	!print *, "x:"
	!print *, x
	!print *, "y:"
	!print *, y
	
	! Deallocate coordinate arrays
	deallocate(s)
	deallocate(x)
	deallocate(y)
	deallocate(z)
end program ropeCoil