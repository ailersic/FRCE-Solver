program ropeCoil
	use arcLength
	use pentadiag
	use zEquation
	use xyEquations
	implicit none
	
	integer								:: i, n
	real*8								:: beta, g, sf = 2, w = 0.01, omega = 6
	real*8, allocatable, dimension(:)	:: s, x, y, z, xt, yt
	
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
	
	! Calculate arc length of fluid rope
	call solveArcLength(sf, beta, g)
	
	! Populate arc length array and initial x vector with linearly spaced nodes
	do i = 1,n
		s(i) = (sf*(i - 1))/(n - 1)
		x(i) = (i - 1)/(abs(omega)*(n - 1))
	end do
	
	! Use analytical solution for z
	!call solveZ(z, s, beta, g)
	
	! Use numerical solution for z
	z = solveZPent(s, beta)
	
	! Solve x and y by successive over-relaxation
	!call solveXY(x, y, s, beta, omega)
	
	! Use numerical solution for x and y
	do i=1,1
		yt = solveYPent(s, x, beta, omega)
		y = (1-w)*y + w*yt
		deallocate(yt)
		
		xt = solveXPent(s, y, beta, omega)
		x = (1-w)*x + w*xt
		deallocate(xt)
	end do
	
	print *, "x:"
	print *, x
	print *, "y:"
	print *, y
	print *, "z:"
	print *, z
	
	! Deallocate coordinate arrays
	deallocate(s)
	deallocate(x)
	deallocate(y)
	deallocate(z)
end program ropeCoil