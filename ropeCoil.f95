program ropeCoil
	use arcLength
	use zEquation
	use xyEquations
	implicit none
	
	integer													:: i, n
	real(SELECTED_REAL_KIND(15))							:: beta, g, sf = 2, w = 1, omega = 6.28
	real(SELECTED_REAL_KIND(15)), allocatable, dimension(:)	:: s, x, y, z
	
	! Request user input
	print *, "Number of nodes:"
	read *, n
	print *, "Beta value:"
	read *, beta
	print *, "Unitless gravity constant:"
	read *, g
	
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
	end do
	
	print *, "s:"
	print *, s
	
	! Use analytical solution for z
	call solveZ(z, s, beta, g)
	
	print *, "z:"
	print *, z
	
	! Solve x and y by successive over-relaxation
	call solveXY(x, y, s, beta, omega)
	
	print *, "x:"
	print *, x
	print *, "y:"
	print *, y
	
	! Deallocate coordinate arrays
	deallocate(s)
	deallocate(x)
	deallocate(y)
	deallocate(z)
end program ropeCoil