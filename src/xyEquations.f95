module xyEquations
contains
	subroutine SORIterX(x, y, s, beta, omega, w)
		implicit none
		integer										:: n, i
		real(SELECTED_REAL_KIND(15))				:: h, am2, am1, a0, ap1, ap2, bm1, bp1
		real(SELECTED_REAL_KIND(15)), intent(in)	:: beta, w, omega, s(:)
		real(SELECTED_REAL_KIND(15)), intent(inout)	:: x(:), y(:)
		
		! Initialize number of nodes and step size
		n = size(s)
		h = s(2) - s(1)
		
		! Calculate coefficients on nodes
		am2 = -beta/(2*h*h*h)
		am1 = 1/(h*h) + beta/(h*h*h)
		a0 = -2/(h*h) - omega*omega
		ap1 = 1/(h*h) - beta/(h*h*h)
		ap2 = beta/(2*h*h*h)
		
		bm1 = -omega/h
		bp1 = omega/h
		
		! Iterate on all internal nodes, with first and last modified for boundary conditions
		x(2) = a0/(a0 + w*am2)*((1-w)*x(2) + w/a0*(bp1*y(3) - bm1*y(1) - ap1*x(3) - ap2*x(4)))
		do i=3,n-2
			x(i) = (1-w)*x(i) + w/a0*(bp1*y(i+1) - bm1*y(i-1) - am2*x(i-2) - am1*x(i-1) - ap1*x(i+1) - ap2*x(i+2))
		end do
		x(n-1) = (1-w)*x(n-1) + w/a0*(bp1*y(n) - bm1*y(n-2) - am2*x(n-3) - am1*x(n-2) - ap1*x(n) - ap2*x(n-1))
	end
	
	subroutine SORIterY(x, y, s, beta, omega, w)
		implicit none
		integer										:: n, i
		real(SELECTED_REAL_KIND(15))				:: h, am2, am1, a0, ap1, ap2, bm1, bp1
		real(SELECTED_REAL_KIND(15)), intent(in)	:: beta, w, omega, s(:)
		real(SELECTED_REAL_KIND(15)), intent(inout)	:: x(:), y(:)
		
		! Initialize number of nodes and step size
		n = size(s)
		h = s(2) - s(1)
		
		! Calculate coefficients on nodes
		am2 = -beta/(2*h*h*h)
		am1 = 1/(h*h) + beta/(h*h*h)
		a0 = -2/(h*h) - omega*omega
		ap1 = 1/(h*h) - beta/(h*h*h)
		ap2 = beta/(2*h*h*h)
		
		bm1 = omega/h
		bp1 = -omega/h
		
		! Iterate on all internal nodes, with first and last modified for boundary conditions
		y(2) = a0/(a0 + w*am2)*((1-w)*y(2) + w/a0*(bp1*x(3) - bm1*x(1) - ap1*y(3) - ap2*y(4)))
		do i=3,n-2
			y(i) = (1-w)*y(i) + w/a0*(bp1*x(i+1) - bm1*x(i-1) - am2*y(i-2) - am1*y(i-1) - ap1*y(i+1) - ap2*y(i+2))
		end do
		y(n-1) = (1-w)*y(n-1) + w/a0*(bp1*x(n) - bm1*x(n-2) - am2*y(n-3) - am1*y(n-2) - ap2*(y(n-1) + 2*h))
	end
	
	subroutine solveXY(x, y, s, beta, omega)
		implicit none
		integer										:: n, i
		real(SELECTED_REAL_KIND(15))				:: beta, w, omega, rf
		real(SELECTED_REAL_KIND(15)), intent(in)	:: s(:)
		real(SELECTED_REAL_KIND(15)), intent(out)	:: x(:), y(:)
		
		! Initialize SOR factor, final radius, and number of nodes
		w = 1.5
		rf = 1/omega
		n = size(s)
		
		! Inialize values of x and y
		do i=1,n
			x(i) = 0
			y(i) = 0
		end do
		x(n) = rf
		
		! Solve x and y by SOR
		do i=1,50 ! Dummy number of iterations
			call SORIterX(x, y, s, beta, omega, w)
			call SORIterY(x, y, s, beta, omega, w)
		end do
	end
end module xyEquations