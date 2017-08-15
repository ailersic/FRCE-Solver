module zEquationSOR
contains
	subroutine SORIterZ(z, s, beta, w)
		implicit none
		integer										:: n, i
		real(SELECTED_REAL_KIND(15))				:: h
		real(SELECTED_REAL_KIND(15)), intent(in)	:: beta, w, s(:)
		real(SELECTED_REAL_KIND(15)), intent(inout)	:: z(:)
		
		! Initialize number of nodes and step size
		n = size(s)
		h = s(2) - s(1)
		
		! Iterate on all internal nodes, with first and last modified for boundary conditions
		z(2) = (1-w)*z(2) + 2*w*h*h*h/(4*h + beta)*(beta/(h*h) + (1/(h*h) - beta/(h*h*h))*z(3) + beta/(2*h*h*h)*z(4))
		do i=3,n-2
			z(i) = (1-w)*z(i) - w*h*h/2*(beta/(2*h*h*h)*z(i-2) - (1/(h*h) + beta/(h*h*h))*z(i-1) &
				- (1/(h*h) - beta/(h*h*h))*z(i+1) - beta/(2*h*h*h)*z(i+2))
		end do
		z(n-1) = (1-w)*z(n-1) + 2*w*h*h*h/(beta - 4*h)*(1/(h*h) - beta/(h*h*h) - beta/(2*h*h*h)*z(n-3) - (1/(h*h) + beta/(h*h*h))*z(n-2))
	end
	
	function residualZ(z, s, beta)
		implicit none
		integer										:: n, i
		real(SELECTED_REAL_KIND(15))				:: h, residualZ
		real(SELECTED_REAL_KIND(15)), intent(in)	:: beta, s(:), z(:)
		
		n = size(s)
		h = s(2) - s(1)
	
		residualZ = 0
		residualZ = residualZ + abs(beta/(h*h) - (2/(h*h) + beta/(2*h*h*h))*z(2) + (1/(h*h) - beta/(h*h*h))*z(3) + beta/(2*h*h*h)*z(4))
		do i=3,n-2
			residualZ = residualZ + abs(-beta/(2*h*h*h)*z(i-2) + (1/(h*h) + beta/(h*h*h))*z(i-1) - 2/(h*h)*z(i) &
				+ (1/(h*h) - beta/(h*h*h))*z(i+1) + beta/(2*h*h*h)*z(i+2))
		end do
		residualZ = residualZ + &
			abs(-beta/(2*h*h*h)*z(n-3) + (1/(h*h) + beta/(h*h*h))*z(n-2) + (beta/(2*h*h*h) - 2/(h*h))*z(n-1) + 1/(h*h) - beta/(h*h*h))
	end
	
	subroutine solveZSOR(z, s, beta)
		use zEquation
		implicit none
		integer										:: n, i
		real(SELECTED_REAL_KIND(15))				:: beta, w, rf, g = 0
		real(SELECTED_REAL_KIND(15)), intent(in)	:: s(:)
		real(SELECTED_REAL_KIND(15)), intent(out)	:: z(:)
		
		! Initialize SOR factor, final radius, and number of nodes
		w = 1
		n = size(s)
		
		! Inialize z with analytical solution
		call solveZ(z, s, beta, g)
		
		print *, 'z:'
		print *, z
		print *, 'residual:'
		print *, residualZ(z, s, beta)
		
		! Solve x and y by SOR
		do i=1,10 ! Dummy number of iterations
			call SORIterZ(z, s, beta, w)
			print *, 'z:'
			print *, z
			print *, 'residual:'
			print *, residualZ(z, s, beta)
		end do
	end
end module zEquationSOR