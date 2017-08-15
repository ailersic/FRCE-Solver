module zEquationSOR
contains
	subroutine solveZSOR(z, s, beta)
		implicit none
		integer										:: n, i, j, iterMax = 10
		real(SELECTED_REAL_KIND(15))				:: zSum1, zSum2, h, res, tol = 1e-8, w = 1
		real(SELECTED_REAL_KIND(15)), intent(in)	:: beta, s(:)
		real(SELECTED_REAL_KIND(15)), intent(inout)	:: z(:)
		
		! Initialize number of nodes, step size, residual, and iteration counter
		n = size(s)
		h = s(2) - s(1)
		res = 1
		j = 0
		
		do while ((res.gt.tol).and.(iterMax.gt.j))
			! Sum values before iteration
			zSum1 = 0
			do i=1,n
				zSum1 = zSum1 + z(i)
			end do
			
			! Iterate on all internal nodes, with first and last modified for boundary conditions
			z(2) = (1-w)*z(2) + 2*w*h*h*h/(4*h + beta)*(beta/(h*h) + (1/(h*h) - beta/(h*h*h))*z(3) + beta/(2*h*h*h)*z(4))
			do i=3,n-2
				z(i) = (1-w)*z(i) - w*h*h/2*(beta/(2*h*h*h)*z(i-2) - (1/(h*h) + beta/(h*h*h))*z(i-1) &
					- (1/(h*h) - beta/(h*h*h))*z(i+1) - beta/(2*h*h*h)*z(i+2))
			end do
			z(n-1) = (1-w)*z(n-1) + 2*w*h*h*h/(beta - 4*h)*(1/(h*h) - beta/(h*h*h) - beta/(2*h*h*h)*z(n-3) - (1/(h*h) + beta/(h*h*h))*z(n-2))
		
			! Sum values after iteration
			zSum2 = 0
			do i=1,n
				zSum2 = zSum2 + z(i)
			end do
			
			! Calculate residual and increment iteration counter
			res = abs(zSum2 - zSum1)/abs(zSum2)
			j = j + 1
		end do
	end
end module zEquationSOR