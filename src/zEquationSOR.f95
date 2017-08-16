module zEquationSOR
contains
	function residualZ(z, s, beta)
		implicit none
		integer										:: n, i
		real(SELECTED_REAL_KIND(15))				:: am2, am1, a0, ap1, ap2, h, residualZ
		real(SELECTED_REAL_KIND(15)), intent(in)	:: beta, s(:), z(:)
		
		n = size(s)
		h = s(2) - s(1)
		
		am2 = -beta/(2*h*h*h)
		am1 = 1/(h*h) + beta/(h*h*h)
		a0 = -2/(h*h)
		ap1 = 1/(h*h) - beta/(h*h*h)
		ap2 = beta/(2*h*h*h)
		
		residualZ = 0
		
		residualZ = residualZ + abs(am2*(z(2) - 2*h) + am1*z(1) + a0*z(2) + ap1*z(3) + ap2*z(4))
		do i=3,n-2
			residualZ = residualZ + abs(am2*z(i-2) + am1*z(i-1) + a0*z(i) + ap1*z(i+1) + ap2*z(i+2))
		end do
		residualZ = residualZ + abs(am2*z(n-3) + am1*z(n-2) + a0*z(n-1) + ap1*z(n) + ap2*z(n-1))
	end

	subroutine solveZSOR(z, s, beta)
		implicit none
		integer										:: n, i, j, iterMax = 1e2
		real(SELECTED_REAL_KIND(15))				:: am2, am1, a0, ap1, ap2, zSum1, zSum2
		real(SELECTED_REAL_KIND(15))				:: h, res, tol = 1e-8, w = 0.1
		real(SELECTED_REAL_KIND(15)), intent(in)	:: beta, s(:)
		real(SELECTED_REAL_KIND(15)), intent(inout)	:: z(:)
		
		! Initialize number of nodes, step size, residual, and iteration counter
		n = size(s)
		h = s(2) - s(1)
		res = 1
		j = 0
		
		am2 = -beta/(2*h*h*h)
		am1 = 1/(h*h) + beta/(h*h*h)
		a0 = -2/(h*h)
		ap1 = 1/(h*h) - beta/(h*h*h)
		ap2 = beta/(2*h*h*h)
		
		do while ((res.gt.tol).and.(iterMax.gt.j))
			! Sum values before iteration
			zSum1 = 0
			do i=1,n
				zSum1 = zSum1 + z(i)
			end do
			
			! Iterate on all internal nodes, with first and last modified for boundary conditions
			z(2) = (1-w)*z(2) + w/(am2 + a0)*(2*h*am2 - am1*z(1) - ap1*z(3) - ap2*z(4))
			do i=3,n-2
				z(i) = (1-w)*z(i) - w/a0*(-am2*z(i-2) - am1*z(i-1) - ap1*z(i+1) - ap2*z(i+2))
			end do
			z(n-1) = (1-w)*z(n-1) + w/(ap2 + a0)*(-am2*z(n-3) - am1*z(n-2) - ap1*z(n))
		
			! Sum values after iteration
			zSum2 = 0
			do i=1,n
				zSum2 = zSum2 + z(i)
			end do
			
			! Calculate residual and increment iteration counter
			res = abs(zSum2 - zSum1)/abs(zSum2)
			j = j + 1
			
			!print *, "ITERATION"
			print *, residualZ(z, s, beta)
		end do
	end
end module zEquationSOR