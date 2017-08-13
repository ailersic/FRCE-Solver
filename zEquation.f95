module zEquation
contains
	subroutine solveZ(z, s, beta, g)
		implicit none
		integer										:: n, i
		real(SELECTED_REAL_KIND(15))				:: beta, g, c1, c2, c3, sf
		real(SELECTED_REAL_KIND(15)), intent(in)	:: s(:)
		real(SELECTED_REAL_KIND(15)), intent(out)	:: z(:)
		
		! Initialize number of nodes
		n = size(s)
		
		! Initialize coefficients in equation
		sf = s(n)
		c1 = (g*sf + 1)/(beta*(1 - exp(-sf/beta)))
		c2 = 1 - (g*sf + 1)/(1 - exp(-sf/beta))
		c3 = beta*(g*sf + 1)/(1 - exp(-sf/beta))
		
		! Solve equation for each node
		do i = 1,n
			z(i) = -c1*beta*beta*exp(-s(i)/beta) + g*s(i)*s(i)/2 + c2*s(i) + c3
		end do
	end
end module zEquation