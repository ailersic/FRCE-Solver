module arcLength
contains
	subroutine solveArcLength(sf, beta, g)
		implicit none
		real(SELECTED_REAL_KIND(15))				:: sf, beta, g, lerr, a = 1, tol = 1e-12

		! While error is above tolerance, do Newton's method
		do while (abs(sfEqn(sf, beta, g)) > tol)
			sf = sf - sfEqn(sf, beta, g)/dsfEqn(sf, beta, g)
		end do
	end

	function sfEqn(sf, beta, g)
		implicit none
		real(SELECTED_REAL_KIND(15)), intent(in)	:: sf, beta, g
		real(SELECTED_REAL_KIND(15))				:: sfEqn
		
		sfEqn = beta*(g*sf + 1) + g*sf*sf/2 + sf*(1 - (g*sf + 1)/(1 - exp(-sf/beta))) - 1
	end

	function dsfEqn(sf, beta, g)
		implicit none
		real(SELECTED_REAL_KIND(15)), intent(in)	:: sf, beta, g
		real(SELECTED_REAL_KIND(15))				:: a, b, dsfEqn
		
		a = sf*exp(sf/beta)*(g*sf + 1) - beta*(exp(sf/beta) - 1)*(2*g*sf*exp(sf/beta) + 1)
		b = beta*(exp(sf/beta) - 1)*(exp(sf/beta) - 1)
		dsfEqn = beta*g + g*sf + a/b
	end
end module arcLength