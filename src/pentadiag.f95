module pentadiag
contains
	subroutine solvePent(u, s, bv, omega, betaV)
		implicit none
		integer								:: i, n
		real*8								:: h, am2, am1, a0, ap1, ap2
		real*8, intent(in)					:: s(:), bv(:), omega, betaV
		real*8, intent(out)					:: u(:)
		real*8, allocatable, dimension(:)	:: alpha, beta, gamma, mu, a, b, c, d, e, v
		
		n = size(s)
		h = s(2) - s(1)
		
		allocate(a(n-1))
		allocate(b(n-2))
		allocate(c(2:n))
		allocate(d(n))
		allocate(e(3:n))
		
		am2 = -betaV/(2*h*h*h)
		am1 = 1/(h*h) + betaV/(h*h*h)
		a0 = -2/(h*h) - omega*omega
		ap1 = 1/(h*h) - betaV/(h*h*h)
		ap2 = betaV/(2*h*h*h)
		
		d(1) = 1;			a(1) = 0;			b(1) = 0
		c(2) = am1;			d(2) = am2 + a0;	a(2) = ap1;			b(2) = ap2
		do i=3,n-2
			e(i) = am2;		c(i) = am1;			d(i) = a0;			a(i) = ap1;			b(i) = ap2
		end do
		e(n-1) = am2;		c(n-1) = am1;		d(n-1) = ap2 + a0;	a(n-1) = ap1
		e(n) = 0;			c(n) = 0;			d(n) = 1
		
		allocate(alpha(n-1))
		allocate(beta(n-2))
		allocate(gamma(n))
		allocate(mu(n))
		allocate(v(n))
		
		do i=1,n
			if (i.eq.1) then
				mu(i)    = d(i)
				alpha(i) = a(i)/mu(i)
				beta(i)  = b(i)/mu(i)
				
				v(i)     = bv(i)
			elseif (i.eq.2) then
				gamma(i) = c(i)
				mu(i)    = d(i) - alpha(i-1)*gamma(i)
				alpha(i) = (a(i) - beta(i-1)*gamma(i))/mu(i)
				beta(i)  = b(i)/mu(i)
				
				v(i)     = (bv(i) - v(i-1)*gamma(i))/mu(i)
			elseif (i.eq.(n-1)) then
				gamma(i) = c(i) - alpha(i-2)*e(i)
				mu(i)    = d(i) - beta(i-2)*e(i) - alpha(i-1)*gamma(i)
				alpha(i) = (a(i) - beta(i-1)*gamma(i))/mu(i)
				
				v(i)     = (bv(i) - v(i-2)*e(i) - v(i-1)*gamma(i))/mu(i)
			elseif (i.eq.n) then
				gamma(i) = c(i) - alpha(i-2)*e(i)
				mu(i)    = d(i) - beta(i-2)*e(i) - alpha(i-1)*gamma(i)
				
				v(i)     = (bv(i) - v(i-2)*e(i) - v(i-1)*gamma(i))/mu(i)
			else
				gamma(i) = c(i) - alpha(i-2)*e(i)
				mu(i)    = d(i) - beta(i-2)*e(i) - alpha(i-1)*gamma(i)
				alpha(i) = (a(i) - beta(i-1)*gamma(i))/mu(i)
				beta(i)  = b(i)/mu(i)
				
				v(i)     = (bv(i) - v(i-2)*e(i) - v(i-1)*gamma(i))/mu(i)
			end if
		end do
		
		do i=n,1,-1
			if (i.eq.n) then
				u(n) = v(n)
			elseif (i.eq.(n-1)) then
				u(n-1) = v(n-1) - alpha(n-1)*u(n)
			else
				u(i) = v(i) - alpha(i)*u(i+1) - beta(i)*u(i+2)
			end if
		end do
		
		deallocate(a)
		deallocate(b)
		deallocate(c)
		deallocate(d)
		deallocate(e)
		
		deallocate(alpha)
		deallocate(beta)
		deallocate(gamma)
		deallocate(mu)
		deallocate(v)
	end
	
	function solveZPent(s, beta)
		implicit none
		integer								:: i, n
		real*8								:: h, omega = 0
		real*8								:: az = 0, azs = 1, bz = 1, bzs = 0
		real*8, intent(in)					:: s(:), beta
		real*8, dimension(:), allocatable	:: bl(:), solveZPent(:)
		
		n = size(s)
		h = s(2) - s(1)
		
		allocate(bl(n))
		allocate(solveZPent(n))
		
		bl(1) = az
		bl(2) = -beta/(h*h)*azs
		do i=3,n-1
			bl(i) = 0
		end do
		bl(n-1) = bzs
		bl(n) = bz
		
		call solvePent(solveZPent, s, bl, omega, beta)
		
		deallocate(bl)
	end
	
	function solveXPent(s, y, beta, omega)
		implicit none
		integer								:: i, n
		real*8								:: h
		real*8, intent(in)					:: s(:), y(:), beta, omega
		real*8								:: ax = 0, axs = 0, bx, bxs = 0
		real*8, dimension(:), allocatable	:: bl(:), solveXPent(:)
		
		n = size(s)
		h = s(2) - s(1)
		
		bx = abs(1/omega)
		
		allocate(bl(n))
		allocate(solveXPent(n))
		
		bl(1) = ax
		do i=2,n-1
			bl(i) = omega*y(i+1)/h - omega*y(i-1)/h
		end do
		bl(2) = bl(2) - beta/(h*h)*axs
		bl(n-1) = bl(n-1) - beta/(h*h)*bxs
		bl(n) = bx
		
		call solvePent(solveXPent, s, bl, omega, beta)
		
		deallocate(bl)
	end
	
	function solveYPent(s, x, beta, omega)
		implicit none
		integer								:: i, n
		real*8								:: h
		real*8, intent(in)					:: s(:), x(:), beta, omega
		real*8								:: ay = 0, ays = 0, by = 0, bys = 1
		real*8, dimension(:), allocatable	:: bl(:), solveYPent(:)
		
		n = size(s)
		h = s(2) - s(1)
		
		allocate(bl(n))
		allocate(solveYPent(n))
		
		bl(1) = ay
		do i=2,n-1
			bl(i) = omega*x(i-1)/h - omega*x(i+1)/h
		end do
		bl(2) = bl(2) - beta/(h*h)*ays
		bl(n-1) = bl(n-1) - beta/(h*h)*bys
		bl(n) = by
		
		call solvePent(solveYPent, s, bl, omega, beta)
		
		deallocate(bl)
	end
end module pentadiag