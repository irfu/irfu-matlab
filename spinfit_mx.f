c Spin fit routime obtained from KHT and matlab interface to it.
c yuri@irfu.se, 2004
c
c $Id$
c


	subroutine mexFunction(nlhs, plhs, nrhs, prhs)
	
c	[bad,x,sigma,iter,lim] = spinfit(fnterms,fitmax,fomega,at,az);
	
	include 'sfit.inc'
	
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on the DEC Alpha
C     64-bit platform
C
	integer plhs(*), prhs(*)
C-----------------------------------------------------------------------
C

	integer nlhs, nrhs
	integer m, n, size
	integer mxGetM, mxGetN, mxGetPr, mxIsNumeric
	real*8 omega
	real*8 tvar
	integer nterms,itmax,iter,lim
	integer nn, ier
	integer pr_omega, pr_nterms, pr_itmax, pr_iter, pr_lim, pr_at, pr_az, pr_bad, pr_x, pr_sigma
c output
	real*8 x(maxterms_fit), sigma	
	
	intrinsic nint
	
c
c check for proper number of arguments
c
	if (nrhs .ne. 5) then
		call mexErrMsgTxt('spinfit requires five input arguments')
	elseif (nlhs .gt. 5) then
		call mexErrMsgTxt('yprime requires five output argument')
	endif

c	nterms	
	m = mxGetM(prhs(1))
	n = mxGetN(prhs(1))
	size = m*n
	if (size .ne. 1) then
		call mexErrMsgTxt('nterms must be scalar')
	endif
	
	pr_nterms = mxGetPr(prhs(1))
	call mxCopyPtrToReal8(pr_nterms, tvar, 1)
	nterms = nint(tvar)
	
c	write(*,*) 'nterms: ',nterms
	
c	itmax	
	m = mxGetM(prhs(2))
	n = mxGetN(prhs(2))
	size = m*n
	if (size .ne. 1) then
		call mexErrMsgTxt('itmax must be scalar')
	endif
	
	pr_itmax = mxGetPr(prhs(2))
	call mxCopyPtrToReal8(pr_itmax, tvar, 1)
	itmax = nint(tvar)	
c	write(*,*) 'itmax: ',itmax
	
c	omega	
	m = mxGetM(prhs(3))
	n = mxGetN(prhs(3))
	size = m*n
	if (size .ne. 1) then
		call mexErrMsgTxt('omega must be scalar')
	endif
	
	pr_omega = mxGetPr(prhs(3))
	call mxCopyPtrToReal8(pr_omega, omega, 1)	

c	write(*,*) 'omega: ', omega

c	time
	pr_at = mxGetPr(prhs(4))
	m = mxGetM(prhs(4))
	n = mxGetN(prhs(4))
	if (m .ne. 1) then
		call mexErrMsgTxt('t must be n x 1 vector')
	endif
	
	nn = n;
	
c	write(*,*) 'nn: ', nn

c	e-field
	m = mxGetM(prhs(5))
	n = mxGetN(prhs(5))
	if (m .ne. 1) then
		call mexErrMsgTxt('e must be n x 1 vector')
	endif
	
	if (n .ne. nn) then
		call mexErrMsgTxt('t and e must be of the same length')
	endif
	
	pr_az = mxGetPr(prhs(5))
	
c	write(*,*) 'All inputs are ready.'
	
	plhs(1) = mxCreateDoubleMatrix(1,nn,0)
	pr_bad = mxGetPr(plhs(1))
	
	plhs(2) = mxCreateDoubleMatrix(1,maxterms_fit,0)
	plhs(3) = mxCreateDoubleMatrix(1,1,0)
	plhs(4) = mxCreateDoubleMatrix(1,1,0)
	plhs(5) = mxCreateDoubleMatrix(1,1,0)
	
c	write(*,*) 'Ready to go.'

	call sfit(nterms,itmax,iter,lim,omega,nn,%val(pr_at),%val(pr_az),%val(pr_bad),x,sigma,ier)
	
	pr_x = mxGetPr(plhs(2))
c10	format(F6.4)
c 	do i = 1,nterms
c		write(*,10) x(i)
c	end do
	call mxCopyReal8ToPtr(x, pr_x, maxterms_fit)
	
	pr_sigma = mxGetPr(plhs(3))
	call mxCopyReal8ToPtr(real(sigma), pr_sigma, 1)
	
	pr_iter = mxGetPr(plhs(4))
	call mxCopyReal8ToPtr(real(iter), pr_iter, 1)
	
	pr_lim = mxGetPr(plhs(5))
	call mxCopyReal8ToPtr(real(lim), pr_lim, 1)
	
	return
	end


c Module name: SFIT
c
c By:
c Bengt H. Nilsson, KTH
c
c Module description:
c
c Modified:
c by Yuri Khotyaintsev for g77 compiler
c 1) replace structures
c 2) use automatic array instead of pointers

	subroutine sfit (fnterms,fitmax,fiter,flim,fomega,nn,at,az,bad,x,sigma,ier)

c Function name: SFIT
c
c Description:
c Fit x(1)+x(2)*cos(omega*t)+x(3)*sin(omega*t)
c         +x(4)*cos(2*omega*t)+x(5)*sin(2*omega*t)+... to data
c
c Constraints:
c
c Interface:
c fit - structure containing fitting control
c nn - number of data points
c at - array of times
c az - array of data
c bad - array indicating bad points
c x - array for resulting coefficients from fit
c sigma - output value
c ier - error indicator
c
c Returns: none

	implicit none

	include 'sfit.inc'
c input
	real*8 fomega
	integer fnterms,fitmax,fiter,flim
	integer nn
	real*8 at(nn),az(nn)
c output
	logical*1 bad(nn)
	real*8 x(maxterms_fit),sigma
	integer ier
c local
	real*8 omega
	real*8 arg,t,y
	real*8 s(maxterms_fit,maxterms_fit+1),w(maxterms_fit+1),
	1	q(maxterms_fit,maxterms_fit+1)
	real*8 diff,ref,const
	integer row,col,i,iter
	logical*1 change

	real*8 adiff(nn)

	real*8 const0,dconst
	data const0/1.4/,dconst/0.4/
                  

	ier = -1
	if (nn .lt. fnterms+1) goto 999
	if (fnterms .gt. maxterms_fit .or.
	1	mod(fnterms,2) .eq. 0) goto 999
	omega = fomega
	
c 20	format(F6.4,' ',F6.4)
c 	do i = 1,nn
c		write(*,20) at(i), az(i)
c	end do

C Build normal equations system
	do row = 1,fnterms
	  do col = 1,fnterms+1
	    s(row,col) = 0.0
	  end do
	end do
C Add to normal equations
	w(1) = 1
	do i = 1,nn                  
	  do row = 2,fnterms,2
	    arg = omega * float(row/2) * at(i)
	    w(row) = cos(arg)
	    w(row+1) = sin(arg)
	  end do
	  w(fnterms+1) = az(i)
	  do row = 1,fnterms
	    do col = row,fnterms+1
	      s(row,col) = s(row,col) + w(row)*w(col)
	    end do
	  end do
	  bad(i) = .false.
	end do
	flim = nn
	const = const0
C Start of iteration loop    
	do iter = 1,fitmax
	  fiter = iter
C Solve normal equations
	  if (flim .lt. fnterms+1) then
	    ier = -1
	    goto 999
	  end if

	  do row = 1,fnterms
	    q(row,row) = s(row,row)
	    do col = row+1,fnterms
	      q(row,col) = s(row,col)
	      q(col,row) = s(row,col)
	    end do                                       
	    q(row,fnterms+1) = s(row,fnterms+1)
	  end do
C Solve
	  call solve (q,fnterms,x,ier)
	  if (ier .ne. 0) then
	    goto 999
	  end if
C Compute sigma
	  sigma = 0.0
	  do i = 1,nn
	    if (.not. bad(i)) then
	      t = at(i)
	      w(1) = 1
	      do row = 2,fnterms,2
	        arg = omega * float(row/2) * t
	        w(row) = cos(arg)
	        w(row+1) = sin(arg)
	      end do
	      y = 0.0
	      do row = 1,fnterms
	        y = y + x(row) * w(row)
	      end do
	      diff = az(i) - y
	      adiff(i) = diff
	      sigma = sigma + diff*diff  
	    end if
	  end do
	  sigma = sqrt(sigma/(flim-fnterms))

	  if (fiter .lt. fitmax) then
	    ref = const*sigma
C Search bad points
	    change = .false.
	    do i = 1,nn
	      if (.not. bad(i)) then
	        if (abs(adiff(i)) .gt. ref) then
C Subtract from normal equations
	          t = at(i)
	          w(1) = 1
	          do row = 2,fnterms,2
	            arg = omega * float(row/2) * t
	            w(row) = cos(arg)
	            w(row+1) = sin(arg)
	          end do
	          w(fnterms+1) = az(i)
	          do row = 1,fnterms
	            do col = row,fnterms+1
	              s(row,col) = s(row,col) - w(row)*w(col)
	            end do
	          end do
	          flim = flim - 1
	          bad(i) = .true.
	          change = .true.
	        end if
	      end if
	    end do                   

	    if (.not. change) go to 999
	    if (flim .le. 1) go to 999
	    const = const + dconst
	  end if
	end do

999	return
	end
	
c Module name: SOLVE
c
c By:
c Bengt H. Nilsson, KTH
c
c Module description:
c
c Modified:

	SUBROUTINE SOLVE (A,N,X,IER)

c Function name: SOLVE
c
c Description:
c equation solver.
c
c Constraints:
c
c Interface:
c A - equation system
c N - number of equations
c X - result
c IER - error indicator
c
c Returns:

	implicit none

	include 'sfit.inc'

	real*8 A(maxterms_fit,maxterms_fit+1),X(maxterms_fit)
	integer n,ier

	integer nm1,np1,j,k,mm,jp1,m,mm1,l
	real*8 work,y             

	IER=0
	NM1=N-1
	NP1=N+1
	DO 70 J=1,NM1
C
C FIND NON ZERO ENTRY IN JTH COLUMN
	DO 10 K=J,N
	MM=K
	IF (ABS(A(K,J)) .GT. 1.E-12) GO TO 20
10	CONTINUE
	IER=J
	RETURN
20	IF (MM .EQ. J) GO TO 40
C
C INTERCHANGE MMTH ROW WITH JTH ROW
	DO 30 K=J,NP1
	WORK=A(MM,K)
	A(MM,K)=A(J,K)
	A(J,K)=WORK
30	CONTINUE
C
C SUBTRACT JTH ROW FROM SUBSEQUENT ROWS
40	JP1=J+1
	DO 60 K=JP1,N
	Y=A(K,J)/A(J,J)
	DO 50 L=J,NP1
	A(K,L)=A(K,L)-Y*A(J,L)
50	CONTINUE
60	CONTINUE
70	CONTINUE
C
C NOW SOLVE
	DO 90 J=1,N
	M=NP1-J  
	X(M)=A(M,NP1)/A(M,M)
	MM1=M-1
	DO 80 K=1,MM1
	A(K,NP1)=A(K,NP1)-X(M)*A(K,M)
80	CONTINUE
90	CONTINUE
	RETURN
	END
