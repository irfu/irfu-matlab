c Spin fit routine obtained from KTH and matlab interface to it.
c yuri@irfu.se, 2007
c
c This is a new 64-bit compatible version
c
c compile with the following flags for the best speed:
c mex -v c_efw_spinfit_mx.F FFLAGS='$FFLAGS -O2 -mtune=opteron -funroll-loops'
c 
c for gfortran under Linux, use the -ffixed-line-length-none flag.
c
c $Id$
c
#include "fintrf.h"
#define MAXTERMS_FIT 9
#define NaN -159e7
#ifndef mwsize
# define mwsize integer*4
#endif


	subroutine mexFunction(nlhs, plhs, nrhs, prhs)
	
	implicit none
c	[ts,sfit,sdev,iter,nout] = ...
c                  c_efw_spinfit_mx(maxit,minpts,nterms,te,data,pha);
c
c maxit  - maximum number of iterations (zero means infinity)
c minpts - minimum fraction of spin to perform fit (1=all points)
c nterms - number of terms to fit
c t0     - start of the first spin
c te     - EFW time in seconds (epoch)
c data   - EFW data from pair in mV/m, should correspond to te
c pha    - phase of the current boom pair in radians, should correspond to te

	
	mwpointer plhs(*), prhs(*)
C-----------------------------------------------------------------------
C

	integer nlhs, nrhs, mxIsNumeric
	mwpointer mxGetPr, mxCreateDoubleMatrix
	mwsize mxGetM, mxGetN
	mwsize m, n, size
	
	mwpointer pr_tmp, pr_te, pr_data, pr_pha
	integer maxit, minpts, nterms
c local
	real*8 tvar
	integer ndata, nspins
	integer*8 tstart
	mwpointer pr_ts, pr_sfit, pr_sdev, pr_iter, pr_nout
	
	intrinsic nint
	
c	write(*,*) 'starting ...'

c
c check for proper number of arguments
c
	if (nrhs .ne. 6) then
	call mexErrMsgTxt('c_efw_spinfit_mx needs six input args')
	elseif (nlhs .gt. 5) then
	call mexErrMsgTxt('c_efw_spinfit_mx needs five output args')
	endif

c	maxit		argument #1
	m = mxGetM(prhs(1))
	n = mxGetN(prhs(1))
	size = m*n
	if (size .ne. 1) then
		call mexErrMsgTxt('MAXIT must be scalar')
	endif
	pr_tmp = mxGetPr(prhs(1))
	call mxCopyPtrToReal8(pr_tmp, tvar, 1)
	maxit = nint(tvar)
	if (maxit .lt. 1) then
		call mexErrMsgTxt('MAXIT must be a positive nonzero number')
	endif
c	write(*,*) 'maxit: ',maxit

c	nterms		argument #3
	m = mxGetM(prhs(3))
	n = mxGetN(prhs(3))
	size = m*n
	if (size .ne. 1) then
		call mexErrMsgTxt('NTERMS must be scalar')
	endif
	pr_tmp = mxGetPr(prhs(3))
	call mxCopyPtrToReal8(pr_tmp, tvar, 1)
	nterms = nint(tvar)
	if (nterms .le. 1 .or. nterms .gt. MAXTERMS_FIT .or.
     +	mod(nterms,2) .eq. 0) then
		call mexErrMsgTxt('NTERMS must be one of 3,5,7,9')
	endif
c	write(*,*) 'nterms: ',nterms

c	minpts		argument #2
	m = mxGetM(prhs(2))
	n = mxGetN(prhs(2))
	size = m*n
	if (size .ne. 1) then
		call mexErrMsgTxt('MINPTS must be scalar')
	endif
	pr_tmp = mxGetPr(prhs(2))
	call mxCopyPtrToReal8(pr_tmp, tvar, 1)
	minpts = nint(tvar)
	if (minpts .le. nterms) then
		call mexErrMsgTxt('MINPTS must be larger than NTERMS')
	endif
c	write(*,*) 'minpts: ',minpts
	
c	te		argument #4
	m = mxGetM(prhs(4))
	ndata = mxGetN(prhs(4))
	if (m .ne. 1) then
		call mexErrMsgTxt('TE must be n x 1 vector')
	endif
	pr_te = mxGetPr(prhs(4))
	
c	e-field		argument #5
	m = mxGetM(prhs(5))
	n = mxGetN(prhs(5))
	if (m .ne. 1) then
		call mexErrMsgTxt('DATA must be n x 1 vector')
	endif
	pr_data = mxGetPr(prhs(5))
	if (n .ne. ndata) then
		call mexErrMsgTxt('TE and DATA must be of the same length')
	endif
	
c	phase		argument #6
	m = mxGetM(prhs(6))
	n = mxGetN(prhs(6))
	if (m .ne. 1) then
		call mexErrMsgTxt('phase must be n x 1 vector')
	endif
	pr_pha = mxGetPr(prhs(6))
	if (n .ne. ndata) then
		call mexErrMsgTxt('TE and PHA must be of the same length')
	endif
	
c	compute number of spins
	call getTStartNSpins(%val(pr_te),ndata,tstart,nspins)
c	write(*,*) 'nspins: ', nspins
c	write(*,*) 'tstart: ', tstart
c	write(*,*) 'All inputs are ready.'

c
c prepare output
c

	plhs(1) = mxCreateDoubleMatrix(1,nspins,0)
	plhs(2) = mxCreateDoubleMatrix(nterms,nspins,0)
	plhs(3) = mxCreateDoubleMatrix(1,nspins,0)
	plhs(4) = mxCreateDoubleMatrix(1,nspins,0)
	plhs(5) = mxCreateDoubleMatrix(1,nspins,0)
	
	pr_ts = mxGetPr(plhs(1))
	pr_sfit = mxGetPr(plhs(2))
	pr_sdev  = mxGetPr(plhs(3))
	pr_iter = mxGetPr(plhs(4))
	pr_nout = mxGetPr(plhs(5))
c	write(*,*) 'Ready to go.'

	call spinfit( maxit,minpts,nterms,tstart,
     + nspins,ndata,%val(pr_te),%val(pr_data),
     + %val(pr_pha),%val(pr_ts),%val(pr_sfit),%val(pr_sdev),
     + %val(pr_iter),%val(pr_nout))
	
c	write(*,*) 'returning...'
	
	return
	end

c
c Compure number of spins
c
	subroutine getTStartNSpins(te,ndata,tstart,nspins)
	implicit none
c in
	integer ndata
	real*8 te(ndata)
c out
	integer nspins
	integer*8 tstart
	
c	write(*,*) 'te(1):', te(1)
c	write(*,*) 'te(end):', te(ndata)
	
	tstart = te(1)/4
	tstart = tstart*4
c	write(*,*) 'tstart:', tstart
	nspins = (te(ndata) - tstart)/4 + 1
	
	return
	end

c
c call onesfit for every spin
c
	subroutine spinfit(maxit,minpts,nterms,t0,nspins,ndata,
     +	te,az,pha,ts,sfit,sdev,iter,nout)
	implicit none
c in
	integer maxit, minpts, nterms, nspins, ndata
	integer*8 t0
	real*8 te(ndata), az(ndata), pha(ndata)
c out
	real*8 ts(nspins),sfit(nterms,nspins),sdev(nspins),iter(nspins)
	real*8 nout(nspins)
c local
	integer idx, nn, i, j, idxs, idxe, ierr, niter, nbad
	integer lim, spinstart
	real*8 x(MAXTERMS_FIT), sigma
	real*8 dt

c initialize arrays
c	write(*,*) 'spinfit: init arrays'
	do i = 1, nspins
c		ts(i) =  t0 +(i-1)*4.0 +2.0
		dt = (i-1)*4.0 +2.0
		ts(i) = (i-1)*4.0 +2.0 + dble(t0)
c		write(*,*) 'i: ',i,' dt:', dt,' ts(i)-t0:',ts(i)-t0
	end do
	do i = 1, nspins
		do j = 1, nterms
			sfit(j,i) = NaN
		end do
	end do
	do i = 1, nspins
		sdev(i) = NaN
	end do
	do i = 1, nspins
		iter(i) = NaN
	end do
	do i = 1, nspins
		nout(i) = NaN
	end do
	
c	check if we have at least one spin of data
	if (ndata .lt. minpts) return
	
	idx = 1
	do i=1, nspins
c		write(*,*) 'spin: ',i
		spinstart =  t0 +(i-1)*4
		idxs = 0
		idxe = 0
		
		do while (idxe .eq. 0)
c			write(*,*) 'spin: ',i,' idx:',idx,' left: ',spinstart+4-te(idx)

c			first point alredy more then one spin later then the first one
			if ( (idxs .eq. 0) .and. (te(idx) .ge. (spinstart+4) ) ) then
				idxe = -1
			else
				if (idxs .eq. 0 .and. te(idx) .ge. spinstart) idxs = idx
				if (idx .eq. ndata .or. te(idx+1) .ge. spinstart+4) idxe = idx
				idx = idx + 1
			endif
		end do
		
c		check number of data points
		nn = idxe-idxs+1
c		write(*,*) 'got: ',nn,' min:',minpts
		if ( nn .gt. minpts) then
c			write(*,*) 'start: ',idxs,' stop:',idxe
			call onesfit(nterms,maxit,niter,lim,nn,pha(idxs),az(idxs),nbad,
     + x,sigma,ierr)
			if (ierr .eq. 0) then
c			write(*,*) '         spin: ok'
				do j=1, nterms
					sfit(j,i) = x(j)
c					write(*,*) 'x(j):',x(j)
				end do
				sdev(i) = sigma
				iter(i) = real(niter)
				nout(i) = nbad
			endif
		endif
	end do
	
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

	subroutine onesfit (fnterms,fitmax,fiter,flim,
     + nn,apha,az,nbad,x,sigma,ier)

c Function name: SFIT
c
c Description:
c Fit x(1)+x(2)*cos(pha)+x(3)*sin(pha)
c         +x(4)*cos(2*pha)+x(5)*sin(2*pha)+... to data
c
c Constraints:
c
c Interface:
c fit - structure containing fitting control
c nn - number of data points
c apha - array of phases
c az - array of data
c bad - array indicating bad points
c x - array for resulting coefficients from fit
c sigma - output value
c ier - error indicator
c
c Returns: none

	implicit none

c input
	real*8 fomega
	integer fnterms,fitmax,fiter,flim
	integer nn
	real*8 apha(nn),az(nn)
c output
	logical*1 bad(nn)
c	real*8 rbad(nn)
	real*8 x(MAXTERMS_FIT),sigma
	integer ier, nbad
c local
	real*8 omega
	real*8 arg,pha,y
	real*8 s(MAXTERMS_FIT,MAXTERMS_FIT+1),w(MAXTERMS_FIT+1),
     + q(MAXTERMS_FIT,MAXTERMS_FIT+1)
	real*8 diff,ref,const
	integer row,col,i,iter
	logical*1 change

	real*8 adiff(nn)

	real*8 const0,dconst
	data const0/1.4/,dconst/0.4/
                  
	nbad = 0
	ier = -1
	if (nn .lt. fnterms+1) goto 999
	if (fnterms .gt. MAXTERMS_FIT .or.
     + mod(fnterms,2) .eq. 0) goto 999

c	write(*,*) 'onesfit E:'
c10	format(F6.4)
c	do i = 1,nn
c		write(*,10) az(i)
c	end do
c	write(*,*) 'onesfit PHA:'
c	do i = 1,nn
c		write(*,*) apha(i)
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
	    arg = float(row/2) * apha(i)
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
	      pha = apha(i)
	      w(1) = 1
	      do row = 2,fnterms,2
	        arg = float(row/2) * pha
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
	  sigma = sqrt(sigma/float(flim-1))

	  if (fiter .lt. fitmax) then
	    ref = const*sigma
C Search bad points
	    change = .false.
	    do i = 1,nn
	      if (.not. bad(i)) then
	        if (abs(adiff(i)) .gt. ref) then
C Subtract from normal equations
	          pha = apha(i)
	          w(1) = 1
	          do row = 2,fnterms,2
	            arg = float(row/2) * pha
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
	
999	do i = 1,nn
		if (bad(i)) then
			nbad = nbad +1
c			rbad(i) = 1.
c			write(*,*) 'bad: ', i
c		else
c			rbad(i) = 0.
		end if
	end do
	
	return
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

	real*8 A(MAXTERMS_FIT,MAXTERMS_FIT+1),X(MAXTERMS_FIT)
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
