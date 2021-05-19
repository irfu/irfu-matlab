C Spin fit routine obtained from KTH and matlab interface to it.
C yuri@irfu.se, 2007
C
C This is a new 64-bit compatible version
C
C compile with the following flags for the best speed:
C mex -v c_efw_spinfit_mx.F \
C        FFLAGS='$FFLAGS -O2 -mtune=opteron -funroll-loops'
C 
C for gfortran under Linux, use the -ffixed-line-length-none flag.
C
C
#include "fintrf.h"
#define MAXTERMS_FIT 9
#define NaN -159e7
#ifndef mwSize
# define mwSize integer*4
#endif


      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

      implicit none
C      [ts,sfit,sdev,iter,nout] = ...
C                  c_efw_spinfit_mx(maxit,minpts,nterms,te,data,pha);
C
C maxit  - maximum number of iterations (zero means infinity)
C minpts - minimum fraction of spin to perform fit (1=all points)
C nterms - number of terms to fit
C t0     - start of the first spin
C te     - EFW time in seconds (epoch)
C data   - EFW data from pair in mV/m, should correspond to te
C pha    - phase of the current boom pair in radians, should correspond
C          to te


      mwPointer plhs(*), prhs(*)
C-----------------------------------------------------------------------
C

      integer nlhs, nrhs, mxIsNumeric
#if MX_HAS_INTERLEAVED_COMPLEX
      mwPointer mxGetDoubles
#else
      mwPointer mxGetPr
#endif
      mwPointer mxCreateDoubleMatrix
      mwSize mxGetM, mxGetN
      mwSize m, n, size

      mwPointer pr_tmp, pr_te, pr_data, pr_pha
      integer maxit, minpts, nterms
C local
      real*8 tvar
      integer ndata, nspins
      integer*8 tstart
      mwPointer pr_ts, pr_sfit, pr_sdev, pr_iter, pr_nout

C      set default rounding be to nearest integer
      intrinsic nint

C      write(*,*) 'starting ...'

C
C check for proper number of arguments
C
      if (nrhs .ne. 6) then
        call mexErrMsgTxt('c_efw_spinfit_mx needs six input args')
      elseif (nlhs .gt. 5) then
        call mexErrMsgTxt('c_efw_spinfit_mx needs five output args')
      endif

C      maxit            argument #1
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      size = m*n
      if (size .ne. 1) then
        call mexErrMsgTxt('MAXIT must be scalar')
      endif
#if MX_HAS_INTERLEAVED_COMPLEX
      pr_tmp = mxGetDoubles(prhs(1))
#else
      pr_tmp = mxGetPr(prhs(1))
#endif
      call mxCopyPtrToReal8(pr_tmp, tvar, 1)
      maxit = nint(tvar)
      if (maxit .lt. 1) then
        call mexErrMsgTxt('MAXIT must be a positive nonzero number')
      endif
C      write(*,*) 'maxit: ',maxit

C      nterms            argument #3
      m = mxGetM(prhs(3))
      n = mxGetN(prhs(3))
      size = m*n
      if (size .ne. 1) then
        call mexErrMsgTxt('NTERMS must be scalar')
      endif
#if MX_HAS_INTERLEAVED_COMPLEX
      pr_tmp = mxGetDoubles(prhs(3))
#else
      pr_tmp = mxGetPr(prhs(3))
#endif
      call mxCopyPtrToReal8(pr_tmp, tvar, 1)
      nterms = nint(tvar)
      if (nterms .le. 1 .or. nterms .gt. MAXTERMS_FIT .or.
     +    mod(nterms,2) .eq. 0) then
        call mexErrMsgTxt('NTERMS must be one of 3,5,7,9')
      endif
C      write(*,*) 'nterms: ',nterms

C      minpts            argument #2
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      size = m*n
      if (size .ne. 1) then
        call mexErrMsgTxt('MINPTS must be scalar')
      endif
#if MX_HAS_INTERLEAVED_COMPLEX
      pr_tmp = mxGetDoubles(prhs(2))
#else
      pr_tmp = mxGetPr(prhs(2))
#endif
      call mxCopyPtrToReal8(pr_tmp, tvar, 1)
      minpts = nint(tvar)
      if (minpts .le. nterms) then
        call mexErrMsgTxt('MINPTS must be larger than NTERMS')
      endif
C      write(*,*) 'minpts: ',minpts

C      te            argument #4
      m = mxGetM(prhs(4))
      ndata = mxGetN(prhs(4))
      if (m .ne. 1) then
        call mexErrMsgTxt('TE must be n x 1 vector')
      endif
#if MX_HAS_INTERLEAVED_COMPLEX
      pr_te = mxGetDoubles(prhs(4))
#else
      pr_te = mxGetPr(prhs(4))
#endif

C      e-field            argument #5
      m = mxGetM(prhs(5))
      n = mxGetN(prhs(5))
      if (m .ne. 1) then
        call mexErrMsgTxt('DATA must be n x 1 vector')
      endif
#if MX_HAS_INTERLEAVED_COMPLEX
      pr_data = mxGetDoubles(prhs(5))
#else
      pr_data = mxGetPr(prhs(5))
#endif
      if (n .ne. ndata) then
        call mexErrMsgTxt('TE and DATA must be of the same length')
      endif

C      phase            argument #6
      m = mxGetM(prhs(6))
      n = mxGetN(prhs(6))
      if (m .ne. 1) then
        call mexErrMsgTxt('phase must be n x 1 vector')
      endif
#if MX_HAS_INTERLEAVED_COMPLEX
      pr_pha = mxGetDoubles(prhs(6))
#else
      pr_pha = mxGetPr(prhs(6))
#endif
      if (n .ne. ndata) then
        call mexErrMsgTxt('TE and PHA must be of the same length')
      endif

C      compute number of spins
      call getTStartNSpins(%val(pr_te),ndata,tstart,nspins)
C      write(*,*) 'nspins: ', nspins
C      write(*,*) 'tstart: ', tstart
C      write(*,*) 'All inputs are ready.'

C
C prepare output
C

      plhs(1) = mxCreateDoubleMatrix(1,nspins,0)
      plhs(2) = mxCreateDoubleMatrix(nterms,nspins,0)
      plhs(3) = mxCreateDoubleMatrix(1,nspins,0)
      plhs(4) = mxCreateDoubleMatrix(1,nspins,0)
      plhs(5) = mxCreateDoubleMatrix(1,nspins,0)

#if MX_HAS_INTERLEAVED_COMPLEX
      pr_ts = mxGetDoubles(plhs(1))
      pr_sfit = mxGetDoubles(plhs(2))
      pr_sdev  = mxGetDoubles(plhs(3))
      pr_iter = mxGetDoubles(plhs(4))
      pr_nout = mxGetDoubles(plhs(5))
#else
      pr_ts = mxGetPr(plhs(1))
      pr_sfit = mxGetPr(plhs(2))
      pr_sdev  = mxGetPr(plhs(3))
      pr_iter = mxGetPr(plhs(4))
      pr_nout = mxGetPr(plhs(5))
#endif
C      write(*,*) 'Ready to go.'

      call spinfit( maxit,minpts,nterms,tstart,
     + nspins,ndata,%val(pr_te),%val(pr_data),
     + %val(pr_pha),%val(pr_ts),%val(pr_sfit),%val(pr_sdev),
     + %val(pr_iter),%val(pr_nout))

C      write(*,*) 'returning...'

      return
      end

C
C Compure number of spins
C
      subroutine getTStartNSpins(te,ndata,tstart,nspins)
      implicit none
C in
      integer ndata
      real*8 te(ndata)
C out
      integer nspins
      integer*8 tstart

C      write(*,*) 'te(1):', te(1)
C      write(*,*) 'te(end):', te(ndata)

      tstart = te(1)/4
      tstart = tstart*4
C      write(*,*) 'tstart:', tstart
      nspins = (te(ndata) - tstart)/4 + 1

      return
      end

C
C call onesfit for every spin
C
      subroutine spinfit(maxit,minpts,nterms,t0,nspins,ndata,
     +  te,az,pha,ts,sfit,sdev,iter,nout)
      implicit none
C in
      integer maxit, minpts, nterms, nspins, ndata
      integer*8 t0
      real*8 te(ndata), az(ndata), pha(ndata)
C out
      real*8 ts(nspins),sfit(nterms,nspins),sdev(nspins),iter(nspins)
      real*8 nout(nspins)
C local
      integer idx, nn, i, j, idxs, idxe, ierr, niter, nbad
      integer lim, spinstart
      real*8 x(MAXTERMS_FIT), sigma
      real*8 dt

C initialize arrays
C      write(*,*) 'spinfit: init arrays'
      do i = 1, nspins
C        ts(i) =  t0 +(i-1)*4.0 +2.0
        dt = (i-1)*4.0 +2.0
        ts(i) = (i-1)*4.0 +2.0 + dble(t0)
C        write(*,*) 'i: ',i,' dt:', dt,' ts(i)-t0:',ts(i)-t0
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

C      check if we have at least one spin of data
      if (ndata .lt. minpts) return

      idx = 1
      do i=1, nspins
C        write(*,*) 'spin: ',i
        spinstart =  t0 +(i-1)*4
        idxs = 0
        idxe = 0

        do while (idxe .eq. 0)
C        write(*,*) 'spin: ',i,' idx:',idx,' left: ',spinstart+4-te(idx)

C        first point alredy more then one spin later then the first one
          if ( (idxs .eq. 0) .and. (te(idx) .ge. (spinstart+4) ) ) then
            idxe = -1
          else
            if (idxs .eq. 0 .and. te(idx) .ge. spinstart) then
              idxs = idx
            endif
            if (idx .eq. ndata .or. te(idx+1) .ge. spinstart+4) then
              idxe = idx
            endif
            idx = idx + 1
          endif
        end do

C        check number of data points
        nn = idxe-idxs+1
C        write(*,*) 'got: ',nn,' min:',minpts
        if ( nn .gt. minpts) then
C          write(*,*) 'start: ',idxs,' stop:',idxe
          call onesfit(nterms,maxit,niter,lim,nn,pha(idxs),
     + az(idxs),nbad,x,sigma,ierr)
          if (ierr .eq. 0) then
C            write(*,*) '         spin: ok'
            do j=1, nterms
              sfit(j,i) = x(j)
C              write(*,*) 'x(j):',x(j)
            end do
            sdev(i) = sigma
            iter(i) = real(niter)
            nout(i) = nbad
          endif
        endif
      end do

      return
      end

C Module name: SFIT
C
C By:
C Bengt H. Nilsson, KTH
C
C Module description:
C
C Modified:
C by Yuri Khotyaintsev for g77 compiler
C 1) replace structures
C 2) use automatic array instead of pointers

      subroutine onesfit (fnterms,fitmax,fiter,flim,
     + nn,apha,az,nbad,x,sigma,ier)

C Function name: SFIT
C
C Description:
C Fit x(1)+x(2)*cos(pha)+x(3)*sin(pha)
C         +x(4)*cos(2*pha)+x(5)*sin(2*pha)+... to data
C
C Constraints:
C
C Interface:
C fit - structure containing fitting control
C nn - number of data points
C apha - array of phases
C az - array of data
C bad - array indicating bad points
C x - array for resulting coefficients from fit
C sigma - output value
C ier - error indicator
C
C Returns: none

      implicit none

C input
      real*8 fomega
      integer fnterms,fitmax,fiter,flim
      integer nn
      real*8 apha(nn),az(nn)
C output
      logical*1 bad(nn)
C      real*8 rbad(nn)
      real*8 x(MAXTERMS_FIT),sigma
      integer ier, nbad
C local
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

C      write(*,*) 'onesfit E:'
C10    format(F6.4)
C      do i = 1,nn
C            write(*,10) az(i)
C      end do
C      write(*,*) 'onesfit PHA:'
C      do i = 1,nn
C            write(*,*) apha(i)
C      end do

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

999   do i = 1,nn
        if (bad(i)) then
          nbad = nbad +1
C          rbad(i) = 1.
C          write(*,*) 'bad: ', i
C        else
C          rbad(i) = 0.
        end if
      end do

      return
      end

C Module name: SOLVE
C
C By:
C Bengt H. Nilsson, KTH
C
C Module description:
C
C Modified:

      SUBROUTINE SOLVE (A,N,X,IER)

C Function name: SOLVE
C
C Description:
C equation solver.
C
C Constraints:
C
C Interface:
C A - equation system
C N - number of equations
C X - result
C IER - error indicator
C
C Returns:

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
10    CONTINUE
      IER=J
      RETURN
20    IF (MM .EQ. J) GO TO 40
C
C INTERCHANGE MMTH ROW WITH JTH ROW
      DO 30 K=J,NP1
      WORK=A(MM,K)
      A(MM,K)=A(J,K)
      A(J,K)=WORK
30    CONTINUE
C
C SUBTRACT JTH ROW FROM SUBSEQUENT ROWS
40    JP1=J+1
      DO 60 K=JP1,N
      Y=A(K,J)/A(J,J)
      DO 50 L=J,NP1
      A(K,L)=A(K,L)-Y*A(J,L)
50    CONTINUE
60    CONTINUE
70    CONTINUE
C
C NOW SOLVE
      DO 90 J=1,N
      M=NP1-J
      X(M)=A(M,NP1)/A(M,M)
      MM1=M-1
      DO 80 K=1,MM1
      A(K,NP1)=A(K,NP1)-X(M)*A(K,M)
80    CONTINUE
90    CONTINUE
      RETURN
      END
