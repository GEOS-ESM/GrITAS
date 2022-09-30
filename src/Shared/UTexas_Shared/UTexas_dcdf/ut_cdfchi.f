      program ut_cdfchi
!     can use this:
!     https://www.fourmilab.ch/rpkp/experiments/analysis/chiCalc.html
!     to corroborate this program
!     Can also compare w/ Matlab function: chi2inv
!     R. Todling, November 2014
      DOUBLE PRECISION bound,df,p,q,x,t
      integer status,which
      integer iarg,argc
      character(len=256)  argv
      argc = iargc()
      if ( argc < 1 ) then
         print *, "Usage: ut_cdfchi.x q df"
         print *, "   "
         print *, " where "
         print *, "  q   is the probability "
         print *, "  df  is the number of degrees of freedom"
         print *, "   "
         call exit(1)
      endif

      iarg=1
      call GetArg ( iarg, argv )
      read(argv,*) q
      iarg=iarg+1
      call GetArg ( iarg, argv )
      read(argv,*) df
!
      which = 2
      p= 1.d0-q
      write(*,'(a,f8.3)') 'With a probability of', real(q,4)
      write(*,'(a,i10)')  'and the number of degrees of freedom: ', nint(df)
      call cdfchi (which,p,q,x,df,status,bound)
      if (status==0) then
          write(*,'(a,f16.7)') 'the chi-square value is: ', x
      else
         call exit(1)
      endif
      write(*,'(a,f8.3)') 'With a probability of', real(q,4)
      write(*,'(a,i10)')  'and the number of degrees of freedom: ', nint(df)
      call cdft   (which,p,q,t,df,status,bound)
      if (status==0) then
          write(*,'(a,f16.7)') 'the t-student value is: ', t
      else
         call exit(1)
      endif
      end
