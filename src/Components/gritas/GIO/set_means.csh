  set Usage     = "set_means.csh SubSetNo [-r result]"
  set Example   = "set_means.csh 2 "

  if ( $#argv < 1 ) goto err
  set SubSetNo  = $argv[1]; shift argv

  if      ( $SubSetNo == 1   ) then
     set Years  = ( 1979 )
#     set Months = ( 01 02 03 04 05 06 07 08 09 10 11 12 )
     set Months = (                      08 09 10 11 12 )
  else if ( $SubSetNo == 2   ) then 
#    set Years  = ( 1980 1981 1982 )
#    set Months = ( 01 02 03 04 05 06 07 08 09 10 11 12 )
     set Years  = (           1982 )
     set Months = (                                  12 )
  else if ( $SubSetNo == 3   ) then
     set Years  = ( 1983 1984 1985 )
     set Months = ( 01 02 03 04 05 06 07 08 09 10 11 12 )
  else if ( $SubSetNo == 3_5a ) then
     set Years  = ( 1986           )
     set Months = ( 01 02 03 04 05 06                   )
  else if ( $SubSetNo == 3_5b ) then
     set Years  = ( 1986           )
     set Months = (                   07 08 09 10 11 12 )
  else if ( $SubSetNo == 3_5c ) then
     set Years  = (      1987      )
     set Months = ( 01 02 03 04 05 06                   )
  else if ( $SubSetNo == 3_5d ) then
     set Years  = (      1987      )
     set Months = (                   07 08 09 10 11 12 )
  else if ( $SubSetNo == 3_5e ) then
     set Years  = (           1988 )
     set Months = ( 01 02 03 04 05 06                   )
  else if ( $SubSetNo == 3_5f ) then
     set Years  = (           1988 )
     set Months = (                   07 08 09 10 11 12 )
  else if ( $SubSetNo == 4   ) then
     set Years  = ( 1989 1990 1991 )
     set Months = ( 01 02 03 04 05 06 07 08 09 10 11 12 )
  else if ( $SubSetNo == 5   ) then
     set Years  = ( 1992 1993 1994 )
     set Months = ( 01 02 03 04 05 06 07 08 09 10 11 12 )
  else if ( $SubSetNo == 6   ) then
     set Years  = ( 1995 1996 1997 )
     set Months = ( 01 02 03 04 05 06 07 08 09 10 11 12 )
  else if ( $SubSetNo == 7   ) then
     set Years  = ( 1998 1999 2000 )
     set Months = ( 01 02 03 04 05 06 07 08 09 10 11 12 )
  else if ( $SubSetNo == 7a  ) then
     set Years  = ( 1998           )
     set Months = (                            10 11 12 )
  else if ( $SubSetNo == 7b  ) then
     set Years  = (      1999 2000 )
     set Months = ( 01 02 03 04 05 06 07 08 09 10 11 12 )
  else if ( $SubSetNo == 8   ) then
     set Years  = ( 2001 2002 2003 )
     set Months = ( 01 02 03 04 05 06 07 08 09 10 11 12 )
  else if ( $SubSetNo == 9   ) then
     set Years  = ( 2004 2005 2006 )
     set Months = ( 01 02 03 04 05 06 07 08 09 10 11 12 )
   else
     echo " Invalid value for for SubSetNo (=$SubSetNo )" 
     exit 1
  endif

  foreach Year  ( $Years  )
  foreach Month ( $Months )
     csh -vx gritas2means.csh ${Year}${Month} -r rms
     if ( $status ) then
        echo " Error status returned from gritas2means (-r = rms)" 
        echo "    Year/Mon = ${Year}${Month}"
        exit 1
     endif
     csh -vx gritas2means.csh ${Year}${Month} -r obrate
     if ( $status ) then
        echo " Error status returned from gritas2means (-r = obrate)" 
        echo "    Year/Mon = ${Year}${Month}"
        exit 1
     endif
  end
  end

# All is well
# ----------- 
  exit 0

# Usage messages
# --------------
  err:
     echo "  usage: $Usage"
     echo "example: $Example"
  exit 1
