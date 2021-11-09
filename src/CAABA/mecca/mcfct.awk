# ----------------------------------------------------------------------------
#
# Author:
#   Rolf Sander, Max-Planck-Institute, Mainz, Germany, 2008-...
#
# Time-stamp: <2013-01-15 14:29:31 sander>
#
# mcfct.awk adds Monte-Carlo factors to rate coefficients in *.eqn file
#
# - usage:
#   gawk -f mcfct.awk mecca.eqn
#
# mz_rs_20130115+
# Note that the input file must be in unicode format because the § has a
# different encoding in unicode (c2a7) than in iso-latin-1 (a7). If
# necessary, convert it with:
#
# iconv -f iso-8859-1 -t utf-8 input.eqn > output.eqn
# mz_rs_20130115-
#
# ----------------------------------------------------------------------------

BEGIN {
  i = 1
}

# ----------------------------------------------------------------------------

# (replacement text "\\1" works only with gensub but not with gsub):

# note that:
#             exp( ln(f) * mcfct ) = f^mcfct

{

  while (match($0, "{§§[^}]*}") != 0) {
    # logarithmic uncertainty (e.g. IUPAC evaluation):
    $0 = gensub("{§§([^}]*)}", "*EXP(\\1*LOG(10.)*mcexp(" i++ "))", "", $0)
  }

  while (match($0, "{§}") != 0) {
    # default uncertainty assumed to be 1.25:
    $0 = gensub("{§}", "*EXP(LOG(1.25)*mcexp(" i++ "))", "", $0)
  }

  while (match($0, "{§[^}]*}") != 0) {
    # uncertainty factor (e.g. JPL evaluation):
    $0 = gensub("{§([^}]*)}", "*EXP(LOG(\\1)*mcexp(" i++ "))", "", $0)
  }

  printf "%s\n", $0

}

# ----------------------------------------------------------------------------

END{
  # make the number of random numbers available to Fortran code via MAX_MCEXP:
  print  "#INLINE F90_GLOBAL"
  print  "  ! from mcfct.awk:"
  printf "  INTEGER, PARAMETER, PUBLIC :: MAX_MCEXP = %d\n", i-1
  print  "!KPPPP_DIRECTIVE vector variable definition start"
  print  "  REAL :: mcexp(MAX_MCEXP) ! Monte-Carlo factor"
  print  "!KPPPP_DIRECTIVE vector variable definition end"
  print  "#ENDINLINE {above lines go to messy_mecca_kpp_global}"
}

# ----------------------------------------------------------------------------
