# ----------------------------------------------------------------------------
#
# Author:
#   Rolf Sander, Max-Planck-Institute, Mainz, Germany, 2006-2014
#
# Time-stamp: <2014-09-04 15:07:40 sander>
#
# substitute.awk substitutes each line in the input file that contains "regexp"
# with the contents of the file "newtext". Output is written to "outfile".
#
# example usage:
# gawk -f substitute.awk -v regexp="REGEXP" -v newtext=NEWTEXT -v outfile=gas2.eqn gas.eqn

# ----------------------------------------------------------------------------

BEGIN {
#system("rm -f " outfile"; touch "outfile) # create empty file
system("echo -n | sed \"s|-n||g\" > " outfile) # create empty file
status = -1
}

# ----------------------------------------------------------------------------

{
if (match($0, regexp) == 0) {
  print >> outfile
} else {
  if (status<0) {
    system("cat " newtext " >> " outfile)
    status = 0 # OK, the reaction number has been found
  }
}
}

# ----------------------------------------------------------------------------

END {
if (status<0) {
  print "ERROR: reaction number has not been found!"
  exit 1
}
}

# ----------------------------------------------------------------------------
