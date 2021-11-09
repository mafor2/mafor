# ----------------------------------------------------------------------------
#
# Author: Rolf Sander, Max-Planck-Institute, Mainz, Germany, 2014
#
# Time-stamp: <2014-09-08 16:35:31 sander>
#
# Implementation of the cpp #include directive via awk
#
# example usage:
# gawk -f cpp_include.awk -v outfile=OUTPUTFILE INPUTFILE

# ----------------------------------------------------------------------------

BEGIN {
if (outfile == "") {outfile = "outfile"}
system("echo -n | sed \"s|-n||g\" > " outfile) # create empty file
status = 0
}

# ----------------------------------------------------------------------------

{
if (match($0, "^#include[ 	]+(.*)", arr) == 0) {
  print >> outfile
} else {
  # include a file
  # first, add prefix "include/" and suffix ".rpl" to filename if necessary:
  includefile = "include/" gensub("\\.rpl$", "", "g", arr[1]) ".rpl"
  if ((getline <includefile) == -1) {
    status = 1
    print "ERROR: cannot find include file " includefile
  } else {
    print "  Including " includefile
    print "// Start of rpl file " includefile >> outfile
    system("cat " includefile " >> " outfile)
    print "// End of rpl file " includefile >> outfile
  } 
  close(includefile)
}
}

# ----------------------------------------------------------------------------

END {
if (status != 0) {
  exit status
}
}

# ----------------------------------------------------------------------------
