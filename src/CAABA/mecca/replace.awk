# ----------------------------------------------------------------------------
#
# Author:
#   Rolf Sander, Max-Planck-Institute, Mainz, Germany, 2006-...
#
# Time-stamp: <2018-01-30 14:55:28 sander>
#
# replace.awk replaces individual reactions in KPP equation files
#
# normally, replace.awk is called via xmecca
#
# example usage:
# gawk -f replace.awk -v infile="gas.eqn" -v outfile="gas2.eqn" example.rpl
#
# ----------------------------------------------------------------------------

BEGIN {
rplfile = "tmp_1.rpl"
eqnfile = "tmp_1.eqn"
system("cp " infile " " eqnfile)
}

# ----------------------------------------------------------------------------

{

if (match($0, "^#REPLACE[ 	]*<([A-Za-z0-9_]*[*]?)>", arr) != 0) {
  # #REPLACE command was found and the reaction number, i.e. the 
  # "([A-Za-z0-9]*[*]?)" part of the above regexp, is stored in arr[1].
  #  system("rm -f " rplfile"; touch "rplfile) # create empty file
  system("echo -n | sed \"s|-n||g\" > " rplfile) # create empty file
  getline
  replacementexists = 0 # 0 = false
  oriarr = arr[1]
  # change the wildcard "*" to a regexp:
  gsub("*", "[^>]*", arr[1])
  # loop until #ENDREPLACE is found:
  while (match($0, "^#ENDREPLACE") == 0) {
    replacementexists = 1 # 1 = true
    if ( arr[1] == oriarr ) {
      # add main reaction number into angle brackets:
      print gensub("<([A-Za-z0-9_]*)>", "<" arr[1] "\\1>", "g") >> rplfile
    } else {
      # leave reaction number unchanged for wildcards:
      print >> rplfile
    }
    getline
  }
  if ( arr[1] == "" ) {
    print "Adding new reaction(s) ..."
    system("cat " rplfile " >> " eqnfile)
  } else {
    if ( replacementexists == 1 ) {
      printf "Replacing reaction %s ...\n", oriarr
    } else {
      printf "Deleting reaction %s ...\n", oriarr
    }
    # insert the new equations into the eqn file:
    command = "gawk -f substitute.awk -v regexp=\"^[ 	]*<" arr[1] ">\" -v newtext=" rplfile " -v outfile=" outfile " " eqnfile
    if (system(command)>0) exit 1
    system("cp -f " outfile " " eqnfile)
  }
} else {
  # #REPLACE command was not found.
  # Empty lines and comments starting with "//" are okay,
  # otherwise print an error message:
  if ( (match($0, "^[ 	]*$") == 0) && (match($0, "^//") == 0) ) {
    printf "ERROR: %s\n", $0
  }
}

}

# ----------------------------------------------------------------------------

END {
system("cp -f " eqnfile " " outfile)
system("rm " rplfile " " eqnfile)
}

# ----------------------------------------------------------------------------
