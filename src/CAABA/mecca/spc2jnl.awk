# ----------------------------------------------------------------------------

# Author: Rolf Sander, Max-Planck-Institute, Mainz, Germany, 2015-2017
# Time-stamp: <2018-12-28 10:31:55 sander>

# spc2jnl.awk transforms species file *.spc into *.jnl file for ferret

# example usage:
# gawk -f spc2jnl.awk -v jnlfile=../jnl/_mecca_spc.jnl mecca.spc

# Normally, however, spc2jnl.awk is called by xmecca and not started
# directly from the command line.

# ----------------------------------------------------------------------------

BEGIN {
  print "running spc2jnl.awk..."
  # remove "../", then replace "." by "_" in input file
  basename = gensub("\\.\\./", "", "g", ARGV[1])
  gsub("\\.", "_", basename)
  # write header line:
  dontedit = "! created automatically by spc2jnl.awk, DO NOT EDIT!"
  print dontedit > jnlfile
  errorstring = ""
}

# ----------------------------------------------------------------------------

{
  # store ferret name for the species like {@H_2SO_4} in ferretname
  if (match($0, "{@([^}]*)}", arr) != 0) {
    ferretname = arr[1]
  } else {
    ferretname = ""
  }
  # convert long sub/superscripts like "_<abc>" into "_a_b_c":
  do {
    oldname = ferretname
    if (match(ferretname, "([_^])<([^>]+)>", arr) != 0) {
      text = gensub("(.)", arr[1] "\\1", "g", arr[2])
      sub("[_^]<[^>]+>", text, ferretname)
      #print oldname, ferretname
    }
  } while (ferretname != oldname)
  # aqueous-phase suffix:
  gsub("(a##)", "($ax)", ferretname)
  # delete all comments {...} from $0:
  gsub("{[^}]*}", "")
  # is current line a line with a species definition?
  if (match($0, "^[ \t]*([A-Za-z][A-Za-z0-9_#]*)[ \t]*=.*;", arr) != 0) {
    kppname = arr[1]
    gsub("_a##", "_($ax)", kppname)
    if (ferretname=="") {
      errorstring = sprintf("%s\nWARNING: species %s has no ferret name",
        errorstring, kppname)
      ferretname = kppname
    }
    if (kppname=="I") {
        kppname = "Iod" # "I" cannot be used in ferret as a variable name
    }
    printf "LET/D %-15s = 0 ; DEF SYMBOL %-15s = %s\n", 
        kppname, kppname, ferretname >> jnlfile
  }
}

# ----------------------------------------------------------------------------

END {
  printf "Input file:  %s\n", ARGV[1]
  printf "Output file: %s\n", jnlfile
  print errorstring
}

# ----------------------------------------------------------------------------
