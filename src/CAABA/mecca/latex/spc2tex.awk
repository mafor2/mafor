# ----------------------------------------------------------------------------

# Author: Rolf Sander, Max-Planck-Institute, Mainz, Germany, 2003-2005
# Time-stamp: <2018-12-28 10:32:46 sander>

# spc2tex.awk transforms species file *.spc into tex file.

# usage:
# gawk -f spc2tex.awk ../mecca.spc

# ----------------------------------------------------------------------------

BEGIN {
  print "running spc2tex.awk..."
  # remove "../", then replace "." by "_" in input file
  basename = gensub("\\.\\./", "", "g", ARGV[1])
  gsub("\\.", "_", basename)
  # define name of output file
  texfile     = basename ".tex"
  # write header line
  dontedit = "% created automatically by spc2tex.awk, DO NOT EDIT!"
  print dontedit > texfile
  print "\\makeatletter" >> texfile
  print "\\def\\defkpp#1#2{\\expandafter\\def\\csname #1\\endcsname{\\chem{#2}}}%" >> texfile
  print "\\def\\kpp#1{\\@ifundefined{#1}{\\errmessage{#1 undefined}}{\\csname #1\\endcsname}}%" >> texfile
  # KPP dummy species hv and PROD:
  print "\\defkpp{hv}{h\\nu}%" >> texfile
  print "\\defkpp{PROD}{products}%" >> texfile
  errorstring = ""
}

# ----------------------------------------------------------------------------

{
  # store LaTeX text for the species like {@H_2SO_4} in latexname
  if (match($0, "{@([^}]*)}", arr) != 0) {
    latexname = arr[1]
  } else {
    latexname = ""
  }
  # replace <> by {}
  gsub("<", "{", latexname)
  gsub(">", "}", latexname)
  # delete all comments {...} from $0
  gsub("{[^}]*}", "")
  # is current line a line with a species definition?
  if (match($0, "^[ \t]*([A-Za-z][A-Za-z0-9_]*)[ \t]*=.*;", arr) != 0) {
    kppname = arr[1]
    if (latexname=="") {
      errorstring = sprintf("%s\nWARNING: species %s has no LaTeX name",
        errorstring, kppname)
      latexname = kppname
    }
    # write a line into the LaTeX table
    printf "\\defkpp{%s}{%s}%%\n", kppname, latexname >> texfile
  }
}

# ----------------------------------------------------------------------------

END {
  print "\\makeatother" >> texfile
  printf "Input file:  %s\n", ARGV[1]
  printf "Output file: %s\n", texfile
  print errorstring
}

# ----------------------------------------------------------------------------
