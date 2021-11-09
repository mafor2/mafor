# ----------------------------------------------------------------------------

# Author: Rolf Sander, Max-Planck-Institute, Mainz, Germany, 2015-2017
# Time-stamp: <2018-12-28 10:31:55 sander>

# spc2mpl.awk transforms species file *.spc into *.py file for matplotlib

# example usage:
# gawk -f spc2mpl.awk -v mplfile=../pycaaba/_mecca_spc.py mecca.spc

# Normally, however, spc2mpl.awk is called by xmecca and not started
# directly from the command line.

# ----------------------------------------------------------------------------

BEGIN {
  print "running spc2mpl.awk..."
  # remove "../", then replace "." by "_" in input file
  basename = gensub("\\.\\./", "", "g", ARGV[1])
  gsub("\\.", "_", basename)
  # write header line:
  dontedit = "# created automatically by spc2mpl.awk, DO NOT EDIT!"
  print dontedit > mplfile
  print "def spc_names():" >> mplfile
  print "    return {"     >> mplfile
  errorstring = ""
}

# ----------------------------------------------------------------------------

{
  # store name for the species like {@H_2SO_4} in texname
  if (match($0, "{@([^}]*)}", arr) != 0) {
    texname = arr[1]
  } else {
    texname = ""
  }
  # convert long sub/superscripts like "_<10>" into "_{10}":
  oldname = texname
  texname = gensub("([_^])<([^>]+)>", "\\1{\\2}", "g", texname)
  # convert \FormatAq<SPECIES><##> to SPECIES(a##):
  texname = gensub("\\\\FormatAq<([^>]+)><([^>]+)>", "\\1(a\\2)", "g", texname)
  # if (texname != oldname) print oldname, texname
  # delete all comments {...} from $0:
  gsub("{[^}]*}", "")
  # is current line a line with a species definition?
  if (match($0, "^[ \t]*([A-Za-z][A-Za-z0-9_#]*)[ \t]*=.*;", arr) != 0) {
    kppname = arr[1]
    if (texname=="") {
      errorstring = sprintf("%s\nWARNING: species %s has no tex name",
        errorstring, kppname)
      texname = kppname
    }
    # define hfill for vertical alignment:
    hfill = substr("               ", 1, 15-length(kppname))
    printf "        '%s'%s: r'%s',\n", kppname, hfill, texname >> mplfile
  }
}

# ----------------------------------------------------------------------------

END {
  # dummy must be added because previous line ends with ",":
  print "        'Dummy'          : r'Dummy'}" >> mplfile
  if (errorstring!="") {
    print errorstring
  }
  printf "Input file:  %s\n", ARGV[1]
  printf "Output file: %s\n", mplfile
}

# ----------------------------------------------------------------------------
