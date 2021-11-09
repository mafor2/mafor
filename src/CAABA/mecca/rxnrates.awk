# ----------------------------------------------------------------------------
#
# Author:
#   Rolf Sander, Max-Planck-Institute, Mainz, Germany, 2009
#
# Time-stamp: <2018-12-28 10:31:55 sander>
#
# rxnrates.awk creates ferret jnl files to plot reaction rates
#
# usage: invoke via xmecca
#
# for debugging, activate lines starting with "#DEBUG"
#
# ----------------------------------------------------------------------------

BEGIN {
  printf "working on %s...\n", ARGV[1]
  logfile  = "rxnrates.log"
  dontedit = "created automatically by rxnrates.awk, DO NOT EDIT!"
  printf "%s\n",   dontedit > logfile
  printf "! %s\n", dontedit > jnlfile1
  printf "! %s\n", dontedit > jnlfile2
  # initialize errorstring:
  errorstring = ""
}

# ----------------------------------------------------------------------------

# analyze one side of the equation:
function analyze(prodloss, oneside, eqnid, equation){
  # split oneside into individual terms:
  # (note "\\" before "+" sign!)
  printf "total %s     = |%s|\n", prodloss, oneside >> logfile
  n = split(oneside, arr3, " *\\+ *")
  for (i = 1; i <= n; i++) {
      term = arr3[i]
      #DEBUG printf "%s %d is |%s|\n", prodloss, i, term >> logfile
      # add space between factor and species:
      sub ("^[0-9.]+", "& ", term)
      #DEBUG printf "%s %d is |%s|  ", prodloss, i, term >> logfile
      # check if there is already a factor:
      if (match(term, "^[A-Za-z]") == 0) {
        #DEBUG printf "(factor yes)\n" >> logfile
      } else {
        #DEBUG printf "(factor no)\n" >> logfile
        # add factor "1":
      sub ("^", "1 ", term)
      }
      #DEBUG printf "%s %d is |%s|\n", prodloss, i, term >> logfile
      # split term into factor and species:
      split(term, arr2, " +")
      factor  = arr2[1]
      species = arr2[2]
      printf "%s factor  %d = |%s|\n", prodloss, i, factor  >> logfile
      printf "%s species %d = |%s|\n", prodloss, i, species >> logfile
      if (prodloss=="loss") {
        sign = "-"
      } else {
        sign = "+"
      }
      # print budget unless species is a "RR*" dummy variable:
      if (match(species, "^RR") == 0) {
        printf "LET/UNITS=\"mol/mol/s\" rate = (%s%s) * RR%s ; ", 
          sign, factor, eqnid >> jnlfile2
        printf "GO _plot_rxnrates_scaled rate \"%s*%s: %s\" \"%s\"\n", 
          factor, eqnid, loss, species >> jnlfile2
      }
  }
}

# ----------------------------------------------------------------------------

{
  printf "\nline: |%s|\n", $0 >> logfile
  # delete all comments {...} from $0:
  gsub("^//.*}", "")
  # store equation ID like G3102b or A8703_a## in eqnid:
  if (match($0, "<([A-Za-z_0-9#]+)>", arr) != 0) {
    eqnid = arr[1]
    printf "eqnid     = |%s|\n", eqnid >> logfile
  } else {
    eqnid = unknown
    printf "no eqnid found\n" >> logfile
  }
  # does current line contain a chemical equation, i.e. something like
  # " = ... : ... ; " ?
  if (match($0, "=.*:.*;") != 0) {
    printf "equation  = |%s|\n", $0 >> logfile
    # check if equation ID exists
    if (eqnid==unknown) {
      errorstring = sprintf("%s\nERROR: This reaction has no eqnid:\n  %s",
        errorstring, $0)
    }
    # delete all comments {...} from $0:
    gsub("{[^}]*}", "")
    printf "equation  = |%s|\n", $0 >> logfile
    # delete eqnid from $0:
    gsub("<" eqnid ">", "")
    printf "equation  = |%s|\n", $0 >> logfile
    # delete rate constant from $0:
    gsub(":.*", "")
    printf "equation  = |%s|\n", $0 >> logfile
    # reduce multiple spaces to one:
    gsub("  +", " ")
    printf "equation  = |%s|\n", $0 >> logfile
    # remove leading spaces:
    gsub("^ +", "")
    printf "equation  = |%s|\n", $0 >> logfile
    # remove trailing spaces:
    gsub(" +$", "")
    printf "equation  = |%s|\n", $0 >> logfile
    # split into loss and prod:
    split($0, arr, " *= *")
    loss = arr[1]
    prod = arr[2]
    analyze("loss", loss, eqnid, $0)
    analyze("prod", prod, eqnid, $0)
    printf "go _plot_rxnrates RR%s \"%s: %s\"\n", eqnid, eqnid, loss >> jnlfile1
  }
}

# ----------------------------------------------------------------------------

END {
  if (errorstring!="") {
    print errorstring
  }
}

# ----------------------------------------------------------------------------
