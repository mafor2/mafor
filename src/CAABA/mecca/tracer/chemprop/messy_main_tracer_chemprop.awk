# ----------------------------------------------------------------------------

# Author: Rolf Sander, Max-Planck-Institute, Mainz, Germany, 2016
# Time-stamp: <2017-02-14 18:05:38 sander>

# messy_main_tracer_chemprop.awk converts chemical properties listed in
# messy_main_tracer_chemprop.tbl into the f90 include file
# messy_main_tracer_chemprop.inc.

# Usage:
# setenv LC_ALL C ; gawk -f messy_main_tracer_chemprop.awk messy_main_tracer_chemprop.tbl

# Normally, messy_main_tracer_chemprop.awk is called via xchemprop and
# not started directly from the command line.

# ----------------------------------------------------------------------------

BEGIN {
  TRUE = 1
  FALSE = 0
  # define name of output file:
  incfile = gensub("\\.tbl", ".inc", "g", ARGV[1])
  # define name of LaTeX definitions file:
  texdefsfile = "tmp_chemprop_defs.tex"
  # write header line
  dontedit = "Created automatically by messy_main_tracer_chemprop.awk, DO NOT EDIT!"
  INT_PUB_PAR = "INTEGER, PUBLIC, PARAMETER :: "
  printf "! -*- f90 -*- %s\n\n", dontedit > incfile
  errorstring = ""
  spcnum = 0
  nprop = 0
  print "------------------------------------"
  print "col chemprop             default"
  print "------------------------------------"
}

# ----------------------------------------------------------------------------

function trim(string)
{
  # remove leading spaces:
  gsub("^ +", "", string)
  # remove trailing spaces:
  gsub(" +$", "", string)
  return string
}

# ----------------------------------------------------------------------------

function add_string(string1, string2)
{
  # Add string2 to string1, using a comma as separator. If necessary,
  # also add a line break.
  separator = ""
  if (length(string1) > 0) {
    separator = ", "
    if (length(gensub(".*\n", "", "g", string1))+length(string2) > 70) {
      separator = ", &\n  "
    }
  }
  return string1 separator string2
}

# ----------------------------------------------------------------------------

{

  # skip empty lines:
  if (match($0, "^ *$", arr) != 0) {
    next
  }

  # skip comment line starting with "#":
  if (match($0, "^#", arr) != 0) {
    next
  }

  # skip line separator:
  if (match($0, "^\\|-", arr) != 0) {
    next
  }

  # LaTeX definitions:
  if (match($0, "^\\\\", arr) != 0) {
    print $0 > texdefsfile
    next
  }

  # define properties:
  if (match($0, "^\\| *<CHEMPROP>", arr) != 0) {
    integernum = 0
    realnum = 0
    stringnum = 0
    NAMES_CASK_I_CHEMPROP = ""
    NAMES_CASK_R_CHEMPROP = ""
    NAMES_CASK_S_CHEMPROP = ""
    # split current line (field separator is "|" ) and save in "chemprop":
    nprop = split($0, chemprop, "|")
    for (i = 3; i < nprop; i++) {
      chemprop[i] = trim(chemprop[i])
      if (match(chemprop[i], "^I_") != 0) {
        integernum += 1
        printf INT_PUB_PAR "%s = %d\n",
          chemprop[i], integernum > incfile
        NAMES_CASK_I_CHEMPROP = add_string(NAMES_CASK_I_CHEMPROP, 
          sprintf("'%-15s'", gensub("^I_", "", "g", chemprop[i])))
      }
      if (match(chemprop[i], "^R_") != 0) {
        realnum += 1
        printf INT_PUB_PAR "%s = %d\n",
          chemprop[i], realnum > incfile
        NAMES_CASK_R_CHEMPROP = add_string(NAMES_CASK_R_CHEMPROP, 
          sprintf("'%-15s'", gensub("^R_", "", "g", chemprop[i])))
      }
      if (match(chemprop[i], "^S_") != 0) {
        stringnum += 1
        printf INT_PUB_PAR "%s = %d\n",
          chemprop[i], stringnum > incfile
        NAMES_CASK_S_CHEMPROP = add_string(NAMES_CASK_S_CHEMPROP, 
          sprintf("'%-15s'", gensub("^S_", "", "g", chemprop[i])))
      }
    }
    # header:
    printf INT_PUB_PAR "MAX_CASK_I_CHEMPROP = %d\n", integernum > incfile
    printf INT_PUB_PAR "MAX_CASK_R_CHEMPROP = %d\n", realnum    > incfile
    printf INT_PUB_PAR "MAX_CASK_S_CHEMPROP = %d\n", stringnum  > incfile
    printf "CHARACTER(LEN=NAMES_CASK_STRLEN), DIMENSION(MAX_CASK_I_CHEMPROP), &\n" > incfile
    printf "  PARAMETER, PUBLIC :: NAMES_CASK_I_CHEMPROP = (/ &\n" > incfile
    printf "  %s /)\n", NAMES_CASK_I_CHEMPROP > incfile
    printf "CHARACTER(LEN=NAMES_CASK_STRLEN), DIMENSION(MAX_CASK_R_CHEMPROP), &\n" > incfile
    printf "  PARAMETER, PUBLIC :: NAMES_CASK_R_CHEMPROP = (/ &\n" > incfile
    printf "  %s /)\n", NAMES_CASK_R_CHEMPROP > incfile
    printf "CHARACTER(LEN=NAMES_CASK_STRLEN), DIMENSION(MAX_CASK_S_CHEMPROP), &\n" > incfile
    printf "  PARAMETER, PUBLIC :: NAMES_CASK_S_CHEMPROP = (/ &\n" > incfile
    printf "  %s /)\n", NAMES_CASK_S_CHEMPROP > incfile
    printf "TYPE, PUBLIC :: t_meta_chemprop\n" > incfile
    printf "  CHARACTER(LEN=20) :: kppname\n" > incfile
    printf "  INTEGER, DIMENSION(MAX_CASK_I_CHEMPROP) :: cask_i\n" > incfile
    printf "  REAL(DP), DIMENSION(MAX_CASK_R_CHEMPROP) :: cask_r\n" > incfile
    printf "  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(MAX_CASK_S_CHEMPROP) :: cask_s\n" > incfile
    printf "END TYPE t_meta_chemprop\n\n" > incfile
    next
  }

  if (nprop == 0) {
    print "ERROR: properties must be defined first!"
    exit 1
  }

  # define LaTeX info:
  if (match($0, "^\\| *<LATEX>", arr) != 0) {
    # split current line (field separator is "|" ) and save in "latex":
    split($0, latex, "|")
    texfilenum = 0
    printlatex = FALSE
    for (i = 3; i < nprop; i++) {
      latex[i] = trim(latex[i])
      latextoken[i] = ""
      if (match(latex[i], "^<")) {
        gsub("^< *", "", latex[i])
        texfile[++texfilenum] = sprintf("tmp_chemprop_%s.tex", chemprop[i])
        latexstring = ""
        printlatex = TRUE
        latextoken[i] = "<"
        columns = "l"
      }
      if (match(latex[i], ">$")) {
        gsub(" *>$", "", latex[i])
        if (latextoken[i] == "<") {
          latextoken[i] = "="
        } else {
          latextoken[i] = ">"
        }
      }
      if (printlatex) {
        columns = columns "ll"
        latexstring = sprintf("%s & %s & Reference ", latexstring, latex[i])
      }
      if (latextoken[i]==">" || latextoken[i]=="=") {
        header = sprintf("KPP name%s\\\\\n\\hline\n", latexstring)
        printf "\\begin{longtable}{%s}\n", columns > texfile[texfilenum]
        printf "\\hline\n%s\\endfirsthead\n", header > texfile[texfilenum]
        printf "\\hline\n%s\\endhead\n", header > texfile[texfilenum]
        printf "\\hline\n\\endfoot\n" > texfile[texfilenum]
        printlatex = FALSE
      }
    }
    ntexfilenum = texfilenum
    next
  }

  # define defaultvalue:
  if (match($0, "^\\| *<DEFAULT>", arr) != 0) {
    # split current line (field separator is "|" ) and save in
    # "defaultvalue":
    n = split($0, defaultvalue, "|")
    for (i = 3; i < n; i++) {
      defaultvalue[i] = trim(defaultvalue[i])
      printf "%3d %-20.20s %s\n", i, chemprop[i], defaultvalue[i]
    }
    next
  }

  # data entry:
  spcnum += 1
  integerdata = ""
  realdata = ""
  stringdata = ""
  # split current line (field separator is "|" ) and save in array items:
  n = split($0, items, "|")
  kppname = trim(items[2])
  if (n != nprop) {
    print "ERROR: ", $0
    exit 1
  }
  texfilenum = 0
  printlatex = FALSE
  # loop over data columns:
  for (i = 3; i < n; i++) {
    item = trim(items[i])
    # if empty, use default:
    if (item == "") {
      item = defaultvalue[i]
      if (item == "") {
        printf "ERROR: %s is missing for %s (there is no default).\n", 
          chemprop[i], kppname
        exit 1
      }
    }
    if (latextoken[i]=="<" || latextoken[i]=="=") {
      printlatex = TRUE
      printf "%s ", kppname > texfile[++texfilenum]
    }
    # store BibTeX references like {&1400} in ref:
    if (match(item, "{&([^}]+)}", arr) != 0) {
      ref = arr[1]
      # delete bibliography info from item:
      gsub(" *{&[^}]*} *", "", item)
      if (ref=="SN") {
        ref = "see notes"
      } else {
        # put each BibTeX reference into \citet{}
        # here, "&" represents the matched regexp:
        gsub("[^, ]+", "\\citet{&}", ref)
      }
    } else {
      ref = "???"
    }
    if (printlatex) {
      latexitem = gensub("_", "\\\\_", "g", item)
      printf "& %s & %s ", latexitem, ref > texfile[texfilenum]
    }
    if (latextoken[i]==">" || latextoken[i]=="=") {
      printf "\\\\\n" > texfile[texfilenum]
      printlatex = FALSE
    }
    # INTEGER:
    if (match(chemprop[i], "^I_") != 0) {
      integerdata = add_string(integerdata, item)
    }
    # REAL:
    if (match(chemprop[i], "^R_") != 0) {
      if (chemprop[i] == "R_molarmass") {
        # convert sum formula into molar mass:
        gsub("[A-Z][a-z]?", "+M&*", item)
        gsub("^+", "", item)
        gsub("\\*\\+", "+", item)
        gsub("*$", "", item)
      }
      # add "_DP" to all REAL numbers:
      realdata = add_string(realdata,
        gensub("\\<([0-9]+\\.?[0-9]*([eE][-+]?[0-9]+)?)", "\\1_DP", "g", item))
    }
    # STRING:
    if (match(chemprop[i], "^S_") != 0) {
      # the number 24 here must be identical to STRLEN_MEDIUM:
      stringdata = add_string(stringdata, sprintf("'%-24s'", item))
    }
  }
  printf "TYPE(t_meta_chemprop), PARAMETER :: spc_%03d = t_meta_chemprop( &\n",
    spcnum > incfile
  printf "  '%s', &\n", kppname > incfile
  printf "  (/ %s /), & ! integer\n", integerdata > incfile
  printf "  (/ %s /), & ! real\n",    realdata    > incfile
  printf "  (/ %s /) ) ! string\n\n",   stringdata  > incfile
}

# ----------------------------------------------------------------------------

END {
  print "------------------------------------"
  # Put parts of LaTeX file together:
  printf "%% %s\n", dontedit > "chemprop.tex"
  system("cat chemprop_begin.tex >> chemprop.tex")
  system("cat " texdefsfile " >> chemprop.tex")
  for (i = 1; i <= ntexfilenum; i++) {
    printf "\\end{longtable}\n\\newpage" > texfile[i]
    system("cat " texfile[i] " >> chemprop.tex")
  }
  system("cat chemprop_end.tex >> chemprop.tex")
  # Create f90 variable "chemprop" in incfile:
  spclist = ""
  printf "! combine all species into one array:\n" > incfile
  printf INT_PUB_PAR "N_CHEMPROP = %d\n", spcnum > incfile
  print "TYPE(t_meta_chemprop), PUBLIC, PARAMETER, DIMENSION(N_CHEMPROP) :: chemprop = &" > incfile
  for (i = 1; i <= spcnum; i++) {
    spclist = add_string(spclist, sprintf("spc_%03d", i))
  }
  printf "  (/ %s /)\n", spclist > incfile
}

# ----------------------------------------------------------------------------
