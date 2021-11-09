# ----------------------------------------------------------------------------
#
# Authors:
#   Rolf Sander, Max-Planck-Institute, Mainz, Germany, 2003-...
#   Astrid Kerkweg, Max-Planck-Institute, Mainz, Germany, 2005
#   Sergey Gromov, Max-Planck-Institute, Mainz, Germany, 2015
#
# Time-stamp: <2018-12-28 10:32:46 sander>
#
# eqn2tex.awk transforms eqn into tex.
#
# usage:
# gawk -f eqn2tex.awk ../mecca.eqn

# ----------------------------------------------------------------------------

BEGIN {
  print "running eqn2tex.awk..."
  # remove "../", then replace "." by "_" in input file:
  basename = gensub("\\.\\./", "", "g", ARGV[1])
  gsub("\\.", "_", basename)
  # define name of output files:
  hfile   = basename "_h.tex"
  afile   = basename "_a.tex"
  gfile   = basename "_g.tex"
  jfile   = basename "_j.tex"
  phfile  = basename "_ph.tex"
  hetfile = basename "_het.tex"
  eqfile  = basename "_eq.tex"
  iexfile = basename "_iex.tex" 
  hnotesfile   = basename "_h_notes.tex"
  anotesfile   = basename "_a_notes.tex"
  gnotesfile   = basename "_g_notes.tex"
  jnotesfile   = basename "_j_notes.tex"
  phnotesfile  = basename "_ph_notes.tex"
  hetnotesfile = basename "_het_notes.tex"
  eqnotesfile  = basename "_eq_notes.tex"
  iexnotesfile = basename "_iex_notes.tex" 
  logfile = "eqn2tex.log"
  # write header line:
  dontedit = "% created automatically by eqn2tex.awk, DO NOT EDIT!"
  print dontedit > hfile
  print dontedit > afile
  print dontedit > gfile
  print dontedit > jfile
  print dontedit > phfile
  print dontedit > hetfile
  print dontedit > eqfile
  print dontedit > iexfile
  print dontedit > hnotesfile
  print dontedit > anotesfile
  print dontedit > gnotesfile
  print dontedit > jnotesfile
  print dontedit > phnotesfile
  print dontedit > hetnotesfile
  print dontedit > eqnotesfile
  print dontedit > iexnotesfile
  print "log file" > logfile
  # initialize some strings:
  errorstring = ""
  unknown = "????"
  graycol = "\\rowcolor[gray]{0.95}"
}

# ----------------------------------------------------------------------------

{
  printf "\nline: %s\n", $0 >> logfile
  # store equation ID like G3102b in eqnid:
  if (match($0, "<([A-Za-z_0-9]+)>", arr) != 0) {
    eqnid = arr[1]
    printf "eqnid     = %s\n", eqnid >> logfile
    baseeqnid = gensub("([A-Z]+[0-9]+).*", "\\1", "g", eqnid)
    printf "baseeqnid = %s\n", baseeqnid >> logfile
  } else {
    eqnid = unknown
    printf "no eqnid found\n" >> logfile
  }
  # does current line contain a chemical equation, i.e. something like
  # " = ... : ... ; " ?
  if (match($0, "=.*:.*;") != 0) {
    printf "equation  = %s\n", $0 >> logfile
    # check if equation ID exists:
    if (eqnid==unknown) {
      errorstring = sprintf("%s\nERROR: This reaction has no eqnid:\n  %s",
        errorstring, $0)
    }
    # store marker like {%StTrGNJ} in marker:
    if (match($0, "{%([^}]*)}", arr) != 0) {
      marker = arr[1]
      printf "marker    = %s\n", marker >> logfile
    } else {
      marker = ""
      printf "no marker found\n" >> logfile
    }
    # store alternative LaTeX text for the rate
    # constant's A-factor like {@} in klatex:
    if (match($0, "{@([^}]*)}", arr) != 0) {
      klatex = arr[1]
      # replace <> by {}:
      gsub("<", "{", klatex)
      gsub(">", "}", klatex)
      printf "klatex    = %s\n", klatex >> logfile
    } else {
      klatex = ""
      printf "no klatex found\n" >> logfile
    }
    # store alternative LaTeX text for the rate constant's exponent
    # like {$} in kexplatex:
    if (match($0, "{\\$([^}]*)}", arr) != 0) {
      kexplatex = arr[1]
      # replace <> by {}:
      gsub("<", "{", kexplatex)
      gsub(">", "}", kexplatex)
      printf "kexplatex = %s\n", kexplatex >> logfile
    } else {
      kexplatex = ""
      printf "no kexplatex found\n" >> logfile
    }
    # store LaTeX note in note:
    if (match($0, "// *(.*)", arr) != 0) {
      note = arr[1]
      printf "note = %s\n", note >> logfile
    } else {
      note = ""
      printf "no note found\n" >> logfile
    }
    # store BibTeX references like {&1400} in ref:
    if (match($0, "{&([^}]+)}", arr) != 0) {
      ref = arr[1]
      if (ref=="SGN") {
        ref = "see general notes$^*$"
      } else {
        # put each BibTeX reference into \citet{}
        # here, "&" represents the matched regexp:
        gsub("[^, ]+", "\\citet{&}", ref)
      }
      # if there is a note, add an asterisk to refer to it:
      if (note!="") {
        ref = ref "$^*$"
      }
      printf "ref       = %s\n", ref >> logfile
    } else {
      if (match(eqnid, "^EQ[0-9]*b") != 0) {
        # equilibrium backward reaction EQ##b doesn't need a reference:
        printf "no ref needed for equilibrium backward reaction %s\n", eqnid >> logfile
      } else if (note!="") {
        ref = "see note$^*$"
        printf "no ref needed for reaction %s because there is a note\n", eqnid >> logfile
      } else {
        ref = "\\rule{8mm}{3mm}"
        errorstring = sprintf("%s\nERROR: reaction %s has no reference",
          errorstring, eqnid)
      }
    }
    # mz_sg_20150506+
    # isotopic composition transfer string
    if (match(marker, "isotrans") != 0) {
      # equnid is not shown, element transfer is passed via marker 
      # [ e.g. {%isotrans:C} ]
      split(marker, arr, ":")
      # remove "//" from $0 and equation:
      gsub("//", "", equation)
      gsub("//", "")
      eqnid = ""
      marker = arr[2]" transfer:"
      ref = ""
      note=""
    }
    # tagging reaction
    if (match(eqnid, "^TAG") != 0) {
      # eqnid is not shown, tagging info is passed via marker
      # [ e.g. {%tag:IC} ]
      # removing "TAG" to fork equation below
      gsub("TAG", "", baseeqnid)
      # adding whitespace after : in marker tag:CFG
      gsub("tag:", "tag: ", marker)
      eqnid = ""
      ref = ""
      note = ""
    }
    # mz_sg_20150506-
    printf "equation  = %s\n", $0 >> logfile
    # delete all comments {...} from $0:
    gsub("{[^}]*}", "")
    # delete all equation IDs <...> from $0:
    gsub("<([A-Za-z_0-9]+)>", "")
    printf "equation  = %s\n", $0 >> logfile
    # turn TABs into spaces:
    gsub("[\t]", " ")
    # reduce multiple spaces to one:
    gsub("  +", " ")
    printf "equation  = %s\n", $0 >> logfile
    # split into equation and rateconst using the kpp-syntax separators ":;":
    split($0, arr, " *[:;] *")
    equation  = arr[1]
    rateconst = "\\code{" arr[2] "}"
    printf "equation  = %s\n", equation >> logfile
    printf "rateconst = %s\n", rateconst >> logfile
    # delete reaction rate pseudo-species RR*:
    gsub("= RR[A-Za-z_0-9]+ \\+", "=", equation)
    printf "equation  = %s\n", equation >> logfile
    # put LaTeX command \kpp{} around each specie in equation
    # here, "&" represents the matched regexp:
    gsub("[A-Za-z][A-Za-z0-9_]*", "\\kpp{&}", equation)
    printf "equation  = %s\n", equation >> logfile
    # create reaction arrow
    gsub("=", "$\\rightarrow$", equation)
    printf "equation  = %s\n", equation >> logfile
    # use alternative LaTeX text for rate constant base?
    if (klatex!="") {
      rateconst = klatex
    }
    # write a line into the specific LaTeX table:
    if (match(baseeqnid, "^G[0-9]") != 0) {
      # G = gas-phase reaction:
      rowcol_g = ( rowcol_g=="" ? graycol : "" ) # alternate row color
      printf "%s\\code{%s} & %s & %s & %s & %s\\\\\n",
        rowcol_g, eqnid, marker, equation, rateconst, ref >> gfile
      # print note:
      if (note!="") {
        printf "\n\\note{%s}{%s}\n", eqnid, note >> gnotesfile
      }
      if (kexplatex!="") {
        errorstring = sprintf("%s\nERROR: reaction %s cannot use {$%s}",
          errorstring, eqnid, kexplatex)
      }
    } else if (match(baseeqnid, "^J[0-9]") != 0) {
      # J = photolysis reaction:
      rowcol_j = ( rowcol_j=="" ? graycol : "" ) # alternate row color
      printf "%s\\code{%s} & %s & %s & %s & %s\\\\\n",
        rowcol_j, eqnid, marker, equation, rateconst, ref >> jfile
      # print note:
      if (note!="") {
        printf "\n\\note{%s}{%s}\n", eqnid, note >> jnotesfile
      }
      if (kexplatex!="") {
        errorstring = sprintf("%s\nERROR: reaction %s cannot use {$%s}",
          errorstring, eqnid, kexplatex)
      }
    } else if (match(baseeqnid, "^PH[0-9]") != 0) {
      # PH = aqueous-phase photolysis reaction:
      rowcol_ph = ( rowcol_ph=="" ? graycol : "" ) # alternate row color
      printf "%s\\code{%s} & %s & %s & %s & %s\\\\\n",
        rowcol_ph, eqnid, marker, equation, rateconst, ref >> phfile
      if (note!="") {
        printf "\n\\note{%s}{%s}\n", eqnid, note >> phnotesfile
      }
      if (kexplatex!="") {
        errorstring = sprintf("%s\nERROR: reaction %s cannot use {$%s}",
          errorstring, eqnid, kexplatex)
      }
    } else if (match(baseeqnid, "^HET[0-9]") != 0) {
      # HET = HET reaction:
      rowcol_het = ( rowcol_het=="" ? graycol : "" ) # alternate row color
      printf "%s\\code{%s} & %s & %s & %s & %s\\\\\n",
        rowcol_het, eqnid, marker, equation, rateconst, ref >> hetfile
      if (note!="") {
        printf "\n\\note{%s}{%s}\n", eqnid, note >> hetnotesfile
      }
      if (kexplatex!="") {
        errorstring = sprintf("%s\nERROR: reaction %s cannot use {$%s}",
          errorstring, eqnid, kexplatex)
      }
    } else if (match(baseeqnid, "^A[0-9]") != 0) {
      # A = aqueous-phase reaction:
      rowcol_a = ( rowcol_a=="" ? graycol : "" ) # alternate row color
      printf "%s\\code{%s} & %s & %s & %s & %s & %s\\\\\n",
        rowcol_a, eqnid, marker, equation, rateconst, kexplatex, ref >> afile
      if (note!="") {
        printf "\n\\note{%s}{%s}\n", eqnid, note >> anotesfile
      }
    } else if (match(eqnid, "^EQ[0-9]*b") != 0) {
      # backward reaction, do nothing
    } else if (match(eqnid, "^EQ[0-9]*f") != 0) {
      # EQ = equilibria (acid-base and others):
      # delete f in equation id:
      gsub("f_", "_", eqnid)
      # replace rightarrow by rightleftharpoons:
      gsub("rightarrow", "rightleftharpoons", equation)
      rowcol_eq = ( rowcol_eq=="" ? graycol : "" ) # alternate row color
      printf "%s\\code{%s} & %s & %s & %s & %s & %s\\\\\n",
        rowcol_eq, eqnid, marker, equation, rateconst, kexplatex, ref >> eqfile
      if (note!="") {
        printf "\n\\note{%s}{%s}\n", eqnid, note >> eqnotesfile
      }
    #  mz_sg_20150505+
    } else if (match(baseeqnid, "^IEX.*") != 0) {
      # IEX = isotope exchange reaction
      # replace rightarrow by rightleftharpoons
      gsub("// ", "", equation)
      gsub("rightarrow", "rightleftharpoons", equation)
      equation = gensub("*", "^{abun}&", 1, equation)
      equation = gensub("*", "^{abun}&", 4, equation)
      equation = gensub("(.+)(.+)(.+)", "^{rare}&", 2, equation)
      equation = gensub("(.+)(.+)(.+)(.+)", "^{rare}&", 3, equation)
      rowcol_iex = ( rowcol_iex=="" ? graycol : "" ) # alternate row color
      printf "%s\\code{%s} & %s & %s & %s & %s\\\\\n",
        rowcol_iex, eqnid, marker, equation, rateconst, ref >> iexfile
      if (note!="") {
        printf "\n\\note{%s}{%s}\n", eqnid, note >> iexnotesfile
      }
    #  mz_sg_20150505-
    } else if (match(baseeqnid, "^H[0-9]") != 0) {
      # H = aqueous-phase heterogenous and Henry reaction:
      rowcol_h = ( rowcol_h=="" ? graycol : "" ) # alternate row color
      printf "%s\\code{%s} & %s & %s & %s & %s\\\\\n",
        rowcol_h, eqnid, marker, equation, rateconst,  ref >> hfile
      if (note!="") {
        printf "\n\\note{%s}{%s}\n", eqnid, note >> hnotesfile
      }
      if (kexplatex!="") {
        errorstring = sprintf("%s\nERROR: reaction %s cannot use {$%s}",
          errorstring, eqnid, kexplatex)
      }
    } else {
      errorstring = sprintf("%s\nERROR: unknown type of equation ID: %s",
        errorstring, eqnid)
      printf "\\code{%s} & %s & %s & %s & %s & %s\\\\\n", eqnid, marker, equation, rateconst, kexplatex, ref >> logfile
    }
  } else {
    # not an equation line
    # check for LaTeX text {@...}:
    if (match($0, "{@([^}]*)}", arr) != 0) {
      klatex = arr[1]
      printf "klatex    = %s\n", klatex >> logfile
      # write a line into the specific LaTeX table:
      if (match(baseeqnid, "^G[0-9]") != 0) {
        # G = gas-phase reaction:
        printf "%s\n", klatex >> gfile
      } else if (match(baseeqnid, "^J[0-9]") != 0) {
        # J = photolysis reaction:
        printf "%s\n", klatex >> jfile
      } else if (match(baseeqnid, "^PH[0-9]") != 0) {
        # PH = aqueous-phase photolysis reaction:
        printf "%s\n", klatex >> phfile
      } else if (match(baseeqnid, "^HET") != 0) {
        # HET = HET reaction:
        printf "%s\n", klatex >> hetfile
      } else if (match(baseeqnid, "^A[0-9]") != 0) {
        # A = aqueous-phase reaction:
        printf "%s\n", klatex >> afile
      } else if (match(baseeqnid, "^EQ[0-9]") != 0) {
        # EQ = equilibria (acid-base and others):
        printf "%s\n", klatex >> eqfile
      } else if (match(baseeqnid, "^H[0-9]") != 0) {
        # H = aqueous-phase heterogenous and Henry reaction:
        printf "%s\n", klatex >> hfile
      #  mz_sg_20150505+
      } else if (match(baseeqnid, "^IEX.*") != 0) {
        # IEX = isotope exchange reaction
        printf "%s\n", klatex >> iexfile
      #  mz_sg_20150505-
      } else {
        errorstring = sprintf("%s\nERROR: unknown type of equation ID: %s",
          errorstring, baseeqnid)
      }
    }
  }
}

# ----------------------------------------------------------------------------

END {
  printf "Input file:   %s\n", ARGV[1]
  printf "Output files: %s_*.tex\n", basename
  print errorstring
  print errorstring >> logfile
  print "(error messages about dummy reactions like D05 can be ignored)"
}

# ----------------------------------------------------------------------------
