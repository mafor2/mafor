# ----------------------------------------------------------------------------
#
# Author:
#   Rolf Sander, Max-Planck-Institute, Mainz, Germany, 2010-...
#
# Time-stamp: <2013-08-19 11:32:42 sander>
#
# findrefs.awk is used by ncv. It provides the input for citefind which
# extracts all necessary references from literat.bib and copies them to
# meccalit.bib.
#
# ----------------------------------------------------------------------------

{
    # store BibTeX references like {&1400} in ref:
    # - either allow only numbers:
    if (match($0, "{&+([0-9]*)}", arr) != 0) {
    # - or allow anything:
    #   if (match($0, "{&+([^}]*)}", arr) != 0) {
      ref = arr[1]
      if (ref!="") {
        # transform comma-separated values to one number per line:
        gsub("[, ]+", "\n", ref)
        # print ref(s) to output file:
        printf "%s\n", ref #>> outfile
      }
    }
}

# ----------------------------------------------------------------------------
