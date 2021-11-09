# ----------------------------------------------------------------------------

# Author: Rolf Sander, Max-Planck-Institute, Mainz, Germany, 2006
# Time-stamp: <2006-02-11 21:00:51 sander>

# aerosol.awk extracts all aerosol species from mecca.spc and creates
# the declarations and initializations of the ind_*_a(:) arrays

# example usage (for 12 aerosol phases):
# gawk -f aerosol.awk -v apn=12 -v declfile=decl.f90 -v initfile=init.f90 mecca.spc

# Normally, however, aerosol.awk is called by xmecca and not started
# directly from the command line.

# ----------------------------------------------------------------------------

BEGIN {
# write header line
infostring = "  ! from aerosol.awk:"
print infostring > declfile
print infostring > initfile
}

# ----------------------------------------------------------------------------

{
  # is current line a definition of an aerosol species like "XYZ_a01 =" ?
  if (match($0, "^ *([A-Za-z0-9]+)_a01 *=", arr) != 0) {
    # define hfill for vertical alignment:
    hfill = substr("         ", 1, 9-length(arr[1]))
    printf "  INTEGER, PUBLIC, DIMENSION(APN) :: ind_%s_a %s= 0\n", 
      arr[1], hfill >> declfile
    for (i=1; i<=apn; i++) {
      printf "  ind_%s_a(%2.2d) %s= ind_%s_a%2.2d\n", 
        arr[1], i, hfill, arr[1], i >> initfile
    }
  }
}

# ----------------------------------------------------------------------------
