# ----------------------------------------------------------------------------

# Author: Rolf Sander, Max-Planck-Institute, Mainz, Germany, 2003-2018
# Time-stamp: <2018-12-28 10:31:55 sander>

# Extract ind_XYZ from messy_mecca_kpp_Parameters.f90 and insert these
# species into messy_mecca_idt_si.inc. For all ind_XYZ that are not zero
# it also inserts the species into messy_mecca_c2mr_si.inc,
# messy_mecca_mr2c_si.inc, and messy_mecca_trac_si.inc, using the tracer
# info in tracdef.tbl.

# 1) messy_mecca_idt_si.inc:  the index declaration "INTEGER idt_XYZ"
# 2) messy_mecca_c2mr_si.inc: the conversion concentration -> mixing ratio
# 3) messy_mecca_mr2c_si.inc: the conversion mixing ratio -> concentration
# 4) messy_mecca_trac_si.inc: the "CALL new_tracer" command including all tracer
#                             properties as defined in tracdef.tbl.

# usage:
# gawk -f tracdef_tbl.awk -v submodel=mecca -v tracdef=tracdef.tbl messy_mecca_kpp_Parameters.f90

# Normally, however, tracdef_tbl.awk is called by xmecca and not started
# directly from the command line.

# ----------------------------------------------------------------------------

BEGIN {
# define names of output files
# (normally, submodel="mecca", however, it can also be "mtchem")
idtfile     = "messy_" submodel "_idt_si.inc"
conc2mrfile = "messy_" submodel "_c2mr_si.inc"
mr2concfile = "messy_" submodel "_mr2c_si.inc"
tracfile    = "messy_" submodel "_trac_si.inc"
# write header line:
dontedit = "! -*- f90 -*- created automatically by xmecca, DO NOT EDIT!"
print dontedit > idtfile
print dontedit > conc2mrfile
print dontedit > mr2concfile
print dontedit > tracfile
print "Executing tracdef_tbl.awk..."
print "List of tracers:"
errorstring = ""
}

# ----------------------------------------------------------------------------

# if prop[i] contains more than just spaces, add an optional parameter
# to "CALL new_tracer" in tracfile. Then, increment i so that the next
# property is check in the folowing call of addproperty

function addproperty(property)
 {
  if (match(prop[i], "^ *$") == 0) 
    printf "  %s = %s, &\n", property, prop[i] >> tracfile
  i++
  }

function makesettracer(property)
 {
     if (match(prop[i], "^ *$") == 0) { 
      printf "CALL set_tracer(status, setname, idt_%s, %s, %s)\n",
       arr[1], property, prop[i] >> tracfile
      printf "CALL tracer_halt(substr, status)\n" >> tracfile
	  }
  i++
  }

# ----------------------------------------------------------------------------

{
# is current line something like "INTEGER ... ind_XYZ = nnn" with any nnn ?
if (match($0, "INTEGER.*ind_([A-Za-z0-9_]+) = ", arr) != 0) {
  # ind_XYZ was found and the species name XYZ, i.e. the "([A-Za-z0-9]+)" part
  # of the above regexp, is stored in arr[1].
  # define hfill for vertical alignment
  hfill = substr("               ", 1, 15-length(arr[1]))
  # 1) add tracer to messy_mecca_idt_si.inc
  printf "INTEGER :: idt_%s%s= 0\n", arr[1], hfill >> idtfile
}
# is current line something like "INTEGER ... ind_XYZ = nnn" with nnn != 0 ?
if (match($0, "INTEGER.*ind_([A-Za-z0-9_]+) = [^0]", arr) != 0) {
  # split arr (field separator is "_" ) and save basename
  # and subname in array basesub
  split(arr[1], basesub, "_")
  # Now try to find this species in the tracdef file. The species name
  # must be the first non-whitespace word in the line.

  if (basesub[2]=="") {
    subhfill = ""
    command = ("grep -E \"^\\| *" basesub[1] " \" " tracdef)
  } else {
    subhfill = substr("            ", 1, 1+length(basesub[2]))
    command = ("grep -E \"^\\| *" basesub[1] " +\\| +" basesub[2] " \" " tracdef)
  }
  if ( ((command | getline) > 0) && ( (arr[1]!="H2O") || (tag=="y") ) ) {
    # 'XYZ' was found and stored in $0 by getline
    printf "%s ", arr[1]
    # ------------------------------------------------------------------------
    # 2) add tracer to messy_mecca_c2mr_si.inc
    printf "zmrac(:,idt_%s)%s= riac(:) * Conc(:,ind_%s)\n",
      arr[1], hfill, arr[1] >> conc2mrfile
    # ------------------------------------------------------------------------
    # 3) add tracer to messy_mecca_mr2c_si.inc
    printf "Conc(:,ind_%s)%s= c_air(:) * zmrbc(:,idt_%s)\n",
      arr[1], hfill, arr[1] >> mr2concfile
    # ------------------------------------------------------------------------
    # 4) add tracer to messy_mecca_trac_si.inc
    # start of CALL new_tracer
    printf "%s%s%s%s%s%s\n", "CALL new_tracer(status, setname, '", 
      basesub[1], "',", hfill, subhfill, "modstr, &" >> tracfile
    printf "  quantity             = AMOUNTFRACTION, &\n" >> tracfile
    # mz_sg_20140407+
    #[per s] units for embudget passive loss/prod tracers XPTL*/XPTP*
    if ( match(basesub[1],"XPT[LP].+") == 1 ) {
      printf "  unit                 = 'mol/mol/s', &\n"  >> tracfile
    } else {
      printf "  unit                 = 'mol/mol', &\n"    >> tracfile
    }
    # mz_sg_20140407-
    # optional parameters for new_tracer if defined in the tracdef file.
    # split $0 (field separator is | surrounded by spaces)
    # save individual tracer properties in elements of array prop
    split($0, prop, " *\\| *")
    if (match(prop[3], "^ *$") == 0) 
      printf "  subname              = '%s', &\n", prop[3] >> tracfile
    # mz_sg_170515+
    # reference for tagged species/tracers
    if (match(prop[4], "^ *$") == 0) 
      printf "  refspec              = '%s', &\n", prop[4] >> tracfile
    # mz_sg_170515-    
    # prop[2] = tracer name, prop[3] = subname, prop[4] = refspec, other
    # properties start at i=5
    i=5
    # the order of the addproperty calls must be the same as the order of
    # the corresponding entries in the tracdef file.
    addproperty("medium              ")
    # last part of CALL new_tracer
    printf "  idx                  = idt_%s)\n", arr[1] >> tracfile
    printf "CALL tracer_halt(substr,status)\n" >> tracfile

    # the order of the makesettracer calls must be the same as the order of
    # the corresponding entries in the tracdef file.
    makesettracer("R_vini              ")
    makesettracer("I_force_init        ")
    makesettracer("I_advect            ")
    makesettracer("I_convect           ")
    makesettracer("I_vdiff             ")
    makesettracer("I_wetdep            ")
    makesettracer("I_drydep            ")
    makesettracer("I_sedi              ")
    makesettracer("I_scav              ")
    makesettracer("I_mix               ")
    makesettracer("I_force_col         ")
    makesettracer("I_integrate         ")
    makesettracer("I_timefilter        ")
    makesettracer("I_aerosol_method    ")
    makesettracer("S_aerosol_model     ")
    makesettracer("I_aerosol_mode      ")
    makesettracer("I_aerosol_nsol      ")

    # ------------------------------------------------------------------------
  } else {
    # print current species in parenthesis
    printf "(%s) ", arr[1]
    # add another line to errorstring
    if (arr[1]=="H2O") {
      # Missing H2O is okay. This tracer is not provided by MECCA:
      errorstring = sprintf("%s\nINFO: MECCA does not produce a H2O tracer.", 
        errorstring)
    } else {
      # Other species should not be missing in process_*.tbl:
      errorstring = sprintf("%s\nWARNING: %s%s not found in process_gas.tbl/process_aqueous.tbl", 
        errorstring, arr[1], hfill)
    }
  }
  # --------------------------------------------------------------------------
  # close grep command to avoid too many open processes
  close(command)
}
}

# ----------------------------------------------------------------------------

END {
print errorstring
printf "(warnings about process_gas.tbl/process_aqueous.tbl are only\n"
printf "for ECHAM, not for the box model)\n\n"
}

# ----------------------------------------------------------------------------
