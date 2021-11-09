### ====================================================================
###  @Awk-file{
###     author          = "Nelson H. F. Beebe",
###     version         = "1.03",
###     date            = "19 February 1999",
###     time            = "12:22:50 MST",
###     filename        = "citefind.awk",
###     address         = "Center for Scientific Computing
###                        University of Utah
###                        Department of Mathematics, 322 INSCC
###                        155 S 1400 E RM 233
###                        Salt Lake City, UT 84112-0090
###                        USA",
###     telephone       = "+1 801 581 5254",
###     FAX             = "+1 801 585 1640, +1 801 581 4148",
###     URL             = "http://www.math.utah.edu/~beebe",
###     checksum        = "35241 163 650 5598",
###     email           = "beebe@math.utah.edu, beebe@acm.org,
###                        beebe@ieee.org (Internet)",
###     codetable       = "ISO/ASCII",
###     keywords        = "BibTeX, bibliography, citation label,
###                        citation tag",
###     supported       = "yes",
###     abstract        = "This program extracts bibliography entries
###                        from one or more bibliography files for
###                        specified citation tags ",
###     docstring       = "*********************************************
###                        This code is hereby placed in the PUBLIC
###                        DOMAIN and may be redistributed without any
###                        restrictions.
###                        *********************************************
###
###                        Read a list of cite tags from the file
###                        given by the first argument, then process
###                        each of the remaining file arguments as
###                        BibTeX bib files, looking up the
###                        bibentries, and output them to stdout.
###                        @preamble{} and @string{} entries are
###                        automatically output as well.
###
###                        Usage:
###                             nawk -f citefind.awk foo.tags bibfile(s) \
###                                 >newbibfile
###
###                        To be recognized, bib entries must look like
###
###                        @keyword{tag,
###                        ...
###                        }
###
###                        where the start @ appears in column 1, and
###                        the complete entry has balanced braces.
###
###                        This program should normally be run via
###                        the separate citefind shell script.
###
###                        The checksum field above contains a CRC-16
###                        checksum as the first value, followed by the
###                        equivalent of the standard UNIX wc (word
###                        count) utility output of lines, words, and
###                        characters.  This is produced by Robert
###                        Solovay's checksum utility.",
###  }
### ====================================================================

### Edit history (reverse chronological order):
### [19-Feb-1999]	1.03	Change BEGIN loop code to work properly.
###
###				Add END{} block with code to check for
###				missing and duplicate entries.
###
###				Add warning() function and three calls
###				to it.
###
###				Remove unneeded statement-terminating
###				semicolons.
###
###				Capitalize global variables.
###
### [30-Oct-1992]	1.02	Update header for new release.
###
### [21-Oct-1992]       1.01    Update for public distribution.
###
### [12-Apr-1989]       1.00    Original version.

BEGIN {
    while ((getline < ARGV[1]) > 0) # process first file only to get tags
	Tagnames[$0] = $0
    ARGV[1] = ""		# the main loop below must ignore this file
}

# @string and @preamble -- collect up to paired closing brace

/^@([Pp][Rr][Ee][Aa][Mm][Bb][Ll][Ee]|[sS][tT][rR][iI][nN][gG]){/        {
    printbraceditem()
    print ""
}

# "@keyword{tag," -- collect up to line starting with right brace
/^@[a-zA-Z0-9]*{/       {
	tag = substr($0,index($0,"{")+1)
	tag = substr(tag,1,index(tag,",")-1)
	if (tag in Tagnames)
	{
	    Count[tag]++
	    printbraceditem(tag)
	    print ""
	}
    }

END {
    # Do a final consistency check that each item was found, and found only once
    for (tag in Tagnames)
    {
	if (!(tag in Count))
	    warning("bibliography entry [" tag "] was not found")
	else if (Count[tag] > 1)
	    warning("bibliography entry [" tag "] was found " Count[tag] " times")
    }
    exit (Error_Count)
}

#=======================================================================

function bracecount(s, k,n)
{
    n = 0
    for (k = 1; k <= length(s); ++k)
    {
	if (substr(s,k,1) == "{")
	    n++
	else if (substr(s,k,1) == "}")
	    n--
    }
    return (n)
}

function printbraceditem(tag, count,last_filename)
{
    # Starting with the current contents of $0, print lines until we
    # reach a zero brace count.
    count = bracecount($0)
    print $0
    last_filename = FILENAME
    while (count != 0)
    {
	if ((last_filename != FILENAME) || (getline <= 0))
	    break
	print $0
	count += bracecount($0)
    }
    if (count != 0)
	warning("unexpected end-of-file in bibliography entry [" tag \
		"] from file [" last_filename "]")
}

function warning(message)
{
    print FILENAME ":" FNR ":%%" message >"/dev/stderr"
    Error_Count++
}

########################################################################
