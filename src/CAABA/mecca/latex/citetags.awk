### ====================================================================
###  @Awk-file{
###     author          = "Nelson H. F. Beebe",
###     version         = "1.03",
###     date            = "19 February 1999",
###     time            = "12:36:04 MST",
###     filename        = "citetags.awk",
###     address         = "Center for Scientific Computing
###                        University of Utah
###                        Department of Mathematics, 322 INSCC
###                        155 S 1400 E RM 233
###                        Salt Lake City, UT 84112-0090
###                        USA",
###     telephone       = "+1 801 581 5254",
###     FAX             = "+1 801 585 1640, +1 801 581 4148",
###     URL             = "http://www.math.utah.edu/~beebe",
###     checksum        = "31982 107 465 4936",
###     email           = "beebe@math.utah.edu, beebe@acm.org,
###                        beebe@ieee.org (Internet)",
###     codetable       = "ISO/ASCII",
###     keywords        = "BibTeX, bibliography",
###     supported       = "yes",
###     docstring       = "*********************************************
###                        This code is hereby placed in the PUBLIC
###                        DOMAIN and may be redistributed without any
###                        restrictions.
###                        *********************************************
###
###                        Read a LaTeX file, or an .aux file, and
###                        extract the citation tags, outputting them
###                        one per line on stdout.
###
###                        The companion program, citefind.awk, can be
###                        used to extract the complete BibTeX entries
###                        for those tags from a bibliography file.
###                        This is useful when a small bibtex file
###                        must be prepared from large citation
###                        databases.
###
###                        Usage:
###                             nawk -f citetags.awk foo.ltx | sort > temp.tags
###                             nawk -f citefind.awk temp.tags bibfile(s) \
###                                   >newbibfile
###                             rm temp.tags
###
###                        The sort step is optional; citetags will
###                        not output duplicate tags, but they will be
###                        in random order.
###
###                        This program should normally be run via
###                        the separate citetags shell script.
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
### [19-Feb-1999]	1.03	Capitalize global variables.
###
### [30-Oct-1992]	1.02	Update header for new release.
###
### [21-Oct-1992]       1.01    Update for public distribution.
###
### [12-Apr-1989]       1.00    Original version.

/(\\cit(e|ation)|\\bibcite)/ {
                    for (;;)    # loop over possible multiple \cite{} commands
                    {
                        $0 = substr($0,index($0,"\\cit"));
                        comment_start = index($0,"%");
                        if (comment_start > 0)
                            $0 = substr($0,1,comment_start-1);
                        line = substr($0,index($0,"{")+1);
                        for (;;)
                        {
                            k = index(line,"}");
                            if (k > 0)
                                break;
                            getline;
                            comment_start = index($0,"%");
                            if (comment_start > 0)
                                $0 = substr($0,1,comment_start-1);
                            line = line $0;
                        }
                        gsub(/[ \t]/,"",line);  # remove white space
                        templine = line;
                        line = substr(templine,1,k-1);
                        # We now have a comma-separated tag word list
                        split(line,taglist,",");
                        for (n in taglist)
                            Tagtable[m++] = taglist[n];
                        $0 = substr(templine,k);
                        if ($0 !~ /\\cite/)
                            break;              # no more \cite{} commands
                    }
                }

# Now print the (unsorted) list of unique tags.
END     {
    for (n in Tagtable)
        print Tagtable[n];
    }
