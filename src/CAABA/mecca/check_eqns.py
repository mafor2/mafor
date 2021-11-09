#!/usr/bin/env python
# -*- coding: utf-8 -*- Time-stamp: <2018-06-14 15:08:41 sander>

# The python script check_eqns.py analyses the reactions present
# in a *.eqn file, assessing the conservation of the chemical elements
# as given in a *.spc file. Reactions found to be non-conserving are
# written to standard output along with a report of the elemental
# imbalance in each case.

# authors:
# original code in perl:            Tim Butler, Mainz (2008)
# rewritten from scratch in python: Rolf Sander, Mainz (2018)

# Usage:
#   - normally via xmecca
#   - manually with:
#     cat gas.spc aqueous.spc > tmp_.spc
#     cat gas.eqn aqueous.eqn > tmp_.eqn
#     ./check_eqns.py tmp_.spc tmp_.eqn

import sys, os
import re
caabadir = os.path.realpath(os.path.dirname(__file__)+'/..')
sys.path.append(os.path.realpath(caabadir+'/pycaaba'))
from pyteetime import tee # from pycaaba

HLINE =  '-' * 78
DEBUG = False #True

##############################################################################

def num_var(chunk):
    # split a chunk into number and variable, e.g.: 0.7H2O2 -> 0.7, H2O2
    search_result = re.search('([0-9.]*)(.*)', chunk)
    number = search_result.group(1)
    if (number):
        number = float(number)
    else:
        number = float(1)
    variable = search_result.group(2)
    return number, variable

##############################################################################

def get_elements_spc(spcfilename):

    elements = {} # create a dictionary
    print >> LOGFILE, HLINE
    SPCFILE = open(spcfilename, 'r')
    for line in iter(SPCFILE): # loop over *.spc file
        line = line.strip() # remove leading and trailing whitespace
        search_result = re.search('^ *([A-z0-9_#]+) *=([A-z0-9+ #]+);', line)
        if (not search_result): # skip lines that do not define a species
            if DEBUG: print 'NO:  |%s|' % (line)
            continue # proceed with next line in *.spc file
        if DEBUG: print HLINE
        species = search_result.group(1)
        elements[species] = {} # create a sub-dictionary for this species
        spc_composition = search_result.group(2).replace(' ','') # rm whitespace
        if DEBUG: print '|%s|' % (line)
        if DEBUG: print '|%s|' % (search_result.group())
        print >> LOGFILE, '%20s = ' % (species),
        if DEBUG: print '|%s|' % (spc_composition)
        for chunk in spc_composition.split('+'): # split into atoms
            count, element = num_var(chunk) # separate stoichiometric factor
            print >> LOGFILE, '%g * %s + ' % (count, element),
            elements[species][element] = count # add element count to dictionary
        print >> LOGFILE
    SPCFILE.close()
    print >> LOGFILE, HLINE
    if DEBUG: print elements
    return elements

##############################################################################

def analyze_eqn(elements, eqnfilename):

    exitstatus = 0
    EQNFILE = open(eqnfilename, 'r')
    regexp = re.compile('^<(.*)>(.*)=(.*):(.*);(.*)$')
    for line0 in iter(EQNFILE): # loop over *.eqn file
        massbal = {} # create a dictionary for mass balance
        # remove comments and leading and trailing whitespace:
        line = re.sub('{[^}]*}','',line0).strip()
        reaction = regexp.search(line) # check if line contains a reaction
        if reaction:
            eqntag = reaction.group(1).strip()        # e.g.: <G2111>
            reactants_str = reaction.group(2).strip() # e.g.: H2O + O1D
            products_str = reaction.group(3).strip()  # e.g.: 2 OH
            rate = reaction.group(4).strip()          # e.g.: 1.63E-10*EXP(60./temp)
            extra = reaction.group(5).strip()         # e.g.: // products assumed
            reactants = reactants_str.split('+')
            products  = products_str.split('+')
            if DEBUG: print '  eqntag:    <%s>' % (eqntag)
            # analyze products first:
            if DEBUG: print '  products:  %s' % (products_str)
            for chunk in products:
                chunk = chunk.replace(' ','') # rm whitespace
                stoic, product = num_var(chunk) # separate stoichiometric factor
                if DEBUG: print '|%s| |%s| |%s|' % (stoic, product, elements[product])
                # add all atoms from current product to mass balance:
                for key, value in elements[product].items():
                    if (key=='Min'): # subtract negative charge from Pls
                        key = 'Pls'
                        value = -value
                    if key in massbal:
                        massbal[key] += stoic * value
                    else:
                        massbal[key] = stoic * value
            # analyze reactants next:
            if DEBUG: print '  reactants: %s' % (reactants_str)
            for chunk in reactants:
                chunk = chunk.replace(' ','') # rm whitespace
                if (chunk=='hv'): continue # ignore pseudo-reactant hv
                stoic, reactant = num_var(chunk) # separate stoichiometric factor
                if DEBUG: print '|%s| |%s| |%s|' % (stoic, reactant, elements[reactant])
                # subtract all atoms from current reactant from mass balance:
                for key, value in elements[reactant].items():
                    if (key=='Min'): # subtract negative charge from Pls
                        key = 'Pls'
                        value = -value
                    if key in massbal:
                        massbal[key] -= stoic * value
                    else:
                        massbal[key] = -stoic * value
            # write errorstring if mass balance is not correct:
            errorstring = ''
            for element, balance in sorted(massbal.items()):
                # activate the following lines to ignore specific elements:
                if (element=='IGNORE'): continue
                if (element=='H'):      continue
                #if (element=='N'):      continue
                if (element=='O'):      continue
                if (abs(balance)>1E-14):
                    errorstring = '%s %+g %s,' % (errorstring, balance, element)
            if (errorstring):
                exitstatus = 4
                print '<%s> %s -> %s' % (eqntag, reactants_str, products_str)
                print '%s\n' % (errorstring.rstrip(','))
            if DEBUG: print '  rate:      %s' % (rate)
            if DEBUG: print '  extra:     %s' % (extra)

        else:
            if DEBUG: print 'NO:  %s' % (line)
    EQNFILE.close()
    return exitstatus

##############################################################################

if __name__ == '__main__':

    LOGFILE = tee.stdout_start('check_eqns.log', append=False) # stdout
    if len(sys.argv) > 2:
        spcfilename = sys.argv[1]
        eqnfilename = sys.argv[2]
    else:
        spcfilename = 'gas.spc'
        eqnfilename = 'gas.eqn'
    print 'spcfile = %s' % (spcfilename)
    print 'eqnfile = %s\n' % (eqnfilename)

    # get elemental composition from *.spc file:
    elements = get_elements_spc(spcfilename)

    # analyze mass balance in *.eqn file:
    exitstatus = analyze_eqn(elements, eqnfilename)
    if (exitstatus): sys.exit(exitstatus)
    tee.stdout_stop()

##############################################################################
