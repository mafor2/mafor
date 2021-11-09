#!/usr/bin/perl -w

# Analyse the reactions present in $ineqnfile, assessing the chemical
# elements as given in $spcfile.

# Author:
#   Rolf Sander, Nov 2015
#   based heavily on check_conservation.pl from Tim Butler

# EXAMPLE USAGE:
#   ./check_eqns.pl
#   ./check_eqns.pl myfile.spc myfile.eqn output_new.eqn output_free.eqn
# CHECK ALL FILES:
#   cd $M2 ; cat gas.eqn aqueous.eqn  rpl/*.rpl rpl/*/*.rpl > ~/tmp/all.eqn
#   cat gas.spc aqueous.spc > ~/tmp/all.spc
#   ./check_eqns.pl ~/tmp/all.spc ~/tmp/all.eqn
#   sort -k2 check_eqns_count.log> ~/tmp/check_eqns_count.log
#   e ~/tmp/check_eqns_count.log

$spcfile     = ($#ARGV>-1) ? "$ARGV[0]" : "gas.spc";
$ineqnfile   = ($#ARGV>0)  ? "$ARGV[1]" : "gas.eqn";
$neweqnfile  = ($#ARGV>1)  ? "$ARGV[2]" : "tmp_gas.eqn_new";
$freeeqnfile = ($#ARGV>2)  ? "$ARGV[3]" : "tmp_gas.eqn_free";
$countfile   = ($#ARGV>3)  ? "$ARGV[4]" : "check_eqns_count.log";
print "spcfile     = $spcfile\n";
print "ineqnfile   = $ineqnfile\n";
print "freeeqnfile = $freeeqnfile\n";
print "neweqnfile  = $neweqnfile\n\n";

# obtain elemental composition of all species:
$elements = &get_elements_spc($spcfile);

# open input and output files:
open (INEQNFILE,   "$ineqnfile")    or die $!;
open (FREEEQNFILE, ">$freeeqnfile") or die $!;
open (NEWEQNFILE,  ">$neweqnfile")  or die $!;
# open variables for printing into it:
open MASSBALANCE_ERROR,  '>', \$massbalance_error;
open CARBONCOUNT_ERROR,  '>', \$carboncount_error;
open ELEMENTLABEL_ERROR, '>', \$elementlabel_error;

# loop over all lines in input equation (*.eqn) file:
foreach $line (<INEQNFILE>) {
    ($number, $reactants, $products, $labels, $rate, $extra) = $line =~ m/
	^\s* <(.*?)>        # reaction number
	(.*?) = (.*?)       # reactants = products
	: \s+ {%(.*?)} \s+  # reaction labels
	(.*?);              # rate expression terminated by ";"
	(.*?)$              # everything else
	/x;
    # skip if current line is not a chemical equation:
    unless (defined $number && defined $reactants && defined $products &&
        defined $labels && defined $rate) {
        print NEWEQNFILE $line;  # print unchanged line
        print FREEEQNFILE $line; # print unchanged line
        next;
    }
    undef $extra;
    $reactants =~ s/{.*?}//g;
    $products =~ s/{.*?}//g;
    $reactants =~ s/^\s*(.*?)\s*$/$1/;
    $products =~ s/^\s*(.*?)\s*$/$1/;
    @reactants = split /\s*\+\s*/, $reactants;
    my @products = split /\s*\+\s*/, $products;

    # 1) MASS BALANCE
    my $imbalance = &stringify(&compare(&sum("reactant",@reactants), &sum("product",@products)));
    printf MASSBALANCE_ERROR "%-13s %-20s %s\n",
        "<$number>:", "$imbalance", $reactants
        if ($imbalance);

    # 2) CARBON COUNT
    $cmax=&cmax(@reactants);
    # 2nd digit of eqntag is 0 for >= 10 C atoms:
    my $cmaxeqntag = ($cmax<10) ? $cmax : 0;
    printf CARBONCOUNT_ERROR "%-10s C%-5s %s\n", "<$number>:", $cmax, $reactants
        if (((substr($number,0,2) eq "G4") || (substr($number,0,2) eq "J4"))
            && (substr($number,2,1) ne $cmaxeqntag));

    # 3) ELEMENT LABELS
    &comparelabels;
    printf ELEMENTLABEL_ERROR "%-10s %-10s -> %-10s %s\n",
        "<$number>:", "{%$labels}", "{%$labels_free}", $reactants
        unless ($labels eq $labels_free);
}

# 4) SPECIES COUNT
open (COUNTFILE, ">$countfile") or die $!;
foreach my $key (keys(%spccount)) {
    printf COUNTFILE "%-20s %4d %4d\n", $key, 
    $spccount{$key}{"reactant"}, $spccount{$key}{"product"};
}
close COUNTFILE;

close INEQNFILE;
close NEWEQNFILE;
close FREEEQNFILE;
close MASSBALANCE_ERROR;
close CARBONCOUNT_ERROR;
close ELEMENTLABEL_ERROR;

my $exitstatus = 0;
# print error messages:
if ($massbalance_error) {
    print "Incorrect mass balance:\n";
    print "$massbalance_error\n";
    $exitstatus += 4; # binary: 100
}
if ($carboncount_error) {
    print "Incorrect carbon count:\n";
    print "$carboncount_error\n";
    $exitstatus += 2; # binary: 010
}
if ($elementlabel_error) {
    print "Element labels are generated automatically,\n";
    print "Do not insert them manually!\n";
    print "$elementlabel_error\n";
    $exitstatus += 1; # binary: 001
}
exit $exitstatus;

sub comparelabels {
    $elementlist = &elementlist(@reactants);
    $labels_free = $labels;
    $labels_free =~ s/([A-Z])/ $1/g;
    $labels_free = "$labels_free "; # add trailing space
    $labels_free =~ s/ Br / /g;
    $labels_free =~ s/ Cl / /g;
    $labels_free =~ s/ C / /g;
    $labels_free =~ s/ F / /g;
    $labels_free =~ s/ H / /g;
    $labels_free =~ s/ Hg / /g;
    $labels_free =~ s/ I / /g;
    $labels_free =~ s/ N / /g;
    $labels_free =~ s/ O / /g;
    $labels_free =~ s/ S / /g;
    $labels_free =~ s/ //g;
    $line_free = $line;
    $line_free=~ s/{%$labels}/{%$labels_free}/g;
    print FREEEQNFILE $line_free; # output current line without element labels
    $line=~ s/{%$labels}/{%$labels_free$elementlist}/g;
    print NEWEQNFILE $line; # output current line with new element labels
}

sub compare {
    my ($lhs, $rhs) = @_;
    my %result;
    $result{$_} += $rhs->{$_} foreach keys %$rhs;
    $result{$_} -= $lhs->{$_} foreach keys %$lhs;
    # ignore some elements in mass balance:
    delete $result{Ignore} if defined $result{Ignore};
    delete $result{H}      if defined $result{H};
    delete $result{O}      if defined $result{O};
    # for charge balance, subtract anions from cations:
    if ((defined $result{Min}) && (defined $result{Pls})) {
        $result{Pls} -= $result{Min};
        delete $result{Min};
    }
    return \%result;
}

sub stringify {
    my $thing = $_[0];
    my $string = "";
    foreach my $element (keys %{ $thing }) {
        my $effect = $thing->{$element};
        next if abs($effect) < 1E-7;
        $effect = round($effect);
        $effect = "+$effect" unless $effect < 0;
        $string .= "$effect $element; ";
    }
    $string =~ s/\s$//; # delete trailing spaces
    return $string;
}

# input: array of strings like "0.7 HCHO" or "CH4"
sub sum {
    my $reacprod = shift; # get first argument ("reactant" or "product")
    my (@chunks) = @_;
    my %tally;
    foreach my $chunk (@chunks) { # loop over reactants or products
        my ($fraction, $species) = &_get_fraction($chunk);
        next if $species eq 'hv';
        die "Don't know anything about $species"
            unless exists $elements->{$species};
        $spccount{$species}{$reacprod} += 1;
        foreach my $element (keys %{ $elements->{$species} }) {
            $tally{$element} += $fraction * $elements->{$species}{$element};
        }
    }
    return \%tally;
}

# input: array of strings like "0.7 HCHO" or "CH4"
sub elementlist {
    my (@chunks) = @_;
    my $elementlist = "";
    my %tally;
    foreach my $chunk (@chunks) { # loop over reactants
        my ($fraction, $species) = &_get_fraction($chunk);
        next if $species eq 'hv';
        die "Don't know anything about $species"
            unless exists $elements->{$species};
        foreach my $element (keys %{ $elements->{$species} }) {
            $tally{$element} += $elements->{$species}{$element};
        }
    }
    foreach my $key (sort (keys(%tally))) {
        $elementlist .= "$key" unless
            (("$key" eq "C") && ($cmax == 1)) || # ignore C1
            ("$key" eq "H") || # ignore H
            ("$key" eq "Pls") || # ignore charge +
            ("$key" eq "Min") || # ignore charge -
            ("$key" eq "Ignore") || # ignore Ignore
            ("$key" eq "O");   # ignore O
    }
    return $elementlist;
}

# input: array of strings like "0.7 HCHO" or "CH4"
# output: max number of C atoms in any of the reactants
sub cmax {
    my (@chunks) = @_;
    my $tally = 0;
    foreach my $chunk (@chunks) { # loop over reactants
        my ($fraction, $species) = &_get_fraction($chunk);
        next if $species eq 'hv';
        die "Don't know anything about $species"
            unless exists $elements->{$species};
        if (exists $elements->{$species}{"C"}) {
            $tmp = $elements->{$species}{"C"};
            $tally = $tmp if $tally < $tmp;
            #print "$tally $elements->{$species}{'C'}\n";
        }
    }
    return $tally;
}

# input: *.spc filename
sub get_elements_spc {
    my $file = $_[0];
    open (SPCFILE, "$file") or die $!;
    my %elements;
    foreach my $line (<SPCFILE>) {
        #  mz_rs_20080628+
        # allow spaces at the start of a line in spc file:
        # my ($species, $elements) = $line =~ /^(\S+)\s*=\s*(.*)\s*;/;
        my ($species, $elements) = $line =~ /^\s*(\S+)\s*=\s*(.*)\s*;/;
        #  mz_rs_20080628-
        next unless defined $species and defined $elements;
        $spccount{$species}{"reactant"} = 0;
        $spccount{$species}{"product"} = 0;
        my @elements = split /\+/, $elements;
        foreach my $chunk (@elements) {
            my ($count, $element) = $chunk =~ /\s*(\d*)\s*([A-Za-z]+)\s*/;
            $count = 1 if $count eq "";
            #print "$species: @elements\n";
            $elements{$species}{$element} += $count;
        }
    }
    close SPCFILE;
    return \%elements;
}

# input is a string, e.g.: "0.7 HCHO" or "CH4"
# output is two strings, e.g.: "0.7", "HCHO" or "1", "CH4"
sub _get_fraction {
    my $string = shift;
    my ($fraction, $species);
    unless (($fraction, $species) = $string =~ /^([.\d]+)\s*(.*)/) {
        ($fraction, $species) = (1, $string);
    }
    return ($fraction, $species);
}

sub round {
    my $value  = $_[0];
    my $digits = (@_>1) ? "$_[1]" : 4; # default = 4
    open MYVAR, '>', \$variable;
    my $format = "%.${digits}g";
    printf MYVAR "$format", "$value";
    close MYVAR;
    #printf "Digits: %d OLD: %10g NEW: %10g\n", "$digits", "$value", "$variable";
    return $variable;
}
