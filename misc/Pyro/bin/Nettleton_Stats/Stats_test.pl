#!/bin/perl
#
# stats_test.pl - testing the Stats.pm module
#
use strict;

use lib "$ENV{HOME}/perllib";
use Stats qw(:DEFAULT factorial);


print "Factorial:\n";
print "factorial(0)\t", &factorial(0), "\n";
print "factorial(4)\t", &factorial(4), "\n";
print "factorial(100)\t", &factorial(100), "\n";
print "factorial(170)\t", &factorial(170), "\n";
print "factorial(171)\t", &factorial(171), "\n";


my @data;
for (1 .. 1000) {
  push @data, rand;
}

my $count = count(\@data);
my $min   = min(\@data);
my $max   = max(\@data);
my $mean  = mean(\@data);
my $var   = var(\@data);
my $std   = std(\@data);

print "\nDescriptive statistics for 1000 random numbers:\n",
  "count\t$count\n",
  "min\t$min\n",
  "max\t$max\n",
  "mean\t$mean\n",
  "var\t$var\n",
  "std\t$std\n\n";

my $N = 50000;
my ($lab, $dat) = &loadKCSdata;
print "NB test statistics for KCS data, $N permutations:\n";
my ($T, $p, $m, $sd) = NB_test($lab, $dat, \&hdist, $N);
print "mean\t$m\nstd\t$sd\nT\t$T\np-value\t$p\n\n";

sub loadKCSdata {

  my @lab = ();
  my @dat = ();

  open KCS, "KCS.txt" or die "Can't open KCS.txt!";
  while (<KCS>) {
    chomp;
    my ($label, @obs) = split / /;
    push @lab, $label;
    push @dat, [@obs];
    }
  close KCS;
  
  return (\@lab, \@dat);
}


print "\nStatistics for 2x2 contingency tables (see Bosch p. 127):\n\n";
my $h = [ [ 'B 9.1', 'Buch', '~Buch' ],
	  [    'TV',  19,        18  ],
	  [   '~TV',  24,        13  ]  ];

my ($chi2, $df, $chi2_P, $Cphi, $smallcount) = Chi2_test($h);
print_contingency_table($h);
print "\nChi square test:\nChi2\t$chi2\ndf\t$df\nP <= \t$chi2_P\nCramer's phi\t$Cphi\nsmall counts\t$smallcount\n";

my ($p1, $p2) = Fexact_test($h);
print "\nFisher's exact test:\np_one_sided\t$p1\np_two_sided\t$p2\n";

my ($Chi2_P, $Fe_P1, $Fe_P2) = test2x2($h);
print "\ntest2x2 (performs both tests if possible):\nChi2_P <= \t$Chi2_P\nFisher's exact 1-sided\t$Fe_P1\nFisher's exact 2-sided\t$Fe_P2\n";

print "\n";
