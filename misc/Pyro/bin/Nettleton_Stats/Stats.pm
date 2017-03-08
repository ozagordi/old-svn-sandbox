#
# Stats.pm - Basic statistics functions
#

package Stats;

use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;

$VERSION = 1.00;
@ISA = qw(Exporter);

# autoexport:
@EXPORT = qw( $pi
	      log2 log10 pow
	      count sum min max mean var std median
	      random rand_perm resample
	      NB_test hdist
	      Fexact_test
	      Chi2_test
	      FDR_control
	      test2x2 print_contingency_table
	      entropy mut_inf
	    );
	     
# export on request:
@EXPORT_OK = qw(factorial factorial_approx n2perm);


my $pi = 3.141592653;

#
# General descriptive statistics
#

sub count {

  # count elements in array

  my $X = shift;

  return (scalar @$X);
}


sub sum {

  # sum of vector components

  my $X = shift;

  my $sum = 0;
  for (@$X) {
    $sum += $_;
  }

  return $sum;
}


sub max {

  # maximum value
  # use either max($x, $y) or max(\@x)

  if (scalar @_ == 2) {
    return max2(@_);
  }
  
  my $X = shift;

  my $n = count($X);
  if ($n < 1) {
    return;
  }
  elsif ($n == 1) {
    return $X->[0];
  }
  
  # else:
  my $i;
  my $max = $X->[0];
  for $i (1 .. $n-1) {
    $max = max2($max, $X->[$i]);
  }
  return $max;
}

sub max2 {

  my $a = shift;
  my $b = shift;

  $a > $b ? return $a : return $b;
}


sub min {

  # minimum value
  # use either min($x, $y) or min(\@x)

  if (scalar @_ == 2) {
    return min2(@_);
  }   
  
  my $X = shift;

  my $n = count($X);
  if ($n < 1) {
    return;
  }
  elsif ($n == 1) {
    return $X->[0];
  }
  
  # else:
  my $i;
  my $min = $X->[0];
  for $i (1 .. $n-1) {
    $min = min2($min, $X->[$i]);
  }
  return $min;
}

sub min2 {

  my $a = shift;
  my $b = shift;

  $a < $b ? return $a : return $b;
}


sub mean {

  # sample mean

  my $X = shift;

  return unless count($X);

  return sum($X) / count($X);
}


sub median {

  # sample median

  my $X = shift;

  my $n = count($X);
  
  return          if ($n <  1);
  return $X->[0]  if ($n == 1);
  
  my @Z = sort {$a <=> $b} @$X;
  if (is_odd($n)) {
    return ($Z[($n-1)/2] + $Z[($n+1)/2]) / 2;
  }
  else {
    return $Z[$n/2];
  }
}


sub is_odd {

  my $a = shift;

  return ($a % 2);
}


sub var {

  # (empirical) variance

  my $X = shift;

  my $n = count($X);
  return if ($n <= 0);
  return 0.0 if ($n == 1);

  my $mean = mean($X);
  my ($i, $s, $ss);

  for (@$X) {
    $s += $_;
    $ss += ($_ * $_);
  }

  return (1 / ($n-1)) * ($ss - (1/$n) * ($s * $s));
}
  
  
sub std {

  # standard deviation

  return sqrt(var(shift));
}

#
# END general descriptive statistics.
#



#
# mjd_permute  (Cookbook p. 126, Example 4-4)
#

# factorial
  
sub factorial {

  my $n = shift;
  
  my $fac = 1;

  for (1 .. $n) {
    $fac *= $_;
  }

  return $fac;
}


sub factorial_approx {

  my $n = shift;

  # approximation for n!
  # see http://mathworld.wolfram.com/StirlingsApproximation.html

  return sqrt((2*$n + (1/3))*$pi) * pow($n,$n) * exp(- $n);
}


sub n2pat {

  my $N   = shift;  # $N-th pattern ...
  my $len = shift;  # of length $len

  my @pat;
  my $i = 1;
  
  while ($i <= $len + 1) {
    push @pat, $N % $i;
    $N = int($N/$i);
    $i++;
  }

  return @pat;
}


sub pat2perm {  # turn pattern into permutation

  my @pat = @_;

  my @source = (0 .. $#pat);
  my @perm;

  push @perm, splice(@source, (pop @pat), 1) while @pat;

  return @perm;
}


sub n2perm {  # n2perm($N, $len) generates the $N-th permutation of $len objects 

  pat2perm(n2pat(@_));
}

#
# END mjd_permute
#



sub random {  # generate an integer random number ...
 
  my $x = shift;  # between $x ...
  my $y = shift;  # and $y (including both).

  return int(rand($y - $x + 1)) + $x;
}


sub rand_perm {

  # returns a random permutation of the $N integers 0, ..., $N-1

  my $N = shift;

  my @list = (0 .. $N-1);
  my @perm;

  while (@list) {
    my $r = random(0, $#list);
    push @perm, splice(@list, $r, 1);
  }

  return @perm;
}


sub resample {

  # resample B samples from {0, ..., N-1}

  my $N = shift;
  my $B = shift;

  my @S;

  for my $i (1 .. $B) {
    push @S, random(0, $N-1);
  } 

  return @S;
}
  

#
# False Discovery Rate (FDR) Controlling Screening of Multiple Tests
# under any dependency structure
#
# returns ($alpha_cutoff) where
#    - $alpha_cutoff is the largest alpha leading to rejection of H0,
#      all p-values larger than $alpha_cutoff lead to acceptance of H0
#      $alpha_cutoff == 0  <==>  no null hypotheses can be rejected
#
# Proceedure by Benjamini and Yekutieli (1997), Working Paper 97-3.
#

sub FDR_control {

  my $P   = shift;          # p-values
  my $FDR = shift || 0.01;  # false discovery rate

  my $m = count($P);
  my @P = sort {$a <=> $b} @$P;  # P[0] <= P[1] <= ... <= P[m-1]

  my (@h, $i);
  for $i (1 .. $m) {
    $h[$i-1] = 1/$i;
  }
  my $q = $FDR/sum(\@h);

  # find k = max{i | $P[i] <= (i/m)*q}
  $i = $m;
  while (($i >= 1) and ($P[$i-1] > ($i/$m) * $q)) {
    $i--;
  }
  
  if ($i == 0) {
    return (0, 0);
  }
  else {
    return ($P[$i-1], $i);
  }
  
  # reject H0_1, ..., H0_k with a FDR <= $FDR
  # i.e reject only those H0 with alpha <= $P[k] 
}

#
# END FDR_control
#
  


#
# Nettelton and Banerjee's proceedure for
#
#    "Testing the equality of distributions
#     of random vectors with categorical components"
#    Dan Nettelton, T. Banerjee
#    Comp. Stat. & Data Anal. 37 (2001) 195-208
#

sub NB_test {
  
  # permutation test
  #
  # - based on Dan Nettelton's S-PLUS implementation
  # - vector components are compared in a supplied function 'dist'
  #
  # returns  * test statistic T for the sample under consideration
  #          * its p-value as estimated by the permutation test
  #          * mean of the permutation test statistic
  #          * std of the permutation test statistic

  my $lab       = shift;  # (column) vector of lables
  my $dat       = shift;  # data matrix (rows = observations, 
                          #              cols = categorical components)
  my $dist_func = shift;  # reference to distance function
  my $N         = shift;  # number of permutations

  return if ($N <= 0);  # do not test

  # get distance matrix
  my $dist = getDistMat($dat, $dist_func);

  # get all nearest neighbors:
  my $NN = getNN($dat, $dist);

  # realized test statistic
  # (number of nearest neighbors from different groups):
  my $T = getT($lab, $NN);

  # permutation test of significance
  my $perm_T = permTest($lab, $NN, $N);
  my $count = 0;
  for (@$perm_T) {
    $count += ($_ <= $T);
  }

  return ($T, $count/$N, mean($perm_T), std($perm_T));
}


sub permTest {

  my $lab = shift;  # labels
  my $NN  = shift;  # nearest neighors
  my $N   = shift;  # number of samples to be drawn

  my $len = count($lab);  # number of elements to be permuted
  
  my ($i, @T);
  for $i (1 .. $N) {
    my @perm_lab = @$lab[rand_perm($len)];  # permute index vector
    push @T, getT(\@perm_lab, $NN);         # calculate T
  }

  return \@T;
}


sub getT {

  my $lab = shift;  # labels
  my $NN  = shift;  # nearest neighbors
  
  my $n = count($lab);
  my ($i, $j);

  my $T = 0;
  for $i (0 .. $n-1) {
    for $j ($i+1 .. $n-1) {
      # count edges (presence indicated in NN) from different groups:
      $T += ($lab->[$i] ne $lab->[$j]) * $NN->[$i][$j];
    }
  }

  return $T;
}


sub getDistMat {  # calculate distance matrix

  my $dat       = shift;  # data (rows = observations, cols = components)
  my $dist_func = shift;  # distance function 

  my $n = count($dat);  # number of observations
  my ($i, $j);
  my @dist = ();
  my $max_dist = 0;  # for diagonal

  # &$dist_func is assumed symmetric:
  for $i (0 .. $n-1) {
    for $j ($i+1 .. $n-1) {
      $dist[$i][$j] = $dist[$j][$i] = &$dist_func($dat->[$i], $dat->[$j]);
      $max_dist = max($max_dist, $dist[$i][$j]);
    }
  }

  # diagonal (this ensures i not to be i's nearest neighbor):
  for $i (0 .. $n-1) {
    $dist[$i][$i] = $max_dist;
  }
  
  return \@dist;
}


sub getNN {  # construct nearest neighbor subgraph

  my $dat  = shift;  # data (rows = observations, cols = components)
  my $dist = shift;  # distance matrix

  my $n = count($dat);  # number of observations
  my ($i, $j);

  # nearest neighbors (NNs):
  my @mindist = ();
  my @NN1 = ();

  for $i (0 .. $n-1) {
    $mindist[$i] = min(\@{$dist->[$i]}); # minimum distance from i to any j
    for $j (0 .. $n-1) {
      $NN1[$i][$j] = ($dist->[$i][$j] == $mindist[$i]);
      # NN1[i][j] == 1 iff j is a NN of i, else 0
    }
  }

  # realize the OR condition in the def. of D in the paper
  # by OR-merging the lower and upper triangle of NN1 to (the triangle) NN:
  my @NN = ();
  for $i (0 .. $n-1) {
    for $j ($i+1 .. $n-1) {
      $NN[$i][$j] = ($NN1[$i][$j] + $NN1[$j][$i] > 0); 
      # now we have NN[i][j] == 1 iff
      #    (j is a NN of i) OR (i is a NN of j), else 0
      # so NN[i][j] indicates if there is an edge between i and j
    }
  }
  
  return \@NN;
}



sub hdist {  # Hammning distance

  my $X = shift;
  my $Y = shift;

  my $n = count($X);
  my ($i);
  
  my $hdist = 0;
  for $i (0 .. $n-1) {
    $hdist += ($X->[$i] ne $Y->[$i]);
  }

  return $hdist;
}

#
# END NB permutation test.
#



#
# Fisher's exact test for 2x2 contingency tables
#

sub Fexact_P {

  # The probability of a given 2x2 table  A B  is:
  #                                       C D
  #     (A+B)!(C+D)!(A+C)!(B+D)!
  # P = ------------------------
  #        A!B!C!D!(A+B+C+D)!

  my ($A, $B, $C, $D) = @_;

  my $N = $A + $B + $C +$D;

  my @num;
  push @num, (1 .. $A + $B);
  push @num, (1 .. $C + $D);
  push @num, (1 .. $A + $C);
  push @num, (1 .. $B + $D);

  my @denom;
  push @denom, (1 .. $A);
  push @denom, (1 .. $B);
  push @denom, (1 .. $C);
  push @denom, (1 .. $D);
  push @denom, (1 .. $N);

  @num   = sort {$a <=> $b} @num;
  @denom = sort {$a <=> $b} @denom;

  my $P = 1.0;
  for my $i (0 .. $#num) {
    $P *= ($num[$i]/$denom[$i]);
  }
  
  return ($P > 0) ? $P : (1/$N);  # if P==0, be conservative
}


sub Fexact_P_old {

  # The probability of a given 2x2 table  A B  is:
  #                                       C D
  #     (A+B)!(C+D)!(A+C)!(B+D)!
  # P = ------------------------
  #        A!B!C!D!(A+B+C+D)!

  my ($A, $B, $C, $D) = @_;

  my $N = $A + $B + $C + $D;
  my $num = factorial($A + $B) * factorial($C + $D) *
    factorial($A + $C) * factorial($B + $D); 
  my $denom = factorial($A) * factorial($B) * 
    factorial($C) * factorial($D) * factorial($N);

  #printf "%g %g\n", $num / $denom, Fexact_P2($A,$B,$C,$D);

  return $num / $denom;
}


sub Fexact_test {

  # Calculate p-values (one and two sided) for Fisher's exact test:
  #
  # sum up all probabilities of configurations with the same row and
  # column sums that are less or equally likely than the given one.
  # Do so in the same direction (one sided) as the given configuration
  # or in both directions (two sided).

  my $h = shift;  # data table (first row and column contains labels)

  my $xa = $h->[1][1];  my $xb = $h->[1][2];
  my $ya = $h->[2][1];  my $yb = $h->[2][2];
  
  my $Pcutoff = Fexact_P($xa, $xb, $ya, $yb);

  # count main diagonal up
  my @right_tail = ( $Pcutoff );
  my $x_a = $xa;  my $x_b = $xb;
  my $y_a = $ya;  my $y_b = $yb;
  while ($x_b>0 and $y_a>0) {
    $x_a++;  $x_b--;
    $y_a--;  $y_b++;
    my $P = Fexact_P($x_a, $x_b, $y_a, $y_b);
    push @right_tail, $P if $P<=$Pcutoff;
  }
  
  # count main diagonal down
  my @left_tail = ( $Pcutoff );
  $x_a = $xa;  $x_b = $xb;
  $y_a = $ya;  $y_b = $yb;
  while ($x_a>0 and $y_b>0) {
    $x_a--;  $x_b++;
    $y_a++;  $y_b--;    
    my $P = Fexact_P($x_a, $x_b, $y_a, $y_b);
    push @left_tail, $P if $P<=$Pcutoff;
  }
  
  # one sided P:
  my $Ponesided = $xa*$yb > $xb*$ya ? sum(\@right_tail) : sum(\@left_tail);

  return ($Ponesided, sum(\@left_tail) + sum(\@right_tail) - $Pcutoff);
}

#
# END Fisher's exact test
#



#
# Pearson's Chi square test for  m x r  contingency tables  
#


sub Chi2_test {

  my $h = shift;  # table
                  # (first row $h[0][] and first column $h[][0] contain labels)
  my $m = shift;  # number of data rows    (optional)
  my $r = shift;  # number of data columns (optional)

  # Performs a chi-square test on the contingency table $h[1-$m][1-$r]
  #
  # see  Bosch, Elementare Einführung in die angewandte Statistik,
  #      Vieweg 1994, S. 128ff.

  ($m, $r) = complete_contingency_table($h, $m, $r);

  return chi_square($h, $m, $r);
}


sub chi_square {
  
  my $h = shift;  # table
                  # (first row $h[0][] and first column $h[][0] contain labels)
  my $m = shift;  # number of data rows  
  my $r = shift;  # number of data columns

  # Calculate Chi**2 test statistic

  my ($i, $k);
  my $chi2 = 0;
  my $sc = 0;  # # of cells with too low a count

  for $i (1 .. $m) {
    for $k (1 .. $r) {
      if ($h->[$i][$r+1] * $h->[$m+1][$k] == 0) {  # division by zero ahead?
	return;
      }
      $chi2 += (($h->[$i][$k] - (($h->[$i][$r+1] * $h->[$m+1][$k]) / $h->[$m+1][$r+1]))**2) / ($h->[$i][$r+1] * $h->[$m+1][$k]);
      if ($h->[$i][$k] < 5) {  # count too small entries
	$sc++;
      }
    }
  }
  $chi2 *= $h->[$m+1][$r+1];

  # degrees of freedom
  my $df = ($r - 1) * ($m - 1);

  # Cramer's phi
  my $Cphi = sqrt($chi2 / ($h->[$m+1][$r+1] * min($m-1, $r-1)));
  
  # p-value:
  my $chi2_P = 1 - chi_square_P($chi2) if ($df == 1);

  return ($chi2, $df, $chi2_P, $Cphi, $sc);
}


sub chi_square_P {

  my $chi2 = shift;

  # determine p-value for a given chi2 statistic (df = 1)

  # chi2 distribution with 1 degree of freedom:
  my %df1 = (  0.00   => 0.025,
	       0.004  => 0.05,
	       0.02   => 0.1,
	       
	       2.71   => 0.9,
	       3.84   => 0.95,
	       5.02   => 0.975,
	       6.63   => 0.99,
	       7.88   => 0.995,
	      10.83   => 0.999,
	      15.1367 => 0.9999,
	      19.5114 => 0.99999,
	      23.9281 => 0.999999,
	      28.3740 => 0.9999999,
	      32.8413 => 0.99999999,
	      37.3249 => 0.999999999,
	      41.8215 => 0.9999999999 );
  
  for (sort {$b <=> $a} keys %df1) {
    if ($chi2 >= $_) {
      return $df1{$_};
    }
  }

  return;  # error: $chi2 < 0
}
  

#
# END Chi square
#



#
# contingency tables
#

sub print_row {

  my $row = shift;
  my $r   = (scalar @{$row}) - 2;  # no. of data columns

  printf "%5s|", $row->[0];

  my $j;
  for $j (1 .. $r-1) {
    printf "%5s ", $row->[$j];
  }

  printf "%5s|%5s\n", $row->[$r], $row->[$r+1];
}


sub print_line {

  my $r = shift;  # no. of data colums

  print "-----+";
  for (1 .. $r-1) {
    print "-----+";
  }
  printf "-----+-----\n";
}


sub print_contingency_table {

  my $h = shift;  # table 

  my $m = (scalar @{$h}) - 2;       # no. of data rows
  my $r = (scalar @{$h->[0]}) - 2;  # no. of data columns
  my $i;

  print_row($h->[0]);
  print_line($r);
  for $i (1 .. $m) {
    print_row($h->[$i]);
  }
  print_line($r);
  print_row($h->[$m+1]);
}


sub complete_contingency_table {

  my $h = shift;  # table
                  # (first row $h[0][] and first column $h[][0] contain labels)
  my $m = shift;  # number of data rows    (optional)
  my $r = shift;  # number of data columns (optional)

  # Calculate row and column sums

  # guess $m and $r if not provided:
  $m = (scalar @{$h}) - 1      unless ($m);
  $r = (scalar @{$h->[0]}) - 1 unless ($r);

  my ($i, $j, $k);

  # row sums
  $h->[0][$r+1] = 'sum';
  for $i (1 .. $m) {
    $h->[$i][$r+1] = 0;
    for $j (1 .. $r) {
      $h->[$i][$r+1] += $h->[$i][$j];
    }
  }

  # column sums
  $h->[$m+1][0] = 'sum';
  for $j (1 .. $r) {
    $h->[$m+1][$j] = 0;
    for $i (1 .. $m) {
      $h->[$m+1][$j] += $h->[$i][$j];
    }
  }

  # total sum (sample size)
  $h->[$m+1][$r+1] = 0;
  for $i (1 .. $m) {
    $h->[$m+1][$r+1] += $h->[$i][$r+1];
  }

  return ($m, $r);
}

#
# END contingency tables
#



#
# 2x2 contingency tables
#

sub test2x2 {
  
  # analyze 2x2 contingency table
  # usage: test2x2($h, $verbose)  OR  test2x2($A, $B, $C, $D, $verbose)
  #    where $h is a reference to the table (first row and column omitted)
  #    or  |$A $B|
  #        |$C $D|  is the table.

  my ($h, $verbose);

  if (scalar @_ <= 2) {
    ($h, $verbose) = @_;
  }
  else {
    ($h->[1][1], $h->[1][2], $h->[2][1], $h->[2][2], $verbose) = @_;
  }
  
  # Chi2 test:
  my ($chi2, $df, $chi2_P, $Cphi, $smallcounts) = Chi2_test($h, 2, 2);
  if ($smallcounts) {
    $chi2_P = undef;
  }
  
  # Fishers exact test:
  my ($P1, $P2);
  #if ($h->[1][3] <= 170 and $h->[2][3] <= 170 and $h->[3][1] <= 170 and $h->[2][2] <= 170) {
    ($P1, $P2) = Fexact_test($h);
#  }

  if ($verbose) {
    print_contingency_table($h);
    if ($verbose >= 2) {
      printf "Chi square P  = %g\n", $chi2_P if ($chi2_P);
      printf "F exact 1-s P = %g\n", $P1     if ($P1);
      printf "F exact 2-s P = %g\n", $P2     if ($P2);
    }
  }
  
  return ($chi2_P, $P1, $P2);
}

#
# END 2x2 tables.
#



#
# Information Theory
#

sub log2 {
  
  my $x = shift;

  if (! $x or $x <= 0) {
    return;
  }
  
  return log($x) / log(2);
}


sub entropy {

  my $X = shift;
  my $Y = shift;  # optional

  # Calculate entropy of @$X 
  # or joint entropy of @$X and @$Y

  my $n = count($X);
  my %freq;
  for my $i (0 .. $n-1) {
    my $key = $Y ? "$X->[$i];$Y->[$i]" : $X->[$i];
    # ';' is a special character!
    $freq{$key}++;
  }
  
  my $H = 0;
  for (values %freq) {
    my $prob = $_ / $n;
    $H -= ($prob * log2($prob));
  }
  
  return $H;
}


sub mut_inf {

  my $X = shift;
  my $Y = shift;

  # Mutual information between @$X and @$Y

  # (this is nice, but not very fast...:)
  return entropy($X) + entropy($Y) - entropy($X, $Y);
}

#
# END Information theory
#



sub print_matrix {

  my $table  = shift;
  my $rowsep = shift;
  my $colsep = shift;

  $rowsep = defined $rowsep ? $rowsep : "\n";
  $colsep = defined $colsep ? $colsep : "\t";

  # prints 'NULL' for undefined values
  for my $i (0 .. $#{$table}) {
    for my $j (0 .. $#{$table->[$i]}) {
      my $v = defined $table->[$i][$j] ? $table->[$i][$j] : 'NULL';
      print $v, $colsep;
    }
    print $rowsep;
  }
}


sub print_right_upper_triangle {

  my $mat = shift;

  my $n = count($mat);
  my ($i, $j);

  for $i (0 .. $n-1) {
    for $j (0 .. $i) {
      print "   ";
    }
    for $j ($i+1 .. $n-1) {
      printf " %2d", $mat->[$i][$j];
    }
    print "\n";
  }
}


sub log10 {
  
  my $x = shift;

  if (! $x or $x <= 0) {
    return;
  }
  
  return log($x) / log(10);
}


sub pow {

  my $a = shift;
  my $b = shift;

  # calculate a to the power of b

  if ($a == 0 or $b == 0) {
    return 1;
  }
  else {
    return exp($b * log($a));
  }
}



1;



  
