
#
# Nettleton and Banerjee's proceedure for
#
#    "Testing the equality of distributions
#     of random vectors with categorical components"
#    Dan Nettleton, T. Banerjee
#    Comp. Stat. & Data Anal. 37 (2001) 195-208
#

sub NB_test {

  # permutation test
  #
  # - based on Dan Nettleton's S-PLUS implementation
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

