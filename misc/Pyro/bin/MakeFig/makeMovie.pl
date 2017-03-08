#!/usr/bin/perl -w
#make a movie...
# but for now, just make one xfig file

#	string	orientation		("Landscape" or "Portrait")
#	string	justification		("Center" or "Flush Left")
#	string	units			("Metric" or "Inches")
#	string	papersize		("Letter", "Legal", "Ledger", "Tabloid",
#					 "A", "B", "C", "D", "E",
#					 "A4",   "A3", "A2", "A1", "A0" and "B5")
#	float	magnification		(export and print magnification, %)
#	string	multiple-page		("Single" or "Multiple" pages)
#	int	transparent color	(color number for transparent color for GIF
#					 export. -3=background, -2=None, -1=Default,
#					 0-31 for standard colors or 32- for user colors)
#	# optional comment		(An optional set of comments may be here,
#					 which are associated with the whole figure)
#	int	resolution coord_system	(Fig units/inch and coordinate system:
#					   1: origin at lower left corner (NOT USED)
#					   2: upper left)

$xfigHead = "#FIG 3.2 
Landscape
Center
Metric
A4
100.00
Single
-2
1200 2";

#$readColor = 2;
# colors for errors
$errorColor = 4;
# the colors for the population.  Should be of same length as @haplo
@haploColor = (8,12,24);

# @haplo contains the haplotypes in the population
#$haplo[0] = 'PQITLWQRPLVTIKIGGQVKEALLQAGADDTVLEEITLPGRWKPKMIGGYGGFIKVKQYDQIPLEICGHK';
#$haplo[1] = 'PQITLWQIPLVTIKIGGQVKEALLQAGADDTVLEEITLPGRWKPKMIGGIGGFIKVKQYDQIPLEICGHK';
#$haplo[2] = 'PQITLWQIPLVTIKIGGQVKEALGQAGADDTVLEEITLPGRWKPKMIGGYGGFIKVKQYDQIPLEICGHK';
#$haplo[0] = 'PQITLWQRPLVTIKIGGQVKEALLQAGADDTVLEEITLPGRWKPKMIGGYGGFIKV';
#$haplo[1] = 'PQITLWAIPLVTIKIGGQVKEALLQAGADDTVLEVITLPGRWKPKMIGVIGGFIKV';
#$haplo[2] = 'PQITLWAIPLVTIKIGGQVKEAGGQAGADDTVLEEITLPGRWKPKMIGGYGGFIKV';
$haplo[0] = 'GATATTGTTATTTATCAGTATATGGATGATTTGTATGTTGGTT';
$haplo[1] = 'GATATTGTTATTTATCAGTATTTTGATGATTTGTATGTTGGAC';



$n = length($haplo[0]);

#parameters for the picture
$errRate = .05;
$readLen = 10;
$numReads = 30;
# have to guess these to put the letters in right pos
$charWidth = 120;
$lineSep = 200;

#output file
$basefile = "test";
$figfile = "$basefile.fig";
$outtype = "pdf";
$outfile = "$basefile.$outtype";

$current = $xfigHead . "\n";
foreach my $j (0..@haplo-1) {
	foreach my $i (1..$n) {
		$current .= textLine(substr($haplo[$j],$i-1,1),$haploColor[$j],($i-1),$j);
	}
}

foreach my $rd (@haplo..$numReads+@haplo) {
	$whichHap = int(rand(@haplo));
	$pos = int(rand($n-3));
	$len = $readLen + int(rand(2)) - 2;
	#print "$pos ", substr($haplo,$pos,$len), "\n";
	$read = substr($haplo[$whichHap],$pos,$len);
	foreach my $j (0..length($read)-1) {
		if (rand() < $errRate) {
			$current .= textLine(drawAlpha(),$errorColor,$pos+$j,$rd);
		} else {
			$current .= textLine(substr($read,$j,1),$haploColor[$whichHap],$pos+$j,$rd);
		}
	}
}

open OUT, ">$figfile";
print OUT $current;
close OUT;
system "fig2dev -L $outtype -F $figfile > $outfile";
system "gv $outfile";

sub drawAlpha {
  #  my @aa = qw(A R N D C Q E G H I L K M F P S T W Y V);
  my @aa = qw(A C G T);
  return $aa[int(rand(scalar(@aa) -1))];
}

## a text line is: 
#   int	object 			(always 4)
#	int	sub_type		(0: Left justified
#					     1: Center justified
#					     2: Right justified)
#	int	color			(enumeration type)
#	int	depth			(enumeration type)
#	int	pen_style		(enumeration , not used)
#	int	font 			(enumeration type)
#	float	font_size 		(font size in points)
#	float	angle			(radians, the angle of the text)
#	int	font_flags		(bit vector)
#	float	height			(Fig units)
#	float	length			(Fig units)
#	int	x, y			(Fig units, coordinate of the origin
#					 of the string.  If sub_type = 0, it is
#					 the lower left corner of the string.
#					 If sub_type = 1, it is the lower
#					 center.  Otherwise it is the lower
#					 right corner of the string.)
#	char	string[]		(ASCII characters; starts after a blank
#					 character following the last number and
#					 ends before the sequence '\001'.  This
#					 sequence is not part of the string.
#					 Characters above octal 177 are
#					 represented by \xxx where xxx is the
#					 octal value.  This permits fig files to
#					 be edited with 7-bit editors and sent
#					 by e-mail without data loss.
#					 Note that the string may contain '\n'.)
# example:
#     *          *
# 4 0 4 50 -1 0 12 0.0000 4 135 525 2565 4455 GGTG\001
sub textLine {
	#input: text, color, xpos, ypos
	my $text = shift;
	my $color = shift;
	my $x = shift;
	my $y = shift;
	$x = $x * $charWidth;
	$y = $y * $lineSep;
	my $height = 120;
	my $length = 120 * length($text);
	my $size = 12;
	#my $font = 12; #courier
	my $font = 5; #latex tt
	my $fontFlag = 0; # 0 = latex fonts, 4=postscript fonts

	#rule for heitht and length:
	#with times new-roman:
	#size = 12 => height = 135 size = 14 => height = 150, 16 => 180

	#with courier:
	#12 => height = 105, length = 120 * length($text)

	my $ans = "4 0 $color 50 -1 $font $size 0.0000 $fontFlag $height $length $x $y $text\\001\n";
	return $ans;
}
