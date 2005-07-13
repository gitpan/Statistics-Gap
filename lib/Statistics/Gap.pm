package Statistics::Gap;

use 5.008005;
use strict;
use warnings;

require Exporter;
use AutoLoader qw(AUTOLOAD);

our @ISA = qw(Exporter);

our @EXPORT = qw( gap );

our $VERSION = '0.02';

# pre-requisites
use GD;
use GD::Text;
use GD::Graph::lines;
use GD::Graph::colour;

# global variable
my $distfuncref;
my @d = ();

# Calculates optimal number of clusters that should be used to 
# cluster the input data using the "Gap Statistics"
sub gap
{
    # input params
    my $prefix = shift;
    my $matrixfile = shift;
    my $distmeasure  = shift;
    my $clustmtd = shift;
    my $K = shift;
    my $B = shift;
    my $i = 0;
    my $j = 0;

    # To optimize the code we have pulled out this if-else loop
    # from the innermost for loop below. Depending upon the
    # distance measure requested we set the function pointer
    # here and thus the if-else loop is executed only once.
    if($distmeasure=~/euclidean/i)    # Euclidean measure 
    {
	$distfuncref = \&dist_euclidean;
    } 
    elsif($distmeasure=~/squared/i)    # Squared Euclidean measure 
    {
	$distfuncref = \&dist_euclidean_sqr;
    }
    elsif($distmeasure=~/manhattan/i)    # Manhattan measure 
    {
	$distfuncref = \&dist_manhattan;
    }
    
    # read the matrix file into a 2 dimensional array.
    my @inpmat = ();
    open(INP,"<$matrixfile") || die "Error opening input matrix file!";

    # extract the number of rows from the first line in the file.
    my $rcnt = 0;
    my $ccnt = 0;
    my $line;

    $line = <INP>;
    chomp($line);
    $line=~s/\s+/ /;
    
    ($rcnt,$ccnt) = split(/\s+/,$line);

    # copy the complete matrix to a 2D array
    while(<INP>)
    {
	# remove the newline at the end of the input line
	chomp;

	# skip empty lines
	if(m/^\s*\s*\s*$/)
	{
	    next;
	}

	# remove leading white spaces
	s/^\s+//;

	# seperate individual values in a line
	my @tmp = ();
	@tmp = split(/\s+/);
	
	# populate them into the 2D matrix
	push @inpmat, [ @tmp ]; 
    }

    close INP;

    my @row1 = ();
    my @row2 = ();

    # calculate the pairwise distances for all possible unique pairs (n(n-1)/2)
    for($i = 0; $i < $rcnt; $i++)
    {
	# for all the rows in the cluster
	for($j = $i+1; $j < $rcnt; $j++)
	{
	    @row1 = @{$inpmat[$i]};
	    
	    @row2 = @{$inpmat[$j]};
	    
	    $d[$i][$j] = &$distfuncref(\@row1, \@row2);
	    
#	    print "distance between $i and $j = $d[$i][$j]\n";
	}
    }
    
    #~~~~~~~~~~~~~~~~~~~~ Step 1: Calculate the Error Measure (Wk) for the observed data ~~~~~~~~~~~~~~~~~~~~~~~~

    #local variables
    my $k = 0;
    my @W = ();
    my @k = ();
    my @tmp_W =();

    # loop through K times
    for($k=1; $k<=$K; $k++)
    {
	my $lineNo = 0;
	my %hash = ();

	# cluster the input matrix(nXm) i.e. cluster n observations into k clusters
	my $out_filename = "tmp.op" . $k . time();
	my $status = 0;
	$status = system("vcluster --clmethod $clustmtd $matrixfile $k >& $out_filename ");
	die "Error running vclusters \n" unless $status==0;
	
	# read the clustering output file
	open(CO,"<$matrixfile.clustering.$k") || die "Error opening clustering output file.";

	my $clust = 0;
	while($clust = <CO>)
	{
	    # hash on the cluster# and append the observation# 
	    chomp($clust);
	    if(exists $hash{$clust})
	    {
		$hash{$clust} .= " $lineNo";
	    }
	    else
	    {
		$hash{$clust} = $lineNo;
	    }

	    # increment the line number
	    $lineNo++;
	}

	close CO;

	# Calculate the "Within Cluster Dispersion Measure / Error Measure" Wk 
	# for given matrix and k value.
	$W[$k] = &error_measure(\@inpmat, \%hash);

	# for the graph plotting
	$tmp_W[$k-1] = $W[$k];
	$k[$k-1] = $k;

	unlink "$out_filename","$matrixfile.clustering.$k";	
    }

    # For plotting of graph of log(Wk) vs. K
    my $my_graph = new GD::Graph::lines(600,480);
    
    $my_graph->set_title_font(gdGiantFont,24);
    $my_graph->set_x_label_font(gdGiantFont,14);
    $my_graph->set_y_label_font(gdGiantFont,14);
    
    $my_graph->set('x_label' => 'number of clusters k',
		   'y_label' => 'log(Wk)',
		   'title' => 'log(Wk) by k',
		   'transparent' => '0',
		   ) or warn $my_graph->error;

    open(FILE,">$prefix.fig2.png");
    print FILE $my_graph->plot([\@k,\@tmp_W])->png;
    close(FILE);

    #~~~~~~~~~~~~~~ End of Step 1 ~~~~~~~~~~~~~~~~~~~~~~~


    #~~~~~~~~~~~~~~~~~~~~ Step 2: Generation of Reference Model  ~~~~~~~~~~~~~~~~~~~~~~~~

    my @refmat = ();
    my @W_M = ();

    my @min = ();
    my @diff = ();

    my $p = 0;

    # Calculate the min and max for each column in the observed
    # matrix and store it in respective arrays

    # for each column (i.e. feature) in the observed data matrix
    for($p = 0; $p < $ccnt; $p++)
    {
	my @inpcol = ();
	
	# extract the column from the observed matrix
	for ($j = 0; $j < $rcnt; $j++) 
	{
	    $inpcol[$j] = $inpmat[$j][$p];
	}
	
	# find the range of the column values (i.e. min amd max)
	my @sorted = ();
	my $max = 0 ;

	# sort numerically ascending
	@sorted = sort {$a <=> $b} @inpcol;

	# store the required values.
	$min[$p] = $sorted[0];
	$max = $sorted[$#sorted];
	$diff[$p] = $max - $min[$p];
	
    }

    # Repeat the complete Reference Distribution generation procedure B times.
    for($i = 1; $i <= $B; $i++)
    {

	# create the reference matrix
	for($p = 0; $p < $ccnt; $p++)
	{
	    # generate n (n=# of rows i.e. observations in the observed matrix) 
	    # random variables over the range of $min and $max and store them as
	    # a column in the reference matrix
	    my $rand;
	    for($j = 0; $j < $rcnt; $j++)
	    {
		# scale the generated random number from the range of (0,1) to 
		# the range of (min,max) found from the observed data.
		$rand = rand();	
		$rand *= $diff[$p];
		$rand += $min[$p];
		$refmat[$j][$p] = $rand;
	    }
	}

	# Calculate Wk* from the generated reference matrix which consists of:
	# 1. Cluster the generated n observations
	# 2. Calculate the Wk*
	
	# Write the matrix to a temporary file
	my $filename = "tmp.ref." . time();
	open(RO,">$filename") || die "Error opening temporary file ($filename) in write mode.\n";
	
	my $ccnt = $#{$refmat[0]} + 1;

	print RO scalar(@refmat) ."	" . $ccnt  . "\n";

	for $j (@refmat)
	{
	    print RO "@$j\n";
	}

	close RO;

	# loop through K times
	for($k=1 ; $k<=$K ; $k++)
	{
	    # cluster the input matrix(nXm) i.e. cluster n observations into k clusters
	    my $out_filename = "tmp.ref.op." . $k . "." . time();
	    my $status = 0;
	    $status = system("vcluster --clmethod $clustmtd $filename $k >& $out_filename");
	    die "Error running vclusters \n" unless $status==0;

	    # read the clustering output file
	    open(CO,"<$filename.clustering.$k") || die "Error opening clustering output file.";
	    
	    # initialize/clean the variables
	    my %hash = ();
	    my $lineNo = 0;
	    my $clust = 0;

	    while($clust = <CO>)
	    {
		# hash on the cluster# and append the observation# 
		chomp($clust);
		if(exists $hash{$clust})
		{
		    $hash{$clust} .= " $lineNo";
		}
		else
		{
		    $hash{$clust} = $lineNo;
		}
		
		# increment the line number
		$lineNo++;
	    }
	    
	    close CO;

	    $W_M[$i][$k] = &error_measure(\@refmat, \%hash);
	    
	    unlink "$out_filename","$filename.clustering.$k";
	}

	unlink "$filename", "$filename.tree";
    }

    my @sum = ();
    my @gap = ();

    # Calculate summationOverB(log(Wkb*))
    my @E_W = ();
    $E_W[0] = 0;
    
    my @tmp_gap = ();
    for($i = 1; $i <= $K; $i++)
    {
	for($j = 1; $j <= $B; $j++)
	{
	    $sum[$i] += $W_M[$j][$i];
	}
	
	# Calculate Gap(k) = 1/B(summationOverB(log(Wkb*))) - log(Wk)
	$gap[$i] = $sum[$i]/$B - $W[$i];
	
	# for the graph
	$tmp_gap[$i-1] = $gap[$i];
	$E_W[$i-1] = $sum[$i]/$B;
    }

    # For plotting of graph of obs and exp log(Wk) vs. K
    my $fig3 = new GD::Graph::lines(600,480);
    
    $fig3->set_title_font(gdGiantFont,24);
    $fig3->set_x_label_font(gdGiantFont,14);
    $fig3->set_y_label_font(gdGiantFont,14);
    
    $fig3->set('x_label' => 'number of clusters k',
		   'y_label' => 'obs and exp log(Wk)',
		   'title' => 'obs and exp by k',
		   'transparent' => '0',
		   ) or warn $my_graph->error;

    $fig3->set_legend( 'Observed', 'Expected' );
    open(FILE,">$prefix.fig3.png");
    print FILE $fig3->plot([\@k,\@tmp_W,\@E_W])->png;
    close(FILE);


    #~~~~~~~~~~~~~~~~~~~~ Step 3: Calculation of Standard Deviation  ~~~~~~~~~~~~~~~~~~~~~~~~
    
    my @sd = ();
    my @s = ();

    # Calculate standard deviation sd for Wk*
    for($i = 1; $i <= $K; $i++)
    {
	for($j = 1; $j <= $B; $j++)
	{
	    $sd[$i] += ($W_M[$j][$i] - $sum[$i]/$B)**2;
	}
	
	$sd[$i] = sprintf("%.4f",sqrt($sd[$i]/$B));
	# Calculate the modified standard deviation to account for 
	# simulation error.
	$s[$i] = $sd[$i] * sqrt(1 + 1/$B);
	$s[$i] = sprintf("%.4f",$s[$i]);
    }

    my $ans = 1;
    my $tmp = 0;

    # Find the optimal #of k
    # Find the smallest possible value of k for which the following 
    # criterion is satisfied:
    # Gap(k) >= Gap(k+1) - s(k+1)
    for($i = 1; $i < $K; $i++)
    {
	if($gap[$i+1] < 0)
	{
	    $tmp = $gap[$i+1] + $s[$i];
	}
	else
	{
	    $tmp = $gap[$i+1] - $s[$i];
	}

	if($gap[$i] >= $tmp)
	{
	    $ans = $i;
	    last;
	}
    }

    # For plotting of graph of Gap vs. K    
    my $fig4 = new GD::Graph::lines(600,480);
    
    $fig4->set_title_font(gdGiantFont,24);
    $fig4->set_x_label_font(gdGiantFont,14);
    $fig4->set_y_label_font(gdGiantFont,14);
    
    $fig4->set('x_label' => 'number of clusters k',
		   'y_label' => 'Gap',
		   'title' => 'Gap by k',
		   'transparent' => '0',
		   ) or warn $my_graph->error;

    open(FILE,">$prefix.fig4.png");
    print FILE $fig4->plot([\@k,\@tmp_gap])->png;
    close(FILE);

    print "$ans\n";
    return $ans;
}


# Calculates the "Within Cluster Dispersion / Error Measure (Wk)"
sub error_measure
{
    # arguments
    my @matrix = @{(shift)};
    my %clustout = %{(shift)};

    #local variables
    my $i; 
    my $j;
    my @rownum;
    my $key;
    my $row1;
    my $row2;
    my $W = 0;
    my @D = ();
    my $tmp;

    # for each cluster
    foreach $key (sort keys %clustout)
    {
	$D[$key] = 0;

	@rownum = split(/\s+/,$clustout{$key});

	# for each instance in the cluster
	for($i = 0; $i < $#rownum; $i++)
	{
	    # for all the rows in the cluster
	    for($j = $i+1; $j <= $#rownum; $j++)
	    {
		# find the distance between the 2 rows of the matrix.
		$row1 = $rownum[$i];
		$row2 = $rownum[$j];

		# store the Dr value
		if(exists $d[$row1][$row2])
		{
		    $D[$key] += $d[$row1][$row2];
		}
		else
		{
		    $D[$key] += $d[$row2][$row1];
		}
	    }
	}

	$W += (1/(2 * ($#rownum + 1))) * $D[$key];
    }

    if($W != 0)
    {
	$W = sprintf("%.4f", log($W)/log(10));
    }

    return $W;
}

# Calculates distance between the given 2 rows
# currently implements using the Euclidean measure.
sub dist_euclidean
{
    # arguments
    my @i = @{(shift)};
    my @j = @{(shift)};

    # local variables
    my $a;
    my $dist = 0;
    my $retvalue = 0;

    # Euclidean measure 
    # square-root(summation on all j (xij - xi'j)^2 where i, i' are the rows indicies)
    for $a (0 .. $#i)
    {
	$dist += (($i[$a] - $j[$a])**2);
    }
    $dist = sqrt($dist);

    $retvalue = sprintf("%.4f",$dist);	
    return $retvalue;
}

sub dist_euclidean_sqr
{
    # arguments
    my @i = @{(shift)};
    my @j = @{(shift)};

    # local variables
    my $a;
    my $dist = 0;
    my $retvalue = 0;

    # Squared Euclidean measure 
    # summation on all j (xij - xi'j)^2 where i, i' are the rows indicies.
    for $a (0 .. $#i)
    {
	$dist += (($i[$a] - $j[$a])**2);
    }

    $retvalue = sprintf("%.4f",$dist);	
    return $retvalue;
}

sub dist_manhattan
{
    # arguments
    my @i = @{(shift)};
    my @j = @{(shift)};

    # local variables
    my $a;
    my $dist = 0;
    my $retvalue = 0;

    # Manhattan measure 
    # summation on all j (|xij - xi'j|) where i, i' are the rows indicies.
    
    for $a (0 .. $#i)
    {
       $dist += abs($i[$a] - $j[$a]);
    }

    $retvalue = sprintf("%.4f",$dist);	
    return $retvalue;
}

1;
__END__

=head1 NAME

Statistics::Gap - Perl extension for the "Gap Statistics"

=head1 SYNOPSIS

  use Statistics::Gap;
  &gap("GapPrefix", "Filename.txt", "manhattan", "agglo", 5, 3);

=head1 DESCRIPTION

    Given a dataset how does one automatically find the optimal number 
    of clusters that the dataset should be grouped into? - is one of the 
    prevailing problems. Statisticians Robert Tibshirani, Guenther Walther 
    and Trevor Hastie  propose a solution for this problem is a Techinal 
    Report named - "Estimating the number of clusters in a dataset via 
    the Gap Statistics". This perl module implements the approach proposed 
    in the above paper.

=head2 EXPORT

 "gap" function by default.

=head1 INPUT

=head2 Prefix
    
    The string that should be used to as a prefix while naming the 
    intermediate files and the .png files (graph files).

=head2 InputFile

    The input dataset is expected in a plain text file where the first
    line in the file gives the dimensions of the dataset and then the 
    dataset in a matrix format should follow. The contexts / observations 
    should be along the rows and the features should be along the column.

=head3 DistanceMeasure

    The Distance Measure that should be used.
    Currrently this module supports the following distance measure:
    1. Manhattan (string that should be used as an argument: "manhattan")
    2. Euclidean (string that should be used as an argument: "euclidean")
    3. Squared Euclidean (string that should be used as an argument: "squared")

=head3 ClusteringAlgorithm

    The Clustering Measures that can be used are:
    1. rb - Repeated Bisections [Default]
    2. rbr - Repeated Bisections for by k-way refinement
    3. direct - Direct k-way clustering
    4. agglo  - Agglomerative clustering
    5. graph  - Graph partitioning-based clustering
    6. bagglo - Partitional biased Agglomerative clustering

=head3 K value
    
    This is an approximate upper bound for the number of clusters that may be
    present in the dataset. Thus for a dataset that you expect to be seperated
    into 3 clusters this value should be set some integer value greater than 3.

=head3 B value
    
    Specifies the number of time the reference distribution should be generated
    Typical value would be 3.

=head1 OUTPUT

    The output returned is a single integer number which indicates the optimal
    number of clusters that the input dataset should be clustered into.

=head1 PRE-REQUISITES

    This module uses suite of C programs called CLUTO for clustering purposes. 
    Thus CLUTO needs to be installed for this module to be functional.
    CLUTO can be downloaded from http://www-users.cs.umn.edu/~karypis/cluto/

=head1 SEE ALSO

    http://citeseer.ist.psu.edu/tibshirani00estimating.html
    http://www-users.cs.umn.edu/~karypis/cluto/

=head1 AUTHOR

    Anagha Kulkarni, University of Minnesota Duluth
    kulka020 <at> d.umn.edu
	
    Guergana Savova, Mayo Clinic
    savova.guergana <at> mayo.edu

=head1 COPYRIGHT AND LICENSE

    Copyright (C) 2005-2006, Guergana Savova and Anagha Kulkarni

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

=cut
