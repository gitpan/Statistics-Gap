package Statistics::Gap;

use 5.008005;
use strict;
use warnings;
use POSIX qw(floor ceil);

require Exporter;
use AutoLoader qw(AUTOLOAD);

our @ISA = qw(Exporter);

our @EXPORT = qw( gap );

our $VERSION = '0.05';

# pre-requisites
use GD;
use GD::Text;
use GD::Graph::lines;
use GD::Graph::colour;
use Data::Dumper;

# global variable
my $distfuncref;
my @d = ();

# Estimates number of clusters that a given data naturally falls into
sub gap
{
    # input params
    my $prefix = shift;
    my $matrixfile = shift;
    my $distmeasure  = shift;
    my $clustmtd = shift;
    my $K = shift;
    my $B = shift;
    my $mtd = shift;
    my $perc = shift;

    my $i = 0;
    my $j = 0;
    my $c = 0;

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

    # Not a valid condition: 
    # If maximum number of clusters requested (k) is greater than the 
    # number of observations.
    if($K > $rcnt)
    {
	print STDERR "The K value ($K) cannot be greater than the number of observations present in the input data ($rcnt). \n";
	exit 1;
    }

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

    my @rmar = ((0) x $rcnt);
    my @cmar = ((0) x $ccnt);
    my $totalElem = 0;

    # calculate the pairwise distances for all possible unique pairs (n(n-1)/2)
    for($i = 0; $i < $rcnt; $i++)
    {
	# for all the rows in the cluster
	for($j = $i+1; $j < $rcnt; $j++)
	{
	    @row1 = @{$inpmat[$i]};
	    
	    @row2 = @{$inpmat[$j]};
	    
	    $d[$i][$j] = &$distfuncref(\@row1, \@row2);
	}

	# calculate the marinals and the total number of non-zero values
	for($c = 0; $c < $ccnt; $c++)
	{
	    if($inpmat[$i][$c] != 0)
	    {
		$rmar[$i] += 1;
		$cmar[$c] += 1;
		$totalElem += 1;
	    }
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
	die "Error running vcluster \n" unless $status==0;
	
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

    my @W_M = ();

    my @min = ();
    my @diff = ();

    my $n = 0;
    my $rand;
    my $tmp_rmar = 0;

    # 2 possible methods of generating reference distribution - uniform & proportional
    if($mtd eq "unif")
    {
	# Repeat the complete Reference Distribution generation procedure B times.
	for($i = 1; $i <= $B; $i++)
	{
	    my @refmat = ();

	    # initialize the reference matrix with "0"s.
	    for($j = 0; $j < $rcnt; $j++)
	    {
		my @tmp_array = ((0) x $ccnt);
		push @refmat, [ @tmp_array ]; 
	    }

	    # create the binary reference matrix using Uniform method
	    for($n = 0; $n < $rcnt; $n++)
	    {
		$tmp_rmar = $rmar[$n];
		while($tmp_rmar)
		{
		    $rand = rand($ccnt);      

		    if(!$refmat[$n][$rand])
		    {
			$refmat[$n][$rand] = 1;
			$tmp_rmar--;
		    }
		}
	    } #for n

	    # Calculate Wk* from the generated reference matrix which consists of:
	    # 1. Cluster the generated n observations
	    # 2. Calculate the Wk*
	    
	    # Write the matrix to a temporary file
	    my $filename = "tmp.ref." . time() . $i;
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
		die "Error running vcluster \n" unless $status==0;
		
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

	} # for $B
    }
    else # proportional
    {
	# calculate the range for each feature
	my @lower = ();
	my @upper = ();

	my $index = 0;
	my $hold = 0;
	my %hash_col = ();

	# sort the array of col marginals in ascending order
	# to minimize the the adjustment of ranges which is 
	# needed after a random number is assigned to a feature
	@cmar = sort {$a <=> $b} (@cmar);

	# create the ranges for all the features
	for($c = 0; $c <= $#cmar; $c++)
	{
	    if($cmar[$c])
	    {
		$lower[$index] = $hold + 1;
		$upper[$index] = $hold + $cmar[$c];
		$hold = $upper[$index];

		# to take care of the situations where a col 
		# marginal for a feature is zero.
		$hash_col{$index} = $c;
		$index++;
	    }
	}

	# Repeat the complete Reference Distribution generation procedure B times. 
	for($i = 1; $i <= $B; $i++)
	{
	    my @refmat = ();

	    # initialize the reference matrix with "0"s.
	    for($j = 0; $j < $rcnt; $j++)
	    {
		my @tmp_array = ((0) x $ccnt);
		push @refmat, [ @tmp_array ]; 
	    }

	    # create the binary reference matrix using Proportional method
	    for($n = 0; $n < $rcnt; $n++)
	    {
		$tmp_rmar = $rmar[$n];

		# copy the ranges to temp array
		# because these temp ranges will
		# be adjusted after every random
		# number is assigned to a feature
		my @tmp_low = @lower;
		my @tmp_upp = @upper;
		my $temp_totalElem = $totalElem;
		my $prev_ans = scalar(@cmar);

		# repeat row marginal times
		while($tmp_rmar)
		{
		    # translate the random generated over [0,$temp_totalElem-1]
		    # to [1,$temp_totalElem]
		    $rand = rand($temp_totalElem);
		    $rand++;

		    # find the feature# to which this random number belongs to.
		    # using binary search
		    my $min = 0;
		    my $max = $#tmp_upp;
		    my $col = 0;
		    
		    my $feat_ind = 0;

		    while($min <= $max)
		    {
			$col = floor(($min + $max)/2);

			if($rand == $tmp_upp[$col])
			{
			    $feat_ind = $hash_col{$col};
			    if($prev_ans <= $feat_ind)
			    {
				$feat_ind++;
			    }
			    $prev_ans = $feat_ind;
			    last;
			}
			elsif($rand < $tmp_upp[$col])
			{
			    if($rand >= $lower[$col])
			    {
				$feat_ind = $hash_col{$col};
				if($prev_ans <= $feat_ind)
				{
				    $feat_ind++;
				}
				$prev_ans = $feat_ind;
				last;
			    }
			    else
			    {
				$max = $col - 1;
			    }
			}
			else
			{
			    $min = $col + 1;
			}
		    } #while

		    $refmat[$n][$feat_ind] = 1;
		    $tmp_rmar--;
		    
		    # adjust the ranges
		    my $tmp_ans = $feat_ind;
		    while($tmp_ans < $#tmp_low )
		    {
			$tmp_low[$tmp_ans] = $tmp_low[$tmp_ans+1] - $cmar[$feat_ind];
			$tmp_upp[$tmp_ans] = $tmp_upp[$tmp_ans+1] - $cmar[$feat_ind];			
			$tmp_ans++;
		    }
		    
		    $temp_totalElem = $temp_totalElem - $cmar[$feat_ind];

		    $tmp_low[$#tmp_low] = undef;
		    $tmp_upp[$#tmp_upp] = undef;
		    $#tmp_low = $#tmp_low - 1;
		    $#tmp_upp = $#tmp_upp - 1;
		} # while
	    } #for n

	    # Calculate Wk* from the generated reference matrix which consists of:
	    # 1. Cluster the generated n observations
	    # 2. Calculate the Wk*
	    
	    # Write the matrix to a temporary file
	    my $filename = "tmp.ref." . time() . $i;
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
		die "Error running vcluster \n" unless $status==0;
		
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

	} # for $B
    }# else "proportional"

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
	
	$sum[$i] = sprintf("%.4f",$sum[$i]/$B) + 0;
	# Calculate Gap(k) = 1/B(summationOverB(log(Wkb*))) - log(Wk)
	$gap[$i] = sprintf("%.4f", $sum[$i] - $W[$i]) + 0;
	
	# for the graph
	$tmp_gap[$i-1] = $gap[$i];
	$E_W[$i-1] = $sum[$i];
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
	    $sd[$i] += ($W_M[$j][$i] - $sum[$i])**2;
	}
	
	$sd[$i] = sprintf("%.4f",sqrt($sd[$i]/$B)) + 0;
	# Calculate the modified standard deviation to account for 
	# simulation error.
	$s[$i] = $sd[$i] * sqrt(1 + 1/$B);
	$s[$i] = sprintf("%.4f",$s[$i]) + 0;
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
	    $tmp = $gap[$i+1] + $s[$i+1];
	}
	else
	{
	    $tmp = $gap[$i+1] - $s[$i+1];
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

    # Printing to the log file.
    open(LO,">$prefix.log") || die "Error opening the Log file ($prefix.log) in write mode.\n";    

    # Print the confidence intervals
    my $alpha = 0;
    my $low_conf_int = 0;
    my $upp_conf_int = 0;
    my @tmp = ();

    printf LO "%3s  %10s  %10s  %10s  %10s  %10s  %30s\n", "K", "Gap(k)", "log(W(k))", "log(W*(k))", "sd(k)", "s(k)", "$perc% Confidence Intervals";   
    printf LO "-" x 95 ."\n";
    for($i = 1; $i <= $K; $i++)
    {
	# Calculate a from 100(1-2a) = %
	$alpha = (1 - $perc/100)/2;

	for($j = 1; $j <= $B; $j++)
	{
	    $tmp[$j-1] = $W_M[$j][$i];
	}
	# sort in the numeric ascending order 
	@tmp = sort {$a <=> $b} (@tmp);    
	
	# Calculate lower bound = average of W*[floor(B*a)] and W*[ceil(B*a)]
	$low_conf_int = ($tmp[floor($B*$alpha)-1] + $tmp[ceil($B*$alpha)-1])/2;
	
	# Calculate upper bound = average of W*[floor(B*(1-a))] and W*[ceil(B*(1-a))]
	$upp_conf_int = ($tmp[floor($B*(1-$alpha))-1] + $tmp[ceil($B*(1-$alpha))-1])/2;
	
	printf LO "%3d  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %30s\n", $i, $gap[$i], $W[$i], $sum[$i], $sd[$i], $s[$i], "$low_conf_int - $upp_conf_int";
    }

    print LO "\nIndividual Dispersion values:\n";
    for($i = 1; $i <= $K; $i++)
    {
	print LO "K=$i\n";
	printf LO "%3s  %10s\n", "B", "log(W*)"; 
	printf LO "-" x 15 . "\n"; 
	for($j = 1; $j <= $B; $j++)
	{
	    $tmp[$j-1] = $W_M[$j][$i];
	    printf LO "%3s  %10s\n", "$j", "$W_M[$j][$i]"; 
	}
	print LO "\n";
    }

    close LO;


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
    foreach $key (keys %clustout)
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

	$W += (1/($#rownum + 1)) * $D[$key];
    }

    if($W != 0)
    {
	$W = sprintf("%.4f", log($W)/log(10)) + 0;
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

    $retvalue = sprintf("%.4f",$dist) + 0;	
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

    $retvalue = sprintf("%.4f",$dist) + 0;	
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

    $retvalue = sprintf("%.4f",$dist) + 0;	
    return $retvalue;
}

1;
__END__

=head1 NAME

 Statistics::Gap - Perl extension for the "Gap Statistic"

=head1 SYNOPSIS

 use Statistics::Gap;
 &gap("GapPrefix", "InputFile", "squared", "agglo", 5, 100, "unif", 90);

 OR

 use Statistics::Gap;
 &gap("GapPrefix", "InputFile", "squared", "agglo", 5, 100, "prop", 90);

 Input file is expected in the "dense" format -
 Sample Input file:
   
 6 5
 1       1       0       0       1
 1       0       0       0       0
 1       1       0       0       1
 1       1       0       0       1
 1       0       0       0       1
 1       1       0       0       1 	  	

=head1 DESCRIPTION

Given a dataset how does one automatically find the optimal number 
of clusters that the dataset should be grouped into? - is one of the 
prevailing problems. Statisticians Robert Tibshirani, Guenther Walther 
and Trevor Hastie  propose a solution for this problem in a Techinal 
Report named - "Estimating the number of clusters in a dataset via 
the Gap Statistic". This perl module implements the approach proposed 
in the above paper. 

If one tries to cluster a dataset (i.e. numerous observations described 
in terms of a feature space) into n groups/clusters and if we plot the
graph of Within Cluster (Dis)Similarity along Y-axis and Number of clusters
along X-axis then this graph generally takes a form of a elbow/knee depending
upon the measure on the Y-axis. The Gap Statistic seeks to locate this
elbow/knee because the value on the X-axis at this elbow is the optimal 
number of clusters for the data.	

NOTE: 
Gap Statistic uses reference distribution in the process of estimating
the number of clusters. The appropriate methodology for generation of this 
reference distribution is dependent on the data to be clustered. 
This module was implemented for data with following characteristics:
1. highly sparse - very few features occur in any given observation.
2. high multivariate dimensionality (i.e. large feature space)
3. binary feature frequency - feature either occurs or does not occur 
in an observation but rarely occurs multiple times in the same observation.

=head2 EXPORT

"gap" function by default.

=head1 INPUT

=head2 Prefix

The string that should be used to as a prefix while naming the 
intermediate files and the .png files (graph files).

=head2 InputFile

The input dataset is expected in "dense" matrix format.
The input dense matrix is expected in a plain text file where the first
line in the file gives the dimensions of the dataset and then the 
dataset in a matrix format should follow. The contexts / observations 
should be along the rows and the features should be along the column.

	eg:
      	6 5
        1       1       0       0       1
        1       0       0       0       0
        1       1       0       0       1
        1       1       0       0       1
        1       0       0       0       1
        1       1       0       0       1 	

The first line (6 5) gives the number of rows (observations) and the 
number of columns (features) present in the following matrix.
Following each line records the frequency of occurrence of the feature
at the column in the given observation. Thus features1 (1st column) occurs
once in the observation1 and infact once in all the other observations too 
while the feature3 does not occur in observation1.

=head2  DistanceMeasure

The Distance Measure that should be used.
Currrently this module supports the following distance measure:
1. Squared Euclidean (string that should be used as an argument: "squared")
2. Manhattan (string that should be used as an argument: "manhattan")
3. Euclidean (string that should be used as an argument: "euclidean")

=head2  ClusteringAlgorithm

The Clustering Measures that can be used are:
1. rb - Repeated Bisections [Default]
2. rbr - Repeated Bisections for by k-way refinement
3. direct - Direct k-way clustering
4. agglo  - Agglomerative clustering
5. graph  - Graph partitioning-based clustering
6. bagglo - Partitional biased Agglomerative clustering

=head2  K value

This is an approximate upper bound for the number of clusters that may be
present in the dataset. Thus for a dataset that you expect to be seperated
into 3 clusters this value should be set some integer value greater than 3.

=head2  B value

Specifies the number of time the reference distribution should be generated.

=head2  ReferenceGenerationMethod

1. Uniform - While generating the reference distribution, all the features
in the feature set have equal probability of being selected for the observation
under consideration.
    
2. Proportional - Each feature is assigned a probability of being selected
depending upon its frequency of occurrence in the observed data. Thus feature
distribution is taken into consideration while selecting the features for the 
reference distribution generation.

=head2  Percentage Confidence Interval

This parameter specifies the percentage confidence to be reported in the log file.
Since Statistics::Gap uses parametric bootstrap method for reference distribution 
generation, it is critical to understand the interval around the sample mean that
could contain the population ("true") mean and with what certainty.

=head1 OUTPUT

1. A single integer number at STDOUT which is the Gap Statistic's 
estimate of number of clusters present in the input dataset.
2. The PREFIX.log file contains the log of various values at different K values.
The first table in the file gives values like Gap(k), log(W(k)) etc. for every K value
experimented with.	
3. The PREFIX.fig*.png files are the graphical representations of different values which
help locate the knee/elbow. 

=head1 PRE-REQUISITES

1. This module uses suite of C programs called CLUTO for clustering purposes. 
Thus CLUTO needs to be installed for this module to be functional.
CLUTO can be downloaded from http://www-users.cs.umn.edu/~karypis/cluto/

2. Following Perl Modules
   1. GD   (http://search.cpan.org/~lds/GD-2.19/GD.pm)
   2. GD::Text     (http://search.cpan.org/~mverb/GDTextUtil-0.86/Text.pm)
   3. GD::Graph::lines     (http://search.cpan.org/~mverb/GDGraph-1.43/)
   4. GD::Graph::colour    (http://search.cpan.org/~mverb/GDGraph-1.43/Graph/colour.pm)

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
