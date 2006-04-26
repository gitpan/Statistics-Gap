# 2_sense.t version 0.02
# (Updated 04/26/2006 -- Anagha)
#
# A script to run tests on the Statistics::Gap module.
# This test cases check for 2 sense input matrix.
# The following are among the tests run by this script:
# 1. Try loading the Statistics::Gap i.e. is it added to the @INC variable
# 2. Compare the answer returned by the Gap Statistics with the actual answer

use strict;
use warnings;

use Test::More tests => 7;

BEGIN { use_ok('Statistics::Gap') };

# optimal number of clusters for the above input data
my $ans = 2;

my $result = 0;
$result = &gap("pre_2", "vec", "t/2_senses", "rb", "e1", 5, 30, "rep", 80, 4);

is($result, $ans, "Comparing Gap Statistics' answer ($result) with the actual optimal number of clusters ($ans) for the input data");

if(-e "pre_2.gap.log")
{
	is("exists","exists");
}
else
{
	is("does not exist pre_2.gap.log","exist");
}

if(-e "pre_2.fig2.dat")
{
	is("exists","exists");
}
else
{
	is("does not exist fig2","exist");
}

if(-e "pre_2.fig3a.dat")
{
	is("exists","exists");
}
else
{
	is("does not exist fig3a","exist");
}

if(-e "pre_2.fig3b.dat")
{
	is("exists","exists");
}
else
{
	is("does not exist fig3b","exist");
}

if(-e "pre_2.fig4.dat")
{
	is("exists","exists");
}
else
{
	is("does not exist fig4","exist");
}

unlink "pre_2.gap.log", "pre_2.fig2.dat","pre_2.fig3a.dat","pre_2.fig3b.dat","pre_2.fig4.dat";

__END__
