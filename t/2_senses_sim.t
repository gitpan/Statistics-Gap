# 2_sense.t version 0.01
# (Updated 08/13/2005 -- Anagha)
#
# A script to run tests on the Statistics::Gap module.
# This test cases check for 2 sense input matrix.
# The following are among the tests run by this script:
# 1. Try loading the Statistics::Gap i.e. is it added to the @INC variable
# 2. Compare the answer returned by the Gap Statistics with the actual andwer

use strict;
use warnings;

use Test::More tests => 7;

BEGIN { use_ok('Statistics::Gap') };

# optimal number of clusters for the above input data
my $ans = 2;

my $result = 0;
$result = &gap("pre_2_sim", "sim", "t/2_senses_sim", "rb", "e1", 10, 15, "rep", 80, 4);

is($result, $ans, "Comparing Gap Statistics' answer ($result) with the actual optimal number of clusters ($ans) for the input data");

if(-e "pre_2_sim.gap.log")
{
	is("exists","exists");
}
else
{
	is("does not exist pre_2_sim.gap.log","exist");
}

if(-e "pre_2_sim.fig2.dat")
{
	is("exists","exists");
}
else
{
	is("does not exist fig2","exist");
}

if(-e "pre_2_sim.fig3a.dat")
{
	is("exists","exists");
}
else
{
	is("does not exist fig3a","exist");
}

if(-e "pre_2_sim.fig3b.dat")
{
	is("exists","exists");
}
else
{
	is("does not exist fig3b","exist");
}

if(-e "pre_2_sim.fig4.dat")
{
	is("exists","exists");
}
else
{
	is("does not exist fig4","exist");
}

unlink "pre_2_sim.gap.log", "pre_2_sim.fig2.dat","pre_2_sim.fig3a.dat","pre_2_sim.fig3b.dat","pre_2_sim.fig4.dat";

__END__
