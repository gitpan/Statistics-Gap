# 2_sense.t version 0.01
# (Updated 07/06/2005 -- Anagha)
#
# A script to run tests on the Statistics::Gap module.
# This test cases check for 2 sense input matrix.
# The following are among the tests run by this script:
# 1. Try loading the Statistics::Gap i.e. is it added to the @INC variable
# 2. Compare the answer returned by the Gap Statistics with the actual andwer

use strict;
use warnings;

use Test::More tests => 5;

BEGIN { use_ok('Statistics::Gap') };

# optimal number of clusters for the above input data
my $ans = 5;

my $result = 0;
$result = &gap("pre_mix", "t/mix", "squared", "agglo", 7, 30, "prop",80);

is($result, $ans, "Comparing Gap Statistics' answer ($result) with the actual optimal number of clusters ($ans) for the input data");

if(-e "pre_mix.fig2.png")
{
	is("exists","exists");
}
else
{
	is("does not exist fig2","exist");
}

if(-e "pre_mix.fig3.png")
{
	is("exists","exists");
}
else
{
	is("does not exist fig3","exist");
}

if(-e "pre_mix.fig4.png")
{
	is("exists","exists");
}
else
{
	is("does not exist fig4","exist");
}

unlink "pre_mix.fig2.png","pre_mix.fig3.png","pre_mix.fig4.png","t/mix.tree";

__END__
