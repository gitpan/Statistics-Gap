use 5.008005;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'Statistics::Gap',
    VERSION	      => '0.02',
    VERSION_FROM      => 'lib/Statistics/Gap.pm',
    PREREQ_PM         => {GD => 2.19,
			  GD::Text => 0.86,
			  GD::Graph::lines => 1.15,
			  GD::Graph::colour => 1.10},
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'lib/Statistics/Gap.pm', # retrieve abstract from module
       AUTHOR         => 'Anagha Kulkarni <kulka020@d.umn.edu>, Guergana Savova <savova.guergana@mayo.edu>') : ()),
);
