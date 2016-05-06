#!/usr/bin/perl

use Getopt::Long qw(:config posix_default no_ignore_case gnu_compat);

#######################################################
# interval2mask.pl
#
#
# [description]
# 
# This program takes sorted list of closed-intervals
# [[l0, r0], [l1, r1], ... , [ln, rn]]
# from stdin and performs the following operation --
# For each interval (li, ri), compute 
#    Li := (li - margin) / resolution
#    Ri := (ri + margin) / resolution
# and output the union of intervals
# [[L0, R0], [L1, R1], ... , [Ln, Rn]]
#
#   
# [usage]
#
# $ cat some_bed_file | grep chr1 | cut -f2,3 \
#   | ./interval2mask.pl --resolution=1000 --margin=1000
#
#######################################################

# set resolution (size of the bin) from command line option
my $res = 0, $margin = 0, $sep = "\n";
GetOptions('res=i' => \$res, 'margin=i' => \$margin);

if ($res > 0){
    # check the option (to avoid zero division)

    my $prev_left = -1, $prev_right = -1;
    # set decoy params
    
    foreach(<STDIN>){

	chomp;
	@entry = split(/\t/, $_);
	$left  = int(($entry[0] - $margin) / $res);
	$right = int(($entry[1] + $margin) / $res);       

	if(($left - 1) <= $prev_right){
	    # there is an overlap between two intervals
	    # or two intervals are adjacent to each other
	    $prev_right = $right;
	    # merge two intervals

	}else{ # there is no overlap	    
	    if($prev_right >= 0){
		print STDOUT $prev_left, $sep, $prev_right, "\n";
		# print the previous interval
	    }
	    $prev_left = $left;
	    $prev_right = $right;

	}
    }
    
    if($prev_left != -1){
	print STDOUT $prev_left, $sep, $prev_right, "\n";
	# print the last interval
    }
}
