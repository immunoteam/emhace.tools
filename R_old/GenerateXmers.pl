#!/usr/bin/perl

# GenerateXmers.pl
# Generates xmers from text read from standard input
# Written to be used by GenerateXmers R function
# standalone usage: cat textfile.txt | GenerateXmers.pl merlength

use strict;
use warnings;
my $merlngth = shift @ARGV;

while(<STDIN>) {
	chomp;
	my $to_pos = length($_);
	foreach my $strt (0..($to_pos-$merlngth)) {
		print(substr($_, $strt, $merlngth), "\n");
	}
	print("--\n")
}
