#!/usr/bin/perl

if(scalar(@ARGV) == 0)
{
    &usage();
}

$file1 = $ARGV[0]; #input fasta file
$file2 = $ARGV[1]; #output fasta file

open (IN, "$file1") or die &usage();
open (OUT, ">$file2") or die &usage();
&initialize ();

while (<IN>)
{
	if (/>/)
	{	chomp;
		if ($. == 1)
		{
			print OUT "$_\n";
		}
		else {print OUT "\n$_\n";}
	}
	else
	{
		chomp;
		print OUT $_;
	}
}
close IN;

sub usage 
{
    print "\n\tUse this script to convert sequences on 60 format to one line format\n\n";
    print "\texample command:\n";
    print "\tParameter Order : input sequence file name, output sequence file name\n\n";
    print "\tperl seq_one_line.pl parameter Order mentioned above\n\n";
    print "\n\tPlease follow the same parameter order as mentioned in above command\n\n";
    exit(1);
}

sub initialize
{
	undef @array;
}
