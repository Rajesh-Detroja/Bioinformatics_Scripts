#!/usr/bin/perl
 # for computing simple statistics about an assembly given just the
# scaffold/contig fastas, just give it the fasta files representing the
# scaffolds as input (but make sure the gapchar is set appropriately for your
# data
use strict;

use Getopt::Long;

my($gapchar, $min_contig_size);
GetOptions (
    "gapchar=s" => \$gapchar,
    "min_contig_size=s" => \$min_contig_size,
            );
if (!defined $gapchar) {
    $gapchar='N';
}
my ($scaffold, $seqtext);
my (@contig_info, @scaffold_info);
my $num_contigs = 0;
my $short_contigs = 0;
my $num_scaffolds = 0;
my $contig_length_sum = 0;
my $scaffold_length_sum = 0;
my $seqtext;
my $max = 0;
my $min = 1000000000;
my $length;

#get the header and newline-less seqtext foreach scaffold and
#handle it
while (<>) {
    if (/^>(\S+)/) {
        $num_scaffolds++;
        if (defined $seqtext) {
                $length = length($seqtext);
                if ($length > $max){
                        $max = $length;
                }
                if ($length < $min){
                        $min = $length;
                }
                &handle_scaffold($scaffold, $seqtext);
        }
        $scaffold = $1;
	$seqtext = "";
    }
    else {
        chomp;
        $seqtext .= $_;
    }
}

$length=length($seqtext);
if ($length > $max){
        $max=$length;
}
if ($length < $min){
        $min = $length;
}

&handle_scaffold($scaffold, $seqtext);

print "num_scaffolds = ", $num_scaffolds, "\n",
"num_contigs = ", $num_contigs,  "\n",
(($short_contigs > 0) ? "num_short_contigs = ". $short_contigs. "\n" : ""),
        "total genome length incl. gaps = ", $scaffold_length_sum,  "\n",
        "total genome length w/o gaps = ", $contig_length_sum,  "\n",
        "avg_scaffold_size incl. gaps = ", $scaffold_length_sum / $num_scaffolds, "\n",
        "avg_scaffold_size w/o gaps = ", $contig_length_sum / $num_scaffolds, "\n",
    "avg_contig_size = ", $contig_length_sum/ $num_contigs, "\n";

@contig_info = sort {$b->{length} <=> $a->{length}} @contig_info;
@scaffold_info = sort {$b->{length} <=> $a->{length}} @scaffold_info;
my $cumlength = 0;
foreach my $contig (@contig_info) {
    $cumlength += $contig->{length};
    #print "ctg length " . $contig->{length} . " cum length " . $cumlength . " contig_length_sum " . $contig_length_sum . "\n";
    if ($cumlength >= $contig_length_sum / 2) {
        print "Contig N50 = " . $contig->{length} . "\n";
        last;
    }
}
$cumlength = 0;
foreach my $scaffold (@scaffold_info) {
    $cumlength += $scaffold->{length};
    if ($cumlength >= $scaffold_length_sum / 2) {
        print "Scaffold N50 = " . $scaffold->{length} . "\n";
        last;
}
}
print "Max scaffold size = $max\n";
print "Min scaffold size = $min\n";

sub handle_scaffold {
    my ($scaffold, $seqtext) = @_;
    my $scaffold_contigs = 0;
    my $scaffold_contig_length_sum = 0;
    while ($seqtext =~ /([^$gapchar]+)/g) {
        my $contig_length = length($1);
        if (defined $min_contig_size && $contig_length < $min_contig_size) {
            $short_contigs++;
            next;
        }
        $scaffold_contigs++;
        $scaffold_contig_length_sum += $contig_length;
        push @contig_info, {
            length => $contig_length
            };
    }
    my $scaffold_length = length($seqtext);
    push @scaffold_info, {
        name => $scaffold,
        length => $scaffold_length,
        num_contigs => $scaffold_contigs,
        contig_length_sum => $scaffold_contig_length_sum,
    };

    #print $scaffold, length($seqtext), $scaffold_contigs, $scaffold_contig_length_sum, "\n";
    $num_contigs += $scaffold_contigs;
    $contig_length_sum += $scaffold_contig_length_sum;
    $scaffold_length_sum += $scaffold_length;
}
