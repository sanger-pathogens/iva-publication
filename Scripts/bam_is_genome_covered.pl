#!/usr/bin/env perl

use strict;
use warnings;

@ARGV == 4 or die "usage: $0 <in.bam> <% genome> <min depth each strand> <min mapped %>
";


my $bam = $ARGV[0];
my $min_percent_cov = $ARGV[1];
my $min_depth = $ARGV[2];
my $min_mapped_percent = $ARGV[3];

my ($reads, $mapped_reads) = reads_mapped($bam);
my $total_length = total_length_from_bam($bam);
my $fwd_cov = bam_to_bases_covered($bam, $min_depth, 0);
my $rev_cov = bam_to_bases_covered($bam, $min_depth, 1);
my $ok_bases = length_intersection($fwd_cov, $rev_cov);
my $depth_ok = 100 * $ok_bases / $total_length >= $min_percent_cov ? 1 : 0;
my $enough_mapped = 100 * $mapped_reads / $reads >= $min_mapped_percent ? 1 : 0;
my $ok = $depth_ok * $enough_mapped;

print "Ref_bases\t$total_length\n";
print "OK_bases\t$ok_bases\n";
print "Reads\t$reads\n";
print "Mapped\t$mapped_reads\n";
print "Pass_QC\t$ok\n";


sub length_intersection {
    my $h1 = shift;
    my $h2 = shift;
    my $total = 0;
    for (keys %$h1) {
        $total += 1 if defined $h2->{$_};
    }
    return $total;
}


sub total_length_from_bam {
    my $bam = shift;
    my $length = 0;
    open F, "samtools view -H $bam |" or die $!;
    while (<F>) {
        if (/^\@SQ.*\tLN:(\d+)\W/) {
            $length += $1;
        }
    }
    close F or die $!;
    $length > 0 or die "Got zero total length from BAM!";
    return $length;
}


sub bam_to_bases_covered {
    my $bam = shift;
    my $min_depth = shift;
    my $reverse_strand = shift;
    my %covered;

    my $flags = $reverse_strand ? '--rf 0x10' : '--ff 0x10';
    open F, "samtools mpileup -B $flags $bam | cut -f 1,2,4 |" or die $!;
    while (<F>) {
       chomp;
       my ($name, $pos, $depth) = split /\t/;
       if ($depth > $min_depth) {
           $covered{"$name.$pos"} = 1;
       }
    }
    close F or die $!;
    return \%covered;
}


sub reads_mapped {
    my $bam = shift;
    my $total_reads = -1;
    my $mapped_reads = -1;
    open F, "samtools flagstat $bam |" or die $!;
    while (<F>) {
        if (/^(\d+) \+ 0 in total/) {
            $total_reads = $1;
        }
        if (/^(\d+) \+ 0 mapped/) {
            $mapped_reads = $1;
        }

    }
    close F or die $!;

    if ($total_reads == -1 or $mapped_reads == -1) {
        die "Error getting read counts from bam $bam";
    }
    return ($total_reads, $mapped_reads);
}
