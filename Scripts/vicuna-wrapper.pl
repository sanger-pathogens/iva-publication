#!/usr/bin/env perl

use strict;
use warnings;
use File::Spec;
use Getopt::Long;

my %options = (
    min_identify => 90,
);

my $options_ok = GetOptions(\%options,
    'help',
    'min_identify=i',
);

if (scalar @ARGV != 5 or !($options_ok)) {
    print STDERR "usage: $0  <reads_1.fq> <reads_2.fq> <outdir> <insert min> <insert max>\n";
    exit(1);
}

my $reads1 = File::Spec->rel2abs($ARGV[0]);
my $reads2 = File::Spec->rel2abs($ARGV[1]);
my $outdir = $ARGV[2];
my $insert_min = $ARGV[3];
my $insert_max = $ARGV[4];
my $reads_dir = File::Spec->rel2abs($outdir) . '/Reads/';

my $config_file = 'config.txt';
mkdir $outdir or die $!;
chdir $outdir or die $!;
mkdir 'Reads' or die $!;
symlink $reads1, 'Reads/reads_1.fastq';
symlink $reads2, 'Reads/reads_2.fastq';

open F, ">$config_file" or die $!;

print F <<TXT;
============================================================================
=				 Configuration file of Diversifier
============================================================================

/* -----------------------------------------------------------------------
 *							Trimmer
 * -----------------------------------------------------------------------
 */

// The path to the Fasta file storing vectors used for trimming; [default ""]

// minimum prefix/suffix length of a read that matches a vector
// for trimming to be applied; [default 7]
//minMSize\t9

// minimum length of an internal sub-read that matches a vector
// for the read to be discarded; has to be shorter than any read
// [default 15]
// minInternalMSize\t15

// a sub-read is considered internal if the distance between either of its
// ends and the ends of the full read exceeding this value; [default 4]
//maxOverhangSize\t2

// minimum length of a read that will be retained after trimming;
// [default 25]
// minReadSize\t25

/* -----------------------------------------------------------------------
 * 							Profiler
 * -----------------------------------------------------------------------
 */

// Fasta file storing Multiple Sequence Alignment (MSA) of HIV
// genomes from database
//MSAFileName\tresult/LANL-HIV-1B-07082011/hiv-1B-cleaned.algn

// number of bins for dividing MSA, in the range of [10, 256] [default 20]
//binNumber\t20

// [default 15]
//kmerLength\t15

// maximum Hamming variation allowed in each kmer [default 1]
//maxHD\t1

// minimum spanning of kmers on a read to call a valid mapping [default 75]
//minSpan\t80

// number of blocks for dividing kmer indices
// has to be >= MaxHD and <= KmerLength, [default 5]
//blockNumber\t5

// output file storing IDs of reads that have been mapped to bins;
// [default ""]
//rMapFileName\toutput/rmap.txt

/* -----------------------------------------------------------------------
 * 							Contiger
 * -----------------------------------------------------------------------
 */

// word size for shingling [default 12]
//w1\t12

// word size for super-shingling [default 5]
//w2\t8

// maxmum % of divergence between read & consensus [default 10]
//Divergence\t8

// number of base pairs that can be ignored towards either end of a read,
// this accounts for insufficient trimming, PCR artifacts,
// sequencing errors, etc. [default 4]
//max_read_overhang\t2


/* The following 3 parameters are used to determine the reliable interval
 * start/end positions of a consensus
 */

// minimum weight of a profile column [default 5]
//min_profile_col_weight\t5

// minimum percentage ratio between the weight of the consensus base and
// the total weight of the profile column [default 85]
//min_consensus_base_ratio\t85

// max length of unreliable region in either end of the consensus [default 10]
//max_contig_overhang\t10

/* The following 2 parameters are used to determine low frequent variants
 * of a contig, which will be removed b4 aligning two contigs
 */
// min frequency of length polymorphic region to be considered for contig
// alignment [default 5]
//min_perc_polymorphism\t5

// maximum length of any variant that will be removed b4 contig alignment
// [default 20]
//max_variant_len\t20

// seed kmer length to computer overlap between two contigs
// has to be in the range of [9, 16] [default 12]
//seed_kmer_len\t12

// minimum length of overlap to merge two contigs in lieu of insufficient
// paired links [default 25]
//min_contig_overlap\t25

// minimum number of paired links to merge two contigs in lieu of
// insufficient overlap: [seed_kmer_len, min_reliable_overlap)
// [default 3]
//min_contig_links\t2

// minimum percent identity to merge two contigs
min_identify\t$options{min_identify}

/* -----------------------------------------------------------------------
 *			General Parameters for Assembly
 * -----------------------------------------------------------------------
 */

// Input folder for paired fastq files, note dir has to end with '/' [default ""]
pFqDir\t$reads_dir

// Input folder for non-paired fastq files [default ""]
//npFqDir\tinput/single/

// Input folder for paired fasta files [default ""]
//pFaDir\t""

// Input folder for non-paired fasta files [default ""]
//npFaDir\t""

// Number of total reads to be processed per batch [default 2M]
//batchSize\t1000000

// [default -1]
LibSizeLowerBound\t$insert_min
// [default -1]
LibSizeUpperBound\t$insert_max

// min length of contig that will be included in the output [default 300]
//min_output_contig_len\t1500

// output directory [default "./"]
//outputDIR\t/seq/viral/analysis/xyang/TMP/454/2599/
TXT


close F;

exec "vicuna $config_file" or die $!;

