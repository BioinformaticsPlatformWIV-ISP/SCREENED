#!/usr/bin/env perl

=head1 USAGE

perl screened.pl

	-ref <input fasta file containing reference genomes>
	-config <input file containing method information> 
	-method <string that contains comma-separated list of method names to be evaluated>
	-cluster <string denoting whether amplicons should be clustered: "no", "greedy" or "exhaustive">
	-out_detail <name of output file containing detailed information>
	-out_summary <name of output file containing summary information>
	-out_cluster <name of output file containing cluster information>
	-deg <string denoting whether degenerate fragments should be considered: "yes" or "no">
	-ext_am <string denoting whether amplicons should be extende: "yes" or "no">
	-ext_frag <string denoting whether fragments should be extende: "yes" or "no">

=head1 VERSION

Current version is 1.0 (first 'stable release')
SCREENED ('polymeraSe Chain Reaction Evaluation through largE-scale miNing of gEnomic Data') Copyright (C) 2017 Kevin Vanneste (Bioinformatics Platform, Scientific Institute of Public Health - Belgium)
http://wiv-isp.be
Contact: bioit@wiv-isp.be

=head1 DESCRIPTION

Script that evaluates the usability of using a PCR "method" (consisting out of an amplicon and "fragments": forward and reverse primer to amplify the amplicon, and probe
to detect the amplicon) for a large set of sequences (typically full reference genomes but can also be smaller sequences that are to be tested).

For each PCR method to be evaluated, the script first searches for the method amplicon in the list of reference genomes using BLAST.
The script currently assumes that the top hit is the correct one. Exhaustive testing has indicated this to be a reasonable assumption.
If no hit is found, a statement "FALSE" is returned for the amplicon sequence and propagated throughout the script and output files.
Even if a small sequence hit is found somewhere in the genome not corresponding to the method amplicon, downstream tests (see further) will ensure
no false positives are propagated in the eventual output. If amplicon extension is used, the found amplicons will be extended in both the 5' and 3'
direction if, and only if, the found amplicon in the reference genome is shorter than the method amplicon. 

For each PCR method to be evaluated, the script then searches in these found amplicons whether the provided forward and reverse primer, plus the probe,
for the method will successfully anneal.
This is done by using BLASTN-SHORT tweaked specifically for these short fragment searches. If degenerate fragments should be considered, all possible
variations are considered. It is assumed that the top hit is the correct one. Exhaustive testing has indicated this to be a reasonable assumption, especially
if all degenerate fragments are considered.
If fragment extension is used, the found forward and reverse primer, and the probe, will be extended in both the 5' end 3' direction if, and only if, either of those found fragments
in the amplicon is shorter than the method fragments.
Three selection criteria are then used to decide whether these top hits for each fragment will lead to 
successful annealing (the fifth criterion is however not used anymore - it is preserved here for legacy reasons):
1. The aligned part for each fragment cannot contain more than 20% mismatches. Mismatches are recalculated after BLAST to account for the redundancy of the
IUPAC nucleotide characters because this is not natively supported by BLAST. If all degenerate fragments are considered, the latter will however not pose problems.
2. The aligned part for each fragment should be at least 80% of the length of the total fragment.
3. For the forward and reverse primers, there can be no mismatch in the last 5 bases at the 3' end.
Only if the 3 fragments pass all these selection criteria, it is assumed that the PCR method works for the reference genome in question.

When more than one method is tested, also referred to as method groups, a method group will be considered successful for a reference genome if at least
one of the methods succeeds.

Additionally, the script can also perform sequence clustering of all found amplicons for each PCR method. If performed, this can be done either
greedy (requires USEARCH, see dependencies) or exhaustive (requires MUSCLE, see dependencies). Exhaustive clustering is guaranteed to find the optimal solution
but can take very long for large sets of amplicons, whereas greedy clustering runs a lot faster and testing has indicated that the most optimal solution
is also almost always found.

=head1 COMMAND LINE OPTIONS

=over

=item -ref <input fasta file containing reference genomes>

A standard fasta file containing a list of reference genomes in nucleotide format (or large sequences to be evaluated otherwise).
Ensure all sequence headers are unique.
Also ensure that the nucleotide sequence only consists out of IUPAC nucleotide characters, excluding 'U'!

=item -config <input file containing method information> 

A tab-delimited config file containing all the required information for the PCR methods.
Each line should contain the information for one method in the following manner: 
	
	The first column contains a string that identifies the method. Avoid special characters!
	The second column contains the sequence of the forward primer. 
	The third column contains the sequence of the reverse primer. 
	The fourth column contains the sequence of the probe. 
	The fifth column contains be the sequence of the amplicon.

As for the reference genomes, the method information should only consist out of IUPAC nucleotide characters, excluding 'U'!
IMPORTANT: for forward primer, reverse primer, probe, and amplicon; all sequence information should be given in 5'->3' direction.
In practice, this means that the forward primer, amplicon, and probe, can be given as is but the reverse primer should be reverse complemented first.
Each PCR method to be evaluated (i.e. the "-method"" argument) should be presented in this file. More methods can be present than are evaluated in any one given run.

=item -method <string that contains comma-separated list of method names to be evaluated>

A string that denotes which PCR methods constitute a method group to be evaluated. Method names should have the corresponding fragment plus amplicon information
present in the config file (i.e. the "-config" argument). For a reference genome (i.e. "-fasta" argument), if any of the methods specified here is successful,
the method group is considered to be successful.

=item -cluster <string denoting whether amplicons should be clustered: "no", "greedy" or "exhaustive">

String that denotes the type of clustering:
"no" for no clustering (an empty output file will still be created).
"greedy" for clustering based on usearch. Usearch should hence be installed and available in the environment (see dependencies). 
Greedy clustering speeds up the computation considerably, but is not guaranteed to deliver an optimal solution. Testing has however indicated
the optimal solution is almost always found in practice.
"exhaustive" for exhaustive clustering. Does an all-versus-all global alignment of found amplicons based on muscle. Muscle should hence be installed and
available in the environment (see dependencies). Exhaustive clustering can take very long for large sets of amplicons, but is guaranteed to deliver an optimal solution.

=item -out_detail <name of output file containing detailed information>

Output file name that contains detailed information for each method (file contains a header line that explains the contents of this file).

=item -out_summary <name of output file containing summary information>

Output file name that contains summary information for the method (pair) to be evaluated
(file contains a header line that explains the contents of this file).

=item -out_cluster <name of output file containing cluster information>

Output file name that contains clustered amplicon information (file contains a header line that explains the contents of this file).
This file will also be created but be empty if no clustering was performed.

=item -deg <string denoting whether degenerate fragments should be considered: "yes" or "no">

String that denotes whether all degenerate fragments should be considered.
"no" means the fragment will be used as is (i.e. forward primer, reverse primer, and probe). "yes" means that all possible versions of the degenerate fragments
will be blasted and only the best one will be reported.
Not considering degenerate fragments can be dangerous if your fragments contain many degenerate IUPAC-characters, as these are considered as mismatches by BLAST.
Because these characters break up the longest possible word size recognized in your fragment, it may be that a random kmer is found someplace else in the amplicon.
Considering degenerate fragments will prevent this but will also increase running time.

=item -ext_am <string denoting whether amplicons should be extended: "yes" or "no">

Stringr that denotes whether amplicons can be extended at either the 5' and/or 3' end. "no" means no extension will be performed.
If set to "yes", found amplicons in the reference genome will be extended if, and only if, the found amplicon is shorter than the PCR method
amplicon. This can be required in cases where either the last or penultimate base of your amplicon is a mismatch. The BLAST algorithm will not consider this as part
of the alignment, and it will not be reported although manual inspection may still reveal it should be included. It is always advised to run the program first
with -ext_am set to "no"" to see if this may be a problem in your dataset first.

=item -ext_frag <string denoting whether fragments should be extended: "yes" or "no"">

String that denotes whether fragments (all 3: forward and reverse primer, and the probe) can be extended at either the 5' and/or 3' end. "no" means no extension will be performed.
If set to "yes", found fragments in the amplicon will be extended if, and only if, the found fragment is shorter than the PCR method
fragment. This can be required in cases where either the last or penultimate base of your fragment is a mismatch. The BLAST algorithm will not consider this as part
of the alignment, and it will not be reported although manual inspection may still reveal it should be included. It is always advised to run the program first
with -ext_frag set to 0 to see if this may be a problem in your dataset first.

=item -mismatch <maximum percentage of allowed mismatches>

Floating point number between 0 and 100 that gives the maximum percentage of allowed mismatches in annealing sites. This is applied to all three fragments (forward primer, reverse primer and probe).

=item -length <minimum percentage of alignment length>

Floating point number between 0 and 100 that gives the minimum percentage of the alignment length in annealing sites. This is applied to all three fragments (forward primer, reverse primer and probe).

=item -end <length of bases at 3' end that cannot have a mismatch>

Integer betweeen 0 and the total fragment length of the smallest method fragment that should be considered to check for mimmatches at the 3' end for the forward and reverse primer. 
This does NOT apply to the probe.

=item -bases <maximum number of mismathces at 3' end>

Integer between 0 and the total fragment length of the smallest method fragment that denotes the number of mismatches which is considered to prevent elongation for the length defined by the '-end' parameter.
This does NOT apply to the probe.

=back

=head1 DEPENDENCIES

=over

=item PERL

Perl should be available in the environment. Version v5.18.2 is tested and supported. Higher versions should also work in theory.

=item BLAST

BLAST should always be available in the environment. Version 2.2.30 is tested and supported. All versions of the new BLAST ("BLAST+")
should work in theory.

=item USEACH

Usearch should be available in the environment if greedy clustering is performed. Rename the binary to 'usearch' as the default binary name typically also
contains version information etc.
Version 8.1.1861 is tested and supported. Higher versions should also work in theory.

=item MUSCLE

Muscle should be available in the environment if exhaustive clustering is performed. Rename the binary to 'muscle' (although this is the default binary name).
Version 3.8.31 is tested and supported.

=back

=head1 BUG REPORTS

No known bugs at the moment.
Please send all bug reports and/or feature requests to bioit@wiv-isp.be.

=head1 COMMENTS

This program started out as a small idea/project to quickly test whether it would be possible to use BLAST to search in large sets of reference genomes whether 
a PCR method works through blasting the PCR method information (amplicon and fragments: forward primer, reverse primer, and probe) against the reference genome sequence.
It turned out to work remarkably well and successive layers of complexity were added (support for clustering, criteria checks, amplicon extension, fragment extension...)
to the extent a fully working program/paper was created. Because there was not a clear-cut design at the beginning of the project, this history of adding successive layers 
of complexity is also present within the code unfortunately. Although code refactoring has been attempted to make the design more straightforward and clear, some
legacy artifacts are still present. Time (and popularity of the tool...) permitting, a novel version might be created in Python communicating directly with an underlying
SQL database to increase flow.

=head1 CITATION

People using this program for their work are asked to cite the corresponding publication. Publication will be added here as soon as it is accepted somewhere.

=cut

use strict;
use warnings;
use Data::Dumper;    # Used for debugging - following up data structures
use Digest::MD5 qw(md5_hex);

#====================================================================
# MAIN PROGRAM BODY
#====================================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Initialization
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Check that arguments are given, if not print help
if ( !$ARGV[0] ) {
    print STDERR `pod2text $0`;
    exit;
}

# Load and make sure all required input parameters are accounted for
# There are also some simple checks but they are not really exhaustive
my %params = @ARGV;
unless ( exists( $params{'-ref'} ) )         { die "Please specify an input fasta file with all your reference genomes!\n" }
unless ( exists( $params{'-config'} ) )      { die "Please specify an input config file!\n" }
unless ( exists( $params{'-out_detail'} ) )  { die "Please specify the name of your detailed output file!\n" }
unless ( exists( $params{'-out_summary'} ) ) { die "Please specify the name of your summary output file!\n" }
unless ( exists( $params{'-out_cluster'} ) ) { die "Please specify the name of your cluster output file!\n" }
unless ( exists( $params{'-method'} ) )      { die "Please specify a comma-separated list with methods to be tested!\n" }
unless ( exists( $params{'-deg'} ) )         { die "Please specify whether all degenerate fragments should be considered (yes or no)!\n" }
unless ( $params{'-deg'} eq "no" || $params{'-deg'} eq "yes" ) { die "Input value for degeneracy not recognized (should be yes or no)!\n" }
unless ( exists( $params{'-cluster'} ) ) { die "Please specify a clustering method (no, greedy, or exhaustive)!\n" }
unless ( $params{'-cluster'} eq "no" || $params{'-cluster'} eq "greedy" || $params{'-cluster'} eq "exhaustive" ) { die "Input value for clustering not recognized (should be no, greedy, or exhaustive)!\n" }
unless ( $params{'-ext_am'} eq "no" || $params{'-ext_am'} eq "yes" )  { die "Please specify whether amplicons can be extended (yes or no)!\n" }
unless ( $params{'-ext_frag'} eq "no" || $params{'-ext_frag'} eq "yes" )  { die "Please specify whether fragments can be extended (yes or no)!\n" }
unless ( exists( $params{'-mismatch'} ) ) { die "Please specifiy a percentage for the maximum number of allowed mismatches in fragment annealing sites (value between 0 and 100)!\n" }
unless ( $params{'-mismatch'} >= 0 && $params{'-mismatch'} <= 100 ) { die "Please specify a value between 0 and 100 for the '-mismatch' parameter!\n" }
unless ( exists( $params{'-length'} ) ) { die "Please specifiy a percentage for the minimum length of the alignment in fragment annealing sites (value between 0 and 100)!\n" }
unless ( $params{'-length'} >= 0 && $params{'-length'} <= 100 ) { die "Please specify a value between 0 and 100 for the '-length' parameter!\n" }
unless ( exists( $params{'-end'} ) ) { die "Please specifiy the number of bases (integer) counting from the 3' end that should be considered for checking mismatches at the 3' end for the forwared and reverse primer (please see manual)!\n" }
unless ( $params{'-end'} =~ /^[+]?\d+$/ ) { die "Please specify a positive integer or 0 for the '-end' parameter!\n" }
unless ( exists( $params{'-bases'} ) ) { die "Please specifiy the number of mismatches (integer) in the region specified by the '-end' parameter that is considered to prevent elongation (please see manual)!\n" }
unless ( $params{'-bases'} =~ /^[+]?\d+$/ ) { die "Please specify a positive integer or 0 for the '-bases' parameter!\n" }
#print Dumper(\%params); exit;#print Dumper(\%params); exit;

# Load IUPAC database first
my %IUPAC = loadIUPAC();

# Process reference genomes and config file
my %refGenomes = &fasta2hash( $params{'-ref'} );
my %info       = &processConfig( $params{'-config'} );
my %methods    = &processMethods( $params{'-method'} );

#print Dumper(\%refGenomes); exit;
#print Dumper(\%info); exit;
#print Dumper(\%methods); exit;

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Find all amplicons in each reference
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Returns all found amplicons
my %foundAmplicons = ();
foreach my $method ( keys %methods ) {
    &blastAmplicon($method);
}

#print Dumper(\%foundAmplicons); #exit;

# First remove all gaps from the found amplicons for sequencing as these gaps can negatively influence both the fragment searches and clustering
# (Although some manual checks indicate there are in fact only seldom gaps within the found amplicons)
# IMPORTANT: foundAmplicons hash is hence updated within for the sequence of the aligned part by removing gaps (can be dangerous if this is not accounted by people modifying this script)!!!
&removeGaps();

#print Dumper(\%foundAmplicons); #exit;

# Extend the amplicons at both ends to account for gaps at the end of the amplicon alignment for which it is impossible to correct using traditional blast
# Note that this makes associated blast parameters invalid (the crucial ones are updated but if you don't check and start checking on others it may lead to erroneous results)!
if ( $params{'-ext_am'} eq "yes" ) {
    &extendAmplicons();
}

#print Dumper(\%foundAmplicons);# exit;

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Cluster the amplicons
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Cluster the amplicons by sequence similarity
my %foundClusters = ();
my %maxClusters   = ();
foreach my $method ( keys %methods ) {
    if ( $params{'-cluster'} eq "no" ) { next }
    elsif ( $params{'-cluster'} eq "greedy" ) {
        &clusterAmpliconGreedy($method);
    } elsif ( $params{'-cluster'} eq "exhaustive" ) {
        &clusterAmpliconExhaustive($method);
    }
}

#print Dumper(\%foundClusters); exit;
#print Dumper(\%maxClusters); #exit;

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Find all fragments in each amplicon
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Returns all found fragments but only for the top hit amplicon in each reference (hard-coded in &blastFragments subroutine!)
my %foundFragments = ();
my @searchFragments      = qw(forward reverse probe);
foreach my $method ( keys %methods ) {
    foreach my $fragment (@searchFragments) {
        &blastFragment( $method, $fragment );
    }
}

#print Dumper(\%foundFragments); exit;

# For all blast results, sort the results for the degenerate fragment searches
my %sortedFragments = ();
&sortVariantScores();

#print Dumper(\%sortedFragments);# exit;

# For all the fragment blast results, allow extension at the 5' end  for the forward and reverse primer
if ( $params{'-ext_frag'} eq "yes" ) {
    &extendFragments();
}

#print Dumper(\%sortedFragments);# exit;

# Recalculate the mismatch scores for the found fragments to account for the degeneracy of the genetic code
# This is still required in version 0.6 with the introduction of searching all possible variants because this does not always necessarily happen (if -deg flag is set to no)
# plus the amplicon sequence itself may also contain degenerate nucleotides!
my %foundFragmentsUpdated = %sortedFragments;
&recalculateMismatch();

#print Dumper(\%foundFragmentsUpdated); exit;

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Determine whether the fragments will anneal properly
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#First create a new hash than will contain the scores assuming everything has score 1
my %scoredFragments = ();
foreach my $ref ( keys %foundFragmentsUpdated ) {
    foreach my $method ( keys %{ $foundFragmentsUpdated{$ref} } ) {
        foreach my $fragment ( keys %{ $foundFragmentsUpdated{$ref}{$method} } ) {
            foreach my $variantHit ( keys %{ $foundFragmentsUpdated{$ref}{$method}{$fragment} } ) {
                foreach my $variant ( keys %{ $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit} } ) {
                    foreach my $blastHit ( keys %{ $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant} } ) {
                        $scoredFragments{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit} = 1;
                    }
                }
            }
        }
    }
}

#print Dumper(\%scoredFragments); exit;

# Criterion 1: Check that the mismatch score for each fragment is not higher than 20%
&checkMismatch();

#print Dumper(\%scoredFragments); exit;

# Criterion 2: Aligned part of the fragment should be at least 90%
&checkLength();

#print Dumper(\%scoredFragments); exit;

# Criterion 3: Forward start site can start before amplicon but not after start amplicon and last 5 bases of 3' end should not contain mismatches
&checkForward();

#print Dumper(\%scoredFragments); exit;

# Criterion 3: Reverse end site can end after amplicon but not before end amplicon and last 5 bases of 3' end at should not contain mismatches
&checkReverse();

#print Dumper(\%scoredFragments); exit;

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Determine whether a method leads to a signal
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# First create a new hash that will contain the scores assuming everything has score 0
# Use the original input hashes to ensure they are all really tested
my %methodResults = ();
foreach my $ref ( keys %refGenomes ) {
    foreach my $method ( keys %methods ) {
        $methodResults{$ref}{$method} = 0;
    }
}

# Determine whether the method leads to a signal
&determineSignal();

#print Dumper(\%methodResults); exit;

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Determine whether a combination of methods leads to at least one signal
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# First create a new hash that will contain the scores assuming everything has score 0
# Use the original input hashes to ensure they are all really tested
my %unionResults = ();
foreach my $ref ( keys %refGenomes ) {
    $unionResults{$ref} = 0;
}

# Determine whether the union of methods leads to a signal
&determineUnion();

#print Dumper(\%unionResults); exit;

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Print to output file(s)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Write output to tab-delimited summary text file that later can be loaded into Excel
&writeSummaryOutput( $params{'-out_summary'} );

# Write output to tab-delimited detailed text file that later can be loaded into Excel
&writeDetailedOutput( $params{'-out_detail'} );

# Write cluster information to tab-delimited detailed text file that can be loaded into Excel
&writeClusterOutput( $params{'-out_cluster'} );

# Bye bye...
exit;

#====================================================================
# SUBROUTINE
#====================================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# generateVariants
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub generateVariants {

    my (@input) = @_;
    if ( $input[0] =~ /(.*)([KDBHRMVYWSN])(.*)/ ) {
        my $head = $1;
        my $tail = $3;
        my @seqs = ();
        foreach my $base ( keys %{ $IUPAC{$2} } ) {
            push( @seqs, generateVariants( $head . $base . $tail ) );
        }
        return (@seqs);
    } else {
        return ( $input[0] );
    }

}    # sub generateVariants

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# loadIUPAC
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub loadIUPAC {

    my %load_IUPAC = (
        'A' => {
            'A' => 1
        },
        'T' => {
            'T' => 1
        },
        'C' => {
            'C' => 1
        },
        'G' => {
            'G' => 1
        },

        q{.} => {
            q{.} => 1,
            q{-} => 1
        },
        q{-} => {
            q{-} => 1,
            q{.} => 1
        },
        'R' => {
            'A' => 1,
            'G' => 1
        },
        'Y' => {
            'C' => 1,
            'T' => 1
        },
        'M' => {
            'C' => 1,
            'A' => 1
        },
        'K' => {
            'T' => 1,
            'G' => 1
        },
        'W' => {
            'T' => 1,
            'A' => 1
        },
        'S' => {
            'C' => 1,
            'G' => 1
        },
        'B' => {
            'C' => 1,
            'T' => 1,
            'G' => 1
        },
        'D' => {
            'A' => 1,
            'T' => 1,
            'G' => 1
        },
        'H' => {
            'A' => 1,
            'T' => 1,
            'C' => 1
        },
        'V' => {
            'A' => 1,
            'C' => 1,
            'G' => 1
        },
        'N' => {
            'A' => 1,
            'C' => 1,
            'G' => 1,
            'T' => 1
        }
    );

    #print Dumper(\%load_IUPAC); exit;

    return (%load_IUPAC);

}    # sub loadIUPAC

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# fasta2hash
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub fasta2hash {

    my ($file) = @_;
    my %id2seq = ();
    my $id     = q{};
    open( my $F, "<", "$file" ) || die "Could not open fasta file $file!:$!\n";
    while (<$F>) {
        chomp;
        if ( $_ =~ /^>(.+)/ ) {
            $id = $1;
        } else {
            my $test = ( ( $_ =~ tr/KDBHRMVYWSN\-\.ATCG//c ) );
            unless ( $test == 0 ) { die "You input in the fasta2hash subroutine contains some non-nucleotide IUPAC characters or U!\n" }
            $id2seq{$id} .= uc($_);
        }
    }
    close $F;
    return (%id2seq);

}    # sub fasta2hash

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# processConfig
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub processConfig {

    my ($file)    = @_;
    my %results   = ();
    my @fragments = qw(-forward -reverse -probe -amplicon);
    open( my $IN, "<", "$file" ) || die "Could not open config file!:$!\n";
    while (<$IN>) {
        chomp( my $row = $_ );
        if ( $row eq q{} || $row =~ /^$/ ) { next }
        my @line = split( m/\t/, $row );
        unless ( defined( $line[0] ) ) { next }    # Added to catch empty lines filled with tabs as sometimes happens

        # First extract the sequence information directly (and make it uppercase)
        my $forward = uc( $line[1] );
        $forward =~ s/\s+//g;
        $results{ $line[0] }{'-forward'} = $forward;
        my $reverse = uc( $line[2] );
        $reverse =~ s/\s+//g;
        my $reverseComp = &reverseComplement($reverse);
        $results{ $line[0] }{'-reverse'} = $reverseComp;
        my $probe = uc( $line[3] );
        $probe =~ s/\s+//g;
        $results{ $line[0] }{'-probe'} = $probe;
        my $amplicon = uc( $line[4] );
        $amplicon =~ s/\s+//g;
        $results{ $line[0] }{'-amplicon'} = $amplicon;

        # Check all the provided sequence information to ensure only nucleotide IUPAC characters (excluding U) are present
        foreach my $i (@fragments) {
            my $seq = $results{ $line[0] }{$i};
            my $test = ( ( $seq =~ tr/KDBHRMVYWSN\-\.ATCG//c ) );
            unless ( $test == 0 ) { die "Your provided input in your config file contains some non-nucleotide IUPAC characters or U!:$seq\n" }
        }

        # For the forwared, reverse, and probe; immediately assign all possible variants accounting for degeneracy genetic code
        @{ $results{ $line[0] }{'-reverseVariants'} } = generateVariants( $results{ $line[0] }{'-reverse'} );
        @{ $results{ $line[0] }{'-probeVariants'} }   = generateVariants( $results{ $line[0] }{'-probe'} );
        @{ $results{ $line[0] }{'-forwardVariants'} } = generateVariants( $results{ $line[0] }{'-forward'} );

    }
    close $IN;
    return (%results);

}    # sub processConfig

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# reverseCompement
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub reverseComplement {

    my ($DNA) = @_;

    # Reverse the DNA sequence
    my $revComp = reverse($DNA);

    # Complement the reversed DNA sequence
    $revComp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    return ($revComp);

}    # sub reverseComplement

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# processMethods
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub processMethods {

    my ($input) = @_;
    my %results = ();
    my @input = split( m/,/, $input );
    foreach (@input) {
        unless ( exists( $info{$_} ) ) { die "You specified a method that is not available in your config file!\n" }
        $results{$_} = 1;
    }
    return (%results);

}    # sub processMethods

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# blastAmplicon
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub blastAmplicon {

    my ($method) = @_;

    # Blast settings
    my $blastString = "7 qseqid qlen sseqid slen qstart qend sstart send qseq sseq evalue bitscore length pident nident mismatch positive gapopen gaps";
    my @blastFields = qw(-qseqid -qlen -sseqid -slen -qstart -qend -sstart -send -qseq -sseq -evalue -bitscore -length -pident -nident -mismatch -positive -gapopen -gaps);

    # Throw warning if amplicon length is too short for conventional blasting
    if ( length( $info{$method}{'-amplicon'} ) <= 30 ) { print STDERR "Your amplicon length is too short for conventional blasting, results may be influenced!" }

    # Make query input fasta file
    my $ampliconFastaFile = "amplicon" . $method . ".fasta";
    open( my $OUT, ">", "$ampliconFastaFile" ) || die "Could not create ouput amplicon fasta file:$!\n";
    print $OUT ">amplicon\n$info{$method}{'-amplicon'}\n";
    close $OUT;

    # Blast amplicon in each reference genome
    foreach my $ref ( keys %refGenomes ) {

        # Make subject fasta file
        my $refFastaFile = "refGenome.fasta";
        open( my $OUT, ">", "$refFastaFile" ) || die "Could not create ouput reference fasta file:$!\n";
        print $OUT ">ref\n$refGenomes{$ref}\n";
        close $OUT;

        # Perform blast, assume one subject but collect all HSPs
        my $blastOut = "ampliconAgainstReferenceBlastOutput";
        system ("blastn -task blastn -query $ampliconFastaFile -outfmt '$blastString' -subject $refFastaFile -out $blastOut -max_target_seqs 1 -strand plus -reward 1 -penalty -1");

        # Collect the output
        my %blastResults = &processBlastOutput( $blastOut, \@blastFields );

        #print Dumper(\%blastResults); exit;
        $foundAmplicons{$ref}{$method}{'-amplicon'} = \%blastResults;

        # Remove extra files
        unlink( $refFastaFile, $blastOut );

    }    # foreach my $ref

    # Delete extra fasta files and return
    unlink($ampliconFastaFile);
    return;

}    # sub blastAmplicon

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# removeGaps
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub removeGaps {

    my $amplicon = "-amplicon";
    my $hit      = "1";
    my $sseq     = "-sseq";
    my $qseq	 = "-qseq";

    foreach my $ref ( keys %foundAmplicons ) {
        foreach my $method ( keys %{ $foundAmplicons{$ref} } ) {
			
			# Remove the gaps in the subject sequence
            my $subjectSequence = $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{$sseq};
            $subjectSequence =~ s/-//g;
            $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{$sseq} = $subjectSequence;         

			# Remove the gaps in the query sequence
            my $querySsequence = $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{$qseq};
            $querySsequence =~ s/-//g;
            $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{$qseq} = $querySsequence; 
        }
    }

    return;

}    # sub removeGaps

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# extendAmplicons
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub extendAmplicons {

    my $amplicon     = "-amplicon";
    my $hit          = "1";
    my $length       = "-length";
    my $maxExtension = $params{'-ext_am'};
    my $sseq         = "-sseq";
	
    foreach my $ref ( keys %foundAmplicons ) {
        foreach my $method ( keys %{ $foundAmplicons{$ref} } ) {

            # If the amplicon is FALSE, do not do any extension
            if ( $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-qseq'} eq "FALSE" ) { next }

            # Needs to be set here or results may be carried over the loop to the next reference genome!
            my ( $forwardExtension, $reverseExtension ) = ( q{}, q{} );

            # Capture if forward extension is required and allowed, and perform if so
            if ( $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-qstart'} != 1 ) {

                # Capture with how many bases the amplicon should be extended, correcting for the fact it can never go further than the genome boundary
                my $extension;
                my $requiredExtension = ( $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-qstart'} - 1 );
                my $possibleExtension = ( $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-sstart'} - 1 );
                if ( $requiredExtension < $possibleExtension ) { $extension = $requiredExtension }
                else { $extension = $possibleExtension }
                
                # Expression for start site has been tested and verified (discrepancy perl and blast in numbering!)
                my $startSite = ( $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-sstart'} - $extension - 1 );
                if ( ( $startSite < 0 ) || ( $extension < 0 ) ) { die "startSite or $extension can never be lower than 0 in extendAmplicons subroutine!\n" }; # Legacy check just to be sure
                $forwardExtension = substr( $refGenomes{$ref}, $startSite, $extension );

                # Start site in query and subject also need to be updated for later criteria checks!
                $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-qstart'} = ( $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-qstart'} - $extension );
            	$foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-sstart'} = ( $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-sstart'} - $extension );

            }

            # Capure if reverse extension is required and allowed, and perform if so, correcting for the fact it can never go further than the genome boundary
            if ( $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-qend'} < $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-qlen'} ) {
                
                # Capture with how many bases the amplicon should be extended
                my $extension;
                my $requiredExtension = ( $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-qlen'} - $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-qend'} );
                my $possibleExtension = ( length( $refGenomes{$ref} ) - $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-send'} );
                if ( $requiredExtension < $possibleExtension ) { $extension = $requiredExtension }
                else { $extension = $possibleExtension }
                
                # Expression for start site has been tested and verified (discrepancy perl and blast in numbering!)
                my $startSite = ( $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-send'} );
                if ( ( $startSite < 0 ) || ( $extension < 0 ) ) { die "startSite or $extension can never be lower than 0 in extendAmplicons subroutine!\n" }; # Legacy check just to be sure
                $reverseExtension = substr( $refGenomes{$ref}, $startSite, $extension );

                # End site in the query and subject also need to be updated for later criteria checks
                $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-qend'} = ( $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-qend'} + $extension );
            	$foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-send'} = ( $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-send'} + $extension );
            
            }

            # Update amplicon with either extension if existing
            $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-sseq'} = $forwardExtension . $foundAmplicons{$ref}{$method}{$amplicon}{$hit}{'-sseq'} . $reverseExtension;

        }
    }

    return;

}    # sub extendAmplicons

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# clusterAmpliconGreedy
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub clusterAmpliconGreedy {

    my ($method) = @_;

    # First make a fasta file dump of all found amplicons for the particular method in question
    my $fasta = "FoundAmplicons" . $method;
    open( my $OUT, ">", "$fasta" ) || die "Could not open output file in &clusterAmpliconGreedy subroutine!:$!\n";
    foreach my $ref ( keys %foundAmplicons ) {
        my $sequence = $foundAmplicons{$ref}{$method}{'-amplicon'}{'1'}{'-sseq'};
        print $OUT ">$ref\n$sequence\n";
    }
    close $OUT;

    # Use usearch to do the clustering
    my $clust = "FoundCLusters" . $method;
    system ("usearch -cluster_fast $fasta -sort length -id 1.0 -uc $clust -notrunclabels -quiet -query_cov 1.0 -target_cov 1.0");

    # Collect output and keep track of highest cluster number
    my $maxCluster = 0;
    open( my $IN, "<", "$clust" ) || die "Could not open input file in &clusterAmpliconGreedy subroutine!:$!\n";
    while (<$IN>) {
        chomp( my $row = $_ );
        if ( $row =~ /^#/ ) { next }    # Usearch comments
        my @line = split( m/\t/, $row );
        $foundClusters{$method}{ $line[1] }{ $line[8] } = 1;
        if ( $line[1] > $maxCluster ) { $maxCluster = $line[1] }    # Keep track of highest cluster number
    }
    close $IN;
    $maxClusters{$method} = $maxCluster;

    # Delete extra files and return
    unlink( $fasta, $clust );
    return;

}    # sub clusterAmpliconGreedy

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# clusterAmpliconExhaustive
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub clusterAmpliconExhaustive {

    my ($method) = @_;
    my %network = ();

    # Perform an all versus all sequence comparison (calculations are done only for one half, other half assigned through symmetry of distance matrix)
    my @references = ( keys %foundAmplicons );
    my $size2      = scalar(@references);
    my $size       = $size2 - 1;

    for ( my $i = 0 ; $i < $size ; $i++ ) {
        for ( my $j = ( $i + 1 ) ; $j < $size2 ; $j++ ) {

            # Initialize
            my $refOne        = $references[$i];
            my $sequenceOne   = $foundAmplicons{$refOne}{$method}{'-amplicon'}{'1'}{'-sseq'};
            my $refTwo        = $references[$j];
            my $sequenceOther = $foundAmplicons{$refTwo}{$method}{'-amplicon'}{'1'}{'-sseq'};

            # Ignore the diagonals of the matrix
            if ( $refOne eq $refTwo ) { next }

            # Ignore amplicons with a "FALSE" sequence (so they are not aligned leading to strange errors)
            if ( $sequenceOne eq "FALSE" || $sequenceOther eq "FALSE" ) { next }

            # Prepare muscle input file
            my $muscleInput = "amplicon" . $method . "MuscleInput.fasta";
            open( my $OUT, ">", "$muscleInput" ) || die "Could not create input fasta file in clusterAmpliconExhaustive subroutine:$!\n";
            print $OUT ">$refOne\n$sequenceOne\n>$refTwo\n$sequenceOther\n";
            close $OUT;

            # Perform a global alignment of both sequences
            my $muscleOutput = "amplicon" . $method . "MuscleOutput.fasta";
            system ("muscle -in $muscleInput -out $muscleOutput -quiet");

            # Read in muscle output and delete it
            my %muscle = &fasta2hash($muscleOutput);
            unlink( $muscleInput, $muscleOutput );

            #print Dumper(\%muscle);

            # Calculate identify score based on global alignment while accounting for degeneracy of IUPAC
            # Computationally not efficient: would there be a better way?
            my $seqOne   = $muscle{$refOne};
            my $seqTwo   = $muscle{$refTwo};
            my $mismatch = 0;
            for ( my $i = 0 ; $i < length($seqOne) ; $i++ ) {    # Global alignment so length should be the same for both seq
                my $baseOne = substr( $seqOne, $i, 1 );
                my $baseTwo = substr( $seqTwo, $i, 1 );
                my @baseOneHits = ( keys %{ $IUPAC{$baseOne} } );
                my @baseTwoHits = ( keys %{ $IUPAC{$baseTwo} } );
                my $hitFlag     = 0;
                foreach my $one (@baseOneHits) {
                    foreach my $two (@baseTwoHits) {
                        if ( $one eq $two ) { $hitFlag = 1; }
                    }
                }
                if    ( $hitFlag == 1 ) { next }
                elsif ( $hitFlag == 0 ) { $mismatch++ }
            }

            # Only connect two sequences if they have 100% identity!
            if ( $mismatch == 0 ) {
                $network{$refTwo}{$refOne} = 1;
                $network{$refOne}{$refTwo} = 1;
            }

        }
    }

    #print Dumper(\%network); exit;

    # Do a recursive search where only connected subcomponens with an identity score equal to 100% are extracted
    my $count   = 0;
    my %visited = ();
    foreach my $ref ( keys %network ) {
        if ( !exists( $visited{$ref} ) ) {
            &explore( $method, $ref, $count, \%visited, \%network );
            $count++;
        }

    }

    # Add all 'singletons', sequences that were not taken up into the network and that hence have not been considered yet
    foreach my $ref ( keys %foundAmplicons ) {
        if ( !exists( $visited{$ref} ) ) {
            $foundClusters{$method}{$count}{$ref} = 1;
            $count++;
        }
    }

    $maxClusters{$method} = --$count;

    return;

}    # sub clusterAmpliconExhaustive

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub explore
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub explore {

    my ( $method, $ref, $count, $refVisited, $refNetwork ) = @_;
    ${$refVisited}{$ref} = 1;
    $foundClusters{$method}{$count}{$ref} = 1;
    foreach my $other ( keys %{ ${$refNetwork}{$ref} } ) {
        if ( !exists( ${$refVisited}{$other} ) ) {
            &explore( $method, $other, $count, $refVisited, $refNetwork );
        }
    }

    return;

}    # sub explore

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# blastFragment
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# WARNING: THIS SUBROUTINE ASSUMES IT SHOULD ONLY BLAST FRAGMENTS AGAINST THE BEST HIT AMPLICON!
sub blastFragment {

    my ( $method, $fragment ) = @_;

    # Blast settings, make sure only to add extra parameters but never omit them as some of them are used throughout the reminder of this script for checking (i.e. -evalue in this subroutine)
    my $blastString = "7 qseqid qlen sseqid slen qstart qend sstart send qseq sseq evalue bitscore length pident nident mismatch positive gapopen gaps";
    my @blastFields = qw(-qseqid -qlen -sseqid -slen -qstart -qend -sstart -send -qseq -sseq -evalue -bitscore -length -pident -nident -mismatch -positive -gapopen -gaps);

    # Loop over the different variants of the fragment, if requested
    my @fragmentVariants = ();
    if ( $params{'-deg'} eq "no" ) { push( @fragmentVariants, $info{$method}{ q{-} . $fragment } ) }
    elsif ( $params{'-deg'} eq "yes" ) { @fragmentVariants = @{ $info{$method}{ q{-} . $fragment . "Variants" } } }
    foreach my $variant (@fragmentVariants) {

        # Make query input fasta file for the variant of the fragment
        my $fragmentFastaFile = $fragment . $method . ".fasta";
        open( my $OUT, ">", "$fragmentFastaFile" ) || die "Could not create ouput fragment fasta file:$!\n";
        print $OUT ">$fragment\n$variant\n";
        close $OUT;

        # Loop over each reference genome (amplicon)
        foreach my $ref ( keys %refGenomes ) {

            # First test if an amplicon was found, if not finding the fragment is impossible, so return FALSE
            if ( $foundAmplicons{$ref}{$method}{'-amplicon'}{"1"}{'-sseq'} eq "FALSE" ) {
                my %falseResults = ();
                for ( my $i = 0 ; $i < scalar(@blastFields) ; $i++ ) {
                    $falseResults{"1"}{ $blastFields[$i] } = "FALSE";    # This mimicks the &processBlast subroutine output
                    $falseResults{"1"}{'-bitscore'} = "0";               # This special case needs to be set because bitscore is used in some subroutines for sorting
                }
                $foundFragments{$ref}{$method}{ q{-} . $fragment }{$variant} = \%falseResults;
                next;
            }

            # Make subject fasta file
            # Here a loop could be introduced to blast against all best amplicons but currently it assumes only top best amplicon should be searhced
            my $refAmpliconFile = "refAmplicon.fasta";
            open( my $OUT, ">", "$refAmpliconFile" ) || die "Could not create ouput reference fasta file:$!\n";
            my $sequence = $foundAmplicons{$ref}{$method}{'-amplicon'}{"1"}{'-sseq'};
            print $OUT ">ref\n$sequence\n";
            close $OUT;

            # Throw warning if fragment length is too long
            # Warning currently disabled because it occurs quite often but does not seem to influence the quality as long as the fragment is not too long
            #if (length($info{$method}{"-".$fragment}) >= 30) {print STDERR "Your fragment length is too long for 'short blasting', results may be influenced!\n"}

            # Perform the blast run, assume one subject but collect all HSPs
            # Blast parameters have been tuned with low penalty and reward scores to account for the short fragment blasting
            my $blastOut = "fragmentAgainstAmpliconBlastOutput";
            system ("blastn -task blastn-short -query $fragmentFastaFile -outfmt '$blastString' -subject $refAmpliconFile -out $blastOut -max_target_seqs 1 -strand plus -penalty -1 -reward 1 -word_size 4");

            # Collect the output
            my %blastResults = &processBlastOutput( $blastOut, \@blastFields );

            # DISABLED: blastn-short is now run anyway with the shortest allowed word length (4)
            # If no output is run, rerun with a lower word size
            #			if ($blastResults{"1"}{"-evalue"} eq "FALSE") { # Checking evalue is a good bet because this value is as good as always collected, might show strange behaviour if for some reasing the evalue is ever ommitted!
            #				my $flag = 1;
            #				my $wordSize = 6; # 7 is the standard word size used in blastn-short so start with 6 here
            #				while ($flag) {
            #					system( "blastn -task blastn-short -query $fragmentFastaFile -outfmt '$blastString' -subject $refAmpliconFile -out $blastOut -max_target_seqs 1 -strand plus -penalty -1 -reward 1 -word_size $wordSize")";
            #					%blastResults = &processBlastOutput($blastOut,\@blastFields);
            #					# Keep rerunning and lowering wordsize until you get results
            #					if ($blastResults{"1"}{"-evalue"} eq "FALSE") {
            #						$wordSize--;
            #					} else {$flag = 0;}
            #					# If wordsize reaches however 3, quit because blastn refuses to run lower or equal to this value
            #					if ($wordSize == 3) {$flag = 0}
            #				}
            #			}
            #			#print Dumper(\%blastResults); exit;
            $foundFragments{$ref}{$method}{ q{-} . $fragment }{$variant} = \%blastResults;

            # Delete extra fasta files
            unlink( $blastOut, $refAmpliconFile );

        }    # foreach my $ref

        # Delete extra fasta file
        unlink($fragmentFastaFile);

    }    #foreach my $variant

    return;

}    # sub blastFragment

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sortVariantScores
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub sortVariantScores {

    # Loop over the hash with the found blasted fragments
    foreach my $ref ( keys %foundFragments ) {
        foreach my $method ( keys %{ $foundFragments{$ref} } ) {
            foreach my $fragment ( keys %{ $foundFragments{$ref}{$method} } ) {

                # Only consider the hash that carries all fragments
                my %HoH = %{ $foundFragments{$ref}{$method}{$fragment} };

                #print "I collected the following information for ref $ref and method $method and fragment $fragment:\n";
                #print Dumper(\%HoH);
                #print "Now I will attempt to sort it!\n";

                # Now sort them on the value of the bitscores at the leaves
                my @sorted = sort { $HoH{$b}{'1'}{'-bitscore'} <=> $HoH{$a}{'1'}{'-bitscore'} } keys %HoH;

                #print Dumper(\@sorted);

                # DEPRECATED: Store ALL results
                #				foreach (my $i=0;$i<scalar(@sorted);$i++) {
                #					my $variant = $sorted[$i];
                #					my $score = $foundFragments{$ref}{$method}{$fragment}{$variant}{'1'}{'-bitscore'};
                #					# Uncomment below to check whether code really worked
                #					# print "Variant $i is $variant with score $score\n";
                #					# Store in the main hash
                #					my $counter = ($i + 1);
                #					%{$sortedFragments{$ref}{$method}{$fragment}{$counter}{$variant}{'1'}} = %{$foundFragments{$ref}{$method}{$fragment}{$variant}{'1'}};
                #				}

                # Only store the best result
                %{ $sortedFragments{$ref}{$method}{$fragment}{'1'}{ $sorted[0] }{'1'} } = %{ $foundFragments{$ref}{$method}{$fragment}{ $sorted[0] }{'1'} };

            }    # foreach my $fragment
        }    # foreach my $method
    }    # foreach my $ref

    return;

}    # sub sortVariantScores

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# extendFragments
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub extendFragments {

    foreach my $ref ( keys %sortedFragments ) {
        foreach my $method ( keys %{ $sortedFragments{$ref} } ) {
            foreach my $fragment ( keys %{ $sortedFragments{$ref}{$method} } ) {
                foreach my $fragmentHit ( keys %{ $sortedFragments{$ref}{$method}{$fragment} } ) {
                    foreach my $variant ( keys %{ $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit} } ) {
                        foreach my $variantHit ( keys %{ $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant} } ) {
												
                            # First catch blast hits that did not have a match ("FALSE")
                            if ( $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-sseq'} eq "FALSE" ) {
                            	 next
                            }
                                                                                   
           					# Needs to be set here or results may be carried over the loop to the next reference genome!
            				my ( $forwardExtensionSubject, $reverseExtensionSubject, $forwardExtensionQuery, $reverseExtensionQuery ) = ( q{}, q{}, q{}, q{} );
                            
                            # For both forward and reverse primer, and amplicon,  extend the LEFT END (5')
         					# Capture if 5' extension is required and allowed, and perform if so
                            if ( $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-qstart'} != 1 ) {
                                                                               
              					# Capture with how many bases the fragment should be extended, correcting for the fact it can never go further than the amplicon boundary
	                            my $extension;
	                            my $requiredExtension = ( $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-qstart'} - 1 );
	                            my $possibleExtension = ( $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-sstart'} - 1 );
                				if ( $requiredExtension < $possibleExtension ) { $extension = $requiredExtension }
                				else { $extension = $possibleExtension }
	
	                            # Update the sequence of the subject (expression for start site has been tested and verified (discrepancy perl and blast in numbering))
	                            my $startSiteSubject = ( $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-sstart'} - $extension - 1 );
	                            if ( ( $startSiteSubject < 0 ) || ( $extension < 0 ) ) { die "startSiteSubject or extension can never be lower than 0 in extendFragments subroutine!\n" }; # Legacy check
	                            $forwardExtensionSubject = substr( $foundAmplicons{$ref}{$method}{'-amplicon'}{'1'}{'-sseq'}, $startSiteSubject, $extension );
	
		                        # Update the sequence of the query (expression for start site has been tested and verified (discrepancy perl and blast in numbering))
	                            my $startSiteQuery = ( $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-qstart'} - $extension - 1 );
	                            if ( ( $startSiteQuery < 0 ) || ( $extension < 0 ) ) { die "startSiteQuery or extension can never be lower than 0 in extendFragments subroutine!\n" }; # Legacy check
	                            $forwardExtensionQuery = substr( $variant, $startSiteQuery, $extension );
	                            
	                            # Update critical blast stats that are used later in the criteria checks
	                            $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-qstart'} = ( $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-qstart'} - $extension );
	                            $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-sstart'} = ( $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-sstart'} - $extension );
	                              
                            }

                            # For both forward and reverse primer, and amplicon, extend the RIGHT END (3')
         					# Capture if 3' extension is required and allowed, and perform if so
                            if ( $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-qend'} < $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-qlen'} ) {
                            	                           	
              					# Capture with how many bases the fragment should be extended, correcting for the fact it can never go further than the amplicon boundary
                           		my $extension;
                            	my $requiredExtension = ( $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-qlen'} - $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-qend'} );
                            	my $possibleExtension = ( length( $foundAmplicons{$ref}{$method}{'-amplicon'}{'1'}{'-sseq'} ) - $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-send'} );
               					if ( $requiredExtension < $possibleExtension ) { $extension = $requiredExtension }
                				else { $extension = $possibleExtension }

	                            # Update the sequence of the subject (expression for start site has been tested and verified (discrepancy perl and blast in numbering))
	                            my $startSiteSubject = ( $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-send'} );
	                            if ( ( $startSiteSubject < 0 ) || ( $extension < 0 ) ) { die "startSiteSubject or extension can never be lower than 0 in extendFragments subroutine!\n" }; # Legacy check
	                            $reverseExtensionSubject = substr( $foundAmplicons{$ref}{$method}{'-amplicon'}{'1'}{'-sseq'}, $startSiteSubject, $extension );

	                            # Update the sequence of the query (expression for start site has been tested and verified (discrepancy perl and blast in numbering))
	                            my $startSiteQuery = ( $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-qend'} );
	                            if ( ( $startSiteQuery < 0 ) || ( $extension < 0 ) ) { die "startSiteQuery or extension can never be lower than 0 in extendFragments subroutine!\n" }; # Legacy check
	                            $reverseExtensionQuery = substr( $variant, $startSiteQuery, $extension );
	                            	                            
	                            # Update critical blast stats
	                            $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-qend'} = $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-qend'} + $extension;
	                            $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-send'} = $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-send'} + $extension;

                            }
 						
	 					    # Update fragment with either extension if existing
	           				$sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-sseq'} = $forwardExtensionSubject . $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-sseq'} . $reverseExtensionSubject;
	           				$sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-qseq'} = $forwardExtensionQuery . $sortedFragments{$ref}{$method}{$fragment}{$fragmentHit}{$variant}{$variantHit}{'-qseq'} . $reverseExtensionQuery;
	           						                     
                        }                         
                    }
                }
            }
        }
    }

    return;

}    # sub extendFragments

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# processBlastOutput
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub processBlastOutput {

    # Initialize
    my ( $input, $blastFieldsRef ) = @_;
    my @output  = ();
    my %results = ();

    # Collect results
    open( my $IN, "<", "$input" ) || die "Could not open blast output file $input:$!\n";
    while (<$IN>) {
        chomp( my $row = $_ );

        # Catch comments
        if ( $row =~ m/^#/ ) {

            # Catch no hits and return directly, mind the scope of %noHits!
            if ( $row =~ m/^# 0 hits found/ ) {
                my %noHits = ();
                for ( my $i = 0 ; $i < scalar(@{$blastFieldsRef}) ; $i++ ) {
                    $noHits{"1"}{ ${$blastFieldsRef}[$i] } = "FALSE";
                }
                $noHits{"1"}{'-bitscore'} = "0";    # This special case needs to be set because bitscore is used in some subroutines for sorting
                return (%noHits);
            } else {
                next;
            }
        }

        # Collect actual data if a hit is found, mind the scope of %hits
        my @line = split( m/\t/, $row );
        my %hits = ();
        for ( my $i = 0 ; $i < scalar(@line) ; $i++ ) {
            $hits{ ${$blastFieldsRef}[$i] } = $line[$i];
        }
        push( @output, \%hits );
    }
    close $IN;

    #print Dumper(\@output); exit;

    # Sort the blast output based on its bit score
    my @sortedOutput = sort { $b->{-bitscore} <=> $a->{-bitscore} } @output;

    #print Dumper(\@sortedOutput);
    my $counter = 1;
    for ( my $i = 0 ; $i < scalar(@sortedOutput) ; $i++ ) {
        my %hash = %{ $sortedOutput[$i] };
        $results{$counter} = \%hash;
        $counter++;
        last;    #DANGEROUS: BRUTE FORCE HACK TO ONLY ALLOW THE BEST RESULT TO BE PRODUCED
    }

    #print Dumper(\%results);

    return (%results);

}    # sub processBlastOutput

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# recalculateMismatch
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub recalculateMismatch {

    # Loop over all found fragments and recalculate the number of mismatches
    foreach my $ref ( keys %foundFragmentsUpdated ) {
        foreach my $method ( keys %{ $foundFragmentsUpdated{$ref} } ) {
            foreach my $fragment ( keys %{ $foundFragmentsUpdated{$ref}{$method} } ) {
                foreach my $variantHit ( keys %{ $foundFragmentsUpdated{$ref}{$method}{$fragment} } ) {
                    foreach my $variant ( keys %{ $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit} } ) {
                        foreach my $blastHit ( keys %{ $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant} } ) {

                            # Extract the sequence of the query and hit aligned part first to make life a bit easier
                            my $sseq = $sortedFragments{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-sseq'};
                            my $qseq = $sortedFragments{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-qseq'};

                            # Both should have the same length
                            unless ( length( $sseq ) == length( $qseq ) ) {
                                print "Error for reference $ref method $method fragment $fragment!\n";
                                print "Your query seq is $qseq\n";
                                print "Your subject seq is $sseq\n";

                                die "In recalculateMismatch subroutine the length of the query and subject sequence should always be the same!\n";
                            }

                            # Filter out "FALSE" (fragments for which no blast hit was found in the amplicon in the first place)
                            if ( $sseq eq "FALSE" ) {
                                $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-mismatch'} = "FALSE";
                                next;
                            }

                            # Loop over them and compare. There is probably a more elegant and computationally less demadning way to do this but tired + head-ache + deadline...
                            my $mismatch = 0;
                            for ( my $i = 0 ; $i < length($sseq) ; $i++ ) {
                                my $sbase = substr( $sseq, $i, 1 );
                                my $qbase = substr( $qseq, $i, 1 );
                                my @sbaseHits = ( keys %{ $IUPAC{$sbase} } );
                                my @qbaseHits = ( keys %{ $IUPAC{$qbase} } );
                                my $hitFlag   = 0;
                                foreach my $q (@qbaseHits) {
                                    foreach my $s (@sbaseHits) {
                                        if ( $q eq $s ) { $hitFlag = 1 }
                                    }
                                }
                                if    ( $hitFlag == 1 ) { next }
                                elsif ( $hitFlag == 0 ) { $mismatch++ }
                            }

                            # Assign novel mismatchscore
                            $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-mismatch'} = $mismatch;

                        }
                    }
                }
            }
        }
    }

    return;

}    # recalculateMismatch

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# checkMismatch
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub checkMismatch {
	
	my $criterion = ( $params{'-mismatch'} / 100 );
	
    foreach my $ref ( keys %foundFragmentsUpdated ) {
        foreach my $method ( keys %{ $foundFragmentsUpdated{$ref} } ) {
            foreach my $fragment ( keys %{ $foundFragmentsUpdated{$ref}{$method} } ) {
                foreach my $variantHit ( keys %{ $foundFragmentsUpdated{$ref}{$method}{$fragment} } ) {
                    foreach my $variant ( keys %{ $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit} } ) {
                        foreach my $blastHit ( keys %{ $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant} } ) {

                            # Define first for simplicity
                            my $mismatch        = $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-mismatch'};
                            my $alignmentLength = length( $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-sseq'} );

                            # Filter out "FALSE" (fragments for which no blast hit was found in the amplicon in the first place)
                            if ( $mismatch eq "FALSE" ) {
                                $scoredFragments{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit} = 0;
                                next;
                            }

                            # Test if mismatch percentage is higher than a specified maximum, if so assign an invalid score
                            if ( ( $mismatch / $alignmentLength ) > $criterion ) {
                                $scoredFragments{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit} = -1;
                            }

                        }
                    }
                }
            }
        }
    }

    return;

}    # sub checkMismatch

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# checkLength
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub checkLength {

	my $criterion = ( $params{'-length'} / 100 );
	
    foreach my $ref ( keys %foundFragmentsUpdated ) {
        foreach my $method ( keys %{ $foundFragmentsUpdated{$ref} } ) {
            foreach my $fragment ( keys %{ $foundFragmentsUpdated{$ref}{$method} } ) {
                foreach my $variantHit ( keys %{ $foundFragmentsUpdated{$ref}{$method}{$fragment} } ) {
                    foreach my $variant ( keys %{ $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit} } ) {
                        foreach my $blastHit ( keys %{ $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant} } ) {

                            # Filter out "FALSE" (fragments for which no blast hit was found in the amplicon in the first place)
                            if ( $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-mismatch'} eq "FALSE" ) {
                                $scoredFragments{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit} = 0;
                                next;
                            }
                            
                            # Define first for simplicity
                            my $alignmentLength = length ( $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-sseq'} );
                            my $fragmentLength  = $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-qlen'};

                            # Test if aligned part of fragment is lower than the specified minimum, if so assign an invalid score
                            if ( ( $alignmentLength / $fragmentLength ) < $criterion ) {
                                $scoredFragments{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit} = -2;
                            }

                        }
                    }
                }
            }
        }
    }

    return;

}    # sub checkLength

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# checkForward
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub checkForward {

    my $regionLength = $params{'-end'} * -1;
    my $criterion = $params{'-bases'};
    
    foreach my $ref ( keys %foundFragmentsUpdated ) {
        foreach my $method ( keys %{ $foundFragmentsUpdated{$ref} } ) {
            my $fragment = "-forward";
            foreach my $variantHit ( keys %{ $foundFragmentsUpdated{$ref}{$method}{$fragment} } ) {
                foreach my $variant ( keys %{ $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit} } ) {
                    foreach my $blastHit ( keys %{ $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant} } ) {

                        # Filter out "FALSE" (fragments for which no blast hit was found in the amplicon in the first place)
                        if ( $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-mismatch'} eq "FALSE" ) {
                            $scoredFragments{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit} = 0;
                            next;
                        }

                        # Do a simple test to see if the region length defined by the user is not larger than the found annealing site in the first place
                        my $totalLength = length( $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-sseq'} );
                        if ( $totalLength < abs( $regionLength) ) { $scoredFragments{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit} = -3 }
                                             
                        # Extract a number of bases at the 3' end as defined by the user
                        my $querySubpiece   = substr( $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-qseq'}, $regionLength );
                        my $subjectSubpiece = substr( $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-sseq'}, $regionLength );

                        # Check for mismatches accounting for degeneracy of IUPAC code
                        my $mismatch = 0;
                        for ( my $i = 0 ; $i < length($querySubpiece) ; $i++ ) {
                            my $sbase = substr( $subjectSubpiece, $i, 1 );
                            my $qbase = substr( $querySubpiece,   $i, 1 );
                            my @sbaseHits = ( keys %{ $IUPAC{$sbase} } );
                            my @qbaseHits = ( keys %{ $IUPAC{$qbase} } );
                            my $hitFlag   = 0;
                            foreach my $q (@qbaseHits) {
                                foreach my $s (@sbaseHits) {
                                    if ( $q eq $s ) { $hitFlag = 1 }
                                }
                            }
                            if    ( $hitFlag == 1 ) { next }
                            elsif ( $hitFlag == 0 ) { $mismatch++ }
                        }
                        if ( $mismatch == 0 ) { next }
                        elsif ( $mismatch >= $criterion ) {
                            $scoredFragments{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit} = -3;
                        }

                    }
                }
            }
        }
    }

    return;

}    # sub checkForward

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# checkReverse
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub checkReverse {

    my $regionLength = $params{'-end'};
    my $criterion = $params{'-bases'};
    
    foreach my $ref ( keys %foundFragmentsUpdated ) {
        foreach my $method ( keys %{ $foundFragmentsUpdated{$ref} } ) {
            my $fragment = "-reverse";
            foreach my $variantHit ( keys %{ $foundFragmentsUpdated{$ref}{$method}{$fragment} } ) {
                foreach my $variant ( keys %{ $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit} } ) {
                    foreach my $blastHit ( keys %{ $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant} } ) {

                        # Filter out "FALSE" (fragments for which no blast hit was found in the amplicon in the first place)
                        if ( $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-mismatch'} eq "FALSE" ) {
                            $scoredFragments{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit} = 0;
                            next;
                        }
                       
                        # Do a simple test to see if the region length defined by the user is not larger than the found annealing site in the first place
                        my $totalLength = length( $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-sseq'} );
                        if ( $totalLength < abs( $regionLength) ) { $scoredFragments{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit} = -3 }
                        
                        # Extract a number of bases at the 3' end as defined by the user but remember that fragments are reverse complemented so to get 3' you need to take the first 5 bases!!!
                        my $querySubpiece   = substr( $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-qseq'}, 0, $regionLength );
                        my $subjectSubpiece = substr( $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-sseq'}, 0, $regionLength );

                        # Check for mismatches accounting for degeneracy of IUPAC code
                        my $mismatch = 0;
                        for ( my $i = 0 ; $i < length($querySubpiece) ; $i++ ) {
                            my $sbase = substr( $subjectSubpiece, $i, 1 );
                            my $qbase = substr( $querySubpiece,   $i, 1 );
                            my @sbaseHits = ( keys %{ $IUPAC{$sbase} } );
                            my @qbaseHits = ( keys %{ $IUPAC{$qbase} } );
                            my $hitFlag   = 0;
                            foreach my $q (@qbaseHits) {
                                foreach my $s (@sbaseHits) {
                                    if ( $q eq $s ) { $hitFlag = 1 }
                                }
                            }
                            if    ( $hitFlag == 1 ) { next }
                            elsif ( $hitFlag == 0 ) { $mismatch++ }
                        }
                        if ( $mismatch == 0 ) { next }
                        elsif ( $mismatch >= $criterion ) {
                            $scoredFragments{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit} = -3;
                        }

                    }
                }
            }
        }
    }

    return;

}    # sub checkReverse

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# determineSignal
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub determineSignal {

    my @fragments  = qw(-forward -reverse -probe);
    my $variantHit = 1;                              # Only the top hit is considered and has been propagated
    my $blastHit   = 1;                              # Only the top hit is considered and has been propagated

    foreach my $ref ( keys %refGenomes ) {
        foreach my $method ( keys %methods ) {

            # Each fragment associated to the method needs to be positive for annealing before the method itself can be considered positive and lead to a signal!
            my $hitCount = 0;
            for ( my $i = 0 ; $i < scalar(@fragments) ; $i++ ) {
                my @variants = ( keys %{ $scoredFragments{$ref}{$method}{ $fragments[$i] }{$variantHit} } );
                if ( $scoredFragments{$ref}{$method}{ $fragments[$i] }{$variantHit}{ $variants[0] }{$blastHit} == 1 ) { $hitCount++ }
            }
            if ( $hitCount == scalar(@fragments) ) { $methodResults{$ref}{$method} = 1 }

        }

    }

    return;

}    # sub determineSignal

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# determineUnion
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub determineUnion {

    foreach my $ref ( keys %refGenomes ) {
        foreach my $method ( keys %methods ) {

            # Only one method needs to succeed to be considered a success
            if ( $methodResults{$ref}{$method} == 1 ) { $unionResults{$ref} = 1 }
        }
    }

    return;

}    # determineUnion

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# writeSummaryOutput
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub writeSummaryOutput {

    my ($input) = @_;

	open( my $OUT_SUMMARY, ">", "$input" ) || die "Could not open summary output file for final results:$!\n";
	
	# Print run information
	&writeInfo($OUT_SUMMARY);

	# Print the actual data	
	print $OUT_SUMMARY "# Reference genome\tPositive signal for combination of methods (o = no, 1 = yes): $params{'-method'}\n";
	foreach my $ref ( keys %refGenomes ) {
	    print $OUT_SUMMARY "$ref\t$unionResults{$ref}\n";
	}
	close $OUT_SUMMARY;
	return;
	
}	# writeSummaryOutput

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# writeDetailedOutput
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub writeDetailedOutput {

    my ($input) = @_;

    my @fragments  = qw(-forward -reverse -probe);
    my $variantHit = 1;                              # Only the top hit is considered and has been propagated
    my $blastHit   = 1;                              # Only the top hit is considered and has been propagated

    open( my $OUT_DETAIL, ">", "$input" ) || die "Could not open detailed output file for final results:$!\n";

	# Print run information
	&writeInfo($OUT_DETAIL);

	# Print data
    print $OUT_DETAIL "# Reference genome tested\tMethod tested\tAmplicon code (see manual!)\tFound amplicon sequence in reference genome\t";
    print $OUT_DETAIL "Forward primer anneals? (1 = yes, < 1 = no (see manual!))\tExact (best) sequence of forward primer used for blasting\tAnnealing sequence of forward primer in amplicon\tAnnealing sequence of forward primer in fragment\tNumber of mismatches in forward primer annealing site\tLength of forward primer annealing site\t";
    print $OUT_DETAIL "Reverse primer anneals? (1 = yes, < 1 = no (see manual!))\tExact (best) sequence of reverse primer used for blasting (reverse complemented!)\tAnnealing sequence of reverse primer in amplicon (reverse complemented!)\tAnnealing sequence of reverse primer in fragment (reverse complemented!)\tNumber of mismatches in reverse primer annealing site\tLength of reverse primer annealing site\t";
    print $OUT_DETAIL "Probe anneals? (1 = yes, < 1 = no (see manual!))\tExact (best) sequence of probe used for blasting\tAnnealing sequence of probe in amplicon\tAnnealing sequence of probe in fragment\tNumber of mismatches in probe annealing site\tLength of probe annealing site\n";

    foreach my $ref ( keys %refGenomes ) {
        foreach my $method ( keys %methods ) {
			
			# Print method information
            print $OUT_DETAIL "$ref\t$method";
			
			# Check amplicon start and end cases
            my ( $ampliconGenomeStart, $ampliconGenomeEnd ) = ( 0, 0 );
            my $foundAmplicon = $foundAmplicons{$ref}{$method}{'-amplicon'}{'1'}{'-sseq'};
            unless ( $foundAmplicon eq "FALSE") {
            	if ( $foundAmplicons{$ref}{$method}{'-amplicon'}{'1'}{'-sstart'} == 1 ) { $ampliconGenomeStart = 1 } 
            	if ( $foundAmplicons{$ref}{$method}{'-amplicon'}{'1'}{'-send'} == length( $refGenomes{$ref} ) ) {$ampliconGenomeEnd = 1 }
            }
            
    		# Print amplicon information
    		my $ampliconFlag = qw{};
    		if ($foundAmplicon eq "FALSE") { $ampliconFlag = 0 }
    		elsif ( $ampliconGenomeStart == 1 ) {$ampliconFlag = -1}        
    		elsif ( $ampliconGenomeEnd == 1 ) {$ampliconFlag = -2}        
    		else { $ampliconFlag = 1 }
    		if ( $ampliconGenomeStart == 1 && $ampliconGenomeEnd == 1 ) { die "Your amplicon cannot be both at the beginning and end of the genome!\n" }        
            print $OUT_DETAIL "\t$ampliconFlag\t$foundAmplicon";
			
			# Print fragment information
            foreach ( my $i = 0 ; $i < scalar(@fragments) ; $i++ ) {
                my $fragment = $fragments[$i];
                my @variants = ( keys %{ $scoredFragments{$ref}{$method}{$fragment}{$variantHit} } );
                my $variant  = $variants[0];
                my $success  = $scoredFragments{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit};
                my $annealAmplicon   = $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-sseq'};
                my $annealFragment   = $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-qseq'};
                my $mismatch = $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-mismatch'};
                my $length;
                my $obtainLength = $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-sseq'};
                if ( $obtainLength eq "FALSE" ) { $length = "FALSE" }
                else { $length = length( $foundFragmentsUpdated{$ref}{$method}{$fragment}{$variantHit}{$variant}{$blastHit}{'-sseq'} ) }
                print $OUT_DETAIL "\t$success\t$variant\t$annealAmplicon\t$annealFragment\t$mismatch\t$length";

            }

            print $OUT_DETAIL "\n";

        }
    }

    close $OUT_DETAIL;
    return;

}    # sub writeDetailedOutput

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# writeClusterOutput
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub writeClusterOutput {

    my ($input) = @_;

	open( my $OUT_CLUSTER, ">", "$input" ) || die "Could not open summary output file for cluster results:$!\n";
	
	# Print run information
	&writeInfo($OUT_CLUSTER);
	
	# Print no data if no clustering was performed
	if ( $params{'-cluster'} eq "no" ) {
	    print $OUT_CLUSTER "No clustering performed...\n";
	} else {
	
	    # Make a nice header line
	    print $OUT_CLUSTER "# Method\tCluster number\tReference\tAmplicon sequence\n";
	    foreach my $method ( keys %foundClusters ) {
	        for ( my $cluster = 0 ; $cluster <= $maxClusters{$method} ; $cluster++ ) {
	            foreach my $ref ( keys %{ $foundClusters{$method}{$cluster} } ) {
	                my $seq = $foundAmplicons{$ref}{$method}{"-amplicon"}{"1"}{"-sseq"};
	                print $OUT_CLUSTER "$method\t$cluster\t$ref\t$seq\n";
	            }
	        }
	    }
	}
	close $OUT_CLUSTER;
	return;
		
} # sub writeClusterOutput

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# writeInfo
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub writeInfo {

	my ($fileHandle) = @_;
	
	# Create checkums for thei input files
	my $ref_md5 = md5_hex( $params{'-ref'} );
	my $config_md5 = md5_hex( $params{'-config'} );
	
	print $fileHandle "# Screened version 1.0\n";
	print $fileHandle "# Input fasta file containing reference genome: $params{'-ref'} (MD5 checksum: $ref_md5)\n";
	print $fileHandle "# Input file containing method information: $params{'-config'} (MD5 checksum: $config_md5)\n";
	print $fileHandle "# Method(s) evaluated: $params{'-method'}\n";
	print $fileHandle "# Clustering method employed: $params{'-cluster'}\n";
	print $fileHandle "# Name of output file containing detailed information: $params{'-out_detail'}\n";
	print $fileHandle "# Name of output file containing summary information: $params{'-out_summary'}\n";
	print $fileHandle "# Name of output file containing cluster information: $params{'-out_cluster'}\n";
	print $fileHandle "# Consider degenerate fragments: $params{'-deg'}\n";
	print $fileHandle "# Perform amplicon extension: $params{'-ext_am'}\n";
	print $fileHandle "# Perform fragment extension: $params{'-ext_frag'}\n";
	print $fileHandle "# Maximum percentage of mismatches allowed in fragment annealing sites: $params{'-mismatch'}\n";
	print $fileHandle "# Minimum percentage of alignment length required in fragment annealing sites: $params{'-length'}\n";
	print $fileHandle "# Number of bases to consider at the 3' end for forward and reverse primer to check if there are mismatches: $params{'-end'}\n";
	print $fileHandle "# Number of mismatches in the region specified by the '-end' parameter that is conisdered to prevent elongation in the forward and reverse primer: $params{'-bases'}\n";
	
	return;
	
} # sub writeInfo