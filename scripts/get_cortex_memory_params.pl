#!/usr/bin/perl -w

use Getopt::Long;

my $cortex_bub = 0;
my $num_kmers;
my $k;
my $quick = 0; 
my $utilisation = 0.75;
my $returnparams = "nbm";
my $size_per_entry;

&GetOptions(
    'cortex_bub'    => \$cortex_bub,
    'numkmers:i'    => \$num_kmers,
    'k:i'           => \$k,
    'utilisation:i' => \$utilisation,
    'return:s'      => \$returnparams
);

die "You must specify a -numkmers parameter\n" if not defined $num_kmers;
die "You must specify a -k parameter\n" if not defined $k;

if ($k <= 31) {
    $size_per_entry = 16;
} elsif ($k <= 63) {
    $size_per_entry = 24;
} elsif ($k <=95) {
    $size_per_entry = 32;
} else {
    die "Invalid kmer size\n";
}

if ($cortex_bub == 1) {
    $size_per_entry += 8;
}

my ($n, $b, $meg) = choose_n_and_b();

my $output = 0;

if ($returnparams =~ /n/) {
    print $output++ > 0 ? "\t$n":"$n";
}

if ($returnparams =~ /b/) {
    print $output++ > 0 ? "\t$b":"$b";
}

if ($returnparams =~ /m/) {
    print $output++ > 0 ? "\t$meg":"$meg";
}

print "\n";

sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1))
}

sub choose_n_and_b {
    my $best_n = 0;
    my $best_b = 0;
    my $best_m = 0;
    my $best_d = 100;

    for (my $b=50; $b<=100; $b+=5) {
        my $kmers_to_allow = ($num_kmers / $utilisation) / $b;
        my $n = roundup(log($kmers_to_allow)/log(2));
        my $entries = (2**$n)*$b;
        my $bytes = $entries * $size_per_entry;
        my $meg = roundup($bytes / (1024 * 1024));
        my $utilised = $num_kmers / $entries;
        my $diff = abs($utilisation - $utilised);
        #print "n=$n b=$b meg=$meg d=$diff\n";

        if ($diff < $best_d) {
            $best_n = $n;
            $best_b = $b;
            $best_d = $diff;
            $best_m = $meg;
        }
    }
    return ($best_n, $best_b, $best_m);
}
