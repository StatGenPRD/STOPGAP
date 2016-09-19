#!/usr/bin/env perl
use Modern::Perl;
use autodie;
use File::HomeDir;
use Text::Autoformat;

# correlation threshold
my $corr = 0;
# significance threshold
my $pval = 1e-5;

# directories
my $original = File::HomeDir->my_home . "/Resources/FunctionalGenomics/data/FANTOM5/original/";
my $processed = File::HomeDir->my_home . "/Resources/FunctionalGenomics/data/FANTOM5/processed/";

# clean up
unlink glob $processed . "CRM/*";

# enhancer to gene associations
my %e2g;

# read enhancer -> gene association file
open my $fh, "<", $original . "Enhancers/enhancer_tss_associations.bed";
while (<$fh>) {
	chomp;
	next if /^track name=/;
	next if /^#/;
	my @flds = split /\t/;
	my @flds2 = split /;/, $flds[3];
	my $enhancer = $flds2[0];
	my $gene = $flds2[2];
	my $r;
	if (defined $flds2[3]) {
		$r = $flds2[3] =~ s/R://r;
	}
	my $fdr;
	if (defined $flds2[4]) {
		$fdr = $flds2[4] =~ s/FDR://r;
	}
	if (defined $gene && defined $r && defined $fdr && $r > $corr && $fdr < $pval) {
		$e2g{$enhancer}{$gene}++;
	}
}
close $fh;

# open enhancer directory
opendir my $dh, $original . "Enhancers/";
# get cell-specific and tissue-specific filenames
my @files = grep { /^(CL|UBERON):\d{7}_\w+_expressed_enhancers\.bed$/ } readdir($dh);
closedir $dh;

# process cell/tissue-specific enhancers
for my $infile (@files) {
	my $name = $infile =~ s/^(CL|UBERON):\d{7}_(\w+)_expressed_enhancers.bed$/$2/r;
	$name = autoformat $name, { case => 'title' }; # this adds newlines apparently
	$name =~ s/_|\n//g;
	$name .= ".bed";
	$infile = $original . "Enhancers/" . $infile;
	my $outfile = $processed . "CRM/" . $name;
	open my $in, "<", $infile;
	open my $out, ">", $outfile;
	while (<$in>) {
		chomp;
		next if /^track name=/;
		next if /^#/;
		my @flds = split /\t/;
		my $chr = $flds[0];
		my $start = $flds[1];
		my $end = $flds[2];
		my $enhancer = $flds[3];
		if (exists $e2g{$enhancer}) {
			my @genes;
			for my $gene (keys %{ $e2g{$enhancer} }) {
				push @genes, $gene;
			}
			print $out $chr . "\t" . $start . "\t" . $end . "\t" . join(",", @genes) . "\n";
		}
		else {
			print $out $chr . "\t" . $start . "\t" . $end . "\t" . "." . "\n";
		}
	}
	close $in;
	close $out;
}
