#!/usr/bin/env perl
use Modern::Perl;
use autodie;
use File::HomeDir;

# directories
my $original = File::HomeDir->my_home . "/Resources/FunctionalGenomics/data/ENCODE/original/";
my $processed = File::HomeDir->my_home . "/Resources/FunctionalGenomics/data/ENCODE/processed/";

# files
opendir my $dh, $original . "DHS/";
my @DHSlinksFiles = grep { /^genomewideCorrs_above0\.75?_promoterPlusMinus500kb_withGeneNames_3[25]celltypeCategories\.bed8$/ } readdir($dh);
closedir($dh);
my $DHSfile = $original . "DHS" . "/" . "files.txt";
my $TFBSfile = $original . "TFBS" . "/" . "files.txt";
my $ChIAfile = $original . "ChIA-PET/files.txt";

# data structures
my %DHSs;
my %d2g; # DHS 2 gene
my %TFBSs;
my %ChIA;

# clean up
unlink glob $processed . "DHS" . "/" . "*";
unlink glob $processed . "TFBS" . "/" . "*";
unlink glob $processed . "ChIA-PET" . "/" . "*";

## read DHSs (keep cells with treatment)
#open my $Dfh, "<", $DHSfile;
#while (<$Dfh>) {
#	chomp;
#	my @flds = split /\s/; 
#	my $filename = $flds[0] =~ s/\.gz//r;
#	my $cell = $flds[6] =~ s/cell=(.+);/$1/r;
#	my $treatment = $flds[7] =~ s/treatment=(.+);/$1/r;
#	my $label;
#	if ($treatment =~ /^None$/) {
#		$label = $cell;
#	}
#	else {
#		$label = $cell . "_". $treatment;
#	}
#	my $oldfile = $original . "DHS" . "/" . $filename;
#	$DHSs{$label} = $oldfile;
#}
#close $Dfh;

# read DHSs (exclude cells with treatment)
open my $Dfh, "<", $DHSfile;
while (<$Dfh>) {
	chomp;
	my @flds = split /\s/; 
	my $filename = $flds[0] =~ s/\.gz//r;
	my $cell = $flds[6] =~ s/cell=(.+);/$1/r;
	my $treatment = $flds[7] =~ s/treatment=(.+);/$1/r;
	my $label;
	if ($treatment =~ /^None$/) {
		$label = $cell;
		my $oldfile = $original . "DHS" . "/" . $filename;
		$DHSs{$label} = $oldfile;
	}
}
close $Dfh;

# read DHS links (distal enhancer DHSs to promoter DHSs associations)
for my $DHSlinksFile (@DHSlinksFiles) { 
	open my $Lfh, "<", $original . "DHS/" . $DHSlinksFile;
	while (<$Lfh>) {
		chomp;
		my @flds = split /\t/;
		my $promoter = $flds[0] . ":" . $flds[1] . "-" . $flds[2];
		my $gene = $flds[3];
		my $enhancer = $flds[4] . ":" . $flds[5] . "-" . $flds[6];
		# enhancers and promoters alike
		$d2g{$promoter}{$gene}++;
		$d2g{$enhancer}{$gene}++;
	}
	close $Lfh;
}

# write DHSs
for my $label (keys %DHSs) {
	my $oldfile = $DHSs{$label};
	my $newfile = $processed . "DHS" . "/" . $label . ".bed";
	open my $in, "<", $oldfile;
	open my $out, ">", $newfile;
	while (<$in>) {
		my @flds = split /\t/;
		my $chr = $flds[0];
		my $start = $flds[1];
		my $end = $flds[2];
		my $dhs = $chr . ":" . $start . "-" . $end;
		if (exists $d2g{$dhs}) {
			my @genes;
			for my $gene (keys %{ $d2g{$dhs} }) {
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

## read TFBSs (keep cells with treatment)
#open my $Tfh, "<", $TFBSfile;
#while (<$Tfh>) {
#	chomp;
#	my @flds = split /\s/; 
#	my $filename = $flds[0] =~ s/\.gz//r;
#	my $cell = $flds[6] =~ s/cell=(.+);/$1/r;
#	my $treatment = $flds[7] =~ s/treatment=(.+);/$1/r;
#	my $label;
#	if ($treatment =~ /^None$/) {
#		$label = $cell;
#	}
#	else {
#		$label = $cell . "_" . $treatment;
#	}
#	my $antibody = $flds[8] =~ s/antibody=(.+?)([_(].+)?;/$1/r;
#	$antibody = "Pol2" if $antibody =~ /^Pol2/;
#	my $oldfile = $original . "TFBS" . "/" . $filename;
#	$TFBSs{$label}{$antibody}{$oldfile}++;
#}
#close $Tfh;

# read TFBSs (exclude cells with treatment)
open my $Tfh, "<", $TFBSfile;
while (<$Tfh>) {
	chomp;
	my @flds = split /\s/; 
	my $filename = $flds[0] =~ s/\.gz//r;
	my $cell = $flds[6] =~ s/cell=(.+);/$1/r;
	my $treatment = $flds[7] =~ s/treatment=(.+);/$1/r;
	my $label;
	if ($treatment =~ /^None$/) {
		$label = $cell;
		my $antibody = $flds[8] =~ s/antibody=(.+?)([_(].+)?;/$1/r;
		$antibody = "Pol2" if $antibody =~ /^Pol2/;
		my $oldfile = $original . "TFBS" . "/" . $filename;
		$TFBSs{$label}{$antibody}{$oldfile}++;
	}
}
close $Tfh;

# write TFBSs
for my $label (keys %TFBSs) {
	my $newfile = $processed . "TFBS" . "/" . $label . ".bed";
	open my $out, ">>", $newfile;
	for my $antibody (keys %{ $TFBSs{$label} }) {
		for my $oldfile (keys %{ $TFBSs{$label}{$antibody} }) {
			open my $in, "<", $oldfile;
			while (<$in>) {
				chomp;
				my @flds = split /\t/;
				my $chr = $flds[0];
				my $start = $flds[1];
				my $end = $flds[2];
				print $out $chr . "\t" . $start . "\t" . $end . "\t" . $antibody . "\n";
			}
			close $in;
		}
	}
	close $out;
}

# read ChIA-PET
open my $Cfh, "<", $ChIAfile;
while (<$Cfh>) {
	my @flds = split /\s/;
	my $filename;
	if ($flds[0] =~ /.+\.bed\.gz$/) {
		$filename = $flds[0] =~ s/\.gz//r;
	}
	else {
		next;
	}
	my $cell = $flds[7] =~ s/cell=(.+);/$1/r;
	my $antibody = $flds[8] =~ s/antibody=(.+?)([_(].+)?;/$1/r;
	$antibody = $antibody . ".";
	my $oldfile = $original . "ChIA-PET" . "/" . $filename;
	my $label = $antibody . $cell;
	$ChIA{$label} = $oldfile;
}
close $Cfh;

# write ChIA-PET
for my $label (keys %ChIA) {
	my $oldfile = $ChIA{$label};
	my $newfile = $processed . "ChIA-PET" . "/" . $label . ".bed";
	open my $in, "<", $oldfile;
	open my $out, ">>", $newfile;
	while (<$in>) {
		my @flds = split /\t/;
		my @flds2 = split /:|\.\.|-|,/, $flds[3];
		my $chr = $flds2[0];
		my $start = $flds2[1];
		my $end = $flds2[2];
		my $name = $flds2[3] . ":" . $flds2[4] . "-" . $flds2[5];
		print $out $chr . "\t" . $start . "\t" . $end . "\t" . $name . "\n";
	}
	close $in;
	close $out;
}

