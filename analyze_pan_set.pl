#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my $species = shift @ARGV;
my $orthologs_file = shift @ARGV;
my $pan_file = shift @ARGV;

my %outgroups = (
    sorghum => [
        {
            label => "Viridiplanteae",
            re => qr/^AT/
        },
        {
            label => "Poaceae",
            re => qr/^Os/
        },
        {
            label => "Andropogoneae",
            re => qr/^Zm/
        }
    ],
    maize => [
        {
            label => "Viridiplanteae",
            re => qr/^AT/
        },
        {
            label => "Poaceae",
            re => qr/^Os/
        },
        {
            label => "Andropogoneae",
            re => qr/^SORBI/
        }
    ],
    rice => [
        {
            label => "Viridiplanteae",
            re => qr/^AT/
        },
        {
            label => "Poaceae",
            re => qr/^Zm/
        },
        {
            label => "Poaceae",
            re => qr/^SORBI/
        }
    ],
    grapevine => [
        {
            label => "Viridiplanteae",
            re => qr/^Zm/
        },
        {
            label => "Viridiplanteae",
            re => qr/^SORBI/
        },
        {
            label => "Viridiplanteae",
            re => qr/^Os/
        },
        {
            label => "Rosids",
            re => qr/^AT/
        }
    ]
);

$species and $outgroups{$species} and -e $orthologs_file and -e $pan_file or die "usage: $0 <species> orthologs.csv pan_set.txt";
print STDERR "reading $pan_file\n";
my %pan;
open (my $fh, "<", $pan_file);
my $header = <$fh>;
chomp $header;
my @cols = split /\t/, $header;
my @taxa = split /,/, $cols[-1];
while (<$fh>) {
    chomp;
    my ($id,$tid,$set_size,$members,$n_taxa,$tax_tally) = split /\t/, $_;
    my @member_list = split /,/, $members;
    my %tt;
    my @tax_tally = split /,/, $tax_tally;
    for(my $i=0;$i<@tax_tally; $i++) {
        $tt{$taxa[$i]} = $tax_tally[$i];
    }
    $pan{$id} = {
        taxon_id => $tid,
        set_size => $set_size,
        members => \@member_list,
        n_taxa => $n_taxa,
        tax_tally => \%tt
    };
    
}
close $fh;

my @tax_order;
my $prev_tax = 0;
print STDERR "reading orthologs\n";
open ($fh, "<", $orthologs_file);
$header = <$fh>;
my %age;
while (<$fh>) {
    s/"//g;
    chomp;
    my ($id,$tid,@orthologs) = split /,/, $_;
    if ($tid != $prev_tax) {
        $prev_tax = $tid;
        push @tax_order, $tid;
    }
    for my $rule (@{$outgroups{$species}}) {
        last if $age{$id};
        my $re = $rule->{re};
        for my $gid (@orthologs) {
            if ($gid =~ m/$re/) {
                $age{$id} = $rule->{label};
                last;
            }
        }
    }
    $age{$id} ||= $species;
    print "AGE\t$id\t$age{$id}\n";
}
close $fh;

print STDERR "munging\n@tax_order\n";
my %pan_tally;
my %genome_tally;
for my $id (keys %pan) {
    $genome_tally{$pan{$id}{n_taxa}}{$age{$id}}++;
    $pan_tally{$pan{$id}{taxon_id}}{$age{$id}}++;
}

my @labels = map {$_->{label}} @{$outgroups{$species}};
push @labels, $species;
print join("\t", "Genome",@labels),"\n";
for my $tid (@tax_order) {
    print join ("\t", $tid, map {$pan_tally{$tid}{$_}} @labels), "\n";
}

print "\n\n", join("\t", "Pan_set_size", @labels),"\n";
for my $pss (sort {$b <=> $a} keys %genome_tally) {
    print join ("\t", $pss, map {$genome_tally{$pss}{$_}} @labels), "\n";
}