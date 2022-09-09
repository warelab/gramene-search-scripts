#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

# input to this script comes from an API call
# try this: curl "https://data.sorghumbase.org/sorghum_v3/search?q=taxonomy__ancestors:4557&fq=biotype:protein_coding&fl=id,taxon_id,homology__all_orthologs&wt=csv&rows=1000000" | awk -F "," '$2!=45770000' > sorghum_v3.orthologs.csv
# also create a blacklist of genes to avoid when choosing representatives
# echo "select g1.stable_id, g2.stable_id from gene_member g1, gene_member g2, homology_member hm, homology_member hm2, homology h, gene_tree_node_attr gtna where gtna.node_type=\"gene_split\" and gtna.node_id = h.gene_tree_node_id and h.homology_id = hm.homology_id and hm.gene_member_id = g1.gene_member_id and hm2.homology_id = h.homology_id and hm2.gene_member_id = g2.gene_member_id and hm.gene_member_id > hm2.gene_member_id" | mysql -h HOSTNAME -u USERNAME -pPASSWORD ensembl_compara_3_87_sorghum1021 | awk 'NR>1{print $1"\n"$2}' |sort | uniq > sorghum_v3.split_genes.txt

my %blacklist;
my %pan; # $pan{id} = rep
my %hasTax;
my %taxa;
my %isRep;
my $bl_file = shift @ARGV;
open(my $fh, "<", $bl_file);
while (<$fh>) {
    chomp;
    $blacklist{$_} = 1;
}
close $fh;

while (<>) {
    chomp;
    s/\"//g;
    my ($id, $taxon_id, @orthologs) = split /,/, $_;
    next if ($id eq "id");
    $taxa{$taxon_id} = 1;
    $hasTax{$id}=$taxon_id;    
    next if $pan{$id};
    if ($blacklist{$id}) {
        print STDERR "$id has been blacklisted\n";
        next;
    }
    $isRep{$id}=\@orthologs;
    for my $o (@orthologs) {
        $pan{$o} = $id;
    }
}
print join ("\t", 'rep_id','taxon_id','set_size','members','n_taxa', join(',',sort keys %taxa)),"\n";
for my $rep (keys %isRep) {
    my @orthologs = grep {$hasTax{$_} and $_ ne $rep} @{$isRep{$rep}};
    my %taxTally;
    for my $tid (keys %taxa) {
        $taxTally{$tid} = 0;
    }
    for my $id (@{$isRep{$rep}}) {
        next unless $hasTax{$id};
        $taxTally{$hasTax{$id}}++;
    }
    if (@orthologs == 0) {
        $taxTally{$hasTax{$rep}}++;
    }
    my $n_taxa=0;
    for my $tid (keys %taxTally) {
        $n_taxa++ if ($taxTally{$tid} > 0);
    }
    print join("\t"
    , $rep
    , $hasTax{$rep}
    , scalar @orthologs + 1
    , join(',',@orthologs)
    , $n_taxa
    , join(',',map {$taxTally{$_}} sort keys %taxa)
    ),"\n";
}