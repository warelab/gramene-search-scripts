#!/usr/bin/env perl

/*
input is a comma separated file with gene_id, taxon_id, and orthologs as it comes out of the gramene api query.
For example:

curl https://data.gramene.org/yeast_v1/search\?q\=system_name:yarr\*\&fl\=id,taxon_id,homology__all_orthologs\&wt\=csv\&rows\=36000 > yali_id.taxon_id.orthologs

Then just read that file on stdin (or pipe it strait out of the API)

The output file has a group identifier and the number of genes from each queried genome in the orthology group.
If a gene has no orthologs because it was not in compara, it reports a 0 in the corresponding taxon_id column
*/
my %pan; # $pan{id} = rep
my %hasTax;
my %isRep;
my %taxa;
my $header = <>;
while (<>) {
    chomp;
    s/\"//g;
    my ($id, $taxon_id, @orthologs) = split /,/, $_;
    $taxa{$taxon_id} = 1;
    $hasTax{$id}=$taxon_id;    
    next if $pan{$id};
    $isRep{$id}=\@orthologs;
    for my $o (@orthologs) {
        $pan{$o} = $id;
    }
}
my $group_num=0;
print join("\t","orthogroup",sort keys %taxa),"\n";
for my $rep (keys %isRep) {
    $group_num++;
    my @orthologs = grep {$hasTax{$_}} @{$isRep{$rep}};
    my %taxTally;
    for my $tid (keys %taxa) {
        $taxTally{$tid} = 0;
    }
    for my $id (@{$isRep{$rep}}) {
        next unless $hasTax{$id};
        $taxTally{$hasTax{$id}}++;
    }
    print join("\t", "GRP_" . $group_num, map {$taxTally{$_}} sort keys %taxa),"\n";
    # print join("\t", $rep, map {$taxTally{$_}} sort keys %taxa),"\n";
}