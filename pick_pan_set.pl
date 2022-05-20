#!/usr/bin/env perl

my %pan; # $pan{id} = rep
my %hasTax;
my %isRep;
while (<>) {
    chomp;
    s/\"//g;
    my ($id, $taxon_id, @orthologs) = split /,/, $_;
    next if ($id eq "id");
    $hasTax{$id}=$taxon_id;    
    next if $pan{$id};
    $isRep{$id}=\@orthologs;
    for my $o (@orthologs) {
        $pan{$o} = $id;
    }
}
for my $rep (keys %isRep) {
    my @orthologs = grep {$hasTax{$_}} @{$isRep{$rep}};
    print join("\t", $rep, $hasTax{$rep}, scalar @orthologs || 1),"\n";
}