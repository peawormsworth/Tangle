#!/usr/bin/perl -wT
#
# 120_unit_icosian.pl - generate the 120 Hamiltonion Quaterions of the 600-cell 
#  then rotate each of its points by 180 around a couple axis
#
# ref: https://en.wikipedia.org/wiki/Icosian
#
use strict;
use lib qw(.);
use CayleyDickson;
use Data::Dumper;

use constant GOLDEN_RATIO     => (sqrt(5)+1)/2;
use constant GOLDEN_CONJUGATE => -1/GOLDEN_RATIO;
use constant GR => GOLDEN_RATIO;
use constant GC => GOLDEN_CONJUGATE;

my @icosians;

foreach my $digit (0,1,2,3) {
   foreach my $polarity (+1, -1) {
      my @digits = ((0) x $digit, 2 * $polarity, (0) x (3 - $digit));
      #push @icosians, CayleyDickson->new((0) x $digit,	$polarity, (0) x (3 - $digit));
      push @icosians, CayleyDickson->new(@digits) / 2;
   }
}


foreach my $p1 (1,-1) {
   foreach my $p2 (1, -1) {
      foreach my $p3 (1,-1) {
         foreach my $p4 (1, -1) {
            push @icosians, CayleyDickson->new( $p1, $p2, $p3, $p4 ) / 2;
         }
      }
   }
}

foreach my $v1 ( 0, 1, GR, GC ) {
   my %set;
   foreach my $s1 ( +1, -1 ) {
      $set{$v1} = $s1;
      foreach my $v2 ( 0, 1, GR, GC ) {
         next if $v2 == $v1;
         foreach my $s2 ( +1, -1 ) {
            $set{$v2} = $s2;
            foreach my $v3 ( 0, 1, GR, GC ) {
               next if ($v3 == $v2 or $v3 == $v1);
               foreach my $s3 ( +1, -1 ) {
                  $set{$v3} = $s3;
                  foreach my $v4 ( 0, 1, GR, GC ) {
                     next if ($v4 == $v3 or $v4 == $v2 or $v4 == $v1);
                     foreach my $s4 ( +1, -1 ) {
                        $set{$v4} = $s4;
                        next if $set{0} == -1;
                        next if $set{GR()} == $set{GC()};
                        my $q = CayleyDickson->new( $s1*$v1, $s2*$v2, $s3*$v3, $s4*$v4 ) / 2;
                        push @icosians, $q;
                     }
                  }
               }
            }
         }
      }
   }
}

print "The 120 Hamiltonian Quaterions representing the vertex points of the 600-cell:\n\n";

my $i = 0;
foreach my $v (@icosians) {
   $i++;
   printf "vertex %d: [%s]\n", $i, join(', ', $v->flat);
}

print "\nNow we rotate the 600-cell:\n\n";

my $rotation = sqrt(1/4) * CayleyDickson->new(0,0,1,-1);
my @rotated;
foreach my $v (@icosians) {
   push @rotated, $rotation * $v * 1/$rotation;
}

$i = 0;
foreach my $v (@rotated) {
   $i++;
   printf "vertex %d: [%s]\n", $i, join(', ', $v->flat);
}



#d(icosians => \@icosians);
printf "Total elements: %d\n", scalar @icosians;

sub d {
   my %a = @_;
   my @k = keys %a;
   my $d = Data::Dumper->new([@a{@k}],[@k]); $d->Purity(1)->Deepcopy(1); print $d->Dump;
}


1;

__END__

