#!/usr/bin/perl -wT
#
# 120_unit_icosian.pl - generate the 120 vertex, 600-cell 
#  using 120 unique quaternions.
#  then rotate each of these points about an axis
#
# ref: https://en.wikipedia.org/wiki/Icosian
#

# make Space is a CayleyDickson class ( hypercomplex calculators)
#
package Space; use base qw(CayleyDickson); 1; 

use strict;
use lib qw(.);
use CayleyDickson;
use Data::Dumper;
use utf8;
use constant DEBUG     => 0;
use constant VERBOSE   => 1;
use constant PRECISION => 0.000000001;



# golden ratio (φ) and golden conjugate (ϕ)
# have a negative inverse relation:
#
#   φ =  –ϕ ⁻¹  and  ϕ =  –φ ⁻¹
#

use constant           POINT =>           0;                                     # zero dimension
use constant            LINE => (     POINT **   POINT          );               #    1 dimension. Depends on CPU philosophy. What is 0⁰? https://youtu.be/r0_mi8ngNnM?t=750
use constant           PLANE => (      LINE  +   LINE           );               #    2 dimension
use constant       SPACETIME => (     PLANE **  PLANE           );               #    4 dimension
use constant     IMAGINATION => (     PLANE **  -LINE           );               #  1/2 dimension. Half rotation between +1 and -1, a distance 1 from 0.
use constant               φ => ((SPACETIME  +   LINE           ) ** IMAGINATION + LINE) * IMAGINATION; # golden ratio
use constant               ϕ => (     -LINE  /   φ              );               # other side of the golden ratio (conjugate)
use constant         SPIN_UP => (      LINE                     );               # spin up
use constant       SPIN_DOWN => (  -SPIN_UP                     );               # spin down
use constant            SPIN => (   SPIN_UP  ,   SPIN_DOWN      );               # polarity (+/-)
use constant       AMPLITUDE => (     POINT  ,   SPIN_UP,  φ , ϕ);               # 4 probability amplitudes to be arranged
use constant        ROTATION => Space->new( POINT, SPIN, POINT  ) * IMAGINATION; # 4d vector representing the rotational axis
use constant      INVOLUTION => SPIN_UP / ROTATION;                              # the inverse of rotation required for completion


my @icosians;

foreach my $dimension (LINE .. SPACETIME) {
   foreach my $polarity (SPIN) {
      my @vector = (POINT) x SPACETIME;
      $vector[$dimension - LINE] =  PLANE * $polarity;
      push @icosians, Space->new(@vector) * IMAGINATION;
   }
}


foreach my $p1 (SPIN) {
   foreach my $p2 (SPIN) {
      foreach my $p3 (SPIN) {
         foreach my $p4 (SPIN) {
            push @icosians, Space->new( $p1, $p2, $p3, $p4 ) * IMAGINATION;
         }
      }
   }
}


foreach my $v1 (AMPLITUDE) {
   foreach my $p1 (SPIN) {
      foreach my $v2 (AMPLITUDE) {
         foreach my $p2 (SPIN) {
            foreach my $v3 (AMPLITUDE) {
               foreach my $p3 (SPIN) {
                  foreach my $v4 (AMPLITUDE) {
                     foreach my $p4 (SPIN) {
                        my %vector;
                        @vector{ $v1, $v2, $v3, $v4 } = ( $p1, $p2, $p3, $p4 );
                        next if keys %vector         != SPACETIME
                             or      $vector{+POINT} == -LINE
                             or      $vector{  +φ  } == $vector{ +ϕ };
                        push @icosians, Space->new( $v1*$p1, $v2*$p2, $v3*$p3, $v4*$p4 ) * IMAGINATION;
                     }
                  }
               }
            }
         }
      }
   }
}

#
# Rotate the object by some amount
#

printf "\nRotate 600-cell by ROTATION: %s\n\n", ROTATION if VERBOSE;

my @rotated;
foreach my $vector (@icosians) {
   push @rotated, ROTATION * $vector * INVOLUTION;
}


#
#  Output ...
#

print "\n\nOriginal 120 Quaterions representing the vertex points of the 600-cell:\n\n";

if (VERBOSE) {
   my $i = 1;
   foreach my $vector (@icosians) {
      printf "vertex %5d: [%s]\n", $i++, basis($vector);
   }
}

printf "\n\n120 Quaterions vertex points rotated by %s\n\n", ROTATION;


if (VERBOSE) {
   my $i = 1;
   foreach my $vector (@rotated) {
      printf "vertex %5d: [%s]\n", $i++, basis($vector);
   }
}

d(icosians => \@icosians) if DEBUG;

printf "Total elements: %d\n", scalar @icosians;


# 
# function tools ...
#

sub basis {
   my $vector = shift;
   join(', ', map(sprintf('  %s%-s', ($_ < 0 ? '-' : '+' ), sprintf('%0.5f', abs($_))), $vector->flat))
}

sub d {
   my %a = @_;
   my @k = keys %a;
   my $d = Data::Dumper->new([@a{@k}],[@k]); $d->Purity(1)->Deepcopy(1); print $d->Dump;
}


1;

__END__
