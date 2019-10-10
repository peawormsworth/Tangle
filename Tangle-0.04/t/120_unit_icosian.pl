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
use constant DEBUG   => 0;
use constant VERBOSE => 1;


# golden ratio (φ) and golden conjugate (ϕ)
# have a negative inverse relation:
#
#   φ =  –ϕ ⁻¹  and  ϕ =  –φ ⁻¹
#

use constant       POINT => 0;                                # zero dimension
use constant        LINE => ( POINT  **   POINT );            #    1 dimension
use constant          UP => +LINE;
use constant        DOWN => -LINE;
use constant        SPIN => ( UP      ,   DOWN  );            # polarity (+/-)
use constant        AREA => ( LINE    +   LINE  );            #    2 dimension
use constant   SPACETIME => ( AREA   **   AREA  );            #    4 dimension
use constant IMAGINATION => ( AREA   **  -LINE  );            #  1/2 dimension. The place rotated half way between +1 and -1, a distance 1 from 0.
use constant           φ => ((SPACETIME + LINE  ) ** IMAGINATION + LINE) * IMAGINATION; # golden ratio
use constant           ϕ => ( -LINE   /   φ     );            # other side of the golden ratio (conjugate)
use constant     ELEMENT => ( POINT   ,   LINE  ,  φ , ϕ );   # combination of 4 rotated elements
use constant    ROTATION => Space->new(   POINT, POINT, SPIN  ) * IMAGINATION;
use constant  INVOLUTION => 1 / ROTATION;


my @ICOSIANS;

foreach my $dimension (LINE .. SPACETIME) {
   foreach my $polarity (SPIN) {
      my @dimensions = (POINT) x SPACETIME;
      $dimensions[$dimension - LINE] = AREA * $polarity;
      push @ICOSIANS, Space->new(@dimensions) * IMAGINATION;
   }
}


foreach my $p1 (SPIN) {
   foreach my $p2 (SPIN) {
      foreach my $p3 (SPIN) {
         foreach my $p4 (SPIN) {
            push @ICOSIANS, Space->new( $p1, $p2, $p3, $p4 ) * IMAGINATION;
         }
      }
   }
}


foreach my $v1 (ELEMENT) {
   my %set;
   foreach my $p1 (SPIN) {
      $set{$v1} = $p1;
      foreach my $v2 (ELEMENT) {
         next if $v2 == $v1;
         foreach my $p2 (SPIN) {
            $set{$v2} = $p2;
            foreach my $v3 (ELEMENT) {
               next if ($v3 == $v2 or $v3 == $v1);
               foreach my $p3 (SPIN) {
                  $set{$v3} = $p3;
                  foreach my $v4 (ELEMENT) {
                     next if ($v4 == $v3 or $v4 == $v2 or $v4 == $v1);
                     foreach my $p4 (SPIN) {
                        $set{$v4} = $p4;
                        next if $set{+POINT} == -LINE;
                        next if $set{+φ } == $set{+ϕ};
                        push @ICOSIANS, Space->new( $p1*$v1, $p2*$v2, $p3*$v3, $p4*$v4 ) * IMAGINATION;
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

my @rotated;
foreach my $v (@ICOSIANS) {
   push @rotated, ROTATION * $v * INVOLUTION;
}


#
# Output results ...
#

print "120 Quaterions representing the vertex points of the 600-cell:\n\n";

if (VERBOSE) {
   my $i = 1;
   foreach my $v (@ICOSIANS) {
      printf "vertex %5d: [%s]\n", $i++, join(', ', map(sprintf('  %s%-s', ($_ < 0 ? '-' : '+' ), sprintf('%0.5f', abs($_))), $v->flat));
   }
}


print "\nRotated 600-cell by ROTATION:\n\n";

if (VERBOSE) {
   my $i = 1;
   foreach my $v (@rotated) {
      printf "vertex %5d: [%s]\n", $i++, join(', ', map(sprintf('  %s%-s', ($_ < 0 ? '-' : '+' ), sprintf('%0.5f', abs($_))), $v->flat));
   }
}

d(icosians => \@ICOSIANS) if DEBUG;
printf "Total elements: %d\n", scalar @ICOSIANS;

sub d {
   my %a = @_;
   my @k = keys %a;
   my $d = Data::Dumper->new([@a{@k}],[@k]); $d->Purity(1)->Deepcopy(1); print $d->Dump;
}


1;

__END__
