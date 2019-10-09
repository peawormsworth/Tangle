package BiComplex;
use base qw(CayleyDickson);

#
# BiComplex.pm - represent and operate on complex numbers
#
# ref: https://en.wikipedia.org/wiki/Bi-complex_numbers
#
# Generally, bi-complex numbers are like complex numbers in 2 dimensions and quaternion numbers in 4 
#

sub conjugate {
   my $m = shift;

   return $m->SUPER::conjugate unless $m->is_quaternion;
   (ref $m)->new( $m->a, -$m->b )
}


1;

__END__

