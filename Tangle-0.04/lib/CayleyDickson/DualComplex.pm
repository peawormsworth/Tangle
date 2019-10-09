package DualComplex;
use base qw(CayleyDickson);

#
# DualComplex.pm - represent and operate on dual complex numbers
#
# ref: https://en.wikipedia.org/wiki/Dual-complex_numbers
#
# Generally, dual complex numbers are like complex numbers in 2 dimensions and dual numbers in 4 
#

sub i_squared { shift->is_quaternion ? 0 : -1 }

1;

__END__

