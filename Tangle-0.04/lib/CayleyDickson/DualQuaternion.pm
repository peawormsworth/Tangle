package DualQuaternion;
use base qw(CayleyDickson);

#
# DualQuaternion.pm - represent and operate on dual quaternion numbers
#
# ref: https://en.wikipedia.org/wiki/Dual_quaternion
#
# Generally, dual quaternion numbers are like complex numbers in 2 dimensions and quaternion in 4 dimensions,
# but in 8 dimensions it is different from Octonions, in that the square of the imaginary is 0 (instead of normal -1).
#

sub i_squared { shift->is_octonion ? 0 : -1 }

1;

__END__

