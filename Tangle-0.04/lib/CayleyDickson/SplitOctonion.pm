package SplitOctonion;
use base qw(CayleyDickson);

#
# SplitOctonion.pm - represent and operate on dual quaternion numbers
#
# ref: https://en.wikipedia.org/wiki/Split_quaternion
#
# Generally, split octonions are like sedenions where each coefficient is a splitquaternion ...
#
#

sub i_squared { shift->is_sedenion ? 0 : 1 }

1;

__END__

