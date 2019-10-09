package Dual;
use base qw(CayleyDickson);

#
# DualComplex.pm - represent and operate on dual numbers
#
# ref: https://en.wikipedia.org/wiki/Dual_number
#
# Generally, dual numbers are like complex numbers which have the difference that
# i^2 = 0 (usually: i^2 = -1)
#

sub i_squared { 0 }

1;

__END__

