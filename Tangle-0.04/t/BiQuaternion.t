#!/usr/bin/perl -wT
use Test::More tests => 45;

use 5.010;
use warnings;
use strict;
use lib qw(./lib/Tangle ./lib/CayleyDickson);

use BiQuaternion;
use Data::Dumper;

use constant DEBUG     => 0;
use constant VERBOSE   => 0;
use constant PRECISION => 0.0001;
use constant PACKAGE   => 'BiQuaternion';

sub d {
   my %a = @_;
   my @k = keys %a;
   my $d = Data::Dumper->new([@a{@k}],[@k]); $d->Purity(1)->Deepcopy(1); print $d->Dump;
}

sub leq {
   my @a = @{(shift)};
   my @b = @{(shift)};
   my $success = 0;
   if ($success = scalar @a == scalar @b) {
      foreach my $i (0 .. @a - 1) {
         #$success = $a->[$i] == $b->[$i];
	 #printf "compare [%s] with [%s]: ", $a[$i], $b[$i];
	 #print('FAILED') unless ($success = ($a[$i] - $b[$i] < PRECISION));
	 #die('FAILED') unless $success = ($a[$i] - $b[$i] < PRECISION);
         last unless $success = (($a[$i] - $b[$i] < PRECISION) ? 1 : 0);
	 #print "works\n";
      }
   }
   #print "success: $success\n";
   $success;
}

#printf "%s\n", leq([1,2],[2,1]) ? 'YES' : 'NO';
#printf "%s\n", leq([1,2],[1,2]) ? 'YES' : 'NO';
#exit;

my ($a,$b,$c,$d,$e,$h,$i,$o,$measures,$t);

print "\n###\n### class method tests ...\n###\n\n" if VERBOSE;
foreach my $method (qw(new add subtract multiply divide conjugate inverse norm tensor)) {
   can_ok(PACKAGE, $method);
}


print "\n###\n### new() tests ...\n###\n\n" if VERBOSE;
foreach my $set (
   [1,0,0,0],
   [0,1,],
   [0,0,1,0],
   [0,0,0,1],
   [0,-1,2,0],
   [0,0,0,3],
   [1,4,2,9,0,3,2,5],
) {
   $a = BiQuaternion->new(@$set);
   ok(leq([$a->flat], $set), sprintf "BiQuaternion->new(%16s) = $a", join(',', @$set));
}

print "\n###\n### tensor() tests ...\n###\n\n" if VERBOSE;
foreach my $set (
   [ [sqrt(1/2),sqrt(1/2)], [1,0], [sqrt(1/2),0,sqrt(1/2),0] ]
) {
   $a = BiQuaternion->new(@{$set->[0]});
   $b = BiQuaternion->new(@{$set->[1]});
   $t = BiQuaternion->new(@{$set->[2]});
   ok(leq([$a->tensor($b)->flat], [$t->flat]), sprintf("[%s] × [%s] = [%s]", otss($a, $b, $t)));
}

$a = BiQuaternion->new(1,0);
$b = BiQuaternion->new(0,-1);
$c = BiQuaternion->new(1,3);
$d = BiQuaternion->new(2,5);

print "\n###\n### add() test ...\n###\n\n" if VERBOSE;

ok(leq([($a+$b)->flat],[1,-1]), sprintf("[%s] + [%s] = [%s]", otss($a,$b,($a+$b))));

print "\n###\n### Quaternion tests ...\n###\n\n" if VERBOSE;

$a = BiQuaternion->new(3.4,5.34,-0.28,1.239);
$b = BiQuaternion->new(7.34,-6.17,-6.11,9.84);
$c = BiQuaternion->new(2,4,8,-2);
$d = BiQuaternion->new(1,0,2,-1);

printf "calculated a+b = %s\n", $a+$b if VERBOSE;
ok(j([($a+$b)->flat]) eq j([10.74,-0.83,-6.39,11.079]), sprintf("[%s] + [%s] = [%s]", otss($a,$b,$a+$b)));

printf "calculated a-b = %s\n", $a-$b if VERBOSE;
ok(j([($a-$b)->flat]) eq j([-3.94,11.51,5.83,-8.601]), sprintf("[%s] - [%s] = [%s]", otss($a,$b,$a+$b)));

printf "calculated a*b = %s\n", $a*$b if VERBOSE;
#ok(($a*$b - BiQuaternion->new(44.00124,23.03269,-83.01943,8.19526))->norm < PRECISION, sprintf("[%s] * [%s] = [%s]", otss($a,$b,$a*$b)));
ok(leq([($a*$b)->flat],[44.00124, 23.03269, -83.01943, 8.19525999999999]), sprintf("[%s] * [%s] = [%s]", otss($a,$b,$a*$b)));

printf "calculated c+d = %s\n", $c/$d if VERBOSE;
ok(leq([($c+$d)->flat],[3,4,10,-3]),       sprintf("[%s] + [%s] = [%s]", otss($c,$d,$c+$d)));

printf "calculated c-d = %s\n", $c/$d if VERBOSE;
ok(leq([($c-$d)->flat],[1,4,6,-1]),        sprintf("[%s] - [%s] = [%s]", otss($c,$d,$c-$d)));

printf "calculated c/d = %s\n", $c/$d if VERBOSE;
ok(leq([($c/$d)->flat],[10/3,4/3,0,-4/3]), sprintf("[%s] / [%s] = [%s]", otss($c,$d,$c/$d)));
#3.33333333333333, 1.33333333333333, 0, -1.33333333333333

printf "calculated c*d = %s\n", $c*$d if VERBOSE;
ok(leq([($c*$d)->flat],[-16, 0, 16, 4]),      sprintf("[%s] * [%s] = [%s]", otss($c,$d,$c*$d)));


print "\n###\n### Octonion Tests ...\n###\n\n" if VERBOSE;

$a = BiQuaternion->new(1,2,3,4,5,6,7,8);
$b = BiQuaternion->new(8,7,6,5,4,3,2,1);

printf "calculated a+b = %s\n", $a+$b if VERBOSE;
printf "calculated a-b = %s\n", $a-$b if VERBOSE;
printf "calculated a/b = %s\n", $a/$b if VERBOSE;

ok(leq([($a*$b)->flat],[-68, -16, -36, 52, 68, 76, -36, 136]), sprintf("[%s] * [%s] = [%s]", otss($a,$b,$a*$b)));
ok(leq([($a+$b)->flat],[9, 9, 9, 9, 9, 9, 9, 9]), sprintf("[%s] + [%s] = [%s]", otss($a,$b,$a+$b)));
ok(leq([($a-$b)->flat],[-7, -5, -3, -1, 1, 3, 5, 7]), sprintf("[%s] - [%s] = [%s]", otss($a,$b,$a-$b)));
ok(leq([($a/$b)->flat],[0.274509803921569, 0.303921568627451, 0.137254901960784, 0.264705882352941, -0.352941176470588, 0.441176470588235, 0.176470588235294, 0.441176470588235]), sprintf("[%s] / [%s] = [%s]", otss($a,$b,$a/$b)));


$b = BiQuaternion->new(8,7,6,5,4,3,2,1);

printf "calculated a+b = %s\n", $a+$b if VERBOSE;
printf "calculated a-b = %s\n", $a-$b if VERBOSE;
printf "calculated a*b = %s\n", $a*$b if VERBOSE;
printf "calculated a/b = %s\n", $a/$b if VERBOSE;

ok(leq([($a + $b)->flat],[9,9,9,9,9,9,9,9          ]), sprintf("[%s] + [%s] = [%s]", otss($a, $b, $a+$b)));
ok(leq([($a - $b)->flat],[-7,-5,-3,-1,1,3,5,7      ]), sprintf("[%s] - [%s] = [%s]", otss($a, $b, $a-$b)));
ok(leq([($a * $b)->flat],[-68, -16, -36, 52, 68, 76, -36, 136]), sprintf("[%s] * [%s] = [%s]", otss($a, $b, $a*$b)));
ok(leq([($a / $b)->flat],[0.274509803921569, 0.303921568627451, 0.137254901960784, 0.264705882352941, -0.352941176470588, 0.441176470588235, 0.176470588235294, 0.441176470588235]), sprintf("[%s] / [%s] = [%s]", otss($a, $b, $a/$b)));


printf "calculated b+a = %s\n", $a+$b if VERBOSE;
printf "calculated b-a = %s\n", $a-$b if VERBOSE;
printf "calculated b*a = %s\n", $a*$b if VERBOSE;
printf "calculated b/a = %s\n", $a/$b if VERBOSE;

ok(leq([($b + $a)->flat],[9,9,9,9,9,9,9,9  ]), sprintf("[%s] + [%s] = [%s]", otss($b, $a, $b+$a)));
ok(leq([($b - $a)->flat],[7,5,3,1,-1,-3,-5,-7]), sprintf("[%s] - [%s] = [%s]", otss($b, $a, $b-$a)));
ok(leq([($b * $a)->flat],[-68, -16, 36, 88, -76, 112, 36, 100]), sprintf("[%s] * [%s] = [%s]", otss($b, $a, $b*$a)));
ok(leq([($b / $a)->flat],[0.274509803921569, 0.303921568627451, -0.215686274509804, -0.264705882352941, 0.352941176470588, -0.441176470588235, -0.176470588235294, -0.441176470588235]), sprintf("[%s] / [%s] = [%s]", otss($b, $a, $b/$a)));

print "\n###\n### Octonion Tests ...\n###\n\n" if VERBOSE;

$a = BiQuaternion->new(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);
$b = BiQuaternion->new(16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1);

printf "calculated a+b = %s\n", $a+$b if VERBOSE;
printf "calculated a-b = %s\n", $a-$b if VERBOSE;
printf "calculated a*b = %s\n", $a*$b if VERBOSE;
printf "calculated a/b = %s\n", $a/$b if VERBOSE;

ok(leq([($a + $b)->flat],[17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17]), sprintf("[%s] + [%s] = [%s]", otss($a, $b, $a+$b)));
ok(leq([($a - $b)->flat],[-15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15]), sprintf("[%s] - [%s] = [%s]", otss($a, $b, $a-$b)));
ok(leq([($a * $b)->flat],[-428, -318, -8, 82, -20, 70, -16, 58, 532, 590, 576, -558, 268, 226, -136, 486]), sprintf("[%s] * [%s] = [%s]", otss($a, $b, $a*$b)));
ok(leq([($a / $b)->flat],[0.267379679144385, 0.275401069518717, 0.149732620320856, -0.0294117647058824, 1.38777878078145e-17, 0.181818181818182, 0, 0.272727272727273, -0.363636363636364, 1.38777878078145e-17, 0.0909090909090909, 0.409090909090909, -0.181818181818182, 0.409090909090909, 0.0909090909090909, 0.318181818181818]), sprintf("[%s] / [%s] = [%s]", otss($a, $b, $a/$b)));

printf "calculated a+b = %s\n", $a+$b if VERBOSE;
printf "calculated a-b = %s\n", $a-$b if VERBOSE;
printf "calculated a*b = %s\n", $a*$b if VERBOSE;
printf "calculated a/b = %s\n", $a/$b if VERBOSE;

ok(leq([($b + $a)->flat],[17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17]), sprintf("[%s] + [%s] = [%s]", otss($b, $a, $b+$a)));
ok(leq([($b - $a)->flat],[15, 13, 11, 9, 7, 5, 3, 1, -1, -3, -5, -7, -9, -11, -13, -15]), sprintf("[%s] - [%s] = [%s]", otss($b, $a, $b-$a)));
ok(leq([($b * $a)->flat],[-428, -318, -8, 218, -20, 342, -16, 466, -556, 46, 168, 598, -276, 634, 136, 486]), sprintf("[%s] * [%s] = [%s]", otss($b, $a, $b*$a)));
ok(leq([($b / $a)->flat],[0.267379679144385, 0.275401069518717, 0.0588235294117647, -0.165775401069519, 1.38777878078145e-17, -0.181818181818182, 0, -0.272727272727273, 0.363636363636364, -1.73472347597681e-17, -0.0909090909090909, -0.409090909090909, 0.181818181818182, -0.409090909090909, -0.0909090909090909, -0.318181818181818]), sprintf("[%s] / [%s] = [%s]", otss($b, $a, $b/$a)));

#printf "1/%s = %s\n", $a, 1/$a;
#$a = BiQuaternion->new(sqrt(1/10),sqrt(2/10),sqrt(3/10),sqrt(4/10));
#printf "1/%s = %s\n", $a, 1/$a;


sub j    { join ', ', @{$_[0]}    }
sub ots  { join ', ', $_[0]->flat }
sub otss { map ots($_), @_        }


1;

__END__

