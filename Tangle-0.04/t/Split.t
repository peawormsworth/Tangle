#!/usr/bin/perl -wT
use Test::More tests => 45;

use 5.010;
use warnings;
use strict;
use lib qw(./lib/Tangle ./lib/CayleyDickson);

use Split;
use Data::Dumper;

use constant DEBUG     => 0;
use constant VERBOSE   => 0;
use constant PRECISION => 0.0001;
use constant PACKAGE   => 'Split';

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
   $a = Split->new(@$set);
   ok(leq([$a->flat], $set), sprintf "Split->new(%16s) = $a", join(',', @$set));
}

print "\n###\n### tensor() tests ...\n###\n\n" if VERBOSE;
foreach my $set (
   [ [sqrt(1/2),sqrt(1/2)], [1,0], [sqrt(1/2),0,sqrt(1/2),0] ]
) {
   $a = Split->new(@{$set->[0]});
   $b = Split->new(@{$set->[1]});
   $t = Split->new(@{$set->[2]});
   ok(leq([$a->tensor($b)->flat], [$t->flat]), sprintf("[%s] × [%s] = [%s]", otss($a, $b, $t)));
}

$a = Split->new(1,0);
$b = Split->new(0,-1);
$c = Split->new(1,3);
$d = Split->new(2,5);

print "\n###\n### add() test ...\n###\n\n" if VERBOSE;

ok(leq([($a+$b)->flat],[1,-1]), sprintf("[%s] + [%s] = [%s]", otss($a,$b,($a+$b))));

print "\n###\n### Quaternion tests ...\n###\n\n" if VERBOSE;

$a = Split->new(3.4,5.34,-0.28,1.239);
$b = Split->new(7.34,-6.17,-6.11,9.84);
$c = Split->new(2,4,8,-2);
$d = Split->new(1,0,2,-1);

printf "calculated a+b = %s\n", $a+$b if VERBOSE;
ok(j([($a+$b)->flat]) eq j([10.74,-0.83,-6.39,11.079]), sprintf("[%s] + [%s] = [%s]", otss($a,$b,$a+$b)));

printf "calculated a-b = %s\n", $a-$b if VERBOSE;
ok(j([($a-$b)->flat]) eq j([-3.94,11.51,5.83,-8.601]), sprintf("[%s] - [%s] = [%s]", otss($a,$b,$a+$b)));

printf "calculated a*b = %s\n", $a*$b if VERBOSE;
ok(leq([($a*$b)->flat],[-7.9918, 18.2176, 37.36103, 8.19525999999999]), sprintf("[%s] * [%s] = [%s]", otss($a,$b,$a*$b)));

printf "calculated c+d = %s\n", $c/$d if VERBOSE;
ok(leq([($c+$d)->flat],[3,4,10,-3]),       sprintf("[%s] + [%s] = [%s]", otss($c,$d,$c+$d)));

printf "calculated c-d = %s\n", $c/$d if VERBOSE;
ok(leq([($c-$d)->flat],[1,4,6,-1]),        sprintf("[%s] - [%s] = [%s]", otss($c,$d,$c-$d)));

printf "calculated c/d = %s\n", $c/$d if VERBOSE;
ok(leq([($c/$d)->flat],[1/3,2/3,4/3,-4/3]), sprintf("[%s] / [%s] = [%s]", otss($c,$d,$c/$d)));

printf "calculated c*d = %s\n", $c*$d if VERBOSE;
ok(leq([($c*$d)->flat],[16, 8, 8, 4]),      sprintf("[%s] * [%s] = [%s]", otss($c,$d,$c*$d)));


print "\n###\n### Octonion Tests ...\n###\n\n" if VERBOSE;

$a = Split->new(1,2,3,4,5,6,7,8);
$b = Split->new(8,7,6,5,4,3,2,1);

printf "calculated a+b = %s\n", $a+$b if VERBOSE;
printf "calculated a-b = %s\n", $a-$b if VERBOSE;
printf "calculated a/b = %s\n", $a/$b if VERBOSE;

ok(leq([($a*$b)->flat],[16, 14, 84, 82, 8, 42, 4, 38]), sprintf("[%s] * [%s] = [%s]", otss($a,$b,$a*$b)));
ok(leq([($a+$b)->flat],[9, 9, 9, 9, 9, 9, 9, 9]), sprintf("[%s] + [%s] = [%s]", otss($a,$b,$a+$b)));
ok(leq([($a-$b)->flat],[-7, -5, -3, -1, 1, 3, 5, 7]), sprintf("[%s] - [%s] = [%s]", otss($a,$b,$a-$b)));
ok(leq([($a/$b)->flat],[0, 0.088235294117647, -0.176470588235294, -0.0882352941176471, 0.352941176470588, 0.264705882352941, 0.529411764705882, 0.441176470588235]), sprintf("[%s] / [%s] = [%s]", otss($a,$b,$a/$b)));

$b = Split->new(8,7,6,5,4,3,2,1);

printf "calculated a+b = %s\n", $a+$b if VERBOSE;
printf "calculated a-b = %s\n", $a-$b if VERBOSE;
printf "calculated a*b = %s\n", $a*$b if VERBOSE;
printf "calculated a/b = %s\n", $a/$b if VERBOSE;

ok(leq([($a + $b)->flat],[9,9,9,9,9,9,9,9          ]), sprintf("[%s] + [%s] = [%s]", otss($a, $b, $a+$b)));
ok(leq([($a - $b)->flat],[-7,-5,-3,-1,1,3,5,7      ]), sprintf("[%s] - [%s] = [%s]", otss($a, $b, $a-$b)));
ok(leq([($a * $b)->flat],[16, 14, 84, 82, 8, 42, 4, 38]), sprintf("[%s] * [%s] = [%s]", otss($a, $b, $a*$b)));
ok(leq([($a / $b)->flat],[0, 0.088235294117647, -0.176470588235294, -0.0882352941176471, 0.352941176470588, 0.264705882352941, 0.529411764705882, 0.441176470588235]), sprintf("[%s] / [%s] = [%s]", otss($a, $b, $a/$b)));


printf "calculated b+a = %s\n", $a+$b if VERBOSE;
printf "calculated b-a = %s\n", $a-$b if VERBOSE;
printf "calculated b*a = %s\n", $a*$b if VERBOSE;
printf "calculated b/a = %s\n", $a/$b if VERBOSE;

ok(leq([($b + $a)->flat],[9,9,9,9,9,9,9,9  ]), sprintf("[%s] + [%s] = [%s]", otss($b, $a, $b+$a)));
ok(leq([($b - $a)->flat],[7,5,3,1,-1,-3,-5,-7]), sprintf("[%s] - [%s] = [%s]", otss($b, $a, $b-$a)));
ok(leq([($b * $a)->flat],[16, 32, -24, -8, 80, 60, 112, 92]), sprintf("[%s] * [%s] = [%s]", otss($b, $a, $b*$a)));
ok(leq([($b / $a)->flat],[0, -0.088235294117647, 0.176470588235294, 0.0882352941176471, -0.352941176470588, -0.264705882352941, -0.529411764705882, -0.441176470588235]), sprintf("[%s] / [%s] = [%s]", otss($b, $a, $b/$a)));

print "\n###\n### Octonion Tests ...\n###\n\n" if VERBOSE;

$a = Split->new(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);
$b = Split->new(16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1);

printf "calculated a+b = %s\n", $a+$b if VERBOSE;
printf "calculated a-b = %s\n", $a-$b if VERBOSE;
printf "calculated a*b = %s\n", $a*$b if VERBOSE;
printf "calculated a/b = %s\n", $a/$b if VERBOSE;

ok(leq([($a + $b)->flat],[17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17]), sprintf("[%s] + [%s] = [%s]", otss($a, $b, $a+$b)));
ok(leq([($a - $b)->flat],[-15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15]), sprintf("[%s] - [%s] = [%s]", otss($a, $b, $a-$b)));
ok(leq([($a * $b)->flat],[32, 98, -108, -42, 296, 226, 428, 358, 16, 14, 284, 282, 8, 142, 4, 138]), sprintf("[%s] * [%s] = [%s]", otss($a, $b, $a*$b)));
ok(leq([($a / $b)->flat],[-3.46944695195361e-18, -0.0227272727272727, 0.136363636363636, 0.113636363636364, -0.0909090909090909, -0.0227272727272727, -0.136363636363636, -0.0681818181818182, 0.181818181818182, 0.204545454545455, 0.0454545454545455, 0.0681818181818182, 0.272727272727273, 0.204545454545454, 0.318181818181818, 0.25]), sprintf("[%s] / [%s] = [%s]", otss($a, $b, $a/$b)));

printf "calculated a+b = %s\n", $a+$b if VERBOSE;
printf "calculated a-b = %s\n", $a-$b if VERBOSE;
printf "calculated a*b = %s\n", $a*$b if VERBOSE;
printf "calculated a/b = %s\n", $a/$b if VERBOSE;

ok(leq([($b + $a)->flat],[17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17]), sprintf("[%s] + [%s] = [%s]", otss($b, $a, $b+$a)));
ok(leq([($b - $a)->flat],[15, 13, 11, 9, 7, 5, 3, 1, -1, -3, -5, -7, -9, -11, -13, -15]), sprintf("[%s] - [%s] = [%s]", otss($b, $a, $b-$a)));
ok(leq([($b * $a)->flat],[32, -4, 232, 196, -112, -12, -184, -84, 288, 320, 80, 112, 416, 312, 480, 376]), sprintf("[%s] * [%s] = [%s]", otss($b, $a, $b*$a)));
ok(leq([($b / $a)->flat],[-3.46944695195361e-18, 0.0227272727272727, -0.136363636363636, -0.113636363636364, 0.0909090909090909, 0.0227272727272727, 0.136363636363636, 0.0681818181818182, -0.181818181818182, -0.204545454545455, -0.0454545454545455, -0.0681818181818182, -0.272727272727273, -0.204545454545454, -0.318181818181818, -0.25]), sprintf("[%s] / [%s] = [%s]", otss($b, $a, $b/$a)));


sub j    { join ', ', @{$_[0]}    }
sub ots  { join ', ', $_[0]->flat }
sub otss { map ots($_), @_        }


1;

__END__

