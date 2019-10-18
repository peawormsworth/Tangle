#!/usr/bin/perl -wT
use Test::More tests => 87;

use 5.010;
use warnings;
use strict;
use lib qw(./lib/CayleyDickson ./lib/Tangle);


use Tangle;
use Data::Dumper;

use constant DEBUG          => 0;
use constant VERBOSE        => 1;
use constant SKIP_CNOT_TEST => 0;
use constant PRECISION      => 10 ** -9;
use constant PACKAGE        => 'Tangle';
use constant METHODS        => qw(new cnot swap x y z i h xx yy zz u d cx cy cz cs not rswap rnot ccnot tips state raw_state measure measures);

sub d {
   my %a = @_;
   my @k = keys %a;
   my $d = Data::Dumper->new([@a{@k}],[@k]); $d->Purity(1)->Deepcopy(1); print $d->Dump;
}

my ($a,$b,$c,$d,$e,$f,$h,$i,$j,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$measures);

diag 'Class method tests ...' if VERBOSE;
foreach my $method (METHODS) {
   can_ok(PACKAGE, $method);
}

diag 'Instantiate tests ...' if VERBOSE;
$a = Tangle->new(1,2);
$e = CayleyDickson->new(1,2);
ok(($a-$e)->norm < PRECISION, "(1,2) : expect: $e, calc: $a");

$a = Tangle->new(3,4);
$e = CayleyDickson->new(3,4);
ok(($a-$e)->norm < PRECISION, "(3,4) : expect: $e, calc: $a");

diag 'Tensor tests ...' if VERBOSE;
$a = Tangle->new(1,2);
$b = Tangle->new(3,4);
$c = $a->tensor($b);
$e = CayleyDickson->new(3,4,6,8);
ok(($c-$e)->norm < PRECISION, " (1,2)⊗ (3,4)                 : expect: $e, calc: $c");


$a = Tangle->new(0,1);
$b = Tangle->new(0,1);
$c = Tangle->new(1,0);
$d = $b->tensor($c);
$e = CayleyDickson->new(0,0,1,0);
ok(($d-$e)->norm < PRECISION, " (0,1)⊗ (1,0)         = |10>  : expect: $e, calc: $d");

$f = $a->tensor($d);
$e = CayleyDickson->new(0,0,0,0,0,0,1,0,0);
ok(($f-$e)->norm < PRECISION, " (0,1)⊗ (0,1)⊗ (1,0)  = |110> : expect: $e, calc: $f");

$a = Tangle->new(1,0);
$b = Tangle->new(1,0);
$c = $b->tensor($a);
$e = CayleyDickson->new(1,0,0,0);
ok(($c-$e)->norm < PRECISION, " (0,1)⊗ (1,0)         = |00>  : expect: $e, calc: $c");

$a = Tangle->new(1,0);
$b = Tangle->new(0,1);
$c = $a->tensor($b);
$e = CayleyDickson->new(0,1,0,0);
ok(($c-$e)->norm < PRECISION, " (1,0)⊗ (1,0)         = |01>  : expect: $e, calc: $c");

$a = Tangle->new(0,1);
$b = Tangle->new(1,0);
$c = $a->tensor($b);
$e = CayleyDickson->new(0,0,1,0);
ok(($c-$e)->norm < PRECISION, " (0,1)⊗ (1,0)         = |10>  : expect: $e, calc: $c");

$a = Tangle->new(0,1);
$b = Tangle->new(0,1);
$c = $a->tensor($b);
$e = CayleyDickson->new(0,0,0,1);
ok(($c-$e)->norm < PRECISION, " (0,1)⊗ (0,1)         = |11>  : expect: $e, calc: $c");

$a = Tangle->new(0,1);
$b = Tangle->new(1,0);
$c = Tangle->new(1,0);
$d = $a->tensor($b)->tensor($c);
$e = CayleyDickson->new(0,0,0,0,1,0,0,0,0);
ok(($d-$e)->norm < PRECISION, " (0,1)⊗ (1,0)⊗ (1,0)  = |100> : expect: $e, calc: $d");

diag 'Gate tests ...' if VERBOSE;

$a = Tangle->new(1,2,3,4);
$a->b->swap;
$e = CayleyDickson->new(1,2,4,3);
#$a->swap;
ok(($a-$e)->norm < PRECISION, "SWAP-B(1,2,3,4)  = (1,2,4,3) : expect: $e, calc: $a");

foreach my $set (
   #['SWAP [1,2,3,4]   = [1,2,4,3]' => [0,1], [0,1],[0,0,1,0],'+1+2i+4j+3k'],
   ['cnot |00> = |00> = [1,0,0,0]' => [1,0], [1,0],[1,0,0,0],'+1'         ],
   ['cnot |01> = |01> = [0,1,0,0]' => [1,0], [0,1],[0,1,0,0],'+0+1i'      ],
   ['cnot |10> = |11> = [0,0,0,1]' => [0,1], [1,0],[0,0,0,1],'+0+1k'      ],
   ['cnot |11> = |10> = [0,0,1,0]' => [0,1], [0,1],[0,0,1,0],'+0+1j'      ],
) {
   my $label = @$set[0];
   my $c = Tangle->new(@{$set->[1]});
   my $t = Tangle->new(@{$set->[2]});
   my $e = CayleyDickson->new(@{$set->[3]});
   d(c => $c) if DEBUG;
   d(t => $t) if DEBUG;
   $c->cnot($t);
   d(c => $c) if DEBUG;
   d(t => $t) if DEBUG;
   ok(($t-$e)->norm < PRECISION, "$label : expect: $e, calc: $t");
}

diag 'Superpositions states...' if VERBOSE;
$a = Tangle->new(1/sqrt(2),1/sqrt(2));
$b = Tangle->new(1/2,sqrt(3)/2);
$d = Tangle->new(1/sqrt(2),-1/sqrt(2));
ok($a->a - $a->b < PRECISION, "(√½,√½)  = $a");
ok($d->a + $d->b < PRECISION, "(√½,-√½) = $d");

diag 'X-gate ...' if VERBOSE;
$c = Tangle->new(0,1);
$c->y_gate;
$e = CayleyDickson->new(-1,0);
ok(($e-$c)->norm < PRECISION, "(-1,0) : expect: $e, clac: $c");

use constant MEASURES => 10000;
use constant MEASURE_PRECISION => 10 ** -1;

###################################################
diag 'Measurements ...' if VERBOSE;
$a = Tangle->new(1/sqrt(2),1/sqrt(2));
$measures = $a->measures(MEASURES);
ok( 
   ( abs($measures->{0}-.5) < MEASURE_PRECISION and abs($measures->{1}-.5) < MEASURE_PRECISION ),
   sprintf('50%%/50%% measure 0/1 on %s runs of (√½,√½) within %s%%', MEASURES, 100 * MEASURE_PRECISION)
);

$a = Tangle->new(1,0);
$measures = $a->measures(MEASURES);
ok( 
   $measures->{0} == 1,
   sprintf('   100%% measure   0 on %s runs of ( 1, 0) within %s%%', MEASURES, 100 * MEASURE_PRECISION)
);

$a = Tangle->new(0,1);
$measures = $a->measures(10000);
ok($measures->{1} == 1, sprintf('   100%% measure   1 on %s runs of ( 0, 1) within %s%%', MEASURES, 100 * MEASURE_PRECISION));

# 50/50 x 50/50 is 1/4,1/4,1/4,1/4
$a = Tangle->new(1/sqrt(2),1/sqrt(2));
$b = Tangle->new(1/sqrt(2),1/sqrt(2));
$h = $a->tensor($b);
$measures = $h->measures(10000);
ok(
   (abs($measures->{0}-.25) < MEASURE_PRECISION and abs($measures->{1}-.25) < MEASURE_PRECISION and abs($measures->{2}-.25) < MEASURE_PRECISION and abs($measures->{3}-.25) < MEASURE_PRECISION), 
   sprintf('  %s measures of (1/2,1/2,1/2,1/2) has equal 25%% probabilities within %s%%', MEASURES, 100 * MEASURE_PRECISION)
);

###################################################
diag 'Hadamard and Y-gate tests ...' if VERBOSE;

$a = Tangle->new(1/sqrt(2),sqrt(3)/2);
$b = Tangle->new(sqrt(3)/2,1/sqrt(2));
$a->swap;
ok(($a->a - $b->a < PRECISION and $a->b - $b->b < PRECISION), 'SWAP(√½,√3/2)            = (√3/2,√½)');
ok(($a->a - $b->a < PRECISION and $a->b - $b->b < PRECISION), 'SWAP(√½,√3/2)            = (√3/2,√½)');

#$a = Tangle->new(1,0);
$a = Tangle->new(1,0);
$a->hadamard;
$e = CayleyDickson->new(sqrt(1/2),sqrt(1/2));
ok(($a - $e)->norm < PRECISION, "Hadamard(|0>)            = |+>  : expect: $e, calc: $a");

$a = Tangle->new(0,1);
$a->hadamard;
$e = CayleyDickson->new(sqrt(1/2),-sqrt(1/2));
ok(($a - $e)->norm < PRECISION, "Hadamard(|1>)            = |->  : expect: $e, calc: $a");

$a = Tangle->new(1,0);
$a->hadamard;
$a->hadamard;
$e = CayleyDickson->new(1,0);
ok(($a - $e)->norm < PRECISION, "Hadamard(Hadamard(|0>))  = |0>  : expect: $e, calc: $a");

$a = Tangle->new(0,1);
$a->hadamard;
$a->hadamard;
$e = CayleyDickson->new(0,1);
ok(($a - $e)->norm < PRECISION, "Hadamard(Hadamard(|1>))  = |1>  : expect: $e, calc: $a");

###################################################
diag 'Hadamard and X-gate tests ...' if VERBOSE;

$a = Tangle->new(1,0);
$a->x_gate;
$a->hadamard;
$e = CayleyDickson->new(sqrt(1/2),-sqrt(1/2));
ok(($a - $e)->norm < PRECISION, "Hadamard(XGate(|0>))     = ( √½,-√½) =  |-> : expect: $e, calc: $a");

$a = Tangle->new(0,1);
$a->x_gate;
$a->hadamard;
$e = CayleyDickson->new(sqrt(1/2),sqrt(1/2));
ok(($a - $e)->norm < PRECISION, "Hadamard(XGate(|1>))     = ( √½, √½) =  |+> : expect: $e, calc: $a");

$a = Tangle->new(1/sqrt(2),1/sqrt(2));
$a->x_gate;
$a->hadamard;
$e = CayleyDickson->new(1,0);
ok(($a - $e)->norm < PRECISION, "H(X(√½,√½)) = H(X(|+>))  = ( √½, √½) =  |+> : expect: $e, calc: $a");

$a = Tangle->new(-1/sqrt(2),1/sqrt(2));
$a->x_gate;
$a->hadamard;
$e = CayleyDickson->new(0,1);
ok(($a - $e)->norm < PRECISION, "H(X(-√½,√½))= H(X(-|+>)) = (  0, 1 ) =  |1> : expect: $e, calc: $a");

$a = Tangle->new(1/sqrt(2),-1/sqrt(2));
$a->x_gate;
$a->hadamard;
$e = CayleyDickson->new(0,-1);
ok(($a - $e)->norm < PRECISION, "H(X(√½,-√½))= H(X(|->))  = (  0,-1 ) = -|0> : expect: $e, calc: $a");

$a = Tangle->new(0,-1);
$a->x_gate;
$a->hadamard;
$e = CayleyDickson->new(-sqrt(1/2),-sqrt(1/2));
ok(($a - $e)->norm < PRECISION, "H(X(0,-1))= H(X(-|1>))   = (-√½,-√½) = -|+> : expect: $e, calc: $a");

$a = Tangle->new(-1,0);
$a->x_gate;
$a->hadamard;
$e = CayleyDickson->new(-sqrt(1/2),sqrt(1/2));
ok(($a - $e)->norm < PRECISION, "H(X(-1,0))= H(X(-|0>))   = (-√½, √½) = -|-> : expect: $e, calc: $a");

$a = Tangle->new(-1/sqrt(2),-1/sqrt(2));
$a->x_gate;
$a->hadamard;
$e = CayleyDickson->new(-1,0);
ok(($a - $e)->norm < PRECISION, "H(X(-√½,-√½))= H(X(-|+>))= ( -1, 0 ) = -|0> : expect: $e, calc: $a");

###################################################
diag '16 step walk about the complex unit circle: X(H(X(H(X(H(X(H(X(H(X(H(X(H(X(H(|0>)))))))) = |0>' if VERBOSE;

$a = Tangle->new(1,0);
$e = CayleyDickson->new(1,0);
ok(($a-$e)->norm < PRECISION                                      , "Step 0: Starting    = (1,0)      : expect: $e, calc: $a");

$a->x_gate;
$e = CayleyDickson->new(+0,1);
ok(($a-$e)->norm < PRECISION                                      , "Step 1: X(1,0)      = (0,1)      : expect: $e, calc: $a");

$a->hadamard;
$e = CayleyDickson->new(+sqrt(1/2),-sqrt(1/2));
ok(($a-$e)->norm < PRECISION                                      , "Step 2: H(0,1)      = (√½,-√½)   : expect: $e, calc: $a");

$a->x_gate;
$e = CayleyDickson->new(-sqrt(1/2),+sqrt(1/2));
ok(($a-$e)->norm < PRECISION                                      , "Step 3: X(√½,-√½)   = (-√½,√½)   : expect: $e, calc: $a");

$a->hadamard;
$e = CayleyDickson->new(+0,-1);
ok(($a-$e)->norm < PRECISION                                      , "Step 4: H(-√½,√½)   = (0,-1)     : expect: $e, calc: $a");

$a->x_gate;
$e = CayleyDickson->new(-1,0);
ok(($a-$e)->norm < PRECISION                                      , "Step 5: X(0,-1)     = (-1,0)     : expect: $e, calc: $a");

$a->hadamard;
$e = CayleyDickson->new(-sqrt(1/2),-sqrt(1/2));
ok(($a-$e)->norm < PRECISION                                      , "Step 6: H(-1,0)     = (-√½,-√½)  : expect: $e, calc: $a");

$a->x_gate;
$e = CayleyDickson->new(-sqrt(1/2),-sqrt(1/2));
ok(($a-$e)->norm < PRECISION                                      , "Step 7: X(-√½,-√½)  = (-√½,-√½)  : expect: $e, calc: $a");

$a->hadamard;
$e = CayleyDickson->new(-1,0);
ok(($a-$e)->norm < PRECISION                                      , "Step 8: H(-√½,-√½)  = (-1,0)     : expect: $e, calc: $a");

$a->x_gate;
$e = CayleyDickson->new(0,-1);
ok(($a-$e)->norm < PRECISION                                      , "Step 9: X(-1,0)     = (0,-1)     : expect: $e, calc: $a");

$a->hadamard;
$e = CayleyDickson->new(-sqrt(1/2),+sqrt(1/2));
ok(($a-$e)->norm < PRECISION                                      , "Step 10: H(0,-1)    = (-√½,√½)   : expect: $e, calc: $a");

$a->x_gate;
$e = CayleyDickson->new(+sqrt(1/2),-sqrt(1/2));
ok(($a-$e)->norm < PRECISION                                      , "Step 11: X(-√½,√½)  = (√½,-√½)   : expect: $e, calc: $a");

$a->hadamard;
$e = CayleyDickson->new(0,1);
ok(($a-$e)->norm < PRECISION                                      , "Step 12: H(√½,-√½)  = (0,1)      : expect: $e, calc: $a");

$a->x_gate;
$e = CayleyDickson->new(1,0);
ok(($a-$e)->norm < PRECISION                                      , "Step 13: X(0,1)     = (1,0)      : expect: $e, calc: $a");

$a->hadamard;
$e = CayleyDickson->new(+sqrt(1/2),+sqrt(1/2));
ok(($a-$e)->norm < PRECISION                                      , "Step 14: H(1,0)     = (√½,√½)    : expect: $e, calc: $a");

$a->x_gate;
$e = CayleyDickson->new(+sqrt(1/2),+sqrt(1/2));
ok(($a-$e)->norm < PRECISION                                      , "Step 15: X(-√½,√½)  = (√½,-√½)   : expect: $e, calc: $a");

$a->hadamard;
$e = CayleyDickson->new(1,0);
ok(($a-$e)->norm < PRECISION                                      , "Step 16: H(√½,√½)   = (1,0)      : expect: $e, calc: $a");


###################################################
diag '4 Black box gate actions: constant 1/0, identity and negate ...' if VERBOSE;

$c = Tangle->new(1,0);
$t = Tangle->new(1,0);
$r = $c->tensor($t);
$e = CayleyDickson->new(1,0,0,0);
ok(($r-$e)->norm < PRECISION, "Constant-0 |00> = |00> : expect: $e, calc: $r");

$c = Tangle->new(0,1);
$t = Tangle->new(1,0);
$r = $c->tensor($t);
$e = CayleyDickson->new(0,0,1,0);
ok(($r-$e)->norm < PRECISION, "Constant-0 |10> = |10> : expect: $e, calc: $r");

$c = Tangle->new(1,0);
$t = Tangle->new(1,0);
$t->x_gate;
$r = $c->tensor($t);
$e = CayleyDickson->new(0,1,0,0);
ok(($r-$e)->norm < PRECISION, "Constant-1 |00> = |01> : expect: $e, calc: $r");

$c = Tangle->new(0,1);
$t = Tangle->new(1,0);
$t->x_gate;
$r = $c->tensor($t);
$e = CayleyDickson->new(0,0,0,1);
ok(($r-$e)->norm < PRECISION, "Constant-1 |10> = |11> : expect: $e, calc: $r");

$c = Tangle->new(1,0);
$t = Tangle->new(1,0);
$r = $c->cnot($t);
$e = CayleyDickson->new(1,0,0,0);
ok(($r-$e)->norm < PRECISION, "Identity   |00> = |00> : expect: $e, calc: $r");

$c = Tangle->new(0,1);
$t = Tangle->new(1,0);
$r = $c->cnot($t);
$e = CayleyDickson->new(0,0,0,1);
ok(($r-$e)->norm < PRECISION, "Identity   |10> = |11> : expect: $e, calc: $r");

$c = Tangle->new(1,0);
$t = Tangle->new(1,0);
$r = $c->cnot($t);
$t->x_gate;
$e = CayleyDickson->new(0,0,1,0);
ok(($r-$e)->norm < PRECISION, "Negate     |00> = |01> : expect: $e, calc: $r");

$c = Tangle->new(0,1);
$t = Tangle->new(1,0);
$r = $c->cnot($t);
$t->x_gate;
$e = CayleyDickson->new(0,-1,0,0);
ok(($r-$e)->norm < PRECISION, "Negate     |10> = |10> : expect: $e, calc: $r");


1;

__END__
