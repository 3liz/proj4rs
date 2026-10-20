<?php

// Run with:
//
//   PROJ4RS_LIB=../../target/release/libproj4rs_c.so php tests/test_all.php
//

declare(strict_types=1);

use Proj4rs\Proj;
use Proj4rs\Proj4rsException;
use Proj4rs\Transform;

spl_autoload_register(function (string $class): void {
    $prefix = 'Proj4rs\\';
    if (str_starts_with($class, $prefix)) {
        require __DIR__ . '/../src/' . substr($class, strlen($prefix)) . '.php';
    }
});

const LAEA = '+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80';
const X_IN = 15.4213696;
const Y_IN = 47.0766716;
const X_OUT = 4732659.007426266;
const Y_OUT = 2677630.7269610995;

function assert_approx(float $expected, float $actual, float $rel = 1e-6): void
{
    if (abs($expected - $actual) > $rel * abs($expected)) {
        throw new AssertionError("Expected $expected, got $actual");
    }
}

function assert_true(bool $cond, string $msg = 'Assertion failed'): void
{
    if (!$cond) {
        throw new AssertionError($msg);
    }
}

$tests = [
    'proj_properties' => function () {
        $src = new Proj('WGS84');
        assert_true($src->projname() === 'longlat', $src->projname());
        assert_true($src->isLatlong());
        assert_true(!$src->isGeocent());
        assert_true($src->axis() === 'enu', $src->axis());
        assert_true($src->isNormalizedAxis());
        assert_true($src->units() === 'degrees', $src->units());

        $dst = new Proj(LAEA);
        assert_true($dst->projname() === 'laea', $dst->projname());
        assert_true(!$dst->isLatlong());
        assert_true($dst->units() === 'm', $dst->units());
        assert_approx(1.0, $dst->toMeter());
    },

    'invalid_proj' => function () {
        try {
            new Proj('+proj=foobar');
        } catch (Proj4rsException $e) {
            assert_true($e->getMessage() !== '', 'Empty error message');
            return;
        }
        throw new AssertionError('Expected Proj4rsException');
    },

    'transform_scalar' => function () {
        $trans = new Transform(new Proj('WGS84'), new Proj(LAEA));
        [$x, $y] = $trans->transform(X_IN, Y_IN);
        assert_approx(X_OUT, $x);
        assert_approx(Y_OUT, $y);
    },

    'transform_scalar_3d' => function () {
        $trans = new Transform('WGS84', LAEA);
        [$x, $y, $z] = $trans->transform(X_IN, Y_IN, 0.0);
        assert_approx(X_OUT, $x);
        assert_approx(Y_OUT, $y);
        assert_true($z === 0.0, "z = $z");
    },

    'transform_arrays' => function () {
        $trans = new Transform('WGS84', LAEA);
        [$xs, $ys] = $trans->transform([X_IN, X_IN], [Y_IN, Y_IN]);
        assert_true(count($xs) === 2 && count($ys) === 2);
        foreach ([0, 1] as $i) {
            assert_approx(X_OUT, $xs[$i]);
            assert_approx(Y_OUT, $ys[$i]);
        }
    },

    'transform_empty' => function () {
        $trans = new Transform('WGS84', LAEA);
        assert_true($trans->transform([], []) === [[], []]);
        assert_true($trans->transformPoints([]) === []);
    },

    'transform_points' => function () {
        $trans = new Transform('WGS84', LAEA);
        $out = $trans->transformPoints([[X_IN, Y_IN], [X_IN, Y_IN]]);
        assert_true(count($out) === 2);
        foreach ($out as [$x, $y]) {
            assert_approx(X_OUT, $x);
            assert_approx(Y_OUT, $y);
        }

        $out = $trans->transformPoints([[X_IN, Y_IN, 10.0], [X_IN, Y_IN, 20.0]]);
        assert_approx(X_OUT, $out[1][0]);
        assert_approx(Y_OUT, $out[1][1]);
        assert_true($out[0][2] === 10.0 && $out[1][2] === 20.0);
    },

    'transform_roundtrip' => function () {
        $trans = new Transform(LAEA, 'WGS84');
        [$x, $y] = $trans->transformPoint([X_OUT, Y_OUT]);
        assert_approx(X_IN, $x);
        assert_approx(Y_IN, $y);
    },

    'mismatched_arrays' => function () {
        $trans = new Transform('WGS84', LAEA);
        try {
            $trans->transform([1.0, 2.0], [1.0]);
        } catch (InvalidArgumentException) {
            return;
        }
        throw new AssertionError('Expected InvalidArgumentException');
    },
];

$failed = 0;
foreach ($tests as $name => $test) {
    try {
        $test();
        echo "PASS $name\n";
    } catch (Throwable $e) {
        $failed++;
        echo "FAIL $name: {$e->getMessage()}\n";
    }
}

echo sprintf("\n%d passed, %d failed\n", count($tests) - $failed, $failed);
exit($failed === 0 ? 0 : 1);
