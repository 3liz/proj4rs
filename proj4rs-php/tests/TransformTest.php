<?php

declare(strict_types=1);

namespace Proj4rs\Tests;

use PHPUnit\Framework\Attributes\DataProvider;
use PHPUnit\Framework\TestCase;
use Proj4rs\Proj;
use Proj4rs\Transform;

final class TransformTest extends TestCase
{
    private const LAEA = '+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80';

    private const LON = 15.4213696;
    private const LAT = 47.0766716;
    private const X = 4732659.007426266;
    private const Y = 2677630.7269610995;

    private const DELTA = 1e-6;

    private Transform $trans;

    protected function setUp(): void
    {
        $this->trans = new Transform('WGS84', self::LAEA);
    }

    public function testConstructFromProjObjects(): void
    {
        $trans = new Transform(new Proj('WGS84'), new Proj(self::LAEA));

        [$x, $y] = $trans->transform(self::LON, self::LAT);
        $this->assertEqualsWithDelta(self::X, $x, self::DELTA);
        $this->assertEqualsWithDelta(self::Y, $y, self::DELTA);
    }

    public function testSourceAndDestination(): void
    {
        $this->assertInstanceOf(Proj::class, $this->trans->source());
        $this->assertSame('longlat', $this->trans->source()->projname());
        $this->assertSame('laea', $this->trans->destination()->projname());
    }

    public function testTransformScalar(): void
    {
        $res = $this->trans->transform(self::LON, self::LAT);

        $this->assertCount(2, $res);
        $this->assertEqualsWithDelta(self::X, $res[0], self::DELTA);
        $this->assertEqualsWithDelta(self::Y, $res[1], self::DELTA);
    }

    public function testTransformScalar3d(): void
    {
        [$x, $y, $z] = $this->trans->transform(self::LON, self::LAT, 10.0);

        $this->assertEqualsWithDelta(self::X, $x, self::DELTA);
        $this->assertEqualsWithDelta(self::Y, $y, self::DELTA);
        $this->assertSame(10.0, $z);
    }

    public function testTransformIntegers(): void
    {
        [$x, $y] = $this->trans->transform(15, 47);

        $this->assertIsFloat($x);
        $this->assertIsFloat($y);
    }

    public function testTransformNullZ(): void
    {
        $this->assertCount(2, $this->trans->transform(self::LON, self::LAT, null));
    }

    public function testTransformArrays(): void
    {
        [$xs, $ys] = $this->trans->transform([self::LON, self::LON, self::LON], [self::LAT, self::LAT, self::LAT]);

        $this->assertCount(3, $xs);
        $this->assertCount(3, $ys);
        foreach ([0, 1, 2] as $i) {
            $this->assertEqualsWithDelta(self::X, $xs[$i], self::DELTA);
            $this->assertEqualsWithDelta(self::Y, $ys[$i], self::DELTA);
        }
    }

    public function testTransformArrays3d(): void
    {
        [$xs, $ys, $zs] = $this->trans->transform([self::LON, self::LON], [self::LAT, self::LAT], [1.0, 2.0]);

        $this->assertEqualsWithDelta(self::X, $xs[1], self::DELTA);
        $this->assertEqualsWithDelta(self::Y, $ys[1], self::DELTA);
        $this->assertSame([1.0, 2.0], $zs);
    }

    public function testTransformArraysIgnoresKeys(): void
    {
        [$xs, $ys] = $this->trans->transform(['a' => self::LON], ['b' => self::LAT]);

        $this->assertSame([0], array_keys($xs));
        $this->assertEqualsWithDelta(self::X, $xs[0], self::DELTA);
        $this->assertEqualsWithDelta(self::Y, $ys[0], self::DELTA);
    }

    public function testTransformEmpty(): void
    {
        $this->assertSame([[], []], $this->trans->transform([], []));
        $this->assertSame([[], [], []], $this->trans->transform([], [], []));
        $this->assertSame([], $this->trans->transformPoints([]));
    }

    public function testTransformWithoutConversion(): void
    {
        [$x, $y] = $this->trans->transform(deg2rad(self::LON), deg2rad(self::LAT), null, false);

        $this->assertEqualsWithDelta(self::X, $x, self::DELTA);
        $this->assertEqualsWithDelta(self::Y, $y, self::DELTA);
    }

    public function testTransformNamedArguments(): void
    {
        [$x, $y] = $this->trans->transform(x: deg2rad(self::LON), y: deg2rad(self::LAT), convert: false);

        $this->assertEqualsWithDelta(self::X, $x, self::DELTA);
        $this->assertEqualsWithDelta(self::Y, $y, self::DELTA);
    }

    public function testConvertDefaultsToTrue(): void
    {
        $params = (new \ReflectionMethod(Transform::class, 'transform'))->getParameters();

        $this->assertSame('convert', $params[3]->getName());
        $this->assertTrue($params[3]->getDefaultValue());
    }

    public function testTransformPoint(): void
    {
        [$x, $y] = $this->trans->transformPoint([self::LON, self::LAT]);

        $this->assertEqualsWithDelta(self::X, $x, self::DELTA);
        $this->assertEqualsWithDelta(self::Y, $y, self::DELTA);
    }

    public function testTransformPoints(): void
    {
        $points = $this->trans->transformPoints([[self::LON, self::LAT], [self::LON, self::LAT]]);

        $this->assertCount(2, $points);
        foreach ($points as [$x, $y]) {
            $this->assertEqualsWithDelta(self::X, $x, self::DELTA);
            $this->assertEqualsWithDelta(self::Y, $y, self::DELTA);
        }
    }

    public function testTransformPoints3d(): void
    {
        $points = $this->trans->transformPoints([[self::LON, self::LAT, 10.0], [self::LON, self::LAT, 20.0]]);

        $this->assertCount(3, $points[0]);
        $this->assertSame(10.0, $points[0][2]);
        $this->assertSame(20.0, $points[1][2]);
    }

    public function testRoundTrip(): void
    {
        $inverse = new Transform(self::LAEA, 'WGS84');

        [$lon, $lat] = $inverse->transform(self::X, self::Y);
        $this->assertEqualsWithDelta(self::LON, $lon, 1e-9);
        $this->assertEqualsWithDelta(self::LAT, $lat, 1e-9);
    }

    public function testGeocentric(): void
    {
        $trans = new Transform('WGS84', '+proj=geocent +ellps=WGS84');

        [$x, $y, $z] = $trans->transform(0.0, 0.0, 0.0);
        $this->assertEqualsWithDelta(6378137.0, $x, self::DELTA);
        $this->assertEqualsWithDelta(0.0, $y, self::DELTA);
        $this->assertEqualsWithDelta(0.0, $z, self::DELTA);
    }

    public function testFailedPointsAreNan(): void
    {
        $trans = new Transform('WGS84', '+proj=merc +ellps=WGS84');

        [$valid, $invalid] = $trans->transformPoints([[0.0, 0.0], [0.0, 95.0]]);
        $this->assertFalse(is_nan($valid[0]));
        $this->assertNan($invalid[0]);
        $this->assertNan($invalid[1]);
    }

    /**
     * @return iterable<string, array{\Closure(Transform): mixed, string}>
     */
    public static function invalidArguments(): iterable
    {
        yield 'invalid source' => [
            fn () => new Transform(42, 'WGS84'),
            'Expecting a Proj4rs\Proj object or a string',
        ];
        yield 'invalid destination' => [
            fn () => new Transform('WGS84', new \stdClass()),
            'Expecting a Proj4rs\Proj object or a string',
        ];
        yield 'mixed scalar and array' => [
            fn (Transform $t) => $t->transform([1.0], 2.0),
            'Coordinates must be all scalars or all arrays',
        ];
        yield 'mismatched lengths' => [
            fn (Transform $t) => $t->transform([1.0, 2.0], [1.0]),
            'Coordinate arrays must have the same length',
        ];
        yield 'non numeric' => [
            fn (Transform $t) => $t->transform(['a'], [1.0]),
            'Expecting numeric coordinates',
        ];
        yield 'point dimension' => [
            fn (Transform $t) => $t->transformPoint([1.0]),
            'Points must have 2 or 3 coordinates',
        ];
        yield 'points dimension' => [
            fn (Transform $t) => $t->transformPoints([[1.0, 2.0, 3.0, 4.0]]),
            'Points must have 2 or 3 coordinates',
        ];
        yield 'mixed dimensions' => [
            fn (Transform $t) => $t->transformPoints([[1.0, 2.0], [1.0, 2.0, 3.0]]),
            'All points must have the same dimension',
        ];
        yield 'point not an array' => [
            fn (Transform $t) => $t->transformPoints([1.0, 2.0]),
            'Points must be arrays',
        ];
    }

    #[DataProvider('invalidArguments')]
    public function testInvalidArgumentsThrow(\Closure $call, string $message): void
    {
        $this->expectException(\InvalidArgumentException::class);
        $this->expectExceptionMessage($message);
        $call($this->trans);
    }

    public function testInvalidProjectionStringThrows(): void
    {
        $this->expectException(\Proj4rs\Proj4rsException::class);
        new Transform('WGS84', '+proj=foobar');
    }
}
