<?php

declare(strict_types=1);

namespace Proj4rs\Tests;

use PHPUnit\Framework\Attributes\DataProvider;
use PHPUnit\Framework\TestCase;
use Proj4rs\Proj;
use Proj4rs\Proj4rsException;

final class ProjTest extends TestCase
{
    public function testClassesAreProvidedByTheExtension(): void
    {
        foreach ([Proj::class, \Proj4rs\Transform::class, Proj4rsException::class] as $class) {
            $this->assertTrue((new \ReflectionClass($class))->isInternal(), $class);
        }
        $this->assertTrue((new \ReflectionClass(Proj::class))->isFinal());
    }

    /**
     * @return iterable<string, array{string, string, bool, bool, string, bool, float, string}>
     */
    public static function definitions(): iterable
    {
        yield 'WGS84' => ['WGS84', 'longlat', true, false, 'enu', true, 1.0, 'degrees'];
        yield 'laea' => [
            '+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80',
            'laea', false, false, 'enu', true, 1.0, 'm',
        ];
        yield 'geocent' => ['+proj=geocent +ellps=WGS84', 'geocent', false, true, 'enu', true, 1.0, 'm'];
        yield 'utm us-ft' => [
            '+proj=utm +zone=32 +ellps=GRS80 +units=us-ft',
            'utm', false, false, 'enu', true, 0.3048006096012192, 'us-ft',
        ];
        yield 'longlat neu' => [
            '+proj=longlat +ellps=WGS84 +axis=neu',
            'longlat', true, false, 'neu', false, 1.0, 'degrees',
        ];
    }

    #[DataProvider('definitions')]
    public function testProperties(
        string $defn,
        string $projname,
        bool $isLatlong,
        bool $isGeocent,
        string $axis,
        bool $isNormalizedAxis,
        float $toMeter,
        string $units,
    ): void {
        $proj = new Proj($defn);

        $this->assertSame($projname, $proj->projname());
        $this->assertSame($isLatlong, $proj->isLatlong());
        $this->assertSame($isGeocent, $proj->isGeocent());
        $this->assertSame($axis, $proj->axis());
        $this->assertSame($isNormalizedAxis, $proj->isNormalizedAxis());
        $this->assertEqualsWithDelta($toMeter, $proj->toMeter(), 1e-12);
        $this->assertSame($units, $proj->units());
    }

    /**
     * @return iterable<string, array{string}>
     */
    public static function invalidDefinitions(): iterable
    {
        yield 'empty' => [''];
        yield 'garbage' => ['not a projection'];
        yield 'unknown projection' => ['+proj=foobar'];
    }

    #[DataProvider('invalidDefinitions')]
    public function testInvalidDefinitionThrows(string $defn): void
    {
        $this->expectException(Proj4rsException::class);
        new Proj($defn);
    }

    public function testExceptionIsARuntimeException(): void
    {
        try {
            new Proj('+proj=foobar');
            $this->fail('Expected Proj4rsException');
        } catch (Proj4rsException $e) {
            $this->assertInstanceOf(\RuntimeException::class, $e);
            $this->assertNotSame('', $e->getMessage());
        }
    }

    public function testConstructorRequiresADefinition(): void
    {
        $this->expectException(\ArgumentCountError::class);
        /** @phpstan-ignore arguments.count */
        new Proj();
    }
}
