<?php

declare(strict_types=1);

namespace Proj4rs;

use FFI;
use FFI\CData;

/**
 * Coordinate transformation between two projections
 *
 * Coordinates that fail to be projected are set to NAN.
 */
final class Transform
{
    private const SIZEOF_DOUBLE = 8;

    private Proj $src;
    private Proj $dst;

    public function __construct(Proj|string $src, Proj|string $dst)
    {
        $this->src = is_string($src) ? new Proj($src) : $src;
        $this->dst = is_string($dst) ? new Proj($dst) : $dst;
    }

    public function source(): Proj
    {
        return $this->src;
    }

    public function destination(): Proj
    {
        return $this->dst;
    }

    /**
     * Transform coordinates given as separate scalars or arrays of values.
     *
     * If `$convert` is true, latlong coordinates are assumed to be in degrees.
     *
     * Returns `[x, y]` or `[x, y, z]`, as floats if the inputs are scalars,
     * as lists of floats otherwise.
     *
     * @param float|float[] $x
     * @param float|float[] $y
     * @param float|float[]|null $z
     * @return array{0: float|float[], 1: float|float[], 2?: float|float[]}
     *
     * @throws Proj4rsException
     */
    public function transform(
        float|array $x,
        float|array $y,
        float|array|null $z = null,
        bool $convert = true,
    ): array {
        $scalar = !is_array($x);
        if ($scalar !== !is_array($y) || ($z !== null && $scalar !== !is_array($z))) {
            throw new \InvalidArgumentException('Coordinates must be all scalars or all arrays');
        }

        $xs = $scalar ? [$x] : array_values($x);
        $ys = $scalar ? [$y] : array_values($y);
        $zs = $z === null ? null : ($scalar ? [$z] : array_values($z));

        $len = count($xs);
        if (count($ys) !== $len || ($zs !== null && count($zs) !== $len)) {
            throw new \InvalidArgumentException('Coordinate arrays must have the same length');
        }
        if ($len === 0) {
            return $zs === null ? [[], []] : [[], [], []];
        }

        $bx = self::toBuffer($xs);
        $by = self::toBuffer($ys);
        $bz = $zs === null ? null : self::toBuffer($zs);

        $this->call(
            FFI::addr($bx[0]),
            FFI::addr($by[0]),
            $bz === null ? null : FFI::addr($bz[0]),
            $len,
            self::SIZEOF_DOUBLE,
            $convert,
        );

        if ($scalar) {
            return $bz === null ? [$bx[0], $by[0]] : [$bx[0], $by[0], $bz[0]];
        }
        $rx = self::fromBuffer($bx, $len);
        $ry = self::fromBuffer($by, $len);
        return $bz === null ? [$rx, $ry] : [$rx, $ry, self::fromBuffer($bz, $len)];
    }

    /**
     * Transform a single point given as `[x, y]` or `[x, y, z]`
     *
     * @param float[] $point
     * @return float[]
     *
     * @throws Proj4rsException
     */
    public function transformPoint(array $point, bool $convert = true): array
    {
        return $this->transformPoints([$point], $convert)[0];
    }

    /**
     * Transform a list of points, each given as `[x, y]` or `[x, y, z]`.
     *
     * All points must have the same dimension.
     *
     * @param array<float[]> $points
     * @return array<float[]>
     *
     * @throws Proj4rsException
     */
    public function transformPoints(array $points, bool $convert = true): array
    {
        $len = count($points);
        if ($len === 0) {
            return [];
        }

        $points = array_values($points);
        $dim = count($points[0]);
        if ($dim !== 2 && $dim !== 3) {
            throw new \InvalidArgumentException('Points must have 2 or 3 coordinates');
        }

        // Interleaved buffer: x0, y0, [z0,] x1, y1, [z1,] ...
        $buf = FFI::new("double[" . ($len * $dim) . "]");
        $i = 0;
        foreach ($points as $p) {
            if (count($p) !== $dim) {
                throw new \InvalidArgumentException('All points must have the same dimension');
            }
            foreach ($p as $v) {
                $buf[$i++] = (float) $v;
            }
        }

        $this->call(
            FFI::addr($buf[0]),
            FFI::addr($buf[1]),
            $dim === 3 ? FFI::addr($buf[2]) : null,
            $len,
            $dim * self::SIZEOF_DOUBLE,
            $convert,
        );

        $out = [];
        for ($k = 0, $i = 0; $k < $len; $k++, $i += $dim) {
            $out[] = $dim === 3
                ? [$buf[$i], $buf[$i + 1], $buf[$i + 2]]
                : [$buf[$i], $buf[$i + 1]];
        }
        return $out;
    }

    private function call(CData $x, CData $y, ?CData $z, int $len, int $stride, bool $convert): void
    {
        $ffi = Library::ffi();
        $res = $ffi->proj4rs_transform(
            $this->src->cdata(),
            $this->dst->cdata(),
            $x,
            $y,
            $z,
            $len,
            $stride,
            $convert,
        );
        if ($res !== 1) {
            throw new Proj4rsException(Library::lastError());
        }
    }

    /** @param float[] $values */
    private static function toBuffer(array $values): CData
    {
        $buf = FFI::new("double[" . count($values) . "]");
        foreach ($values as $i => $v) {
            $buf[$i] = (float) $v;
        }
        return $buf;
    }

    /** @return float[] */
    private static function fromBuffer(CData $buf, int $len): array
    {
        $out = [];
        for ($i = 0; $i < $len; $i++) {
            $out[] = $buf[$i];
        }
        return $out;
    }
}
