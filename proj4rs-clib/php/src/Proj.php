<?php

declare(strict_types=1);

namespace Proj4rs;

use FFI;
use FFI\CData;

/**
 * A projection definition
 */
final class Proj
{
    private CData $cdata;

    /**
     * `$defn` may be:
     *
     * - A projstring
     * - "WGS84" - equivalent to the projstring "+proj=longlat +ellps=WGS84"
     * - An EPSG code (i.e "EPSG:3035") if the library is built
     *   with the `crs-definitions` feature
     *
     * @throws Proj4rsException
     */
    public function __construct(string $defn)
    {
        $ffi = Library::ffi();
        $cdata = $ffi->proj4rs_proj_new($defn);
        if ($cdata === null || FFI::isNull($cdata)) {
            throw new Proj4rsException(Library::lastError());
        }
        $this->cdata = $cdata;
    }

    public function __destruct()
    {
        Library::ffi()->proj4rs_proj_delete($this->cdata);
    }

    public function __clone()
    {
        throw new \LogicException('Proj objects cannot be cloned');
    }

    /** @internal */
    public function cdata(): CData
    {
        return $this->cdata;
    }

    public function projname(): string
    {
        return Library::ffi()->proj4rs_proj_projname($this->cdata);
    }

    public function isLatlong(): bool
    {
        return Library::ffi()->proj4rs_proj_is_latlong($this->cdata);
    }

    public function isGeocent(): bool
    {
        return Library::ffi()->proj4rs_proj_is_geocent($this->cdata);
    }

    /**
     * Return the axis directions as a 3 chars string (i.e "enu")
     *
     * - x: 'e' (East) or 'w' (West)
     * - y: 'n' (North) or 's' (South)
     * - z: 'u' (Up) or 'd' (Down)
     */
    public function axis(): string
    {
        return FFI::string(Library::ffi()->proj4rs_proj_axis($this->cdata), 3);
    }

    public function isNormalizedAxis(): bool
    {
        return Library::ffi()->proj4rs_proj_is_normalized_axis($this->cdata);
    }

    public function toMeter(): float
    {
        return Library::ffi()->proj4rs_proj_to_meter($this->cdata);
    }

    /**
     * Return units of the projection (i.e "degrees", "m", "km", ...)
     */
    public function units(): string
    {
        return Library::ffi()->proj4rs_proj_units($this->cdata);
    }
}
