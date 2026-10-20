<?php

declare(strict_types=1);

namespace Proj4rs;

use FFI;

/**
 * Holds the FFI handle to the proj4rs C library.
 *
 * The library is resolved in this order:
 *
 * 1. The path given to `Library::load()`
 * 2. The `PROJ4RS` FFI scope if `proj4rs.h` was preloaded with `FFI::load()`
 *    from an opcache.preload script
 * 3. The `PROJ4RS_LIB` environment variable
 * 4. The library bundled in the package `lib/` directory
 * 5. The platform default library name, resolved by the dynamic linker
 */
final class Library
{
    public const SCOPE = 'PROJ4RS';

    private static ?FFI $ffi = null;

    /**
     * Explicitly load the C library from `$path`
     */
    public static function load(?string $path = null): FFI
    {
        $path = $path ?: (getenv('PROJ4RS_LIB') ?: (self::bundledLib() ?? self::defaultLibName()));
        try {
            self::$ffi = FFI::cdef(self::cdef(), $path);
        } catch (FFI\Exception $e) {
            throw new Proj4rsException("Failed to load proj4rs library '$path': {$e->getMessage()}", 0, $e);
        }
        return self::$ffi;
    }

    public static function ffi(): FFI
    {
        if (self::$ffi === null) {
            try {
                self::$ffi = FFI::scope(self::SCOPE);
            } catch (FFI\Exception) {
                self::load();
            }
        }
        return self::$ffi;
    }

    /**
     * Return the last error message from the library
     */
    public static function lastError(): string
    {
        // FFI converts `const char *` return values to PHP strings
        return self::ffi()->proj4rs_last_error();
    }

    private static function cdef(): string
    {
        $header = file_get_contents(__DIR__ . '/../proj4rs.h');
        // FFI::cdef() does not support preprocessor directives
        return preg_replace('/^\s*#.*$/m', '', $header);
    }

    private static function bundledLib(): ?string
    {
        $path = __DIR__ . '/../lib/' . self::defaultLibName();
        return is_file($path) ? $path : null;
    }

    private static function defaultLibName(): string
    {
        return match (PHP_OS_FAMILY) {
            'Windows' => 'proj4rs_c.dll',
            'Darwin' => 'libproj4rs_c.dylib',
            default => 'libproj4rs_c.so',
        };
    }
}
