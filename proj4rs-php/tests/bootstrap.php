<?php

declare(strict_types=1);

require __DIR__ . '/../vendor/autoload.php';

if (!extension_loaded('proj4rs-php')) {
    fwrite(STDERR, "The proj4rs-php extension is not loaded, run the tests with:\n\n"
        . "  php -d extension=../target/debug/libproj4rs_php.so vendor/bin/phpunit\n\n"
        . "or `cargo make php.unittest`\n");
    exit(1);
}
