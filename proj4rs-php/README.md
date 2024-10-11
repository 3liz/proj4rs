# Minimal binding for php

This crate implement a minimal binding for php.

Do not expect good performances for batch processing because of the nature
of php arrays and constante marshalling for points.

If you are seeking for performance, you would rather go with the C bindings
and php [FFI module](https://www.php.net/manual/en/book.ffi.php).
