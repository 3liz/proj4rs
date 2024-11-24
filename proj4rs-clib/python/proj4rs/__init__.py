from array import array
from collections import abc
from typing import Any, Tuple, TypeVar, Union, overload

from ._proj4rs import ffi, lib


class Proj:

    def __init__(self, defn: str):
        _defn = ffi.new("char[]", defn.encode())
        self._cdata = lib.proj4rs_proj_new(_defn)

    def __del__(self):
        lib.proj4rs_proj_delete(self._cdata)

    @property
    def projname(self) -> str:
        _rv = lib.proj4rs_proj_projname(self._cdata)
        return ffi.string(_rv).decode()

    @property
    def is_latlong(self) -> bool:
        return lib.proj4rs_proj_is_latlong(self._cdata)

    @property
    def is_geocent(self) -> bool:
        return lib.proj4rs_proj_is_geocent(self._cdata)

    @property
    def axis(self) -> bytes:
        _rv = lib.proj4rs_proj_axis(self._cdata)
        return bytes(ffi.cast("uint8_t[3]", _rv))

    @property
    def is_normalized_axis(self) -> bool:
        return lib.proj4rs_proj_is_normalized_axis(self._cdata)

    @property
    def to_meter(self) -> float:
        return lib.proj4rs_proj_to_meter(self._cdata)

    @property
    def units(self) -> str:
        _rv = lib.proj4rs_proj_units(self._cdata)
        return ffi.string(_rv).decode()


def _scalar_to_buffer(x):
    return array("d", (float(x),))


def _copy_buffer(x, inplace):
    match x:
        case array():
            if not inplace or x.typecode != 'd':
                x = array("d", x)
        case memoryview():
            # Ensure 1 dimensional data
            if x.ndim != 1:
                raise ValueError("Expecting 1 dimensional array")
            if not inplace or x.format != 'd':
                x = array("d", x)
        case abc.Sequence():
            x = array("d", x)
        case _:
            raise ValueError("Invalid buffer type")
    return x


SIZEOF_DOUBLE = ffi.sizeof("double")


T = TypeVar('T')


class Transform:

    def __init__(self, src: Proj | str, dst: Proj | str):
        self._from = Proj(src) if isinstance(src, str) else src
        self._to = Proj(dst) if isinstance(dst, str) else dst

    @property
    def source(self) -> Proj:
        return self._from

    @property
    def destination(self) -> Proj:
        return self._to

    @overload
    def transform(
        self,
        x: abc.Buffer,
        *,
        convert: bool = True,
        inplace: bool = False,
    ) -> Union[
        Tuple[abc.Buffer, abc.Buffer],
        Tuple[abc.Buffer, abc.Buffer, abc.Buffer],
    ]: ...

    @overload
    def transform(
        self,
        x: float | int,
        y: float | int,
        *,
        convert: bool = True,
        inplace: bool = False,
    ) -> Union[
        Tuple[float, float],
    ]: ...

    @overload
    def transform(
        self,
        x: float | int,
        y: float | int,
        z: float | int,
        *,
        convert: bool = True,
        inplace: bool = False,
    ) -> Union[
        Tuple[float, float, float],
    ]: ...

    @overload
    def transform(
        self,
        x: list | tuple,
        y: list | tuple,
        *,
        convert: bool = True,
        inplace: bool = False,
    ) -> Union[
        Tuple[array, array],
    ]: ...

    @overload
    def transform(
        self,
        x: list | tuple,
        y: list | tuple,
        z: list | tuple,
        *,
        convert: bool = True,
        inplace: bool = False,
    ) -> Union[
        Tuple[array, array, array],
    ]: ...

    @overload
    def transform(
        self,
        x: abc.Buffer,
        y: abc.Buffer,
        *,
        convert: bool = True,
        inplace: bool = False,
    ) -> Union[
        Tuple[abc.Buffer, abc.Buffer],
    ]: ...

    @overload
    def transform(
        self,
        x: abc.Buffer,
        y: abc.Buffer,
        z: abc.Buffer,
        *,
        convert: bool = True,
        inplace: bool = False,
    ) -> Union[
        Tuple[abc.Buffer, abc.Buffer, abc.Buffer],
    ]: ...

    def transform(
        self,
        x: T,
        y: T | None = None,
        z: T | None = None,
        *,
        convert: bool = True,
        inplace: bool = False,
    ) -> Union[
        Tuple[Any, Any],
        Tuple[Any, Any, Any],
    ]:
        """ Transform coordinates

            Parameters
            ----------

            xx: scalar or sequence, input x coordinate(s)
            yy: scalar or sequence, optional, input x coordinate(s)
            zz: scalar or sequence, optional, input x coordinate(s)

            convert: if true, assume that coordinates are  in degrees and the transformation
                     will convert data accordingly
            inplace: if true, convert data inplace if the input data implement the Buffer
                     protocol. The buffer must be writable

            Returns
            -------
            A tuple of buffer objects in the case the input is a Sequence,
            a tuple of float otherwise.
        """
        match (x, y, z):
            case (abc.Buffer(), None, None):
                scalar = False
                m = memoryview(x)
                if m.ndim != 2:
                    raise ValueError("Expecting two-dimensional buffer")
                if m.shape is None:
                    raise ValueError("Invalid buffer shape (None)")
                size, dim = m.shape
                if dim != 2 and dim != 3:
                    raise ValueError(f"Expecting geometry dimensions of 2 or 3, found {dim}")
                # Flatten buffer
                flatten = m.cast('b').cast(m.format)  # type: ignore [call-overload]
                if not inplace or m.format != 'd' or not m.c_contiguous:
                    _x = array('d', flatten[0::dim])
                    _y = array('d', flatten[1::dim])
                    _z = array('d', flatten[2::dim]) if dim > 2 else None
                else:
                    _x = flatten[0::dim]
                    _y = flatten[1::dim]
                    _z = flatten[2::dim] if dim > 2 else None

                stride = dim * SIZEOF_DOUBLE

            case (abc.Sequence(), abc.Sequence(), _):
                scalar = False
                if len(y) != len(x) and (not z or len(z) != len(x)):  # type: ignore [arg-type]
                    raise ValueError("Arrays must have the same length")
                _x = _copy_buffer(x, inplace)
                _y = _copy_buffer(y, inplace)
                _z = _copy_buffer(z, inplace) if z else None
                size = len(_x)
                stride = SIZEOF_DOUBLE
            case (abc.Buffer(), abc.Buffer(), _):
                scalar = False
                mx = memoryview(x)
                my = memoryview(y)
                mz = memoryview(z) if z else None    # type: ignore [arg-type]
                if len(my) != len(mx) and (not mz or len(mz) != len(mx)):  #
                    raise ValueError("Buffers must have same length")
                _x = _copy_buffer(mx, inplace)
                _y = _copy_buffer(my, inplace)
                _z = _copy_buffer(mz, inplace) if mz else None
                size = len(_x)
                stride = SIZEOF_DOUBLE
            case _:
                scalar = True
                _x = _scalar_to_buffer(x)
                _y = _scalar_to_buffer(y)
                _z = _scalar_to_buffer(z) if z else None
                size = 1
                stride = SIZEOF_DOUBLE

        _t = "double[]"

        _xx = ffi.from_buffer(_t, _x, require_writable=True)
        _yy = ffi.from_buffer(_t, _y, require_writable=True)
        _zz = ffi.from_buffer(_t, _z, require_writable=True) if z else ffi.NULL
        res = lib.proj4rs_transform(
            self._from._cdata,
            self._to._cdata,
            _xx,
            _yy,
            _zz,
            size,
            stride,
            convert,
        )
        if res != 1:
            error = lib.proj4rs_last_error()
            raise RuntimeError(ffi.string(error).decode())

        if scalar:
            return (_x[0], _y[0], _z[0]) if _z else (_x[0], _y[0])
        else:
            return (_x, _y, _z) if _z else (_x, _y)
