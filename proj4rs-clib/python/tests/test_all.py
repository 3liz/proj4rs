import proj4rs
import pytest


def test_transform_sequence():

    src = proj4rs.Proj("WGS84")
    dst = proj4rs.Proj("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80")

    print("Src:", src.projname)
    print("Dst:", dst.projname)

    x = [15.4213696]
    y = [47.0766716]

    trans = proj4rs.Transform(src, dst)

    x, y = trans.transform(x, y)

    print("x =", x)  # Should be 4732659.007426
    print("y =", y)  # Should be 2677630.726961

    assert x[0] == pytest.approx(4732659.007426266, 1e-6)
    assert y[0] == pytest.approx(2677630.7269610995, 1e-6)



def test_transform_scalar():

    src = proj4rs.Proj("WGS84")
    dst = proj4rs.Proj("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80")

    print("Src:", src.projname)
    print("Dst:", dst.projname)

    x = 15.4213696
    y = 47.0766716

    print("Transform")
    trans = proj4rs.Transform(src, dst)

    x, y = trans.transform(x, y)

    print("x =", x)  # Should be 4732659.007426
    print("y =", y)  # Should be 2677630.726961

    assert x == pytest.approx(4732659.007426266, 1e-6)
    assert y == pytest.approx(2677630.7269610995, 1e-6)


def test_transform_buffer_inplace():
    from array import array

    src = proj4rs.Proj("WGS84")
    dst = proj4rs.Proj("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80")

    print("Src:", src.projname)
    print("Dst:", dst.projname)

    x = array('d', [15.4213696])
    y = array('d', [47.0766716])

    trans = proj4rs.Transform(src, dst)

    xx, yy = trans.transform(x, y, inplace=True)

    print("x =", x)  # Should be 4732659.007426
    print("y =", y)  # Should be 2677630.726961

    assert xx is x
    assert yy is y

    assert xx[0] == pytest.approx(4732659.007426266, 1e-6)
    assert yy[0] == pytest.approx(2677630.7269610995, 1e-6)


def test_transform_invalid_buffer():
    from array import array

    src = proj4rs.Proj("WGS84")
    dst = proj4rs.Proj("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80")

    print("Src:", src.projname)
    print("Dst:", dst.projname)

    x = array('d', [15.4213696])

    trans = proj4rs.Transform(src, dst)

    with pytest.raises(ValueError, match="Expecting two-dimensional buffer"):
        trans.transform(x, inplace=True)


def test_transform_buffer_2d():
    from array import array

    x = array('d', [15.4213696, 47.0766716])

    # Reshape to a two dimensionnal array
    m = memoryview(x).cast('b').cast('d', shape=(1,2))
    print("* shape =", m.shape, "ndim", m.ndim)

    transform = proj4rs.Transform(
        "WGS84",
        "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80",
    ).transform

    xx, yy = transform(m, inplace=True)

    print("xx =", list(xx))
    print("yy =", list(yy))

    assert x[0] == pytest.approx(4732659.007426266, 1e-6)
    assert x[1] == pytest.approx(2677630.7269610995, 1e-6)



