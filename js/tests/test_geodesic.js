import { Proj } from '../proj4.js';

// Geodesic: checks the `geodesic` namespace against the reference values
// used in the proj4rs-geodesic crate (from GeographicLib).
function test_geodesic() {
    const { Geodesic, Polygon, Mask, Flags } = Proj.geodesic;

    let failed = 0;
    let total = 0;

    function check(name, ok) {
        total += 1;
        if (!ok) {
            failed += 1;
        }
        console.log(`${ok ? 'ok' : 'FAIL'} ${name}`);
    }

    function close(name, actual, expected, tol) {
        const diff = Math.abs(actual - expected);
        check(`${name}: ${actual} (expected ${expected}, diff ${diff.toExponential(2)})`, diff <= tol);
    }

    const g = Geodesic.wgs84();
    close('wgs84 a', g.a, 6378137.0, 0);
    close('wgs84 f', g.f, 1.0 / 298.257223563, 0);

    const g2 = new Geodesic(6378145.0, 1.0 / 298.25);
    close('custom a', g2.a, 6378145.0, 0);
    close('custom f', g2.f, 0.003352891869237217, 1e-18);

    // Inverse: JFK to Singapore Changi Airport
    let inv = g.inverse(40.64, -73.78, 1.36, 103.99);
    close('inverse JFK-SIN s12', inv.s12, 15347512.94051294, 1e-6);
    close('inverse JFK-SIN azi1', inv.azi1, 3.3057734780176125, 1e-12);
    close('inverse JFK-SIN azi2', inv.azi2, 177.48784020815515, 1e-12);

    // Inverse: Berkeley to Port Moresby
    inv = g.inverse(37.87622, -122.23558, -9.4047, 147.1597);
    close('inverse Berkeley-Port Moresby s12', inv.s12, 10700471.955233702, 1e-6);
    close('inverse Berkeley-Port Moresby azi1', inv.azi1, -96.91639942294974, 1e-12);
    close('inverse Berkeley-Port Moresby azi2', inv.azi2, -127.32548874543627, 1e-12);

    // Inverse along the equator
    inv = g.inverse(0.0, 0.0, 0.0, 10.0);
    close('inverse equator s12', inv.s12, 1113194.9079327357, 1e-5);
    close('inverse equator azi1', inv.azi1, 90.0, 0);

    // Direct: 10000 km NE of JFK
    const dir = g.direct(40.64, -73.78, 45.0, 10e6);
    close('direct lat2', dir.lat2, 32.621100463725796, 1e-12);
    close('direct lon2', dir.lon2, 49.05248709295982, 1e-12);
    close('direct azi2', dir.azi2, 140.40598587680074, 1e-12);

    // General direct problem, round trip through the inverse
    const pos = g.genDirect(40.64, -73.78, 45.0, Flags.None, 10e6, Mask.All);
    close('genDirect lat2', pos.lat2, dir.lat2, 1e-12);
    close('genDirect s12', pos.s12, 10e6, 1e-6);
    check('genDirect computes all quantities',
        [pos.a12, pos.m12, pos.scale12, pos.scale21, pos.area12].every(Number.isFinite));
    const arc = g.genDirect(40.64, -73.78, 45.0, Flags.ArcMode, pos.a12, Mask.Latitude | Mask.Longitude);
    close('genDirect arc mode lat2', arc.lat2, dir.lat2, 1e-9);
    close('genDirect arc mode lon2', arc.lon2, dir.lon2, 1e-9);

    // General inverse problem: masked out quantities are NaN
    let ginv = g.genInverse(40.64, -73.78, 1.36, 103.99, Mask.Distance);
    close('genInverse s12', ginv.s12, 15347512.94051294, 1e-6);
    check('genInverse masked out quantities are NaN',
        [ginv.m12, ginv.scale12, ginv.scale21, ginv.area12].every(Number.isNaN));
    ginv = g.genInverse(40.64, -73.78, 1.36, 103.99, Mask.All);
    check('genInverse Mask.All computes all quantities',
        [ginv.a12, ginv.s12, ginv.m12, ginv.scale12, ginv.scale21, ginv.area12].every(Number.isFinite));

    // Geodesic line: way points between JFK and Singapore Changi Airport
    const line = g.inverseLine(40.64, -73.78, 1.36, 103.99, Mask.None);
    close('inverseLine s13', line.s13, 15347512.94051294, 1e-6);
    const end = line.position(line.s13);
    close('inverseLine end lat', end.lat2, 1.36, 1e-12);
    close('inverseLine end lon', end.lon2, 103.99, 1e-12);
    const mid = line.genPosition(Flags.ArcMode, line.a13 * 0.5, Mask.Latitude | Mask.Longitude);
    const half = g.inverse(40.64, -73.78, mid.lat2, mid.lon2);
    close('inverseLine midpoint distance', half.s12, line.s13 * 0.5, 1e-6);

    const dline = g.directLine(40.64, -73.78, 45.0, 10e6, Mask.Latitude | Mask.Longitude);
    close('directLine s13', dline.s13, 10e6, 1e-6);
    close('directLine end lat', dline.position(dline.s13).lat2, dir.lat2, 1e-12);

    // Polygon: triangle (0°N,0°E), (0°N,90°E), (90°N,0°E)
    const p = new Polygon(g, false);
    p.addPoint(0, 0);
    p.addPoint(0, 90);
    const tested = p.testPoint(90, 0, false, true);
    p.addPoint(90, 0);
    let r = p.compute(false, true);
    check(`polygon num: ${p.num}`, p.num === 3 && r.num === 3);
    close('polygon perimeter', r.perimeter, 30022685.630020067, 1e-6);
    close('polygon area', r.area, 63758202715511.055, 1);
    close('polygon testPoint area', tested.area, r.area, 1e-3);
    p.clear();
    check('polygon clear', p.num === 0);

    // Polyline: only the length is computed
    const pl = g.polygon(true);
    pl.addPoint(0, 0);
    pl.addEdge(90, 1113194.9079327357);
    r = pl.compute(false, true);
    close('polyline length', r.perimeter, 1113194.9079327357, 1e-6);
    check('polyline area is NaN', Number.isNaN(r.area));

    // Polygon area: Antarctica
    const lats = [-72.9, -71.9, -74.9, -74.3, -77.5, -77.4, -71.7, -65.9, -65.7,
        -66.6, -66.9, -69.8, -70.0, -71.0, -77.3, -77.9, -74.7];
    const lons = [-74., -102., -102., -131., -163., 163., 172., 140., 113.,
        88., 59., 25., -4., -14., -33., -46., -61.];
    r = g.polygonArea(new Float64Array(lats), new Float64Array(lons));
    close('Antarctica area', r.area, 13376856682207.375, 1);
    close('Antarctica perimeter', r.perimeter, 14710425.40697364, 1e-6);

    let thrown = false;
    try {
        g.polygonArea(new Float64Array([0, 0]), new Float64Array([0]));
    } catch (e) {
        thrown = true;
    }
    check('polygonArea refuses mismatched lengths', thrown);

    console.log(`test_geodesic: ${total - failed} of ${total} passed`);
    return failed === 0;
}

export const geodesic_ok = test_geodesic();
