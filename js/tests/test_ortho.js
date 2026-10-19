import { proj4 } from '../proj4.js';

// Orthographic: forward and inverse against reference values from
// PROJ 9.8.1 (through pyproj 3.8.0), sphere and ellipsoid, all aspects.
function test_ortho() {
    const cases = [
        ['+proj=longlat +R=6370997', '+proj=ortho +R=6370997 +lat_0=-42 +lon_0=147', [
            [[150, -30], [288760.773147241, 1319545.149494376]],
            [[100, -60], [-2329726.119477315, -2646571.464291994]],
        ]],
        ['+proj=longlat +R=6370997', '+proj=ortho +R=6370997 +lat_0=90', [
            [[10, 60], [553156.009485740, -3137103.620008760]],
        ]],
        ['+proj=longlat +ellps=WGS84', '+proj=ortho +ellps=WGS84 +lat_0=-42 +lon_0=147', [
            [[150, -30], [289326.597573795, 1316745.581434381]],
            [[100, -60], [-2338214.274110457, -2649929.082205364]],
            [[-133, -37], [5022512.762328568, -2265554.666728071]],
        ]],
        ['+proj=longlat +ellps=WGS84', '+proj=ortho +ellps=WGS84 +lat_0=0 +lon_0=0', [
            [[10, 20], [1041166.061652678, 2167696.787828757]],
        ]],
        ['+proj=longlat +ellps=WGS84', '+proj=ortho +ellps=WGS84 +lat_0=-90 +lon_0=30', [
            [[10, -60], [-1093474.169046881, 3004295.588213213]],
        ]],
        ['+proj=longlat +ellps=WGS84',
         '+proj=ortho +ellps=WGS84 +lat_0=40 +lon_0=-100 +alpha=30 +k_0=0.9996 +x_0=1000 +y_0=-500', [
            [[-90, 45], [380806.930194124, 909970.369811893]],
        ]],
    ];

    let failed = 0;
    let total = 0;
    for (const [src, dst, points] of cases) {
        const t = proj4(src, dst);
        for (const [lonlat, expected] of points) {
            total += 1;
            const xy = t.forward(lonlat.slice());
            const back = t.inverse(xy.slice());
            const dxy = Math.hypot(xy[0] - expected[0], xy[1] - expected[1]);
            const dll = Math.hypot(back[0] - lonlat[0], back[1] - lonlat[1]);
            const ok = dxy < 1e-6 && dll < 1e-9;
            if (!ok) {
                failed += 1;
            }
            console.log(`${ok ? 'ok' : 'FAIL'} ${dst} ${lonlat} => ${xy} (diff ${dxy.toExponential(2)} m, round trip ${dll.toExponential(2)} deg)`);
        }
    }

    // A point on the far hemisphere must be refused: the wasm module
    // returns NaN for points it cannot transform (or throws)
    total += 1;
    let refused = false;
    try {
        const xy = proj4('+proj=longlat +ellps=WGS84', '+proj=ortho +ellps=WGS84 +lat_0=-42 +lon_0=147', [-33, 42]);
        refused = !(Number.isFinite(xy[0]) && Number.isFinite(xy[1]));
    } catch (e) {
        refused = true;
    }
    if (!refused) {
        failed += 1;
    }
    console.log(`${refused ? 'ok' : 'FAIL'} far hemisphere point refused`);

    console.log(`test_ortho: ${total - failed} of ${total} passed`);
    return failed === 0;
}

export const ortho_ok = test_ortho();
