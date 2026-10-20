import { proj4 } from '../proj4.js';

// Cassini, equidistant conic, Robinson and Equal Earth: forward and inverse
// against reference values from PROJ 9.8.1 (through pyproj 3.8.0).
function test_four() {
    const cases = [
        ['+proj=longlat +R=6370997', '+proj=cass +R=6370997 +lat_0=-42 +lon_0=147', [
            [[150, -30], [288859.731012250, 1330553.868361174]],
        ], 1e-6],
        ['+proj=longlat +ellps=WGS84', '+proj=cass +ellps=WGS84 +lat_0=-42 +lon_0=147 +x_0=500000 +y_0=10000000', [
            [[150, -30], [789425.749226268, 11327731.424806831]],
        ], 1e-4],
        ['+proj=longlat +R=6370997', '+proj=eqdc +R=6370997 +lat_1=-20 +lat_2=-40 +lon_0=135', [
            [[150, -30], [1418501.215579857, -3428347.067879791]],
        ], 1e-6],
        ['+proj=longlat +ellps=GRS80', '+proj=eqdc +ellps=GRS80 +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60', [
            [[-120, 30], [-2193479.005798156, -818253.437336989]],
        ], 1e-4],
        ['+proj=longlat +R=6370997', '+proj=robin +R=6370997', [
            [[-100, 60], [-7536475.401958603, 6328945.808890665]],
        ], 1e-6],
        ['+proj=longlat +ellps=WGS84', '+proj=eqearth +ellps=WGS84', [
            [[-100, 60], [-7225128.360577337, 6924035.137552897]],
        ], 1e-6],
    ];

    let failed = 0;
    let total = 0;
    for (const [src, dst, points, tol] of cases) {
        const t = proj4(src, dst);
        for (const [lonlat, expected] of points) {
            total += 1;
            const xy = t.forward(lonlat.slice());
            const back = t.inverse(xy.slice());
            const dxy = Math.hypot(xy[0] - expected[0], xy[1] - expected[1]);
            const dll = Math.hypot(back[0] - lonlat[0], back[1] - lonlat[1]);
            const ok = dxy < tol && dll < 1e-7;
            if (!ok) {
                failed += 1;
            }
            console.log(`${ok ? 'ok' : 'FAIL'} ${dst} ${lonlat} => ${xy} (diff ${dxy.toExponential(2)} m, round trip ${dll.toExponential(2)} deg)`);
        }
    }

    console.log(`test_four: ${total - failed} of ${total} passed`);
    return failed === 0;
}

export const four_ok = test_four();
