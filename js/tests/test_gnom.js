import { proj4 } from '../proj4.js';

// Gnomonic: forward and inverse against reference values from
// PROJ 9.8.1 (through pyproj 3.8.0), sphere and ellipsoid, all aspects.
function test_gnom() {
    const cases = [
        ['+proj=longlat +R=6370997', '+proj=gnom +R=6370997 +lat_0=-42 +lon_0=147', [
            [[150, -30], [295478.296841106, 1350242.102235714]],
            [[100, -60], [-2797139.686113018, -3177553.796137307]],
            [[-170, -30], [4672935.288165815, 413276.766432783]],
        ]],
        ['+proj=longlat +R=6370997', '+proj=gnom +R=6370997 +lat_0=0 +lon_0=0', [
            [[10, 20], [1123378.665112688, 2354625.320024798]],
        ]],
        ['+proj=longlat +R=6370997', '+proj=gnom +R=6370997 +lat_0=90', [
            [[10, 60], [638729.541960903, -3622415.238975614]],
        ]],
        ['+proj=longlat +R=6370997', '+proj=gnom +R=6370997 +lat_0=-90 +lon_0=30', [
            [[10, -60], [-1258051.610002067, 3456468.390989212]],
        ]],
        ['+proj=longlat +ellps=WGS84', '+proj=gnom +ellps=WGS84 +lat_0=-42 +lon_0=147', [
            [[150, -30], [296026.117484902, 1347275.543766303]],
            [[100, -60], [-2807859.531655353, -3182754.956885839]],
            [[-170, -30], [4683812.890027311, 407828.140511455]],
        ]],
        ['+proj=longlat +ellps=WGS84', '+proj=gnom +ellps=WGS84 +lat_0=0 +lon_0=0', [
            [[10, 20], [1124713.736254139, 2342041.514079806]],
        ]],
        ['+proj=longlat +ellps=WGS84', '+proj=gnom +ellps=WGS84 +lat_0=90', [
            [[10, 60], [640943.892175542, -3634973.443090167]],
        ]],
        ['+proj=longlat +ellps=GRS80',
         '+proj=gnom +ellps=GRS80 +lat_0=40 +lon_0=-100 +k_0=0.9996 +x_0=1000 +y_0=-500', [
            [[-90, 45], [795033.992842076, 605635.247739543]],
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

    // A point 90 degrees or more from the centre must be refused: the wasm
    // module returns NaN for points it cannot transform (or throws)
    total += 1;
    let refused = false;
    try {
        for (const [src, dst] of [
            ['+proj=longlat +R=6370997', '+proj=gnom +R=6370997 +lat_0=-42 +lon_0=147'],
            ['+proj=longlat +ellps=WGS84', '+proj=gnom +ellps=WGS84 +lat_0=-42 +lon_0=147'],
        ]) {
            const xy = proj4(src, dst, [-33, 42]);
            refused = !(Number.isFinite(xy[0]) && Number.isFinite(xy[1]));
            if (!refused) {
                break;
            }
        }
    } catch (e) {
        refused = true;
    }
    if (!refused) {
        failed += 1;
    }
    console.log(`${refused ? 'ok' : 'FAIL'} far hemisphere point refused`);

    console.log(`test_gnom: ${total - failed} of ${total} passed`);
    return failed === 0;
}

export const gnom_ok = test_gnom();
