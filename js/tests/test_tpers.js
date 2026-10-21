import { proj4 } from '../proj4.js';

// Tilted and near-sided perspective: forward and inverse against reference
// values from PROJ 9.8.1 (through pyproj 3.8.0).
function test_tpers() {
    const cases = [
        ['+proj=longlat +R=6370997', '+proj=tpers +R=6370997 +lat_0=-42 +lon_0=147 +h=3000000 +tilt=30 +azi=45', [
            [[150, -30], [-575211.027126064, 1036328.413044423]],
            [[140, -50], [344750.161414741, -1371938.406334842]],
        ]],
        ['+proj=longlat +ellps=WGS84', '+proj=tpers +ellps=WGS84 +lat_0=90 +h=1000000 +tilt=20 +azi=-30', [
            [[10, 80], [-527515.170085313, -1542351.175448220]],
        ]],
        ['+proj=longlat +ellps=WGS84', '+proj=tpers +ellps=WGS84 +lat_0=0 +lon_0=0 +h=35786000', [
            [[10, 20], [1027105.010111250, 2152833.713367413]],
        ]],
        ['+proj=longlat +ellps=WGS84', '+proj=nsper +ellps=WGS84 +lat_0=40 +lon_0=-100 +h=3000000', [
            [[-90, 45], [763619.441029603, 584966.429458347]],
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

    // A point beyond the horizon must be refused: the wasm module returns
    // NaN for points it cannot transform (or throws)
    total += 1;
    let refused = false;
    try {
        const xy = proj4('+proj=longlat +R=6370997', '+proj=tpers +R=6370997 +lat_0=-42 +lon_0=147 +h=3000000 +tilt=30 +azi=45', [-33, 42]);
        refused = !(Number.isFinite(xy[0]) && Number.isFinite(xy[1]));
    } catch (e) {
        refused = true;
    }
    if (!refused) {
        failed += 1;
    }
    console.log(`${refused ? 'ok' : 'FAIL'} point beyond the horizon refused`);

    console.log(`test_tpers: ${total - failed} of ${total} passed`);
    return failed === 0;
}

export const tpers_ok = test_tpers();
