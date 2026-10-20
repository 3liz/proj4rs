import { proj4 } from '../proj4.js';

// Oblique Mercator: forward and inverse against reference values from
// PROJ 9.8.1 (through pyproj 3.8.0), alpha/gamma, gamma only, two points.
function test_omerc() {
    const cases = [
        ['+proj=longlat +ellps=GRS80',
         '+proj=omerc +lat_0=4 +lonc=102.25 +alpha=323.0257964666666 +k=0.99984 +x_0=804671 +y_0=0 +no_uoff +gamma=323.1301023611111 +ellps=GRS80', [
            [[101.5, 3.1], [389287.174191490, 343128.619662765]],
            [[104, 1.5], [667104.332086940, 165868.754500688]],
        ]],
        ['+proj=longlat +ellps=evrstSS',
         '+proj=omerc +lat_0=4 +lonc=115 +alpha=53.31582047222222 +k=0.99984 +x_0=590476.87 +y_0=442857.65 +gamma=53.13010236111111 +ellps=evrstSS', [
            [[117, 6], [811179.866049288, 664991.329329078]],
        ]],
        ['+proj=longlat +ellps=GRS80',
         '+proj=omerc +lat_0=57 +lonc=-133.6666666666667 +alpha=323.1301023611111 +k=0.9999 +x_0=5000000 +y_0=-5000000 +no_uoff +gamma=323.1301023611111 +ellps=GRS80', [
            [[-135, 58], [739841.600324485, 687227.401595449]],
        ]],
        ['+proj=longlat +ellps=WGS84', '+proj=omerc +lat_0=40 +lonc=10 +gamma=20 +ellps=WGS84', [
            [[12, 42], [139309.854872632, 241309.802316528]],
        ]],
        ['+proj=longlat +ellps=WGS84',
         '+proj=omerc +lat_1=45 +lon_1=-100 +lat_2=50 +lon_2=-80 +lat_0=47 +ellps=WGS84', [
            [[-80, 50], [993675.885257817, 423831.381319845]],
        ]],
        ['+proj=longlat +R=6370997', '+proj=omerc +R=6370997 +lat_0=-42 +lonc=147 +alpha=30', [
            [[150, -30], [294093.952081180, 1329914.798619571]],
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

    console.log(`test_omerc: ${total - failed} of ${total} passed`);
    return failed === 0;
}

export const omerc_ok = test_omerc();
