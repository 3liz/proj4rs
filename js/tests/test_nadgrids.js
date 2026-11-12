import { proj4, Proj } from '../proj4.js';

// Nadgrids loaded from JS with `add_nadgrid`.
//
// Grid files are fetched from `./fixtures` (a link to `proj4rs/fixtures`).
//
// Reference values were computed with the native (non-wasm) proj4rs build,
// reading the same grid files from disk: ED50 (Catalonia) to ETRS89.
const DST = '+proj=longlat +ellps=GRS80 +towgs84=0,0,0';

const POINTS = [
    [[2.0, 41.5], [1.998844997221, 41.498877805604]],
    [[1.5, 42.0], [1.498826419380, 41.998891316652]],
    [[0.8, 41.0], [0.798813254939, 40.998857449426]],
    [[2.8, 42.3], [2.798862839990, 42.298906350562]],
];

async function fetchGrid(name) {
    const response = await fetch(new URL(`./fixtures/${name}`, import.meta.url));
    if (!response.ok) {
        throw new Error(`Failed to fetch ${name}: ${response.status}`);
    }
    return new Uint8Array(await response.arrayBuffer());
}

function isUnknownFormat(e) {
    return String(e?.message ?? e).includes('Unknown grid format');
}

async function test_nadgrids() {
    let failed = 0;
    let total = 0;

    const check = (ok, msg) => {
        total += 1;
        if (!ok) {
            failed += 1;
        }
        console.log(`${ok ? 'ok' : 'FAIL'} ${msg}`);
    };

    const grids = [
        // [file, requires the 'tiff' feature]
        ['100800401.gsb', false],
        ['es_cat_icgc_100800401.tif', true],
        ['es_cat_icgc_100800401_pixel.tif', true],
    ];

    for (const [name, isTiff] of grids) {
        const data = await fetchGrid(name);
        try {
            Proj.add_nadgrid(name, data);
        } catch (e) {
            if (isTiff && isUnknownFormat(e)) {
                // The wasm module is not built with the 'tiff' feature
                console.log(`skip ${name}: tiff support not enabled`);
                continue;
            }
            check(false, `add_nadgrid ${name}: ${e}`);
            continue;
        }
        check(true, `add_nadgrid ${name} (${data.length} bytes)`);

        const src = `+proj=longlat +ellps=intl +nadgrids=${name}`;
        const t = proj4(src, DST);
        for (const [lonlat, expected] of POINTS) {
            const xy = t.forward(lonlat.slice());
            const back = t.inverse(xy.slice());
            const dxy = Math.hypot(xy[0] - expected[0], xy[1] - expected[1]);
            const dll = Math.hypot(back[0] - lonlat[0], back[1] - lonlat[1]);
            check(
                dxy < 1e-9 && dll < 1e-9,
                `${name} ${lonlat} => ${xy} (diff ${dxy.toExponential(2)} deg, round trip ${dll.toExponential(2)} deg)`,
            );
        }
    }

    // Invalid data must be rejected
    try {
        Proj.add_nadgrid('invalid.gsb', new Uint8Array([1, 2, 3, 4]));
        check(false, 'invalid grid data rejected');
    } catch (e) {
        check(true, `invalid grid data rejected (${e.message ?? e})`);
    }

    // Truncated data must be rejected (reading past the end of the array)
    try {
        const data = await fetchGrid('100800401.gsb');
        Proj.add_nadgrid('truncated.gsb', data.subarray(0, 1024));
        check(false, 'truncated grid data rejected');
    } catch (e) {
        check(true, `truncated grid data rejected (${e.message ?? e})`);
    }

    console.log(`test_nadgrids: ${total - failed} of ${total} passed`);
    return failed === 0;
}

export const nadgrids_ok = await test_nadgrids();
