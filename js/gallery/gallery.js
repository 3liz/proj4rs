// proj4rs projection gallery
//
// Every map is drawn by proj4rs (wasm, through the proj4js-compatible
// wrapper): Natural Earth coastline in grey, graticule in light blue.
// The red rings are the graticule intersections as computed by PROJ
// (reference.json, from make_reference.py). Where proj4rs agrees with PROJ
// the rings sit exactly on the line crossings. Each panel also reports the
// largest forward difference (projected units) and the largest inverse
// difference (degrees, proj4rs inverse of the PROJ coordinates).

import { proj4 } from '../proj4.js';

const SIZE = 320;
const PAD = 12;

async function load(url) {
    const r = await fetch(url);
    if (!r.ok) {
        throw new Error(`${url}: ${r.status}`);
    }
    return r.json();
}

// The source lon/lat system uses the same earth model as the projection so
// that no datum shift is involved.
function longlatFor(proj) {
    const keep = proj.split(/\s+/).filter(
        (t) => /^\+(ellps|R|a|b|rf|f|es|e)=/.test(t));
    return ['+proj=longlat', ...keep].join(' ');
}

function makeForward(t) {
    return (lon, lat) => {
        try {
            const [x, y] = t.forward([lon, lat]);
            return Number.isFinite(x) && Number.isFinite(y) ? [x, y] : null;
        } catch (e) {
            return null;
        }
    };
}

function makeInverse(t) {
    return (x, y) => {
        try {
            const [lon, lat] = t.inverse([x, y]);
            return Number.isFinite(lon) && Number.isFinite(lat) ? [lon, lat] : null;
        } catch (e) {
            return null;
        }
    };
}

function inWindow(lon, lat, e) {
    return lon >= e.lon[0] && lon <= e.lon[1] && lat >= e.lat[0] && lat <= e.lat[1];
}

// Shift a longitude into the window when the window extends past 180.
function wrapLon(lon, e) {
    if (lon < e.lon[0]) {
        lon += 360;
    }
    if (lon > e.lon[1]) {
        lon -= 360;
    }
    return lon;
}

function graticule(e) {
    const lines = [];
    const d = 0.5;
    for (let lon = e.lon[0]; lon <= e.lon[1] + 1e-9; lon += e.step) {
        const line = [];
        for (let lat = e.lat[0]; lat <= e.lat[1] + 1e-9; lat += d) {
            line.push([lon, lat]);
        }
        lines.push(line);
    }
    for (let lat = e.lat[0]; lat <= e.lat[1] + 1e-9; lat += e.step) {
        const line = [];
        for (let lon = e.lon[0]; lon <= e.lon[1] + 1e-9; lon += d) {
            line.push([lon, lat]);
        }
        lines.push(line);
    }
    return lines;
}

function coastline(geojson, e) {
    const lines = [];
    for (const f of geojson.features) {
        const g = f.geometry;
        const parts = g.type === 'LineString' ? [g.coordinates] : g.coordinates;
        for (const part of parts) {
            let line = [];
            for (const [lon0, lat] of part) {
                const lon = wrapLon(lon0, e);
                if (inWindow(lon, lat, e)) {
                    line.push([lon, lat]);
                } else if (line.length) {
                    lines.push(line);
                    line = [];
                }
            }
            if (line.length) {
                lines.push(line);
            }
        }
    }
    return lines;
}

// Project lon/lat polylines, splitting where a point does not project.
function projectLines(lines, fwd) {
    const out = [];
    for (const line of lines) {
        let cur = [];
        for (const [lon, lat] of line) {
            const p = fwd(lon, lat);
            if (p) {
                cur.push(p);
            } else if (cur.length) {
                out.push(cur);
                cur = [];
            }
        }
        if (cur.length > 1) {
            out.push(cur);
        }
    }
    return out;
}

function bounds(lines, points) {
    let [x0, y0, x1, y1] = [Infinity, Infinity, -Infinity, -Infinity];
    const add = ([x, y]) => {
        x0 = Math.min(x0, x); y0 = Math.min(y0, y);
        x1 = Math.max(x1, x); y1 = Math.max(y1, y);
    };
    lines.forEach((l) => l.forEach(add));
    points.forEach(add);
    return [x0, y0, x1, y1];
}

function drawLines(ctx, lines, toPx, jump) {
    ctx.beginPath();
    for (const line of lines) {
        let prev = null;
        for (const p of line) {
            const q = toPx(p);
            // Split at wrap-arounds (a segment jumping across the map).
            if (prev && Math.hypot(q[0] - prev[0], q[1] - prev[1]) < jump) {
                ctx.lineTo(q[0], q[1]);
            } else {
                ctx.moveTo(q[0], q[1]);
            }
            prev = q;
        }
    }
    ctx.stroke();
}

function fmt(v) {
    if (v === 0) {
        return '0';
    }
    return v.toExponential(1);
}

function compare(ref, fwd, inv) {
    const r = {
        points: ref.length, compared: 0, visibility: 0, invFailed: 0,
        dxy: 0, dll: 0, worstFwd: null, worstInv: null,
    };
    for (const [lon, lat, x, y] of ref) {
        const p = fwd(lon, lat);
        if ((x === null) !== (p === null)) {
            r.visibility += 1;
            continue;
        }
        if (p === null) {
            continue;
        }
        r.compared += 1;
        const d = Math.hypot(p[0] - x, p[1] - y);
        if (d > r.dxy) {
            r.dxy = d;
            r.worstFwd = [lon, lat];
        }
        const q = inv(x, y);
        if (q === null) {
            r.invFailed += 1;
            r.worstInv = [lon, lat];
            continue;
        }
        // Longitude is undefined at the poles.
        const dlon = Math.abs(lat) > 89.999999 ? 0
            : (((q[0] - lon) % 360) + 540) % 360 - 180;
        const dl = Math.hypot(dlon * Math.cos(lat * Math.PI / 180), q[1] - lat);
        if (dl > r.dll) {
            r.dll = dl;
            if (!r.invFailed) {
                r.worstInv = [lon, lat];
            }
        }
    }
    return r;
}

function verdict(r, e) {
    const tol = e.tol ?? 1e-3;
    if (r.compared === 0 || r.visibility || r.invFailed || r.dxy > tol || r.dll > 1e-7) {
        return e.expect === 'differs' ? 'known' : 'differs';
    }
    if (e.expect === 'differs') {
        return 'unexpected match';
    }
    return r.dxy === 0 ? 'identical' : 'matches';
}

function panel(e, ref, coast) {
    const src = longlatFor(e.proj);
    const t = proj4(src, e.proj);
    const fwd = makeForward(t);
    const inv = makeInverse(t);

    const div = document.createElement('figure');
    div.className = 'panel';
    const canvas = document.createElement('canvas');
    const dpr = window.devicePixelRatio || 1;
    canvas.width = SIZE * dpr;
    canvas.height = SIZE * dpr;
    canvas.style.width = `${SIZE}px`;
    canvas.style.height = `${SIZE}px`;
    div.appendChild(canvas);

    const result = { id: e.id, proj: e.proj, ...compare(ref, fwd, inv) };
    result.verdict = verdict(result, e);

    const grat = projectLines(graticule(e), fwd);
    const coastXY = projectLines(coast, fwd);
    const refXY = ref.filter((r) => r[2] !== null).map((r) => [r[2], r[3]]);
    const [x0, y0, x1, y1] = bounds(grat, refXY);
    const scale = (SIZE - 2 * PAD) / Math.max(x1 - x0, y1 - y0);
    const ox = (SIZE - (x1 - x0) * scale) / 2;
    const oy = (SIZE - (y1 - y0) * scale) / 2;
    const toPx = ([x, y]) => [ox + (x - x0) * scale, SIZE - oy - (y - y0) * scale];

    const ctx = canvas.getContext('2d');
    ctx.scale(dpr, dpr);
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, SIZE, SIZE);
    ctx.lineWidth = 0.6;
    ctx.strokeStyle = '#9ec5e8';
    drawLines(ctx, grat, toPx, SIZE / 4);
    ctx.lineWidth = 0.8;
    ctx.strokeStyle = '#555555';
    drawLines(ctx, coastXY, toPx, SIZE / 4);
    ctx.lineWidth = 1;
    ctx.strokeStyle = '#d62728';
    for (const p of refXY) {
        const [px, py] = toPx(p);
        ctx.beginPath();
        ctx.arc(px, py, 3, 0, 2 * Math.PI);
        ctx.stroke();
    }

    const cap = document.createElement('figcaption');
    cap.innerHTML = `<b>${e.id}</b> <span class="${result.verdict.replace(' ', '-')}">${result.verdict}</span><br>`
        + `<code>${e.proj}</code><br>`
        + `${result.compared}/${result.points} points, fwd max ${fmt(result.dxy)}, `
        + `inv max ${fmt(result.dll)} deg`
        + (result.visibility ? `, <b>${result.visibility} visibility mismatches</b>` : '')
        + (result.invFailed ? `, <b>${result.invFailed} inverse failures</b>` : '')
        + (e.note ? `<br><i>${e.note}</i>` : '');
    div.appendChild(cap);
    return [div, result];
}

async function main() {
    const [entries, reference, coast] = await Promise.all([
        load('projections.json'),
        load('reference.json'),
        load('data/ne_110m_coastline.geojson'),
    ]);
    const only = new URLSearchParams(location.search).get('only');
    document.getElementById('meta').textContent =
        `Reference: PROJ ${reference.proj_version} (pyproj ${reference.pyproj_version}).`;
    const grid = document.getElementById('grid');
    const results = [];
    for (const raw of entries) {
        if (only && !only.split(',').includes(raw.id)) {
            continue;
        }
        const e = { ...reference.defaults, ...raw };
        const ref = reference.entries[e.id];
        if (!ref) {
            throw new Error(`${e.id}: not in reference.json, rerun make_reference.py`);
        }
        const [div, result] = panel(e, ref, coastline(coast, e));
        grid.appendChild(div);
        results.push(result);
        console.log(`${result.verdict} ${e.id}: ${result.compared}/${result.points} fwd ${fmt(result.dxy)} inv ${fmt(result.dll)} vis ${result.visibility} invfail ${result.invFailed}`);
    }
    window.galleryResults = results;
}

main().catch((err) => {
    console.error(err);
    window.galleryResults = { error: String(err) };
});
