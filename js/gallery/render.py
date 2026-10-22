#!/usr/bin/env python3
"""Render the proj4rs projection gallery headlessly and report the results.

Serves the js/ directory, opens gallery/index.html in Chromium through
Playwright, waits for every panel, then writes:

  out/gallery.png        the whole page
  out/<id>.png           one image per projection
  out/results.json       per-projection comparison with PROJ

and prints one line per projection. Exits 1 if any projection differs
from PROJ beyond its tolerance, unless the difference
is recorded as known in projections.json.

Usage: python3 render.py [--only id1,id2] [--out DIR]
Needs: the wasm module built into js/pkg (cargo make wasm in proj4rs/),
reference.json (make_reference.py) and the playwright python package.
"""
import argparse
import functools
import http.server
import json
import os
import re
import sys
import threading

from playwright.sync_api import sync_playwright

HERE = os.path.dirname(os.path.abspath(__file__))
JS_DIR = os.path.dirname(HERE)


def serve():
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=JS_DIR)
    handler.log_message = lambda *a: None
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    return httpd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="comma separated projection ids")
    ap.add_argument("--out", default=os.path.join(HERE, "out"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    httpd = serve()
    url = f"http://127.0.0.1:{httpd.server_port}/gallery/index.html"
    if args.only:
        url += "?only=" + args.only

    with sync_playwright() as pw:
        exe = os.environ.get("CHROMIUM")
        browser = pw.chromium.launch(executable_path=exe) if exe else pw.chromium.launch()
        page = browser.new_page(viewport={"width": 1360, "height": 900}, device_scale_factor=2)
        page.on("pageerror", lambda e: print("PAGEERROR", e, file=sys.stderr))
        page.goto(url)
        page.wait_for_function("window.galleryResults !== undefined", timeout=120000)
        results = page.evaluate("window.galleryResults")
        if isinstance(results, dict) and "error" in results:
            sys.exit(f"gallery failed: {results['error']}")
        page.screenshot(path=os.path.join(args.out, "gallery.png"), full_page=True)
        for i, el in enumerate(page.query_selector_all(".panel")):
            name = re.sub(r"[^A-Za-z0-9]+", "-", results[i]["id"]).strip("-")
            el.screenshot(path=os.path.join(args.out, f"{name}.png"))
        ua = page.evaluate("navigator.userAgent")
        browser.close()
    httpd.shutdown()

    with open(os.path.join(args.out, "results.json"), "w") as f:
        json.dump(results, f, indent=1)
    bad = 0
    for r in results:
        bad += r["verdict"] in ("differs", "unexpected match")
        print(f"{r['verdict']:9} {r['id']:16} {r['compared']:4}/{r['points']:<4} "
              f"fwd {r['dxy']:.1e} at {r['worstFwd']}  inv {r['dll']:.1e} deg at {r['worstInv']}  "
              f"visibility {r['visibility']}  inverse failures {r['invFailed']}")
    known = sum(r["verdict"] == "known" for r in results)
    print(f"{len(results) - bad - known} of {len(results)} match PROJ, {known} known differences, {bad} unexpected ({ua})")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
