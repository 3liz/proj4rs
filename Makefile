#
# Proj4rs Makefile
#
SHELL = bash
.ONESHELL:

main:
	@echo "Actions: "
	echo "  npm-package: Create wasm npm bundle"


npm-package:
	@rm -rf js/pkg-bundle
	cargo make --cwd proj4rs  wasm_bundle
	cp README.md js/proj4.js js/pkg-bundle
	python3  << 'EOF'
	import json
	js = json.load(open("js/pkg-bundle/package.json"))
	js["files"].extend(("proj4rs_bg.wasm.d.ts", "proj4.js"))
	json.dump(js, open("js/pkg-bundle/package.json", "w"), indent=4)
	EOF

