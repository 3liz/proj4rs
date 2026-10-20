#
# Proj4rs Makefile
# Used in gitlab ci 
#
SHELL = bash
.ONESHELL:

clean:
	rm -rf pkg *.tgz

# Build for web
build:
	{\
	  set -e; \
	  wasm-pack build \
	  --target=web \
	  --features=logging \
	  --features=proj4js-compat \
	  --features=with-wasm-entrypoint; \
	  cd pkg && ln -sf ../../js/proj4.js; \
	}

pkg-%:
	@cat pkg/package.json | jq -r .$*

CI_COMMIT_SHORT_SHA ?= $(shell git rev-parse --short HEAD)

manifest:
	@{ \
	  buildid=$$(date +"%Y%m%d%H%M"); \
	  echo "commit=$(CI_COMMIT_SHORT_SHA)" > pkg/manifest; \
	  echo echo "buildid=$$buildid" >> pkg/manifest; \
	  cat pkg/manifest; \
	}

package: manifest
	{ \
	  ls -al;
	  name=$$($(MAKE) -s pkg-name); \
	  version=$$($(MAKE) -s pkg-version); \
	  package="$$name-$$version.tgz"; \
	  echo "Building package $$package"; \
	  tar -czvf $$package \
	  --dereference \
	  --transform=s/pkg/$$name/ pkg; \
	}

