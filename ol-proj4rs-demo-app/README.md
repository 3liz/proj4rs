# OpenLayers + Proj4rs with Vite

This example demonstrates how `OpenLayers` can be used with `Proj4rs`. It is based on `OpenLayers` + [Vite](https://vitejs.dev/).

To get started, run the following (requires Node 14+):

    cd ol-proj4rs-demo-app
    npm update
    cp ../js/proj4.js assets/js/
    cp ../pkg/proj4rs_bg.wasm assets/pkg
    cp ../pkg/proj4rs_bg.wasm.d.ts assets/pkg
    cp ../pkg/proj4rs.js assets/pkg
    cp ../pkg/proj4rs.d.ts assets/pkg
    npm start

Then go to http://localhost:5173 with your browser, you have started start a development server.

To generate a build ready for production:

    npm run build

Then deploy the contents of the `dist` directory to your server.  You can also run `npm run serve` to serve the results of the `dist` directory for preview.
