# OpenLayers + Proj4rs with Vite

This example demonstrates how `OpenLayers` can be used with `Proj4rs`. It is based on `OpenLayers` + [Vite](https://vitejs.dev/).

To get started, run the following (requires Node 14+):

    cd ol-proj4rs-demo-app
    npm update
    npm start

Then go to http://localhost:5173 with your browser, you have started start a development server.

To generate a build ready for production:

    npm run build

Then deploy the contents of the `dist` directory to your server.  You can also run `npm run serve` to serve the results of the `dist` directory for preview.


## Running the demo from Docker container

From the root of this repository, run the cli command 

```
.docker/ol-run.sh
```

This will build the nodejs image and run the application. Once the application
is started, navigate no  http://localhost:5173/.



