# Running the WASM examples

## Running the js tests

There is a `index.html` file for testing the WASM module in a navigator.

For security reasons, you need to run it from a server.
You can start a Python server with the following command:

```bash
python3 -m http.server
```

The server will automatically serve the `index.html` file in the current directory

## Running the OpenLayers demo from Docker container

```bash
.docker/ol-run.sh
```

This will build the Node.js image and run the application. Once the application
is started, navigate to http://localhost:5173/.
