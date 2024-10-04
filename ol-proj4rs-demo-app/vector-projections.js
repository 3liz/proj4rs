import Map from 'ol/Map.js';
import TileGrid from 'ol/tilegrid/TileGrid.js';
import TileLayer from 'ol/layer/Tile.js';
import View from 'ol/View.js';
import WMTS, {optionsFromCapabilities} from 'ol/source/WMTS.js';
import WMTSCapabilities from 'ol/format/WMTSCapabilities.js';
import {proj4} from 'proj4rs/proj4.js';
import {OSM, TileImage, TileWMS} from 'ol/source.js';
import {createXYZ} from 'ol/tilegrid.js';
import {getCenter, getWidth} from 'ol/extent.js';
import {get as getProjection, transformExtent} from 'ol/proj.js';
import {register} from 'ol/proj/proj4.js';

/*
const layers = {};

layers['osm'] = new TileLayer({
  source: new OSM(),
});

layers['wms4326'] = new TileLayer({
  source: new TileWMS({
    url: 'https://ahocevar.com/geoserver/gwc/service/wms',
    crossOrigin: '',
    params: {
      'LAYERS': 'ne:NE1_HR_LC_SR_W_DR',
      'TILED': true,
      'VERSION': '1.1.1',
    },
    projection: 'EPSG:4326',
    // Source tile grid (before reprojection)
    tileGrid: createXYZ({
      extent: [-180, -90, 180, 90],
      maxResolution: 360 / 512,
      maxZoom: 10,
    }),
    // Accept a reprojection error of 2 pixels
    reprojectionErrorThreshold: 2,
  }),
});

layers['wms21781'] = new TileLayer({
  source: new TileWMS({
    attributions:
      '© <a href="https://shop.swisstopo.admin.ch/en/products/maps/national/lk1000"' +
      'target="_blank">Pixelmap 1:1000000 / geo.admin.ch</a>',
    crossOrigin: 'anonymous',
    params: {
      'LAYERS': 'ch.swisstopo.pixelkarte-farbe-pk1000.noscale',
      'FORMAT': 'image/jpeg',
    },
    url: 'https://wms.geo.admin.ch/',
    projection: 'EPSG:21781',
  }),
});

const parser = new WMTSCapabilities();

layers['wmts3413'] = new TileLayer();
const urlA =
  'https://map1.vis.earthdata.nasa.gov/wmts-arctic/' +
  'wmts.cgi?SERVICE=WMTS&request=GetCapabilities';
fetch(urlA)
  .then(function (response) {
    return response.text();
  })
  .then(function (text) {
    const result = parser.read(text);
    const options = optionsFromCapabilities(result, {
      layer: 'OSM_Land_Mask',
      matrixSet: 'EPSG3413_250m',
    });
    options.crossOrigin = '';
    options.projection = 'EPSG:3413';
    options.wrapX = false;
    layers['wmts3413'].setSource(new WMTS(options));
  });

layers['bng'] = new TileLayer();
const urlB =
  'https://tiles.arcgis.com/tiles/qHLhLQrcvEnxjtPr/arcgis/rest/services/OS_Open_Raster/MapServer/WMTS';
fetch(urlB)
  .then(function (response) {
    return response.text();
  })
  .then(function (text) {
    const result = parser.read(text);
    const options = optionsFromCapabilities(result, {
      layer: 'OS_Open_Raster',
    });
    options.attributions =
      'Contains OS data © Crown Copyright and database right ' +
      new Date().getFullYear();
    options.crossOrigin = '';
    options.projection = 'EPSG:27700';
    options.wrapX = false;
    layers['bng'].setSource(new WMTS(options));
  });

const startResolution = getWidth(getProjection('EPSG:3857').getExtent()) / 256;
const resolutions = new Array(22);
for (let i = 0, ii = resolutions.length; i < ii; ++i) {
  resolutions[i] = startResolution / Math.pow(2, i);
}

layers['states'] = new TileLayer({
  source: new TileWMS({
    url: 'https://ahocevar.com/geoserver/wms',
    crossOrigin: '',
    params: {'LAYERS': 'topp:states'},
    serverType: 'geoserver',
    tileGrid: new TileGrid({
      extent: [-13884991, 2870341, -7455066, 6338219],
      resolutions: resolutions,
      tileSize: [512, 256],
    }),
    projection: 'EPSG:3857',
  }),
});

const map = new Map({
  layers: [layers['osm'], layers['bng']],
  target: 'map',
  view: new View({
    projection: 'EPSG:3857',
    center: [0, 0],
    zoom: 2,
  }),
});
 */

import GeoJSON from 'ol/format/GeoJSON.js';
import Graticule from 'ol/layer/Graticule.js';
import Projection from 'ol/proj/Projection.js';
import VectorLayer from 'ol/layer/Vector.js';
import VectorSource from 'ol/source/Vector.js';
import {Fill, Style} from 'ol/style.js';


proj4.defs(
  'EPSG:27700',
  '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 ' +
  '+x_0=400000 +y_0=-100000 +ellps=airy ' +
  '+towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489 ' +
  '+units=m +no_defs'
);
proj4.defs(
  'EPSG:23032',
  '+proj=utm +zone=32 +ellps=intl ' +
  '+towgs84=-87,-98,-121,0,0,0,0 +units=m +no_defs'
);
proj4.defs(
  'EPSG:5479',
  '+proj=lcc +lat_1=-76.66666666666667 +lat_2=' +
  '-79.33333333333333 +lat_0=-78 +lon_0=163 +x_0=7000000 +y_0=5000000 ' +
  '+ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
);
proj4.defs(
  'EPSG:3413',
  '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 ' +
  '+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
);
proj4.defs(
  'EPSG:2163',
  '+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 ' +
  '+a=6370997 +b=6370997 +units=m +no_defs'
);
proj4.defs(
  'ESRI:54009',
  '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 ' + '+units=m +no_defs'
);
register(proj4);

const proj27700 = getProjection('EPSG:27700');
proj27700.setExtent([-650000, -150000, 1350000, 1450000]);
proj27700.setWorldExtent([-65, -15, 135, 145]);

const proj23032 = getProjection('EPSG:23032');
proj23032.setExtent([-1206118.71, 4021309.92, 1295389.0, 8051813.28]);
proj23032.setWorldExtent([-121, 20, 130, 75]);

const proj5479 = getProjection('EPSG:5479');
proj5479.setExtent([6825737.53, 4189159.8, 9633741.96, 5782472.71]);
proj5479.setWorldExtent([0,0,0,0]);

const proj3413 = getProjection('EPSG:3413');
proj3413.setExtent([-4194304, -4194304, 4194304, 4194304]);
proj3413.setWorldExtent([-179.99, -40, 179.99, 84]);

const proj2163 = getProjection('EPSG:2163');
proj2163.setExtent([-8040784.5135, -2577524.921, 3668901.4484, 4785105.1096]);
proj2163.setWorldExtent([-180, 10, -10, 84]);

const proj54009 = getProjection('ESRI:54009');
proj54009.setExtent([-18019909.21177587, -9009954.605703328, 18019909.21177587, 9009954.605703328]);
proj54009.setWorldExtent([-179, -89.99, 179, 89.99]);

// Configure the Sphere Mollweide projection object with an extent,
// and a world extent. These are required for the Graticule.
/*const sphereMollweideProjection = new Projection({
  code: 'ESRI:54009',
  extent: [
    -18019909.21177587, -9009954.605703328, 18019909.21177587,
    9009954.605703328,
  ],
  worldExtent: [-179, -89.99, 179, 89.99],
});

 */

const style = new Style({
  fill: new Fill({
    color: '#eeeeee',
  }),
});

let vectorMap = new VectorLayer({
  source: new VectorSource({
    url: 'https://openlayers.org/data/vector/ecoregions.json',
    format: new GeoJSON(),
  }),
  style: function (feature) {
    const color = feature.get('COLOR_BIO') || '#eeeeee';
    style.getFill().setColor(color);
    return style;
  },
})

const map = new Map({
  keyboardEventTarget: document,
  layers: [
    vectorMap,
    new Graticule(),
  ],
  target: 'map',
  view: new View({
    center: [0, 0],
    projection: 'EPSG:3857',
    zoom: 2,
  }),
});

const viewProjSelect = document.getElementById('view-projection');

function updateViewProjection() {
  const newProj = getProjection(viewProjSelect.value);
  console.log(newProj)
  const newProjExtent = newProj.getExtent();
  console.log(newProjExtent)
  const newView = new View({
    projection: newProj,
    center: getCenter(newProjExtent || [0, 0, 0, 0]),
    zoom: 0,
    extent: newProjExtent || undefined,
  });
  map.setView(newView);

  console.log(map.getLayers().getArray()[0])
}

/**
 * Handle change event.
 */
viewProjSelect.onchange = function () {
  updateViewProjection();
};

updateViewProjection();
