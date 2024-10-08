import Map from 'ol/Map.js';
import View from 'ol/View.js';
import {proj4} from './assets/js/proj4.js';
import {getCenter} from 'ol/extent.js';
import {get as getProjection} from 'ol/proj.js';
import {register} from 'ol/proj/proj4.js';
import GeoJSON from 'ol/format/GeoJSON.js';
import Graticule from 'ol/layer/Graticule.js';
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
proj5479.setWorldExtent([0, 0, 0, 0]);

const proj3413 = getProjection('EPSG:3413');
proj3413.setExtent([-4194304, -4194304, 4194304, 4194304]);
proj3413.setWorldExtent([-179.99, -40, 179.99, 84]);

const proj2163 = getProjection('EPSG:2163');
proj2163.setExtent([-8040784.5135, -2577524.921, 3668901.4484, 4785105.1096]);
proj2163.setWorldExtent([-180, 10, -10, 84]);

const proj54009 = getProjection('ESRI:54009');
proj54009.setExtent([-18019909.21177587, -9009954.605703328, 18019909.21177587, 9009954.605703328]);
proj54009.setWorldExtent([-179, -89.99, 179, 89.99]);

const viewProjSelect = document.getElementById('view-projection');

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
    projection: getProjection(viewProjSelect.value),
    center: getCenter(getProjection(viewProjSelect.value).getExtent() || [0, 0, 0, 0]),
    zoom: 0,
    extent: getProjection(viewProjSelect.value).getExtent() || undefined,
  }),
});

function updateViewProjection() {
  const newProj = getProjection(viewProjSelect.value);
  const newProjExtent = newProj.getExtent();
  const newView = new View({
    projection: newProj,
    center: getCenter(newProjExtent || [0, 0, 0, 0]),
    zoom: 0,
    extent: newProjExtent || undefined,
  });
  updateMapVar(newView)
}

function updateMapVar(view) {
  map.setView(view);

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

  map.setLayers([vectorMap, new Graticule()]);
}

/**
 * Handle change event.
 */
viewProjSelect.onchange = function () {
  updateViewProjection();
};
