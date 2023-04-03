import { proj4 } from '../proj4.js';

var defs = proj4.defs;

function test_defs() {
   defs(
        'EPSG:27700',
        '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 ' +
            '+x_0=400000 +y_0=-100000 +ellps=airy ' +
            '+towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489 ' +
            '+units=m +no_defs'
  );
  defs(
    'EPSG:23032',
    '+proj=utm +zone=32 +ellps=intl ' +
        '+towgs84=-87,-98,-121,0,0,0,0 +units=m +no_defs'
  );
  defs(
    'EPSG:5479',
    '+proj=lcc +lat_1=-76.66666666666667 +lat_2=' +
        '-79.33333333333333 +lat_0=-78 +lon_0=163 +x_0=7000000 +y_0=5000000 ' +
        '+ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
  );
  /* somerc projection not supported yet
  defs(
    'EPSG:21781',
    '+proj=somerc +lat_0=46.95240555555556 ' +
        '+lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel ' +
        '+towgs84=674.4,15.1,405.3,0,0,0,0 +units=m +no_defs'
  );*/
  defs(
    'EPSG:3413',
    '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 ' +
        '+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
  );
  /* laea projection not supported yet
  defs(
    'EPSG:2163',
    '+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 ' +
        '+a=6370997 +b=6370997 +units=m +no_defs'
  );*/
  /* moll projection not supported yet
  defs(
    'ESRI:54009',
    '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 ' + '+units=m +no_defs'
  );*/
  console.log(defs('EPSG:27700'))
  console.log(defs['EPSG:27700'])
  console.log(defs['EPSG:4326']);
  console.log(defs['EPSG:4269']);
  console.log(defs['EPSG:3857']);
  console.log(defs.WGS84 == defs('EPSG:4326'));
  console.log(defs['EPSG:3785'] == defs('EPSG:3785'));
  console.log(defs.GOOGLE == defs('EPSG:3785'));
  console.log(defs['EPSG:900913'] == defs('EPSG:3785'));
  console.log(defs['EPSG:102113'] == defs('EPSG:3785'));
}

test_defs();

