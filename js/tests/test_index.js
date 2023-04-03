import { proj4 } from '../proj4.js';

function test_index() {
        proj4.defs('EPSG:3006', '+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs');
        proj4.defs('EPSG:3021', '+lon_0=15.808277777799999 +lat_0=0.0 +k=1.0 +x_0=1500000.0 +y_0=0.0 +proj=tmerc +ellps=bessel +units=m +towgs84=414.1,41.3,603.1,-0.855,2.141,-7.023,0 +no_defs');
        console.log(proj4(proj4.defs('EPSG:3006'), proj4.defs('EPSG:3021'), [319180, 6399862]));
        console.log(proj4(proj4.defs('EPSG:3006'), proj4.defs('EPSG:3021')).forward([319180, 6399862]));

        proj4.defs('EPSG:2154', '+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs');
        proj4.defs('EPSG:3857', '+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs +type=crs');
        console.log(proj4(proj4.defs('EPSG:2154'), proj4.defs('EPSG:3857')).forward([489353.59, 6587552.2]));
}

test_index();
