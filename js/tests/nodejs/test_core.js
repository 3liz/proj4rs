let Proj = await import("../../pkg-node/proj4rs.js");

function test_core() {
    //EPSG:3006
    var sweref99tm = '+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs';
    // EPSG:3021
    var rt90 = '+lon_0=15.808277777799999 +lat_0=0.0 +k=1.0 +x_0=1500000.0 +y_0=0.0 +proj=tmerc +ellps=bessel +units=m +towgs84=414.1,41.3,603.1,-0.855,2.141,-7.023,0 +no_defs';
    let from = new Proj.Projection(sweref99tm);
    let to = new Proj.Projection(rt90);
    let point = new Proj.Point(319180, 6399862, 0.0);
    Proj.transform(from, to, point);
    console.log(`=> ${point.x} ${point.y}`);

    var epsg2154 = '+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs';
    var epsg3857 = '+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs +type=crs';
    from = new Proj.Projection(epsg2154);
    to = new Proj.Projection(epsg3857);
    point = new Proj.Point(489353.59, 6587552.2, 0.0);
    Proj.transform(from, to, point);
    console.log(`=> ${point.x} ${point.y}`);
}

test_core();
