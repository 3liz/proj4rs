//
// Proj4js compatibility wrapper
//
// This module is a wrapper around proj4rs in order to provide the same interface as
// Proj4js
//

// Initialize the proj4rs module
import init, * as Proj from "../pkg/proj4rs.js";
await init();

/// Used as factory for creating projection alias
function defs(name) {
    var that = this;
    if (arguments.length === 2) {
        var def = arguments[1];
        if (typeof def === 'string') {
            if (def.charAt(0) === '+') {
                defs[name] = new Proj.Projection(def);
            }
            else {
                // WKT not supported
                throw new Error('WKT format is not supported!');
            }
        } else if (def instanceof Proj.Projection){
            defs[name] = def;
        } else {
            // Unknown def
            throw new Error('Unknown projection object!');
        }
    }
    else if (arguments.length === 1) {
        if (Array.isArray(name)) {
            return name.map(function(v) {
                if (Array.isArray(v)) {
                    defs.apply(that, v);
                }
                else {
                    defs(v);
                }
            });
        }
        else if (typeof name === 'string') {
            if (name in defs) {
                return defs[name];
            }
        }
        else {
            throw new Error('Unknown projection name! '+name);
        }
    }
}

// Transform array or object to Proj point
function toPoint(coords) {
    if (Array.isArray(coords) && coords.length>1) {
        return new Proj.Point(coords[0], coords[1], (coords.length > 2 ? coords[2] : 0.0));
    } else if (typeof coords === 'object' &&  coords.x !== undefined && coords.y !== undefined) {
        return new Proj.Point(coords.x, coords.y, (coords.z !== undefined ? coords.z : 0.0));
    } else {
        throw new Error(`Cannot convert ${typeof coords} to point.`)
    }
}

// Transfrom any point-like object or array with z and/or m coordinates.
// Handle denormalized axis.
function transform(source, dest, point) {
    let projPoint = toPoint(point);
    var hasZ = point.z !== undefined;
    var hasM = point.m !== undefined;
    Proj.transform(source, dest, projPoint);
    let newPoint = {
      x: projPoint.x,
      y: projPoint.y
    };
    if (hasZ) {
      newPoint.z = projPoint.z;
    }
    if (hasM) {
      newPoint.m = point.m;
    }
    return newPoint;
}


// From proj4js
// Transform and return the same structural object: return an array in case input is an array
// or an object in case input is an object,
// by recopying extra fields from input to output
function transformer(from, to, coords) {
    let transformed = transform(from, to, coords);
    let geocent = from.isGeocentric || to.isGeocentric;
    if (Array.isArray(coords)) {
      if (coords.length > 2) {
          if (typeof transformed.z === 'number' && geocent) {
            return [transformed.x, transformed.y, transformed.z].concat(coords.splice(3));
          } else {
            return [transformed.x, transformed.y].concat(coords.splice(2));
          }
      } else {
        return [transformed.x, transformed.y];
      }
    } else {
      keys = Object.keys(coords);
      if (keys.length > 2) {
        keys.forEach(function (key) {
          if (key !== 'x' && key !== 'y' && (key !== 'z' || geocent)) {
            transformed[key] = coords[key];
          }
        });
      }
      return transformed;
    }
}

var wgs84 = new Proj.Projection('+proj=longlat +ellps=WGS84 +datum=WGS84 +units=degrees');


function testDef(code){
  return code in defs;
}

function getProj(item) {
    if (typeof item === 'string') {
        if (testDef(item)) {
            return defs[item];
        }
        return new Proj.Projection(item);
    }
    if (item.oProj) {
      return item.oProj;
    }
    return item
}

// This method allows for variadic arguments:
// The first argument must always be a projection
// Number of arguments:
// * 1 argument: Transformation from/to WGS84 to given projection
// * 2 arguments: two cases
//     * The second argument is a projection: return a transform object
//       from projection in arg1 to projection in arg2
//     * The second argument is coordinates (object coordinate's like or array:
//       return the result of projection from WGS84 to given projection applied
//       to the input coordinates in arg2.
// * 3 arguments:
//    * First and second arguments are respectively source and destination projections.
//      Third argument is coordinates on which to apply transformation, returns
//      the results of the transformation.
// The projection arguments could be:
// * The projection code stored in proj4.defs
// * The proj4 string
// * proj4.Proj object
//
function proj4(fromProj, toProj, coords) {
    fromProj = getProj(fromProj);
    var single = false;
    var obj;
    if (typeof toProj === 'undefined') {
      toProj = fromProj;
      fromProj = wgs84;
      single = true;
    } else if (typeof toProj.x !== 'undefined' || Array.isArray(toProj)) {
      coords = toProj;
      toProj = fromProj;
      fromProj = wgs84;
      single = true;
    }
    toProj = getProj(toProj);
    if (coords) {
      return transformer(fromProj, toProj, coords);
    } else {
      obj = {
        forward: function (coords) {
          return transformer(fromProj, toProj, coords);
        },
        inverse: function (coords) {
          return transformer(toProj, fromProj, coords);
        }
      };
      if (single) {
        obj.oProj = toProj;
      }
      return obj;
    }
}

// Init globals definition

function make_globals() {
    defs('EPSG:4326', "+title=WGS 84 (long/lat) +proj=longlat +ellps=WGS84 +datum=WGS84 +units=degrees");
    defs('EPSG:4269', "+title=NAD83 (long/lat) +proj=longlat +a=6378137.0 +b=6356752.31414036 +ellps=GRS80 +datum=NAD83 +units=degrees");
    defs('EPSG:3857', "+title=WGS 84 / Pseudo-Mercator +proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs");

    defs.WGS84 = defs['EPSG:4326'];
    defs['EPSG:3785'] = defs['EPSG:3857']; // maintain backward compat, official code is 3857
    defs.GOOGLE = defs['EPSG:3857'];
    defs['EPSG:900913'] = defs['EPSG:3857'];
    defs['EPSG:102113'] = defs['EPSG:3857'];
}

make_globals();

proj4.defs = defs;
proj4.defaultDatum = 'WGS84';
proj4.WGS84 = defs.WGS84;
proj4.toPoint = toPoint;
proj4.defs = defs;
proj4.transform = transform;

// The only method exported
export { proj4, Proj };
