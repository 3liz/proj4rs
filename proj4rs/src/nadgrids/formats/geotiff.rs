//!
//! Read grid from GeoTiff
//!
//! Inspired from the original PROJ tiff grid support handling
//!
//! Original copyright notice and permission:
//!
//! Copyright (c) 2000, Frank Warmerdam <warmerdam@pobox.com>
//! Copyright (c) 2019, Even Rouault, <even.rouault at spatialys.com>
//!
//! Permission is hereby granted, free of charge, to any person obtaining a
//! copy of this software and associated documentation files (the "Software"),
//! to deal in the Software without restriction, including without limitation
//! the rights to use, copy, modify, merge, publish, distribute, sublicense,
//! and/or sell copies of the Software, and to permit persons to whom the
//! Software is furnished to do so, subject to the following conditions:
//!

use std::io::{Read, Seek};

use tiff::decoder::{Decoder, DecodingResult, DecodingSampleType};
use tiff::tags::{self, Tag};

use crate::errors::{Error, Result};
use crate::log;
use crate::math::consts::SEC_TO_RAD;
use crate::nadgrids::grid::{Grid, GridId, Lp, REL_TOLERANCE_HGRIDSHIFT};
use crate::nadgrids::{Catalog, header::Header};

const DEG_TO_RAD: f64 = 1.0f64.to_radians();

pub(super) fn is_tiff<const N: usize>(size: usize, hdr: Header<N>) -> bool {
    // Test combinations of signature for ClassicTIFF/BigTIFF little/big endian
    size >= 4
        && (((hdr.get_u8(0) == b'I' && hdr.get_u8(1) == b'I')
            || (hdr.get_u8(0) == b'M' && hdr.get_u8(1) == b'M'))
            && ((hdr.get_u8(2) == 0x2a && hdr.get_u8(3) == 0)
                || (hdr.get_u8(3) == 0x2a && hdr.get_u8(2) == 0)
                || (hdr.get_u8(2) == 0x2b && hdr.get_u8(3) == 0)
                || (hdr.get_u8(3) == 0x2b && hdr.get_u8(2) == 0)))
}

const TIFFTAG_GEOPIXELSCALE: Tag = Tag::Unknown(33_550);
const TIFFTAG_GEOTIEPOINTS: Tag = Tag::Unknown(33_922);
const TIFFTAG_GEOTRANSMATRIX: Tag = Tag::Unknown(34_264);

/// GeoTiff reader
pub(super) fn read_tiff<R: Read + Seek>(catalog: &Catalog, key: &str, read: &mut R) -> Result<()> {
    let mut reader = Decoder::new(read)?;

    for ifd_index in 0u32.. {
        read_tiff_grid(catalog, key, &mut reader, ifd_index)?;
        if !reader.more_images() {
            break;
        }
        reader.next_image()?;
    }
    Ok(())
}

fn read_tiff_grid<R: Read + Seek>(
    catalog: &Catalog,
    key: &str,
    reader: &mut Decoder<R>,
    ifd_index: u32,
) -> Result<()> {
    log::debug!("Reading TIFF Grid index {ifd_index} a {key}");

    // Check accepted photometric value
    if let Some(tag) = reader.find_tag(Tag::PhotometricInterpretation)?
        && tag.into_u16()? != tags::PhotometricInterpretation::BlackIsZero.to_u16()
    {
        return Err(Error::InvalidTiffGridFormat(
            "Unsupported Photometric value",
        ));
    }

    // Wed don't want old JPEG
    if let Some(tag) = reader.find_tag(Tag::Compression)?
        && tag.into_u16()? == tags::CompressionMethod::JPEG.to_u16()
    {
        return Err(Error::InvalidTiffGridFormat(
            "Unsupported compression method",
        ));
    }

    // Get geokeys
    let geokeys = match reader.get_tag_u16_vec(Tag::GeoKeyDirectoryTag) {
        Err(_) => {
            return Err(Error::InvalidTiffGridFormat("No GeoKeys tag"));
        }
        Ok(v) if v.len() < 4 || v.len() % 4 != 0 => {
            return Err(Error::InvalidTiffGridFormat(
                "Wrong number of values in Geokeys tag",
            ));
        }
        Ok(v) => v,
    };

    if geokeys[0] != 1 {
        return Err(Error::InvalidTiffGridFormat(
            "Unsupported GeoTIFF major version",
        ));
    }

    // We only support GeoTIFF 1.0 and 1.1 atm
    if geokeys[1] != 1 || geokeys[2] > 1 {
        log::warn!("GeoTIFF {}.{} possibly not handled", geokeys[1], geokeys[2]);
    }

    let mut is_geographic = true;
    let mut pixel_is_area = false;

    let mut i = 4;
    while i + 3 < geokeys.len() {
        const GT_MODEL_TYPE_GEO_KEY: u16 = 1024;
        const MODEL_TYPE_PROJECTED: u16 = 1;
        const MODEL_TYPE_GEOGRAPHIC: u16 = 2;
        const GT_RASTER_TYPE_GEO_KEY: u16 = 1025;
        const RASTER_PIXEL_IS_AREA: u16 = 1;

        match geokeys[i] {
            GT_MODEL_TYPE_GEO_KEY => match geokeys[i + 3] {
                MODEL_TYPE_GEOGRAPHIC => {}
                MODEL_TYPE_PROJECTED => {
                    is_geographic = false;
                }
                _ => {
                    return Err(Error::InvalidTiffGridFormat(
                        "Unsupported GTModelTypeGeoKey",
                    ));
                }
            },
            GT_RASTER_TYPE_GEO_KEY if geokeys[i + 3] == RASTER_PIXEL_IS_AREA => {
                pixel_is_area = true;
            }
            _ => {}
        }
        i += 4;
    }

    let mut hres;
    let mut vres;
    let mut west;
    let mut north;

    if let Ok(matrix) = reader.get_tag_f64_vec(TIFFTAG_GEOTRANSMATRIX)
        && matrix.len() == 16
    {
        if matrix[1] != 0. || matrix[3] != 0. {
            return Err(Error::InvalidTiffGridFormat(
                "Rotational terms not supported in GeoTransformationMatrix",
            ));
        }
        hres = matrix[0];
        vres = -matrix[5]; // negation to simulate GeoPixelScale convention
        west = matrix[3];
        north = matrix[7];
    } else {
        match reader.get_tag_f64_vec(TIFFTAG_GEOPIXELSCALE) {
            Err(_) => {
                return Err(Error::InvalidTiffGridFormat("No GeoPixelScale tag"));
            }
            Ok(v) if v.len() != 3 => {
                return Err(Error::InvalidTiffGridFormat("Invalid GeoPixelScale size"));
            }
            Ok(v) => {
                hres = v[0];
                vres = v[1];
            }
        }

        match reader.get_tag_f64_vec(TIFFTAG_GEOTIEPOINTS) {
            Err(_) => {
                return Err(Error::InvalidTiffGridFormat("No GeoTiePoints tag"));
            }
            Ok(v) if v.len() != 6 => {
                return Err(Error::InvalidTiffGridFormat("Invalid GeoTiePoints size"));
            }
            Ok(v) => {
                west = v[3] - v[0] * hres;
                north = v[4] + v[1] * vres;
            }
        }
    }

    let (rowsize, nrows) = reader.dimensions()?;
    let (rowsize, nrows) = (rowsize as usize, nrows as usize);

    if pixel_is_area {
        west += 0.5 * hres;
        north -= 0.5 * vres;
    }

    let mut south = north - vres * ((nrows - 1) as f64);
    let mut east = west + hres * ((rowsize - 1) as f64);
    let bottom_up = vres < 0.;

    if bottom_up {
        (north, south) = (south, north);
    }

    vres = vres.abs();

    if is_geographic {
        west = west.to_radians();
        north = north.to_radians();
        east = east.to_radians();
        south = south.to_radians();
        hres = hres.to_radians();
        vres = vres.to_radians();
    }

    let layout = reader.image_buffer_layout()?;
    // Make sure that we have at least two planes
    if layout.planes < 2 {
        if ifd_index == 0 {
            return Err(Error::InvalidTiffGridFormat("No enough samples"));
        } else {
            // Skip that image
            log::debug!("Skipping IFD {ifd_index} because it has no enough samples");
            return Ok(());
        }
    }

    let mut result = layout
        .sample_type
        .and_then(|sample_type| match sample_type {
            DecodingSampleType::I16 => Some(DecodingResult::I16(vec![])),
            DecodingSampleType::U16 => Some(DecodingResult::U16(vec![])),
            DecodingSampleType::I32 => Some(DecodingResult::I32(vec![])),
            DecodingSampleType::U32 => Some(DecodingResult::U32(vec![])),
            DecodingSampleType::F32 => Some(DecodingResult::F32(vec![])),
            DecodingSampleType::F64 => Some(DecodingResult::F64(vec![])),
            _ => None,
        })
        .ok_or(Error::InvalidTiffGridFormat("Unsupported sample type"))?;

    // Read GDAL metadata
    let metadata = Metadata::read(reader)?;

    // Read grid data
    let _ = reader.read_image_to_buffer(&mut result)?;

    // Ensure everything was read
    if result.as_buffer(0).byte_len() != layout.complete_len {
        return Err(Error::InvalidTiffGridFormat("Grid size too big"));
    }

    fn to_cvs<T: Into<f64> + Copy>(
        v: Vec<T>,
        nrows: usize,
        rowsize: usize,
        bottom_up: bool,
        metadata: &Metadata,
    ) -> Vec<Lp> {
        // All non-TIFF grids have the first rows in the file being the one
        // corresponding to the southern-most row. In GeoTIFF, the convention is
        // *generally* different (when m_bottomUp == false), TIFF being an
        // image-oriented image. If m_bottomUp == true, then we had GeoTIFF hints
        // that the first row of the image is the southern-most
        let gcount = nrows * rowsize;
        let lat_offset = (metadata.lat_sample_idx as usize) * gcount;
        let lon_offset = (metadata.lon_sample_idx as usize) * gcount;

        let unit_factor = metadata.unit_factor;

        if !bottom_up {
            // Follow NTv2 convention (southern-most is first row)
            // NOTE: rows are reversed
            (0..gcount)
                .rev()
                .map(|i| Lp {
                    lam: v[i + lon_offset].into() * unit_factor,
                    phi: v[i + lat_offset].into() * unit_factor,
                })
                .collect()
        } else {
            (0..gcount)
                .map(|i| Lp {
                    lam: v[i + lon_offset].into() * unit_factor,
                    phi: v[i + lat_offset].into() * unit_factor,
                })
                .collect()
        }
    }

    let mut cvs = match result {
        DecodingResult::I16(v) => to_cvs(v, nrows, rowsize, bottom_up, &metadata),
        DecodingResult::U16(v) => to_cvs(v, nrows, rowsize, bottom_up, &metadata),
        DecodingResult::I32(v) => to_cvs(v, nrows, rowsize, bottom_up, &metadata),
        DecodingResult::U32(v) => to_cvs(v, nrows, rowsize, bottom_up, &metadata),
        DecodingResult::F32(v) => to_cvs(v, nrows, rowsize, bottom_up, &metadata),
        DecodingResult::F64(v) => to_cvs(v, nrows, rowsize, bottom_up, &metadata),
        _ => unreachable!(),
    };

    // Apply scale

    if let Some((scale, offset)) = metadata.lon_adf_scale {
        let offset = offset * metadata.unit_factor;
        for v in cvs.iter_mut() {
            v.lam = v.lam * scale + offset;
        }
    }

    if let Some((scale, offset)) = metadata.lat_adf_scale {
        let offset = offset * metadata.unit_factor;
        for v in cvs.iter_mut() {
            v.phi = v.phi * scale + offset;
        }
    }

    // Positive east
    if !metadata.positive_east {
        for v in cvs.iter_mut() {
            v.lam = -v.lam;
        }
    }

    // In proj4, rows are stored in reverse order
    if !bottom_up {
        for i in 0..nrows {
            let offs = i * rowsize;
            cvs[offs..(offs + rowsize)].reverse();
        }
    }

    // Lower left corner
    let ll = Lp {
        lam: west,
        phi: south,
    };

    // Upper rigth corner
    let ur = Lp {
        lam: east,
        phi: north,
    };

    let del = Lp {
        lam: hres,
        phi: vres,
    };

    // Size limit
    let lim = Lp {
        lam: rowsize as f64,
        phi: nrows as f64,
    };

    let epsilon = (del.lam.abs() + del.phi.abs()) * REL_TOLERANCE_HGRIDSHIFT;
    let lineage = GridId::root();

    println!("{:?}", metadata.grid_name);

    let id = metadata
        .grid_name
        .unwrap_or_else(|| GridId::from((0u32, ifd_index)));

    catalog.add_grid(
        key.into(),
        Grid {
            id,
            lineage,
            ll,
            ur,
            del,
            lim,
            epsilon,
            cvs: cvs.into_boxed_slice(),
            is_geographic,
        },
    )
}

// =================
// GDAL Metadata
// =================

struct SampleMeta {
    t: Descr,
    offset: u16,
    adf_scale: Option<f64>,
    adf_offset: Option<f64>,
    positive_east: bool,
}

impl SampleMeta {
    fn new(t: Descr, offset: u16) -> Self {
        Self {
            t,
            offset,
            adf_scale: None,
            adf_offset: None,
            positive_east: true,
        }
    }

    fn adf_scale(&self) -> Option<(f64, f64)> {
        match (self.adf_scale, self.adf_offset) {
            (None, None) => None,
            (Some(scale), None) => Some((scale, 0.)),
            (None, Some(offset)) => Some((1.0, offset)),
            (Some(scale), Some(offset)) => Some((scale, offset)),
        }
    }
}

#[derive(Default)]
struct Metadata {
    grid_name: Option<GridId>,
    unit_factor: f64,
    lon_sample_idx: u16,
    lat_sample_idx: u16,
    lon_adf_scale: Option<(f64, f64)>,
    lat_adf_scale: Option<(f64, f64)>,
    positive_east: bool,
}

#[derive(PartialEq, Copy, Clone)]
enum UnitType {
    ArcSecond,
    Radian,
    Degree,
}

#[derive(PartialEq)]
enum Descr {
    Lat,
    Lon,
}

impl UnitType {
    fn from_str(s: &str) -> Result<Self> {
        match s {
            "arc-second" | "arc-second per year" => Ok(Self::ArcSecond),
            "radian" => Ok(Self::Radian),
            "degree" => Ok(Self::Degree),
            _ => Err(Error::InvalidTiffGridFormat("Unsupported unit")),
        }
    }
}

impl Metadata {
    const TIFFTAG_GDAL_METADATA: Tag = Tag::Unknown(42112);

    fn parse_attributs(attr_str: &str) -> Option<(&str, Option<&str>, Option<&str>)> {
        let mut sample = None;
        let mut role = None;
        let mut name = None;
        for attr in attr_str.split_ascii_whitespace() {
            if let Some((attr, value)) = attr.split_once("=") {
                let value = value.trim_matches('\"');
                match attr {
                    "name" => name = Some(value),
                    "sample" => sample = Some(value),
                    "role" => role = Some(value),
                    _ => continue,
                }
            }
        }
        name.map(|name| (name, sample, role))
    }

    fn read<R: Read + Seek>(reader: &mut Decoder<R>) -> Result<Self> {
        Self::parse(
            &reader
                .get_tag_ascii_string(Self::TIFFTAG_GDAL_METADATA)
                .unwrap_or_default()
                .trim(),
        )
    }

    fn parse(md: &str) -> Result<Self> {
        let mut units: Option<UnitType> = None;
        let mut grid_name = None;
        let mut sample0 = SampleMeta::new(Descr::Lat, 0);
        let mut sample1 = SampleMeta::new(Descr::Lon, 1);

        fn parse_f64(s: &str) -> Result<f64> {
            use crate::parse::FromStr;
            f64::from_str(s).map_err(|_| Error::InvalidTiffGridFormat("f64 Parse error"))
        }

        if !md.is_empty() {
            log::debug!("Found GDAL Metadata '{md}'");

            // GDAL metadata is XML format
            md.strip_prefix("<GDALMetadata>")
                .and_then(|s| s.strip_suffix("</GDALMetadata>"))
                .unwrap_or_default()
                .split("<Item ")
                .try_for_each(|item| {
                    if let Some((attrs, value)) = item
                        .trim()
                        .strip_suffix("</Item>")
                        .and_then(|s| s.split_once(">"))
                    {
                        let value = value.trim();
                        match Self::parse_attributs(attrs) {
                            Some(("UNITTYPE", Some(sample), Some("unittype"))) => match sample {
                                "0" | "1" => {
                                    let unittype = UnitType::from_str(value)?;
                                    match units {
                                        Some(units) => {
                                            if units != unittype {
                                                return Err(Error::InvalidTiffGridFormat(
                                                    "Samples have different units",
                                                ));
                                            }
                                        }
                                        None => units = Some(unittype),
                                    }
                                }
                                _ => {}
                            },
                            Some(("DESCRIPTION", Some(sample), Some("description"))) => {
                                match sample {
                                    "0" => match value {
                                        "latitude_offset" => sample0.t = Descr::Lat,
                                        "longitude_offset" => sample0.t = Descr::Lon,
                                        _ => {
                                            return Err(Error::InvalidTiffGridFormat(
                                                "Unexpected offset type for sample 0",
                                            ));
                                        }
                                    },
                                    "1" => match value {
                                        "latitude_offset" => sample1.t = Descr::Lat,
                                        "longitude_offset" => sample1.t = Descr::Lon,
                                        _ => {
                                            return Err(Error::InvalidTiffGridFormat(
                                                "Unexpected offset type for sample 1",
                                            ));
                                        }
                                    },
                                    _ => {}
                                }
                            }
                            Some(("positive_value", Some(sample), _)) => match sample {
                                "0" => sample0.positive_east = value == "east",
                                "1" => sample1.positive_east = value == "east",
                                _ => {}
                            },
                            Some((_, Some(sample), Some("scale"))) => match sample {
                                "0" => sample0.adf_scale = Some(parse_f64(value)?),
                                "1" => sample1.adf_scale = Some(parse_f64(value)?),
                                _ => {}
                            },
                            Some((_, Some(sample), Some("offset"))) => match sample {
                                "0" => sample0.adf_offset = Some(parse_f64(value)?),
                                "1" => sample1.adf_offset = Some(parse_f64(value)?),
                                _ => {}
                            },
                            Some(("grid_name", _, _)) => {
                                grid_name = Some(GridId::from(value.as_bytes()));
                            }
                            _ => {}
                        }
                    }
                    Ok(())
                })?
        }

        if sample0.t == sample1.t {
            return Err(Error::InvalidTiffGridFormat("Incoherent sample indexes"));
        }

        let (lon_sample, lat_sample) = if sample0.t == Descr::Lon {
            (sample0, sample1)
        } else {
            (sample1, sample0)
        };

        Ok(Self {
            grid_name,
            unit_factor: units
                .map(|units| match units {
                    UnitType::ArcSecond => SEC_TO_RAD,
                    UnitType::Radian => 1.,
                    UnitType::Degree => DEG_TO_RAD,
                })
                .unwrap_or(SEC_TO_RAD),
            lon_sample_idx: lon_sample.offset,
            lat_sample_idx: lat_sample.offset,
            lon_adf_scale: lon_sample.adf_scale(),
            lat_adf_scale: lat_sample.adf_scale(),
            positive_east: lon_sample.positive_east,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nadgrids::Catalog;
    use crate::tests::setup;
    use std::env;
    use std::fs::File;
    use std::io::BufReader;
    use std::path::Path;

    macro_rules! fixture {
        ($name:expr) => {
            Path::new(&env::var("CARGO_MANIFEST_DIR").unwrap())
                .join("fixtures")
                .as_path()
                .join($name)
                .as_path()
        };
    }

    macro_rules! load_tiff {
        ($cat:expr, $name:expr) => {
            // Use a BufReader or efficiency
            let file = File::open(fixture!($name)).unwrap();
            let mut read = BufReader::new(file);
            read_tiff($cat, $name, &mut read).unwrap();
        };
    }

    #[test]
    fn tiff_gdal_metadata() {
        let md = Metadata::parse(concat!(
            r#"<GDALMetadata>"#,
            r#"<Item name="area_of_use">Spain - Catalonia</Item>"#,
            r#"<Item name="grid_name">0INT2GRS</Item>"#,
            r#"<Item name="target_crs_epsg_code">4258</Item>"#,
            r#"<Item name="TYPE">HORIZONTAL_OFFSET</Item>"#,
            r#"<Item name="UNITTYPE" sample="0" role="unittype">arc-second</Item>"#,
            r#"<Item name="DESCRIPTION" sample="0" role="description">latitude_offset</Item>"#,
            r#"<Item name="positive_value" sample="1">east</Item>"#,
            r#"<Item name="UNITTYPE" sample="1" role="unittype">arc-second</Item>"#,
            r#"<Item name="DESCRIPTION" sample="1" role="description">longitude_offset</Item>"#,
            r#"</GDALMetadata>"#,
        ))
        .expect("Failed to parse GDAL metadata");

        assert_eq!(md.grid_name, Some(GridId::from("0INT2GRS".as_bytes())));
        assert_eq!(md.grid_name.unwrap().as_str(), "0INT2GRS");
        assert_eq!(md.unit_factor, SEC_TO_RAD);
        assert_eq!(md.lon_sample_idx, 1);
        assert_eq!(md.lat_sample_idx, 0);
        assert_eq!(md.lon_adf_scale, None);
        assert_eq!(md.lat_adf_scale, None);
        assert!(md.positive_east);
    }

    #[test]
    fn tiff_100800401() {
        setup();

        let catalog = Catalog::default();
        load_tiff!(&catalog, "es_cat_icgc_100800401.tif");

        let grids = catalog
            .find("es_cat_icgc_100800401.tif")
            .unwrap()
            .collect::<Vec<_>>();
        assert_eq!(grids.len(), 1);

        let grid = grids[0];
        assert!(grid.is_root());
        //assert_eq!(grid.id.as_str(), "0INT2GRS");
        assert_eq!(grid.cvs.len(), 1591);
    }
}
