#define FFI_SCOPE "PROJ4RS"
#define FFI_LIB "libproj4rs_c.so"

// FFI declarations for proj4rs, usable with FFI::load() and opcache.preload.
// The FFI_* defines must stay at the very top of the file.
// Keep in sync with proj4rs-clib/src/lib.rs

typedef struct Proj4rs Proj4rs;

const char *proj4rs_last_error(void);
Proj4rs *proj4rs_proj_new(const char *c_defn);
void proj4rs_proj_delete(Proj4rs *c_ptr);
const char *proj4rs_proj_projname(const Proj4rs *c_ptr);
bool proj4rs_proj_is_latlong(const Proj4rs *c_ptr);
bool proj4rs_proj_is_geocent(const Proj4rs *c_ptr);
const uint8_t *proj4rs_proj_axis(const Proj4rs *c_ptr);
bool proj4rs_proj_is_normalized_axis(const Proj4rs *c_ptr);
double proj4rs_proj_to_meter(const Proj4rs *c_ptr);
const char *proj4rs_proj_units(const Proj4rs *c_ptr);
int proj4rs_transform(const Proj4rs *src_ptr,
                      const Proj4rs *dst_ptr,
                      double *x,
                      double *y,
                      double *z,
                      ssize_t len,
                      ssize_t stride,
                      bool convert);
