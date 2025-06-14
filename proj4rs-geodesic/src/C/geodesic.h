/**
 * \file geodesic.h
 * \brief API for the geodesic routines in C
 *
 * These routines are a simple transcription of the corresponding C++ classes
 * in <a href="https://geographiclib.sourceforge.io"> GeographicLib</a>.  The
 * "class data" is represented by the structs geod_geodesic, geod_geodesicline,
 * geod_polygon and pointers to these objects are passed as initial arguments
 * to the member functions.  Most of the internal comments have been retained.
 * However, in the process of transcription some documentation has been lost
 * and the documentation for the C++ classes, GeographicLib::Geodesic,
 * GeographicLib::GeodesicLine, and GeographicLib::PolygonAreaT, should be
 * consulted.  The C++ code remains the "reference implementation".  Think
 * twice about restructuring the internals of the C code since this may make
 * porting fixes from the C++ code more difficult.
 *
 * Copyright (c) Charles Karney (2012-2022) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEODESIC_H)
#define GEODESIC_H 1

/**
 * The major version of the geodesic library.  (This tracks the version of
 * GeographicLib.)
 **********************************************************************/
#define GEODESIC_VERSION_MAJOR 2
/**
 * The minor version of the geodesic library.  (This tracks the version of
 * GeographicLib.)
 **********************************************************************/
#define GEODESIC_VERSION_MINOR 1
/**
 * The patch level of the geodesic library.  (This tracks the version of
 * GeographicLib.)
 **********************************************************************/
#define GEODESIC_VERSION_PATCH 0

/**
 * Pack the version components into a single integer.  Users should not rely on
 * this particular packing of the components of the version number; see the
 * documentation for ::GEODESIC_VERSION, below.
 **********************************************************************/
#define GEODESIC_VERSION_NUM(a,b,c) ((((a) * 10000 + (b)) * 100) + (c))

/**
 * The version of the geodesic library as a single integer, packed as MMmmmmpp
 * where MM is the major version, mmmm is the minor version, and pp is the
 * patch level.  Users should not rely on this particular packing of the
 * components of the version number.  Instead they should use a test such as
 * @code{.c}
   #if GEODESIC_VERSION >= GEODESIC_VERSION_NUM(1,40,0)
   ...
   #endif
 * @endcode
 **********************************************************************/
#define GEODESIC_VERSION \
 GEODESIC_VERSION_NUM(GEODESIC_VERSION_MAJOR, \
                      GEODESIC_VERSION_MINOR, \
                      GEODESIC_VERSION_PATCH)

#define GEOD_DLL

  /**
   * The struct containing information about the ellipsoid.  This must be
   * initialized by geod_init() before use.
   **********************************************************************/
  struct geod_geodesic {
    double a;                   /**< the equatorial radius */
    double f;                   /**< the flattening */
    /**< @cond SKIP */
    double f1, e2, ep2, n, b, c2, etol2;
    double A3x[6], C3x[15], C4x[21];
    /**< @endcond */
  };

  /**
   * The struct containing information about a single geodesic.  This must be
   * initialized by geod_lineinit(), geod_directline(), geod_gendirectline(),
   * or geod_inverseline() before use.
   **********************************************************************/
  struct geod_geodesicline {
    double lat1;                /**< the starting latitude */
    double lon1;                /**< the starting longitude */
    double azi1;                /**< the starting azimuth */
    double a;                   /**< the equatorial radius */
    double f;                   /**< the flattening */
    double salp1;               /**< sine of \e azi1 */
    double calp1;               /**< cosine of \e azi1 */
    double a13;                 /**< arc length to reference point */
    double s13;                 /**< distance to reference point */
    /**< @cond SKIP */
    double b, c2, f1, salp0, calp0, k2,
      ssig1, csig1, dn1, stau1, ctau1, somg1, comg1,
      A1m1, A2m1, A3c, B11, B21, B31, A4, B41;
    double C1a[6+1], C1pa[6+1], C2a[6+1], C3a[6], C4a[6];
    /**< @endcond */
    unsigned caps;              /**< the capabilities */
  };

  /**
   * The struct for accumulating information about a geodesic polygon.  This is
   * used for computing the perimeter and area of a polygon.  This must be
   * initialized by geod_polygon_init() before use.
   **********************************************************************/
  struct geod_polygon {
    double lat;                 /**< the current latitude */
    double lon;                 /**< the current longitude */
    /**< @cond SKIP */
    double lat0;
    double lon0;
    double A[2];
    double P[2];
    int polyline;
    int crossings;
    /**< @endcond */
    unsigned num;               /**< the number of points so far */
  };

  /**
   * Initialize the library
   */
  void GEOD_DLL Init();

  /**
   * Initialize a geod_geodesic object.
   *
   * @param[out] g a pointer to the object to be initialized.
   * @param[in] a the equatorial radius (meters).
   * @param[in] f the flattening.
   **********************************************************************/
  bool GEOD_DLL geod_init(struct geod_geodesic* g, double a, double f);

  /**
   * Solve the direct geodesic problem.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] lat1 latitude of point 1 (degrees).
   * @param[in] lon1 longitude of point 1 (degrees).
   * @param[in] azi1 azimuth at point 1 (degrees).
   * @param[in] s12 distance from point 1 to point 2 (meters); it can be
   *   negative.
   * @param[out] plat2 pointer to the latitude of point 2 (degrees).
   * @param[out] plon2 pointer to the longitude of point 2 (degrees).
   * @param[out] pazi2 pointer to the (forward) azimuth at point 2 (degrees).
   *
   * \e g must have been initialized with a call to geod_init().  \e lat1
   * should be in the range [&minus;90&deg;, 90&deg;].  The values of \e lon2
   * and \e azi2 returned are in the range [&minus;180&deg;, 180&deg;].  Any of
   * the "return" arguments \e plat2, etc., may be replaced by 0, if you do not
   * need some quantities computed.
   *
   * If either point is at a pole, the azimuth is defined by keeping the
   * longitude fixed, writing \e lat = &plusmn;(90&deg; &minus; &epsilon;), and
   * taking the limit &epsilon; &rarr; 0+.  An arc length greater that 180&deg;
   * signifies a geodesic which is not a shortest path.  (For a prolate
   * ellipsoid, an additional condition is necessary for a shortest path: the
   * longitudinal extent must not exceed of 180&deg;.)
   *
   * Example, determine the point 10000 km NE of JFK:
   @code{.c}
   struct geod_geodesic g;
   double lat, lon;
   geod_init(&g, 6378137, 1/298.257223563);
   geod_direct(&g, 40.64, -73.78, 45.0, 10e6, &lat, &lon, 0);
   printf("%.5f %.5f\n", lat, lon);
   @endcode
   **********************************************************************/
  void GEOD_DLL geod_direct(const struct geod_geodesic* g,
                            double lat1, double lon1, double azi1, double s12,
                            double* plat2, double* plon2, double* pazi2);

  /**
   * The general direct geodesic problem.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] lat1 latitude of point 1 (degrees).
   * @param[in] lon1 longitude of point 1 (degrees).
   * @param[in] azi1 azimuth at point 1 (degrees).
   * @param[in] flags bitor'ed combination of ::geod_flags; \e flags &
   *   ::GEOD_ARCMODE determines the meaning of \e s12_a12 and \e flags &
   *   ::GEOD_LONG_UNROLL "unrolls" \e lon2.
   * @param[in] s12_a12 if \e flags & ::GEOD_ARCMODE is 0, this is the distance
   *   from point 1 to point 2 (meters); otherwise it is the arc length
   *   from point 1 to point 2 (degrees); it can be negative.
   * @param[out] plat2 pointer to the latitude of point 2 (degrees).
   * @param[out] plon2 pointer to the longitude of point 2 (degrees).
   * @param[out] pazi2 pointer to the (forward) azimuth at point 2 (degrees).
   * @param[out] ps12 pointer to the distance from point 1 to point 2
   *   (meters).
   * @param[out] pm12 pointer to the reduced length of geodesic (meters).
   * @param[out] pM12 pointer to the geodesic scale of point 2 relative to
   *   point 1 (dimensionless).
   * @param[out] pM21 pointer to the geodesic scale of point 1 relative to
   *   point 2 (dimensionless).
   * @param[out] pS12 pointer to the area under the geodesic
   *   (meters<sup>2</sup>).
   * @return \e a12 arc length from point 1 to point 2 (degrees).
   *
   * \e g must have been initialized with a call to geod_init().  \e lat1
   * should be in the range [&minus;90&deg;, 90&deg;].  The function value \e
   * a12 equals \e s12_a12 if \e flags & ::GEOD_ARCMODE.  Any of the "return"
   * arguments, \e plat2, etc., may be replaced by 0, if you do not need some
   * quantities computed.
   *
   * With \e flags & ::GEOD_LONG_UNROLL bit set, the longitude is "unrolled" so
   * that the quantity \e lon2 &minus; \e lon1 indicates how many times and in
   * what sense the geodesic encircles the ellipsoid.
   **********************************************************************/
  double GEOD_DLL geod_gendirect(const struct geod_geodesic* g,
                                 double lat1, double lon1, double azi1,
                                 unsigned flags, double s12_a12,
                                 double* plat2, double* plon2, double* pazi2,
                                 double* ps12, double* pm12,
                                 double* pM12, double* pM21,
                                 double* pS12);

  /**
   * Solve the inverse geodesic problem.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] lat1 latitude of point 1 (degrees).
   * @param[in] lon1 longitude of point 1 (degrees).
   * @param[in] lat2 latitude of point 2 (degrees).
   * @param[in] lon2 longitude of point 2 (degrees).
   * @param[out] ps12 pointer to the distance from point 1 to point 2
   *   (meters).
   * @param[out] pazi1 pointer to the azimuth at point 1 (degrees).
   * @param[out] pazi2 pointer to the (forward) azimuth at point 2 (degrees).
   *
   * @return \e a12 arc length from point 1 to point 2 (degrees).
   *
   * \e g must have been initialized with a call to geod_init().  \e lat1 and
   * \e lat2 should be in the range [&minus;90&deg;, 90&deg;].  The values of
   * \e azi1 and \e azi2 returned are in the range [&minus;180&deg;, 180&deg;].
   * Any of the "return" arguments, \e ps12, etc., may be replaced by 0, if you
   * do not need some quantities computed.
   *
   * If either point is at a pole, the azimuth is defined by keeping the
   * longitude fixed, writing \e lat = &plusmn;(90&deg; &minus; &epsilon;), and
   * taking the limit &epsilon; &rarr; 0+.
   *
   * The solution to the inverse problem is found using Newton's method.  If
   * this fails to converge (this is very unlikely in geodetic applications
   * but does occur for very eccentric ellipsoids), then the bisection method
   * is used to refine the solution.
   *
   * Example, determine the distance between JFK and Singapore Changi Airport:
   @code{.c}
   struct geod_geodesic g;
   double s12;
   geod_init(&g, 6378137, 1/298.257223563);
   geod_inverse(&g, 40.64, -73.78, 1.36, 103.99, &s12, 0, 0);
   printf("%.3f\n", s12);
   @endcode
   **********************************************************************/
   void GEOD_DLL geod_inverse(const struct geod_geodesic* g,
                              double lat1, double lon1,
                              double lat2, double lon2,
                              double* ps12, double* pazi1, double* pazi2);

  /**
   * The general inverse geodesic calculation.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] lat1 latitude of point 1 (degrees).
   * @param[in] lon1 longitude of point 1 (degrees).
   * @param[in] lat2 latitude of point 2 (degrees).
   * @param[in] lon2 longitude of point 2 (degrees).
   * @param[out] ps12 pointer to the distance from point 1 to point 2
   *  (meters).
   * @param[out] pazi1 pointer to the azimuth at point 1 (degrees).
   * @param[out] pazi2 pointer to the (forward) azimuth at point 2 (degrees).
   * @param[out] pm12 pointer to the reduced length of geodesic (meters).
   * @param[out] pM12 pointer to the geodesic scale of point 2 relative to
   *   point 1 (dimensionless).
   * @param[out] pM21 pointer to the geodesic scale of point 1 relative to
   *   point 2 (dimensionless).
   * @param[out] pS12 pointer to the area under the geodesic
   *   (meters<sup>2</sup>).
   * @return \e a12 arc length from point 1 to point 2 (degrees).
   *
   * \e g must have been initialized with a call to geod_init().  \e lat1 and
   * \e lat2 should be in the range [&minus;90&deg;, 90&deg;].  Any of the
   * "return" arguments \e ps12, etc., may be replaced by 0, if you do not need
   * some quantities computed.
   **********************************************************************/
  double GEOD_DLL geod_geninverse(const struct geod_geodesic* g,
                                  double lat1, double lon1,
                                  double lat2, double lon2,
                                  double* ps12, double* pazi1, double* pazi2,
                                  double* pm12, double* pM12, double* pM21,
                                  double* pS12);

  /**
   * Initialize a geod_geodesicline object.
   *
   * @param[out] l a pointer to the object to be initialized.
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] lat1 latitude of point 1 (degrees).
   * @param[in] lon1 longitude of point 1 (degrees).
   * @param[in] azi1 azimuth at point 1 (degrees).
   * @param[in] caps bitor'ed combination of ::geod_mask values specifying the
   *   capabilities the geod_geodesicline object should possess, i.e., which
   *   quantities can be returned in calls to geod_position() and
   *   geod_genposition().
   *
   * \e g must have been initialized with a call to geod_init().  \e lat1
   * should be in the range [&minus;90&deg;, 90&deg;].
   *
   * The ::geod_mask values are:
   * - \e caps |= ::GEOD_LATITUDE for the latitude \e lat2; this is
   *   added automatically,
   * - \e caps |= ::GEOD_LONGITUDE for the latitude \e lon2,
   * - \e caps |= ::GEOD_AZIMUTH for the latitude \e azi2; this is
   *   added automatically,
   * - \e caps |= ::GEOD_DISTANCE for the distance \e s12,
   * - \e caps |= ::GEOD_REDUCEDLENGTH for the reduced length \e m12,
   * - \e caps |= ::GEOD_GEODESICSCALE for the geodesic scales \e M12
   *   and \e M21,
   * - \e caps |= ::GEOD_AREA for the area \e S12,
   * - \e caps |= ::GEOD_DISTANCE_IN permits the length of the
   *   geodesic to be given in terms of \e s12; without this capability the
   *   length can only be specified in terms of arc length.
   * .
   * A value of \e caps = 0 is treated as ::GEOD_LATITUDE | ::GEOD_LONGITUDE |
   * ::GEOD_AZIMUTH | ::GEOD_DISTANCE_IN (to support the solution of the
   * "standard" direct problem).
   *
   * When initialized by this function, point 3 is undefined (l->s13 = l->a13 =
   * NaN).
   **********************************************************************/
  void GEOD_DLL geod_lineinit(struct geod_geodesicline* l,
                              const struct geod_geodesic* g,
                              double lat1, double lon1, double azi1,
                              unsigned caps);

  /**
   * Initialize a geod_geodesicline object in terms of the direct geodesic
   * problem.
   *
   * @param[out] l a pointer to the object to be initialized.
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] lat1 latitude of point 1 (degrees).
   * @param[in] lon1 longitude of point 1 (degrees).
   * @param[in] azi1 azimuth at point 1 (degrees).
   * @param[in] s12 distance from point 1 to point 2 (meters); it can be
   *   negative.
   * @param[in] caps bitor'ed combination of ::geod_mask values specifying the
   *   capabilities the geod_geodesicline object should possess, i.e., which
   *   quantities can be returned in calls to geod_position() and
   *   geod_genposition().
   *
   * This function sets point 3 of the geod_geodesicline to correspond to point
   * 2 of the direct geodesic problem.  See geod_lineinit() for more
   * information.
   **********************************************************************/
  void GEOD_DLL geod_directline(struct geod_geodesicline* l,
                                const struct geod_geodesic* g,
                                double lat1, double lon1,
                                double azi1, double s12,
                                unsigned caps);

  /**
   * Initialize a geod_geodesicline object in terms of the direct geodesic
   * problem specified in terms of either distance or arc length.
   *
   * @param[out] l a pointer to the object to be initialized.
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] lat1 latitude of point 1 (degrees).
   * @param[in] lon1 longitude of point 1 (degrees).
   * @param[in] azi1 azimuth at point 1 (degrees).
   * @param[in] flags either ::GEOD_NOFLAGS or ::GEOD_ARCMODE to determining
   *   the meaning of the \e s12_a12.
   * @param[in] s12_a12 if \e flags = ::GEOD_NOFLAGS, this is the distance
   *   from point 1 to point 2 (meters); if \e flags = ::GEOD_ARCMODE, it is
   *   the arc length from point 1 to point 2 (degrees); it can be
   *   negative.
   * @param[in] caps bitor'ed combination of ::geod_mask values specifying the
   *   capabilities the geod_geodesicline object should possess, i.e., which
   *   quantities can be returned in calls to geod_position() and
   *   geod_genposition().
   *
   * This function sets point 3 of the geod_geodesicline to correspond to point
   * 2 of the direct geodesic problem.  See geod_lineinit() for more
   * information.
   **********************************************************************/
  void GEOD_DLL geod_gendirectline(struct geod_geodesicline* l,
                                   const struct geod_geodesic* g,
                                   double lat1, double lon1, double azi1,
                                   unsigned flags, double s12_a12,
                                   unsigned caps);

  /**
   * Initialize a geod_geodesicline object in terms of the inverse geodesic
   * problem.
   *
   * @param[out] l a pointer to the object to be initialized.
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] lat1 latitude of point 1 (degrees).
   * @param[in] lon1 longitude of point 1 (degrees).
   * @param[in] lat2 latitude of point 2 (degrees).
   * @param[in] lon2 longitude of point 2 (degrees).
   * @param[in] caps bitor'ed combination of ::geod_mask values specifying the
   *   capabilities the geod_geodesicline object should possess, i.e., which
   *   quantities can be returned in calls to geod_position() and
   *   geod_genposition().
   *
   * This function sets point 3 of the geod_geodesicline to correspond to point
   * 2 of the inverse geodesic problem.  See geod_lineinit() for more
   * information.
   **********************************************************************/
  void GEOD_DLL geod_inverseline(struct geod_geodesicline* l,
                                 const struct geod_geodesic* g,
                                 double lat1, double lon1,
                                 double lat2, double lon2,
                                 unsigned caps);

  /**
   * Compute the position along a geod_geodesicline.
   *
   * @param[in] l a pointer to the geod_geodesicline object specifying the
   *   geodesic line.
   * @param[in] s12 distance from point 1 to point 2 (meters); it can be
   *   negative.
   * @param[out] plat2 pointer to the latitude of point 2 (degrees).
   * @param[out] plon2 pointer to the longitude of point 2 (degrees); requires
   *   that \e l was initialized with \e caps |= ::GEOD_LONGITUDE.
   * @param[out] pazi2 pointer to the (forward) azimuth at point 2 (degrees).
   *
   * \e l must have been initialized with a call, e.g., to geod_lineinit(),
   * with \e caps |= ::GEOD_DISTANCE_IN (or \e caps = 0).  The values of \e
   * lon2 and \e azi2 returned are in the range [&minus;180&deg;, 180&deg;].
   * Any of the "return" arguments \e plat2, etc., may be replaced by 0, if you
   * do not need some quantities computed.
   *
   * Example, compute way points between JFK and Singapore Changi Airport
   * the "obvious" way using geod_direct():
   @code{.c}
   struct geod_geodesic g;
   double s12, azi1, lat[101], lon[101];
   int i;
   geod_init(&g, 6378137, 1/298.257223563);
   geod_inverse(&g, 40.64, -73.78, 1.36, 103.99, &s12, &azi1, 0);
   for (i = 0; i < 101; ++i) {
     geod_direct(&g, 40.64, -73.78, azi1, i * s12 * 0.01, lat + i, lon + i, 0);
     printf("%.5f %.5f\n", lat[i], lon[i]);
   }
   @endcode
   * A faster way using geod_position():
   @code{.c}
   struct geod_geodesic g;
   struct geod_geodesicline l;
   double lat[101], lon[101];
   int i;
   geod_init(&g, 6378137, 1/298.257223563);
   geod_inverseline(&l, &g, 40.64, -73.78, 1.36, 103.99, 0);
   for (i = 0; i <= 100; ++i) {
     geod_position(&l, i * l.s13 * 0.01, lat + i, lon + i, 0);
     printf("%.5f %.5f\n", lat[i], lon[i]);
   }
   @endcode
   **********************************************************************/
  void GEOD_DLL geod_position(const struct geod_geodesicline* l, double s12,
                              double* plat2, double* plon2, double* pazi2);

  /**
   * The general position function.
   *
   * @param[in] l a pointer to the geod_geodesicline object specifying the
   *   geodesic line.
   * @param[in] flags bitor'ed combination of ::geod_flags; \e flags &
   *   ::GEOD_ARCMODE determines the meaning of \e s12_a12 and \e flags &
   *   ::GEOD_LONG_UNROLL "unrolls" \e lon2; if \e flags & ::GEOD_ARCMODE is 0,
   *   then \e l must have been initialized with \e caps |= ::GEOD_DISTANCE_IN.
   * @param[in] s12_a12 if \e flags & ::GEOD_ARCMODE is 0, this is the
   *   distance from point 1 to point 2 (meters); otherwise it is the
   *   arc length from point 1 to point 2 (degrees); it can be
   *   negative.
   * @param[out] plat2 pointer to the latitude of point 2 (degrees).
   * @param[out] plon2 pointer to the longitude of point 2 (degrees); requires
   *   that \e l was initialized with \e caps |= ::GEOD_LONGITUDE.
   * @param[out] pazi2 pointer to the (forward) azimuth at point 2 (degrees).
   * @param[out] ps12 pointer to the distance from point 1 to point 2
   *   (meters); requires that \e l was initialized with \e caps |=
   *   ::GEOD_DISTANCE.
   * @param[out] pm12 pointer to the reduced length of geodesic (meters);
   *   requires that \e l was initialized with \e caps |= ::GEOD_REDUCEDLENGTH.
   * @param[out] pM12 pointer to the geodesic scale of point 2 relative to
   *   point 1 (dimensionless); requires that \e l was initialized with \e caps
   *   |= ::GEOD_GEODESICSCALE.
   * @param[out] pM21 pointer to the geodesic scale of point 1 relative to
   *   point 2 (dimensionless); requires that \e l was initialized with \e caps
   *   |= ::GEOD_GEODESICSCALE.
   * @param[out] pS12 pointer to the area under the geodesic
   *   (meters<sup>2</sup>); requires that \e l was initialized with \e caps |=
   *   ::GEOD_AREA.
   * @return \e a12 arc length from point 1 to point 2 (degrees).
   *
   * \e l must have been initialized with a call to geod_lineinit() with \e
   * caps |= ::GEOD_DISTANCE_IN.  The value \e azi2 returned is in the range
   * [&minus;180&deg;, 180&deg;].  Any of the "return" arguments \e plat2,
   * etc., may be replaced by 0, if you do not need some quantities
   * computed.  Requesting a value which \e l is not capable of computing
   * is not an error; the corresponding argument will not be altered.
   *
   * With \e flags & ::GEOD_LONG_UNROLL bit set, the longitude is "unrolled" so
   * that the quantity \e lon2 &minus; \e lon1 indicates how many times and in
   * what sense the geodesic encircles the ellipsoid.
   *
   * Example, compute way points between JFK and Singapore Changi Airport using
   * geod_genposition().  In this example, the points are evenly spaced in arc
   * length (and so only approximately equally spaced in distance).  This is
   * faster than using geod_position() and would be appropriate if drawing the
   * path on a map.
   @code{.c}
   struct geod_geodesic g;
   struct geod_geodesicline l;
   double lat[101], lon[101];
   int i;
   geod_init(&g, 6378137, 1/298.257223563);
   geod_inverseline(&l, &g, 40.64, -73.78, 1.36, 103.99,
                    GEOD_LATITUDE | GEOD_LONGITUDE);
   for (i = 0; i <= 100; ++i) {
     geod_genposition(&l, GEOD_ARCMODE, i * l.a13 * 0.01,
                      lat + i, lon + i, 0, 0, 0, 0, 0, 0);
     printf("%.5f %.5f\n", lat[i], lon[i]);
   }
   @endcode
   **********************************************************************/
  double GEOD_DLL geod_genposition(const struct geod_geodesicline* l,
                                   unsigned flags, double s12_a12,
                                   double* plat2, double* plon2, double* pazi2,
                                   double* ps12, double* pm12,
                                   double* pM12, double* pM21,
                                   double* pS12);

  /**
   * Specify position of point 3 in terms of distance.
   *
   * @param[in,out] l a pointer to the geod_geodesicline object.
   * @param[in] s13 the distance from point 1 to point 3 (meters); it
   *   can be negative.
   *
   * This is only useful if the geod_geodesicline object has been constructed
   * with \e caps |= ::GEOD_DISTANCE_IN.
   **********************************************************************/
  void GEOD_DLL geod_setdistance(struct geod_geodesicline* l, double s13);

  /**
   * Specify position of point 3 in terms of either distance or arc length.
   *
   * @param[in,out] l a pointer to the geod_geodesicline object.
   * @param[in] flags either ::GEOD_NOFLAGS or ::GEOD_ARCMODE to determining
   *   the meaning of the \e s13_a13.
   * @param[in] s13_a13 if \e flags = ::GEOD_NOFLAGS, this is the distance
   *   from point 1 to point 3 (meters); if \e flags = ::GEOD_ARCMODE, it is
   *   the arc length from point 1 to point 3 (degrees); it can be
   *   negative.
   *
   * If flags = ::GEOD_NOFLAGS, this calls geod_setdistance().  If flags =
   * ::GEOD_ARCMODE, the \e s13 is only set if the geod_geodesicline object has
   * been constructed with \e caps |= ::GEOD_DISTANCE.
   **********************************************************************/
  void GEOD_DLL geod_gensetdistance(struct geod_geodesicline* l,
                                    unsigned flags, double s13_a13);

  /**
   * Initialize a geod_polygon object.
   *
   * @param[out] p a pointer to the object to be initialized.
   * @param[in] polylinep non-zero if a polyline instead of a polygon.
   *
   * If \e polylinep is zero, then the sequence of vertices and edges added by
   * geod_polygon_addpoint() and geod_polygon_addedge() define a polygon and
   * the perimeter and area are returned by geod_polygon_compute().  If \e
   * polylinep is non-zero, then the vertices and edges define a polyline and
   * only the perimeter is returned by geod_polygon_compute().
   *
   * The area and perimeter are accumulated at two times the standard floating
   * point precision to guard against the loss of accuracy with many-sided
   * polygons.  At any point you can ask for the perimeter and area so far.
   *
   * An example of the use of this function is given in the documentation for
   * geod_polygon_compute().
   **********************************************************************/
  void GEOD_DLL geod_polygon_init(struct geod_polygon* p, int polylinep);

  /**
   * Clear the polygon, allowing a new polygon to be started.
   *
   * @param[in,out] p a pointer to the object to be cleared.
   **********************************************************************/
  void GEOD_DLL geod_polygon_clear(struct geod_polygon* p);

  /**
   * Add a point to the polygon or polyline.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in,out] p a pointer to the geod_polygon object specifying the
   *   polygon.
   * @param[in] lat the latitude of the point (degrees).
   * @param[in] lon the longitude of the point (degrees).
   *
   * \e g and \e p must have been initialized with calls to geod_init() and
   * geod_polygon_init(), respectively.  The same \e g must be used for all the
   * points and edges in a polygon.  \e lat should be in the range
   * [&minus;90&deg;, 90&deg;].
   *
   * An example of the use of this function is given in the documentation for
   * geod_polygon_compute().
   **********************************************************************/
  void GEOD_DLL geod_polygon_addpoint(const struct geod_geodesic* g,
                                      struct geod_polygon* p,
                                      double lat, double lon);

  /**
   * Add an edge to the polygon or polyline.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in,out] p a pointer to the geod_polygon object specifying the
   *   polygon.
   * @param[in] azi azimuth at current point (degrees).
   * @param[in] s distance from current point to next point (meters).
   *
   * \e g and \e p must have been initialized with calls to geod_init() and
   * geod_polygon_init(), respectively.  The same \e g must be used for all the
   * points and edges in a polygon.  This does nothing if no points have been
   * added yet.  The \e lat and \e lon fields of \e p give the location of the
   * new vertex.
   **********************************************************************/
  void GEOD_DLL geod_polygon_addedge(const struct geod_geodesic* g,
                                     struct geod_polygon* p,
                                     double azi, double s);

  /**
   * Return the results for a polygon.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] p a pointer to the geod_polygon object specifying the polygon.
   * @param[in] reverse if non-zero then clockwise (instead of
   *   counter-clockwise) traversal counts as a positive area.
   * @param[in] sign if non-zero then return a signed result for the area if
   *   the polygon is traversed in the "wrong" direction instead of returning
   *   the area for the rest of the earth.
   * @param[out] pA pointer to the area of the polygon (meters<sup>2</sup>);
   *   only set if \e polyline is non-zero in the call to geod_polygon_init().
   * @param[out] pP pointer to the perimeter of the polygon or length of the
   *   polyline (meters).
   * @return the number of points.
   *
   * The area and perimeter are accumulated at two times the standard floating
   * point precision to guard against the loss of accuracy with many-sided
   * polygons.  Arbitrarily complex polygons are allowed.  In the case of
   * self-intersecting polygons the area is accumulated "algebraically", e.g.,
   * the areas of the 2 loops in a figure-8 polygon will partially cancel.
   * There's no need to "close" the polygon by repeating the first vertex.  Set
   * \e pA or \e pP to zero, if you do not want the corresponding quantity
   * returned.
   *
   * More points can be added to the polygon after this call.
   *
   * Example, compute the perimeter and area of the geodesic triangle with
   * vertices (0&deg;N,0&deg;E), (0&deg;N,90&deg;E), (90&deg;N,0&deg;E).
   @code{.c}
   double A, P;
   int n;
   struct geod_geodesic g;
   struct geod_polygon p;
   geod_init(&g, 6378137, 1/298.257223563);
   geod_polygon_init(&p, 0);

   geod_polygon_addpoint(&g, &p,  0,  0);
   geod_polygon_addpoint(&g, &p,  0, 90);
   geod_polygon_addpoint(&g, &p, 90,  0);
   n = geod_polygon_compute(&g, &p, 0, 1, &A, &P);
   printf("%d %.8f %.3f\n", n, P, A);
   @endcode
   **********************************************************************/
  unsigned GEOD_DLL geod_polygon_compute(const struct geod_geodesic* g,
                                         const struct geod_polygon* p,
                                         int reverse, int sign,
                                         double* pA, double* pP);

  /**
   * Return the results assuming a tentative final test point is added;
   * however, the data for the test point is not saved.  This lets you report a
   * running result for the perimeter and area as the user moves the mouse
   * cursor.  Ordinary floating point arithmetic is used to accumulate the data
   * for the test point; thus the area and perimeter returned are less accurate
   * than if geod_polygon_addpoint() and geod_polygon_compute() are used.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] p a pointer to the geod_polygon object specifying the polygon.
   * @param[in] lat the latitude of the test point (degrees).
   * @param[in] lon the longitude of the test point (degrees).
   * @param[in] reverse if non-zero then clockwise (instead of
   *   counter-clockwise) traversal counts as a positive area.
   * @param[in] sign if non-zero then return a signed result for the area if
   *   the polygon is traversed in the "wrong" direction instead of returning
   *   the area for the rest of the earth.
   * @param[out] pA pointer to the area of the polygon (meters<sup>2</sup>);
   *   only set if \e polyline is non-zero in the call to geod_polygon_init().
   * @param[out] pP pointer to the perimeter of the polygon or length of the
   *   polyline (meters).
   * @return the number of points.
   *
   * \e lat should be in the range [&minus;90&deg;, 90&deg;].
   **********************************************************************/
  unsigned GEOD_DLL geod_polygon_testpoint(const struct geod_geodesic* g,
                                           const struct geod_polygon* p,
                                           double lat, double lon,
                                           int reverse, int sign,
                                           double* pA, double* pP);

  /**
   * Return the results assuming a tentative final test point is added via an
   * azimuth and distance; however, the data for the test point is not saved.
   * This lets you report a running result for the perimeter and area as the
   * user moves the mouse cursor.  Ordinary floating point arithmetic is used
   * to accumulate the data for the test point; thus the area and perimeter
   * returned are less accurate than if geod_polygon_addedge() and
   * geod_polygon_compute() are used.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] p a pointer to the geod_polygon object specifying the polygon.
   * @param[in] azi azimuth at current point (degrees).
   * @param[in] s distance from current point to final test point (meters).
   * @param[in] reverse if non-zero then clockwise (instead of
   *   counter-clockwise) traversal counts as a positive area.
   * @param[in] sign if non-zero then return a signed result for the area if
   *   the polygon is traversed in the "wrong" direction instead of returning
   *   the area for the rest of the earth.
   * @param[out] pA pointer to the area of the polygon (meters<sup>2</sup>);
   *   only set if \e polyline is non-zero in the call to geod_polygon_init().
   * @param[out] pP pointer to the perimeter of the polygon or length of the
   *   polyline (meters).
   * @return the number of points.
   **********************************************************************/
  unsigned GEOD_DLL geod_polygon_testedge(const struct geod_geodesic* g,
                                          const struct geod_polygon* p,
                                          double azi, double s,
                                          int reverse, int sign,
                                          double* pA, double* pP);

  /**
   * A simple interface for computing the area of a geodesic polygon.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] lats an array of latitudes of the polygon vertices (degrees).
   * @param[in] lons an array of longitudes of the polygon vertices (degrees).
   * @param[in] n the number of vertices.
   * @param[out] pA pointer to the area of the polygon (meters<sup>2</sup>).
   * @param[out] pP pointer to the perimeter of the polygon (meters).
   *
   * \e lats should be in the range [&minus;90&deg;, 90&deg;].
   *
   * Arbitrarily complex polygons are allowed.  In the case self-intersecting
   * of polygons the area is accumulated "algebraically", e.g., the areas of
   * the 2 loops in a figure-8 polygon will partially cancel.  There's no need
   * to "close" the polygon by repeating the first vertex.  The area returned
   * is signed with counter-clockwise traversal being treated as positive.
   *
   * Example, compute the area of Antarctica:
   @code{.c}
   double
     lats[] = {-72.9, -71.9, -74.9, -74.3, -77.5, -77.4, -71.7, -65.9, -65.7,
               -66.6, -66.9, -69.8, -70.0, -71.0, -77.3, -77.9, -74.7},
     lons[] = {-74, -102, -102, -131, -163, 163, 172, 140, 113,
                88, 59, 25, -4, -14, -33, -46, -61};
   struct geod_geodesic g;
   double A, P;
   geod_init(&g, 6378137, 1/298.257223563);
   geod_polygonarea(&g, lats, lons, (sizeof lats) / (sizeof lats[0]), &A, &P);
   printf("%.0f %.2f\n", A, P);
   @endcode
   **********************************************************************/
  void GEOD_DLL geod_polygonarea(const struct geod_geodesic* g,
                                 double lats[], double lons[], int n,
                                 double* pA, double* pP);

  /**
   * mask values for the \e caps argument to geod_lineinit().
   **********************************************************************/
  enum geod_mask {
    GEOD_NONE         = 0U,                    /**< Calculate nothing */
    GEOD_LATITUDE     = 1U<<7  | 0U,           /**< Calculate latitude */
    GEOD_LONGITUDE    = 1U<<8  | 1U<<3,        /**< Calculate longitude */
    GEOD_AZIMUTH      = 1U<<9  | 0U,           /**< Calculate azimuth */
    GEOD_DISTANCE     = 1U<<10 | 1U<<0,        /**< Calculate distance */
    GEOD_DISTANCE_IN  = 1U<<11 | 1U<<0 | 1U<<1,/**< Allow distance as input  */
    GEOD_REDUCEDLENGTH= 1U<<12 | 1U<<0 | 1U<<2,/**< Calculate reduced length */
    GEOD_GEODESICSCALE= 1U<<13 | 1U<<0 | 1U<<2,/**< Calculate geodesic scale */
    GEOD_AREA         = 1U<<14 | 1U<<4,        /**< Calculate reduced length */
    GEOD_ALL          = 0x7F80U| 0x1FU         /**< Calculate everything */
  };

  /**
   * flag values for the \e flags argument to geod_gendirect() and
   * geod_genposition()
   **********************************************************************/
  enum geod_flags {
    GEOD_NOFLAGS      = 0U,     /**< No flags */
    GEOD_ARCMODE      = 1U<<0,  /**< Position given in terms of arc distance */
    GEOD_LONG_UNROLL  = 1U<<15  /**< Unroll the longitude */
  };

#endif
