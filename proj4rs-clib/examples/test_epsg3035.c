
#include <stdlib.h>
#include <stdio.h>

// Sample library usage

#include "proj4rs.h"


int main(void) {

    printf("Initializing\n");
    Proj4rs* src = proj4rs_proj_new("EPSG:4326");
    Proj4rs* dst = proj4rs_proj_new("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80");

    printf("Src: %s\n", proj4rs_proj_projname(src));
    printf("Dst: %s\n", proj4rs_proj_projname(dst));

    double x = 15.4213696; 
    double y = 47.0766716;

    printf("Transform\n");
    int res = proj4rs_transform(src, dst, &x , &y, NULL, 1, sizeof(double), true);
    if ( res != 1 ) {
        printf("Error:\n");
        printf("%s\n", proj4rs_last_error());
        return 1;
    } else {
        printf("x = %f\n", x); // Should be 4732659.007426
        printf("x = %f\n", y); // Should be 2677630.726961
    }

    printf("Deleting\n");
    proj4rs_proj_delete(src);
    proj4rs_proj_delete(dst);
    return 0;
}
