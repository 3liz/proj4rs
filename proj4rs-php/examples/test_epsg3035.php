<?php

$src = new Projection("WGS84");
$dst = new Projection("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80");

var_dump($src->projName);
var_dump($dst->projName);

$point = new Point(15.4213696, 47.0766716, 0.);

var_dump($point);

print "==> Transform <==\n";
transform_point($src, $dst, $point, true);

var_dump($point);
