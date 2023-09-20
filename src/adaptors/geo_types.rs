use geo_types::Point;

use crate::transform::Transform;

impl Transform for Point {
    fn transform_coordinates<F>(&mut self, f: &mut F) -> crate::errors::Result<()>
    where
        F: FnMut(f64, f64, f64) -> crate::errors::Result<(f64, f64, f64)>,
    {
        let mut xy = (self.0.x, self.0.y);
        (&mut xy).transform_coordinates(f)?;
        self.set_x(xy.0);
        self.set_y(xy.1);

        Ok(())
    }
}
