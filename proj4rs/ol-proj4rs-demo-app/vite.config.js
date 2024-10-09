export default {
  base: './',
  build: {
    sourcemap: true,
    target: 'esnext',
    rollupOptions: {
      input: {
        index: 'index.html',
        'reprojection-image': 'reprojection-image.html',
        'reprojection': 'reprojection.html',
        'sphere-mollweide': 'sphere-mollweide.html',
        'wms-image-custom-proj': 'wms-image-custom-proj.html',
        'vector-projections': 'vector-projections.html',
      },
    },
  },
  server: {
    fs: {
        allow: ['/..']
    }
  }
}
