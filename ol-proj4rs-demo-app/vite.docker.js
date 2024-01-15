export default {
  build: {
    sourcemap: true,
  },
  server: {
    fs: {
        allow: ['/src']
    }
  }
}
