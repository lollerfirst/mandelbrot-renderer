# Test 1
# MacBookPro12,1 (8GB) (Tile size: 1048576)     1m45s
# MacBookPro12,1 (8GB) (Tile size: 1024)        2m25s
time build/mandelbrot \
    -w 1500 \
    -h 1500 \
    --palette palettes/palette3.csv \
    --maxiter 500 \
    --rmin -0.5811427303467864211317884580 \
    --rmax -0.5811427303466769052682115420 \
    --imin 0.6532705281013165171682115420 \
    --imax 0.6532705281014260330317884580 \
    -o mandelbrot.png

# Test 2
# MacBookPro12,1 (8GB) (Tile size: 1024)        5m58s
time build/mandelbrot \
    -w 1500 \
    -h 1500 \
    --palette palettes/palette2.csv \
    --maxiter 1200 \
    --rmin -0.8748992984154357691300006 \
    --rmax -0.8748992970805744040499994 \
    --imin 0.2320902470861389459899994 \
    --imax 0.2320902484210003110700006 \
    -o mandelbrot.png

# Test 3
# MacBookPro12,1 (8GB) (Tile size: 1024)        52s
time ./mandelbrot --rmin -1.28 --rmax -1.23 --imin 0.01 --imax 0.06

