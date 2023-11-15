# MacBookPro12,1 (8GB) (Tile size: 1048576)     1m45.643s
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
open mandelbrot.png
