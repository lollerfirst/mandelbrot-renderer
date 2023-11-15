#pragma once
#include <boost/multi_array.hpp>
#include <color.h>

typedef boost::multi_array<float, 2> grid_t;
typedef boost::multi_array<Color, 2> image_t;
image_t to_image(grid_t grid, Palette &palette);
void save_png(std::string filename, image_t image);