#include <boost/multi_array.hpp>
#include <color.h>

typedef boost::multi_array<float, 2> grid_t;
typedef boost::multi_array<Color, 2> image_t;
image_t toImage(grid_t grid, Palette &palette);
void savePNG(std::string filename, image_t image);