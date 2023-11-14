#include <png.h>
#include <string.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <boost/multi_array.hpp>
#include <omp.h>
#include <getopt.h>
#include <utility>
#include <algorithm>
#include <array>
#include <color.h>

#define THRESHOLD (32UL * 32UL)

struct collision_t
{
    bool active;
    bool D[4]; // E-S-W-N
};

typedef std::pair<std::size_t, std::size_t> bounds_t;
typedef boost::multi_array<float, 2> grid_t;

static double rmin = -1.5;
static double rmax =  0.5;
static double imin = -1.0;
static double imax =  1.0;
static float maxdist = 1000.0;
static int maxit = 100;
static char filepath[512] = "mb.png";

std::complex<double> scale(grid_t &grid, int x, int y)
{
    int width = grid.shape()[0];
    int height = grid.shape()[1];
    double r = rmin + (rmax - rmin) * (double(x) / double(width));
    double i = imin + (imax - imin) * (double(y) / double(height));
    return std::complex<double>(r, i);
}

collision_t mendelbrot(grid_t &grid, bounds_t x_bounds, bounds_t y_bounds)
{
    collision_t col;
#pragma omp parallel for shared(grid, col)
    for (std::size_t y = y_bounds.first; y <= y_bounds.second; y++)
    {
    	int it;
        for (std::size_t x = x_bounds.first; x <= x_bounds.second; x++)
        {
            std::complex<double> c = scale(grid, x, y);
            std::complex<double> z = c;
            for (it = 0; it < maxit; it++)
            {
                z = z * z + c;
                if (std::abs(z) > maxdist)
                {
                    break;
                }
            }
            float di = float(it);
                // Point is in Mandelbrot set.
                if (it == maxit)
                {
                    if (x == x_bounds.second)
                    {
                        col.D[0] = true;
                    }
                    else if (x == x_bounds.first)
                    {
                    	col.D[2] = true;
                    }
                }
            	else
            	{
                	di -= std::log2(std::log(std::abs(z)) / std::log(maxdist));
            	}
            	
            grid[x][y] = di / float(maxit);
        }
        
        if (it == maxit)
        {
        	if (y == y_bounds.first)
        	{
        		col.D[3] = true;
        	}
        	else if (y == y_bounds.second)
        	{
        		col.D[1] = true;
        	}
        }
    }
    return col;
}

collision_t recursion(grid_t &grid, bounds_t x_bounds, bounds_t y_bounds, int starting_quadrant = 0)
{
    // Check if bounds are above base_grid limit
    if ((x_bounds.second - x_bounds.first + 1) * (y_bounds.second - y_bounds.first + 1) < THRESHOLD)
    {
        return mendelbrot(grid, x_bounds, y_bounds);
    }
    std::array<std::pair<bounds_t, bounds_t>, 4> quadrants;

    for (int i = 0; i < 4; ++i)
    {
        auto &q = quadrants[i];

        switch (i)
        {
        case 0:
            q.first = {x_bounds.first, (x_bounds.first + x_bounds.second) >> 1};
            q.second = {y_bounds.first, (y_bounds.first + y_bounds.second) >> 1};
            break;

        case 1:
            q.first = {((x_bounds.first + x_bounds.second) >> 1) + 1, x_bounds.second};
            q.second = {y_bounds.first, (y_bounds.first + y_bounds.second) >> 1};
            break;

        case 2:
            q.first = {((x_bounds.first + x_bounds.second) >> 1) + 1, x_bounds.second};
            q.second = {((y_bounds.first + y_bounds.second) >> 1) + 1, y_bounds.second};
            break;

        case 3:
            q.first = {x_bounds.first, (x_bounds.first + x_bounds.second) >> 1};
            q.second = {((y_bounds.first + y_bounds.second) >> 1) + 1, y_bounds.second};
        default:
            break;
        }
    }

    collision_t b[4];
    collision_t r;
    bool b_calcd[4];

    memset(&b_calcd, 0, sizeof(b_calcd));
    memset(&b, 0, sizeof(b));
    memset(&r, 0, sizeof(r));

    // Look for a quadrant that is active
    for (int j = 0, i = starting_quadrant; j < 4; ++j, i = (i+1)&3)
    {

        b[i] = recursion(grid, quadrants[i].first, quadrants[i].second, i);
        b_calcd[i] = true;

        if (b[i].active)
        {
            r.active = true;
            r.D[(2 + i) & 3] |= b[i].D[(2 + i) & 3]; // W + i
            r.D[(3 + i) & 3] |= b[i].D[(3 + i) & 3]; // N + i

            if (b[i].D[i] && !b_calcd[(i + 1) & 3]) // LOOK RIGHT + i
            {
	      b[(i + 1) & 3] = recursion(grid, quadrants[(i + 1) & 3].first, quadrants[(i + 1) & 3].second, i);
                b_calcd[(i + 1) & 3] = true;
                r.D[(3 + i) & 3] |= b[(i + 1) & 3].D[(3 + i) & 3];
                r.D[i] |= b[(i + 1) & 3].D[i];

                if (b[(i + 1) & 3].D[(i + 1) & 3] && !b_calcd[(i + 2) & 3]) // LOOK DOWN + i
                {
		  b[(i + 2) & 3] = recursion(grid, quadrants[(i + 2) & 3].first, quadrants[(i + 2) & 3].second, i);
                    b_calcd[(i + 2) & 3] = true;
                    r.D[i] |= b[(i + 2) & 3].D[i];                     // E + i
                    r.D[(i + 1) & 3] |= b[(i + 2) & 3].D[(i + 1) & 3]; // S + i

                    if (b[(i + 2) & 3].D[(i + 2) & 3] && !b_calcd[(i + 3) & 3])
                    {
		      b[(i + 3) & 3] = recursion(grid, quadrants[(i + 3) & 3].first, quadrants[(i + 3) & 3].second, i);
                        b_calcd[(i + 3) & 3] = true;

                        r.D[(i + 1) & 3] |= b[(i + 3) & 3].D[(i + 1) & 3]; // S + i
                        r.D[(i + 2) & 3] |= b[(i + 3) & 3].D[(i + 2) & 3]; // W + i
                    }
                }
            }

            if (b[i].D[(i + 1) & 3] && !b_calcd[(i + 3) & 3]) // LOOK DOWN + i IF NOT YET CALCULATED
            {
	      b[(i + 3) & 3] = recursion(grid, quadrants[(i + 3) & 3].first, quadrants[(i + 3) & 3].second, i);
                r.D[(i + 1) & 3] |= b[(i + 3) & 3].D[(i + 1) & 3]; // S + i
                r.D[(i + 2) & 3] |= b[(i + 3) & 3].D[(i + 2) & 3]; // W + i

                if (b[(i + 3) & 3].D[i] && !b_calcd[(i + 2) & 3]) // LOOK RIGHT + i IF NOT YET CALCULATED (b[(i+2)&3])
                {
		  b[(i + 2) & 3] = recursion(grid, quadrants[(i + 2) & 3].first, quadrants[(i + 2) & 3].second, i);
                    r.D[i] |= b[(i + 2) & 3].D[i];                     // E + i
                    r.D[(i + 1) & 3] |= b[(i + 2) & 3].D[(i + 1) & 3]; // S + i

                    if (b[(i + 2) & 3].D[(i + 3) & 3] && !b_calcd[(i + 1) & 3])
                    {
		      b[(i + 1) & 3] = recursion(grid, quadrants[(i + 3) & 3].first, quadrants[(i + 3) & 3].second, i);
                        b_calcd[(i + 1) & 3] = true;

                        r.D[(i + 3) & 3] |= b[(i + 1) & 3].D[(i + 3) & 3]; // N + i
                        r.D[i] |= b[(i + 1) & 3].D[i];                     // E + i
                    }
                }
            }

            break; // exit the for loop
        }
    }
    return r;
}

using Image = boost::multi_array<Color, 2>;

void SaveImage(std::string filename, Image image)
{
    int width = image.shape()[0];
    int height = image.shape()[1];
    //int channels = 3;
    int bit_depth = 8;
    FILE *fp = fopen(filename.c_str(), "wb+");
    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_init_io(png, fp);
    png_infop info = png_create_info_struct(png);
    png_set_IHDR(
        png, info, width, height, bit_depth, PNG_COLOR_TYPE_RGB,
        PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_bytep *row_pointers = (png_bytep *)malloc(sizeof(png_bytep) * height);
    for (int y = 0; y < height; y++)
    {
        row_pointers[y] = (png_byte *)malloc(png_get_rowbytes(png, info));
        for (int x = 0; x < width; x++)
        {
            auto px = &(row_pointers[y][x * 3]);
            px[0] = image[x][y].r;
            px[1] = image[x][y].g;
            px[2] = image[x][y].b;
        }
    }
    png_write_info(png, info);
    png_write_image(png, row_pointers);
    png_write_end(png, NULL);
    fclose(fp);
    for (int y = 0; y < height; y++)
    {
        free(row_pointers[y]);
    }
    png_destroy_write_struct(&png, &info);
    free(row_pointers);
}

Image toImage(grid_t grid)
{
    std::size_t width = grid.shape()[0];
    std::size_t height = grid.shape()[1];
    Image image = Image{boost::extents[width][height]};
    for (std::size_t y = 0; y < height; y++)
    {
        for (std::size_t x = 0; x < width; x++)
        {
            auto h = grid[x][y] * 360.0;
            float r, g, b;
            HSVtoRGB(h, 1.0f, 1.0f, r, g, b);
            image[x][y] = Color{
                int(255.0 * r),
                int(255.0 * g),
                int(255.0 * b)};
        }
    }
    return image;
}

void parse_arguments(int argc, char *argv[], std::size_t &width, std::size_t &height) {
  struct option long_options[] = {
    {"output", required_argument, NULL, 'o'},
    {"width", required_argument, NULL, 'w'},
    {"height", required_argument, NULL, 'h'},
    {"maxiter", required_argument, NULL, '1'},
    {"rmin", required_argument, NULL, '2'},
    {"rmax", required_argument, NULL, '3'},
    {"imin", required_argument, NULL, '4'},
    {"imax", required_argument, NULL, '5'},
    {0, 0, 0, 0}
  };

  int opt;
  int option_index = 0;
  
  while ((opt = getopt_long(argc, argv, "o:w:h:", long_options, &option_index)) != -1) {
        switch (opt) {
        case '2':
            rmin = std::stod(optarg);
            break;
        case '4':
            imin = std::stod(optarg);
            break;
        case '3':
            rmax = std::stod(optarg);
            break;
        case '5':
            imax = std::stod(optarg);
            break;
        case '1':
            maxit = std::stoi(optarg);
            break;
        case 'w':
            width = std::stol(optarg);
            break;
        case 'h':
            height = std::stol(optarg);
            break;
	case 'o':
	    strncpy(filepath, optarg, 511);
	    break;
        default:
            std::cerr << "Usage: mandelbrot [-rmin <value>] [-rmax <value>] [-imin <value>] [-imax <value>] "
              << "[-maxdist <value>] [-maxit <value>] [-width <value>] [-height <value>]" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

int main(int argc, char** argv)
{
    std::size_t width = 2000, height = 2000;
    parse_arguments(argc, argv, width, height);
 
    grid_t grid = grid_t(boost::extents[width][height]);
    bounds_t x_bounds(0, width - 1), y_bounds(0, height - 1);

    std::cout << "Image size: "<< width <<" x " << height << "\nTile size: "<< THRESHOLD << "\nomp threads = " << omp_get_max_threads() << "\nOutput file: " << filepath << '\n';
    
    recursion(grid, x_bounds, y_bounds);
    Image image = toImage(grid);
    SaveImage(filepath, image);

    return 0;
}
