#include <algorithm>
#include <array>
#include <boost/multi_array.hpp>
#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <getopt.h>
#include <omp.h>
#include <string.h>
#include <utility>
#include <limits>

#include <color.h>
#include <image.h>
#include <mandelbrot.h>

#define THRESHOLD (32UL * 32UL)

struct collision_t
{
    bool active;
    bool D[4]; // E-S-W-N
};

MandelGrid::MandelGrid(
    int width, int height,
    double rmin, double rmax,
    double imin, double imax,
    long maxit, double maxdist,
    bool do_folding) : 
    grid(boost::multi_array<float, 2>{boost::extents[width][height]}),
    itmin(std::numeric_limits<float>::max())
{
    this->width = width;
    this->height = height;
    this->rmin = rmin;
    this->rmax = rmax;
    this->imin = imin;
    this->imax = imax;
    this->maxit = maxit;
    this->maxdist = maxdist;
    this->do_folding = do_folding;
}

boost::detail::multi_array::sub_array<float, 1UL> 
MandelGrid::operator[](int x)
{
    return this->grid[x];
}

typedef std::pair<std::size_t, std::size_t> bounds_t;
static char filepath[512] = "mandelbrot.png";
static char palettepath[512] = "palette.csv";

/// @brief Calculates the complex number corresponding to a
/// position on a discrete grid, using the global ranges
/// [rmin..rmax] and [imin..imax] to render.
/// @param g A discrete grid matching image resolution.
/// @param x The axis mapped to the real number range.
/// @param y The axis mapped to the imaginary number range.
/// @return The complex number corresponding to the position.
std::complex<double> scale(grid_t &g, int x, int y)
{
    int width = g.width;
    int height = g.height;
    double r = g.rmin + (g.rmax - g.rmin) * (double(x) / double(width));
    double i = g.imin + (g.imax - g.imin) * (double(height - y) / double(height));
    return std::complex<double>(r, i);
}

/// @brief Iterates a sub-range of complex numbers on a grid and
/// saves the iteration count, with fractional scaling by distance,
/// needed for escape outside of the bounds of that distance.
/// @param grid Grid of complex numbers.
/// @param x_bounds Sub-range on real axis.
/// @param y_bounds Sub-range on imaginary axis.
/// @return Information about borders of this sub-range overlapping
/// the Mandelbrot set, for pruning of regions.
collision_t mandelbrot(grid_t &grid, bounds_t x_bounds, bounds_t y_bounds)
{
    collision_t col;
    memset(&col, 0, sizeof(collision_t));

#pragma omp parallel for shared(grid, col)
    for (std::size_t y = y_bounds.first; y <= y_bounds.second; y++)
    {
        int it;
        for (std::size_t x = x_bounds.first; x <= x_bounds.second; x++)
        {
            std::complex<double> c = scale(grid, x, y);
            std::complex<double> z = c;
            for (it = 0; it < grid.maxit; it++)
            {
                z = z * z + c;
                if (std::abs(z) > grid.maxdist)
                {
                    break;
                }
                if (grid.do_folding && it == int(grid.maxit / 2))
                {
                    c += z;
                }
            }
            float di = float(it);
            // Point is in Mandelbrot set.
            if (it == grid.maxit)
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
                di -= std::log2(std::log(std::abs(z)) / std::log(grid.maxdist));
            }
            di /= float(grid.maxit);
            grid.itmin = std::min(grid.itmin, di);
            grid.itmax = std::max(grid.itmax, di);
            grid[x][y] = di;
        }

        if (it == grid.maxit)
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

/// @brief Prunes regions of the complex number space to render and
/// calls mandelbrot function to calculate iteration counts.
/// @param g Grid to track iteration counts.
/// @param xb X bounds.
/// @param yb Y bounds.
/// @param starting_quadrant
/// @return Overlap with Mandelbrot set.
collision_t recursion(grid_t &g, bounds_t xb, bounds_t yb, int starting_quadrant = 0)
{
    // Check if bounds are above base_grid limit
    if ((xb.second - xb.first + 1) * (yb.second - yb.first + 1) < THRESHOLD)
    {
        return mandelbrot(g, xb, yb);
    }
    std::array<std::pair<bounds_t, bounds_t>, 4> quadrants;

    for (int i = 0; i < 4; ++i)
    {
        auto &q = quadrants[i];

        switch (i)
        {
        case 0:
            q.first = {xb.first, (xb.first + xb.second) >> 1};
            q.second = {yb.first, (yb.first + yb.second) >> 1};
            break;

        case 1:
            q.first = {((xb.first + xb.second) >> 1) + 1, xb.second};
            q.second = {yb.first, (yb.first + yb.second) >> 1};
            break;

        case 2:
            q.first = {((xb.first + xb.second) >> 1) + 1, xb.second};
            q.second = {((yb.first + yb.second) >> 1) + 1, yb.second};
            break;

        case 3:
            q.first = {xb.first, (xb.first + xb.second) >> 1};
            q.second = {((yb.first + yb.second) >> 1) + 1, yb.second};
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
    for (int j = 0, i = starting_quadrant; j < 4; ++j, i = (i + 1) & 3)
    {
        b[i] = recursion(g, quadrants[i].first, quadrants[i].second, i);
        b_calcd[i] = true;

        if (b[i].active)
        {
            r.active = true;
            r.D[(2 + i) & 3] |= b[i].D[(2 + i) & 3]; // W + i
            r.D[(3 + i) & 3] |= b[i].D[(3 + i) & 3]; // N + i

            if (b[i].D[i] && !b_calcd[(i + 1) & 3]) // LOOK RIGHT + i
            {
                b[(i + 1) & 3] = recursion(
                    g, 
                    quadrants[(i + 1) & 3].first, 
                    quadrants[(i + 1) & 3].second, 
                    i);
                b_calcd[(i + 1) & 3] = true;
                r.D[(3 + i) & 3] |= b[(i + 1) & 3].D[(3 + i) & 3];
                r.D[i] |= b[(i + 1) & 3].D[i];

                if (b[(i + 1) & 3].D[(i + 1) & 3] && !b_calcd[(i + 2) & 3]) // LOOK DOWN + i
                {
                    b[(i + 2) & 3] = recursion(
                        g, 
                        quadrants[(i + 2) & 3].first, 
                        quadrants[(i + 2) & 3].second, 
                        i);
                    b_calcd[(i + 2) & 3] = true;
                    r.D[i] |= b[(i + 2) & 3].D[i];                     // E + i
                    r.D[(i + 1) & 3] |= b[(i + 2) & 3].D[(i + 1) & 3]; // S + i

                    if (b[(i + 2) & 3].D[(i + 2) & 3] && !b_calcd[(i + 3) & 3])
                    {
                        b[(i + 3) & 3] = recursion(
                            g, 
                            quadrants[(i + 3) & 3].first, 
                            quadrants[(i + 3) & 3].second, 
                            i);
                        b_calcd[(i + 3) & 3] = true;

                        r.D[(i + 1) & 3] |= b[(i + 3) & 3].D[(i + 1) & 3]; // S + i
                        r.D[(i + 2) & 3] |= b[(i + 3) & 3].D[(i + 2) & 3]; // W + i
                    }
                }
            }

            if (b[i].D[(i + 1) & 3] && !b_calcd[(i + 3) & 3]) // LOOK DOWN + i IF NOT YET CALCULATED
            {
                b[(i + 3) & 3] = recursion(
                    g, 
                    quadrants[(i + 3) & 3].first, 
                    quadrants[(i + 3) & 3].second, 
                    i);
                r.D[(i + 1) & 3] |= b[(i + 3) & 3].D[(i + 1) & 3]; // S + i
                r.D[(i + 2) & 3] |= b[(i + 3) & 3].D[(i + 2) & 3]; // W + i

                if (b[(i + 3) & 3].D[i] && !b_calcd[(i + 2) & 3]) // LOOK RIGHT + i IF NOT YET CALCULATED (b[(i+2)&3])
                {
                    b[(i + 2) & 3] = recursion(
                        g, 
                        quadrants[(i + 2) & 3].first, 
                        quadrants[(i + 2) & 3].second, 
                        i);
                    r.D[i] |= b[(i + 2) & 3].D[i];                     // E + i
                    r.D[(i + 1) & 3] |= b[(i + 2) & 3].D[(i + 1) & 3]; // S + i

                    if (b[(i + 2) & 3].D[(i + 3) & 3] && !b_calcd[(i + 1) & 3])
                    {
                        b[(i + 1) & 3] = recursion(
                            g, 
                            quadrants[(i + 3) & 3].first, 
                            quadrants[(i + 3) & 3].second, 
                            i);
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

/// @brief Print usage and exit.
void print_usage() 
{
    std::cerr << "Usage: mandelbrot "
              << "[--rmin=<value>] "
              << "[--rmax=<value>] "
              << "[--imin=<value>] "
              << "[--imax=<value>] "
              << "[--maxdist=<value>] "
              << "[--maxit=<value>] "
              << "[--width=<value>] "
              << "[--height=<value>] "
              << "[--output=<path>] "
              << "[--palette=<path>] "
              << "[-f]" 
              << std::endl;
    exit(EXIT_FAILURE);
}

/// @brief Parse command line arguments and set globals and references as appropriate.
/// @param argc
/// @param argv
/// @param width
/// @param height
void parse_arguments(int argc, char *argv[], grid_t &grid)
{
    struct option long_options[] = {
        {"output", required_argument, NULL, 'o'},
        {"palette", required_argument, NULL, 'p'},
        {"width", required_argument, NULL, 'w'},
        {"height", required_argument, NULL, 'h'},
        {"maxiter", required_argument, NULL, '1'},
        {"rmin", required_argument, NULL, '2'},
        {"rmax", required_argument, NULL, '3'},
        {"imin", required_argument, NULL, '4'},
        {"imax", required_argument, NULL, '5'},
        {"fold", no_argument, NULL, 'f'},
        {0, 0, 0, 0}};

    int opt;
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "o:w:h:f", long_options, &option_index)) != -1)
    {
        try
        {
            switch (opt)
            {
            case '2':
                grid.rmin = std::stod(optarg);
                break;
            case '4':
                grid.imin = std::stod(optarg);
                break;
            case '3':
                grid.rmax = std::stod(optarg);
                break;
            case '5':
                grid.imax = std::stod(optarg);
                break;
            case '1':
                grid.maxit = std::stoi(optarg);
                break;
            case 'w':
                grid.width = std::stol(optarg);
                break;
            case 'h':
                grid.height = std::stol(optarg);
                break;
            case 'o':
                strncpy(filepath, optarg, 511);
                break;
            case 'p':
                strncpy(palettepath, optarg, 511);
                break;
            case 'f':
                grid.do_folding = true;
                break;
            default:
                print_usage();
            }
        }
        catch (std::invalid_argument const &ex)
        {
            print_usage();
        }
    }
}

int main(int argc, char **argv)
{
    grid_t g;
    parse_arguments(argc, argv, g);
    std::cout << "Image size:\t\t" << g.width << "x" << g.height 
              << "\nTile size:\t\t" << THRESHOLD
              << "\nReal range:\t\t" << g.rmin << ".." << g.rmax
              << "\nImaginary range:\t" << g.imin << ".." << g.imax
              << "\nIteration bound:\t" << g.maxit
              << "\nEscape distance:\t" << g.maxdist
              << "\nFolding:\t\t" << g.do_folding
              << "\nOMP threads:\t\t" << omp_get_max_threads() 
              << "\nOutput file:\t\t" << filepath << '\n';
    bounds_t x_bounds(0, g.width - 1), y_bounds(0, g.height - 1);
    recursion(g, x_bounds, y_bounds);
    CSVPalette csv = CSVPalette();
    if (csv.read(palettepath))
    {
        image_t image = to_image(g, csv);
        save_png(filepath, image);
    }
    else
    {
        HSVPalette hsv = HSVPalette();
        image_t image = to_image(g, hsv);
        save_png(filepath, image);
    }
    return 0;
}
