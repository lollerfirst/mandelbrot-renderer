#pragma once
#include <boost/multi_array.hpp>

class MandelGrid 
{
    boost::multi_array<float, 2> grid;

public:
    double rmin, rmax;
    double imin, imax;
    long maxit;
    double maxdist;
    bool do_folding;
    float itmin, itmax;
    
    MandelGrid(
        std::size_t width=1000, 
        std::size_t height=1000,
        double rmin=-1.5,
        double rmax=0.5,
        double imin=-1.0, 
        double imax=1.0,
        long maxit=100, 
        double maxdist=1000.0,
        bool do_folding=false);

    boost::detail::multi_array::sub_array<float, 1>
    operator[](int x);

    void width(std::size_t width);
    void height(std::size_t height);
    std::size_t width();
    std::size_t height();
};
typedef MandelGrid grid_t;

void print_usage();

