#include <png.h>
#include <image.h>

/// @brief Convert a grid with Mandelbrot iteration counts 
/// to an image of colors.
/// @param grid 
/// @param palette 
/// @return 
image_t to_image(grid_t grid, Palette &palette)
{
    std::size_t width = grid.shape()[0];
    std::size_t height = grid.shape()[1];
    image_t image = image_t{boost::extents[width][height]};
    for (std::size_t y = 0; y < height; y++)
    {
        for (std::size_t x = 0; x < width; x++)
        {
            float index = (grid[x][y] - itmin) / (itmax - itmin);
            image[x][y] = palette.get(index);
        }
    }
    return image;
}

/// @brief Save an image as PNG.
/// @param filename The path to save the file as.
/// @param image The image with color values.
void save_png(std::string filename, image_t image)
{
    int width = image.shape()[0];
    int height = image.shape()[1];
    int channels = 3;
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
            auto px = &(row_pointers[y][x * channels]);
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
