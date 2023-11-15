#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "color.h"

/// @brief Read a color palette from a CSV file. 
/// Each row starts with a float designating the fractional
/// range that color occupies: .5 is half the range, .25 a
/// quarter etc. The rest of the lines represent values in 
/// RGB ranging from 0 to 255.
/// Example:
/// 0.0, 0,255,0
/// 0.5, 255,0,0
/// 1.0, 0,0,255
/// @param file_path 
/// @return 
bool CSVPalette::read(std::string file_path)
{
    std::ifstream file(file_path);
    if (!file.is_open()) 
    {
        std::cerr << "Cannot read palette: " << file_path << std::endl;
        return false;
    }
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;

        std::vector<std::string> tokens;
        while (std::getline(iss, token, ',')) {
            tokens.push_back(token);
        }
        std::pair<float, Color> color;
        color.first = std::stod(tokens[0]);
        color.second = Color
        {                
            .r = std::stoi(tokens[1]),
            .g = std::stoi(tokens[2]),
            .b = std::stoi(tokens[3]),
        };
        this->colors.push_back(color);
    }
    file.close();
    return true;
}

void CSVPalette::print()
{
    for (const auto& color : this->colors) {
        std::cout << color.first << "\t"
                  << color.second.r << ","
                  << color.second.g << ","
                  << color.second.b << std::endl;
    }
}

/// @brief 
/// @param index 
/// @return 
Color CSVPalette::get(float index) 
{
    Color c1, c2;
    float i1, i2;
    for (auto &c : this->colors) 
    {
        i2 = c.first;
        c2 = c.second;
        if (c.first < index)
        {
            i1 = c.first;
            c1 = c.second;
        } 
        else 
        {
            break;
        }
    }

    if (i1 == i2)
    {
        return c1;
    }
    float d = (index - i1) / (i2 - i1);
    return Color
    {
        .r = (int)lerp(c1.r, c2.r, d),
        .g = (int)lerp(c1.g, c2.g, d),
        .b = (int)lerp(c1.b, c2.b, d),
    };
}

/// @brief Interpolate between two values x1 and x2.
/// @param x1 
/// @param x2 
/// @param d 
/// @return 
float lerp(float x1, float x2, float d) 
{
    return x1 + (x2 - x1) * d;
}

/// @brief Converts an RGB color to HSV.
/// @param h Hue.
/// @param s Saturation.
/// @param v Value.
/// @param r Red.
/// @param g Green.
/// @param b Blue.
void hsv_to_rgb(float h, float s, float v, float &r, float &g, float &b)
{
    float c = v * s; // Chroma
    float hprime = std::fmod(h / 60.0, 6);
    float x = c * (1 - std::fabs(std::fmod(hprime, 2) - 1));
    float m = v - c;

    if (0 <= hprime && hprime < 1)
    {
        r = c;
        g = x;
        b = 0;
    }
    else if (1 <= hprime && hprime < 2)
    {
        r = x;
        g = c;
        b = 0;
    }
    else if (2 <= hprime && hprime < 3)
    {
        r = 0;
        g = c;
        b = x;
    }
    else if (3 <= hprime && hprime < 4)
    {
        r = 0;
        g = x;
        b = c;
    }
    else if (4 <= hprime && hprime < 5)
    {
        r = x;
        g = 0;
        b = c;
    }
    else if (5 <= hprime && hprime < 6)
    {
        r = c;
        g = 0;
        b = x;
    }
    else
    {
        r = 0;
        g = 0;
        b = 0;
    }
    r += m;
    g += m;
    b += m;
}

/// @brief Maps iterations evenly across HSV spectrum.
/// @param iter Mandelbrot iteration count.
/// @return HSV color.
Color HSVPalette::get(float iter)
{
    float r, g, b;
    float h = iter * 360.0;
    hsv_to_rgb(h, 1.0f, 1.0f, r, g, b);
    return Color
    {
        int(255.0 * r),
        int(255.0 * g),
        int(255.0 * b)
    };
}
