#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "color.h"

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

/// @brief 
/// @param fH 
/// @param fS 
/// @param fV 
/// @param fR 
/// @param fG 
/// @param fB 
void HSVtoRGB(float fH, float fS, float fV, float &fR, float &fG, float &fB)
{
    float fC = fV * fS; // Chroma
    float fHPrime = std::fmod(fH / 60.0, 6);
    float fX = fC * (1 - std::fabs(std::fmod(fHPrime, 2) - 1));
    float fM = fV - fC;

    if (0 <= fHPrime && fHPrime < 1)
    {
        fR = fC;
        fG = fX;
        fB = 0;
    }
    else if (1 <= fHPrime && fHPrime < 2)
    {
        fR = fX;
        fG = fC;
        fB = 0;
    }
    else if (2 <= fHPrime && fHPrime < 3)
    {
        fR = 0;
        fG = fC;
        fB = fX;
    }
    else if (3 <= fHPrime && fHPrime < 4)
    {
        fR = 0;
        fG = fX;
        fB = fC;
    }
    else if (4 <= fHPrime && fHPrime < 5)
    {
        fR = fX;
        fG = 0;
        fB = fC;
    }
    else if (5 <= fHPrime && fHPrime < 6)
    {
        fR = fC;
        fG = 0;
        fB = fX;
    }
    else
    {
        fR = 0;
        fG = 0;
        fB = 0;
    }
    fR += fM;
    fG += fM;
    fB += fM;
}

/// @brief Maps iterations evenly across HSV spectrum.
/// @param iter Mandelbrot iteration count.
/// @return HSV color.
Color HSVPalette::get(float iter)
{
    float r, g, b;
    float h = iter * 360.0;
    HSVtoRGB(h, 1.0f, 1.0f, r, g, b);
    return Color
    {
        int(255.0 * r),
        int(255.0 * g),
        int(255.0 * b)
    };
}
