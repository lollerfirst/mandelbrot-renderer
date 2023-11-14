#pragma once
#include <utility>

struct Color
{
    int r, g, b;
};

/* 
Example:
    CSV: 
    0.0, 0,255,0
    0.5, 255,0,0
    1.0, 0,0,255

    Results:
    Palette p = Palette("palette.csv");
    p.get(.75) -> 127 0 127
    p.get(.1) -> 51 204 0
    p.get(.9313123) -> 35 0 219
    p.get(1.9313123) -> 0 0 255
*/
class Palette 
{
    std::vector<std::pair<float, Color> > colors;
    void read(std::string file_path);
public:
    Palette(std::string path);
    Color get(float index);
    void print();
};

float lerp(float x1, float x2, float d);
void HSVtoRGB(float fH, float fS, float fV, float &fR, float &fG, float &fB);
