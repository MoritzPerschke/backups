#define _USE_MATH_DEFINES
#include <cmath>  
#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <vector>

const char* inf = "[*] ";
const char* ok  = "[+] ";
const char* err = "[-] ";

struct Vector
{        
    double x, y, z;           /* Position XYZ or color RGB */

    Vector(const Vector &b) : x(b.x), y(b.y), z(b.z) {}
    Vector(double x_=0, double y_=0, double z_=0) : x(x_), y(y_), z(z_) {}
    
    Vector operator+(const Vector &b) const 
    {
        return Vector(x + b.x, y + b.y, z + b.z);
    }

    Vector operator-(const Vector &b) const
    {
        return Vector(x - b.x, y - b.y, z - b.z);
    }

    Vector operator/(double c) const 
    {
        return Vector(x / c, y / c, z / c);
    }

    Vector operator*(double c) const 
    {
        return Vector(x * c, y * c, z * c);
    }

    friend Vector operator*(double c, const Vector &b) 
    { 
        return b * c; 
    }

    const double LengthSquared() const 
    { 
        return x * x + y * y + z * z; 
    }

    const double Length() const 
    { 
        return sqrt(LengthSquared()); 
    }

    Vector MultComponents(const Vector &b) const
    {
        return Vector(x * b.x, y * b.y, z * b.z);
    }

    const Vector Normalized() const
    {
        return Vector(x, y, z) / sqrt(x*x + y*y + z*z);
    }

    const double Dot(const Vector &b) const 
    {
        return x * b.x + y * b.y + z * b.z;
    }

    const Vector Cross(const Vector &b) const
    {
        return Vector((y * b.z) - (z * b.y), 
                      (z * b.x) - (x * b.z), 
                      (x * b.y) - (y * b.x));
    }
     
    const double Max() 
    {
        return fmax(x, fmax(x, y));
    }

    Vector& clamp() 
    {
        x = x<0 ? 0.0 : x>1.0 ? 1.0 : x;
        y = y<0 ? 0.0 : y>1.0 ? 1.0 : y;
        z = z<0 ? 0.0 : z>1.0 ? 1.0 : z;
        return *this;   
    }
};

Vector sub_abs(const Vector& a, const Vector& b){
    return Vector(
        abs(a.x - b.x),
        abs(a.y - b.y),
        abs(a.z - b.z)
    );
}

/*------------------------------------------------------------------
| Struct for standard Vector operations in 3D 
| (used for points, vectors, and colors)
------------------------------------------------------------------*/
typedef Vector Color;

struct Image 
{
    int width, height;
    Color *pixels;

    Image(int _w, int _h) : width(_w), height(_h) 
    {
        pixels = new Color[width * height];        
    }

    Image(std::string _path){

        std::ifstream image(_path);
        
        char p[2]; // see "magic number" https://de.wikipedia.org/wiki/Portable_Anymap 
        int max;
        image >> p[0];
        image >> p[1];
        image >> width;
        image >> height;
        image >> max;
        
        pixels = new Color[width * height];
        std::cout << inf << "P value: " << p << std::endl;
        std::cout << inf << "Dimensions of image: " << width << "x" << height << std::endl;
        std::cout << inf << "Max Brightness: " << max << std::endl << std::endl;

        int x, y, z;
        for (int i = 0; i < width*height ; i++){
            image >> x >> y >> z;
            pixels[i] = Color(Vector(x, y, z));
        }
        image.close();
    }

    Color getColor(int x, int y) 
    {
        int image_index = (height-y-1) * width + x; 
        return pixels[image_index];
    }

    void setColor(int x, int y, const Color &c) 
    {
        int image_index = (height-y-1) * width + x; 
        if(c.x > 0 && c.y > 0 && c.z > 0){
            pixels[image_index] = c;
        }
        else{
            pixels[image_index] = Color(Vector(0., 0., 0.));
        }
    }

    void addColor(int x, int y, const Color &c) 
    {
        int image_index = (height-y-1) * width + x; 
        pixels[image_index] = pixels[image_index] + c;
    }

    int toInteger(double x)
    { 
        /* Clamp to [0,1] */
        if (x<0.0)
            x = 0.0;        
        
        if (x>1.0)
            x = 1.0;             

        /* Apply gamma correction and convert to integer */
        return int(pow(x,1/2.2)*255+.5); 
    }

    void Save(const std::string &filename) 
    {
        /* Save image in PPM format */
        FILE *f = fopen(filename.c_str(), "wb");
        fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
        for (int i = 0; i < width * height; i++)
            fprintf(f,"%d %d %d ", toInteger(pixels[i].x), 
                                   toInteger(pixels[i].y), 
                                   toInteger(pixels[i].z));
        fclose(f);
    }

    void SaveRaw(const std::string &filename) 
    {
        /* Save image in PPM format */
        FILE *f = fopen(filename.c_str(), "wb");
        fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
        for (int i = 0; i < width * height; i++)
            fprintf(f,"%d %d %d ", (int)pixels[i].x, 
                                   (int)pixels[i].y, 
                                   (int)pixels[i].z);
        fclose(f);
    }
};

enum Rand_Gen_t { DRAND48, HALTON_17_19, STRATIFIED, BLUE_NOISE }; 
std::string to_string(Rand_Gen_t gen){
    switch(gen){
        case DRAND48:
            return "DRAND48";
        case HALTON_17_19:
            return "HALTON_17_19";
        case STRATIFIED:
            return "STRATIFIED";
        case BLUE_NOISE:
            return "BLUE_NOISE";
        default:
            break;
    }
    throw("ToString error");
}