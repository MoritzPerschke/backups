/* Standard includes */
#define _USE_MATH_DEFINES
#include <math.h>   
#include <cstdlib> 
#include <iostream>
#include <fstream>

#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
// Optional. define TINYOBJLOADER_USE_MAPBOX_EARCUT gives robust trinagulation. Requires C++11
//#define TINYOBJLOADER_USE_MAPBOX_EARCUT
#include "external/tinyobjloader/tiny_obj_loader.h"

#ifdef _MSC_VER

double drand48()
{
    return (double)rand() / RAND_MAX;
}

#endif

using namespace std;

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

/*------------------------------------------------------------------
| Struct for standard Vector operations in 3D 
| (used for points, vectors, and colors)
------------------------------------------------------------------*/
typedef Vector Color;
const Color BackgroundColor(0.0, 0.0, 0.0);

/*------------------------------------------------------------------
| Structure for rays (e.g. viewing ray, path tracing)
------------------------------------------------------------------*/
struct Ray 
{
    Vector org, dir;    /* Origin and direction */
    Ray(const Vector org_, const Vector &dir_) : org(org_), dir(dir_) {}
};

/*------------------------------------------------------------------
| Struct holds pixels/colors of rendered image
------------------------------------------------------------------*/
struct Image 
{
    int width, height;
    Color *pixels;

    Image(int _w, int _h) : width(_w), height(_h) 
    {
        pixels = new Color[width * height];        
    }

    Color getColor(int x, int y) 
    {
        int image_index = (height-y-1) * width + x; 
        return pixels[image_index];
    }

    void setColor(int x, int y, const Color &c) 
    {
        int image_index = (height-y-1) * width + x; 
        pixels[image_index] = c;
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

    void Save(const string &filename) 
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
};

/*------------------------------------------------------------------
| Scene objects are spheres and triangles; material either perfectly 
| diffuse, specular (mirror reflection), transparent (refraction/
| reflection), glossy or translucent.
| (DIFFuse, SPECular, REFRactive, GLOSsy, TRANslucent)
------------------------------------------------------------------*/
enum Refl_t { DIFF, SPEC, REFR, GLOS, TRAN }; 

struct Primitive 
{
    Color emission, color;      
    Refl_t refl;     
    
    Primitive(Vector emission_, Vector color_, Refl_t refl_):
           emission(emission_), color(color_), refl(refl_) {}

    virtual double Intersect(const Ray &ray) const = 0;
    virtual void Print() = 0;
    virtual Vector getNormal(Vector hitpoint) const = 0;
    virtual Vector getPosition() const = 0;
    virtual double getRadius() const = 0;
};

/*------------------------------------------------------------------
| Basic geometric element of scene description;
| Triangles are subdivided into smaller patches for radiosity
| computation (subdivision equal for all rectangles)
------------------------------------------------------------------*/
struct Triangle : public Primitive {
    Vector p0;
    Vector edge_a, edge_b;
    Vector normal;

    double a_len, b_len; /* Edge lengths */

    Triangle(const Vector p0_, const Vector a_, const Vector b_,
        const Color emission_, const Color color_, const Refl_t refl_)
        : Primitive(emission_, color_, refl_), p0(p0_), edge_a(a_), edge_b(b_) {
        normal = edge_a.Cross(edge_b);
        normal = normal.Normalized();
        a_len = edge_a.Length();
        b_len = edge_b.Length();
    }

    void Print() {
        cout << "Triangle - p0 = " << p0.x << ", " << p0.y << ", " << p0.z << ") and "
            << edge_a.x << ", " << edge_a.y << ", " << edge_a.z << ") and "
            << edge_b.x << ", " << edge_b.y << ", " << edge_b.z << ") and "
            << "emission = " << emission.x << ", " << emission.y << ", " << emission.z << ")"
            << "color = " << color.x << ", " << color.y << ", " << color.z << ")"
            << "normal = " << normal.x << ", " << normal.y << ", " << normal.z << ")"
            << "refl = " << refl << "\n";
    }

    Vector getNormal(Vector hitpoint) const {
        // For the triangle, the hitpoint does not matter. The normal is the same everywhere on the triangle.
        return normal;
    }

    Vector getPosition() const {
        /*
            This function returns the position of the triangle, for which we simply take its centroid (= Schwerpunkt).
            However, this function is more a dummy function, as it is never called in the actual program.
        */
        Vector p1 = p0 + edge_a;
        Vector p2 = p0 + edge_b;

        return (1 / 3) * (p0 + p1 + p2);
    }

    double getRadius() const {
        /*
            This function returns the circumradius (= Umkreis-Radius) of the triangle.
            However, this function is more a dummy function, as it is never called in the actual program.
        */
        Vector edge_c = edge_a - edge_b;
        double a = edge_a.Length();
        double b = edge_b.Length();
        double c = edge_c.Length();

        double s = (a + b + c) / 2;     // This parameter is called the semiperimeter.
        double area = sqrt(s * (s - a) * (s - b) * (s - c));       // Area of the triangle by Heron's Formula.
        double circumradius = (a * b * c) / (4 * area);

        return circumradius;
    }

    /* Triangle-ray intersection */
    double Intersect(const Ray &ray) const {
        /* Check for plane-ray intersection first */
        const double t = (p0 - ray.org).Dot(normal) / ray.dir.Dot(normal);
        if (t <= 0.00001)
            return 0.0;

        /* Determine if intersection is within triangle */
        Vector p = ray.org + ray.dir * t;
        // Vector d = p - p0;  // Vector from point p0 of triangle to point p in
        // question.

        // Triangle inside-outside testing:
        /*
            In particular, p lies inside of triangle if it lies on the inside of all
            of the lines determined by edges of the triangle.

            Idea from:
            https://courses.cs.washington.edu/courses/csep557/10au/lectures/triangle_intersection.pdf
        */

        // Triangle consists of 3 points:
        Vector point_A = p0;
        Vector point_B = point_A + edge_a;
        Vector point_C = point_A + edge_b;

        // p needs to lie inside all 3 edges of the triangle to lie on the inside:
        const bool inside_edge_a =
            ((point_B - point_A).Cross(p - point_A)).Dot(normal) >= 0;
        const bool inside_edge_b =
            ((point_C - point_B).Cross(p - point_B)).Dot(normal) >= 0;
        const bool inside_edge_c =
            ((point_A - point_C).Cross(p - point_C)).Dot(normal) >= 0;

        if (inside_edge_a && inside_edge_b && inside_edge_c) {
        // p lies on inside of the triangle.
            return t;
        } else {
        // p lies on outside of the triangle.
            return 0.0;
        }
    }
};

struct Sphere : public Primitive {
    double radius;       
    Vector position; 
    
    Sphere(double radius_, Vector position_, Vector emission_, 
           Vector color_, Refl_t refl_):
           Primitive(emission_, color_, refl_), radius(radius_), position(position_) {}

    void Print() {
        cout << "Sphere - radius = " << radius << " and pos = " << position.x << ", " << position.y << ", " << position.z << ") and "
        << "emission = " << emission.x << ", " << emission.y << ", " << emission.z << ")"
        << "color = " << color.x << ", " << color.y << ", " << color.z << ")"
        << "refl = " << refl << "\n";
    }

    Vector getNormal(Vector hitpoint) const {
        // For the sphere, the normal depends on the hitpoint. The normal is NOT the same everywhere on the sphere.
        Vector normal = (hitpoint - position).Normalized();  /* Normal at intersection */ 
        return normal;
    }

    Vector getPosition() const {
        return position;
    }

    double getRadius() const {
        return radius;
    }

    double Intersect(const Ray &ray) const 
    { 
        /* Check for ray-sphere intersection by solving for t:
            t^2*d.d + 2*t*(o-p).d + (o-p).(o-p) - R^2 = 0 */
        Vector op = position - ray.org; 
        double eps = 1e-4;
        double b = op.Dot(ray.dir);
        double radicant = b*b - op.Dot(op) + radius*radius;
        if (radicant < 0.0) 
            return 0.0;      /* No intersection */
        else   
            radicant = sqrt(radicant);
        
        double t;
        t = b - radicant;    /* Check smaller root first */
        if(t > eps)
            return t;
        
        t = b + radicant;
        if(t > eps)          /* Check second root */
            return t;
        
        return 0.0;          /* No intersection in ray direction */  
    }
};