/******************************************************************
*
* PathTracing.cpp
*
* Description: This program demonstrates global illumination rendering
* based on the path tracing method. The intergral in the rendering
* equation is approximated via Monte-Carlo integration; explicit 
* direct lighting is included to improve quality; the rendered image 
* is saved in PPM format.
*
* The code is largely based on the software smallpt by Kevin Beason,
* released under the MIT License.
*
* Advanced Computer Graphics Proseminar WS 2019
* 
* Interactive Graphics and Simulation Group
* Department of Computer Science
* University of Innsbruck
*
*******************************************************************/

/* Standard includes */
#define _USE_MATH_DEFINES
#include <math.h>   
#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <cmath>    // For floor().
#include "lutLDBN.h"    // For Blue Noise sampling.
#include <unistd.h>     // TODO: REMOVE IT.

#ifdef _MSC_VER

double drand48()
{
    return (double)rand() / RAND_MAX;
}

#endif

using namespace std;

/*------------------------------------------------------------------
| Struct for standard Vector operations in 3D 
| (used for points, vectors, and colors)
------------------------------------------------------------------*/
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
| Scene objects are spheres; material either perfectly diffuse, 
| specular (mirror reflection) or transparent (refraction/reflection)
| (DIFFuse, SPECular, REFRactive)
------------------------------------------------------------------*/
enum Refl_t { DIFF, SPEC, REFR }; 

struct Tuple {
    double x;
    double y;

    Tuple(double x_, double y_): x(x_), y(y_) {}
};

struct Sphere 
{
    double radius;       
    Vector position; 
    Color emission, color;      
    Refl_t refl;     
    
    Sphere(double radius_, Vector position_, Vector emission_, 
           Vector color_, Refl_t refl_):
           radius(radius_), position(position_), emission(emission_), 
           color(color_), refl(refl_) {}

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

/******************************************************************
* Hard-coded scene definition: the geometry is composed of spheres
* (i.e. Cornell box walls are part of very large spheres). 
* These are defined by:
* - radius, center 
* - emitted light (light sources), surface reflectivity (~color), 
*   material
*******************************************************************/
Sphere spheres[] = 
{
    Sphere( 1e5, Vector( 1e5  +1,      40.8,      81.6),  Vector(), Vector(.75,.25,.25), DIFF), /* Left wall */
    Sphere( 1e5, Vector(-1e5 +99,      40.8,      81.6),  Vector(), Vector(.25,.25,.75), DIFF), /* Rght wall */
    Sphere( 1e5, Vector(      50,      40.8,       1e5),  Vector(), Vector(.75,.75,.75), DIFF), /* Back wall */
    Sphere( 1e5, Vector(      50,      40.8, -1e5 +170),  Vector(), Vector(),            DIFF), /* Front wall */
    Sphere( 1e5, Vector(      50,       1e5,      81.6),  Vector(), Vector(.75,.75,.75), DIFF), /* Floor */
    Sphere( 1e5, Vector(      50,-1e5 +81.6,      81.6),  Vector(), Vector(.75,.75,.75), DIFF), /* Ceiling */

    Sphere(16.5, Vector(27, 16.5, 47), Vector(), Vector(1,1,1)*.999,  SPEC), /* Mirror sphere */
    Sphere(16.5, Vector(73, 16.5, 78), Vector(), Vector(1,1,1)*.999,  REFR), /* Glas sphere */

    Sphere( 1.5, Vector(50, 81.6-16.5, 81.6), Vector(4,4,4)*100, Vector(), DIFF), /* Light */
};




enum Rand_Gen_t { DRAND48, HALTON_17_19, STRATIFIED, BLUE_NOISE }; 
Rand_Gen_t rand_generator_type = BLUE_NOISE;

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

// Global variables for halton sequence sampling:
int halton_seed = 0;

// Global variables for stratified sampling:
const int num_intervals_stratified = 100;
double stratified_sampling_arr[num_intervals_stratified];
int used_samples = num_intervals_stratified;

// Global variables for blue noise sampling:
std::vector<Point> blue_noise_samples;
const int num_blue_noise_samples = 1000000;

double halton(int index, int base) {
    /*
        The nth number of this sequence is the number n written in binary representation, inverted, and written after 
        the decimal point.
        Code inspired by: https://en.wikipedia.org/wiki/Halton_sequence and https://observablehq.com/@jrus/halton
    */

    double fraction = 1.0;
    double result = 0.0;

    while (index > 0) {
        fraction /= base;
        result += fraction * (index % base);
        index = floor(index / base);
    }
    return result;
}

void populate_arr() {
    /*
        This function populates the array stratified_sampling_arr with one random sample per interval.
    */
    for (int i = 0; i < num_intervals_stratified; i++) {
        // Generate random sample in respective interval:
        double rand_number;
        double lower_bound = ((double) i) / ((double) num_intervals_stratified);
        double upper_bound = ((double) (i + 1)) / ((double) num_intervals_stratified);
        //do {
        //    rand_number = drand48();
        //} while (!(rand_number >= lower_bound && rand_number < upper_bound));
        rand_number = lower_bound + drand48() * (upper_bound - lower_bound);

        stratified_sampling_arr[i] = rand_number;
        //stratified_sampling_arr[i] = (((double) i) / ((double) num_intervals_stratified)) + (drand48() / ((double) num_intervals_stratified));
        //cout << "stratified_sampling_arr[i] = " << stratified_sampling_arr[i] << std::endl;
        //cout << "(((double) i) / ((double) num_intervals_stratified)) = " << (((double) i) / ((double) num_intervals_stratified)) << std::endl;
        //cout << "(drand48() / ((double) num_intervals_stratified)) = " << (drand48() / ((double) num_intervals_stratified));
    }
}

double stratified() {
    if (used_samples == num_intervals_stratified) {
        // Populated array is empty, so not populated anymore.
        populate_arr();
        used_samples = 0;
    }

    // Try random interval until an interval still contains a number (so is not yet set to -1.0).
    int index;
    do {
        index = rand() % num_intervals_stratified;
    } while (stratified_sampling_arr[index] < 0.0);

    used_samples++;

    double rand_value = stratified_sampling_arr[index];

    stratified_sampling_arr[index] = -1.0;  // Invalidate this number.

    return rand_value;
}

Tuple blue_noise() {
    /*
        Using the library, num_blue_noise_samples random Blue Noise samples have been generated in the main() and
        placed inside blue_noise_samples.
        We take a random number from these Blue Noise samples.
    */
   
    //for(auto &s: samples)
    //    std::cout<< s[0]<<" "<< s[1]<<std::endl;

    // Random sampling of Blue noise needed, as close samples are very similar in blue_noise_samples.
    int position = rand() % num_blue_noise_samples;
    Tuple rand_tuple = Tuple((blue_noise_samples.at(position))[0], (blue_noise_samples.at(position))[1]);

    return rand_tuple;
}


Tuple random_generator() {
    if (rand_generator_type == DRAND48) {
        return Tuple(drand48(), drand48());
    } else if (rand_generator_type == HALTON_17_19) {
        halton_seed++;
        //cout << "(" << halton(halton_seed, 2) << ", " << halton(halton_seed, 3) << ")" << std::endl;
        //sleep(1);
        return Tuple(halton(halton_seed, 17), halton(halton_seed, 19));
    } else if (rand_generator_type == STRATIFIED) {
        return Tuple(stratified(), stratified());
    } else if (rand_generator_type == BLUE_NOISE) {
        return blue_noise();
    }
    throw("random_generator reached end of control");
}






/******************************************************************
* Check for closest intersection of a ray with the scene;
* returns true if intersection is found, as well as ray parameter
* of intersection and id of intersected object
*******************************************************************/
bool Intersect(const Ray &ray, double &t, int &id)
{
    const int n = int(sizeof(spheres) / sizeof(Sphere));
    t = 1e20;

    for (int i = 0; i < n; i ++) 
    {
        double d = spheres[i].Intersect(ray);
        if (d > 0.0  && d < t) 
        {
            t = d;
            id = i;
        }
    }
    return t < 1e20;
}


/******************************************************************
* Recursive path tracing for computing radiance via Monte-Carlo
* integration; only considers perfectly diffuse, specular or 
* transparent materials; 
* after 5 bounces Russian Roulette is used to possibly terminate rays;  
* emitted light from light source only included on first direct hit 
* (possibly via specular reflection, refraction), controlled by 
* parameter E = 0/1;  
* on diffuse surfaces light sources are explicitely sampled;
* for transparent objects, Schlick�s approximation is employed;
* for first 3 bounces obtain reflected and refracted component,
* afterwards one of the two is chosen randomly   
*******************************************************************/
Color Radiance(const Ray &ray, int depth, int E)
{
    depth++;

    int numSpheres = int(sizeof(spheres) / sizeof(Sphere));

    double t;                               
    int id = 0;  
                             
    if (!Intersect(ray, t, id))   /* No intersection with scene */
        return BackgroundColor; 

    const Sphere &obj = spheres[id];     

    Vector hitpoint = ray.org + ray.dir * t;    /* Intersection point */
    Vector normal = (hitpoint - obj.position).Normalized();  /* Normal at intersection */ 
    Vector nl = normal;

    /* Obtain flipped normal, if object hit from inside */
    if (normal.Dot(ray.dir) >= 0) 
        nl = nl*-1.0;

    Color col = obj.color; 

    /* Maximum RGB reflectivity for Russian Roulette */
    double p = col.Max();

    if (depth > 5 || !p)   /* After 5 bounces or if max reflectivity is zero */
    {
        if (random_generator().x < p)            /* Russian Roulette */
            col = col * (1/p);        /* Scale estimator to remain unbiased */
        else 
            return obj.emission * E;  /* No further bounces, only return potential emission */
    }

    if (obj.refl == DIFF)
    {                
        Tuple tuple = random_generator();
        /* Compute random reflection vector on hemisphere */
        double r1 = 2.0 * M_PI * tuple.x; 
        double r2 = tuple.y; 
        double r2s = sqrt(r2); 
        
        /* Set up local orthogonal coordinate system u,v,w on surface */
        Vector w = nl; 
        Vector u;
        
        if(fabs(w.x) > .1)
            u = Vector(0.0, 1.0, 0.0);
        else
            u = (Vector(1.0, 0.0, 0.0).Cross(w)).Normalized(); 

        Vector v = w.Cross(u);  

        /* Random reflection vector d */
        Vector d = (u * cos(r1) * r2s + 
                    v * sin(r1) * r2s + 
                    w * sqrt(1 - r2)).Normalized();  

        /* Explicit computation of direct lighting */
        Vector e;
        for (int i = 0; i < numSpheres; i ++)
        {
            const Sphere &sphere = spheres[i];
            if (sphere.emission.x <= 0 && sphere.emission.y <= 0 && sphere.emission.z <= 0) 
                continue; /* Skip objects that are not light sources */
      
            /* Randomly sample spherical light source from surface intersection */

            /* Set up local orthogonal coordinate system su,sv,sw towards light source */
            Vector sw = sphere.position - hitpoint;
            Vector su;
            
            if(fabs(sw.x) > 0.1)
                su = Vector(0.0, 1.0, 0.0);
            else
                su = Vector(1.0, 0.0, 0.0);

            su = (su.Cross(w)).Normalized();
            Vector sv = sw.Cross(su);

            /* Create random sample direction l towards spherical light source */
            Tuple tuple = random_generator();
            double cos_a_max = sqrt(1.0 - sphere.radius * sphere.radius / 
                                    (hitpoint - sphere.position).Dot(hitpoint-sphere.position));
            double eps1 = tuple.x;
            double eps2 = tuple.y;
            double cos_a = 1.0 - eps1 + eps1 * cos_a_max;
            double sin_a = sqrt(1.0 - cos_a * cos_a);
            double phi = 2.0*M_PI * eps2;
            Vector l = su * cos(phi) * sin_a + 
                       sv * sin(phi) * sin_a + 
                       sw * cos_a;
            l = l.Normalized();

            /* Shoot shadow ray, check if intersection is with light source */
            if (Intersect(Ray(hitpoint,l), t, id) && id==i)
            {  
                double omega = 2*M_PI * (1 - cos_a_max);

                /* Add diffusely reflected light from light source; note constant BRDF 1/PI */
                e = e + col.MultComponents(sphere.emission * l.Dot(nl) * omega) * M_1_PI; 
            }
        }
   
        /* Return potential light emission, direct lighting, and indirect lighting (via
           recursive call for Monte-Carlo integration */      
        return obj.emission * E + e + col.MultComponents(Radiance(Ray(hitpoint,d), depth, 0));
    } 
    else if (obj.refl == SPEC) 
    {  
        /* Return light emission mirror reflection (via recursive call using perfect
           reflection vector) */
        return obj.emission + 
            col.MultComponents(Radiance(Ray(hitpoint, ray.dir - normal * 2 * normal.Dot(ray.dir)),
                               depth, 1));
    }

    /* Otherwise object transparent, i.e. assumed dielectric glass material */
    Ray reflRay (hitpoint, ray.dir - normal * 2 * normal.Dot(ray.dir));  /* Prefect reflection */  
    bool into = normal.Dot(nl) > 0;       /* Bool for checking if ray from outside going in */
    double nc = 1;                        /* Index of refraction of air (approximately) */  
    double nt = 1.5;                      /* Index of refraction of glass (approximately) */
    double nnt;

    if(into)      /* Set ratio depending on hit from inside or outside */
        nnt = nc/nt;
    else
        nnt = nt/nc;

    double ddn = ray.dir.Dot(nl);
    double cos2t = 1 - nnt * nnt * (1 - ddn*ddn);

    /* Check for total internal reflection, if so only reflect */
    if (cos2t < 0)  
        return obj.emission + col.MultComponents( Radiance(reflRay, depth, 1));

    /* Otherwise reflection and/or refraction occurs */
    Vector tdir;

    /* Determine transmitted ray direction for refraction */
    if(into)
        tdir = (ray.dir * nnt - normal * (ddn * nnt + sqrt(cos2t))).Normalized();
    else
        tdir = (ray.dir * nnt + normal * (ddn * nnt + sqrt(cos2t))).Normalized();

    /* Determine R0 for Schlick�s approximation */
    double a = nt - nc;
    double b = nt + nc;
    double R0 = a*a / (b*b);
  
    /* Cosine of correct angle depending on outside/inside */
    double c;
    if(into)
        c = 1 + ddn;
    else
        c = 1 - tdir.Dot(normal);

    /* Compute Schlick�s approximation of Fresnel equation */ 
    double Re = R0 + (1 - R0) *c*c*c*c*c;   /* Reflectance */
    double Tr = 1 - Re;                     /* Transmittance */

    /* Probability for selecting reflectance or transmittance */
    double P = .25 + .5 * Re;
    double RP = Re / P;         /* Scaling factors for unbiased estimator */
    double TP = Tr / (1 - P);

    if (depth < 3)   /* Initially both reflection and trasmission */
        return obj.emission + col.MultComponents(Radiance(reflRay, depth, 1) * Re + 
                                                 Radiance(Ray(hitpoint, tdir), depth, 1) * Tr);
    else             /* Russian Roulette */ 
        if (random_generator().x < P)
            return obj.emission + col.MultComponents(Radiance(reflRay, depth, 1) * RP);
        else
            return obj.emission + col.MultComponents(Radiance(Ray(hitpoint,tdir), depth, 1) * TP);
}

/******************************************************************
* Main routine: Computation of path tracing image (2x2 subpixels)
* Key parameters
* - Image dimensions: width, height 
* - Number of samples per subpixel (non-uniform filtering): samples 
* Rendered result saved as PPM image file
*******************************************************************/

int main(int argc, char *argv[]) 
{
    int width = 1024;
    int height = 768;
    int samples = 1;

    if(argc == 3){
        samples = atoi(argv[1]);
        rand_generator_type = Rand_Gen_t(atoi(argv[2]));
    }else{
        std::cout << "[-] " << "Usage: <number of samples> <random generator type >" << std::endl
        << "Options are:\n0: drand48\n1: Halton_17_19\n2: Stratified\n3: Blue Noise" << std::endl;
    }

    srand (time(NULL));

    // Using the library, num_blue_noise_samples random Blue Noise samples are generated and placed inside blue_noise_samples:
    initSamplers();
    ldbnBNOT(num_blue_noise_samples, blue_noise_samples);

    /* Set camera origin and viewing direction (negative z direction) */
    Ray camera(Vector(50.0, 52.0, 295.6), Vector(0.0, -0.042612, -1.0).Normalized());

    /* Image edge vectors for pixel sampling */
    Vector cx = Vector(width * 0.5135 / height);
    Vector cy = (cx.Cross(camera.dir)).Normalized() * 0.5135;

    /*
        Render only the window region from (1 / 4) * width to (3 / 4) * width and
         (1 / 4) * height to (3 / 4) * height.
    */
    const int new_width_start = int((1.0 / 4.0) * width) - 50;       // 0 for full rendering.
    const int new_width_end = int((3.0 / 4.0) * width) + 50;         // width for full rendering.
    const int new_height_start = int((1.0 / 4.0) * height) - 180;     // 0 for full rendering.
    const int new_height_end = int((3.0 / 4.0) * height) - 80;       // height for full rendering.

    // cout << "new_width_start = " << new_width_start << std::endl;
    // cout << "new_width_end = " << new_width_end << std::endl;
    // cout << "new_height_start = " << new_height_start << std::endl;
    // cout << "new_height_end = " << new_height_end << std::endl;

    // for (int i = 1; i < 10; i++) {
    //     cout << "halton(i, 2) = " << halton(i, 2) << std::endl;
    // }
    // for (int i = 1; i < 10; i++) {
    //     cout << "halton(i, 3) = " << halton(i, 3) << std::endl;
    // }
    // for (int i = 1; i < 10; i++) {
    //     cout << "stratified(i) = " << stratified() << std::endl;
    // }

    /* Final rendering */
    Image img(new_width_end - new_width_start, new_height_end - new_height_start);

    /* Loop over image rows */
    for (int y = new_height_start; y < new_height_end; y ++) 
    { 
        cout << "\rRendering (" << samples * 4 << " spp) " << (100.0 * y / (new_height_end - 1)) << "%     ";
        srand(y * y * y);
 
        /* Loop over row pixels */
        for (int x = new_width_start; x < new_width_end; x ++)  
        {
            img.setColor(x - new_width_start, y - new_height_start, Color());
 
            /* 2x2 subsampling per pixel */
            for (int sy = 0; sy < 2; sy ++) 
            {
                for (int sx = 0; sx < 2; sx ++) 
                {
                    Color accumulated_radiance = Color();

                    /* Compute radiance at subpixel using multiple samples */
                    for (int s = 0; s < samples; s ++) 
                    {
                        Tuple tuple = random_generator();
                        const double r1 = 2.0 * tuple.x;
                        const double r2 = 2.0 * tuple.y;

                        /* Transform uniform into non-uniform filter samples */
                        double dx;
                        if (r1 < 1.0)
                            dx = sqrt(r1) - 1.0;
                        else
                            dx = 1.0 - sqrt(2.0 - r1);

                        double dy;
                        if (r2 < 1.0)
                            dy = sqrt(r2) - 1.0;
                        else
                            dy = 1.0 - sqrt(2.0 - r2);
        
                        /* Ray direction into scene from camera through sample */
                        Vector dir = cx * ((x + (sx + 0.5 + dx) / 2.0) / width - 0.5) +
                                     cy * ((y + (sy + 0.5 + dy) / 2.0) / height - 0.5) + 
                                     camera.dir;
                        
                        /* Extend camera ray to start inside box */
                        Vector start = camera.org + dir * 130.0;

                        dir = dir.Normalized();

                        /* Accumulate radiance */
                        accumulated_radiance = accumulated_radiance + 
                            Radiance( Ray(start, dir), 0, 1) / samples;
                    } 
                    
                    accumulated_radiance = accumulated_radiance.clamp() * 0.25;

                    img.addColor(x - new_width_start, y - new_height_start, accumulated_radiance);
                }
            }
        }
    }
    cout << endl;
    img.Save(string("images/" + to_string(rand_generator_type) + "/" + to_string(samples) + ".ppm"));
}