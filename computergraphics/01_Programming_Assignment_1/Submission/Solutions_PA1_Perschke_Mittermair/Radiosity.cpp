/*==================================================================*/
//           Thomas Mittermair, Moritz Perschke                     //
/*==================================================================*/

/******************************************************************
 *
 * Radiosity.cpp
 *
 * Description: This file demonstrates global illumination rendering
 * based on the radiosity method. The geometry is divided into patches
 * for which the form factors are determined employing Monte Carlo
 * integration. Radiosity values for the patches are computed with
 * an iterative solver. The final image (i.e. radiance) is obtained
 * via tracing rays into the scene. Two output files are saved -
 * one with constant shading of patches and one with bicubic color
 * interpolation.
 *
 * The code is extended from software by user Hole and Kevin Beason,
 * released under the MIT License.
 *
 * http://kagamin.net/hole/license.txt
 * http://kagamin.net/hole/smallpt-license.txt
 *
 * Advanced Computer Graphics Practical Session SS 2022
 *
 * Interactive Graphics and Simulation Group
 * Institute of Computer Science
 * University of Innsbruck
 * adjusted to enable for MSVC compile
 *******************************************************************/

/* Standard includes */
#define _USE_MATH_DEFINES
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <random>
#include <vector>

using namespace std;

const double Over_M_PI = 1.0 / M_PI;

static double* form_factor;
static int patch_num = 0;

#ifdef _MSC_VER

double drand48() { return (double)rand() / RAND_MAX; }

#endif

/*------------------------------------------------------------------
| Struct for standard vector operations in 3D
| (used for points, vectors, and colors)
------------------------------------------------------------------*/
struct Vector {
  double x, y, z; /* Position XYZ or color RGB */

  Vector(const Vector &b) : x(b.x), y(b.y), z(b.z) {}
  Vector(double x_ = 0, double y_ = 0, double z_ = 0) : x(x_), y(y_), z(z_) {}

  Vector operator+(const Vector &b) const {
    return Vector(x + b.x, y + b.y, z + b.z);
  }

  Vector operator-(const Vector &b) const {
    return Vector(x - b.x, y - b.y, z - b.z);
  }

  // Unary minus operator:
  Vector operator-() const { return Vector(-x, -y, -z); }

  Vector operator/(double c) const { return Vector(x / c, y / c, z / c); }

  Vector operator*(double c) const { return Vector(x * c, y * c, z * c); }

  friend Vector operator*(double c, const Vector &b) { return b * c; }

  Vector MultComponents(const Vector &b) const {
    return Vector(x * b.x, y * b.y, z * b.z);
  }

  const double LengthSquared() const { return x * x + y * y + z * z; }

  const double Length() const { return sqrt(LengthSquared()); }

  const Vector Normalized() const {
    return Vector(x, y, z) / sqrt(x * x + y * y + z * z);
  }

  const double Dot(const Vector &b) const {
    return x * b.x + y * b.y + z * b.z;
  }

  const Vector Cross(const Vector &b) const {
    return Vector((y * b.z) - (z * b.y), (z * b.x) - (x * b.z),
                  (x * b.y) - (y * b.x));
  }
};

typedef Vector Color;
const Color BackgroundColor(0.0, 0.0, 0.0);

/*------------------------------------------------------------------
| Structure for rays (e.g. viewing ray, ray tracing)
------------------------------------------------------------------*/
struct Ray {
  Vector org, dir; /* Origin and direction */
  Ray(const Vector org_, const Vector &dir_) : org(org_), dir(dir_) {}
};

/*------------------------------------------------------------------
| Struct holds pixels/colors of rendered image
------------------------------------------------------------------*/
struct Image {
  int width, height;
  Color *pixels;

  Image(int _w, int _h) : width(_w), height(_h) {
    pixels = new Color[width * height];
  }

  Color getColor(int x, int y) {
    int image_index = (height - y - 1) * width + x;
    return pixels[image_index];
  }

  void setColor(int x, int y, const Color &c) {
    int image_index = (height - y - 1) * width + x;
    pixels[image_index] = c;
  }

  void addColor(int x, int y, const Color &c) {
    int image_index = (height - y - 1) * width + x;
    pixels[image_index] = pixels[image_index] + c;
  }

  int toInteger(double x) {
    /* Clamp to [0,1] */
    if (x < 0.0)
      x = 0.0;

    if (x > 1.0)
      x = 1.0;

    /* Apply gamma correction and convert to integer */
    return int(pow(x, 1 / 2.2) * 255 + .5);
  }

  void Save(const string &filename) {
    /* Save image in PPM format */
    FILE *f = fopen(filename.c_str(), "wb");
    fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
    for (int i = 0; i < width * height; i++)
      fprintf(f, "%d %d %d ", toInteger(pixels[i].x), toInteger(pixels[i].y),
              toInteger(pixels[i].z));
    fclose(f);
  }
};

/*------------------------------------------------------------------
| Basic geometric element of scene description;
| Triangles are subdivided into smaller patches for radiosity
| computation (subdivision equal for all rectangles)
------------------------------------------------------------------*/
struct Triangle {
  Vector p0;
  Vector edge_a, edge_b;
  Color emission, color;
  Vector normal;

  vector<Color> patch;
  int a_num, b_num;    /* Number of patches/subdivision of edges */
  double a_len, b_len; /* Edge lengths */

  Triangle(const Vector p0_, const Vector &a_, const Vector &b_,
           const Color &emission_, const Color &color_)
      : p0(p0_), edge_a(a_), edge_b(b_), emission(emission_), color(color_) {
    normal = edge_a.Cross(edge_b);
    normal = normal.Normalized();
    a_len = edge_a.Length();
    b_len = edge_b.Length();
  }

  Color sample_patch(int ia, int ib) const {
    if (ia < 0)
      ia = 0;
    if (ia >= a_num)
      ia = a_num - 1;
    if (ib < 0)
      ib = 0;
    if (ib >= b_num)
      ib = b_num - 1;
    return patch[ia * b_num + ib];
  }

  void init_patchs(const int a_num_, const int b_num_) {
    a_num = a_num_;
    b_num = b_num_;
    patch.clear();
    patch.resize(a_num * b_num);
  }

  /* Triangle-ray intersection */
  const double intersect(const Ray &ray) {
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

    // const double ddota = d.Dot(edge_a);
    // if (ddota < 0.0 || ddota > edge_a.LengthSquared())
    //     return 0.0;

    // const double ddotb = d.Dot(edge_b);
    // if (ddotb < 0.0 || ddotb > edge_b.LengthSquared())
    //     return 0.0;

    // return t;
  }
};

/******************************************************************
 * Hard-coded scene definition: the geometry is composed of triangles.
 * These are defined by:
 * - vector to corner(origin), edge a, edge b
 * - emitted light energy (light sources), surface reflectivity (~color)
 *******************************************************************/
Triangle tris_of_scene[] = {
    /*
        Cornell Box walls
        Each wall consists of two triangles.
    */
    Triangle(Vector(0.0, 0.0, 0.0), Vector(100.0, 0.0, 0.0),
             Vector(0.0, 80.0, 0.0), Vector(),
             Color(0.75, 0.75, 0.75)), /* Back */
    Triangle(Vector(0.0, 0.0, 0.0) + Vector(100.0, 0.0, 0.0) +
                 Vector(0.0, 80.0, 0.0),
             -Vector(100.0, 0.0, 0.0), -Vector(0.0, 80.0, 0.0), Vector(),
             Color(0.75, 0.75, 0.75)), /* Back */

    Triangle(Vector(0.0, 0.0, 170.0), Vector(100.0, 0.0, 0.0),
             Vector(0.0, 0.0, -170.0), Vector(),
             Color(0.75, 0.75, 0.75)), /* Bottom */
    Triangle(Vector(0.0, 0.0, 170.0) + Vector(100.0, 0.0, 0.0) +
                 Vector(0.0, 0.0, -170.0),
             -Vector(100.0, 0.0, 0.0), -Vector(0.0, 0.0, -170.0), Vector(),
             Color(0.75, 0.75, 0.75)), /* Bottom */

    Triangle(Vector(0.0, 80.0, 0.0), Vector(100.0, 0.0, 0.0),
             Vector(0.0, 0.0, 170.0), Vector(),
             Color(0.75, 0.75, 0.75)), /* Top */
    Triangle(Vector(0.0, 80.0, 0.0) + Vector(100.0, 0.0, 0.0) +
                 Vector(0.0, 0.0, 170.0),
             -Vector(100.0, 0.0, 0.0), -Vector(0.0, 0.0, 170.0), Vector(),
             Color(0.75, 0.75, 0.75)), /* Top */

    Triangle(Vector(0.0, 0.0, 170.0), Vector(0.0, 0.0, -170.0),
             Vector(0.0, 80.0, 0.0), Vector(),
             Color(0.75, 0.25, 0.25)), /* Left */
    Triangle(Vector(0.0, 0.0, 170.0) + Vector(0.0, 0.0, -170.0) +
                 Vector(0.0, 80.0, 0.0),
             -Vector(0.0, 0.0, -170.0), -Vector(0.0, 80.0, 0.0), Vector(),
             Color(0.75, 0.25, 0.25)), /* Left */

    Triangle(Vector(100.0, 0.0, 0.0), Vector(0.0, 0.0, 170.0),
             Vector(0.0, 80.0, 0.0), Vector(),
             Color(0.25, 0.25, 0.75)), /* Right */
    Triangle(Vector(100.0, 0.0, 0.0) + Vector(0.0, 0.0, 170.0) +
                 Vector(0.0, 80.0, 0.0),
             -Vector(0.0, 0.0, 170.0), -Vector(0.0, 80.0, 0.0), Vector(),
             Color(0.25, 0.25, 0.75)), /* Right */

    Triangle(Vector(100.0, 0.0, 170.0), Vector(-100.0, 0.0, 0.0),
             Vector(0.0, -80.0, 0.0), Vector(),
             Color(0, 1, 0)), /* Front (not visible) */
    Triangle(Vector(100.0, 0.0, 170.0) + Vector(-100.0, 0.0, 0.0) +
                 Vector(0.0, -80.0, 0.0),
             -Vector(-100.0, 0.0, 0.0), -Vector(0.0, -80.0, 0.0), Vector(),
             Color(0, 1, 0)), /* Front (not visible) */

    /* Area light source on top */
    Triangle(Vector(40.0, 79.99, 65.0), Vector(20.0, 0.0, 0.0),
             Vector(0.0, 0.0, 20.0), Vector(12, 12, 12),
             Color(0.75, 0.75, 0.75)),
    Triangle(Vector(40.0, 79.99, 65.0) + Vector(20.0, 0.0, 0.0) +
                 Vector(0.0, 0.0, 20.0),
             -Vector(20.0, 0.0, 0.0), -Vector(0.0, 0.0, 20.0),
             Vector(12, 12, 12), Color(0.75, 0.75, 0.75)),

    /* Cuboid in room */
    Triangle(Vector(30.0, 0.0, 100.0), Vector(0.0, 0.0, -20.0),
             Vector(0.0, 40.0, 0.0), Vector(),
             Color(0.75, 0.75, 0.75)), /* Right */
    Triangle(Vector(30.0, 0.0, 100.0) + Vector(0.0, 0.0, -20.0) +
                 Vector(0.0, 40.0, 0.0),
             -Vector(0.0, 0.0, -20.0), -Vector(0.0, 40.0, 0.0), Vector(),
             Color(0.75, 0.75, 0.75)), /* Right */

    Triangle(Vector(10.0, 0.0, 80.0), Vector(0.0, 0.0, 20.0),
             Vector(0.0, 40.0, 0.0), Vector(),
             Color(0.75, 0.75, 0.75)), /* Left */
    Triangle(Vector(10.0, 0.0, 80.0) + Vector(0.0, 0.0, 20.0) +
                 Vector(0.0, 40.0, 0.0),
             -Vector(0.0, 0.0, 20.0), -Vector(0.0, 40.0, 0.0), Vector(),
             Color(0.75, 0.75, 0.75)), /* Left */

    Triangle(Vector(10.0, 0.0, 100.0), Vector(20.0, 0.0, 0.0),
             Vector(0.0, 40.0, 0.0), Vector(),
             Color(0.75, 0.75, 0.75)), /* Front */
    Triangle(Vector(10.0, 0.0, 100.0) + Vector(20.0, 0.0, 0.0) +
                 Vector(0.0, 40.0, 0.0),
             -Vector(20.0, 0.0, 0.0), -Vector(0.0, 40.0, 0.0), Vector(),
             Color(0.75, 0.75, 0.75)), /* Front */

    Triangle(Vector(30.0, 0.0, 80.0), Vector(-20.0, 0.0, 0.0),
             Vector(0.0, -40.0, 0.0), Vector(),
             Color(0.75, 0.75, 0.75)), /* Back */
    Triangle(Vector(30.0, 0.0, 80.0) + Vector(-20.0, 0.0, 0.0) +
                 Vector(0.0, -40.0, 0.0),
             -Vector(-20.0, 0.0, 0.0), -Vector(0.0, -40.0, 0.0), Vector(),
             Color(0.75, 0.75, 0.75)), /* Back */

    Triangle(Vector(10.0, 40.0, 100.0), Vector(20.0, 0.0, 0.0),
             Vector(0.0, 0.0, -20.0), Vector(),
             Color(0.75, 0.75, 0.75)), /* Top */
    Triangle(Vector(10.0, 40.0, 100.0) + Vector(20.0, 0.0, 0.0) +
                 Vector(0.0, 0.0, -20.0),
             -Vector(20.0, 0.0, 0.0), -Vector(0.0, 0.0, -20.0), Vector(),
             Color(0.75, 0.75, 0.75)), /* Top */
};

/*
    The tris-array will contain all triangles of the scene, subdivided to
   patches (which are themselves triangles, so instances of the Triangle-struct)
   up to a certain depth.
*/
Triangle *tris = tris_of_scene;
int num_tris = int(sizeof(tris_of_scene) /
                   sizeof(Triangle)); // Number of triangles in the tris-array.

Triangle *Split_Triangles(Triangle *tris) {
  /*
      This function splits every triangle in tris into its four sub-triangles.
      All subtriangles are then returned.
  */
  std::vector<Triangle> *splitted_triangles =
      new vector<Triangle>; // Allocate vector on the heap, to use it after
                            // the function.

  for (int i = 0; i < num_tris; i++) { // Iterate over all triangles.
    // Create 4 new small triangles out of our big old triangle.
    Triangle current_triangle = tris[i];

    Triangle t1 = Triangle(current_triangle.p0, current_triangle.edge_a / 2,
                           current_triangle.edge_b / 2,
                           current_triangle.emission, current_triangle.color);
    Triangle t2 =
        Triangle(current_triangle.p0 + current_triangle.edge_a / 2,
                 current_triangle.edge_a / 2, current_triangle.edge_b / 2,
                 current_triangle.emission, current_triangle.color);
    Triangle t3 =
        Triangle(current_triangle.p0 + current_triangle.edge_b / 2,
                 current_triangle.edge_a / 2, current_triangle.edge_b / 2,
                 current_triangle.emission, current_triangle.color);
    Triangle t4 =
        Triangle(current_triangle.p0 + (current_triangle.edge_a / 2 +
                                        current_triangle.edge_b / 2),
                 -current_triangle.edge_a / 2, -current_triangle.edge_b / 2,
                 current_triangle.emission, current_triangle.color);

    // Add the 4 sub-triangles to our array:
    splitted_triangles->push_back(t1);
    splitted_triangles->push_back(t2);
    splitted_triangles->push_back(t3);
    splitted_triangles->push_back(t4);
  }

  num_tris = 4 * num_tris;

  return &(*splitted_triangles)[0]; // Convert vector<> to array and return it.
}

/******************************************************************
 * Check for closest intersection of a ray with the scene;
 * Returns true if intersection is found, as well as ray parameter
 * of intersection and id of intersected object
 *******************************************************************/
bool Intersect_Scene(const Ray &ray, double *t, int *id, Vector *normal) {
  const int n = num_tris;
  *t = 1e20;
  *id = -1;

  for (int i = 0; i < n; i++) {
    double d = tris[i].intersect(ray);
    if (d > 0.0 && d < *t) {
      *t = d;
      *id = i;
      *normal = tris[i].normal;
    }
  }
  return *t < 1e20;
}

/******************************************************************
 * Determine all form factors for all pairs of patches (of all
 * triangles);
 * Evaluation of integrals in form factor equation is done via
 * Monte Carlo integration; samples are uniformly distributed and
 * equally weighted;
 * Computation accelerated by exploiting symmetries of form factor
 * estimation;
 *******************************************************************/
void Calculate_Form_Factors(const int a_div_num, const int b_div_num,
                            const int mc_sample) {
  /* Total number of patches in scene */
  const int n = num_tris;
  for (int i = 0; i < n; i++) {
    tris[i].init_patchs(a_div_num, b_div_num);
    patch_num += tris[i].a_num * tris[i].b_num;
  }

  std::cout << "Number of triangles: " << n << endl;
  cout << "Number of patches: " << patch_num << endl;
  int form_factor_num = patch_num * patch_num;
  cout << "Number of form factors: " << form_factor_num << endl;

  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dis(0.0, 1.0);

  /* 1D-array to hold form factor pairs */
  form_factor = new double[form_factor_num];
  memset(form_factor, 0.0, sizeof(double) * form_factor_num);

  /* 1D-array with patch areas */
  double *patch_area = new double[patch_num];
  memset(patch_area, 0.0, sizeof(double) * patch_num);

  /* Precompute patch areas, assuming same size for each triangle */
  for (int i = 0; i < n; i++) {
    int patch_i = 0;
    for (int k = 0; k < i; k++)
      patch_i += tris[k].a_num * tris[k].b_num;

    for (int ia = 0; ia < tris[i].a_num; ia++) {
      for (int ib = 0; ib < tris[i].b_num; ib++) {
        patch_area[patch_i + ia * tris[i].b_num + ib] =
            ((tris[i].edge_a / tris[i].a_num)
                 .Cross((tris[i].edge_b / tris[i].b_num)))
                .Length() /
            2;
      }
    }
  }

  /* Offsets for indexing of patches in 1D-array */
  int *offset = new int[n];

  for (int i = 0; i < n; i++) {
    offset[i] = 0;
    for (int k = 0; k < i; k++)
      offset[i] += tris[k].a_num * tris[k].b_num;
  }

  /* Loop over all triangles in scene */
  // #pragma omp parallel for
  for (int i = 0; i < n; i++) {
    int patch_i = offset[i];

    cout << i << " ";

    /* Loop over all patches in triangle i */
    for (int ia = 0; ia < tris[i].a_num; ia++) {
      // cout << "*" << flush;
      for (int ib = 0; ib < tris[i].b_num; ib++) {
        const Vector normal_i = tris[i].normal;

        int patch_j = 0;

        /* Loop over all triangles in scene for triangle i */
        for (int j = 0; j < n; j++) {
          const Vector normal_j = tris[j].normal;

          /* Loop over all patches in triangle j */
          for (int ja = 0; ja < tris[j].a_num; ja++) {
            for (int jb = 0; jb < tris[j].b_num; jb++) {
              /* Do not compute form factors for patches on same triangle;
                 also exploit symmetry to reduce computation;
                 intemediate values; will be divided by patch area below */
              if (i < j) {
                double F = 0;

                /* Monte Carlo integration of form factor double integral */
                const int Ni = mc_sample, Nj = mc_sample;

                /* Uniform PDF for Monte Carlo (1/Ai)x(1/Aj) */
                const double pdf =
                    (1.0 / patch_area[offset[i] + ia * tris[i].b_num + ib]) *
                    (1.0 / patch_area[offset[j] + ja * tris[j].b_num + jb]);

                Vector ip_zero, ip_one, ip_two, jp_zero, jp_one, jp_two;

                ip_zero = tris[i].p0;
                ip_one = tris[i].p0 + tris[i].edge_a;
                ip_two = tris[i].p0 + tris[i].edge_b;

                jp_zero = tris[j].p0;
                jp_one = tris[j].p0 + tris[j].edge_a;
                jp_two = tris[j].p0 + tris[j].edge_b;

                /* Determine rays of NixNi uniform samples of patch
                   on i to NjxNj uniform samples of patch on j */
                for (int ia = 0; ia < Ni; ia++) {
                  double xi_zero_i = dis(gen);
                  double xi_one_i = dis(gen);

                  double iu0, iu1, iu2;
                  iu0 = 1. - sqrt(xi_zero_i);
                  iu1 = xi_one_i * sqrt(xi_zero_i);
                  iu2 = 1 - iu0 - iu1;

                  Vector xi = ip_zero * iu0 + ip_one * iu1 + ip_two * iu2;

                  for (int jb = 0; jb < Nj; jb++) {
                    /* Determine sample points xi, xj on both patches */
                    double xi_zero_j = dis(gen);
                    double xi_one_j = dis(gen);

                    double ju0, ju1, ju2;
                    ju0 = 1. - sqrt(xi_zero_j);
                    ju1 = xi_one_j * sqrt(xi_zero_j);
                    ju2 = 1. - ju0 - ju1;

                    /* Calculate sampling points */
                    Vector xj = jp_zero * ju0 + jp_one * ju1 + jp_two * ju2;

                    /* Check for visibility between sample points */
                    const Vector ij = (xj - xi).Normalized();

                    double t;
                    int id;
                    Vector normal;
                    if (Intersect_Scene(Ray(xi, ij), &t, &id, &normal) &&
                        id != j) {
                      continue; /* If intersection with other triangle */
                    }

                    /* Cosines of angles beteen normals and ray inbetween */
                    const double d0 = normal_i.Dot(ij);
                    const double d1 = normal_j.Dot(-1.0 * ij);

                    /* Continue if patches facing each other */
                    if (d0 > 0.0 && d1 > 0.0) {
                      /* Sample form factor */
                      const double K =
                          d0 * d1 / (M_PI * (xj - xi).LengthSquared());

                      /* Add weighted sample to estimate */
                      F += K / pdf;
                    }
                  }
                }

                /* Divide by number of samples */
                F /= mc_sample * mc_sample;

                form_factor[patch_i * patch_num + patch_j] = F;
              }
              patch_j++;
            }
          }
        }
        patch_i++;
      }
    }

    cout << endl;
  }

  /* Copy upper to lower triangular values */
  for (int i = 0; i < patch_num - 1; i++) {
    for (int j = i + 1; j < patch_num; j++) {
      form_factor[j * patch_num + i] = form_factor[i * patch_num + j];
    }
  }

  /* Divide by area to get final form factors */
  for (int i = 0; i < patch_num; i++) {
    for (int j = 0; j < patch_num; j++) {
      form_factor[i * patch_num + j] /= patch_area[i];

      /* Clamp to [0,1] */
      if (form_factor[i * patch_num + j] > 1.0)
        form_factor[i * patch_num + j] = 1.0;
    }
  }
}

/******************************************************************
 * Iterative computation of radiosity via Gathering; i.e. solution
 * using Gauss-Seidel iteration - reuse already computed values;
 * run-time O(n^2)
 *******************************************************************/

void Calculate_Radiosity(const int iteration) {
  const int n = num_tris;
  int patch_i = 0;

  for (int i = 0; i < n; i++) {
    for (int ia = 0; ia < tris[i].a_num; ia++) {
      for (int ib = 0; ib < tris[i].b_num; ib++) {
        Color B;

        int patch_j = 0;
        for (int j = 0; j < n; j++) {
          for (int ja = 0; ja < tris[j].a_num; ja++) {
            for (int jb = 0; jb < tris[j].b_num; jb++) {
              const double Fij = form_factor[patch_i * patch_num + patch_j];

              /* Add form factor multiplied with radiosity of previous step */
              if (Fij > 0.0)
                B = B + Fij * tris[j].patch[ja * tris[j].b_num + jb];

              patch_j++;
            }
          }
        }
        /* Multiply sum with color of patch and add emission */
        B = tris[i].color.MultComponents(B) + tris[i].emission;

        /* Store overall patch radiosity of current iteration */
        tris[i].patch[ia * tris[i].b_num + ib] = B;
        patch_i++;
      }
    }
  }
}

/******************************************************************
 * Helper functions for smooth bicubic (Catmull-Rom) interpolation
 * using 4x4 color patches;
 * First interpolate in y, followed by interpolation of results in x
 *******************************************************************/

Color cubicInterpolate(Color p[4], double x) {
  return p[1] + 0.5 * x *
                    (p[2] - p[0] +
                     x * (2.0 * p[0] - 5.0 * p[1] + 4.0 * p[2] - p[3] +
                          x * (3.0 * (p[1] - p[2]) + p[3] - p[0])));
}

Color bicubicInterpolate(Color p[4][4], double x, double y) {
  Color arr[4];

  arr[0] = cubicInterpolate(p[0], y);
  arr[1] = cubicInterpolate(p[1], y);
  arr[2] = cubicInterpolate(p[2], y);
  arr[3] = cubicInterpolate(p[3], y);

  return cubicInterpolate(arr, x);
}

/******************************************************************
 * Compute radiance from radiosity by shooting rays into the scene;
 * Radiance directly proportional to radiosity for assumed diffuse
 * emitters/surfaces (multiply by PI);
 * At intersections either constant patch color is returned or a
 * smoothly interpolated color of 4x4 neighboring patches
 *******************************************************************/

Color Radiance(const Ray &ray, const int depth, bool interpolation = true) {
  double t;
  int id;
  Vector normal;

  /* Find intersected triangle */
  if (!Intersect_Scene(ray, &t, &id, &normal)) {
    return BackgroundColor;
  }

  /* Determine intersection point on triangle */
  const Triangle &obj = tris[id];
  const Vector hitpoint = ray.org + t * ray.dir;

  /* Determine intersected patch */
  const Vector v = hitpoint - obj.p0;
  const double a_len = v.Dot(obj.edge_a.Normalized());
  const double b_len = v.Dot(obj.edge_b.Normalized());

  double da = obj.a_num * a_len / obj.a_len;
  double db = obj.b_num * b_len / obj.b_len;

  int ia = int(da);
  if (ia >= obj.a_num)
    ia--;
  int ib = int(db);
  if (ib >= obj.b_num)
    ib--;

  /* Bicubic interpolation for smooth image */
  if (interpolation) {
    Color c[4][4];

    int ia = int(da - 0.5);
    int ib = int(db - 0.5);

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        c[i][j] = obj.sample_patch(ia + i - 1, ib + j - 1);
      }
    }

    int ia0 = int(da - 0.5);
    int ib0 = int(db - 0.5);
    double dx = da - ia0 - 0.5;
    double dy = db - ib0 - 0.5;

    if (dx < 0.0)
      dx = 0.0;
    if (dx >= 1.0)
      dx = 1.0;
    if (dy < 0.0)
      dy = 0.0;
    if (dy >= 1.0)
      dy = 1.0;

    return bicubicInterpolate(c, dx, dy) * Over_M_PI;
  } else {
    return obj.patch[ia * obj.b_num + ib] * Over_M_PI;
  }
}

/******************************************************************
 * Main routine: Computation of radiosity image
 * Key parameters
 * - Image dimensions: width, height
 * - Number of samples for antialiasing (non-uniform filter): samples
 * - Number of patches along edges a,b: patches_a, patches_b
 * - Number of uniform samples per patch edge: MC_samples
 * - Number of iterations for iterative solver: iterations
 * Rendered result saved as PPM image file
 *******************************************************************/

int main(int argc, char **argv) {
  int width = 1280;
  int height = 720;
  int samples = 4;

  /* Set camera origin and viewing direction (negative z direction) */
  Ray camera(Vector(50.0, 52.0, 295.6),
             Vector(0.0, -0.042612, -1.0).Normalized());

  /* Image edge vectors for pixel sampling */
  Vector cx = Vector(width * 0.5135 / height);
  Vector cy = (cx.Cross(camera.dir)).Normalized() * 0.5135;

  /* Two final renderings; one with constant, one with interpolated patch colors
   */
  Image img(width, height);
  Image img_interpolated(width, height);

  cout << "Calculating form factors" << endl;

  /*
      The depth decides how many triangle-patches each triangle will be split
     into in the final rendering. A depth of n corresponds to every triangle
     being split into 4^n triangle-patches.
  */
  int depth = 0;

  for (int i = 1; i <= depth; i++) {
    tris = Split_Triangles(tris);
  }

  int patches_a = 1; // Should remain 1 now.
  int patches_b = 1; // Should remain 1 now.
  int MC_samples = 5;

  Calculate_Form_Factors(patches_a, patches_b, MC_samples);

  /* Iterative solution of radiosity linear system */
  cout << "Calculating radiosity" << endl;
  int iterations = 40;
  for (int i = 0; i < iterations; i++) {
    cout << i << " ";
    Calculate_Radiosity(i);
  }
  cout << endl;

  /* Loop over image rows */
  for (int y = 0; y < height; y++) {
    cout << "\rRendering (" << samples * 4 << " spp) "
         << (100.0 * y / (height - 1)) << "%     ";
    srand(y * y * y);

    /* Loop over row pixels */
    for (int x = 0; x < width; x++) {
      img.setColor(x, y, Color());
      img_interpolated.setColor(x, y, Color());

      /* 2x2 subsampling per pixel */
      for (int sy = 0; sy < 2; sy++) {
        for (int sx = 0; sx < 2; sx++) {
          Color accumulated_radiance = Color();
          Color accumulated_radiance2 = Color();

          /* Computes radiance at subpixel using multiple samples */
          for (int s = 0; s < samples; s++) {
            const double r1 = 2.0 * drand48();
            const double r2 = 2.0 * drand48();

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

            /* Determine constant radiance */
            accumulated_radiance =
                accumulated_radiance +
                Radiance(Ray(start, dir.Normalized()), 0, false) / samples;

            /* Determine interpolated radiance */
            accumulated_radiance2 =
                accumulated_radiance2 +
                Radiance(Ray(start, dir.Normalized()), 0, true) / samples;
          }

          img.addColor(x, y, accumulated_radiance);
          img_interpolated.addColor(x, y, accumulated_radiance2);
        }
      }
    }
  }

  cout << endl;

  img.Save(string("image_patches.ppm"));
  img_interpolated.Save(string("image_smooth.ppm"));
}
