/******************************************************************
*
****************************************************
* Team Members: Thomas Mittermair, Moritz Perschke *
****************************************************
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

#include "ThinLens.hpp"

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
    // this blocks all light when using a lens 
    // Sphere( 1e5, Vector(      50,      40.8, -1e5 +170),  Vector(), Vector(),            DIFF), /* Front wall */
    Sphere( 1e5, Vector(      50,       1e5,      81.6),  Vector(), Vector(.75,.75,.75), DIFF), /* Floor */
    Sphere( 1e5, Vector(      50,-1e5 +81.6,      81.6),  Vector(), Vector(.75,.75,.75), DIFF), /* Ceiling */

    Sphere(16.5, Vector(20, 16.5, 60), Vector(), Vector(1,1,1)*.999,  REFR), /* Mirror sphere */
    Sphere(16.5, Vector(75, 16.5, 60), Vector(), Vector(1,1,1)*.999,  SPEC), /* Glas sphere */
    Sphere(16.5, Vector(20, 16.5 + 2 * 16.5, 60), Vector(), Vector(1,1,1)*.999,  TRAN), /* Mirror sphere */
    Sphere(16.5, Vector(75, 16.5 + 2 * 16.5, 60), Vector(), Vector(1,1,1)*.999,  GLOS), /* Glas sphere */

    Sphere( 1.5, Vector(50, 81.6-16.5, 81.6), Vector(4,4,4)*100, Vector(), DIFF), /* Light */
    Sphere( 0.5, Vector(50, 17, 110), Vector(4,4,4)*100, Vector(), DIFF), /* New Light */
};


Primitive **primitives;      // Array of pointer to Primitive-objects.
int num_primitives = 0;

Primitive **Create_Primitives() {
    std::vector<Primitive *> *all_primitives =
      new vector<Primitive *>; // Allocate vector on the heap, to use it after
                            // the function.

    const int num_spheres = int(sizeof(spheres) / sizeof(Sphere));

    for (int i = 0; i < num_spheres; i ++) 
    {
        all_primitives->push_back(&(spheres[i]));
    }

    // Add all triangles to the list:

    // Code from official GitHub: https://github.com/tinyobjloader/tinyobjloader
    //std::string inputfile = "../../assets/jeep_lowpoly/jeep.obj";
    std::string inputfile = "../../assets/cornell-box/CornellBox-Original.obj";
    tinyobj::ObjReaderConfig reader_config;
    reader_config.triangulate = true;
    reader_config.mtl_search_path = "../../assets/cornell-box/"; // Path to material files

    tinyobj::ObjReader reader;

    if (!reader.ParseFromFile(inputfile, reader_config)) {
        if (!reader.Error().empty()) {
        std::cerr << "TinyObjReader: " << reader.Error();
        }
        exit(1);
    }

    if (!reader.Warning().empty()) {
        std::cout << "TinyObjReader: " << reader.Warning();
    }

    auto& attrib = reader.GetAttrib();
    auto& shapes = reader.GetShapes();
    auto& materials = reader.GetMaterials();

    // Loop over shapes
    for (size_t s = 0; s < shapes.size(); s++) {
        // Loop over faces(polygon)
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);

            std::vector<Vector> *points = new std::vector<Vector>;      // Points of the triangle.

            // Loop over vertices in the face.
            for (size_t v = 0; v < fv; v++) {
                // access to vertex
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                tinyobj::real_t vx = attrib.vertices[3*size_t(idx.vertex_index)+0];
                tinyobj::real_t vy = attrib.vertices[3*size_t(idx.vertex_index)+1];
                tinyobj::real_t vz = attrib.vertices[3*size_t(idx.vertex_index)+2];

                Vector vec = Vector(vx, vy, vz);

                /*
                    We now apply transformations (scaling and translation) onto the vector. This way, the final mesh will be
                    at the correct position in the scene.
                */
                Vector transformedVec = 10 * vec + Vector(50, 0, 110);     // (right, up, front)

                points->push_back(transformedVec);

                /*
                    // Check if `normal_index` is zero or positive. negative = no normal data
                    if (idx.normal_index >= 0) {
                        tinyobj::real_t nx = attrib.normals[3*size_t(idx.normal_index)+0];
                        tinyobj::real_t ny = attrib.normals[3*size_t(idx.normal_index)+1];
                        tinyobj::real_t nz = attrib.normals[3*size_t(idx.normal_index)+2];
                    }

                    // Check if `texcoord_index` is zero or positive. negative = no texcoord data
                    if (idx.texcoord_index >= 0) {
                        tinyobj::real_t tx = attrib.texcoords[2*size_t(idx.texcoord_index)+0];
                        tinyobj::real_t ty = attrib.texcoords[2*size_t(idx.texcoord_index)+1];
                    }

                    // Optional: vertex colors
                    tinyobj::real_t red   = attrib.colors[3*size_t(idx.vertex_index)+0];
                    tinyobj::real_t green = attrib.colors[3*size_t(idx.vertex_index)+1];
                    tinyobj::real_t blue  = attrib.colors[3*size_t(idx.vertex_index)+2];
                */
            }

            index_offset += fv;

            // per-face material
            shapes[s].mesh.material_ids[f];

            unsigned int materialId = shapes[s].mesh.material_ids[f];
            Vector faceColor = Vector(materials[materialId].diffuse[0], materials[materialId].diffuse[1], materials[materialId].diffuse[2]);

            Vector emissionColor = Vector(materials[materialId].emission[0], materials[materialId].emission[1], materials[materialId].emission[2]);

            /*
                At this point, points holds the 3 vertices of the triangle. From these 3 points, the normal is automatically computed.
                faceColor holds the color.
            */

            Triangle *t = new Triangle(points->at(0), points->at(1) - points->at(0), points->at(2) - points->at(0),
                emissionColor, faceColor, DIFF);

            all_primitives->push_back(t);
        }
    }

    num_primitives = all_primitives->size();

    return &(*all_primitives)[0]; // Convert vector<> to array and return it.
}


/******************************************************************
* Check for closest intersection of a ray with the scene;
* returns true if intersection is found, as well as ray parameter
* of intersection and id of intersected object
*******************************************************************/
bool Intersect(const Ray &ray, double &t, int &id)
{
    //const int n = int(sizeof(spheres) / sizeof(Sphere));
    const int n = num_primitives;
    t = 1e20;

    for (int i = 0; i < n; i ++) 
    {
        double d = primitives[i]->Intersect(ray);
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

    int numPrimitives = num_primitives;

    double t;                               
    int id = 0;  
                             
    if (!Intersect(ray, t, id))   /* No intersection with scene */
        return BackgroundColor; 

    const Primitive &obj = *(primitives[id]);    

    Vector hitpoint = ray.org + ray.dir * t;    /* Intersection point */

    Vector normal;
    normal = obj.getNormal(hitpoint);

    Vector nl = normal;

    /* Obtain flipped normal, if object hit from inside */
    if (normal.Dot(ray.dir) >= 0) 
        nl = nl*-1.0;

    Color col = obj.color; 

    /* Maximum RGB reflectivity for Russian Roulette */
    double p = col.Max();

    if (depth > 5 || !p)   /* After 5 bounces or if max reflectivity is zero */
    {
        if (drand48() < p)            /* Russian Roulette */
            col = col * (1/p);        /* Scale estimator to remain unbiased */
        else 
            return obj.emission * E;  /* No further bounces, only return potential emission */
    }

    // enum Refl_t { DIFF, SPEC, REFR, GLOS, TRAN }; 
    if (obj.refl == SPEC || obj.refl == GLOS) {
        /* Return light emission mirror reflection (via recursive call using perfect
        reflection vector) */

        Vector perfectReflDir = ray.dir - normal * 2 * normal.Dot(ray.dir);
        Vector reflDir;

        if (obj.refl == SPEC) {
            reflDir = perfectReflDir;
        } else {
            Vector sw = perfectReflDir.Normalized();
            Vector su;
            
            if(fabs(sw.x) > 0.1)
                su = Vector(0.0, 1.0, 0.0);
            else
                su = Vector(1.0, 0.0, 0.0);

            su = (su.Cross(nl)).Normalized();
            Vector sv = sw.Cross(su);

            /* Create random sample direction l towards spherical light source */
            double cos_a_max = 0.975;     // Parameter that controlls glossyness.
            double eps1 = drand48();
            double eps2 = drand48();
            double cos_a = 1.0 - eps1 + eps1 * cos_a_max;
            double sin_a = sqrt(1.0 - cos_a * cos_a);
            double phi = 2.0*M_PI * eps2;
            Vector randomReflDir = su * cos(phi) * sin_a + 
                        sv * sin(phi) * sin_a + 
                        sw * cos_a;
            randomReflDir = randomReflDir.Normalized();
            // randomReflDir is a random vector in a cone around the direction of the mirror direction.
            reflDir = randomReflDir;
        }
        return obj.emission + 
            col.MultComponents(Radiance(Ray(hitpoint, reflDir),
                            depth, 1));
    } else if (obj.refl == REFR || obj.refl == TRAN) {
        /* Otherwise object transparent, i.e. assumed dielectric glass material */
        // Use original normal for the ideal mirror reflection:
        Ray reflRay (hitpoint, ray.dir - normal * 2 * normal.Dot(ray.dir));  /* Prefect reflection */ 

        if (obj.refl == TRAN) { 
            /*
                Instead of perturbing the refraction vector, we perturb the normal:
                The normal is set to a random vector in a cone around the direction of the real normal.
            */
            Vector sw = normal.Normalized();
            Vector su;
            
            if(fabs(sw.x) > 0.1)
                su = Vector(0.0, 1.0, 0.0);
            else
                su = Vector(1.0, 0.0, 0.0);

            su = (su.Cross(nl)).Normalized();
            Vector sv = sw.Cross(su);

            // Create random sample direction l towards spherical light source 
            double cos_a_max = 0.975;     // Parameter that controlls glossyness.
            double eps1 = drand48();
            double eps2 = drand48();
            double cos_a = 1.0 - eps1 + eps1 * cos_a_max;
            double sin_a = sqrt(1.0 - cos_a * cos_a);
            double phi = 2.0*M_PI * eps2;
            normal = su * cos(phi) * sin_a + 
                        sv * sin(phi) * sin_a + 
                        sw * cos_a;
            normal = normal.Normalized();
            // randomReflDir is a random vector in a cone around the direction of the mirror direction.

            nl = normal;

            if (normal.Dot(ray.dir) >= 0) 
                nl = nl*-1.0;
        }

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
            if (drand48() < P)
                return obj.emission + col.MultComponents(Radiance(reflRay, depth, 1) * RP);
            else
                return obj.emission + col.MultComponents(Radiance(Ray(hitpoint,tdir), depth, 1) * TP);
    } else {    // obj.refl == DIFF
        /* Compute random reflection vector on hemisphere */
        double r1 = 2.0 * M_PI * drand48(); 
        double r2 = drand48(); 
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
        for (int i = 0; i < numPrimitives; i ++)
        {
            //const Sphere &sphere = spheres[i];
            const Primitive &primitive = *(primitives[i]);

            if (primitive.emission.x <= 0 && primitive.emission.y <= 0 && primitive.emission.z <= 0) 
                continue; /* Skip objects that are not light sources */
      
            /* Randomly sample spherical light source from surface intersection */

            /* Set up local orthogonal coordinate system su,sv,sw towards light source */
            Vector sw = primitive.getPosition() - hitpoint;
            Vector su;
            
            if(fabs(sw.x) > 0.1)
                su = Vector(0.0, 1.0, 0.0);
            else
                su = Vector(1.0, 0.0, 0.0);

            su = (su.Cross(w)).Normalized();
            Vector sv = sw.Cross(su);

            /* Create random sample direction l towards spherical light source */
            double cos_a_max = sqrt(1.0 - primitive.getRadius() * primitive.getRadius() / 
                                    (hitpoint - primitive.getPosition()).Dot(hitpoint-primitive.getPosition()));
            double eps1 = drand48();
            double eps2 = drand48();
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
                e = e + col.MultComponents(primitive.emission * l.Dot(nl) * omega) * M_1_PI; 
            }
        }
   
        /* Return potential light emission, direct lighting, and indirect lighting (via
           recursive call for Monte-Carlo integration */      
        return obj.emission * E + e + col.MultComponents(Radiance(Ray(hitpoint,d), depth, 0));
    }
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
    double aperture = 5;

    /* Set camera origin and viewing direction (negative z direction) */
    Vector camera_org = Vector(50.0, 52.0, 295.6);
    Vector camera_dir = Vector(0.0, -0.042612, -1.0);

    /* center beween refractive and glossy sphere in focus*/
    Vector center = (Vector(20, 16.5, 60) + Vector(75, 16.5 + 2 * 16.5, 60)) / 2;

    double f = (center - camera_org).Length();

    switch (argc)
    {
    case 1:
        std::cout << "Change default amount samples, aperture and focal length with './Pathtracing <samples> <aperture> <focal length>'" << std::endl;
        break;
    case 2:
        samples = atoi(argv[1]);
        break;
    case 3:
        samples = atoi(argv[1]);
        aperture = atof(argv[2]);
    case 4:
        samples = atoi(argv[1]);
        aperture = atof(argv[2]);
        f = atof(argv[3]);
    default:
        break;
    }
     
    /* glass sphere in focus */
    // double f = (Vector(73, 16.5, 78) - camera_org).Length();
    
    /* mirror sphere in focus */
    // double f = (Vector(27, 16.5, 47) - camera_org).Length();

    /* translucent sphere in focus*/
    // double f = Vector(20, 16.5 + 2 * 16.5, 60) - camera_org).Length();
    
    /* neither in focus*/
    // double f = 251.;

    ThinLens lens = ThinLens(camera_org, aperture, f);
    Camera camera = Camera(camera_org, camera_dir);

    /* Image edge vectors for pixel sampling */
    Vector cx = Vector(width * 0.5135 / height);
    Vector cy = (cx.Cross(camera.direction)).Normalized() * 0.5135;

    /* Final rendering */
    Image img(width, height);
    primitives = Create_Primitives();

    /* Loop over image rows */
    for (int y = 0; y < height; y ++) 
    { 
        cout << "\rRendering (" << samples * 4 << " spp) " << (100.0 * y / (height - 1)) << "%     ";
        srand(y * y * y);
 
        /* Loop over row pixels */
        for (int x = 0; x < width; x ++)  
        {
            img.setColor(x, y, Color());
            
            /* 2x2 subsampling per pixel */
            for (int sy = 0; sy < 2; sy ++) 
            {
                for (int sx = 0; sx < 2; sx ++) 
                {
                    Color accumulated_radiance = Color();

                    /******************************************************************
                     * According to lecture, steps for thin lens "distribution ray    *
                     * tracing":                                                      *
                     * 1. Cast Primary ray through image plane from lens center       *
                     * 2. Find focus point as intersection with focal plane           *
                     * 3. sample positions on lens and cast multiple rays             *
                     ******************************************************************/
                    
                    // 1.
                    Vector dir = cx * ((x + (sx + 0.5) / 2.0) / width - 0.5) +
                                    cy * ((y + (sy + 0.5) / 2.0) / height - 0.5) +
                                    camera.direction;
                    // 2.
                    Vector focalPoint = lens.position + lens.focalLength * dir;

                    // 3.
                    /* Compute radiance at subpixel using multiple samples */
                    for (int s = 0; s < samples; s ++) {
                        // 3.
                        Vector orgp = lens.position + lens.radius * lens.sample();
                        Vector dirp = (focalPoint - orgp).Normalized();

                        /* Accumulate radiance */
                        accumulated_radiance = accumulated_radiance + 
                            Radiance( Ray(orgp, dirp), 0, 1) / samples;
                    } 
                    
                    accumulated_radiance = accumulated_radiance.clamp() * 0.25;

                    img.addColor(x, y, accumulated_radiance);
                }
            }
        }
    }
    cout << endl;

    img.Save(string("image.ppm"));
}