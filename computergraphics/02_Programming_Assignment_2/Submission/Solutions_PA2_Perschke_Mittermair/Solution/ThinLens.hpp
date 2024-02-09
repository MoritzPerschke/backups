#include "PathTracing.hpp"

class ThinLens{
    private:

    public:
        Vector position;
        double radius;
        double focalLength;

        ThinLens(const Vector& position, double radius, double f)
            : position(position), radius(radius), focalLength(f){}
            
        Vector sample(){
            double xi1 = 2.0 * M_PI * drand48();
            double xi2 = sqrt(drand48());

            Vector offset(xi2 * cos(xi1), + xi2 * sin(xi1), 0.0);
            return offset;
        }
};

class Camera{
    private:

    public:
        Vector position;
        Vector direction;

        Camera(const Vector& position, const Vector& direction)
            :position(position), direction(direction.Normalized()){};
};