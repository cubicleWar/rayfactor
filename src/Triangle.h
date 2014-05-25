#ifndef RayFactor_Triangle_h
#define RayFactor_Triangle_h

#include "Primitive.h"

using namespace std;

class Triangle : public Primitive
{
private:
    
public:
    float Ax, Ay, Az, Aw;       // Note Aw, ABw and BCw are spare floats not used
    float ABx, ABy, ABz, ABw;
    float BCx, BCy, BCz, BCw;
    //float cx, cy, cz, cd;
    float nx, ny, nz, nd;
    float ux, uy, uz, ud;
    float vx, vy, vz, vd;
    
    Triangle(int obID) : Primitive(obID) {};
    void hit(const Rayx4& r, __m128& tbest, __m128& objectNo);
    void traceFactors(Primitive *head, VFMatrix* vfm);
    
};

#endif
