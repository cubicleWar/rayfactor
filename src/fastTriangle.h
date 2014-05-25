

#ifndef RayFactor_fastTriangle_h
#define RayFactor_fastTriangle_h

#include "Primitive.h"

using namespace std;

class fastTriangle : public Primitive
{
private:
    
public:

    float cx, cy, cz, cw;
    float nx, ny, nz, nd;
    float ux, uy, uz, ud;
    float vx, vy, vz, vd;
    
    fastTriangle(int obID) : Primitive(obID) {};
    void hit(const Rayx4& r, __m128& tbest, __m128& objectNo);
    void traceFactors(Primitive *head, VFMatrix* vfm);
    
};
#endif
