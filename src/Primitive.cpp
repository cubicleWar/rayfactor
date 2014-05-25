#include "Primitive.h"

int Primitive::numThreads;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Primitive::Primitive
//
//	Comments : The constructor
//
//	Arguments:	obID is the ID of the object being created
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
Primitive::Primitive( int primitiveID )
{
    this->primitiveID = primitiveID;
    this->iden = _mm_castsi128_ps(_mm_set1_epi32(primitiveID)); // Make primitive id greater than 0
	isBoundingElement = false;
	rayDensity = kDefaultRayDensity;
	next = NULL;
	scaleVector = Vector(1.0f,1.0f,1.0f, 0.0f);
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Primitive::getID
//
//	Comments : Returns the objects ID
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int Primitive::getID()
{
    return primitiveID;
}

void Primitive::setNext(Primitive &n)
{
    next = &n;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Primitive::setBoundingElement
//
//	Comments : Set whether the object is a bounding object
//
//	Arguments: bE is a boolean value indicating whether the object should be set as a bounding object
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void Primitive::setIsBounding(bool bE)
{
	isBoundingElement = bE;
}

void Primitive::setRayDensity(int rD)
{
	rayDensity = rD;
}

void Primitive::setWillAnalyse( bool willAnalyse )
{
	this->willAnalysePrimitive = willAnalyse;
}





/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Primitive::rotate
//
//	Comments : Adds the a rotation affine transformation to the affine transformation stack of this 
//			   object
//
//	Arguments: angle is the angle to rotate this object
//			   u is a vector containing the axis to rotate the object about eg. (0,1,0) for a y rotation
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void Primitive::rotate(float angle, Vector u)
{
    affineTransformation::rotate(affine, angle, u);
    affineTransformation::invRotate(invAffine, angle, u);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Primitive::Scale
//
//	Comments : Adds a scale affine transformation to the affine transformation stack
//
//	Arguments: sx the scale percentage along the x-axis
//			   sy the scale percentage along the y-axis
//			   sx the scale percentage along the z-axis
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void Primitive::scale(float sx, float sy, float sz)
{
    affineTransformation::scale(affine, sx, sy, sz);
    affineTransformation::invScale(invAffine, sx, sy, sz);

	scaleVector.set(sx,sy,sz);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Primitive::translate
//
//	Comments : Adds a translation affine transformation to the affine transformation stack
//
//	Arguments: d is a vector by which this object should be translated
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void Primitive::translate(Vector d)
{
    affineTransformation::translate(affine, d);
    affineTransformation::invTranslate(invAffine, d);
    //aS->translate(d);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Primitive::getFirstHit
//
//	Comments : Determines the object which is hit first by a ray travelling from this objects surface
//			   into the surrounding environment
//
//	Arguments: head is the start of a the list of scene objects in the current scene
//			   ray is the to be traced
//			   Intersection is a intersection record for this ray
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void Primitive::getFirstHit(Primitive* head, Rayx4& ray, __m128 &objectNo)
{
    //float tbest = inf;
    //float t = 0.0f;
    //__m128 t = *(__m128*)_ps_inf;
    //__m128 tbest = *(__m128*)_ps_inf;
    __m128 tbest = _mm_set1_ps(100000.f);
    
    //__m128 mask;
    objectNo = _mm_castsi128_ps(_mm_set1_epi32((int)-1));

    for(Primitive* pobj = head; pobj != NULL; pobj = pobj->next)
    {
        pobj->hit(ray, tbest, objectNo);
        //mask = _mm_cmplt_ps(t, tbest);
        
        //tbest = _mm_blendv_ps(tbest, t,mask);
        
        //objectNo = _mm_blendv_ps(objectNo,pobj->iden,mask);
    }
}



