#include "Scene.h"



/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Scene::Scene
//
//	Comments : Default Constructor
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
Scene::Scene()
{
    head = NULL;
    tail = NULL;
	numPrimitives = 0;
}

//Note put a deconstructor: see free Scene

#pragma mark readScene

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Scene::readScene
//
//	Comments : The function for reading an xml NDRay project file and creating the scene
//
//	Arguments : filename is a c string that holds the file name of the xml project file
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	23/01/10	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void Scene::readScene(const char* filename)
{
	TiXmlDocument doc(filename);
	bool didLoadOk = doc.LoadFile();
	if(didLoadOk) {
		TiXmlHandle docHandle( &doc );
		
		TiXmlElement *settings = docHandle.FirstChild( "project" ).FirstChild( "settings" ).ToElement();
		TiXmlElement *geometry = docHandle.FirstChild( "project" ).FirstChild( "geometry" ).ToElement();

		this->loadSettings( settings );
		this->loadGeometry( geometry );
	} else {
		printf("Failed to load file \"%s\"\n", filename);
	}
	
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Scene::loadSettings
//
//	Comments : Parses the project settings from the <settings> xml block
//
//	Arguments : settings is a pointer to a TiXmlElement representing the <settings> block.
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	23/01/10	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Scene::loadSettings( TiXmlElement *settings )
{
	TiXmlElement *child = settings->FirstChild()->ToElement();
	
	//For each element in the <settings> block
	for( child; child; child=child->NextSiblingElement())
	{
		string typeStr = child->ValueStr();
		if(typeStr.compare(kGlobalRayDensity) == 0) 
		{
			child->QueryValueAttribute(kValue,&globalRayDensity);
		} 
		else if(typeStr.compare(kDescription) == 0) 
		{
			sceneDescription = string(child->GetText());
		}
	}
	
	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Scene::loadGeometry
//
//	Comments : Is the higher level handler for the primitives stored in the <geometry> xml block. It 
//			   will determine what type of primitive is specified, create the primitive object, extract
//			   any unique attributes than pass it on for parsing of the generic primitive properties
//
//	Arguments : geometry is pointer to a TiXmlElement which contains the <geometry> block.
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	23/01/10	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Scene::loadGeometry( TiXmlElement *geometry )
{
	TiXmlElement *primitive = geometry->FirstChild()->ToElement();
	int objectID = 0;
	
	//For each <primitive> block
	for( primitive; primitive; primitive=primitive->NextSiblingElement() )
	{
		string prType = "";
		
		bool willAnalyse = true;
		bool isBounding = false;
		
		if(primitive->Attribute( "type") != NULL) {
			prType = primitive->Attribute("type");
		}
		
		if(primitive->Attribute("isBounding") != NULL) {
			string isBoundingStr = primitive->Attribute("isBounding");
			isBounding = (isBoundingStr.compare("true") == 0);
		}
		
		if(primitive->Attribute("analyse") != NULL) {
			string analyseStr = primitive->Attribute("analyse");
			willAnalyse = (analyseStr.compare("true") == 0);
		}
		
		if(prType.compare("cylinderSurface") == 0 || prType.compare("taperedCylinderSurface") == 0) {
			Primitive *cylinder;
			float smallRadius;
			primitive->FirstChildElement( "smallRadius" )->QueryValueAttribute( "value", &smallRadius );
            
            // Handle whether to use a faster primitive based on its type i.e cone, cylinder, frustum
            
			if(smallRadius == 1.0) {
                cylinder = new CylinderSurface(objectID);
            } else {
                cylinder = new TaperedCylinderSurface(objectID, smallRadius);
            }
            
			
			cylinder->setWillAnalyse(willAnalyse);
			cylinder->setIsBounding(isBounding);
			
			this->initPrimitive( primitive, *cylinder);

			this->addObject( *cylinder );
			
		} else if(prType.compare( "annulus" ) == 0) {
			
			float smallRadius, largeRadius;
			primitive->FirstChildElement( "smallRadius" )->QueryValueAttribute( "value", &smallRadius );
			primitive->FirstChildElement( "largeRadius" )->QueryValueAttribute( "value", &largeRadius );
			
			//Assumes largeRaidus and smallRadius are defined and correct
			Annulus *annulus = new Annulus(objectID, largeRadius, smallRadius);
			
			annulus->setWillAnalyse(willAnalyse);
			annulus->setIsBounding(isBounding);

			this->initPrimitive( primitive, *annulus);
			
			this->addObject( *annulus );
			
		}
		else if(prType.compare( "rectangle" ) == 0) {
			
			Rectangle *rectangle = new Rectangle(objectID);
			
			rectangle->setWillAnalyse(willAnalyse);

			this->initPrimitive( primitive, *rectangle);
			
			this->addObject( *rectangle );

		}
		else if(prType.compare( "disc" ) == 0) {
			
			Disc *disc = new Disc(objectID);
			
			disc->setWillAnalyse(willAnalyse);

			this->initPrimitive( primitive, *disc);
			
			this->addObject( *disc );
			
		} else if(prType.compare( "sphere" ) == 0) {
			
			Sphere *sphere = new Sphere(objectID);
			
			sphere->setWillAnalyse(willAnalyse);
			sphere->setIsBounding(isBounding);
			
			this->initPrimitive( primitive, *sphere);
			
			this->addObject( *sphere );
			
		} else if(prType.compare( "triangle" ) == 0) {
            Triangle *triangle = new Triangle(objectID);
            
            triangle->setWillAnalyse(willAnalyse);

            // Get the verticies
            /*
            AB = B - A     (3D Point)
            AC = C - A     (3D Point)
             
             N = AB x AC
             nd =  A.N
             
             U = AC x N / |N|^2
             Ud = -U.A
             
             V = N x AB / |N|^2
             Vd = -V.A
            
            */
            float Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz;
            
            
            TiXmlElement *A = primitive->FirstChild( "A" )->ToElement();
            
            A->QueryValueAttribute( "x", &Ax );
            A->QueryValueAttribute( "y", &Ay );
            A->QueryValueAttribute( "z", &Az );

            TiXmlElement *B = primitive->FirstChild( "B" )->ToElement();
            
            B->QueryValueAttribute( "x", &Bx );
            B->QueryValueAttribute( "y", &By );
            B->QueryValueAttribute( "z", &Bz );

            TiXmlElement *C = primitive->FirstChild( "C" )->ToElement();
            
            C->QueryValueAttribute( "x", &Cx );
            C->QueryValueAttribute( "y", &Cy );
            C->QueryValueAttribute( "z", &Cz );

            
            triangle->Ax = Ax;
            triangle->Ay = Ay;
            triangle->Az = Az;
            
            float ABx = Bx - Ax;
            float ABy = By - Ay;
            float ABz = Bz - Az;
            
            
            float ACx = Cx - Ax;
            float ACy = Cy - Ay;
            float ACz = Cz - Az;
            
            triangle->ABx = ABx;
            triangle->ABy = ABy;
            triangle->ABz = ABz;
            
            triangle->BCx = Cx - Bx;
            triangle->BCy = Cy - By;
            triangle->BCz = Cz - Bz;
            
            triangle->nx = ABy*ACz - ABz*ACy;
            triangle->ny = ABz*ACx - ABx*ACz;
            triangle->nz = ABx*ACy - ABy*ACx;
            triangle->nd = triangle->nx*Ax + triangle->ny*Ay + triangle->nz*Az;
            
            float mgN = sqrtf(triangle->nx*triangle->nx + triangle->ny*triangle->ny + triangle->nz*triangle->nz);
            float mgNsq = mgN*mgN;
        
            triangle->Aw = 1.0f/mgN;
            
            triangle->ux = (ACy*triangle->nz - ACz*triangle->ny)/mgNsq;
            triangle->uy = (ACz*triangle->nx - ACx*triangle->nz)/mgNsq;
            triangle->uz = (ACx*triangle->ny - ACy*triangle->nx)/mgNsq;
            triangle->ud = -triangle->ux*Ax - triangle->uy*Ay - triangle->uz*Az;
            
            triangle->vx = (triangle->ny*ABz -  triangle->nz*ABy)/mgNsq;
            triangle->vy = (triangle->nz*ABx - triangle->nx*ABz)/mgNsq;
            triangle->vz = (triangle->nx*ABy - triangle->ny*ABx)/mgNsq;
            triangle->vd = -triangle->vx*Ax - triangle->vy*Ay - triangle->vz*Az;

            
            this->initPrimitive( primitive, *triangle);

			this->addObject( *triangle );
        } else if(prType.compare( "fasttriangle" ) == 0) {
            fastTriangle *fTriangle = new fastTriangle(objectID);
            
            fTriangle->setWillAnalyse(willAnalyse);
            
            float Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz;
            
            
            TiXmlElement *A = primitive->FirstChild( "A" )->ToElement();
            
            A->QueryValueAttribute( "x", &Ax );
            A->QueryValueAttribute( "y", &Ay );
            A->QueryValueAttribute( "z", &Az );
            
            TiXmlElement *B = primitive->FirstChild( "B" )->ToElement();
            
            B->QueryValueAttribute( "x", &Bx );
            B->QueryValueAttribute( "y", &By );
            B->QueryValueAttribute( "z", &Bz );
            
            TiXmlElement *C = primitive->FirstChild( "C" )->ToElement();
            
            C->QueryValueAttribute( "x", &Cx );
            C->QueryValueAttribute( "y", &Cy );
            C->QueryValueAttribute( "z", &Cz );
            
            
            float ACx = Cx - Ax;
            float ACy = Cy - Ay;
            float ACz = Cz - Az;
            
            float ABx = Bx - Ax;
            float ABy = By - Ay;
            float ABz = Bz - Az;
            
            fTriangle->cx = (Ax + Bx + Cx)/3.0f;
            fTriangle->cy = (Ay + By + Cy)/3.0f;
            fTriangle->cz = (Az + Bz + Cz)/3.0f;
            
            fTriangle->nx = ABy*ACz - ABz*ACy;
            fTriangle->ny = ABz*ACx - ABx*ACz;
            fTriangle->nz = ABx*ACy - ABy*ACx;
            fTriangle->nd = fTriangle->nx*Ax + fTriangle->ny*Ay + fTriangle->nz*Az;
            
            float mgN = sqrtf(fTriangle->nx*fTriangle->nx + fTriangle->ny*fTriangle->ny + fTriangle->nz*fTriangle->nz);
            float mgNsq = mgN*mgN;
            
            fTriangle->cw = 1.0f/mgN;
            
            fTriangle->ux = (ACy*fTriangle->nz - ACz*fTriangle->ny)/mgNsq;
            fTriangle->uy = (ACz*fTriangle->nx - ACx*fTriangle->nz)/mgNsq;
            fTriangle->uz = (ACx*fTriangle->ny - ACy*fTriangle->nx)/mgNsq;
            fTriangle->ud = -fTriangle->ux*Ax - fTriangle->uy*Ay - fTriangle->uz*Az;
            
            fTriangle->vx = (fTriangle->ny*ABz - fTriangle->nz*ABy)/mgNsq;
            fTriangle->vy = (fTriangle->nz*ABx - fTriangle->nx*ABz)/mgNsq;
            fTriangle->vz = (fTriangle->nx*ABy - fTriangle->ny*ABx)/mgNsq;
            fTriangle->vd = -fTriangle->vx*Ax - fTriangle->vy*Ay - fTriangle->vz*Az;
            
            this->initPrimitive( primitive, *fTriangle);
            
			this->addObject( *fTriangle );

        }
		

		objectID++;
	}
	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Scene::initPrimitive
//
//	Comments : Parsers all the attributes common to all primitives and applies them to the primitive 
//			   object
//
//	Arguments : xmlElement is a pointer to a TiXmlElement representing the <primitive> block.
//				primative is a reference to a ScnObject to which the parsed attributes will be applied.
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	23/01/10	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Scene::initPrimitive( TiXmlElement *xmlElement, Primitive &primitive )
{
	if(!(xmlElement->FirstChild())) {
		return true;
	}
	TiXmlElement *child = xmlElement->FirstChild()->ToElement();
	bool foundRayDensity = false;
	// For each child in a <primitive> block
	for( child; child; child=child->NextSiblingElement())
	{
		string typeStr = child->ValueStr();

		//Should be done in the order of scale, rotate, translate
		if(typeStr.compare("rayDensity") == 0) 
		{
			int rayDensity;
			child->Attribute( "value", &rayDensity);
			primitive.setRayDensity(rayDensity);
			foundRayDensity = true;	
		}
		else if(typeStr.compare("translation") == 0) 
		{
			float xTrans, yTrans, zTrans;
			
			child->QueryValueAttribute( "x", &xTrans );
			child->QueryValueAttribute( "y", &yTrans );
			child->QueryValueAttribute( "z", &zTrans );
			
			primitive.translate(Vector(xTrans,yTrans,zTrans));
		} 
		else if(typeStr.compare("scale") == 0) 
		{
			float xScale, yScale, zScale;
			
			child->QueryValueAttribute( "x", &xScale );
			child->QueryValueAttribute( "y", &yScale );
			child->QueryValueAttribute( "z", &zScale );
			
			primitive.scale(xScale,yScale,zScale);
		} 
		else if(typeStr.compare("rotate") == 0) 
		{
			float degRot, xAxis, yAxis, zAxis;
			
			child->QueryValueAttribute( "degrees", &degRot );
			child->QueryValueAttribute( "x", &xAxis );
			child->QueryValueAttribute( "y", &yAxis );
			child->QueryValueAttribute( "z", &zAxis );
			
			primitive.rotate(degRot, Vector(xAxis,yAxis,zAxis));
		}
	}
	
	//Set the scene defaults if no specific attributes were found
	if(!foundRayDensity) {
		primitive.setRayDensity(globalRayDensity);
	}

	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Scene::addObject
//
//	Comments : Adds an object to the system to be analysed
//
//	Arguments: so is the object to be added to the current system
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Scene::addObject(Primitive &so)
{

    if(head == NULL)
    {
        head = &so;
        tail = &so;
		numPrimitives++;
        return true;
    }
    else
    {
        tail->setNext(so);
        tail = &so;
		numPrimitives++;
        return true;
    }
	
	return false;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Scene::getObject
//
//	Comments : Returns a pointer to a scene object with a given ID
//
//	Arguments: iden is the ID of the scene object to return
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
Primitive* Scene::getObject(int obID) //is this function used?
{
    for(Primitive *so = head; so != NULL; so = so->next)
    {
        if(so->getID() == obID)
            return so;
    }
	return NULL;
}

float Scene::getNumPrimitives()
{
	return numPrimitives;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Scene::findViewFactors
//
//	Comments : Driver routine to find to process each object and calculate the view factor matrix
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void Scene::findViewFactors()
{
	int currSurf = 1;

	VFMatrix *vfm = new VFMatrix(numPrimitives);
	
	//Seed the random number generator
	//srand((unsigned)time(NULL));
	
	
	struct timeval startTime, endTime;
	
	for(Primitive* pobj = head; pobj != NULL; pobj = pobj->next)
	{
		cout << "Processing Object " << currSurf << " of " << numPrimitives << endl;
		//pobj->dirAffine->print();
		gettimeofday(&startTime, NULL);
		
		if(pobj->willAnalyse()) {
			pobj->traceFactors(head, vfm);
		}
		
		
		gettimeofday(&endTime, NULL);
		float t = (1000*(endTime.tv_sec-startTime.tv_sec)+(endTime.tv_usec-startTime.tv_usec)/1000);
		
		std::cout << "Time " << t << " ms" << endl;
		
		currSurf++;
	}
	
	vfm->print();
	
    return;
}
