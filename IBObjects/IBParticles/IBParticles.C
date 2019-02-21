/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBParticles class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "IBParticles.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"

// * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IBParticles, 0);
    addToRunTimeSelectionTable(IBObject, IBParticles, dictionary);


// * * * * * * * * * * * * Private Member Fuctions * * * * * * * * * * * * * //

void IBParticles::readTXTFile(word fileName, point center, scalar R)
{
	ifstream datafile(fileName.c_str());
    if (datafile.is_open())
    {
        //-Read points
        datafile>>nPoints_;
        lPoints_.setSize(nPoints_, vector::zero);
        scalar x;
        scalar y;
        scalar z;
        for (int pointI=0;pointI<nPoints_;pointI++)
        {
            datafile>>x;
            datafile>>y;
            datafile>>z;

            lPoints_[pointI].x()=R*x+center.x();
            lPoints_[pointI].y()=R*y+center.y();
            lPoints_[pointI].z()=R*z+center.z();
        }

        //-Read faces connectivity:
        datafile>>nFaces_;
        nPointsOfFaces_.setSize(nFaces_,3);
        pointOfFace_.setSize(nFaces_);
        label point1;
        label point2;
        label point3;
        for (int i=0;i<nFaces_;i++)
        {
            labelList PointsOfFaceI(0);
            datafile>>point1;
            datafile>>point2;
            datafile>>point3;

            PointsOfFaceI.append(point1); 
            PointsOfFaceI.append(point2);
            PointsOfFaceI.append(point3);
            pointOfFace_[i] = PointsOfFaceI;
            nPointsOfFaces_[i] = pointOfFace_[i].size();
        }
    }
    datafile.close();
}

void IBParticles::createParticle2D(const dictionary& dict)
{
    Info<< setw(4) <<"  2D particle with parameters:"<<endl;
	center_ = dict.lookup("CG");
	R_ = readScalar(dict.lookup("radius"));
	nPoints_ = readScalar(dict.lookup("nPoints"));
	rho_ = readScalar(dict.lookup("rho"));
	V_ = PI*R_*R_;
    Ip_.x() = 0.0;
    Ip_.y() = 0.0;
	Ip_.z() = rho_*V_*R_*R_/2.0;

    lPoints_.setSize(nPoints_, vector::zero);

    scalar dAlpha = 2.0*PI/nPoints_;
    forAll(lPoints_, i)
    {
        lPoints_[i] = point (   center_.x() + R_*Foam::cos(dAlpha*i),
                                center_.y() + R_*Foam::sin(dAlpha*i),
                                center_.z()
                            );
    }
    nFaces_ = nPoints_;
    nPointsOfFaces_.setSize(nFaces_,3);
    pointOfFace_.setSize(nFaces_);
    for(int i=0;i<nFaces_;i++)
    {
        pointOfFace_[i].setSize(3);
        pointOfFace_[i][0] = i;
        pointOfFace_[i][1] = (i+1)%nPoints_;
        pointOfFace_[i][2] = nPoints_;
    }
}

void IBParticles::createParticle3D(const dictionary& dict)
{
    Info<< setw(4) <<"  3D particle with parameters:"<<endl;
	center_ = dict.lookup("CG");
	R_ = readScalar(dict.lookup("radius"));
    nPoints_ = readLabel(dict.lookup("nPoints"));
	rho_ = readScalar(dict.lookup("rho"));
	V_ = 4.0/3.0*PI*R_*R_*R_;
    Ip_.x() = 2.0/5.0*rho_*V_*R_*R_;    
    Ip_.y() = 2.0/5.0*rho_*V_*R_*R_;    
    Ip_.z() = 2.0/5.0*rho_*V_*R_*R_;	

    //- Distribute nPoints points on surface of unit sphere
    lPoints_ = createUnitSphereEqAreaPartition(nPoints_);
    
    //- Translate and scale unit sphere according to center_ and R_
    //  to obtain desired 3D sphere
    for (int pointI=0;pointI<nPoints_;pointI++)
    {
        lPoints_[pointI].x()=R_*lPoints_[pointI].x()+center_.x();
        lPoints_[pointI].y()=R_*lPoints_[pointI].y()+center_.y();
        lPoints_[pointI].z()=R_*lPoints_[pointI].z()+center_.z();
    }

}

pointField IBParticles::createUnitSphereEqAreaPartition(label nPoints)
{
    //- Paul Leopardi algorithm 2007
    //  Evenly distribute nPoints points on a unit sphere
    //- 1. Determine colatitudes of polar caps

        //- Area of a region
        scalar Sr = 4.0*PI/nPoints;
        
        //- Collatitude of the North collar cap.
        //  Calculate form area of the cap: Sr = 4*PI*[sin(phiNC/2)]^2
        scalar phiNC = 2.0*asin(sqrt(Sr/(4.0*PI)));
    
    //- 2. Determine an ideal collar angle
    //  (when number of points reach infinity)
       
        scalar deltaI = pow(Sr, 0.5);
    //- 3. Determine an ideal number of collars

        scalar nI = (PI-2.0*phiNC) / deltaI;

    //- 4. Determine the actual number of collars
       
        label n = static_cast<int>(nI+0.5);

    //- 5. Create a list of the ideal number of regions in each collar
    
        //- Fitting collar angle
        scalar deltaF = (PI - 2.0*phiNC)/n;

        //- List of colatitudes of caps
        List<scalar> colatitudes(n+2, 0.0);
        
        forAll(colatitudes, i)
        {
            if (i == 0)
                colatitudes[i] = phiNC;
            else if( i != 0 && i < n+1 )
                colatitudes[i] = phiNC + i*deltaF;
            else
                colatitudes[i] = PI;    
        }

        //- Ideal number of region in each collar
        List<scalar> idealRegionsPerCollar(n+2, 0.0);
        
        forAll(idealRegionsPerCollar, i)
        {
            if (i==0 || i==n+1)
            {
                idealRegionsPerCollar[i] = 1;
            }
            else 
            {
                idealRegionsPerCollar[i] 
                = 
                    (capArea(colatitudes[i]) - capArea(colatitudes[i-1]))/Sr;
            }
        }

    //- 6. Createa a list of actial number of regions in each collar
        
        List<label> actualRegionsPerCollar(n+2,0);
        
        forAll(actualRegionsPerCollar, i)
        {
            if (i==0 || i==n+1)
            {
                actualRegionsPerCollar[i] = 1;
            }
            else if (i == 1)
            {
                actualRegionsPerCollar[i] = static_cast<int>(idealRegionsPerCollar[i]);
            }
            else 
            {
                scalar ai(0.0);
                for(int j=0; j<i; j++)
                {
                    ai += (idealRegionsPerCollar[j] - actualRegionsPerCollar[j]) ;
                } 

                actualRegionsPerCollar[i] = static_cast<int>(idealRegionsPerCollar[i]+ai);
            }
        }

    //- 7. Create a list of colatitudes of points on each collar
        List<List<scalar>> pointLatitudes(n+2);
        List<List<scalar>> pointLongtitudes(n+2);

        for(int i=0; i<n+2; i++)
        {
            if (i==0)
            {
                pointLatitudes[i].setSize(1);
                pointLongtitudes[i].setSize(1);
                forAll(pointLatitudes[i], j)
                {
                    pointLatitudes[i][j] = 0;
                    pointLongtitudes[i][j] = 0;
                }
            }
            else if ( i!=0 && i<n+1)
            {
                pointLatitudes[i].setSize(actualRegionsPerCollar[i]);
                pointLongtitudes[i].setSize(actualRegionsPerCollar[i]);

                forAll(pointLatitudes[i], j)
                {
                    pointLatitudes[i][j] = colatitudes[i-1] + 0.5*(colatitudes[i] - colatitudes[i-1]);
                    pointLongtitudes[i][j] = j*2.0*PI/actualRegionsPerCollar[i];
                }

            }
            else 
            {
                pointLatitudes[i].setSize(1);
                pointLongtitudes[i].setSize(1);
                forAll(pointLatitudes[i], j)
                {
                    pointLatitudes[i][j] = PI;
                    pointLongtitudes[i][j] = 0;
                }
            }
        }

    //8. Create list of points
        List<List<point> > points(n+2);
        pointField allPoints;

        for(int i=0; i<n+2; i++)
        {   
            points[i].setSize(actualRegionsPerCollar[i]);
            forAll(pointLatitudes[i], j)
            {
                point p;
                p.x() = sin(pointLatitudes[i][j])*cos(pointLongtitudes[i][j]);
                p.y() = sin(pointLatitudes[i][j])*sin(pointLongtitudes[i][j]); 
                p.z() = cos(pointLatitudes[i][j]);
                points[i][j] = p;
                allPoints.append(p);
            }
        }
    
    //9. Write VTU
        // nFaces_ = 2*nPoints_ - 4;
        // nPointsOfFaces_.setSize(nFaces_, 3);
        // pointOfFace_.setSize(nFaces_);
        // for (int i=0; i<n+2; i++)
        // {
        //     labelList pofi(3);
        //     forAll(pointLatitudes[i], j)
        //     {
        //         if (i==0)
        //         {
        //             continue;
        //         }
        //         if (i==1)
        //         {
        //             pofi[0] = points[i][j];
        //             pofi[1] = points[i][(j+1)%pointLatitudes[i].size()];
        //             pofi[2] = points[i-1];
        //         }
        //     }
        // }
    return allPoints;
}

scalar IBParticles::capArea(scalar angle)
{
    return (4.0*PI*pow(sin(angle/2.0), 2));
}

void IBParticles::addMotions(const dictionary& dict)
{
    if(movable_)
    {
        const dictionary& motionDict = dict.subDict("motions");
        
        label i=0;
        motions_.setSize(motionDict.size());
        forAllConstIter(IDLList<entry>, motionDict, iter)
        {
            if(iter().isDict())
            {
                motions_.set
                (
                    i++,
                    IBMotion::New 
                    (
                        iter().keyword(),
                        *this,
                        iter().dict()
                    )
                );
                Info<<"  - Motion"<<i<<": "<<iter().keyword()<<endl;
                Info<< setw(15) <<"U = "<<motions_[i-1].UTranslate()<<endl;
                Info<< setw(15) <<"w = "<<motions_[i-1].URotate()<<endl;
            }
        }
        motions_.setSize(i);
    }
    else
    {
        Info<< setw(11) <<"  - Motions: N/A"<<nl<<endl;
    }
}
//--------------------------------Constructors-------------------------------//

IBParticles::IBParticles
(	
	const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
	IBObject(name, mesh, dict),
	PI(3.14159265359),
    mesh_(mesh),
    name_(name),
    objectType_(dict.lookup("type")),
    center_(point::zero),
    R_(0.0),
    rho_(0.0),
    V_(0.0),
    Ip_(vector::zero),
    nPoints_(0),   
	lPoints_(),
    movable_(dict.lookupOrDefault<Switch>("movable", false)),
    motions_(),
	nFaces_(0),
	nPointsOfFaces_(),
	pointOfFace_()
{
    if(mesh.nGeometricD() == 2)
        createParticle2D(dict);
    else
        createParticle3D(dict);

    Info <<"  - Center : "<< center_ <<endl;
    Info <<"  - Radius : "<< R_ <<endl;
    Info <<"  - nPoints: "<< nPoints_ <<endl;
    Info <<"  - Density: "<< rho_ <<endl;

    addMotions(dict);
}

//--------------------------------Member Functions---------------------------//
word IBParticles::name()
{
    return name_;
}

word IBParticles::objectType() 
{
    return objectType_;
}

PtrList<IBMotion>& IBParticles::motions() 
{
    return motions_;
}

scalar IBParticles::rho() 
{
    return rho_;
}

scalar IBParticles::V() 
{
    return V_;
}

vector IBParticles::Ip() 
{
    return Ip_;
}

scalar IBParticles::R() 
{
    return R_;
}

point IBParticles::CG() 
{
    return center_;
}

bool IBParticles::movable() 
{
    return movable_;
}

label IBParticles::nPoints() 
{
    return nPoints_;
}

pointField IBParticles::lPoints() 
{
    return lPoints_;
}

vector IBParticles::translationalVelocity()
{
    vector UT = vector::zero;
    if (movable_)
    {
        forAll(motions_, i)
        {
            UT += motions_[i].UTranslate();
        }
    }

    return UT;
}

vector IBParticles::rotationalVelocity()
{
    vector UR = vector::zero;
    if (movable_)
    {
        forAll(motions_, i)
        {
            UR += motions_[i].URotate();
        }
    }

    return UR;
}

label IBParticles::nFaces()
{
    return nFaces_;
}

labelList IBParticles::nPointsOfFaces()
{
    return nPointsOfFaces_;
}

labelListList IBParticles::pointOfFace()
{
    return pointOfFace_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //