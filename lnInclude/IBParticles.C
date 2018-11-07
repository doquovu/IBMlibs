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
	Ip_ = rho_*V_*R_*R_/2.0;

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
    nPointsOfFaces_.setSize(nFaces_);
    pointOfFace_.setSize(nFaces_);
    for(int i=0;i<nFaces_;i++)
    {
        nPointsOfFaces_[i] = 3;
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
	rho_ = readScalar(dict.lookup("rho"));
	V_ = 4.0/3.0*PI*R_*R_*R_;
	Ip_ = 2.0/5.0*rho_*V_*R_*R_;	

	word filename = dict.lookup("fileName");
	readTXTFile(filename, center_, R_);
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
	PI(3.141593),
    mesh_(mesh),
    name_(name),
    objectType_(dict.lookup("type")),
    center_(point::zero),
    R_(0.0),
    rho_(0.0),
    V_(0.0),
    Ip_(0.0),
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

scalar IBParticles::Ip() 
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