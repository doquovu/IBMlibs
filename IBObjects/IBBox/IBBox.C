/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBBox class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "IBBox.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"

// * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IBBox, 0);
    addToRunTimeSelectionTable(IBObject, IBBox, dictionary);


// * * * * * * * * * * * * Private Member Fuctions * * * * * * * * * * * * * //

void IBBox::readTXTFile(word fileName, point center, scalar R)
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

void IBBox::createBox2D(const dictionary& dict)
{
    Info<< setw(4) <<"  2D Box with parameters:"<<endl;
    max_ = dict.lookup("max");
	min_ = dict.lookup("min");
    scalar dx = max_.x() - min_.x();
    scalar dy = max_.y() - min_.y();
    scalar dz = max_.z() - min_.z();
	R_ = dx;
	nPoints_ = readScalar(dict.lookup("nPoints"));
	rho_ = readScalar(dict.lookup("rho"));
	V_ = dx*dy;
	Ipx_ = rho_*V_(dx*dx + dy*dy)/12.0;
    Ipy_ = Ipx_;
    Ipz_ = Ipx_;
}

void IBBox::createBox3D(const dictionary& dict)
{
    Info<< setw(4) <<"  3D Box with parameters:"<<endl;
	max_ = dict.lookup("max");
    min_ = dict.lookup("min");
    scalar dx = max_.x() - min_.x();
    scalar dy = max_.y() - min_.y();
    scalar dz = max_.z() - min_.z();
    R_ = dx;
    nPoints_ = readScalar(dict.lookup("nPoints"));
    rho_ = readScalar(dict.lookup("rho"));
    V_ = dx*dy*dz;
    Ipx_ = rho_*V_(dy*dy + dz*dz)/12.0;
    Ipy_ = rho_*V_(dz*dz + dx*dx)/12.0;
    Ipz_ = rho_*V_(dx*dx + dy*dy)/12.0;
}

void IBBox::addMotions(const dictionary& dict)
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

IBBox::IBBox
(	
	const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
	IBObject(name, mesh, dict),
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
    Fk_(),
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
word IBBox::name()
{
    return name_;
}

word IBBox::objectType() 
{
    return objectType_;
}

PtrList<IBMotion>& IBBox::motions() 
{
    return motions_;
}

scalar IBBox::rho() 
{
    return rho_;
}

scalar IBBox::V() 
{
    return V_;
}

scalar IBBox::Ip() 
{
    return Ip_;
}

scalar IBBox::R() 
{
    return R_;
}

point IBBox::CG() 
{
    return center_;
}

bool IBBox::movable() 
{
    return movable_;
}

label IBBox::nPoints() 
{
    return nPoints_;
}

pointField IBBox::lPoints() 
{
    return lPoints_;
}

vector IBBox::translationalVelocity()
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

vector IBBox::rotationalVelocity()
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

label IBBox::nFaces()
{
    return nFaces_;
}

labelList IBBox::nPointsOfFaces()
{
    return nPointsOfFaces_;
}

labelListList IBBox::pointOfFace()
{
    return pointOfFace_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //