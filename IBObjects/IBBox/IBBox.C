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
void IBBox::createBox2D(const dictionary& dict)
{
    Info<< setw(4) <<"  2D Box with parameters:"<<endl;
    max_ = dict.lookup("max");
	min_ = dict.lookup("min");
    scalar lx = max_.x() - min_.x();
    scalar ly = max_.y() - min_.y();
    center_ = point( min_.x() + 0.5*lx, 
                     min_.y() + 0.5*ly,
                     min_.z());
	R_ = 0.5*sqrt(lx*lx + ly*ly);
    rho_ = readScalar(dict.lookup("rho"));
    V_ = lx*ly;
    Ip_.x() = 0.0;
    Ip_.y() = 0.0;
    Ip_.z() = rho_*V_*(lx*lx + ly*ly)/12.0;

    //- NOTE: number of interval, not number of points
    label nx = readScalar(dict.lookup("nx"));  
    label ny = readScalar(dict.lookup("ny")); 
    nPoints_ = 2 * (nx + ny);
    lPoints_.setSize(nPoints_);

    scalar dx = lx/nx;
    scalar dy = ly/ny;
    forAll(lPoints_, i)
    {
        if (i == 0)
            lPoints_[i] = min_;
        if ( (i!=0) and (i<=nx) )
        {
            lPoints_[i].x() = lPoints_[i-1].x() + dx; 
            lPoints_[i].y() = lPoints_[i-1].y(); 
            lPoints_[i].z() = lPoints_[i-1].z();
        }
        if ( (i>nx) and (i<=nx+ny) )
        {
            lPoints_[i].x() = lPoints_[i-1].x(); 
            lPoints_[i].y() = lPoints_[i-1].y() + dy; 
            lPoints_[i].z() = lPoints_[i-1].z();
        }
        if ( (i>nx+ny) and (i<=nx+ny+nx) )
        {
            lPoints_[i].x() = lPoints_[i-1].x() - dx; 
            lPoints_[i].y() = lPoints_[i-1].y(); 
            lPoints_[i].z() = lPoints_[i-1].z();
        }
        if ( (i>nx+ny+nx) and (i<nPoints_))
        {
            lPoints_[i].x() = lPoints_[i-1].x(); 
            lPoints_[i].y() = lPoints_[i-1].y() - dy; 
            lPoints_[i].z() = lPoints_[i-1].z();
        }
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

void IBBox::createBox3D(const dictionary& dict)
{
    Info<< setw(4) <<"  3D Box with parameters:"<<endl;
	max_ = dict.lookup("max");
    min_ = dict.lookup("min");
    scalar lx = max_.x() - min_.x();
    scalar ly = max_.y() - min_.y();
    scalar lz = max_.z() - min_.z();
    center_ = point( min_.x() + 0.5*lx, 
                     min_.y() + 0.5*ly,
                     min_.z() + 0.5*lz);
    R_ = 0.5*sqrt(lx*lx + ly*ly + lz*lz);
    rho_ = readScalar(dict.lookup("rho"));
    V_ = lx*ly*lz;
    Ip_.x() = rho_*V_*(ly*ly + lz*lz)/12.0;
    Ip_.y() = rho_*V_*(lz*lz + lx*lx)/12.0;
    Ip_.z() = rho_*V_*(lx*lx + ly*ly)/12.0;

    //- NOTE: number of interval, not number of points
    label nx = readScalar(dict.lookup("nx"));  
    label ny = readScalar(dict.lookup("ny")); 
    label nz = readScalar(dict.lookup("nz")); 
    nPoints_ = 2*(nx + ny )*(nz + 1);
    lPoints_.setSize(nPoints_);

    // scalar dx = lx/nx;
    // scalar dy = ly/ny;
    // scalar dz = lz/nz;

    // forAll(lPoints_, i)
    // {
    //     scalar zInstant(0.0);
    //     label nxy = 2*(nx + ny);
    //     if ((i%nxy)==0) 
    //         zInstant = min_.z() + i*dz/nxy;
    // }
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
    nPoints_(0),   
	lPoints_(),
    movable_(dict.lookupOrDefault<Switch>("movable", false)),
    motions_(),
	nFaces_(0),
	nPointsOfFaces_(),
	pointOfFace_()
{
    if(mesh.nGeometricD() == 2)
        createBox2D(dict);
    else
        createBox3D(dict);

    Info <<"  - Center    : "<< center_ <<endl;
    Info <<"  - Hypotenuse: "<< 2.0*R_ <<endl;
    Info <<"  - nPoints   : "<< nPoints_ <<endl;
    Info <<"  - Density   : "<< rho_ <<endl;

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

vector IBBox::Ip() 
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