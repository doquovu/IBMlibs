/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBSTL class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "IBSTL.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"

// * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IBSTL, 0);
    addToRunTimeSelectionTable(IBObject, IBSTL, dictionary);

// * * * * * * * * * * * * Private Member Fuctions * * * * * * * * * * * * * //
void IBSTL::readSTL(const dictionary& dict)
{
    word fileName = word(dict.lookup("fileName"));
    
    Info<<"  Reading "<<fileName<<nl<<endl;
    
    triSurfaceMesh triSurf_ 
    (
        IOobject
        (
            fileName,
            mesh_.time().constant(),
            "triSurface",
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    lPoints_ = triSurf_.localPoints();
    nPoints_ = lPoints_.size();
    nFaces_  = triSurf_.faceCentres().size();
    nPointsOfFaces_.setSize(nFaces_, 3); 
    pointOfFace_.setSize(nFaces_);

    forAll(nPointsOfFaces_, fi)
    {
        pointOfFace_[fi].setSize(nPointsOfFaces_[fi]);

        forAll(pointOfFace_[fi], pi)
        {
            pointOfFace_[fi][pi] = triSurf_.localFaces()[fi][pi];
        }
    }
    Ip_ = dict.lookup("Ip");   
}

void IBSTL::addMotions(const dictionary& dict)
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

IBSTL::IBSTL
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
    readSTL(dict);
    addMotions(dict);
}

//--------------------------------Member Functions---------------------------//
word IBSTL::name()
{
    return name_;
}

word IBSTL::objectType() 
{
    return objectType_;
}

PtrList<IBMotion>& IBSTL::motions() 
{
    return motions_;
}

scalar IBSTL::rho() 
{
    return rho_;
}

scalar IBSTL::V() 
{
    return V_;
}

vector IBSTL::Ip() 
{
    return Ip_;
}

point IBSTL::CG() 
{
    return center_;
}

scalar IBSTL::R()
{
    return R_;
}

bool IBSTL::movable() 
{
    return movable_;
}

label IBSTL::nPoints() 
{
    return nPoints_;
}

pointField IBSTL::lPoints() 
{
    return lPoints_;
}

vector IBSTL::translationalVelocity()
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

vector IBSTL::rotationalVelocity()
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

label IBSTL::nFaces()
{
    return nFaces_;
}

labelList IBSTL::nPointsOfFaces()
{
    return nPointsOfFaces_;
}

labelListList IBSTL::pointOfFace()
{
    return pointOfFace_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //