/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBWall class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "IBWall.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"

// * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IBWall, 0);
    addToRunTimeSelectionTable(IBObject, IBWall, dictionary);


// * * * * * * * * * * * * Private Member Fuctions * * * * * * * * * * * * * //
void IBWall::createWall2D(const dictionary& dict)
{
    Info<< setw(4) <<"  2D Wall with parameters:"<<endl;
    max_ = dict.lookup("max");
	min_ = dict.lookup("min");
    scalar length = mag(max_ - min_);

    label nInterval = readLabel(dict.lookup("nInterval"));
    nPoints_ = nInterval + 1;
    lPoints_.setSize(nPoints_);

    scalar delta = length/nInterval;
    scalar sinAlpha = (max_.y() - min_.y())/length;
    scalar cosAlpha = (max_.x() - min_.x())/length;
    
    forAll(lPoints_, i)
    {
        if (i == 0)
        {
            lPoints_[i] = min_;
        }
        else if (i == lPoints_.size()-1)
        {
            lPoints_[i] = max_;
        }
        else 
        {
            lPoints_[i] 
            = 
                point 
                (
                    lPoints_[i-1].x() + delta*cosAlpha,
                    lPoints_[i-1].y() + delta*sinAlpha,
                    lPoints_[i-1].z()
                );
        }
    }
}

//--------------------------------Constructors-------------------------------//

IBWall::IBWall
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
    nPoints_(),   
	lPoints_(),
    motions_()
{
    createWall2D(dict);

    Info <<"  - nPoints   : "<< nPoints_ <<endl;
    Info <<"  - Max bound : "<<max_<<endl;
    Info <<"  - Min bound : "<<max_<<endl;
}

//--------------------------------Member Functions---------------------------//
word IBWall::name()
{
    return name_;
}

word IBWall::objectType() 
{
    return objectType_;
}

PtrList<IBMotion>& IBWall::motions() 
{
    return motions_;
}

scalar IBWall::rho() 
{
    return 0.0;
}

scalar IBWall::V() 
{
    return 0.0;
}

vector IBWall::Ip() 
{
    return vector::zero;
}

scalar IBWall::R() 
{
    return 0.0;
}

point IBWall::CG() 
{
    return (min_ + 0.5*(max_-min_));
}

bool IBWall::movable() 
{
    return false;
}

label IBWall::nPoints() 
{
    return nPoints_;
}

pointField IBWall::lPoints() 
{
    return lPoints_;
}

vector IBWall::translationalVelocity()
{
    return vector::zero;
}

vector IBWall::rotationalVelocity()
{
    return vector::zero;
}

label IBWall::nFaces()
{
    return 0;
}

labelList IBWall::nPointsOfFaces()
{
    labelList nullList(0);
    return nullList;
}

labelListList IBWall::pointOfFace()
{
    labelListList nullList(0);
    return nullList;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //