/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBObjectRegistry class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "IBObjectRegistry.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IBObjectRegistry, 0);
}
//-----------------------------Protected Functions---------------------------//
Foam::scalar Foam::IBObjectRegistry::cellSize(label cellID)
{
    scalar delta = 0.0;
    if (mesh_.nGeometricD() == 3)
    {
        delta = pow(mesh_.V().field()[cellID], 1.0/3.0);
    }
    else
    {
        scalar thickness = 0.0;
        Vector<label> directions = mesh_.geometricD();
        for (direction dir = 0; dir < directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = mesh_.bounds().span()[dir];
                break;
            }
        }

        delta = Foam::sqrt(mesh_.V().field()[cellID]/thickness);
    }

    return delta;
}

Foam::labelListList Foam::IBObjectRegistry::findNeiCells
(
    const pointField& points
)
{
    labelListList neighbour(points.size());

    forAll(points, pointI)
    {
        forAll(mesh_.C(), cellI)
        { 
            scalar cellToLpoint = mag(mesh_.C()[cellI] - points[pointI]);
            
            scalar span = 2.0*h_;
            
            if (cellToLpoint <= span)
            {
                neighbour[pointI].append(cellI);
            }
        }
    }

    return neighbour;
}

Foam::labelList Foam::IBObjectRegistry::findSolidCells
(
    label objectID,
    point C
)
{
    labelList slc;
    scalar span = R_[objectID]+h_/sqrt(2.0);
    forAll(mesh_.C(), cellI)
    { 
        scalar dR = mag(mesh_.C()[cellI] - C);
        //- not very effective!
        
        if (dR <= span)
        {
            slc.append(cellI);
        }
    }
    return slc;
}

Foam::labelList Foam::IBObjectRegistry::findSolidCellsInt
(
    label objectID,
    point C
)
{
    labelList slcI;
    scalar R = R_[objectID];

    forAll(mesh_.C(), cellI)
    {
        scalar dR = mag(mesh_.C()[cellI] - C);
        if (dR <= R)
        {
            slcI.append(cellI);
        }
    }
    return slcI;
}

bool Foam::IBObjectRegistry::isInsideRegion1(const label objectID, const point p)
{
    const boundBox& box = mesh_.bounds();
    const point& min = box.min();
    const point& max = box.max();
    scalar R = R_[objectID] + 2.0*h_;
    if(    mesh_.nGeometricD() == 2 
        && p.x() > (min.x()+R) && p.x() < (max.x()-R) 
        && p.y() > (min.y()+R) && p.y() < (max.y()-R)
        && p.z() > min.z()     && p.z() < max.z() )
        return true;
    else if (    mesh_.nGeometricD() == 3 
              && p.x() > (min.x()+R) && p.x() < (max.x()-R) 
              && p.y() > (min.y()+R) && p.y() < (max.y()-R)
              && p.z() > (min.z()+R) && p.z() < (max.z()-R) )
        return true;
    else if (    p.x() == (min.x()+R) || p.x() == (max.x()-R) 
              || p.y() == (min.y()+R) || p.y() == (max.y()-R)
              || p.z() ==  min.z()    || p.z() ==  max.z()   
              || p.z() == (min.z()+R) || p.z() == (max.z()-R) )  
        return true;
    else
        return false;
} 

bool Foam::IBObjectRegistry::isInsideRegion2(const label objectID, const point p)
{
    const boundBox& box = mesh_.bounds();
    const point& min = box.min();
    const point& max = box.max();
    if (   p.x() > (min.x()-SMALL) && p.x() < (max.x()+SMALL) 
        && p.y() > (min.y()-SMALL) && p.y() < (max.y()+SMALL)
        && p.z() > (min.z()-SMALL) && p.z() < (max.z()+SMALL) )
        return true;
    else if (   p.x() == min.x() || p.x() == max.x() 
             || p.y() == min.y() || p.y() == max.y()
             || p.z() == min.z() || p.z() == max.z() )
        return true;
    else
        return false;
} 

bool Foam::IBObjectRegistry::isInsideRegion3(const label objectID, const point p)
{
    const boundBox& box = mesh_.bounds();
    const point& min = box.min();
    const point& max = box.max();
    scalar R = R_[objectID] + 2.0*h_ +SMALL;
    if (   mesh_.nGeometricD() == 2 
        && p.x() > (min.x()-R) && p.x() < (max.x()+R) 
        && p.y() > (min.y()-R) && p.y() < (max.y()+R)
        && p.z() >  min.z()    && p.z() < max.z() )
        return true;
    else if(   mesh_.nGeometricD() == 3 
            && p.x() > (min.x()-R) && p.x() < (max.x()+R) 
            && p.y() > (min.y()-R) && p.y() < (max.y()+R)
            && p.z() > (min.z()-R) && p.z() < (max.z()+R) )
        return true;
    else if (    p.x() == (min.x()+R) || p.x() == (max.x()+R) 
              || p.y() == (min.y()-R) || p.y() == (max.y()+R)
              || p.z() ==  min.z()    || p.z() ==  max.z()   
              || p.z() == (min.z()-R) || p.z() == (max.z()+R) ) 
        return true;
    else
        return false;
} 

void Foam::IBObjectRegistry::cartesianGridSize()
{
	scalarField cellsSize(mesh_.C().size(), 0.0);

    forAll(mesh_.C(), cellI)
    {
        cellsSize[cellI] = cellSize(cellI);
    }
    
    h_ = min(cellsSize);

    if (mesh_.nGeometricD() == 2)
    {
    	dV_ = h_*h_;
    }
    else 
    {
    	dV_ = h_*h_*h_;
    }
}

void Foam::IBObjectRegistry::detectCyclicBoundary()
{
    cyclicDistance_ = vector::zero;
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];
        if (isA<cyclicPolyPatch>(pp))
        {
            periodicBC_ = true;
            
            Info<< "\nFOUND CYCLIC BOUNDARY!"<<endl;
            
            const cyclicPolyPatch& cpp1 = refCast<const cyclicPolyPatch>
            (
                mesh_.boundaryMesh()[patchI]
            );
            //if (cpp1.parallel())
            {
                const boundBox& box = mesh_.bounds();
                const vector patchNorm = cpp1.faceNormals()[0];
                const vector maxToMin = box.max() - box.min();
                cyclicDistance_ = (maxToMin & patchNorm)*patchNorm;
            }
            Info<< "  cyc1 to cyc2 = "<<cyclicDistance_<<nl<<endl;
            return;
        }
    }
}

void Foam::IBObjectRegistry::createObjects
(
    const dictionary& dict
)
{
    if (dict.found("IBObjects"))
    {
        const dictionary& IBObjectDict = dict.subDict("IBObjects");

        label i=0;
        nObjects_ = IBObjectDict.size();
        objects_.setSize(nObjects_);
        forAllConstIter(IDLList<entry>, IBObjectDict, iter)
        {
            if (iter().isDict())
            {
                Info<<"Creating Immersed Object "<<i+1
                    <<" named "<<iter().keyword()<<":"<<endl;
                objects_.set
                (
                    i++,
                    IBObject::New 
                    (
                        iter().keyword(),
                        mesh_,
                        iter().dict()
                    )
                );
            }
        }
        objects_.setSize(i);
    }
}

void Foam::IBObjectRegistry::createShadows()
{
    if (!periodicBC_)
    {
        Shd1CoG_.setSize(0);
        Shd2CoG_.setSize(0);
        Shd1LPoints_.setSize(0);
        Shd2LPoints_.setSize(0);
        Shd1NeiCells_.setSize(0);
        Shd2NeiCells_.setSize(0);
        Shd1SolidCells_.setSize(0);
        Shd2SolidCells_.setSize(0);
        enableShadows_.setSize(nObjects_);
        forAll(enableShadows_, i)
            enableShadows_[i] = false;
        return;
    }
    else
    {
        enableShadows_.setSize(nObjects_);
        Shd1CoG_.setSize(nObjects_);
        Shd2CoG_.setSize(nObjects_);
        Shd1LPoints_.setSize(nObjects_);
        Shd2LPoints_.setSize(nObjects_);
        Shd1NeiCells_.setSize(nObjects_);
        Shd2NeiCells_.setSize(nObjects_);
        Shd1SolidCells_.setSize(nObjects_);
        Shd2SolidCells_.setSize(nObjects_);
        
        for(int i=0; i<nObjects_; i++)
        {
            if (!objects_[i].movable())
            {
                enableShadows_[i] = false;
                Shd1LPoints_[i].setSize(0);
                Shd2LPoints_[i].setSize(0);
                Shd1NeiCells_[i].setSize(0);
                Shd2NeiCells_[i].setSize(0);
                Shd1SolidCells_[i].setSize(0);
                Shd2SolidCells_[i].setSize(0);
            }
            else
            {
                enableShadows_[i] = true;
                Shd1CoG_[i] = CoG_[i] - cyclicDistance_;
                Shd2CoG_[i] = CoG_[i] + cyclicDistance_;

                Shd1LPoints_[i].setSize(LPoints_[i].size());
                Shd2LPoints_[i].setSize(LPoints_[i].size());
                forAll(LPoints_[i], pointI)
                {
                    Shd1LPoints_[i][pointI] = LPoints_[i][pointI] - cyclicDistance_;
                    Shd2LPoints_[i][pointI] = LPoints_[i][pointI] + cyclicDistance_;
                }

                if (isInsideRegion1(i, CoG_[i]))
                {
                    Shd1NeiCells_[i].setSize(0);
                    Shd2NeiCells_[i].setSize(0);
                    Shd1SolidCells_[i].setSize(0);
                    Shd2SolidCells_[i].setSize(0);
                }
                else if (isInsideRegion3(i, Shd1CoG_[i]))
                {
                    Shd1NeiCells_[i] = findNeiCells(Shd1LPoints_[i]);
                    Shd2NeiCells_[i].setSize(0);
                    Shd1SolidCells_[i] = findSolidCells(i, Shd1CoG_[i]);
                    Shd2SolidCells_[i].setSize(0);
                }
                else 
                {
                    Shd1NeiCells_[i].setSize(0);
                    Shd2NeiCells_[i] = findNeiCells(Shd2LPoints_[i]);
                    Shd1SolidCells_[i].setSize(0);
                    Shd2SolidCells_[i] = findSolidCells(i, Shd2CoG_[i]);
                }
            }
        }
    }
}

void Foam::IBObjectRegistry::swap(point& p1, point& p2)
{
    point dummyP = p1;

    p1 = p2;

    p2 = dummyP;
}

void Foam::IBObjectRegistry::getObjectsData()
{
    R_.setSize(nObjects_);
    CoG_.setSize(nObjects_);
    LPoints_.setSize(nObjects_);
    neiCells_.setSize(nObjects_);
    solidCells_.setSize(nObjects_);
    solidCellsInt_.setSize(nObjects_);
    UTranslate_.setSize(nObjects_);
    URotate_.setSize(nObjects_);

    nFaces_.setSize(nObjects_);
    nPointsOfFaces_.setSize(nObjects_);
    pointOfFace_.setSize(nObjects_);

    forAll(objects_, i)
    {
        R_[i]           = objects_[i].R();
        CoG_[i]         = objects_[i].CG();
        nFaces_[i]      = objects_[i].nFaces();
        LPoints_[i]     = objects_[i].lPoints();
        neiCells_[i]    = findNeiCells(LPoints_[i]);
        solidCells_[i]  = findSolidCells(i,CoG_[i]);
        solidCellsInt_[i]  = findSolidCellsInt(i,CoG_[i]);
        URotate_[i]         = objects_[i].rotationalVelocity();
        UTranslate_[i]      = objects_[i].translationalVelocity();
        nPointsOfFaces_[i]  = objects_[i].nPointsOfFaces();
        pointOfFace_[i]     = objects_[i].pointOfFace();
    }
}

//---------------------------------Constructors------------------------------//
Foam::IBObjectRegistry::IBObjectRegistry
(
	const fvMesh& mesh,
	const dictionary& dict
)
:
	mesh_(mesh),
    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
	h_(),
	dV_(),
    rhoF_(readScalar(transportProperties_.lookup("rhoF"))),
    g_(dimensionedVector(transportProperties_.lookup("g")).value()),
	nObjects_(),
	objects_(), 
    LPoints_(),
    R_(),
    CoG_(),
	neiCells_(),
    UTranslate_(),
    URotate_(),
    nFaces_(),
    nPointsOfFaces_(),
    pointOfFace_(),
    periodicBC_(false),
    enableShadows_(),
    cyclicDistance_(),
    Shd1LPoints_(),
    Shd2LPoints_(),
    Shd1CoG_(),
    Shd2CoG_(),
    Shd1NeiCells_(),
    Shd2NeiCells_()
{
	cartesianGridSize();
	createObjects(dict);
	getObjectsData();
    detectCyclicBoundary();
    createShadows();
}

//-------------------------------Member functions----------------------------//

Foam::scalar Foam::IBObjectRegistry::volFraction
(
    label objID, 
    label cellID, 
    point c
)
{
    const faceList& ff = mesh_.faces();
    const pointField& pp = mesh_.points();
    const cell& cc = mesh_.cells()[cellID];
    pointField cellVertices = cc.points(ff, pp);

    scalar alphaIJK(0.0);
    scalar sumPhi(0.0);
    forAll(cellVertices, pointI)
    {
        scalar phi_m = LevelSetFunc(objID, cellVertices[pointI], c);
        alphaIJK += -phi_m * HeavisideFunc(-phi_m);
        sumPhi += mag(phi_m);
    }

    return (alphaIJK/sumPhi);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
