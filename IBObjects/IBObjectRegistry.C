/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBObjectRegistry class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "IBObjectRegistry.H"
#include "labelList.H"
#include "IBBox.H"
#include "vectorTools.H"
#include "IBSTL.H"
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

Foam::vector Foam::IBObjectRegistry::rotate(vector V, point root, scalar angle)
{
    point tip = point(root.x()+V.x(), root.y()+V.y(), root.z()+V.z() );
    point root_Vrot(0,0,0);
    
    if(root.x() == root.y() && root.x() == 0.0)
        root_Vrot = root;
    else
        root_Vrot = rotate(root, angle);

    point tip_Vrot = rotate(tip, angle);

    return (tip_Vrot - root_Vrot);
}
Foam::point Foam::IBObjectRegistry::rotate(point p, scalar angle)
{
    scalar p_angle = getAngle(p);
    scalar rot_angle = p_angle + angle;

    //- Assume spiral channel center is always (0 0 0)
    scalar R = mag(p);
    point P_rot = point (R*Foam::cos(rot_angle), R*Foam::sin(rot_angle), p.z());

    return P_rot;
}
Foam::scalar Foam::IBObjectRegistry::getAngle(point p)
{
    vector pProjectToXOY = vector(p.x(), p.y(), 0);
    vector R = pProjectToXOY/(mag(pProjectToXOY));
    vector i(1, 0, 0);
    scalar alpha = Foam::vectorTools::radAngleBetween(R, i);
    if (p.y() >= 0)
        return alpha;
    else 
        return (-alpha);
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
    Info<<"IBM: Finding solid cells for object "<<objects_[objectID].name()<<endl;
    labelList slc;
    if (objects_[objectID].objectType() == "IBParticles")
    {
        // scalar span = R_[objectID]+h_/sqrt(2.0);
        scalar span = R_[objectID] ;//+ 1.5*h_;
        forAll(mesh_.C(), cellI)
        { 
            scalar dR = mag(mesh_.C()[cellI] - C);
            //- not very effective!
            
            if (dR <= span)
            {
                slc.append(cellI);
            }
        }
    }
    if (objects_[objectID].objectType() == "IBBox")
    {
        IBBox& box = refCast<IBBox>
        (
            objects_[objectID]
        );
        forAll(mesh_.C(), cellI)
        {
            if 
            (
                mesh_.nGeometricD() == 2
             && mesh_.C()[cellI].x() > box.min().x() 
             && mesh_.C()[cellI].y() > box.min().y()
             && mesh_.C()[cellI].x() < box.max().x()
             && mesh_.C()[cellI].y() < box.max().y()
            )
            {
                slc.append(cellI);
            }
        }
    }
    if (objects_[objectID].objectType() == "IBWall")
    {
        slc.setSize(0);
    }
    return slc;
}

void Foam::IBObjectRegistry::readEnvironmentVariables()
{
    if (transportProperties_.found("phases"))
    {
        Pair<word> phases = transportProperties_.lookup("phases");
        dimensionedScalar rho1("rho", dimDensity, transportProperties_.subDict(phases[0]).lookup("rho"));
        dimensionedScalar rho2("rho", dimDensity, transportProperties_.subDict(phases[1]).lookup("rho"));
        if (rho1.value() > rho2.value())
            rhoF_ = rho1.value();
        else 
            rhoF_ = rho2.value();

        uniformDimensionedVectorField gravity
        (
            IOobject
            (
                "g",
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        g_ = gravity.value();
    }
    else
    {
        rhoF_ = transportProperties_.lookupOrDefault("rhoF",1000.0);
        g_ = transportProperties_.lookupOrDefault("g", vector::zero);
    }
}
bool Foam::IBObjectRegistry::isInsideRegion1(const label objectID, const point p)
{
    const boundBox& box = mesh_.bounds();
    const point& min = box.min();
    const point& max = box.max();
    if (bendingAngle_ == 0)
    {
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
    else
    {
        scalar anglePosition = getAngle(p);
        scalar R = R_[objectID] + 2.0*h_;
        scalar dAngle = R/innerRadius_;
        if (    anglePosition > dAngle 
             && anglePosition < (bendingAngle_ - dAngle) 
             && p.z() > min.z()     
             && p.z() < max.z() )
            return true;
        else 
            return false;
    }
} 

bool Foam::IBObjectRegistry::isInsideRegion2(const label objectID, const point p)
{
    const boundBox& box = mesh_.bounds();
    const point& min = box.min();
    const point& max = box.max();

    if (bendingAngle_ == 0)
    {
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
    else
    {
        scalar anglePosition = getAngle(p);
        if (    anglePosition > 0 
             && anglePosition < bendingAngle_ 
             && p.z() > min.z()     
             && p.z() < max.z() )
            return true;
        else 
            return false;
    }
} 

bool Foam::IBObjectRegistry::isInsideRegion3(const label objectID, const point p)
{
    const boundBox& box = mesh_.bounds();
    const point& min = box.min();
    const point& max = box.max();
    scalar R = R_[objectID] + 2.0*h_ +SMALL;

    if ( bendingAngle_ == 0)
    {
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
    else
    {
        scalar anglePosition = getAngle(p);
        scalar R = R_[objectID] + 2.0*h_;
        scalar dAngle = R/innerRadius_;

        if (    anglePosition > -dAngle 
             && anglePosition < (bendingAngle_ + dAngle) 
             && p.z() > min.z()     
             && p.z() < max.z() )
            return true;
        else 
            return false;
    }
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

            if (cpp1.parallel())
            {
                const boundBox& box = mesh_.bounds();
                const vector patchNorm = cpp1.faceNormals()[patchI];
                const vector maxToMin = box.max() - box.min();
                cyclicDistance_ = (maxToMin & patchNorm)*patchNorm;
                Info<< "  Two cyclic planes are parallel, distance = "
                    <<cyclicDistance_<<nl<<endl;
            }
            else
            {
                if (!Pstream::parRun())
                {
                    const cyclicPolyPatch& cpp2 = cpp1.neighbPatch();
                    const vector n1 = cpp1.faceNormals()[0];
                    const vector n2 = -cpp2.faceNormals()[0];
                    bendingAngle_ = Foam::vectorTools::radAngleBetween(n1,n2);
                    Info<< "  Curved channel with bending angle = "
                        <<bendingAngle_*180.0/3.1415<<" deg ("
                        <<bendingAngle_<<" rad)"<<nl<<endl;
                }
                else
                {
                    dictionary blockMeshDict
                    (
                        IOdictionary
                        (
                            IOobject
                            (
                                "blockMeshDict",
                                mesh_.time().system(),
                                mesh_,
                                IOobject::MUST_READ,
                                IOobject::NO_WRITE,
                                false
                            )
                        )
                    );
                    scalar bendingAngleDeg = readScalar(blockMeshDict.subDict("mesh").lookup("bendingAngle"));
                    bendingAngle_ = bendingAngleDeg*3.1415/180;
                    // bendingAngle_ = 3*3.1415/180;
                    Info<< "  Curved channel with bending angle = "
                        <<bendingAngleDeg<<" deg ("
                        <<bendingAngle_<<" rad)"<<nl<<endl;
                }
            
            }
            return;
        }
    }
}

void Foam::IBObjectRegistry::createObjects
(
    const dictionary& dict
)
{
    if (bendingAngle_ != 0)
    {
        innerRadius_ = readScalar(dict.lookup("innerRadius"));
    }

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
               
                if (bendingAngle_ == 0)
                {
                    //- Straight channel
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
                else
                {
                    //- Spiral channel
                    Shd1CoG_[i] = rotate(CoG_[i], -bendingAngle_);
                    Shd2CoG_[i] = rotate(CoG_[i],  bendingAngle_);
                    
                    Shd1LPoints_[i].setSize(LPoints_[i].size());
                    Shd2LPoints_[i].setSize(LPoints_[i].size());
                    forAll(LPoints_[i], pointI)
                    {
                        Shd1LPoints_[i][pointI] = rotate(LPoints_[i][pointI], -bendingAngle_);
                        Shd2LPoints_[i][pointI] = rotate(LPoints_[i][pointI],  bendingAngle_);
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
        if (objects_[i].objectType() == "IBSTL")
        {
            IBSTL& ibstl_ = refCast<IBSTL>
			(
			    objects_[i]
			);
            solidCells_[i] = ibstl_.solidCells();    
        }
        else
            solidCells_[i]  = findSolidCells(i,CoG_[i]);
        URotate_[i]         = objects_[i].rotationalVelocity();
        UTranslate_[i]      = objects_[i].translationalVelocity();
        nPointsOfFaces_[i]  = objects_[i].nPointsOfFaces();
        pointOfFace_[i]     = objects_[i].pointOfFace();
    }
}

template<class T>
Foam::List<T> Foam::IBObjectRegistry::make1DList(const List<List<T>>& twoDList)
{
    List<T> oneDList;
    if (twoDList.empty())
    {
        oneDList.setSize(0);
        return oneDList;
    }
    else
    {
        forAll(twoDList, i)
        {
            if (twoDList[i].empty())
                continue;
            else
            {
                forAll(twoDList[i], j)
                    oneDList.append(twoDList[i][j]);
            }
        }
        return oneDList;
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
    rhoF_(),
    g_(vector::zero),
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
    Shd2NeiCells_(),
    bendingAngle_(0.0),
    innerRadius_(GREAT)
{
    readEnvironmentVariables();
    detectCyclicBoundary();
	cartesianGridSize();
	createObjects(dict);
	getObjectsData();
    createShadows();
    writeNeighbourCells();
    
    dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).optionalSubDict(typeName + "Coeffs")
    );
    maxRefinement_ = readLabel(refineDict.lookup("maxRefinement"));

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

void Foam::IBObjectRegistry::writeNeighbourCells()
{
    volScalarField neighbourCells
    (
        IOobject
        (
            "neighbourCells",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("neighbourCells", dimless, 0)
    ); 

    
    for(int objI=0; objI< nObjects(); objI++)
    {
        forAll(LPoints_[objI], pointI)
        {
            forAll(neiCells_[objI][pointI], cellI)
            {
                neighbourCells[neiCells_[objI][pointI][cellI]] = 1.0;
            }
        }
        if (!Shd1NeiCells_.empty())
            forAll(Shd1NeiCells_[objI], pointI)
                if (!Shd1NeiCells_[objI][pointI].empty())
                    forAll(Shd1NeiCells_[objI][pointI], cellI)
                    {
                        neighbourCells[Shd1NeiCells_[objI][pointI][cellI]] = 1.0;
                    }


        if (!Shd2NeiCells_.empty())
            forAll(Shd2NeiCells_[objI], pointI)
                if (!Shd2NeiCells_[objI][pointI].empty())
                    forAll(Shd2NeiCells_[objI][pointI], cellI)
                    {
                        neighbourCells[Shd2NeiCells_[objI][pointI][cellI]] = 1.0;
                    }
    }
    
    cellsToRefine_.clear();
    forAll(neighbourCells, cellI)
    {
        if (neighbourCells[cellI] > 0)
            cellsToRefine_.append(cellI);
    }

    volScalarField sldCells
    (
        IOobject
        (
            "sldCells",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("sldCells", dimless, 0)
    ); 
    for(int objI=0; objI< nObjects(); objI++)
    {
        forAll(solidCells_[objI], cellI)
        {
            sldCells[solidCells_[objI][cellI]] = 1.0;
        }
    }
    if (mesh_.time().outputTime() || mesh_.time().timeIndex() == 0)
    {
        neighbourCells.write();
        sldCells.write();
    }
}
void Foam::IBObjectRegistry::updateCartesianGridSize()
{
    if (mesh_.time().timeIndex() == 1)
    {
        for (int i=0; i<maxRefinement_; i++)
            h_ /= 2.0;
        
        if (mesh_.nGeometricD() == 2)
            dV_ = h_*h_;
        else 
            dV_ = h_*h_*h_;
        
        Info<<"IBM: Update value:  h = "<<h_<<"; dV_ = "<<dV_<<endl;
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
