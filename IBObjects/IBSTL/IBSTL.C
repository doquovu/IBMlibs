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
#include "pointMesh.H"
#include "pointPatchField.H"
#include "meshSearch.H"
// * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IBSTL, 0);
    addToRunTimeSelectionTable(IBObject, IBSTL, dictionary);

// * * * * * * * * * * * * Private Member Fuctions * * * * * * * * * * * * * //
scalar IBSTL::cellSize(label cellID)
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

void IBSTL::readSTL(const dictionary& dict)
{
    Info<<"  Reading "<<STLFile_<<nl<<endl;

    //- For running ship model only
    forAll(triSurf_.localPoints(), i)
    {
        if (triSurf_.localPoints()[i].y() <= 0)
        {
            lPoints_.append(triSurf_.localPoints()[i]);
        }
    }
    // lPoints_ = triSurf_.localPoints();
    nPoints_ = lPoints_.size();
    nFaces_  = triSurf_.faceCentres().size();
    nPointsOfFaces_.setSize(triSurf_.localPoints().size(), 3); 
    pointOfFace_.setSize(triSurf_.localPoints().size());

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

void IBSTL::makeInterpolationStencil()
{
    //- Find solid cells and mark them with gamma = 1
        scalarField& gammaI = gamma_.ref();

        const triSurfaceSearch& tss(triSurf_);
        
        boolList inside = tss.calcInside(mesh_.cellCentres());
        
        forAll(gammaI, cellI)
        {
            if (inside[cellI])
            {
                gammaI[cellI] = 1;
                solidCells_.append(cellI);
            }
        }
    
    //- Find ib cells
        labelHashSet ibCellSet;
        const unallocLabelList& owner = mesh_.owner();
        const unallocLabelList& neighbour = mesh_.neighbour();

        forAll(neighbour, faceI)
        {
            scalar dGamma 
            = gammaI[neighbour[faceI]] - gammaI[owner[faceI]];
            if (mag(dGamma) > SMALL)
            {
                if (gammaI[owner[faceI]] > SMALL)
                {
                    if (!ibCellSet.found(neighbour[faceI]))
                    {
                        ibCellSet.insert(neighbour[faceI]);
                    }
                }
                else
                {
                    if (!ibCellSet.found(owner[faceI]))
                    {
                        ibCellSet.insert(owner[faceI]);
                    }
                }
            }
        }
        ibCells_ = labelList(ibCellSet.toc());
        
    //- Find interface points correspond to ibCells centre
        ifPts_.setSize(ibCells_.size());
        const pointField ibcc(mesh_.C(), ibCells_);
        List<pointIndexHit> nearInfo;
        triSurf_.findNearest(ibcc, scalarField(ibCells_.size(), GREAT), nearInfo);
        label nMiss = 0;

        forAll(nearInfo, i)
        {
            const pointIndexHit& pi = nearInfo[i];
            if (pi.hit())
            {
                ifPts_[i] = pi.hitPoint();
            }
            else 
            {
                nMiss++;
                ifPts_[i] = ifPts_[i];
            }

            if (nMiss > 0)
            {
                FatalErrorIn
                (
                    "IBSTL::makeInterpolationStencilNew()"
                )   << "Error projecting ibCell centre to surface:"
                    << nMiss <<"faces could not be determined"<<nl
                    <<abort(FatalError);
            }
        }

    //- Find virtual points corresponding to each ibCells centre
        vtPts_.setSize(ibCells_.size());
        
        forAll(vtPts_, i)
        {
            vector h1 = mesh_.C()[ibCells_[i]] - ifPts_[i];
            vtPts_[i] = ifPts_[i] + 2.0*h1;
        }

        Pout<<"@@@@@@@@@ found "<<vtPts_.size()<<" virtual points."<<endl;

        // List<pointField> procVtPts(Pstream::nProcs());
        // procVtPts[Pstream::myProcNo()] = vtPts_;

        // Pstream::gatherList(procVtPts);
        // Pstream::scatterList(procVtPts);

    //- Find cell enclosing each virtual point
        vtPointCells_.setSize(vtPts_.size());
        meshSearch ms(mesh_);

        forAll(vtPointCells_, pointI)
        {
            label cellID = ms.findCell(vtPts_[pointI]);
            if (cellID != -1)
            {
                vtPointCells_[pointI] = cellID;
            }
            else 
            {
                FatalErrorIn
                (
                    "IBSTL::makeInterpolationStencil()"
                )   << "Error: cannot find cell enclosing virtual point "<<pointI
                    << ": "<<vtPts_[pointI]<<nl
                    <<abort(FatalError);
            }

        }

        // List<List<labelList>> procVtPtsCell(Pstream::nProcs());
        // procVtPtsCell[Pstream::myProcNo()] = vtPointCells_;

        // Pstream::gatherList(procVtPtsCell);
        // Pstream::scatterList(procVtPtsCell);

        // for (label procI=0; procI<Pstream::nProcs(); procI++)
        // {
        //     if (procI != Pstream::myProcNo())
        //     {
        //         forAll(procVtPtsCell[procI], pI)
        //         {
        //             if (procVtPtsCell[procI][pI][0] == -1)
        //             {
        //                 point curPoint = procVtPts[procI][pI];
        //                 label cellID = ms.findCell(curPoint);
        //                 if (cellID != -1)
        //                 {
        //                     procVtPtsCell[procI][pI][0] = Pstream::myProcNo();
        //                     procVtPtsCell[procI][pI][1] = cellID;
        //                 }
        //             }
        //             reduce(procVtPtsCell[procI][pI], maxOp<labelList>());
        //         }
        //     }
        // }
        // vtPointCells_ = procVtPtsCell[Pstream::myProcNo()];


    // - Find neighbour cells of ibCells
        const labelListList& cellPoints = mesh_.cellPoints();
        const labelListList& pointCells = mesh_.pointCells();
        ibNeiCells_.setSize(ibCells_.size());
        
        forAll(ibCells_, i)
        {
            labelHashSet cellSet;

            const labelList& curCellPoints = cellPoints[ibCells_[i]];

            forAll (curCellPoints, pointI)
            {
                label curPoint = curCellPoints[pointI];
                const labelList& curPointCells = pointCells[curPoint];

                forAll (curPointCells, cI)
                {
                    if ( gammaI[curPointCells[cI]] > SMALL)
                    {
                        continue;
                    }
                    else
                    {
                        if (!cellSet.found(curPointCells[cI]))
                        {
                            cellSet.insert(curPointCells[cI]);
                        }
                    }
                }
            }

            cellSet.erase(ibCells_[i]);
            ibNeiCells_[i] = cellSet.toc();
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
    STLFile_(word(dict.lookup("fileName"))),
    triSurf_ 
    (
        IOobject
        (
            STLFile_,
            mesh_.time().constant(),
            "triSurface",
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
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
    gamma_
    (
        IOobject
        (
            "gamma",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0)
    )
{
    readSTL(dict);
    Info <<"  - Center : "<< center_ <<endl;
    Info <<"  - Radius : "<< R_ <<endl;
    Info <<"  - nPoints: "<< nPoints_ <<endl;
    Info <<"  - Density: "<< rho_ <<endl;
    addMotions(dict);
    makeInterpolationStencil();
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