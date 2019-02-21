/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - directForcing class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "directForcing.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(directForcing, 0);
	addToRunTimeSelectionTable(IBTechnique, directForcing, dictionary);
}

//---------------------------------Constructors------------------------------//

Foam::directForcing::directForcing
(
	// const word& name,
	IBdynamicFvMesh& mesh,
	const dictionary& dict
)
:	
	IBTechnique(mesh, dict),
	mesh_(mesh),
	nMDF_(readScalar(dict.lookup("multiDirForcingIter")))
{
    Fk_.setSize(nObjects());
    forAll(Fk_, i)
    {
        Fk_[i].setSize(objects()[i].nPoints(), vector::zero);
    }
    Info<<"Cartesian grid size :  h = "<<h()<<nl
        <<"Cartesian volumetric: dV = "<<dV()<<nl<<endl;
}

// -------------------------------Member Functions----------------------------//

Foam::volVectorField Foam::directForcing::ibForce(volVectorField& U)
{
    Info<< "IBM: Calculating ibForce ..."<<endl;
    
    volVectorField ibForce_
    (
        IOobject
        (
            "ibForce",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("ibForce",
                dimensionSet(0,1,-2,0,0,0,0), vector(0,0,0))
    );   
    if ( nObjects() != 0)
    {
        const scalar dT = mesh_.time().deltaTValue();
        const List<pointField>& LPoints_ = LPoints();
        const List<labelListList>& neiCells_ = neiCells();

        const List<pointField>& Shd1LPoints_ = Shd1LPoints();
        const List<labelListList>& Shd1NeiCells_ = Shd1NeiCells();

        const List<pointField>& Shd2LPoints_ = Shd2LPoints();
        const List<labelListList>& Shd2NeiCells_ = Shd2NeiCells();

        for(int objI=0; objI< nObjects(); objI++)
        {
            vectorField& Ub = UBoundary()[objI];

            vectorField ULagr(LPoints_[objI].size(), vector::zero);        
            vectorField FLagr(LPoints_[objI].size(), vector::zero); 
            forAll(LPoints_[objI], pointI)
            {

            //- Interpolating Velocity at Lagragian points
                forAll(neiCells_[objI][pointI], cellI)
                {
                    scalar deltaFuncValue = 
                        deltaFunc
                        (
                            mesh_.C()[neiCells_[objI][pointI][cellI]], 
                            LPoints_[objI][pointI]
                        );
                    
                    ULagr[pointI] 
                        += U[neiCells_[objI][pointI][cellI]] * deltaFuncValue * dV();
                    
                }
                if (enableShadows()[objI])
                {   
                    if (!Shd1NeiCells_[objI].empty())
                    {
                        forAll(Shd1NeiCells_[objI][pointI], cellI)
                        {
                            scalar deltaFuncValue = 
                                deltaFunc
                                (
                                    mesh_.C()[Shd1NeiCells_[objI][pointI][cellI]], 
                                    Shd1LPoints_[objI][pointI]
                                );
                            
                            ULagr[pointI] 
                                += U[Shd1NeiCells_[objI][pointI][cellI]] * deltaFuncValue * dV();
                        }
                    }
                    if (!Shd2NeiCells_[objI].empty())
                    {
                        forAll(Shd2NeiCells_[objI][pointI], cellI)
                        {
                            scalar deltaFuncValue = 
                                deltaFunc
                                (
                                    mesh_.C()[Shd2NeiCells_[objI][pointI][cellI]], 
                                    Shd2LPoints_[objI][pointI]
                                );
                            
                            ULagr[pointI] 
                                += U[Shd2NeiCells_[objI][pointI][cellI]] * deltaFuncValue * dV();
                        }
                    }
                }
                if (mesh_.nGeometricD() == 2)
                {
                    ULagr[pointI].z() = 0;
                }

                reduce(ULagr[pointI], sumOp<vector>());

            //- Calculate Force at Lagrang points
                FLagr[pointI] = ( Ub[pointI] - ULagr[pointI] ) / dT;

            //- Spread Force from Lagrang points to Euler cells.
                forAll(neiCells_[objI][pointI], cellI)
                {
                    scalar deltaFuncValue = 
                        deltaFunc
                        (
                            mesh_.C()[neiCells_[objI][pointI][cellI]], 
                            LPoints_[objI][pointI]
                        );
                    
                    ibForce_[neiCells_[objI][pointI][cellI]] 
                        += FLagr[pointI] * deltaFuncValue * dV();
                }
                if (enableShadows()[objI])
                {   
                    if (!Shd1NeiCells_[objI].empty())
                    {
                        forAll(Shd1NeiCells_[objI][pointI], cellI)
                        {
                            scalar deltaFuncValue = 
                                deltaFunc(mesh_.C()[Shd1NeiCells_[objI][pointI][cellI]], 
                                        Shd1LPoints_[objI][pointI]);
                            
                            ibForce_[Shd1NeiCells_[objI][pointI][cellI]] 
                                += FLagr[pointI] * deltaFuncValue * dV();
                        }
                    }
                    if (!Shd2NeiCells_[objI].empty())
                    {
                        forAll(Shd2NeiCells_[objI][pointI], cellI)
                        {
                            scalar deltaFuncValue = 
                                deltaFunc(mesh_.C()[Shd2NeiCells_[objI][pointI][cellI]], 
                                        Shd2LPoints_[objI][pointI]);
                            
                            ibForce_[Shd2NeiCells_[objI][pointI][cellI]] 
                                += FLagr[pointI] * deltaFuncValue * dV();
                        }
                    }
                }
            }

            //- Save Lagrangian force for later use
            Fk_[objI] = FLagr;
        }
    }   
	return ibForce_;
}
Foam::volVectorField Foam::directForcing::ibForceInt()
{
    Info<< "IBM: Calculating ibForceInt..."<<endl;
    
    volVectorField ibForce_
    (
        IOobject
        (
            "ibForce",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("ibForce",
                dimensionSet(1,-2,-2,0,0,0,0), vector(0,0,0))
    );   
    volScalarField TSolidCells
    (
        IOobject
        (
            "TSolidCells",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("TSolidCells", dimless, 0)
    );

    const scalar dT = mesh_.time().deltaTValue();
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");

    for(int objI=0; objI< nObjects(); objI++)
    {
        labelList& solidCells_ = solidCells()[objI];
        vector US = UTranslate()[objI];
        forAll(solidCells_, cellI)
        {
            ibForce_[solidCells_[cellI]] = rho[solidCells_[cellI]]*(US - U[solidCells_[cellI]])/dT;
            TSolidCells[solidCells_[cellI]] = 1;
        }
    }
    
    // Fk_[objI] += FLagr;

    if (mesh_.time().outputTime())
    {
        TSolidCells.write();
    }

    return ibForce_;
}
void Foam::directForcing::multiDirectForcing
(
    volVectorField& u,
    volVectorField& ibForce_
)
{
    if(nMDF_ > 0)
    {
        dimensionedScalar dT("dT",dimTime,mesh_.time().deltaTValue());;

        for(int i=0; i<nMDF_; i++)
        {
            Info<< "IBM: Multi-direct Forcing Iteration "<<i+1<<endl;
            volVectorField f = ibForce(u);

            u += dT*f;
        }
    }
}

void Foam::directForcing::update()
{
    //- Calculate velocity of Lagrang points, 
    //  move points accordingly, and update
    //  neighbour cells
    if (mesh_.time().timeIndex() > 0)
    {
        for(int objI=0; objI<nObjects(); objI++)
        {
            updateObjectMotions(objI, Fk_[objI]);
        }
    
        movePoints();
    }

    writeNeighbourCells();
    // //- Refine mesh following neighbour cells indicator
    // mesh_.update(cellsToRefine());

    // //- Update neighbour cells after mesh refinement
    // for(int i=0; i<nObjects(); i++)
    // {
    //     neiCells()[i] = findNeiCells(LPoints()[i]);
    // }

    // //- Update h and dV (only do 1 time at timeIndex = 1)
    // updateCartesianGridSize();

    // //- Write neighbour cells indicator
    // writeNeighbourCells();
}

void Foam::directForcing::write()
{
    for(int objI=0; objI<nObjects(); objI++)
    {
        writeLagrPoints(objI, LPoints()[objI]);
        writeLagrForces(objI, LPoints()[objI], Fk_[objI]);
        writeIBForces(objI, Fk_[objI], dV());
        writeObjData(objI, CoG()[objI], UTranslate()[objI], URotate()[objI]);
    
        writeObjVTU
        (
            objI, 
            LPoints()[objI], 
            CoG()[objI], 
            nFaces()[objI], 
            nPointsOfFaces()[objI], 
            pointOfFace()[objI]
        );
        if (enableShadows()[objI])
        {
            writeShadPoints(objI, Shd1LPoints()[objI], Shd2LPoints()[objI]);
        }
        // Fk_[objI].clear();
        // Fk_[objI].setSize(objects()[objI].nPoints(), vector::zero);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
