/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - twoPhaseDirectForcing class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "twoPhaseDirectForcing.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"
#include "pointPatchField.H"
#include "volPointInterpolation.H"
#include "pointVolInterpolation.H"
#include "IBSTL.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(twoPhaseDirectForcing, 0);
	addToRunTimeSelectionTable(IBTechnique, twoPhaseDirectForcing, dictionary);
}

//--------------------------Private Member Functions-------------------------//
void Foam::twoPhaseDirectForcing::makeIbForceUhlmann
(
    const volVectorField& U, 
    volVectorField& f
)
{
    if ( nObjects() != 0)
    {
        const scalar dT = mesh_.time().deltaTValue();
        const List<pointField>& LPoints_ = LPoints();
        const List<labelListList>& neiCells_ = neiCells();

        const List<pointField>& Shd1LPoints_ = Shd1LPoints();
        const List<labelListList>& Shd1NeiCells_ = Shd1NeiCells();

        const List<pointField>& Shd2LPoints_ = Shd2LPoints();
        const List<labelListList>& Shd2NeiCells_ = Shd2NeiCells();

        const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
        // const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

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
                    
                    f[neiCells_[objI][pointI][cellI]] 
                    += 
                        rho[neiCells_[objI][pointI][cellI]]
                    * FLagr[pointI] * deltaFuncValue * dV();
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
                            
                            f[Shd1NeiCells_[objI][pointI][cellI]] 
                            += 
                                rho[Shd1NeiCells_[objI][pointI][cellI]]
                            * FLagr[pointI] * deltaFuncValue * dV();
                        }
                    }
                    if (!Shd2NeiCells_[objI].empty())
                    {
                        forAll(Shd2NeiCells_[objI][pointI], cellI)
                        {
                            scalar deltaFuncValue = 
                                deltaFunc(mesh_.C()[Shd2NeiCells_[objI][pointI][cellI]], 
                                        Shd2LPoints_[objI][pointI]);
                            
                            f[Shd2NeiCells_[objI][pointI][cellI]] 
                            +=
                                rho[Shd2NeiCells_[objI][pointI][cellI]] 
                            * FLagr[pointI] * deltaFuncValue * dV();
                        }
                    }
                }
            }

            //- Save Lagrangian force for later use
            Fk_[objI] = FLagr;
        }
    }
}

void Foam::twoPhaseDirectForcing::makeIbForce_forcingPoints
(
    volVectorField& U, 
    volVectorField& f
)
{
    dimensionedScalar dT = mesh_.time().deltaT();

    //- Interpolate value of field at point mesh
    pointMesh pMesh(mesh_);
    volPointInterpolation pointData(mesh_);

    //- Interpolate velocity to pointMesh
    pointVectorField Up = pointData.interpolate(U);
    pointVectorField Updesire(Up);
    
    forAll(objects(), i)
    {
        if (objects()[i].objectType() == "IBSTL")
        {
            IBSTL& stlObject = refCast<IBSTL>(objects()[i]);

            //- Set velocity at solid points to be zero
            const labelList& sldPts = stlObject.solidPoints();
            
            forAll(sldPts, pointI)
            {
                Updesire[sldPts[pointI]] = vector::zero;
            }

            //- Interpolate velocity at virtual points and ibPoints
            const labelList& ibPts = stlObject.forcingPoints();
            const pointField& vtPts = stlObject.virtualPoints();
            const labelListList& neiPts = stlObject.neighbourPoints();

            forAll(vtPts, pointI)
            {
                vector Uvt = trilinearInterpolate(Up, vtPts[pointI], neiPts[pointI]);
                vector Uif = vector::zero;
                Updesire[ibPts[pointI]] = 0.5*(Uvt + Uif);
            }
        }
    }

    //- Interpolate velocity back to volMesh
    pointVolInterpolation volData(pMesh, mesh_);
    volVectorField Udesire = volData.interpolate(Updesire);

    //- Caluclate ibForce using desired velocity and predicted velocity
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    f = rho*(Udesire - U)/dT;
}

void Foam::twoPhaseDirectForcing::makeIbForce_forcingCells
(
    volVectorField& U, 
    volVectorField& f
)
{
    dimensionedScalar dT = mesh_.time().deltaT();
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho"); 
    volVectorField Udesire(U);

    forAll(objects(), i)
    {
        if (objects()[i].objectType() == "IBSTL")
        {
            IBSTL& stlObject = refCast<IBSTL>(objects()[i]);

            //- Set velocity at solid cells to be zero
            const labelList& slc = stlObject.solidCells();
            forAll(slc, cellI)
            {
                Udesire[slc[cellI]] = vector::zero;
            }

            const labelList& ibc = stlObject.ibCells();
            const labelListList& ibnc = stlObject.ibNeiCells();
            const labelList& vtpc = stlObject.virtualPointCells();
            const pointField& vtp = stlObject.virtualPoints();
            const pointField& ifp = stlObject.interfacePoints();

            //- Interpolate velocity at pointMesh
            volPointInterpolation pointData(mesh_);
            pointVectorField Up = pointData.interpolate(U);
            
            forAll(vtp, i)
            {
                const labelList& curCellPoint = mesh_.cellPoints()[vtpc[i]];

                //- Calculate velocity at virtual points
                vector Uvt = trilinearInterpolate(Up, vtp[i], curCellPoint);

                //- Calculate velocity at object's surface
                vector Uib = UTranslate()[i];

                //- Calculate velocity at ibCells by linear interpolation
                scalar h1 = mag(ifp[i] - mesh_.C()[ibc[i]]);
                scalar h2 = mag(mesh_.C()[ibc[i]] - vtp[i]);
                scalar h = h1 + h2;
                Udesire[ibc[i]] = h1/h*Uvt + h2/h*Uib;
            }
            
            //- Create force field
            f = rho*(Udesire - U)/dT;
            
            //- Spread force from ibCells to neighbour cells
            forAll(ibc, i)
            {
                forAll(ibnc[i], neiI)
                {
                    scalar deltaFuncValue = 
                        deltaFunc
                        (
                            mesh_.C()[ibnc[i][neiI]], 
                            mesh_.C()[ibc[i]]
                        );

                    f[ibnc[i][neiI]] += f[ibc[i]]*deltaFuncValue*dV();
                }
            }
        }
    }
}
//---------------------------------Constructors------------------------------//
Foam::twoPhaseDirectForcing::twoPhaseDirectForcing
(
	// const word& name,
	IBdynamicFvMesh& mesh,
	const dictionary& dict
)
:	
	IBTechnique(mesh, dict),
	mesh_(mesh),
	nMDF_(readScalar(dict.lookup("multiDirForcingIter"))),
    ibMethod_(dict.lookup("method"))
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

Foam::volVectorField Foam::twoPhaseDirectForcing::ibForce(volVectorField& U)
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
                dimensionSet(1,-2,-2,0,0,0,0), vector(0,0,0))
    );   
    
    if (ibMethod_ == "Uhlmann")
    {
        makeIbForceUhlmann(U, ibForce_);
    }
    else if (ibMethod_ == "forcingPoints")
    {
        makeIbForce_forcingPoints(U, ibForce_);
    }
    else if (ibMethod_ == "forcingCells")
    {
        makeIbForce_forcingCells(U, ibForce_);
    }
    else 
    {
        FatalErrorIn
        (
            "twoPhaseDirectForcing::ibForce()"
        )   << "Unknown twoPhaseDirectForcing method: "
            << ibMethod_<<nl<<"Valid methods are:"<<nl
            << "("<<nl<<"Uhlmann"<<nl<<"forcingCells"<<nl
            << "forcingPoints"<<nl<<")"<<nl
            <<abort(FatalError);
    }

	return ibForce_;
}

void Foam::twoPhaseDirectForcing::multiDirectForcing
(
    volVectorField& u,
    volVectorField& ibForce_
)
{
    // const labelList& solidCells_ = solidCells()[0];
    // Info<<"####################################################### BEFORE MULTI FORCING"<<endl;
    if(nMDF_ > 0)
    {
        dimensionedScalar dT("dT",dimTime,mesh_.time().deltaTValue());;
        const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    
        for(int i=0; i<nMDF_; i++)
        {
            Info<< "IBM: Multi-direct Forcing Iteration "<<i+1<<endl;
            volVectorField f = ibForce(u);

            u += dT*f/rho;
            ibForce_ += f;
        }
    }
    // Info<<"####################################################### After MULTI FORCING"<<endl;
}

void Foam::twoPhaseDirectForcing::update()
{
    // const labelList& solidCells_ = solidCells()[0];
    // const volVectorField& u = mesh_.lookupObject<volVectorField>("U");
    // Info<<"####################################################### AFTER pEqn.H"<<endl;
    // forAll(solidCells_, cellI)
    //     Info<<u[solidCells_[cellI]]<<endl;
    for(int objI=0; objI<nObjects(); objI++)
    {
        updateObjectMotions(objI, Fk_[objI]);
    }

    movePoints();

    writeNeighbourCells();
    // //- Refine mesh following neighbour cells indicator
    // mesh_.update(cellsToRefine());

    //- Update neighbour cells after mesh refinement
    // for(int i=0; i<nObjects(); i++)
    // {
        // neiCells()[i] = findNeiCells(LPoints()[i]);
    // }

    // //- Update h and dV (only do 1 time at timeIndex = 1)
    // updateCartesianGridSize();

    //- Write neighbour cells indicator
    // writeNeighbourCells();
}

void Foam::twoPhaseDirectForcing::write()
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
        Fk_[objI].clear();
        Fk_[objI].setSize(objects()[objI].nPoints(), vector::zero);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
