/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - twoPhaseDirectForcingBallaras class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "twoPhaseDirectForcingBallaras.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"
#include "pointPatchField.H"
#include "volPointInterpolation.H"
#include "pointVolInterpolation.H"
#include "IBSTL.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(twoPhaseDirectForcingBallaras, 0);
	addToRunTimeSelectionTable(IBTechnique, twoPhaseDirectForcingBallaras, dictionary);
}

//---------------------------Private Member Functions------------------------//
Foam::vector Foam::twoPhaseDirectForcingBallaras::linearInterpolate
(
    const point p1,
    const point p2,
    const point target,
    const vector v1,
    const vector v2
)
{
    vector targetValue(0,0,0);
    scalar h1 = mag(p1 - target);
    scalar h2 = mag(p2 - target);
    scalar h = h1 + h2;

    targetValue = h2/h*v1 + h1/h*v2;

    return targetValue;
}

Foam::vector Foam::twoPhaseDirectForcingBallaras::bilinearInterpolate
(
    const pointVectorField& phi,
    const point& target,
    const labelList& sourceLabel
)
{
    //- Assuming the interpolation is perfom in the XOY plane
    //- Check number of points
    if (sourceLabel.size() != 4)
    {
        FatalErrorIn
        (
            "twoPhaseDirectForcingBallaras::bilinearInterpolate()"
        )   << "Error: bilinearInterpolate needs 4 points, "
            << sourceLabel.size() <<" points provided"<<nl
            <<abort(FatalError);
    }
    //- Find intermediate points
    const pointField source(mesh_.points(), sourceLabel);

    Pair<label> pair1(0,0);
    Pair<label> pair2(0,0);

    pair1.first() = sourceLabel[0];

    for (int i=1; i<source.size(); i++)
    {
        if (source[i].y() == source[0].y())
        {
            pair1.second() = sourceLabel[i];

            if (i == 1)
            {
                pair2.first() = sourceLabel[2];
                pair2.second() = sourceLabel[3];
            }
            else if (i == 2)
            {
                pair2.first() = sourceLabel[1];
                pair2.second() = sourceLabel[3];
            }
            else 
            {
                pair2.first() = sourceLabel[1];
                pair2.second() = sourceLabel[2];
            }
            break;
        }
    }

    point p1(target.x(), mesh_.points()[pair1.first()].y(), target.z());
    point p2(target.x(), mesh_.points()[pair2.first()].y(), target.z());
    vector v1 
    = 
        linearInterpolate
        (
            mesh_.points()[pair1.first()],
            mesh_.points()[pair1.second()],
            p1,
            phi[pair1.first()],
            phi[pair1.second()]
        );
    vector v2
    = 
        linearInterpolate
        (
            mesh_.points()[pair2.first()],
            mesh_.points()[pair2.second()],
            p2,
            phi[pair2.first()],
            phi[pair2.second()]
        );

    return (linearInterpolate(p1, p2, target, v1, v2));
}

Foam::vector Foam::twoPhaseDirectForcingBallaras::trilinearInterpolate
(
    const pointVectorField& phi,
    const point& target,
    const labelList& sourceLabel
)
{
    //- Check number of points
    if (sourceLabel.size() != 8)
    {
        FatalErrorIn
        (
            "twoPhaseDirectForcingBallaras::trilinearInterpolate()"
        )   << "Error: trilinearInterpolate needs 8 points, "
            << sourceLabel.size() <<" points provided"<<nl
            <<abort(FatalError);
    }
    const pointField source(mesh_.points(), sourceLabel);
    labelList pointStage1;
    labelList pointStage2;
    
    pointStage1.append(sourceLabel[0]);

    for (int i=1; i<source.size(); i++)
    {
        if (source[i].z() == source[0].z())
        {
            pointStage1.append(sourceLabel[i]);
        }
        else 
        {
            pointStage2.append(sourceLabel[i]);
        }
    }

    point p1(target.x(), target.y(), mesh_.points()[pointStage1[0]].z());
    point p2(target.x(), target.y(), mesh_.points()[pointStage2[0]].z());
    vector v1 = bilinearInterpolate(phi, p1, pointStage1);
    vector v2 = bilinearInterpolate(phi, p2, pointStage2);
    
    return (linearInterpolate(p1, p2, target, v1, v2));
}

void Foam::twoPhaseDirectForcingBallaras::makeIbForce
(
    volVectorField& U, 
    volVectorField& f
)
{
    dimensionedScalar dT = mesh_.time().deltaT();

    //- Interpolate value of field at point mesh
    pointMesh pMesh(mesh_);
    volPointInterpolation pointData(mesh_);

    pointVectorField Up = pointData.interpolate(U);
    pointVectorField Updesire(Up);
    
    forAll(objects(), i)
    {
        if (objects()[i].objectType() == "IBSTL")
        {
            IBSTL& stlObject = refCast<IBSTL>(objects()[i]);

            //- Set velocity at solid points to be zero
            const labelList& sldPts = stlObject.solidPoints();
            const labelList& ibPts = stlObject.forcingPoints();
            // const pointField& ifPts = stlObject.interfacePoints();
            const pointField& vtPts = stlObject.virtualPoints();
            const labelListList& neiPts = stlObject.neighbourPoints();

            forAll(sldPts, pointI)
            {
                Updesire[sldPts[pointI]] = vector::zero;
            }

            //- Interpolate velocity at virtual points and ibPoints
            forAll(vtPts, pointI)
            {
                vector Uvt = trilinearInterpolate(Up, vtPts[pointI], neiPts[pointI]);
                vector Uif = vector::zero;
                Updesire[ibPts[pointI]] = 0.5*(Uvt + Uif);
            }
        }
    }

    pointVolInterpolation volData(pMesh, mesh_);
    volVectorField Udesire = volData.interpolate(Updesire);

    f = (Udesire - U)/dT;
}

void Foam::twoPhaseDirectForcingBallaras::makeIbForceNew
(
    volVectorField& U, 
    volVectorField& f
)
{
    dimensionedScalar dT = mesh_.time().deltaT();
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho"); 
    volVectorField Udesire(U);
    //- Make desire velocity field
    forAll(objects(), i)
    {
        if (objects()[i].objectType() == "IBSTL")
        {
            IBSTL& stlObject = refCast<IBSTL>(objects()[i]);

            const labelList& slc = stlObject.solidCells();
            const labelList& ibc = stlObject.ibCells();
            const labelListList& ibnc = stlObject.ibNeiCells();
            const labelList& vtpc = stlObject.virtualPointCells();
            const pointField& vtp = stlObject.virtualPoints();
            const pointField& ifp = stlObject.interfacePoints();

            //- Set velocity at solid cells to be zero
            forAll(slc, cellI)
            {
                Udesire[slc[cellI]] = vector::zero;
            }

            //- Interpolate velocity at virtual points
            //  then calculate desired velocity at ibCells
            volPointInterpolation pointData(mesh_);
            pointVectorField Up = pointData.interpolate(U);
            forAll(vtp, i)
            {
                const labelList& curCellPoint = mesh_.cellPoints()[vtpc[i]];

                vector Uvt = trilinearInterpolate(Up, vtp[i], curCellPoint);
                vector Uib = UTranslate()[i];

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
Foam::twoPhaseDirectForcingBallaras::twoPhaseDirectForcingBallaras
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

Foam::volVectorField Foam::twoPhaseDirectForcingBallaras::ibForce(volVectorField& U)
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
        dimensionedVector("zero",
                dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
    );   
    if ( nObjects() != 0)
    {
        // makeIbForce(U, ibForce_);
        makeIbForceNew(U, ibForce_);
    }   

    if (mesh_.time().outputTime())
    {
        ibForce_.write();
    }
	return ibForce_;
}

Foam::volVectorField Foam::twoPhaseDirectForcingBallaras::ibForceInt()
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
Foam::volVectorField Foam::twoPhaseDirectForcingBallaras::ibForceInt(const volVectorField& rhs)
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
    
    return ibForce_;
}

void Foam::twoPhaseDirectForcingBallaras::multiDirectForcing
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

void Foam::twoPhaseDirectForcingBallaras::update()
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
    mesh_.update(cellsToRefine());

    //- Update neighbour cells after mesh refinement
    for(int i=0; i<nObjects(); i++)
    {
        neiCells()[i] = findNeiCells(LPoints()[i]);
    }

    // //- Update h and dV (only do 1 time at timeIndex = 1)
    updateCartesianGridSize();

    //- Write neighbour cells indicator
    writeNeighbourCells();
}

void Foam::twoPhaseDirectForcingBallaras::write()
{
    for(int objI=0; objI<nObjects(); objI++)
    {
        writeLagrPoints(objI, LPoints()[objI]);
        writeLagrForces(objI, LPoints()[objI], Fk_[objI]);
        writeObjData(objI, CoG()[objI], UTranslate()[objI], URotate()[objI]);

        if (bendingAngle() == 0)
        {
            writeIBForces(objI, Fk_[objI], dV());
        }
        else
        {
            writeIBForces(objI, Fk_[objI], dV(), CoG()[objI]);
        }

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
        if (objects()[objI].objectType() == "IBSTL")
        {
            IBSTL& stlObject = refCast<IBSTL>(objects()[objI]);
            const pointField fp(mesh_.points(), stlObject.forcingPoints());
            const pointField sp(mesh_.points(), stlObject.solidPoints());
            const pointField ip = stlObject.interfacePoints();
            const pointField vp = stlObject.virtualPoints();

            // writePointSet(fp, "forcingPoints");
            // writePointSet(sp, "solidPoints");
            if (mesh_.time().outputTime())
            {
                writePointSet(ip, "interfacePoints");
                writePointSet(vp, "virtualPoints");
                stlObject.gamma().write();
            }

            // List<pointField> interpolateStencil;
            // interpolateStencil.setSize(fp.size());
            // forAll(interpolateStencil, ii)
            // {
            //     interpolateStencil[ii].append(fp[ii]);
            //     interpolateStencil[ii].append(ip[ii]);
            //     forAll(stlObject.neighbourPoints()[ii], nii)
            //     {
            //         interpolateStencil[ii].append(mesh_.points()[stlObject.neighbourPoints()[ii][nii]]);
            //     }
            //     word name = "point"+Foam::name(ii);
            //     writePointSet(interpolateStencil[ii], name);
            // }

        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
