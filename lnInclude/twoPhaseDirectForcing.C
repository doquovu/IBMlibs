/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - twoPhaseDirectForcing class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "twoPhaseDirectForcing.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(twoPhaseDirectForcing, 0);
	addToRunTimeSelectionTable(IBTechnique, twoPhaseDirectForcing, dictionary);
}

//---------------------------------Constructors------------------------------//

Foam::twoPhaseDirectForcing::twoPhaseDirectForcing
(
	// const word& name,
	const fvMesh& mesh,
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

Foam::volVectorField Foam::twoPhaseDirectForcing::ibForce()
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
    volScalarField TneiCells
    (
        IOobject
        (
            "TneiCells",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("TneiCells", dimless, 0)
    );

    const scalar dT = mesh_.time().deltaTValue();
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");

    for(int objI=0; objI< nObjects(); objI++)
    {
        labelList& solidCells_ = solidCellsInt()[objI];
        vector US = UTranslate()[objI];
        forAll(solidCells_, cellI)
        {
            ibForce_[solidCells_[cellI]] = rho[solidCells_[cellI]]*(US - U[solidCells_[cellI]])/dT;
        }
    }
    
    // Fk_[objI] += FLagr;

    if (mesh_.time().outputTime())
    {
        TneiCells.write();
    }

    return ibForce_;
}

void Foam::twoPhaseDirectForcing::multiDirectForcing
(
    volVectorField& u,
    volVectorField& ibForce_
)
{
const labelList& solidCellsInt_ = solidCellsInt()[0];
    Info<<"####################################################### BEFORE MULTI FORCING"<<endl;
forAll(solidCellsInt_, cellI)
Info<<u[solidCellsInt_[cellI]]<<endl;
    if(nMDF_ > 0)
    {
        dimensionedScalar dT("dT",dimTime,mesh_.time().deltaTValue());;
        const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    
        for(int i=0; i<nMDF_; i++)
        {
            Info<< "IBM: Multi-direct Forcing Iteration "<<i+1<<endl;
            volVectorField f = ibForce();

            u += dT*f/rho;
            ibForce_ += f;
        }
    }
Info<<"####################################################### After MULTI FORCING"<<endl;
forAll(solidCellsInt_, cellI)
Info<<u[solidCellsInt_[cellI]]<<endl;
}

void Foam::twoPhaseDirectForcing::update()
{
    const labelList& solidCellsInt_ = solidCellsInt()[0];
    const volVectorField& u = mesh_.lookupObject<volVectorField>("U");
    Info<<"####################################################### AFTER pEqn.H"<<endl;

    forAll(solidCellsInt_, cellI)
        Info<<u[solidCellsInt_[cellI]]<<endl;
    for(int objI=0; objI<nObjects(); objI++)
    {
        updateObjectMotions(objI, Fk_[objI]);
    }

    movePoints();
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
