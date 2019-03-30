/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBModel class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "IBModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(IBModel, 0);
}
//---------------------------------Constructors------------------------------//

Foam::IBModel::IBModel
(
	IBdynamicFvMesh& mesh
)
:	
	IOdictionary
	(
		IOobject
		(
			"IBMProperties",
			mesh.time().constant(),
			mesh,
        	IOobject::MUST_READ,
        	IOobject::NO_WRITE
		)
	),
	techniquePtr_(IBTechnique::New(mesh, subDict("parameters")))
{}

// -------------------------------Member Functions----------------------------//
Foam::volVectorField Foam::IBModel::pressureGradientField()
{
	return techniquePtr_->pressureGradField();
}

Foam::volVectorField Foam::IBModel::ibForce(volVectorField& U)
{
   return techniquePtr_->ibForce(U);
}

void Foam::IBModel::multiDirectForcing
(
	volVectorField& u,
	volVectorField& ibForce
)
{
	return techniquePtr_->multiDirectForcing(u, ibForce);
}

void Foam::IBModel::update()
{
	techniquePtr_->update();
}

void Foam::IBModel::write()
{
	techniquePtr_->write();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
