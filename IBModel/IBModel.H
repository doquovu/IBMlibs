/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBModel class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef IBModel_H
#define IBModel_H

#include "autoPtr.H"
#include "IBTechnique.H"
#include "IBdynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
							Class IBModel declaration
\*---------------------------------------------------------------------------*/

class IBModel
:
	IOdictionary
{
	//- Pointer to chose approaching method
	autoPtr<IBTechnique> techniquePtr_;

	// //- Disallow default bitwise copy construct
	// IBModel(const IBModel&);		
	
	// //- Disallow default bitwise assignment
	// void operator=(const IBModel&);
	
public:
	
	TypeName("IBModel");
	
	IBModel(IBdynamicFvMesh& mesh);

	virtual ~IBModel()
	{}

	virtual volVectorField pressureGradientField();

	virtual volVectorField ibForce(volVectorField& U);

	virtual volVectorField ibForceInt();

	virtual volVectorField ibForceInt(const volVectorField& rhs);

	virtual void multiDirectForcing
	(
		volVectorField& u,
		volVectorField& ibForce
	);

	virtual void update();

	virtual void write();

}; //End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
