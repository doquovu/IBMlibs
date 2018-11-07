/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - twoPhaseDirectForcing class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef twoPhaseDirectForcing_H
#define twoPhaseDirectForcing_H

#include "IBTechnique.H"
#include "IBObject.H"
#include "IBObjectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
					 Class twoPhaseDirectForcing declaration
\*---------------------------------------------------------------------------*/

class twoPhaseDirectForcing
: 
	public IBTechnique
{
	//- Private data

		//- Reference to fvMesh
		const fvMesh& mesh_;
		
		//- Number of multi-direct forcing iteration
		label nMDF_;
		
		List<vectorField> Fk_;

	//- Disallow default bitwise copy construct
	twoPhaseDirectForcing(const twoPhaseDirectForcing&);
	
	//- Disallow default bitwis assignment
	void operator=(const twoPhaseDirectForcing&);

	
public:

	//- Runtime type information
	TypeName("twoPhaseDirectForcing");

	//- Constructor
		twoPhaseDirectForcing
		(
			// const word& name,
			const fvMesh& mesh,
			const dictionary& dict
		);

	//- Destructor
		virtual ~twoPhaseDirectForcing()
		{}

	//- Member functions

		//- Return force field for momentum equation
		virtual volVectorField ibForce();

		//- Perform multi-Firect forcing
		virtual void multiDirectForcing
		(
			volVectorField& u,
			volVectorField& ibForce
		);
		
		//- Update object's velocity and position 
		virtual void update();

		//- Write data for post-processing
		virtual void write();

}; //End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
