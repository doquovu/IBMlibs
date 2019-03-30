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
		IBdynamicFvMesh& mesh_;
		
		//- Number of multi-direct forcing iteration
		label nMDF_;
		
		//- Lagrangian force field
		List<vectorField> Fk_;

		//- Direct forcing method
		word ibMethod_;
		
		void makeIbForceUhlmann
		(
			const volVectorField& U, 
			volVectorField& f
		);

		void makeIbForce_forcingCells
		(
			volVectorField& U, 
			volVectorField& f
		);

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
			IBdynamicFvMesh& mesh,
			const dictionary& dict
		);

	//- Destructor
		virtual ~twoPhaseDirectForcing()
		{}

	//- Member functions

		//- Return force field for momentum equation
		virtual volVectorField ibForce(volVectorField& U);

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
