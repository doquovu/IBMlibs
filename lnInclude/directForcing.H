/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - directForcing class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef directForcing_H
#define directForcing_H

#include "IBTechnique.H"
#include "IBObject.H"
#include "IBObjectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
							Class IBModel declaration
\*---------------------------------------------------------------------------*/

class directForcing
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
	directForcing(const directForcing&);
	
	//- Disallow default bitwis assignment
	void operator=(const directForcing&);

	
public:

	//- Runtime type information
	TypeName("directForcing");

	//- Constructor
		directForcing
		(
			// const word& name,
			const fvMesh& mesh,
			const dictionary& dict
		);

	//- Destructor
		virtual ~directForcing()
		{}

	//- Member functions

		inline List<vectorField> Fk()
		{
			return Fk_;
		}

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
