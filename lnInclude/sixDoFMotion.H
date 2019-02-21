/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - sixDoFMotion class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef sixDoFMotion_H
#define sixDoFMotion_H

#include "IBMotion.H"
#include "IBObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
						Class sixDoFMotion declaration
\*---------------------------------------------------------------------------*/
class sixDoFMotion
:
	public IBMotion
{
	//- Private data

		//- Reference to Immersed object
		IBObject& obj_;

		//- Solver to calculate velocities
		word solver_;
		
		//- Translational velocity
		vector UTranslate_;

		//- Rotational velocity
		vector URotate_;

		//- Constraint
		vector constraintDir_;

		void initialise();

		void readConstraint(const dictionary& dict);

public:

	//- Runtime type information
	TypeName("sixDoFMotion");

	//- Constructor
		sixDoFMotion
		(
			const word& name,
			IBObject& obj,
			const dictionary& dict
		);

	//- Destructor
		virtual ~sixDoFMotion()
		{}

	//- Member functions

		//- Solve Newton equation of motion
		void updateMotion
		(
			vector& uTransl,
			vector& uRotate,
			const fvMesh& mesh,
			const vectorField& Fk, 
			vector repusliveForce,
			point center,
			pointField lPoints,
			scalar rhoF,
			scalar dV,
			vector g
		);
		void updateMotion
		(
			vector& uTransl,
			vector& uRotate,
			const fvMesh& mesh,
			const vectorField& Fk, 
			vector repulsiveForce,
			point center,
			pointField lPoints,
			scalar rhoF,
			scalar dV,
			vector g,
			vector volIntegralU,
			vector volIntegralRxU
		);

		//- Return translational velocity
		virtual vector UTranslate();

		//- Return rotational velocity
		virtual vector URotate();

}; // End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //