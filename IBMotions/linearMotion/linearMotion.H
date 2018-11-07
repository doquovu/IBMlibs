/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - linearMotion class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef linearMotion_H
#define linearMotion_H

#include "IBMotion.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

/*---------------------------------------------------------------------------*\
						Class linearMotion declaration
\*---------------------------------------------------------------------------*/
class linearMotion
:
	public IBMotion
{
	//- Private data

		//- Translational velocity
		vector UTranslate_;

		//- Rotational velocity
		vector URotate_;

public:

	//- Runtime type information
	TypeName("linearMotion");

	//- Constructor
		linearMotion
		(
			const word& name,
			IBObject& obj,
			const dictionary& dict
		);

	//- Destructor
		virtual ~linearMotion()
		{}

	//- Member functions

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