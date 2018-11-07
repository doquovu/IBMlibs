/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - rotatingMotion class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef rotatingMotion_H
#define rotatingMotion_H

#include "IBMotion.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

/*---------------------------------------------------------------------------*\
						Class rotatingMotion declaration
\*---------------------------------------------------------------------------*/
class rotatingMotion
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
	TypeName("rotatingMotion");

	//- Constructor
		rotatingMotion
		(
			const word& name,
			IBObject& obj,
			const dictionary& dict
		);

	//- Destructor
		virtual ~rotatingMotion()
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
