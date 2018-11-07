/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBMotion class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef IBMotion_H
#define IBMotion_H

#include "autoPtr.H"
#include "IBObject.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        	Class  IBMotion Declaration
\*---------------------------------------------------------------------------*/

class IBMotion
{

		const word motionType_;

public:
	//- Runtime type information
	TypeName("IBMotion");

	//- Declare runtime constructor selection table
		declareRunTimeSelectionTable
		(
			autoPtr,
			IBMotion,
			dictionary,
			(
				const word& name,
				IBObject& obj,
				const dictionary& dict
			),
			(name, obj, dict)
		);

	//- Constructor
		IBMotion
		(
			const word& name,
			IBObject& obj,
			const dictionary& dict
		);

	//- Selectors
		static autoPtr<IBMotion> New 
		(
			const word& name,
			IBObject& obj,
			const dictionary& dict
		);

	//- Destructor
		virtual ~IBMotion()
		{}

	//- Member functions

		//- Return object's motion type
		word motionType()
		{
			return motionType_;
		}

		//- Return translational velocity
		virtual vector UTranslate() = 0; 
		
		//- Return rotational velocity
		virtual vector URotate() = 0;

}; // End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
