/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBObject class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef IBObject_H
#define IBObject_H

#include "point.H"
#include "fvMesh.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Class forward declaration
class IBMotion;

/*---------------------------------------------------------------------------*\
                        Class  IBObjectSelection Declaration
\*---------------------------------------------------------------------------*/

class IBObject
{

public:
	//- Runtime type information
	TypeName("IBObject");

	//- Declare runtime constructor selection table
		declareRunTimeSelectionTable
		(
			autoPtr,
			IBObject,
			dictionary,
			(
				const word& name,
				const fvMesh& mesh,
				const dictionary& dict
			),
			(name, mesh, dict)
		);

	//- Constructor
		IBObject
		(
			const word& name,
			const fvMesh& mesh,
			const dictionary& dict
		);

	//- Selectors
		static autoPtr<IBObject> New 
		(
			const word& name,
			const fvMesh& mesh,
			const dictionary& dict
		);

	//- Destructor
		virtual ~IBObject()
		{}

	//- Member functions

		//- Access

			//- Return object's name
			virtual word name() = 0;

			virtual word objectType() = 0; 

			virtual PtrList<IBMotion>& motions() = 0;

			//- Return object's density
			virtual scalar rho() = 0;
			
			//- Return object's volume
			virtual scalar V() = 0;
			
			//- Return Moment of Inertia
			virtual vector Ip() = 0;
			
			//- Return the radious that enclosed objects
			virtual scalar R() = 0;

			//- Return the center of gravity
			virtual point CG() = 0;
		
			//- Return true if object moves
			virtual bool movable() = 0;
			
			//- Return number of Lagrangian points
			virtual label nPoints() = 0;

			//- Return Lagrang points
			virtual pointField lPoints() = 0;

			//- Return translational velocity
			virtual vector translationalVelocity() = 0;

			//- Return rotational velocity
			virtual vector rotationalVelocity() = 0;

			//- Return number of faces to write object
			virtual label nFaces() = 0;

			//- Return number of points of each object's face
			virtual labelList nPointsOfFaces() = 0;

			//- Return coordinate of points of each object's face
			virtual labelListList pointOfFace() = 0;

}; // End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
