/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBWall class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef IBWall_H
#define IBWall_H

#include "fvCFD.H"
#include "IBObject.H"
#include "IBMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
							Class IBWall declaration
\*---------------------------------------------------------------------------*/

class IBWall
:
	public IBObject
{

		const fvMesh& mesh_;

	//- Particle parameters
		word name_;
		
		word objectType_;

		point max_;

		point min_;

		label nPoints_;

		pointField lPoints_;

		PtrList<IBMotion> motions_;

	//- Private functions
		
		//- Create 2D particles automatically
		void createWall2D(const dictionary& dict);
		
public:

	//- Runtime type information
	TypeName("IBWall");

	//- Constructors
		IBWall
		(	
			const word& name,
			const fvMesh& mesh,
			const dictionary& dict
		);
	
	//- Destructor
		virtual ~IBWall()
		{};

	//- Member functions
		point min() {return min_;}
		point max() {return max_;}

		virtual word name();

		virtual word objectType();

		virtual PtrList<IBMotion>& motions();
		
		virtual scalar rho();

		virtual scalar V();

		virtual vector Ip();

		virtual scalar R();

		virtual point CG();
		
		virtual bool movable();

		//- Return number of Lagrangian points
		virtual label nPoints();

		//- Return Lagrang points
		virtual pointField lPoints();

		//- Return translational velocity
		virtual vector translationalVelocity();

		//- Return rotational velocity
		virtual vector rotationalVelocity();

		//- Return number of faces to write object
		virtual label nFaces();

		//- Return number of points of each object's face
		virtual labelList nPointsOfFaces();

		//- Return coordinate of points of each object's face
		virtual labelListList pointOfFace();

}; // End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //