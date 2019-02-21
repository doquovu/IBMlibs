/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBSTL class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef IBSTL_H
#define IBSTL_H

#include "fvCFD.H"
#include "IBObject.H"
#include "IBMotion.H"
#include "triSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
							Class IBSTL declaration
\*---------------------------------------------------------------------------*/

class IBSTL
:
	public IBObject
{

		const fvMesh& mesh_;

		word name_;
		
		word objectType_;

		point max_;

		point min_;

		point center_;

		scalar R_;

		vector Ip_;
		
		scalar rho_;

		scalar V_;

		label nPoints_;

		pointField lPoints_;

		bool movable_;

		PtrList<IBMotion> motions_;

		label nFaces_;
		
		labelList nPointsOfFaces_;
		
		labelListList pointOfFace_;

	//- Private functions
		
		//- Read STL file
		void readSTL(const dictionary& dict);

		//- Create motion
		void addMotions(const dictionary& dict);

public:

	//- Runtime type information
	TypeName("IBSTL");

	//- Constructors
		IBSTL
		(	
			const word& name,
			const fvMesh& mesh,
			const dictionary& dict
		);
	
	//- Destructor
		virtual ~IBSTL()
		{};

	//- Member functions
		virtual word name();

		virtual word objectType();

		virtual PtrList<IBMotion>& motions();
		
		virtual scalar rho();

		virtual scalar V();

		virtual vector Ip();

		virtual point CG();

		virtual scalar R();

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