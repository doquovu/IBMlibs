/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBBox class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef IBBox_H
#define IBBox_H

#include "fvCFD.H"
#include "IBObject.H"
#include "IBMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
							Class IBBox declaration
\*---------------------------------------------------------------------------*/

class IBBox
:
	public IBObject
{

		const fvMesh& mesh_;

	//- Particle parameters
		word name_;
		
		word objectType_;

		point max_;

		point min_;

		scalar R_;

		scalar rho_;

		scalar V_;

		scalar Ipx_;

		scalar Ipy_;

		scalar Ipz_;

		label nPoints_;

		pointField lPoints_;

		bool movable_;

		PtrList<IBMotion> motions_;

		label nFaces_;
		
		labelList nPointsOfFaces_;
		
		labelListList pointOfFace_;

	//- Private functions
		
		//- Read input TXT file for 3D particles
		void readTXTFile(word fileName, point center, scalar R);
		
		//- Create 2D particles automatically
		void createBox2D(const dictionary& dict);
		
		//- Create 3D particles from input TXT file
		void createBox3D(const dictionary& dict);

		//- Create motion
		void addMotions(const dictionary& dict);

public:

	//- Runtime type information
	TypeName("IBBox");

	//- Constructors
		IBBox
		(	
			const word& name,
			const fvMesh& mesh,
			const dictionary& dict
		);
	
	//- Destructor
		virtual ~IBBox()
		{};

	//- Member functions
		virtual word name();

		virtual word objectType();

		virtual PtrList<IBMotion>& motions();
		
		virtual scalar rho();

		virtual scalar V();

		virtual scalar Ip();

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