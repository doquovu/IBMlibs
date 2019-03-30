/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBParticles class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef IBParticles_H
#define IBParticles_H

#include "fvCFD.H"
#include "IBObject.H"
#include "IBMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
							Class IBParticles declaration
\*---------------------------------------------------------------------------*/

class IBParticles
:
	public IBObject
{

		const scalar PI;

		const fvMesh& mesh_;

	//- Particle parameters
		word name_;
		word objectType_;

		//- Center
		point center_;
		//- Radius
		scalar R_;
		//- Density 
		scalar rho_;
		//- Particle volume
		scalar V_;
		//- Moment of Inertia
		vector Ip_;
		//- Number of Lagrang points
		label nPoints_;
		//- Lagrang points
		pointField lPoints_;
		//- Do object move?
		bool movable_;
		//- Pointer to object's motion
		PtrList<IBMotion> motions_;

	//- Points connectivity for writing geometry
		label nFaces_;
		labelList nPointsOfFaces_;
		labelListList pointOfFace_;

	//- Private functions
		
		//- Distribute evenly Lagrangian points on sphere surface
		pointField createUnitSphereEqAreaPartition(label nPoints);
		scalar capArea(scalar angle);

		//- Create 2D particles automatically
		void createParticle2D(const dictionary& dict);
		
		//- Create 3D particles from input TXT file
		void createParticle3D(const dictionary& dict);

		//- Create motion
		void addMotions(const dictionary& dict);

public:

	//- Runtime type information
	TypeName("IBParticles");

	//- Constructors
		IBParticles
		(	
			const word& name,
			const fvMesh& mesh,
			const dictionary& dict
		);
	
	//- Destructor
		virtual ~IBParticles()
		{};

	//- Member functions
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