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
#include "triSurfaceSearch.H"

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
		fileName STLFile_;
		triSurfaceMesh triSurf_;

		word name_;
		word objectType_;
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
		
	//- Data for new ibm method: forcing cells
		
		//- Interface points
		pointField ifPts_;
		//- Virtual points
		pointField vtPts_;
		//- Fluid cells indicator
		volScalarField gamma_;

		labelList solidCells_;
		labelList ibCells_;
		labelListList ibNeiCells_;
		labelList vtPointCells_;

	//- Private functions
		scalar cellSize(label cellID);

		//- Read STL file
		void readSTL(const dictionary& dict);

		//- Create motion
		void addMotions(const dictionary& dict);

		//- 27-3: Vu, New approach
		void makeInterpolationStencil();

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

		const volScalarField& gamma() {return gamma_;}
		const pointField& interfacePoints() {return ifPts_;}
		const pointField& virtualPoints() {return vtPts_;}
		const labelList& solidCells() {return solidCells_;}
		const labelList& ibCells() {return ibCells_;}
		const labelListList& ibNeiCells() {return ibNeiCells_;}
		const labelList& virtualPointCells() {return vtPointCells_;}


		virtual word name();
		virtual word objectType();
		virtual PtrList<IBMotion>& motions();
		virtual scalar rho();
		virtual scalar V();
		virtual vector Ip();
		virtual point CG();
		virtual scalar R();
		virtual bool movable();
		virtual label nPoints();
		virtual pointField lPoints();
		virtual vector translationalVelocity();
		virtual vector rotationalVelocity();
		virtual label nFaces();
		virtual labelList nPointsOfFaces();
		virtual labelListList pointOfFace();

}; // End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //