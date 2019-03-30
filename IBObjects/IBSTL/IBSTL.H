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
		

		//- Data for new ibm method: Yang, Ballaras
			
			labelList solidCells_;
			labelList ibCells_;
			labelListList ibNeiCells_;
			labelListList vtNeiCells_;
			labelList vtPointCells_;

			//- forcing points
				labelList ibPts_;
			//- Interface points
				pointField ifPts_;
			//- Virtual points
				pointField vtPts_;
			//- neighbour points of virtual points
				labelListList neiPts_;
			//- solidPoints
				labelList sldPts_;
			//- Fluid cells indicator
				volScalarField gamma_;

	//- Private functions
		scalar cellSize(label cellID);

		//- Read STL file
		void readSTL(const dictionary& dict);

		//- Create motion
		void addMotions(const dictionary& dict);

		//- Interpolation stencil
		void makeInterpolationStencil();
		
		//- 27-3: Vu, New approach
		void makeInterpolationStencilNew();

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

		const labelList& solidCells() {return solidCells_;}
		const labelList& ibCells() {return ibCells_;}
		const labelListList& ibNeiCells() {return ibNeiCells_;}
		const labelListList& vtNeiCells() {return vtNeiCells_;}
		const labelList& virtualPointCells() {return vtPointCells_;}

		const labelList& forcingPoints() {return ibPts_;}
		const pointField& interfacePoints() {return ifPts_;}
		const pointField& virtualPoints() {return vtPts_;}
		const labelListList& neighbourPoints() {return neiPts_;}
		const labelList& solidPoints() {return sldPts_;}
		const volScalarField& gamma() {return gamma_;}

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