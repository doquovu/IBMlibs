/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBObjectRegistry class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef IBObjectRegistry_H
#define IBObjectRegistry_H

#include "fvMesh.H"
#include "volFields.H"
#include "IBObject.H"
#include "pointField.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

/*---------------------------------------------------------------------------*\
						Class IBOBjectRegistry declaration
\*---------------------------------------------------------------------------*/
class IBObjectRegistry
{
	//- Domain data

		const fvMesh& mesh_;

		const IOdictionary transportProperties_;

		//- Cartesian grid size
		scalar h_;

		//- Cartesian volume
		scalar dV_;

		//- Density of fluid
		scalar rhoF_;

		//- Gravitational acceleration
		vector g_;

	//- Objects data
		
		label nObjects_;

		PtrList<IBObject> objects_;

		List<pointField> LPoints_;

		List<scalar> R_;

		List<point> CoG_;

		List<labelListList> neiCells_;
		
		List<labelList> solidCells_;

		List<labelList> solidCellsExt_;

		List<vector> UTranslate_;

		List<vector> URotate_;

		List<label> nFaces_;

		List<labelList> nPointsOfFaces_;

		List<labelListList> pointOfFace_;

	//- Shadows data

		bool periodicBC_;

		List<bool> enableShadows_;

		vector cyclicDistance_;

		List<pointField> Shd1LPoints_;
		List<pointField> Shd2LPoints_;

		List<point> Shd1CoG_;
		List<point> Shd2CoG_;

		List<labelListList> Shd1NeiCells_;
		List<labelListList> Shd2NeiCells_;
		List<labelList> Shd1SolidCells_;
		List<labelList> Shd2SolidCells_;

		scalar bendingAngle_;
		scalar innerRadius_;
	//- Data for AMR
		labelList cellsToRefine_;
		label maxRefinement_;

protected:

		scalar cellSize(label cellID);
		
		vector rotate(vector v, point root, scalar angle);
		point rotate(point p, scalar angle);
		scalar getAngle(point p);

		labelListList findNeiCells(const pointField& points);
		
		labelList findSolidCells(label objectID, point C);

		void readEnvironmentVariables();

		bool isInsideRegion1(const label objectID, const point p);
		
		bool isInsideRegion2(const label objectID, const point p);

		bool isInsideRegion3(const label objectID, const point p);
		
		void cartesianGridSize();
		
		void detectCyclicBoundary();
		
		void createObjects(const dictionary& dict);

		void createShadows();

		void swap(point& p1, point& p2);
		
		void getObjectsData();

		template<class T>
    	List<T> make1DList(const List<List<T>>& twoDList);

public:

	TypeName("IBObjectRegistry");

	//- Constructor
		IBObjectRegistry
		(
			const fvMesh& mesh,
			const dictionary& dict
		);

	//- Destructor
		virtual ~IBObjectRegistry()
		{}

	//- Member functions

		inline scalar bendingAngle();
		inline scalar h();
		inline scalar dV();
		inline scalar rhoF();
		inline vector g();
		inline label nObjects();
		inline PtrList<IBObject>& objects();
		inline List<pointField>& LPoints();
		inline List<pointField>& Shd1LPoints();
		inline List<pointField>& Shd2LPoints();
		inline List<point>& CoG();
		inline List<point>& Shd1CoG();
		inline List<point>& Shd2CoG();
		inline List<labelListList>& neiCells();
		inline List<labelList>& solidCells();
		inline List<labelListList>& Shd1NeiCells();
		inline List<labelListList>& Shd2NeiCells();
		inline List<labelList>& Shd1SolidCells();
		inline List<labelList>& Shd2SolidCells();
		inline List<vector>& UTranslate();
		inline List<vector>& URotate();
		inline List<bool> enableShadows();
		inline vector cyclicDistance();
		inline List<label>& nFaces();
		inline List<labelList>& nPointsOfFaces();
		inline List<labelListList>& pointOfFace();
		inline scalar HeavisideFunc(scalar phi);
		inline scalar LevelSetFunc(	label objID, point p, point c);
		inline const labelList& cellsToRefine();

		scalar volFraction(label objID, label cellID, point c);
		void writeNeighbourCells();
		void updateCartesianGridSize();

}; // End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IBObjectRegistryI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //