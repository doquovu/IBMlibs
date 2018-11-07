/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - movingIBObjects class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef IBMotionRegistry_H
#define IBMotionRegistry_H

#include "IOmanip.H"
#include "autoPtr.H"
#include "IBMotion.H"
#include "IBObjectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
							Class movingIBObjects declaration
\*---------------------------------------------------------------------------*/

class movingIBObjects
:
	public IBObjectRegistry
{

	const fvMesh& mesh_;
	
	//- Velocity at Lagrangian points on object's surface
	List<vectorField> uBoundary_;

	bool calcRepulsive_;

	scalar epsilonP_;

	scalar epsilonW_;

public:
	
	TypeName("movingIBObjects");
	
	movingIBObjects(const fvMesh& mesh, const dictionary& dict);

	~movingIBObjects()
	{}
	
	void calcUBoundary(); 

	void movePoints();

	void updateObjectMotions
	(
		label objID,
		const vectorField& ForceLagrange
	);

	vector volIntegralU
	(
	    label objID,
	    const volVectorField& U
	);
	
	vector volIntegralRxU
	(
	    label objID,
	    const volVectorField& U
	);

	List<vector> calcRepulsiveForce();

	vector objMutualRepulsive(label obj1ID, label obj2ID);

	vector objWallRepulsive(label objID);

	List<vectorField>& UBoundary()
	{
		return uBoundary_;
	}

}; //End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //