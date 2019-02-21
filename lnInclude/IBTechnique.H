/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBTechnique class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef IBTechnique_H
#define IBTechnique_H

#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "volFields.H"
#include "IBOStream.H"
#include "IBObjectRegistry.H"
#include "movingIBObjects.H"
#include "dictionary.H"
#include "IBdynamicFvMesh.H"
// #include "IBAdaptiveRefineFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
							Class IBModel declaration
\*---------------------------------------------------------------------------*/


class IBTechnique
:
	public IBOStream,
	public movingIBObjects
	// public IBAdaptiveRefineFvMesh
{
		
	//- Private data
	IBdynamicFvMesh& mesh_;
	
	const dictionary& dict_;

	//- Disallow default bitwise copy construct
	IBTechnique(const IBTechnique&);
	
	//- Disallow default bitwise assignment
	void operator=(const IBTechnique&);
		
public:		
	//- Runtime type information
	TypeName("IBTechnique");

		declareRunTimeSelectionTable
		(
			autoPtr,
			IBTechnique,
			dictionary,
			(
				IBdynamicFvMesh& mesh,
				const dictionary& dict
			),
			(mesh, dict)
		);
		
	//- Constructor
		IBTechnique
		(
			IBdynamicFvMesh& mesh,
			const dictionary& dict
		);
		
	//- Selector
	static autoPtr<IBTechnique> New
	(
		IBdynamicFvMesh& mesh,
		const dictionary& dict
	);
	//- Destructor
		virtual ~IBTechnique()
		{}

	//- Member functions
		//- Return delta function for interpolation
		scalar deltaFunc(point pE, point pL);

		//- Return 1D delta function
		scalar deltaFunc1D(scalar xE, scalar xL);

		//- Return pressure gradient field for periodic problems
		volVectorField pressureGradField();

		virtual volVectorField ibForce(volVectorField& U) = 0;

		virtual volVectorField ibForceInt() = 0;

		//- Perform multi-Firect forcing
		virtual void multiDirectForcing
		(
			volVectorField& u,
			volVectorField& ibForce
		) = 0;

		virtual void update() = 0;

		virtual void write() = 0;

}; //End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
