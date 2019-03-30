/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - directForcingBallaras class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef directForcingBallaras_H
#define directForcingBallaras_H

#include "IBTechnique.H"
#include "IBObject.H"
#include "IBObjectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
							Class IBModel declaration
\*---------------------------------------------------------------------------*/

class directForcingBallaras
: 
	public IBTechnique
{
	//- Private data

		//- Reference to fvMesh
		IBdynamicFvMesh& mesh_;
		
		//- Number of multi-direct forcing iteration
		label nMDF_;
		
		List<vectorField> Fk_;

	vector linearInterpolate
	(
		const point p1,
		const point p2,
		const point target,
		const vector v1,
		const vector v2
	);

	vector bilinearInterpolate
	(
		const pointVectorField& phi,
		const point& target,
		const labelList& sourceLabel
	);

	vector trilinearInterpolate
	(
		const pointVectorField& phi,
		const point& target,
		const labelList& sourceLabel
	);
	
	void makeIbForce
	(
		volVectorField& U, 
		volVectorField& f
	);

	void makeIbForceNew
	(
		volVectorField& U, 
		volVectorField& f
	);

	//- Disallow default bitwise copy construct
	directForcingBallaras(const directForcingBallaras&);
	
	//- Disallow default bitwis assignment
	void operator=(const directForcingBallaras&);

	
public:

	//- Runtime type information
	TypeName("directForcingBallaras");

	//- Constructor
		directForcingBallaras
		(
			// const word& name,
			IBdynamicFvMesh& mesh,
			const dictionary& dict
		);

	//- Destructor
		virtual ~directForcingBallaras()
		{}

	//- Member functions

		inline List<vectorField> Fk()
		{
			return Fk_;
		}

		//- Return force field for momentum equation
		virtual volVectorField ibForce(volVectorField& U);
		
		virtual volVectorField ibForceInt();

		// 12/3/2019 
		virtual volVectorField ibForceInt(const volVectorField& rhs);
		
		//- Perform multi-Firect forcing
		virtual void multiDirectForcing
		(
			volVectorField& u,
			volVectorField& ibForce
		);
		
		//- Update object's velocity and position 
		virtual void update();

		//- Write data for post-processing
		virtual void write();

}; //End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
