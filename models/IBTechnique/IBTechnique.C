/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBTechnique class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "IBTechnique.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(IBTechnique, 0);
	defineRunTimeSelectionTable(IBTechnique, dictionary);


//---------------------------------Constructors------------------------------//

IBTechnique::IBTechnique
(
	IBdynamicFvMesh& mesh,
	const dictionary& dict
)
:	
	IBOStream(mesh),
	movingIBObjects(mesh,dict),
	mesh_(mesh),
	dict_(dict)
{}

// -------------------------------Member Functions----------------------------//
scalar IBTechnique::deltaFunc(point p_Eul, point p_Lag)
{
    scalar deltaX = deltaFunc1D(p_Eul.x(), p_Lag.x());
    scalar deltaY = deltaFunc1D(p_Eul.y(), p_Lag.y());
    scalar deltaZ = deltaFunc1D(p_Eul.z(), p_Lag.z());

    if (mesh_.nGeometricD() == 3)
        return (deltaX*deltaY*deltaZ);
    else
        return (deltaX*deltaY);
}

scalar IBTechnique::deltaFunc1D(scalar x_Eul, scalar x_Lag)
{
    const scalar r = mag((x_Eul - x_Lag) / h());
    scalar phiR = 0.0;
    if (r < 0.5)
    {
        phiR = 1.0/3.0 * (1.0 + Foam::sqrt(-3.0 * (r*r) + 1.0));
    }
    else if (r <= 1.5)
    {
        phiR = 1.0/6.0 * (5.0 - 3.0*r - 
            Foam::sqrt(-3.0* (1.0-r)*(1.0-r) + 1.0));    
    }
    else
    {
        phiR = 0.0;
    }

    scalar delta1D = phiR / h();

    return delta1D;
}

volVectorField IBTechnique::pressureGradField()
{
	scalar dP = dict_.lookupOrDefault("gradPField", 0.0);
	word dPDir = dict_.lookup("gradPDir");
	Info<<"IBM: Adding pressure gradient field: "<< dP <<nl<<endl;
	volVectorField gradP
	(
	    IOobject
	    (
	        "gradP",
	        mesh_.time().timeName(),
	        mesh_,
	        IOobject::NO_READ,
	        IOobject::AUTO_WRITE
	    ),
	    mesh_,
	    dimensionedVector("gradP", dimensionSet(0,1,-2,0,0,0,0), vector::zero)
	);
	forAll(gradP, celli)
	{
		if (dPDir == "horizontal")
			gradP[celli] = vector(dP, 0, 0);
		else if (dPDir == "vertical")
			gradP[celli] = vector(0, dP, 0);
		else if (dPDir == "curvedChannel")
		{
			//- Assuming center of curve channel is (0 0 0)
			vector position = mesh_.C()[celli];
			vector flowDir = vector(-position.y(), position.x(), 0);
			flowDir /= mag(flowDir);
			gradP[celli] = dP * flowDir;
		}
		else
		FatalErrorIn
		(
			"IBTechnique::pressureGradField()"
		)	<< "Unknown direction of external pressure gradient field "<< dPDir <<endl
			<< "Valid approaching method are : " <<endl
			<< "(" <<endl
			<< "horizontal" <<endl 
			<< "vertical" <<endl
			<< "curvedChannel"<<endl
			<<")"<<endl
			<< exit(FatalError);
	}

	return gradP;
}

autoPtr<IBTechnique> IBTechnique::New
(
	IBdynamicFvMesh& mesh,
	const dictionary& dict
)
{
	word typeName = dict.lookup("IBtechnique");
	
	Info<< "Selecting Immersed Boundary approaching method: "<<typeName<< endl;
	
	dictionaryConstructorTable::iterator cstrIter = 
		dictionaryConstructorTablePtr_->find(typeName);
		
	if (cstrIter == dictionaryConstructorTablePtr_->end())
	{
		FatalErrorIn
		(
			"IBTechnique::New(const word& name, const IBdynamicFvMesh& mesh,"
			"const dictionary& dict)"
		)	<< "Unknown approaching method "<< typeName <<endl <<endl
			<< "Valid approaching method are : " <<endl
			<< dictionaryConstructorTablePtr_->toc()
			<< exit(FatalError);
	}
	
	return autoPtr<IBTechnique>(cstrIter()(mesh, dict));
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} //- End namespace Foam

// ************************************************************************* //
