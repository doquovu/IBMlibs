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
	const fvMesh& mesh,
	const dictionary& dict
)
:	
	IBOStream(mesh),
	movingIBObjects(mesh,dict),
	mesh_(mesh),
	pressureGradField_(dict.lookup("gradP"))
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
	Info<<"IBM: Adding pressure gradient field:"<<pressureGradField_<<nl<<endl;
	volVectorField gradP
	(
	    IOobject
	    (
	        "gradP",
	        mesh_.time().timeName(),
	        mesh_,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
	    mesh_,
	    dimensionedVector("gradP", dimensionSet(0,1,-2,0,0,0,0), pressureGradField_)
	);

	return gradP;
}

autoPtr<IBTechnique> IBTechnique::New
(
	const fvMesh& mesh,
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
			"IBTechnique::New(const word& name, const fvMesh& mesh,"
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
