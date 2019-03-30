/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBTechnique class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "IBTechnique.H"
#include "pointMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(IBTechnique, 0);
	defineRunTimeSelectionTable(IBTechnique, dictionary);

//--------------------------Protected Members Functions-----------------------//
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

vector IBTechnique::linearInterpolate
(
    const point p1,
    const point p2,
    const point target,
    const vector v1,
    const vector v2
)
{
    vector targetValue(0,0,0);
    scalar h1 = mag(p1 - target);
    scalar h2 = mag(p2 - target);
    scalar h = h1 + h2;

    targetValue = h2/h*v1 + h1/h*v2;

    return targetValue;
}

vector IBTechnique::bilinearInterpolate
(
    const pointVectorField& phi,
    const point& target,
    const labelList& sourceLabel
)
{
    //- Assuming the interpolation is perfom in the XOY plane
    //- Check number of points
    if (sourceLabel.size() != 4)
    {
        FatalErrorIn
        (
            "directForcingBallaras::bilinearInterpolate()"
        )   << "Error: bilinearInterpolate needs 4 points, "
            << sourceLabel.size() <<" points provided"<<nl
            <<abort(FatalError);
    }
    //- Find intermediate points
    const pointField source(mesh_.points(), sourceLabel);

    Pair<label> pair1(0,0);
    Pair<label> pair2(0,0);

    pair1.first() = sourceLabel[0];

    for (int i=1; i<source.size(); i++)
    {
        if (source[i].y() == source[0].y())
        {
            pair1.second() = sourceLabel[i];

            if (i == 1)
            {
                pair2.first() = sourceLabel[2];
                pair2.second() = sourceLabel[3];
            }
            else if (i == 2)
            {
                pair2.first() = sourceLabel[1];
                pair2.second() = sourceLabel[3];
            }
            else 
            {
                pair2.first() = sourceLabel[1];
                pair2.second() = sourceLabel[2];
            }
            break;
        }
    }

    point p1(target.x(), mesh_.points()[pair1.first()].y(), target.z());
    point p2(target.x(), mesh_.points()[pair2.first()].y(), target.z());
    vector v1 
    = 
        linearInterpolate
        (
            mesh_.points()[pair1.first()],
            mesh_.points()[pair1.second()],
            p1,
            phi[pair1.first()],
            phi[pair1.second()]
        );
    vector v2
    = 
        linearInterpolate
        (
            mesh_.points()[pair2.first()],
            mesh_.points()[pair2.second()],
            p2,
            phi[pair2.first()],
            phi[pair2.second()]
        );

    return (linearInterpolate(p1, p2, target, v1, v2));
}

vector IBTechnique::trilinearInterpolate
(
    const pointVectorField& phi,
    const point& target,
    const labelList& sourceLabel
)
{
    //- Check number of points
    if (sourceLabel.size() != 8)
    {
        FatalErrorIn
        (
            "directForcingBallaras::trilinearInterpolate()"
        )   << "Error: trilinearInterpolate needs 8 points, "
            << sourceLabel.size() <<" points provided"<<nl
            <<abort(FatalError);
    }
    const pointField source(mesh_.points(), sourceLabel);
    labelList pointStage1;
    labelList pointStage2;
    
    pointStage1.append(sourceLabel[0]);

    for (int i=1; i<source.size(); i++)
    {
        if (source[i].z() == source[0].z())
        {
            pointStage1.append(sourceLabel[i]);
        }
        else 
        {
            pointStage2.append(sourceLabel[i]);
        }
    }

    point p1(target.x(), target.y(), mesh_.points()[pointStage1[0]].z());
    point p2(target.x(), target.y(), mesh_.points()[pointStage2[0]].z());
    vector v1 = bilinearInterpolate(phi, p1, pointStage1);
    vector v2 = bilinearInterpolate(phi, p2, pointStage2);
    
    return (linearInterpolate(p1, p2, target, v1, v2));
}
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
