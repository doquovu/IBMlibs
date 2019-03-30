/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBModel class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "manageFlowRate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam 
{
    defineTypeNameAndDebug(manageFlowRate, 0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Foam::scalar Foam::manageFlowRate::initializeDP()
{

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::manageFlowRate::manageFlowRate
(
    IBdynamicFvMesh& mesh,
    const dictionary& dict 
)
:
    mesh_(mesh),
    dict_(dict),
    addPressureGradientField_(dict.lookup("addPressureGradField")),
    direction_(dict.lookup("direction")),
    type_(dict.lookup("type")),
    dP_(0.0),
    flowRate_(0.0)
{
    if (type_ == "fixedValue")
    {
        dP_ = readScalar(dict.lookup("value"));
    }
    else if (type_ = "fixedFlowRate")
    {
        flowRate_ = readScalar(dict.lookup("flowRate"));
        dP_ = initializeDP();
    }
    else 
        FatalErrorIn
		(
			"manageFlowRate::manageFlowRate(IBdynamicFvMesh&, "
			"const dictionary&)"
		)	<< "Unknown pressure gradient type "<< type_ 
			<< endl
			<< "Valid types are: "<<endl
			<< "(" <<endl<<"fixedValue"<<endl<<"calculated"<<endl<<")"<<endl
			<<exit(FatalError);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::manageFlowRate::correctPressureGradient()
{

}

Foam::volVectorField Foam::manageFlowRate::pressureField()
{
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

    if (addPressureGradientField_)
    {

        forAll(gradP, celli)
        {
            if (direction_ == "horizontal")
                gradP[celli] = vector(dP_, 0, 0);
            else if (direction_ == "vertical")
                gradP[celli] = vector(0, dP_, 0);
            else if (direction_ == "curvedChannel")
            {
                //- Assuming center of curve channel is (0 0 0)
                vector position = mesh_.C()[celli];
                vector flowDir = vector(-position.y(), position.x(), 0);
                flowDir /= mag(flowDir);
                gradP[celli] = dP_ * flowDir;
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
    }

  	return gradP;
}
