/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - rotatingMotion class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "rotatingMotion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rotatingMotion, 0);
    addToRunTimeSelectionTable(IBMotion, rotatingMotion, dictionary);

//---------------------------------Constructors------------------------------//
rotatingMotion::rotatingMotion
(
	const word& name,
	IBObject& obj,
	const dictionary& dict
)
:
	IBMotion(name,obj, dict),
	UTranslate_(vector::zero),
	URotate_(dict.lookup("omega"))
{}

//-------------------------------Member functions----------------------------//
vector rotatingMotion::UTranslate()
{
	return UTranslate_;
}

vector rotatingMotion::URotate()
{
	return URotate_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
