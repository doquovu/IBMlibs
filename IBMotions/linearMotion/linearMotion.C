/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - linearMotion class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "linearMotion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearMotion, 0);
    addToRunTimeSelectionTable(IBMotion, linearMotion, dictionary);

//---------------------------------Constructors------------------------------//
linearMotion::linearMotion
(
	const word& name,
	IBObject& obj,
	const dictionary& dict
)
:
	IBMotion(name, obj, dict),
	UTranslate_(dict.lookup("u")),
	URotate_(vector::zero)
{}

//-------------------------------Member functions----------------------------//
vector linearMotion::UTranslate()
{
	return UTranslate_;
}

vector linearMotion::URotate()
{
	return URotate_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
