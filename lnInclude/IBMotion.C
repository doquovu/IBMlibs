/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBMotion class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "IBMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IBMotion, 0);
    defineRunTimeSelectionTable(IBMotion, dictionary);


//---------------------------------Constructors------------------------------//
IBMotion::IBMotion
(
	const word& name,
	IBObject& obj,
	const dictionary& dict
)
:
	motionType_(name)
{}

//------------------------------------Selectors------------------------------//

autoPtr<IBMotion> IBMotion::New
(
	const word& name,
	IBObject& obj,
	const dictionary& dict
)
{
	word typeName = name;//dict.lookup("type");
	dictionaryConstructorTable::iterator dictCstrIter = 
		dictionaryConstructorTablePtr_->find(typeName);

	if(dictCstrIter == dictionaryConstructorTablePtr_->end())
	{
		FatalErrorIn
		(
			"IBMotion::New(const fvMesh& mesh, "
			"const word& name, const dictionary& dict)"
		)	<< "Unknown object type "<< typeName 
			<< endl
			<< "Valid object type are: "<<endl
			<<dictionaryConstructorTablePtr_->toc()
			<<exit(FatalError);
	}

	return autoPtr<IBMotion>(dictCstrIter()(name, obj, dict));
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} //- End namespace Foam

// ************************************************************************* //
