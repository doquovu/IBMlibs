/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBObject class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/
#include "error.H"
#include "IBObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IBObject, 0);
    defineRunTimeSelectionTable(IBObject, dictionary);


//---------------------------------Constructors------------------------------//
IBObject::IBObject
(
	const word& name,
	const fvMesh& mesh,
	const dictionary& dict
)
{}

//------------------------------------Selectors------------------------------//

autoPtr<IBObject> IBObject::New
(
	const word& name,
	const fvMesh& mesh,
	const dictionary& dict
)
{
	word typeName = dict.lookup("type");
	dictionaryConstructorTable::iterator dictCstrIter = 
		dictionaryConstructorTablePtr_->find(typeName);

	if(dictCstrIter == dictionaryConstructorTablePtr_->end())
	{
		FatalErrorIn
		(
			"IBObject::New(const fvMesh& mesh, "
			"const word& name, const dictionary& dict)"
		)	<< "Unknown object type "<< typeName 
			<< endl
			<< "Valid object type are: "<<endl
			<<dictionaryConstructorTablePtr_->toc()
			<<exit(FatalError);
	}

	return autoPtr<IBObject>(dictCstrIter()(name, mesh, dict));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} //- End namespace Foam

// ************************************************************************* //
