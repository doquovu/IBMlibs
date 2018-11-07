/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - sixDoFMotion class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "sixDoFMotion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sixDoFMotion, 0);
    addToRunTimeSelectionTable(IBMotion, sixDoFMotion, dictionary);

//--------------------------------Private Members----------------------------//
void sixDoFMotion::initialise()
{
	UTranslate_ = vector::zero;
	URotate_ = vector::zero;
}

//---------------------------------Constructors------------------------------//
sixDoFMotion::sixDoFMotion
(
	const word& name,
	IBObject& obj,
	const dictionary& dict
)
:
	IBMotion(name, obj, dict),
	obj_(obj),
	UTranslate_(),
	URotate_()
{
	initialise();
}
//-------------------------------Member functions----------------------------//
void sixDoFMotion::updateMotion
(
	vector& uTransl,
	vector& uRotate,
	const fvMesh& mesh,
	const vectorField& Fk, 
	vector repulsiveForce,
	point center,
	pointField lPoints,
	scalar rhoF,
	scalar dV,
	vector g
)
{
	//- M. Uhlmann 2005
	const scalar dT = mesh.time().deltaTValue();  

	vector force_ = vector::zero;
	forAll(Fk, i)
	    force_ += Fk[i]*dV;

	vector uTranslNew = uTransl + dT*( - rhoF*force_/(obj_.V()*(obj_.rho()-rhoF)) 
	                		 + repulsiveForce /(obj_.V()*(obj_.rho() - rhoF))
							 + g
	                		); 
	uTransl = 0.5*(uTranslNew+uTransl);

	vector f = vector::zero;
	forAll(lPoints, pointI)
	{
		vector dr = lPoints[pointI] - center;
	    f += (dr ^ Fk[pointI])*dV;  
	}

	vector uRotateNew = uRotate - dT*obj_.rho()*rhoF*f/(obj_.Ip()*(obj_.rho()-rhoF));
	uRotate = 0.5*(uRotateNew+uRotate);
}

void sixDoFMotion::updateMotion
(
	vector& uTransl,
	vector& uRotate,
	const fvMesh& mesh,
	const vectorField& Fk, 
	vector repulsiveForce,
	point center,
	pointField lPoints,
	scalar rhoF,
	scalar dV,
	vector g,
	vector volIntegralU,
	vector volIntegralRxU
)
{
	//- Tobias Kempe & J.Frohlich 2012
	const scalar dT = mesh.time().deltaTValue();  

	vector f = vector::zero;
	forAll(Fk, i)
	    f += Fk[i]*dV;

	//vector uTranslNew =
	uTransl = 
	  uTransl + dT/(obj_.V()*obj_.rho())*( - rhoF*f 
	        					   + rhoF*volIntegralU
								   + obj_.V()*(obj_.rho()-rhoF)*g
								   + repulsiveForce ); 
	
	//uTransl = 0.5*(uTranslNew+uTransl);

	vector fxR = vector::zero;
	forAll(lPoints, pointI)
	{
		vector r = lPoints[pointI] - center;
	    fxR += (r ^ Fk[pointI])*dV;  
	}

	//vector uRotateNew = 
		uRotate = uRotate + dT*rhoF/obj_.Ip()*( -fxR + volIntegralRxU );

	//uRotate = 0.5*(uRotateNew+uRotate);

	if (mesh.nGeometricD() == 2) 
	{
		uTransl.z() = 0.0;
		uRotate.x() = 0.0;
		uRotate.y() = 0.0;
	}
}

vector sixDoFMotion::UTranslate()
{
	return UTranslate_;
}

vector sixDoFMotion::URotate()
{
	return URotate_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
