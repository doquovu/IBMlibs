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

void sixDoFMotion::readConstraint(const dictionary& dict)
{
	if (dict.found("constraint"))
	{
		constraintDir_ = dict.subDict("constraint").lookup("axis");
		Info<<"FOUND CONSTRAINT MOTION: "<<constraintDir_<<endl;
	}
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
	solver_(dict.lookup("solver")),
	UTranslate_(),
	URotate_(),
	constraintDir_(vector::zero)
{
	initialise();
	readConstraint(dict);
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
	Info<<"IBM: Calculating volume integral following"
                        <<" M. Ulhmann assumption"<<endl;
	const scalar dT = mesh.time().deltaTValue();  

	vector force_ = vector::zero;
	forAll(Fk, i)
	    force_ += Fk[i]*dV;

	uTransl
	= 
		  uTransl 
		+ dT*( - rhoF*force_/(obj_.V()*(obj_.rho()-rhoF)) 
		+ repulsiveForce /(obj_.V()*(obj_.rho() - rhoF))
		+ g );

	vector fxR = vector::zero;
	forAll(lPoints, pointI)
	{
		vector dr = lPoints[pointI] - center;
	    fxR += (dr ^ Fk[pointI])*dV;  
	}

	if (mesh.nGeometricD() == 2) 
	{
		if (mesh.geometricD().x() < 0)
		{
			uTransl.x() = 0.0;
			uRotate.y() = 0.0;
			uRotate.z() = 0.0;
			uRotate.x() 
			= 
				uRotate.x() 
			  - dT*obj_.rho()*rhoF*fxR.x()/(obj_.Ip().x()*(obj_.rho()-rhoF));
		}
		else if (mesh.geometricD().y() < 0)
		{
			uTransl.y() = 0.0;
			uRotate.z() = 0.0;
			uRotate.x() = 0.0;
			uRotate.y() 
			= 
				uRotate.y() 
			  - dT*obj_.rho()*rhoF*fxR.y()/(obj_.Ip().y()*(obj_.rho()-rhoF));
		}
		else 
		{
			uTransl.z() = 0.0;
			uRotate.x() = 0.0;
			uRotate.y() = 0.0;
			uRotate.z() 
			= 
				uRotate.z() 
			  - dT*obj_.rho()*rhoF*fxR.z()/(obj_.Ip().z()*(obj_.rho()-rhoF));
		}
	}
	else if (mesh.nGeometricD() == 3)
	{
		uRotate.x() 
		= 
			uRotate.x() 
		  - dT*obj_.rho()*rhoF*fxR.x()/(obj_.Ip().x()*(obj_.rho()-rhoF));
		uRotate.y() 
		= 
			uRotate.y() 
		  - dT*obj_.rho()*rhoF*fxR.y()/(obj_.Ip().y()*(obj_.rho()-rhoF));
		uRotate.z() 
		= 
			uRotate.z() 
		  - dT*obj_.rho()*rhoF*fxR.z()/(obj_.Ip().z()*(obj_.rho()-rhoF));
	}
    else
    {
    	FatalErrorIn
    	(
    		"sixDoFMotion::updateMotion(vector& uTransl, "
			"vector& uRotate, const fvMesh& mesh," 
			"const vectorField& Fk, vector repulsiveForce," 
			"point center, pointField lPoints, scalar rhoF, "
			"scalar dV, vector g)"
    	)	<< "The mesh is neither 2D nor 3D " << endl
    		<< "Please check again"<<endl
    		<<exit(FatalError);
    }

	if (constraintDir_ != vector::zero)
	{
		uTransl = (uTransl & constraintDir_)*constraintDir_;
	}
	
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
	Info<<"IBM: Calculating volume integral following"
                        <<" Tobias K. formula"<<endl;

	const scalar dT = mesh.time().deltaTValue();  

	vector force_ = vector::zero;
	forAll(Fk, i)
	    force_ += Fk[i]*dV;

	uTransl 
	=
	  	uTransl 
	  + dT/(obj_.V()*obj_.rho())
	  * ( 
	  		- rhoF*force_ 
	  		+ rhoF*volIntegralU 
	  		+ obj_.V()*(obj_.rho()-rhoF)*g 
	  		+ repulsiveForce 
	  	);
	
	vector fxR = vector::zero;
	forAll(lPoints, pointI)
	{
		vector r = lPoints[pointI] - center;
	    fxR += (r ^ Fk[pointI])*dV;  
	}

    if (mesh.nGeometricD() == 2) 
    {
    	if (mesh.geometricD().x() < 0)
    	{
    		uTransl.x() = 0.0;
    		uRotate.y() = 0.0;
    		uRotate.z() = 0.0;
    		uRotate.x() 
    		= 
    			uRotate.x() 
    		  + dT*rhoF/obj_.Ip().x()*( -fxR.x() + volIntegralRxU.x());
    	}
    	else if (mesh.geometricD().y() < 0)
    	{
    		uTransl.y() = 0.0;
    		uRotate.z() = 0.0;
    		uRotate.x() = 0.0;
    		uRotate.y() 
    		= 
    			uRotate.y() 
    		  + dT*rhoF/obj_.Ip().y()*( -fxR.y() + volIntegralRxU.y());
    	}
    	else 
    	{
    		uTransl.z() = 0.0;
    		uRotate.x() = 0.0;
    		uRotate.y() = 0.0;
    		uRotate.z() 
    		= 
    			uRotate.z() 
    		  + dT*rhoF/obj_.Ip().z()*( -fxR.z() + volIntegralRxU.z());
    	}
    }
    else if (mesh.nGeometricD() == 3)
    {
    	uRotate.x()
    	= 
    		uRotate.x() 
    	  + dT*rhoF/obj_.Ip().x()*( -fxR.x() + volIntegralRxU.x() );
    	uRotate.y() 
    	= 
    		uRotate.y() 
    	  + dT*rhoF/obj_.Ip().y()*( -fxR.y() + volIntegralRxU.y() );
    	uRotate.z() 
    	= 
    		uRotate.z() 
    	  + dT*rhoF/obj_.Ip().z()*( -fxR.z() + volIntegralRxU.z() );
    }
    else
    {
    	FatalErrorIn
    	(
    		"sixDoFMotion::updateMotion(vector& uTransl, "
			"vector& uRotate, const fvMesh& mesh," 
			"const vectorField& Fk, vector repulsiveForce," 
			"point center, pointField lPoints, scalar rhoF, "
			"scalar dV, vector g, vector volIntegralU, "
			"vector volIntegralRxU )"
    	)	<< "The mesh is neither 2D nor 3D " << endl
    		<< "Please check again"<<endl
    		<<exit(FatalError);
    }

	if (constraintDir_ != vector::zero)
	{
		uTransl = (uTransl & constraintDir_)*constraintDir_;
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
