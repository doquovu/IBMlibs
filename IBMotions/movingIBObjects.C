/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - movingIBObjects class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "movingIBObjects.H"
#include "sixDoFMotion.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(movingIBObjects, 0);
}

//---------------------------------Constructors------------------------------//
Foam::movingIBObjects::movingIBObjects
(
	const fvMesh& mesh,
	const dictionary& dict
)
:
	IBObjectRegistry(mesh,dict),
	mesh_(mesh),
	uBoundary_(),
	calcRepulsive_(dict.lookupOrDefault<Switch>("calcRepulsive", false)),
	epsilonP_(dict.lookupOrDefault<scalar>("epsilonP",1e-07)),
	epsilonW_(dict.lookupOrDefault<scalar>("epsilonW",1e-07))
{
	calcUBoundary();
}

// -------------------------------Member Functions----------------------------//

void Foam::movingIBObjects::calcUBoundary()
{
	uBoundary_.setSize(nObjects());
	for(int i=0; i<nObjects(); i++)
	{	
		uBoundary_[i].setSize(LPoints()[i].size());
	    forAll(uBoundary_[i], pI)
	    {
	        uBoundary_[i][pI] = 
	        	UTranslate()[i] 
	            + (URotate()[i] ^ (LPoints()[i][pI] - CoG()[i]));
	    }
	}
}

void  Foam::movingIBObjects::movePoints()
{
    for(int i=0; i<nObjects(); i++)
    {
        const scalar dT = mesh_.time().deltaTValue();
        if (objects()[i].movable())
        {
            Info<< "IBM: Moving immersed object "<<i+1
                <<" named "<<objects()[i].name()<<endl;

            if (bendingAngle() == 0)
            {
                //- Move Lagrang point according to its velocity
                CoG()[i] += dT*UTranslate()[i];
                forAll(LPoints()[i], pI)
                {
                    LPoints()[i][pI] += dT*uBoundary_[i][pI];
                }
                Info<<"    New CoG     = "<<CoG()[i]<<endl;
                Info<<"    Ut Object "<<i+1<<" = "<<UTranslate()[i]<<endl;
                Info<<"    Ur Object "<<i+1<<" = "<<URotate()[i]<<endl;

                //- HANLE MOTION THROUGH PERIODIC BOUNDARY
                if (enableShadows()[i])
                {
                    Info<< "    Moving shadow of "<<objects()[i].name()<<nl<<endl;
                    Shd1CoG()[i] += dT*UTranslate()[i];
                    Shd2CoG()[i] += dT*UTranslate()[i];
                    forAll(Shd1LPoints()[i], pI)
                    {
                        Shd1LPoints()[i][pI] += dT*uBoundary_[i][pI];
                        Shd2LPoints()[i][pI] += dT*uBoundary_[i][pI];
                    }

                    if (isInsideRegion1(i, CoG()[i]))
                    {
                        Shd1NeiCells()[i].clear();
                        Shd2NeiCells()[i].clear();
                        Shd1SolidCells()[i].clear();
                        Shd2SolidCells()[i].clear();
                        neiCells()[i] = findNeiCells(LPoints()[i]);
                        solidCells()[i] = findSolidCells(i, CoG()[i]);
                        continue;
                    }
                    else if (isInsideRegion2(i, CoG()[i]))
                    {
                        if (isInsideRegion3(i, Shd1CoG()[i]))
                        {
                            Shd1NeiCells()[i] = findNeiCells(Shd1LPoints()[i]);
                            Shd2NeiCells()[i].clear();
                            Shd1SolidCells()[i] = findSolidCells(i, Shd1CoG()[i]);
                            Shd2SolidCells()[i].clear();
                            neiCells()[i] = findNeiCells(LPoints()[i]);
                            solidCells()[i] = findSolidCells(i, CoG()[i]);
                            continue;
                        }
                        else
                        {
                            Shd1NeiCells()[i].clear();
                            Shd2NeiCells()[i] = findNeiCells(Shd2LPoints()[i]);
                            Shd1SolidCells()[i].clear();
                            Shd2SolidCells()[i] = findSolidCells(i, Shd2CoG()[i]);
                            neiCells()[i] = findNeiCells(LPoints()[i]);
                            solidCells()[i] = findSolidCells(i, CoG()[i]);
                            continue;
                        }
                    }
                    else
                    {
                        if (isInsideRegion2(i, Shd1CoG()[i]))
                        {
                            swap(CoG()[i], Shd1CoG()[i]);
                            swap(Shd1CoG()[i], Shd2CoG()[i]);
                            Shd1CoG()[i] = CoG()[i] - cyclicDistance();

                            forAll(LPoints()[i], pI)
                            {
                                swap(LPoints()[i][pI], Shd1LPoints()[i][pI]);
                                swap(Shd1LPoints()[i][pI], Shd2LPoints()[i][pI]);
                                Shd1LPoints()[i][pI] = LPoints()[i][pI] - cyclicDistance();
                            }
                                
                            Shd1NeiCells()[i].clear();
                            Shd2NeiCells()[i] = findNeiCells(Shd2LPoints()[i]);
                            Shd1SolidCells()[i].clear();
                            Shd2SolidCells()[i] = findSolidCells(i, Shd2CoG()[i]);
                            neiCells()[i] = findNeiCells(LPoints()[i]);
                            solidCells()[i] = findSolidCells(i, CoG()[i]);
                            continue;
                        }
                        else if (isInsideRegion2(i, Shd2CoG()[i]))
                        {
                            swap(CoG()[i], Shd2CoG()[i]);
                            swap(Shd2CoG()[i], Shd1CoG()[i]);
                            Shd2CoG()[i] = CoG()[i] + cyclicDistance();

                            forAll(LPoints()[i], pI)
                            {
                                swap(LPoints()[i][pI], Shd2LPoints()[i][pI]);
                                swap(Shd2LPoints()[i][pI], Shd1LPoints()[i][pI]);
                                Shd2LPoints()[i][pI] = LPoints()[i][pI] + cyclicDistance();
                            }

                            Shd1NeiCells()[i] = findNeiCells(Shd1LPoints()[i]);
                            Shd2NeiCells()[i].clear();
                            Shd1SolidCells()[i] = findSolidCells(i, Shd1CoG()[i]);
                            Shd2SolidCells()[i].clear();
                            neiCells()[i] = findNeiCells(LPoints()[i]);
                            solidCells()[i] = findSolidCells(i, CoG()[i]);
                            continue;
                        }
                    }
                }
            }
            else
            {
                //- Move Lagrang point according to its velocity
                CoG()[i] += dT*UTranslate()[i];
                forAll(LPoints()[i], pI)
                {
                    LPoints()[i][pI] += dT*uBoundary_[i][pI];
                }
                Info<<"    New CoG     = "<<CoG()[i]<<endl;
                Info<<"    Ut Object "<<i+1<<" = "<<UTranslate()[i]<<endl;
                Info<<"    Ur Object "<<i+1<<" = "<<URotate()[i]<<endl;

                //- HANLE MOTION THROUGH PERIODIC BOUNDARY
                if (enableShadows()[i])
                {
                    Info<< "    Moving shadow of "<<objects()[i].name()<<nl<<endl;
                   
                    // Shd1CoG()[i] += dT*rotate(UTranslate()[i], CoG()[i], -bendingAngle());
                    // Shd2CoG()[i] += dT*rotate(UTranslate()[i], CoG()[i], bendingAngle());
                    // forAll(Shd1LPoints()[i], pI)
                    // {
                    //     Shd1LPoints()[i][pI] += dT*rotate(uBoundary_[i][pI], LPoints()[i][pI], -bendingAngle());
                    //     Shd2LPoints()[i][pI] += dT*rotate(uBoundary_[i][pI], LPoints()[i][pI], bendingAngle());
                    // }
                    Shd1CoG()[i] = rotate(CoG()[i], -bendingAngle());
                    Shd2CoG()[i] = rotate(CoG()[i],  bendingAngle());

                    forAll(Shd1LPoints()[i], pI)
                    {
                        Shd1LPoints()[i][pI] = rotate(LPoints()[i][pI], -bendingAngle());
                        Shd2LPoints()[i][pI] = rotate(LPoints()[i][pI],  bendingAngle());
                    }

                    if (isInsideRegion1(i, CoG()[i]))
                    {
                        Shd1NeiCells()[i].clear();
                        Shd2NeiCells()[i].clear();
                        Shd1SolidCells()[i].clear();
                        Shd2SolidCells()[i].clear();
                        neiCells()[i] = findNeiCells(LPoints()[i]);
                        solidCells()[i] = findSolidCells(i, CoG()[i]);
                        continue;
                    }
                    else if (isInsideRegion2(i, CoG()[i]))
                    {
                        if (isInsideRegion3(i, Shd1CoG()[i]))
                        {
                            Shd1NeiCells()[i] = findNeiCells(Shd1LPoints()[i]);
                            Shd2NeiCells()[i].clear();
                            Shd1SolidCells()[i] = findSolidCells(i, Shd1CoG()[i]);
                            Shd2SolidCells()[i].clear();
                            neiCells()[i] = findNeiCells(LPoints()[i]);
                            solidCells()[i] = findSolidCells(i, CoG()[i]);
                            continue;
                        }
                        else
                        {
                            Shd1NeiCells()[i].clear();
                            Shd2NeiCells()[i] = findNeiCells(Shd2LPoints()[i]);
                            Shd1SolidCells()[i].clear();
                            Shd2SolidCells()[i] = findSolidCells(i, Shd2CoG()[i]);
                            neiCells()[i] = findNeiCells(LPoints()[i]);
                            solidCells()[i] = findSolidCells(i, CoG()[i]);
                            continue;
                        }
                    }
                    else
                    {
                        if (isInsideRegion2(i, Shd1CoG()[i]))
                        {
                            swap(CoG()[i], Shd1CoG()[i]);
                            swap(Shd1CoG()[i], Shd2CoG()[i]);
                            Shd1CoG()[i] = rotate(CoG()[i], -bendingAngle());

                            for (int pI=0; pI<LPoints()[i].size(); pI++)
                            {
                                swap(LPoints()[i][pI], Shd1LPoints()[i][pI]);
                                swap(Shd1LPoints()[i][pI], Shd2LPoints()[i][pI]);
                                Shd1LPoints()[i][pI] = rotate(LPoints()[i][pI], -bendingAngle());
                            }
                                
                            Shd1NeiCells()[i].clear();
                            Shd2NeiCells()[i] = findNeiCells(Shd2LPoints()[i]);
                            Shd1SolidCells()[i].clear();
                            Shd2SolidCells()[i] = findSolidCells(i, Shd2CoG()[i]);
                            neiCells()[i] = findNeiCells(LPoints()[i]);
                            solidCells()[i] = findSolidCells(i, CoG()[i]);
                            continue;
                        }
                        else if (isInsideRegion2(i, Shd2CoG()[i]))
                        {
                            swap(CoG()[i], Shd2CoG()[i]);
                            swap(Shd2CoG()[i], Shd1CoG()[i]);
                            Shd2CoG()[i] = rotate(CoG()[i], bendingAngle());

                            for (int pI=0; pI<LPoints()[i].size(); pI++)
                            {
                                swap(LPoints()[i][pI], Shd2LPoints()[i][pI]);
                                swap(Shd2LPoints()[i][pI], Shd1LPoints()[i][pI]);
                                Shd2LPoints()[i][pI] = rotate(LPoints()[i][pI], bendingAngle());
                            }

                            Shd1NeiCells()[i] = findNeiCells(Shd1LPoints()[i]);
                            Shd2NeiCells()[i].clear();
                            Shd1SolidCells()[i] = findSolidCells(i, Shd1CoG()[i]);
                            Shd2SolidCells()[i].clear();
                            neiCells()[i] = findNeiCells(LPoints()[i]);
                            solidCells()[i] = findSolidCells(i, CoG()[i]);
                            continue;
                        }
                    }
                }
            }

            neiCells()[i] = findNeiCells(LPoints()[i]);
            solidCells()[i] = findSolidCells(i, CoG()[i]);
        }
    }
}
void Foam::movingIBObjects::updateObjectMotions
(
	label objID,
	const vectorField& ForceLagrange
)
{
	IBObject& objectI = objects()[objID];
		
	if (!objectI.movable())
	{
		return;
	}
	else
	{
		forAll(objectI.motions(), motionI)
		{
			if (   objectI.motions()[motionI].motionType() == "linearMotion"
				|| objectI.motions()[motionI].motionType() == "rotatingMotion" )
			{
				continue;
			}
			else if (objectI.motions()[motionI].motionType() == "sixDoFMotion")
			{
				Info<<"IBM: Recalculate 6DoF Motion for Object "<<objID+1<<endl;
				sixDoFMotion& sDoF = refCast<sixDoFMotion>
				(
				    objectI.motions()[motionI]
				);
                vector repulsiveForce = objWallRepulsive(objID);
                for(int objI=0; objI<nObjects(); objI++)
                {
                    if (objID == objI)
                        continue;
                    else
                    {
                        repulsiveForce += objMutualRepulsive(objID, objI);
                    }
                }

                if (sDoF.solver() == "uhlmann")
                {
                    sDoF.updateMotion
                    (
                        UTranslate()[objID],
                        URotate()[objID],
                        mesh_,
                        ForceLagrange,
                        repulsiveForce,
                        CoG()[objID],
                        LPoints()[objID],
                        rhoF(),
                        dV(),
                        g()
                    );
                }
                else if (sDoF.solver() == "tobias")
                {
                    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
                    const volVectorField& Uold = U.oldTime();
                    const scalar dT = mesh_.time().deltaTValue();

                    vector volIntegralU_ 
                    = 
                    (
                        volIntegralU(objID, U) 
                      - volIntegralU(objID, Uold)
                    ) / dT;

                    vector volIntegralRxU_ 
                    = 
                    (
                        volIntegralRxU(objID, U) 
                      - volIntegralRxU(objID, Uold)
                    ) / dT;
                    sDoF.updateMotion
                    (
                        UTranslate()[objID],
                        URotate()[objID],
                        mesh_,
                        ForceLagrange,
                        repulsiveForce,
                        CoG()[objID],
                        LPoints()[objID],
                        rhoF(),
                        dV(),
                        g(),
                        volIntegralU_,
                        volIntegralRxU_
                    );   
                }
                else
                {
                    FatalErrorIn
                    (
                        "movingIBObjects::updateObjectMotions(label, "
                        "const vectorField&)" 
                    )	<< "Unknown solver for sixDoF motion " << endl
                        << "Valid solver are: 2( tobias, uhlmann )"<<endl
                        <<exit(FatalError);
                }
                forAll(uBoundary_[objID], pI)
				{
				    uBoundary_[objID][pI] 
                    = 
				    	UTranslate()[objID] 
				      + (URotate()[objID]
                      ^ (LPoints()[objID][pI]-CoG()[objID]) );
				}
			}
		}
	}
}

Foam::vector Foam::movingIBObjects::volIntegralU
(
    label objID,
    const volVectorField& uField
)
{   
    vector IU(vector::zero);
    volScalarField volumeFraction
    (
        IOobject
        (
            "volumeFraction",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("volumeFraction", dimless, 0)
    ); 
    forAll(solidCells()[objID], cellI)
    {
        scalar vf = volFraction(objID, solidCells()[objID][cellI], CoG()[objID]);
        volumeFraction[solidCells()[objID][cellI]] = vf;
        IU += uField[solidCells()[objID][cellI]]*dV()*vf;
    }

    if (enableShadows()[objID])
    {
        if (!Shd1SolidCells()[objID].empty())
        {
            forAll(Shd1SolidCells()[objID], cellI)
            {
                scalar vf = 
                    volFraction
                    (
                        objID,
                        Shd1SolidCells()[objID][cellI], 
                        Shd1CoG()[objID]
                    );
                volumeFraction[Shd1SolidCells()[objID][cellI]] = vf;
                IU += uField[Shd1SolidCells()[objID][cellI]]*dV()*vf;
            }
        }
        if (!Shd2SolidCells()[objID].empty())
        {
            forAll(Shd2SolidCells()[objID], cellI)
            {
                scalar vf = 
                    volFraction
                    (
                        objID,
                        Shd2SolidCells()[objID][cellI], 
                        Shd2CoG()[objID]
                    );
                volumeFraction[Shd2SolidCells()[objID][cellI]] = vf;
                IU += uField[Shd2SolidCells()[objID][cellI]]*dV()*vf;
            }
        }
    }
    if (mesh_.time().outputTime() || mesh_.time().timeIndex() == 0)
    {
        volumeFraction.write();
    }

    return IU;
}

Foam::vector Foam::movingIBObjects::volIntegralRxU
(
    label objID,
    const volVectorField& U
)
{
    const volVectorField& cc = mesh_.C();

    vector IRxU(vector::zero);
    forAll(solidCells()[objID], cellI)
    {
        scalar vf = volFraction(objID, solidCells()[objID][cellI], CoG()[objID]);
        IRxU +=
        ( 
            (cc[solidCells()[objID][cellI]] - CoG()[objID]) ^ U[solidCells()[objID][cellI]]
          * dV()
          * vf
        );
    } 

    if (enableShadows()[objID])
    {
        if (!Shd1SolidCells()[objID].empty())
        {
            forAll(Shd1SolidCells()[objID], cellI)
            {
                scalar vf = 
                    volFraction
                    (
                        objID,
                        Shd1SolidCells()[objID][cellI], 
                        Shd1CoG()[objID]
                    );
                IRxU +=
                ( 
                    (
                        (cc[Shd1SolidCells()[objID][cellI]] - Shd1CoG()[objID])
                      ^ U[Shd1SolidCells()[objID][cellI]]
                    )
                  * dV()
                  * vf
                );
            }
        }
        if (!Shd2SolidCells()[objID].empty())
        {
            forAll(Shd2SolidCells()[objID], cellI)
            {
                scalar vf = 
                    volFraction
                    (
                        objID,
                        Shd2SolidCells()[objID][cellI], 
                        Shd2CoG()[objID]
                    );
                IRxU +=
                ( 
                    (
                        (cc[Shd2SolidCells()[objID][cellI]] - Shd2CoG()[objID]) 
                      ^ U[Shd2SolidCells()[objID][cellI]]
                    )
                  * dV()
                  * vf
                );
            }
        }
    }

    return IRxU;
}

Foam::vector Foam::movingIBObjects::objMutualRepulsive
(
    label obj1ID,   //- the object is considered
    label obj2ID    
)
{
    //- Calculate particle-particle interaction
    //  Wan, Turek 2007
    vector force;
	IBObject& object1 = objects()[obj1ID];
	IBObject& object2 = objects()[obj2ID];
	scalar R1 = object1.R();
	scalar R2 = object2.R();
	point& C1 = CoG()[obj1ID];
	point& C2 = CoG()[obj2ID];

    scalar d12 = mag(C1 - C2);
    scalar xi = 2.0*h();
    scalar epsilonPP = epsilonP_;
    if (d12 <= R1+R2)
    {
        force = ((C1-C2)*(R1+R2-d12)/epsilonP_);
    }
    if ( (d12 > R1+R2) && (d12 <= R1+R2+xi) )
    {
        force = ((C1-C2)*Foam::pow(R1+R2+xi-d12, 2.0)/epsilonPP);
    }
    if ( d12 > R1+R2+xi ) 
    {
        force = vector::zero;
    }
    Info<<"IBM: Calculating repulsive force of object "<<obj1ID<<" and "<<obj2ID
                            <<" : "<<force<<endl;
    return force;
}

Foam::vector Foam::movingIBObjects::objWallRepulsive(label objID)
{
    //- Calculate particle-wall interaction
    //  Wan, Turek 2007
	IBObject& object = objects()[objID];
	scalar R = object.R();
	point C = CoG()[objID];

	vector Fwr(vector::zero);

    if(calcRepulsive_)
    {
        Info<<"IBM: Calculating wall repulsive force for "<<object.name()<<endl;
        forAll(mesh_.boundary(), patchI)
        {
            const polyPatch& pp = mesh_.boundaryMesh()[patchI];
            if(mesh_.boundary()[patchI].Cf().size() == 0)
            {
                continue;
            }
            else if(isA<wallPolyPatch>(pp))
            {
                const vectorField fc = pp.faceCentres();
                    
                //- Find face closest to C
                scalar d(GREAT);
                label nearestFace(0);
                forAll(fc, facei)
                {
                    scalar CToFace = mag(C - fc[facei]);
                    if (CToFace <= d)
                    {
                        d = CToFace;
                        nearestFace = facei;
                    }
                }

                //- Find projection of C on nearest face
                //  **projection = p - ((p-c)*n)*n
                const vector fcn = fc[nearestFace];

                //- Error: cannot access mesh_.boundary().nf()[nearestFace]???
                vector Sf_ = mesh_.boundary()[patchI].Sf()[nearestFace];
                vector nf_ = Sf_/mag(Sf_);
                vector p_fcn = C - fcn;
                scalar distance = p_fcn & nf_;
                point Cj = C - distance*nf_;

                scalar di = 2.0*mag(C - Cj);
                scalar xi = 2.0*h();
                scalar epsilonWW = epsilonW_;
                //- Theoretical value of stiffness values:

                if (di <= 2.0*R)
                {
                    Fwr += 2.0*(C-Cj)*(2.0*R-di)/epsilonW_;
                }
                else if ((di > 2.0*R) && (di <= 2.0*R+xi))
                {
                    Fwr += 2.0*(C-Cj)*Foam::pow(2.0*R+xi-di, 2.0)/epsilonWW;
                }
                else 
                    Fwr += vector::zero;
            }
            else
                continue;
        }
    }

    return Fwr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
