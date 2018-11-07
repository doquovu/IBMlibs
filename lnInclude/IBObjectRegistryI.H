/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBObjectRegistry class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * Inline Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::IBObjectRegistry::h()
{
	return h_;
}

inline Foam::scalar Foam::IBObjectRegistry::dV()
{
	return dV_;
}

inline Foam::scalar Foam::IBObjectRegistry::rhoF()
{
	return rhoF_;
}

inline Foam::vector Foam::IBObjectRegistry::g()
{
	return g_;
}

inline Foam::label Foam::IBObjectRegistry::nObjects()
{
	return nObjects_;
}

inline Foam::PtrList<IBObject>& Foam::IBObjectRegistry::objects()
{
	return objects_;
}

inline Foam::List<pointField>& Foam::IBObjectRegistry::LPoints()
{
	return LPoints_;
}

inline Foam::List<pointField>& Foam::IBObjectRegistry::Shd1LPoints()
{
	return Shd1LPoints_;
}
inline Foam::List<pointField>& Foam::IBObjectRegistry::Shd2LPoints()
{
	return Shd2LPoints_;
}

inline Foam::List<point>& Foam::IBObjectRegistry::CoG()
{
	return CoG_;
}

inline Foam::List<point>& Foam::IBObjectRegistry::Shd1CoG()
{
	return Shd1CoG_;
}
inline Foam::List<point>& Foam::IBObjectRegistry::Shd2CoG()
{
	return Shd2CoG_;
}

inline Foam::List<labelListList>& Foam::IBObjectRegistry::neiCells()
{
	return neiCells_;
}

inline Foam::List<labelList>& Foam::IBObjectRegistry::solidCells()
{
	return solidCells_;
}

inline Foam::List<labelList>& Foam::IBObjectRegistry::solidCellsInt()
{
	return solidCellsInt_;
}

inline Foam::List<labelListList>& Foam::IBObjectRegistry::Shd1NeiCells()
{
	return Shd1NeiCells_;
}
inline Foam::List<labelListList>& Foam::IBObjectRegistry::Shd2NeiCells()
{
	return Shd2NeiCells_;
}

inline Foam::List<labelList>& Foam::IBObjectRegistry::Shd1SolidCells()
{
	return Shd1SolidCells_;
}
inline Foam::List<labelList>& Foam::IBObjectRegistry::Shd2SolidCells()
{
	return Shd2SolidCells_;
}

inline Foam::List<vector>& Foam::IBObjectRegistry::UTranslate()
{
	return UTranslate_;
}

inline Foam::List<vector>& Foam::IBObjectRegistry::URotate()
{
	return URotate_;
}

inline Foam::List<bool> Foam::IBObjectRegistry::enableShadows()
{
	return enableShadows_;
}

inline Foam::vector Foam::IBObjectRegistry::cyclicDistance()
{
	return cyclicDistance_;
}

inline Foam::List<label>& Foam::IBObjectRegistry::nFaces()
{
	return nFaces_;
}

inline Foam::List<labelList>& Foam::IBObjectRegistry::nPointsOfFaces()
{
	return nPointsOfFaces_;
}

inline Foam::List<labelListList>& Foam::IBObjectRegistry::pointOfFace()
{
	return pointOfFace_;
}

inline Foam::scalar Foam::IBObjectRegistry::HeavisideFunc(scalar phi)
{
	if (phi<=0) return 0;
	if (phi>0) return 1;
	return 0;
}

inline Foam::scalar Foam::IBObjectRegistry::LevelSetFunc
(
	label objID, 
	point p, 
	point c
)
{
	if(mesh_.nGeometricD() == 2)
	{
	    return 
	    (
	    	sqrt(pow(p.x()-c.x(),2) + pow(p.y()-c.y(),2)) 
	      - R_[objID]
	    );
	}
	else 
	{
		return ( mag(p-CoG_[objID]) - R_[objID]);
	}
}
// ************************************************************************* //