/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBModel class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef manageFlowRate_H
#define manageFlowRate_H

#include "dictionary.H"
#include "IBdynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam 
{

class manageFlowRate
{
    IBdynamicFvMesh& mesh_;

    const dictionary& dict_;

    bool addPressureGradientField_;
    
    word direction_;
    
    word type_;
    
    scalar dP_;

    scalar flowRate_;

    scalar initializeDP();

public:
    TypeName("manageFlowRate");

    manageFlowRate(IBdynamicFvMesh& mesh, const dictionary& dict);

    ~manageFlowRate(){}

    void correctPressureGradient();

    volVectorField pressureField();

}; // End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //