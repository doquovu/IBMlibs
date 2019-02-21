/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::IBstaticFvMesh

Description
    Foam::IBstaticFvMesh

SourceFiles
    IBstaticFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef IBstaticFvMesh_H
#define IBstaticFvMesh_H

#include "IBdynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IBstaticFvMesh Declaration
\*---------------------------------------------------------------------------*/

class IBstaticFvMesh
:
    public IBdynamicFvMesh
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        IBstaticFvMesh(const IBstaticFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const IBstaticFvMesh&);


public:

    //- Runtime type information
    TypeName("IBstaticFvMesh");


    // Constructors

        //- Construct from IOobject
        IBstaticFvMesh(const IOobject& io);


    //- Destructor
    ~IBstaticFvMesh();


    // Member Functions

        //- Is mesh dynamic
        virtual bool dynamic() const
        {
            return false;
        }

        //- Dummy update function which does not change the mesh
        virtual bool update(const labelList& neighbourCells);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
