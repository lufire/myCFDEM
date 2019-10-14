/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

    cloud class managing DEM data for CFD-DEM coupling

Class
    Foam::myCfdemCloud

SourceFiles
    myCfdemCloud.C
    myCfdemCloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef myCfdemCloud_H
#define myCfdemCloud_H

#include "cfdemCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class myCfdemCloud Declaration
\*---------------------------------------------------------------------------*/

class myCfdemCloud 
:
    public cfdemCloud
{

// protected data
protected:

    mutable double **convCds_;

    mutable double **convDEMForces_;
    
    mutable double **convFluidVel_;

    const word unitConversion_;

// Protected member functions
    virtual void getDEMdata();

    virtual void giveDEMdata();

public:

    //- Construct from mesh and a list of particles
    myCfdemCloud
    (
            const fvMesh& mesh
    );

    //- Destructor
    virtual ~myCfdemCloud();

// public Member Functions

    // Access
    
       vector position(int) const;

       inline scalar radius(int) const;

    // Write

      // write cfdemCloud internal data
        virtual bool reAllocArrays() const;

        virtual bool reAllocArrays(int nP, bool forceRealloc) const;  //force number of particles during reallocation
        inline const word& unitConversion() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "myCfdemCloudI.H"

#endif

// ************************************************************************* //
