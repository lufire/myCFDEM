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
\*---------------------------------------------------------------------------*/

#include "myCfdemCloud.H"
#include "forceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::myCfdemCloud::myCfdemCloud
(
    const fvMesh& mesh
)
:
    cfdemCloud(mesh),
    convCds_(NULL),
    convDEMForces_(NULL),
    convFluidVel_(NULL),
    unitConversion_(couplingProperties_.lookup("unitConversion"))
{}

// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //
Foam::myCfdemCloud::~myCfdemCloud()
{
    dataExchangeM().destroy(convCds_,1);
    dataExchangeM().destroy(convDEMForces_,3);
    dataExchangeM().destroy(convFluidVel_,3);
}

// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //
void Foam::myCfdemCloud::getDEMdata()
{
    if (unitConversion_ == "SI-SI" || unitConversion_ == "CGS-CGS")
    {
        cfdemCloud::getDEMdata();
    }
    else if (unitConversion_ == "SI-CGS")
    {
        dataExchangeM().getData("radius","scalar-atom",radii_);
        dataExchangeM().getData("x","vector-atom",positions_);
        dataExchangeM().getData("v","vector-atom",velocities_);
        // Convert from cgs- to si-unit system
        for(int particleI=0; particleI<numberOfParticles(); ++particleI)
        {
            for(int coordI=0; coordI<1; coordI++)
            {
                radii_[particleI][coordI] = radii_[particleI][coordI]*1e-2;
            }
            for(int coordI=0; coordI<3; ++coordI)
            {
                positions_[particleI][coordI] = 
                    positions_[particleI][coordI]*1e-2;
                velocities_[particleI][coordI] = 
                    velocities_[particleI][coordI]*1e-2;
            }
        }
        if(impDEMdragAcc_)
        {
            // array is used twice - might be necessary to clean it first
            dataExchangeM().getData("dragAcc","vector-atom",fAcc_); 
            // Convert from cgs- to si-unit system
            for(int particleI=0; particleI<numberOfParticles(); ++particleI)
            {
                for(int coordI=0; coordI<3; ++coordI)
                {
                    fAcc_[particleI][coordI] = 
                        fAcc_[particleI][coordI]*1e-5;
                }
            }
        }
    }
    else
    {

        FatalError
            << "myCfdemCloud::getDEMdata(): " << endl
            << "    unknown unitConversion type "
            << unitConversion_ << "." << endl
            << "    Valid locateModel types are :" << endl
            << "SI-SI" << endl
            << "SI-CGS" << endl
            << "CGS-CGS" << endl
            << abort(FatalError);
    }
}

void Foam::myCfdemCloud::giveDEMdata()
{
    if(forceM(0).coupleForce())
    {
        if (unitConversion_ == "SI-SI" || unitConversion_ == "CGS-CGS")
        {
            cfdemCloud::giveDEMdata();
        }
        else if (unitConversion_ == "SI-CGS")
        {
            // Conversion from si- to cgs-unit system
            for(int particleI=0; particleI<numberOfParticles(); ++particleI)
            {
                for(int coordI=0; coordI<3; coordI++)
                {
                    convDEMForces_[particleI][coordI] = DEMForces_[particleI][coordI]*1e5;
                }
            }
            dataExchangeM().giveData("dragforce","vector-atom", convDEMForces_);

            if(impDEMdrag_)
            {
                if(verbose_) Info << "sending Ksl and uf" << endl;
                
                // Conversion from si- to cgs-unit system
                for(int particleI=0; particleI<numberOfParticles(); ++particleI)
                {
                    for(int coordI=0; coordI<1; coordI++)
                    {
                        convCds_[particleI][coordI] = Cds_[particleI][coordI]*1e3;
                    }
                    for(int coordI=0; coordI<3; coordI++)
                    {
                        convFluidVel_[particleI][coordI] = fluidVel_[particleI][coordI]*1e2;
                    }
                }
                if(verbose_) Info << "sending Ksl and uf" << endl;
                dataExchangeM().giveData("Ksl","scalar-atom", convCds_);
                dataExchangeM().giveData("uf","vector-atom",convFluidVel_);
            }
            if(verbose_) Info << "giveDEMdata done." << endl;
        }
        else
        {

            FatalError
                << "myCfdemCloud::giveDEMdata(): " << endl
                << "    unknown unitConversion type "
                << unitConversion_ << "." << endl
                << "    Valid locateModel types are :" << endl
                << "SI-SI" << endl
                << "SI-CGS" << endl
                << "CGS-CGS" << endl
                << abort(FatalError);
        }
    }
}

// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //
//
// * * * * * * * * * * * * * * * ACCESS  * * * * * * * * * * * * * //

vector Foam::myCfdemCloud::position(int index) const
{
    vector pos;
    for(int i=0;i<3;i++) pos[i] = positions()[index][i];
    return pos;
}


// * * * * * * * * * * * * * * * WRITE  * * * * * * * * * * * * * //

// * * *   write cfdemCloud internal data   * * * //
//
bool Foam::myCfdemCloud::reAllocArrays() const
{

    if(numberOfParticlesChanged_ && !arraysReallocated_)
    {
        // get arrays of new length
        dataExchangeM().allocateArray(convCds_,0.,1);
        dataExchangeM().allocateArray(convDEMForces_,0.,3);
        dataExchangeM().allocateArray(convFluidVel_,0.,3);
        cfdemCloud::reAllocArrays();
        return true;
    }
    return false;
}

bool Foam::myCfdemCloud::reAllocArrays(int nP, bool forceRealloc) const
{
    if( (numberOfParticlesChanged_ && !arraysReallocated_) || forceRealloc)
    {
        // get arrays of new length
        dataExchangeM().allocateArray(convCds_,0.,1,nP);
        dataExchangeM().allocateArray(convDEMForces_,0.,3,nP);
        dataExchangeM().allocateArray(convFluidVel_,0.,3,nP);
        cfdemCloud::reAllocArrays(nP,forceRealloc);
        return true;
    }
    return false;
}
