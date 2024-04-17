/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "codedFvModelTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "fvMatrix.H"

//{{{ begin codeInclude
#line 29 "/home/mlouis9/PythonProjects/senior-design/cfd/3Dpincell/shortFuel/constant/fvOptions/fissionHeatSource"
#include "fvCFD.H"
//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = ba3f8712cb66180f7252d01150943bead3a78fa2
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void fissionHeatSource_ba3f8712cb66180f7252d01150943bead3a78fa2(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(fissionHeatSourceFvModelscalarSource, 0);

addRemovableToRunTimeSelectionTable
(
    fvModel,
    fissionHeatSourceFvModelscalarSource,
    dictionary
);


const char* const fissionHeatSourceFvModelscalarSource::SHA1sum =
    "ba3f8712cb66180f7252d01150943bead3a78fa2";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fissionHeatSourceFvModelscalarSource::
fissionHeatSourceFvModelscalarSource
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    set_(mesh, coeffs())
{
    if (false)
    {
        Info<<"construct fissionHeatSource sha1: ba3f8712cb66180f7252d01150943bead3a78fa2"
            " from components\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fissionHeatSourceFvModelscalarSource::
~fissionHeatSourceFvModelscalarSource()
{
    if (false)
    {
        Info<<"destroy fissionHeatSource sha1: ba3f8712cb66180f7252d01150943bead3a78fa2\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fissionHeatSourceFvModelscalarSource::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (false)
    {
        Info<<"fissionHeatSourceFvModelscalarSource::addSup()\n";
    }

//{{{ begin code
    #line 38 "/home/mlouis9/PythonProjects/senior-design/cfd/3Dpincell/shortFuel/constant/fvOptions/fissionHeatSource"
const dimensionedScalar heatSource("heatSource", dimPower/dimVolume, 100000000);
        eqn += heatSource*mesh.V();
//}}} end code
}


void fissionHeatSourceFvModelscalarSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (false)
    {
        Info<<"fissionHeatSourceFvModelscalarSource::addSup()\n";
    }

//{{{ begin code
    
//}}} end code
}


void fissionHeatSourceFvModelscalarSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (false)
    {
        Info<<"fissionHeatSourceFvModelscalarSource::addSup()\n";
    }

//{{{ begin code
    
//}}} end code
}


bool fissionHeatSourceFvModelscalarSource::movePoints()
{
    set_.movePoints();
    return true;
}


void fissionHeatSourceFvModelscalarSource::topoChange(const polyTopoChangeMap& map)
{
    set_.topoChange(map);
}


void fissionHeatSourceFvModelscalarSource::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void fissionHeatSourceFvModelscalarSource::distribute
(
    const polyDistributionMap& map
)
{
    set_.distribute(map);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace fv

// ************************************************************************* //

