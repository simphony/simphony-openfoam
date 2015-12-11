/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | Unsupported Contributions for OpenFOAM
 \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2015 SimPhoNy -project
   \\/     M anipulation  |
-------------------------------------------------------------------------------


License
    This file is a derivative work of OpenFOAM.

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

SourceFiles
    TimeMod.C

\*---------------------------------------------------------------------------*/

#ifndef TimeMod_H
#define TimeMod_H

#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class TimeMod Declaration
\*---------------------------------------------------------------------------*/

class TimeMod
:
    public Time
{
    
public:

    //TypeName("time");

    //- The default control dictionary name (normally "controlDict")
    static word controlDictName;


    // Constructors

        //- Construct given name of dictionary to read and argument list
        TimeMod
        (
            const word& name,
            const argList& args,
            const word& systemName = "system",
            const word& constantName = "constant"
        ):
        Time(name,args,systemName,constantName){};

        //- Construct given name of dictionary to read, rootPath and casePath
        TimeMod
        (
            const word& name,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant",
            const bool enableFunctionObjects = true
        ):
        Time(name,rootPath,caseName,systemName,constantName,enableFunctionObjects){};

        //- Construct given dictionary, rootPath and casePath
        TimeMod
        (
            const dictionary& dict,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant",
            const bool enableFunctionObjects = true
        ):
        Time(dict,rootPath,caseName,systemName,constantName,enableFunctionObjects){};

        //- Construct given endTime, rootPath and casePath
        TimeMod
        (
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant",
            const bool enableFunctionObjects = true
        ):
        Time(rootPath,caseName,systemName,constantName,enableFunctionObjects){};


    //- Destructor
    virtual ~TimeMod();

    bool update();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
