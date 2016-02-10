blockMeshDict = """
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.2.1                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

vertices
(
    (0 0 -0.1)
    (8.65 0 -0.1)
    (8.65 0.1 -0.1)
    (0 0.1 -0.1)
    (8.65 1 -0.1)
    (0 1 -0.1)
    (0 0 0.1)
    (8.65 0 0.1)
    (8.65 0.1 0.1)
    (0 0.1 0.1)
    (8.65 1 0.1)
    (0 1 0.1)
);

blocks
(
    hex (0 1 2 3 6 7 8 9) (200 4 1) simpleGrading (1 1 1)
    hex (3 2 4 5 9 8 10 11) (200 36 1) simpleGrading (1 1 1)
);

edges
(
);

// boundary0 -> inlet
// boundary1 -> outlet
// boundary2 -> bottomWall
// boundary3 -> endWall
// boundary4 -> top
// boundary5 -> frontAndBack

boundary
(
    inlet
    {
        type patch;
        faces
        (
            (0 6 9 3)
            (3 9 11 5)
        );
    }
    outlet
    {
        type patch;
        faces
        (
            (1 2 8 7)
        );
    }
    bottomWall
    {
        type wall;
        faces
        (
            (0 1 7 6)
        );
    }
    endWall
    {
        type wall;
        faces
        (
            (2 4 10 8)
        );
    }
    top
    {
        type patch;
        faces
        (
            (5 11 10 4)
        );
    }
    frontAndBack
    {
        type empty;
        faces
        (
            (0 3 2 1)
            (6 7 8 9)
            (3 5 4 2)
            (9 8 10 11)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //


"""
