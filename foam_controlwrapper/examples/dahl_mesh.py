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
    hex (0 1 2 3 6 7 8 9) (2 4 1) simpleGrading (1 1 1)
    hex (3 2 4 5 9 8 10 11) (2 36 1) simpleGrading (1 1 1)
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
    boundary0
    {
        type patch;
        faces
        (
            (0 6 9 3)
            (3 9 11 5)
        );
    }
    boundary1
    {
        type patch;
        faces
        (
            (1 2 8 7)
        );
    }
    boundary2
    {
        type wall;
        faces
        (
            (0 1 7 6)
        );
    }
    boundary3
    {
        type wall;
        faces
        (
            (2 4 10 8)
        );
    }
    boundary4
    {
        type patch;
        faces
        (
            (5 11 10 4)
        );
    }
    boundary5
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
