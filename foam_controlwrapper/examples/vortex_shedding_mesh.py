blockMeshDict = """
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.1.x                                 |
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
    (-23 -25 0)			//0
    ( 23 -25 0)			//1
    (40 -25 0)			//2
    (-0.7071 -0.7071 0)		//3
    ( 0.7071 -0.7071 0)		//4
    (-0.7071  0.7071 0)		//5
    ( 0.7071  0.7071 0)		//6
    (-23 25 0)			//7
    ( 23 25 0)			//8
    (40 25 0)			//9


    (-23 -25 1)
    ( 23 -25 1)
    (40 -25 1)
    (-0.7071 -0.7071 1)
    ( 0.7071 -0.7071 1)
    (-0.7071  0.7071 1)
    ( 0.7071  0.7071 1)
    (-23 25 1)
    ( 23 25 1)
    (40 25 1)
);

blocks
(
    hex (0 3 5 7 10 13 15 17) (50 20 1) simpleGrading (0.20 1 1)
    hex (1 4 3 0 11 14 13 10) (50 20 1) simpleGrading (0.20 1 1)
    hex (4 1 8 6 14 11 18 16) (50 30 1) simpleGrading (5 1 1)
    hex (6 8 7 5 16 18 17 15) (50 10 1) simpleGrading (5 1 1)
    hex (1 2 9 8 11 12 19 18) (30 30 1) simpleGrading (1 1 1)
);

edges
(
   arc 3 5 (-1 0 0)
   arc 13 15 (-1 0 1)
   arc 4 6 (1 0 0)
   arc 14 16 (1 0 1)
   arc 5 6 (0 1 0)
   arc 15 16 (0 1 1)
   arc 3 4 (0 -1 0)
   arc 13 14(0 -1 1)
);

// inlet    boundary0
// outlet   boundary1
// top      boundary2
// bottom   boundary3
// cylinder boundary4
// sides    boundary5

boundary
(
    boundary0
    {
        type wall;
        faces
        (
            (17 7 0 10)
        );
    }

    boundary1
    {
        type wall;
        faces
        (
            (9 19 12 2)
        );
    }

    boundary2
    {
        type wall;
        faces
        (
            (8 7 17 18)
            (9 8 18 19)
         );
    }

    boundary3
    {
        type wall;
        faces
        (
            (0 1 11 10)
            (1 2 12 11)
         );
    }

    boundary4
    {
        type wall;
        faces
        (
            (5 6 16 15)
            (16 6 4 14)
            (4 3 13 14)
            (3 5 15 13)
         );
    }

    boundary5
    {
        type empty;
        faces
        (
            (7 8 6 5)
            (7 5 3 0)
            (3 4 1 0)
            (8 1 4 6)
            (8 9 2 1)

            (18 17 15 16)
            (15 17 10 13)
            (14 13 10 11)
            (18 16 14 11)
            (19 18 11 12)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //

"""
