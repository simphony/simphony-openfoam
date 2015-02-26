""" foam_files module

Module for writing OpenFOAM initial files 

"""
import os
from foam_templates import head, dictionaryTemplates, scalarTemplates, vectorTemplates


class FoamFiles():

    def __init__(self):
        pass

    def createFileContent(self, solver):
        """ create content mapping to files

        """   
        version = '2.2'
        foamClass = 'dictionary'
        fileContent = {}
        for solver in dictionaryTemplates:
            for foamFile in dictionaryTemplates[solver]:
                foamClass = 'dictionary'
                location = '\"' + os.path.dirname(foamFile) + '\"'
                foamObject = os.path.basename(foamFile)
                heading = head.format(version=version, foamclass=foamClass,
                                      location=location, foamobject=foamObject)
                fileContent[foamFile] = heading +\
                    dictionaryTemplates[solver][foamFile]

        for solver in scalarTemplates:
            for foamFile in scalarTemplates[solver]:
                foamClass = 'volScalarField'
                location = '\"' + os.path.dirname(foamFile) + '\"'
                foamObject = os.path.basename(foamFile)
                heading = head.format(version=version, foamclass=foamClass,
                                      location=location, foamobject=foamObject)
                fileContent[foamFile] = heading +\
                    scalarTemplates[solver][foamFile]

        for solver in vectorTemplates:
            for foamFile in vectorTemplates[solver]:
                foamClass = 'volVectorField'
                location = '\"' + os.path.dirname(foamFile) + '\"'
                foamObject = os.path.basename(foamFile)
                heading = head.format(version=version, foamclass=foamClass,
                                      location=location, foamobject=foamObject)
                fileContent[foamFile] = heading +\
                    vectorTemplates[solver][foamFile]

        return fileContent

    def createDirectories(self, caseDirectory):
        """ create default directories

        """   
        directories = ('constant', 'system', '0',
                       os.path.join('constant', 'polyMesh'))
        for dir in directories:
            directory = os.path.join(caseDirectory, dir)
            if not os.path.exists(directory):
                os.makedirs(directory)

    def writeFoamFiles(self, caseDirectory, solver):
        """ write OpenFOAm -files base on solver attribute to given directory
        
        Parameters
        ----------
        caseDirectory : caseDirectory
            directory to write files on.
        solver : solver
            name of the solver


        Raises
        ------
        Exception when IOError occurs.

        """   

        self.createDirectories(caseDirectory)
        fileContent = self.createFileContent(solver)
        for file in fileContent:
            try:
                f = open(os.path.join(caseDirectory, file), 'w')
                f.write(fileContent[file])
                f.close()
            except IOError:
                error_str = "Can't write file: {}"
                raise ValueError(error_str.format(os.path.join(caseDirectory,
                                                               file)))
