from distutils.core import setup, Extension
import os

openfoam_src_dir = os.environ['FOAM_SRC']
openfoam_libbin = os.environ['FOAM_LIBBIN']
openfoam_user_libbin = os.environ['FOAM_USER_LIBBIN']
openfoam_user_src_dir = '../libs'
module = Extension('simphonyfoaminterface',
                   include_dirs=[openfoam_src_dir +
                                 '/finiteVolume/lnInclude',
                                 'lnInclude',
                                 openfoam_src_dir +
                                 '/OpenFOAM/lnInclude',
                                 openfoam_src_dir +
                                 '/meshTools/lnInclude',
                                 openfoam_src_dir +
                                 '/transportModels',
                                 openfoam_src_dir +
                                 '/transportModels/incompressible \
                                 /singlePhaseTransportModel',
                                 openfoam_src_dir +
                                 '/turbulenceModels',
                                 openfoam_user_src_dir +
                                 'incompressibleturbulenceModels/RAS/RASModel',
                                 openfoam_src_dir +
                                 '/OSspecific/POSIX/lnInclude'
                                 ],
                   libraries=['foaminterface',
                              'finiteVolume',
                              'incompressibleTurbulenceModelSimphony',
                              'incompressibleRASModelsSimphony',
                              'incompressibleTransportModels',
                              'genericPatchFields'],
                   library_dirs=[openfoam_user_libbin, openfoam_libbin],
                   extra_compile_args=['-Dlinux64 -DWM_DP -DNoRepository'],
                   sources=['libpythonfoam.cpp'])
setup(name='FoamInterface',
      version='0.0.1',
      description='Foam interface to Python',
      author='Juan Marcelo Gimenez and \
      Santiago Marquez Damian and Norberto Nigro',
      author_email='jmarcelogimenez@gmail.com',
      url='http://www.cimec.org.ar',
      ext_modules=[module])
