from distutils.core import setup, Extension
import os

openfoam_src_dir = os.environ['FOAM_SRC']
openfoam_libbin = os.environ['FOAM_LIBBIN']
openfoam_user_libbin = os.environ['FOAM_USER_LIBBIN']
module = Extension('simphonyfoaminterface',
                   include_dirs=[openfoam_src_dir +
                                 '/finiteVolume/lnInclude',
                                 'lnInclude',
                                 openfoam_src_dir +
                                 '/OpenFOAM/lnInclude',
                                 openfoam_src_dir +
                                 '/OSspecific/POSIX/lnInclude'],
                   libraries=['foaminterface',
                              'finiteVolume',
                              'genericPatchFields'],
                   library_dirs=[openfoam_user_libbin, openfoam_libbin],
                   extra_compile_args=['-Dlinux64 -DWM_DP -DNoRepository'],
                   sources=['libpythonfoam.cpp'])
setup(name='FoamInterface',
      version='0.0.1',
      description='Foam interface to Python',
      author='Kai Hiltunen',
      author_email='kai.hiltunen@numerola.fi',
      url='http://www.numerola.fi',
      ext_modules=[module])
