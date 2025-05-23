from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'Ligand and receptor analysis package'
LONG_DESCRIPTION = 'Package with useful functions for ligand-receptor analysis'

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="mellon", 
        version=VERSION,
        author="Melissa Grant-Peters",
        author_email="<melissagrantpeters@gmail.com>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=[], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'
        
        keywords=['python', 'bioinformatics'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Research",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
        ]
)