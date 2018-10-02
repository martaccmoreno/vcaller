from setuptools import setup, find_packages
import os

# Utility function to read the README file.
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='vcaller',
    author='Marta Moreno',
    version='0.1',
    author_email='mccm20@gmail.com',
    description='A wrapper CLI for variant calling pipelines',
    long_description=read('README.md'),
    #license = "BSD",
    keywords='variant calling pipeline',
    url = 'https://github.com/martaccmoreno/vcaller',
    include_package_data=True,
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research'
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Operating System :: POSIX :: Linux'
       # "License :: OSI Approved :: BSD License",
    ],
    py_modules=['vcaller'],
    install_requires=[
        'click',
        'progress'
    ],
    entry_points='''
        [console_scripts]
        vcaller=vcaller:cli
    ''',
)