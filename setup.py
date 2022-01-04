from setuptools import setup, find_packages

with open('README.rst', 'r') as f_readme:
    long_description = f_readme.read()

setup(
    name='adcpyproc',
    version='0.0.1',
    description=('Post-processing and basic analysis of ocean ADCP'
                 ' data using python.'),
    py_modules=['rdi_adcp', 'rdi_toolbox', 'variables_dictionary',
                'internal_functions'],
    package_dir={'':'src'},
    packages=['adcpyproc'],
    package_data={'adcpyproc': ['WMM_COF_files/*.COF']},
    setup_requires=['numpy'],
    install_requires=[
        'numpy',
        'scipy',
        'datetime',
        ],
    extras_require={
        'save_to_pickle':['pickle'],
        'magdec_correction':['geomag']
    },
    license='MIT', 
    author='Øyvind Lundesgaard',
    author_email='oyvind.lundesgaard@npolar.no',
    classifiers=[
        'Development status :: 3 - Alpha',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Operating System :: OS Independent',
    ],
    long_description=long_description,
    long_description_content_type='text/x-rst',
)