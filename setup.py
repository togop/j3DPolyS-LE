from setuptools import setup, find_packages

setup(name='py3dpolys_le',
      version='2021.5.12',
      description='3D Polymer Simulations - Loop Extrusion model',
      url='https://gitlab.com/togop/3DPolyS-LE',
      author='Daniel Jost, Todor Gitchev',
      author_email='daniel.jost@ens-lyon.fr, todor.gitchev@gmail.com',
      license='MIT',
      classifiers=[  # https://pypi.org/classifiers/
          'Development Status :: 4 - Beta',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Fortran',
          'Programming Language :: Python :: 3.8',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          ],
      keywords='polymer chromosome 3D dynamics simulations loop-extrusion',
      python_requires='>3.8.0',
      packages=find_packages(exclude=['test']),
      include_package_data=True,
      # package_data={'bin': ['3dpolys_le'']},
      # package_data={'data': ['input.dat']},
      # package_dir={'py3dpolys_le': 'py3dpolys_le'},
      # setup_requires=['cython', 'numpy'],
      install_requires=[
          'numpy',
          'pandas',
          'matplotlib',
          'scipy',
          'filelock',
          'h5py',
          'pyranges',
          'dask',
          'cooler'],
      entry_points={
          'console_scripts': [
              'plot_hic = py3dpolys_le.plot_hic:main',
              '3dpolys_le_runner = py3dpolys_le.3dpolys_le_runner:main',
              '3dpolys_le_stats = py3dpolys_le.3dpolys_le_stats:main',
              # 'hdf5_to_cool = py3dpolys_le.hic_converters:hdf5_to_cool', TODO was not ready yet
          ],
      },
      zip_safe=False)
