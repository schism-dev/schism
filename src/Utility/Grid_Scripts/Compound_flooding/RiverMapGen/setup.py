import setuptools
import io

with io.open('README.md','r', encoding='utf8') as fh:
  long_description = fh.read()

  setuptools.setup(
  name='RiverMapGen',
  version='0.0.1',
  author='Fei Ye',
  author_email='feiye@vims.edu',
  description='Python tools for generating watershed river arcs for meshing',
  long_description=long_description,
  long_description_content_type="text/markdown",
  url='',
  project_urls = {
    "Issues": ""
  },
  license='MIT',
  packages=[
    'RiverMapGen',
  ],
  package_data={'RiverMapGen': ['Datafiles/*']},
  install_requires=[
    'numpy',
    'pandas',
    'gdal',
    'shapely',
    'pyshp',
    'geopandas'
  ],
)
