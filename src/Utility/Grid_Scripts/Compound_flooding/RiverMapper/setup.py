import setuptools
import io

with io.open('README.md','r', encoding='utf8') as fh:
  long_description = fh.read()

  setuptools.setup(
  name='RiverMapper',
  version='0.0.4',
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
    'RiverMapper',
  ],
  package_data={'RiverMapper': ['Datafiles/*']},
  install_requires=[
    'numpy',
    'pandas',
    'gdal>=3.6.0',
    'shapely>=2.0.0',
    'pyshp>=2.3.0',
    'geopandas>=0.12.0',
  ],
)
