import setuptools
import io

with io.open('README.md','r', encoding='utf8') as fh:
  long_description = fh.read()

  setuptools.setup(
  name='stofs3d_atlantic_preproc',
  version='0.0.1',
  description='Python tools for generating inputs for STOFS3D-ATL',
  long_description=long_description,
  long_description_content_type="text/markdown",
  url='',
  project_urls = {
    "Issues": ""
  },
  license='MIT',
  packages=[
    '',
  ],
  package_data={},
  install_requires=[
    'numpy',
    'pandas',
    'gdal>=3.6.0',
    'shapely>=2.0.0',
    'geopandas>=0.12.0',
  ],
)
