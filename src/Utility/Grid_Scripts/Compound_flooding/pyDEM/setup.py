import setuptools
import io

with io.open('README.md','r', encoding='utf8') as fh:
  long_description = fh.read()

  setuptools.setup(
  name='pyDEM',
  version='0.0.1',
  author='Linlin Cui',
  author_email='lcui@vims.edu',
  description='Python tools for extracting thalwegs from DEMs',
  long_description=long_description,
  long_description_content_type="text/markdown",
  url='',
  project_urls = {
    "Issues": ""
  },
  license='MIT',
  packages=[
    'pyDEM',
  ],
  package_data={'pyDEM': ['Datafiles/*']},
  install_requires=[
    'numpy',
    'gdal',
    'shapely',
    'richdem',
    'fiona',
    'geopandas'
  ],
)
