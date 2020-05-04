from setuptools import setup

setup(name='pyliger',
      version='0.1',
      description='The Python version of LIGER package.',
      url='https://github.com/welch-lab/pyliger',
      author='',
      author_email='',
      license='MIT',
      packages=['pyliger'],
      keywords='LIGER',
      install_requires=[
              'pandas',
              'numpy',
              'anndata'
              ],
      zip_safe=False)