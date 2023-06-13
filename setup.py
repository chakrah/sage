from setuptools import setup, find_packages

setup(name="SAGE",
      version="1.0.0",
      description="Stellar Activity Grid for Exoplanets",
      author="Hritam Chakraborty",
      author_email="Hritam.Chakraborty@unige.ch",
      packages=setuptools.find_packages(),
      install_requires= ['numpy', 'matplotlib', 'scipy', 'astropy']
      )