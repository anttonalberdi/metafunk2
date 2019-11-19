import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='metafunk2',  
     version='2.0',
     scripts=['metafunk2'] ,
     author="Antton Alberdi",
     author_email="antton.alberdi@gbio.ku.dk",
     description="Taxonomic and functional metagenomics pipeline",
     long_description=long_description,
   long_description_content_type="text/markdown",
     url="https://github.com/javatechy/dokr",
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: GPL-3",
         "Operating System :: OS Independent",
     ],
 )
