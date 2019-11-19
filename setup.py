import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='metafunk2',
     version='2.0',
     author="Antton Alberdi",
     author_email="antton.alberdi@gbio.ku.dk",
     packages=['metafunk2', 'metafunk2.test'],
     scripts=['bin/quality_filtering','bin/duplicate_removal'],
     description="Taxonomic and functional metagenomics pipeline",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/anttonalberdi/metafunk2",
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: GPL-3",
         "Operating System :: OS Independent",
     ],
 )
