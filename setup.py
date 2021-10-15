import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="tess-point",
    version="0.6.2",
    author="Christopher J. Burke",
    author_email="cjburke@mit.edu",
    description="Determine pixel coordinates for TESS targets",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/christopherburke/tess-point",
    license="MIT",
    install_requires=[
        'numpy',
        'astropy',
        'scipy'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy"
    ],
    py_modules=['tess_stars2px'],
)
