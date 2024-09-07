import setuptools

setuptools.setup(
    name="optimo",
    version="0.1",
    author="Gokul Krishnan",
    author_email="gk7huki@gmail.com",
    description="Algorithms for Multi-Objective optimization.",
    long_description="",
    long_description_content_type="text/markdown",
    url="https://github.com/gk7huki/optimo",
    packages=setuptools.find_packages(include=['optimo']),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 2 - Pre-Alpha"
    ],
    python_requires='>=3.6',
    install_requires=['matplotlib>=3'],
    platforms='any'
)
