from setuptools import setup, find_packages

setup(
    name="hubmap-dbgap",
    version="2025.01",
    description="Generates submission for dbGaP",
    url="https://github.com/hubmapconsortium/py-hubmap-dbgap",
    author="Ivan Cao-Berg, Gesina Phillips",
    author_email="icaoberg@psc.edu",
    license="MIT",  # Add a license (you can change this if needed)
    install_requires=[
        "pandas",
        "numpy",
        "tabulate",
        "tqdm",
    ],
    packages=find_packages(),  # This will automatically find your package directories
    classifiers=[  # These help users find your package based on its usage
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  # Make sure this matches your license
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",  # Adjust according to your code's Python version compatibility
    long_description_content_type="text/markdown",  # If you use a README.md file
    long_description=open('README.md').read(),  # This reads your README file (make sure you have one)
)

