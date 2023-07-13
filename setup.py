from setuptools import setup

setup(
    name="hubmap-dbgap",
    version="1.0",
    description="Generates submission for dbGaP",
    url="https://github.com/hubmapconsortium/py-hubmap-dbgap",
    author="Ivan Cao-Berg, Gesina Phillips",
    author_email="icaoberg@psc.edu",
    install_requires=[
        "pandas",
        "numpy",
        "tabulate",
        "tqdm",
    ],
    packages=["hubmapdbgap"],
)
