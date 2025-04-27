from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='PDBPrepper',
    version='0.2.0',
    packages=find_packages(),
    author="Abhishek Kognole",
    author_email="aakognole@gmail.com",
    description="Simple program to quickly prepare PDB",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/aakognole/pdbprepper",
    project_urls={
        "Bug Tracker": "https://github.com/aakognole/pdbprepper/issues",
    },
    install_requires=[
        'pdbfixer',
        'openmm',
    ],
    python_requires=">=3.7",
    zip_safe=False,
    package_data={"PDBPrepper.data": ["*.dat"]},
    include_package_data=True,
)
