from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='PDBPrepper',
    version='0.1.0',
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
    entry_points={
        'console_scripts': [
            'pdbprepper=pdbprepper.core:main',
        ],
    },
    python_requires=">=3.6",
    zip_safe=False,
)
