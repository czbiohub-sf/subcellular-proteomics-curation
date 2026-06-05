from setuptools import setup

with open("requirements.txt") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.lstrip().startswith("#")]

setup(
    name="cellxgene-schema",
    version="7.0.0",
    url="https://github.com/czbiohub-sf/subcellular-proteomics-curation",
    license="MIT",
    author="Chan Zuckerberg Initiative",
    author_email="cellxgene@chanzuckerberg.com",
    description="Adjusted cellxgene-schema for validating and labeling subcellular proteomics datasets",
    long_description="Adjusted cellxgene-schema for validating and labeling subcellular proteomics datasets",
    install_requires=requirements,
    python_requires=">=3.10",
    packages=["cellxgene_schema"],
    package_dir={"cellxgene_schema": "cellxgene_schema"},
    package_data={"cellxgene_schema": ["uniprot_files/*.tsv.gz", "uniprot_files/*.yml", "schema_definitions/*yaml"]},
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
