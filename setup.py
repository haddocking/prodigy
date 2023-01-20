import codecs
import os

from setuptools import setup


def read(fname):
    return codecs.open(
        os.path.join(os.path.dirname(__file__), fname), encoding="utf-8"
    ).read()


with open("requirements.txt") as f:
    required = f.read().splitlines()

setup(
    name="prodigy",
    version="2.1.1",
    description=("PROtein binDIng enerGY prediction."),
    url="http://github.com/haddocking/prodigy",
    author="Computational Structural Biology Group @ Utrecht University",
    author_email="prodigy.bonvinlab@gmail.com",
    license="Apache 2.0",
    packages=["prodigy", "prodigy.lib"],
    package_dir={"prodigy": "prodigy"},
    package_data={
        "prodigy": [
            "naccess.config",
        ]
    },
    long_description=read("README.md"),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10"
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.7, <3.11",
    install_requires=required,
    entry_points={
        "console_scripts": [
            "prodigy = prodigy.predict_IC:main",
        ]
    },
    zip_safe=True,
)
