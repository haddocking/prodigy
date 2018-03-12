import os
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

requirements = [
    "numpy",
    "biopython",
]

setup(
    name='prodigy',
    version = "2.1.0",
    description=("PROtein binDIng enerGY prediction."),
    url='http://github.com/haddocking/prodigy',
    author='Computational Structural Biology Group @ Utrecht University',
    author_email='prodigy.bonvinlab@gmail.com',
    license='Apache 2.0',
    packages=['prodigy', 'prodigy.lib'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: Apache Software License",
    ],
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'prodigy = prodigy.predict_IC:main',
        ]
    },
    zip_safe=False
)
