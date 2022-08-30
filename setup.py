#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = [
    "numpy",
    "matplotlib",
]

test_requirements = [ ]

setup(
    author="Skandan Chandrasekar",
    author_email='s39chand@uwaterloo.ca',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Simple 2D FDTD Solver written in Python",
    entry_points={
        'console_scripts': [
            'pyfdtd=pyfdtd.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    include_package_data=True,
    keywords='pyfdtd',
    name='pyfdtd',
    packages=find_packages(include=['pyfdtd', 'pyfdtd.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/SkandanC/pyfdtd',
    version='0.1.0',
    zip_safe=False,
)
