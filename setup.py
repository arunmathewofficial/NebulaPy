from setuptools import setup, find_packages

with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()

setup(
    name='NebulaPy',
    description='',
    long_description=long_description,
    version='0.0.1',
    author='Arun Mathew',
    author_email='arun@cp.dias.ie',
    url='',
    packages=find_packages(),
    classifiers=[
        'Development Status :: Planning',
        'Intended Audience :: Science/Research',
        'Intended Audience :: End Users/Desktop',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    entry_points={
        'console_scripts': [
            'download_database=NebulaPy.src.database:download_database',
        ],
    },
)
