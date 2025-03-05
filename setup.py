from setuptools import setup, find_packages

with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()

setup(
    name='NebulaPy',
    description='',
    long_description=long_description,
    long_description_content_type="text/markdown",
    version='1.0.0-beta',
    author='Arun Mathew',
    author_email='arun@cp.dias.ie',
    url='',
    packages=find_packages(),

    include_package_data=True,  # This is crucial
    package_data={
        'NebulaPy': ['data/PoWR.tar.xz', 'data/CMFGEN.tar.xz', 'data/Chianti.tar.xz'],
    },

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: End Users/Desktop',
        "License :: OSI Approved :: MIT License",
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    entry_points={
        'console_scripts': [
            'download-database=NebulaPy.src.Database:DownloadDatabase.run',
        ],
    },
)
