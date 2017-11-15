import os
from setuptools import setup, find_packages

def setup_package():
    setup(
        name = "themisra",
        version = "0.1.1",
        author = "Jay Laura",
        author_email = "jlaura@usgs.gov",
        description = ("Thermal modeling using THEMIS temperature data."),
        license = "Public Domain",
        keywords = "geophysical, modeling",
        url = "",
        packages=find_packages(),
        include_package_data=True,
        zip_safe=False,
        install_requires=['pysis'],
        entry_points={'console_scripts': ['interpolate_krc=krc.interpolate_krc_cluster:main',
                                          'themis_to_brightness=krc.themis_to_brightness:main']},
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Topic :: Utilities",
            "License :: Public Domain",
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
        ],
    )

if __name__ == '__main__':
    setup_package()
