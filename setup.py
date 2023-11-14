import os.path

from setuptools import setup

root_path = os.path.abspath(os.path.join(__file__, '..'))

setup(
	name='conc-plots',
	packages=['conc_plots'],
	package_dir={'conc_plots': root_path},
	python_requires='>2.6,<3',
	install_requires=[
        'numpy>=1.14.2',
        'pandas>=0.22.0',
        'matplotlib>=2.2.2',
    ],
)
