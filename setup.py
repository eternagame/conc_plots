import os.path

from setuptools import setup

root_path = os.path.abspath(os.path.join(__file__, '..'))

setup(
	name='conc-plots',
	packages=['conc_plots'],
	package_dir={'conc_plots': root_path},
	python_requires='>3.10,<4',
	install_requires=[
        'numpy>=1.14,<2',
        'pandas>=0.22.0,<3',
        'matplotlib>=3.1,<4',
    ],
)
