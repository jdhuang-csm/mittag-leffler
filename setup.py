from setuptools import setup

with open("README.md", "r",encoding="utf8") as fh:
    long_description = fh.read()

setup(name='Mittag-Leffler',
	version='0.1',
	description='Python implementation of the Mifflag-Leffler function',
	long_description=long_description,
    long_description_content_type="text/markdown",
	url='https://github.com/jdhuang-csm/mittag-leffler',
	author='Jake Huang',
	author_email='jdhuang@mines.edu',
	license='BSD 3-clause',
	packages=['mitlef'],
	install_requires=[
		'numpy',
		'scipy'
		],
	include_package_data=True
	)