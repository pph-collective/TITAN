import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open('requirements.txt') as f:
    required = f.read().splitlines()

setuptools.setup(
    name="titan", # Replace with your own username
    version="2.0.0",
    # author="Example Author",
    # author_email="author@example.com",
    description="TITAN Agent Based Model",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/marshall-lab/TITAN",
    packages=setuptools.find_packages(),
    install_requires=required,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    package_data={'titan': ['params/*.yml', 'settings/*/*.yml']},
    entry_points={
        'console_scripts': ['run_titan=titan:script_init']
    }
)
