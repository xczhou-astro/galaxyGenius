from setuptools import setup, find_packages

with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()

setup(
    name="GalaxyGenius",
    version="1.0.0",
    description="A Mock Galaxy Image Generator for Various Telescopes from Hydrodynamical Simulations",
    author="Xingchen Zhou et al.",
    author_email="xczhou95@gmail.com",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/xczhou-astro/galaxyGenius",
    packages=find_packages(),
    python_requires=">=3.11.10",
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3.11.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
