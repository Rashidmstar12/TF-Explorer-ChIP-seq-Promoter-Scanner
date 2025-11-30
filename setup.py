from setuptools import setup, find_packages

setup(
    name="tf-explorer",
    version="1.1.0",
    author="Rashidmstar12",
    author_email="rashidmstar@pondiuni.ac.in",
    description="A tool to explore TF binding sites for any gene using ENCODE and JASPAR.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Rashidmstar12/TF-Explorer-ChIP-seq-Promoter-Scanner",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    install_requires=[
        "requests",
        "pandas",
        "biopython",
        "matplotlib",
        "seaborn",
        "mygene",
        "pyyaml",
        "streamlit"
    ],
    entry_points={
        "console_scripts": [
            "tf-explorer=tf_explorer.cli:main",
        ],
    },
)
