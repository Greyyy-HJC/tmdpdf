# TMDPDF Documentation

This is TMDPDF project from LPC group, SJTU.
Documentation of Transverse Momentum Dependent Parton Distribution Function (TMDPDF).

## Installation

install for development

```shell
git clone `url to this repo`
cd tmdpdf
pip install -r requirements.txt
pip install -e .
```

or install for usage

```shell
pip install git+`url to this repo`
```

## Run Example Scripts

```shell
python3 scripts/main.py
```

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.
    libtmdpdf/
        ...
    scripts/      # entrypoints

## Mkdocs Commands

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs -h` - Print help message and exit.
