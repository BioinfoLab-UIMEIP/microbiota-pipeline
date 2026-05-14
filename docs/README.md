# Documentation

This directory contains the GitHub Pages documentation for the 16S Microbiota Analysis Pipeline.

The documentation is built with Jekyll (minima theme) and is published at:

```
https://<your-github-username>.github.io/<repository-name>/
```

To browse the documentation online, visit the GitHub Pages URL above.

## Local preview

To preview the site locally, install Jekyll and run:

```bash
cd docs/
bundle exec jekyll serve
```

Then open http://localhost:4000 in a browser.

## Pages

- `index.md` — Landing page: overview, quick start, navigation
- `installation.md` — Installation guide (conda, R packages, databases)
- `configuration.md` — Configuration file reference
- `usage.md` — Running the pipeline, output directory structure, flags TSV format
- `methods.md` — Analytical framework description for each step
- `troubleshooting.md` — Common problems and solutions
- `_config.yml` — Jekyll site configuration
