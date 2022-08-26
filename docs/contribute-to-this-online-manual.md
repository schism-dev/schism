# Overview
All SCHISM developers are encouraged to keep this online manual up to date with their code development.
This section serves as a brief guide for contributing to this manual.  
<br />

# Directory structure
All files used to build this manual are included in [SCHISM's Git Repository](https://github.com/schism-dev/schism).
As of April 9, 2021, You automatically get a copy of them after you download the SCHISM code:
```
git clone https://github.com/schism-dev/schism
```

For you to update this online manual, the relavant files are

- The mkdoc configuration file:
```
$your_dir/schism/mkdocs.yml
```

- The "docs" folder:
```
$your_dir/schism/docs/
```
, which mainly contains markdown files:
```
docs
├── schism
│   ├── barotropic-solver.md
│   ├── eulerian-lagrangian-method.md
│   ├── geometry-discretization.md
│   ├── momentum-equation.md
│   ├── overview.md
│   ├── physical-formulation.md
│   ...
└── verification-tests.md
    ...
```
and assets files (such as \*.png files)
```
docs
├── assets
│   ├── bay-delta-schism-grid.png
│   ├── bnd-elem.png
│   ├── case-study-guam.png
│   ...
```
<br />

# Updating the manual

## Minor edits

If you are only making minor changes on an existing manual page, the easiest way is probably directly editing the \*.md file on Github.
Simply click the "edit" button at the top-right corner of any manual page:

![](./assets/contribute_to_this-edit.png)

It will take you to the Github page, where you can directly edit the \*.md file then commit changes.
We do ask you to provide a brief and sensible commit message.


## Adding a short content
Follow these steps if your content is short enough to fit in a single page:

- Decide where your new page fits in this manual and specify it in the last section of "mkdocs.yml".
For example, see how "your-new-page" is inserted in "mkdocs.yml" below:
```
...
nav:
  - Home: index.md
  ...
  ...
  - Modules:
    - Overview: modules/overview.md
    - Generic tracer module: modules/generic-tracer.md
    - AGE: modules/age.md
    ...
    ...
  - Title of your new page: your-new-page.md
    ...
    ...
  - Contribute to this online manual: contribute-to-this-online-manual.md
```
Please give a meaningful name to the \*.md file associated with your new page.

- Put your \*.md file directly under "docs/", for example:
```
$your_dir/schism/docs/your-new-page.md
```
and put any additional materials in
```
$your_dir/schism/docs/assets/
```
Please give a sensible name to any additional materials, for example: "your-new-page-some-figure.png"
<br />

## Adding a content with multiple subsections
- Decide where your new content fits in this manual and specify it in the last section of "mkdocs.yml".
For example, see how "your-multi-page-content" is inserted in "mkdocs.yml" below:
```
...
nav:
  - Home: index.md
  ...
  ...
  - Modules:
    - Overview: modules/overview.md
    - Generic tracer module: modules/generic-tracer.md
    - AGE: modules/age.md
    ...
    ...
  - Title of your-multi-page-content:
    - Overview: your-multi-page-content/overview.md
    - Some topic: your-multi-page-content/some-topic.md
    ...
  - Contribute to this online manual: contribute-to-this-online-manual.md
```

- Make a subfolder under "docs", which will contain all new \*.md files
```
$your_dir/schism/docs/your-multi-page-content/
```

- For each subsection, create one \*.md file:
```
$your_dir/schism/docs/your-multi-page-content/overview.md
$your_dir/schism/docs/your-multi-page-content/some-topic.md
...
```

Although other directory structures can be used to make multi-page contents, we kindly ask you to follow our convention.
<br />
<br />

# Preview your edits
If you are making non-trivial changes, you may want to preview your edits before committing them to Git.
To do this, you will have to install "mkdocs" on your local machine.
This can be done system-wide using the system python installation (and associated pip command):
```bash
pip install mkdocs mkdocs_material mkdocs-with-pdf
```
Then, under
```
$your_dir/schism/
```
do
```
mkdocs serve
```
If there are no errors, it will return an IP address of the manual under editing at the end of the initial command line message:
```
INFO     -  Building documentation...
WARNING  -  without generate PDF(set environment variable ENABLE_PDF_EXPORT to 1 to enable)
INFO     -  Cleaning site directory
WARNING  -  Both index.md and readme.md found. Skipping readme.md from /media/feiye/My Book/schism/docs
WARNING  -  Documentation file 'contribute-to-this-online-manual.md' contains a link to '{{ site.baseurl }}{% link docs/compound-flood/compound-flood.md %}' which is not found in the documentation files.
WARNING  -  Documentation file 'contribute-to-this-online-manual.md' contains a link to '{{ site.baseurl }}{% link docs/acetool-tips/acetool-tips.md %}' which is not found in the documentation files.
INFO     -  Documentation built in 0.77 seconds
INFO     -  [22:21:41] Watching paths for changes: 'docs', 'mkdocs.yml'
INFO     -  [22:21:41] Serving on http://127.0.0.1:8000/
```
You can open it with a web browser and any changes in the source files will be reflected immediately.

See more instructions [here](https://github.com/schism-dev/schism/blob/master/docs/readme.md).
<br />
<br />

# Markdown
A quick reference of the [Markdown syntax](https://www.markdownguide.org/cheat-sheet/).
<br />
<br />

#Cross-referencing
Most equations, figures and papers in this document are labeled, so these can be referenced in another place in the
 same or different .md files. The general syntax is [] followed by (#<label>). Some examples are: 

- `[Zhang et al. 2016](#zhang2016)   (#zhang2016 in same .md)`
- `(Figure [5ab](./barotropic-solver.md#figure05))   (refer to figure05 inside barotropic-solver.md in the same dir)`

# More features
If you need advanced features for your tutorial page, you can use some HTML syntax in your \*.md files.
However, we recommend that you keep the page layout and syntax as simple as possible.


