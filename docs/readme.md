# Documenting on mkdocs

## Installation (optional; for local development)
The basic install only requires `mkdocs` and `mkdocs_materials` to be installed. This can be done system-wide using the system python installation (and associated `pip` command). 

```bash
pip install mkdocs mkdocs_material mkdocs-with-pdf
```

To check where `pip` is located, use `which pip` in linux/unix system.

If you are a conda/anaconda user, you can also create a new environment from the provided `environment.yml` file

## Writing markdown document
**Basic text**

Sections, subsections: `#`, `##` etc.

Bold: `**Text**` $\rightarrow$ **Text**

Italic: `_Text_` $\rightarrow$ _Text_

**Image** 

```markdown
<figure markdown id='figureXX'>
![alt-text](relative/path/to/figure.png)
<figcaption>A good figure caption</figcaption>
</figure>
```

To reference this figure in the text, create a link `[Figure XX](path/to/page.md#figureXX)`. If the figure is referenced in the same page, the `path/to/page.md` can be omitted.

**Table**

```markdown
| header 1| header 2|
|---------|---------|
| item 1 | item 2|
```
will produce - 

| header 1| header 2|
|---------|---------|
| item 1 | item 2|

Add a `span` before for a referencable caption

```
<span id="table01">Table caption</span>
```

and refer using link notation `Table [1](#table01)`.

**Equations**

Equations can be added using standard latex `\begin{equation}\label{eqxx}\end{equation}` block, or inline with `$...$` block. Other latex command works inside these blocks. The equation inside `equation` block with label `eqxx` can be referred inline using `$\ref{eqxx}$`.

**References**

References can be added as link. For example `[Zhang et al. 2016](#zhang2016)`. An entry for this reference cab be added with html `<span>` in the end of each page as `<span id="zhang2016">Zhang, Y., Stanev, E.V. and S. Grashorn (2016) Unstructured-grid model for the North Sea and Baltic Sea: validation against observations, Ocean Modelling, 97, 91-108.</span>`. Here `id` is what makes the link works.

## Build and deploy site
```bash
mkdocs serve # to see locally
mkdocs build # build site
mkdocs gh-deploy --force # serve site to github on gh-page branch
```

You will need github pages fully setup. See [Github Pages](https://pages.github.com/) for more information.
