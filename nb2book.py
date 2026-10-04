#!/usr/bin/env python3
"""
Build the book PDF from the notebooks in notebooks/.

Adapted from the CS2101 production line (~/cs2101/jupyter/nb2book.py).

For each chapter notebook:
  1. Execute all cells in memory (Julia kernel named in the notebook).
     The notebooks on disk are never modified: they are kept without outputs.
  2. Turn the title cell into a plain chapter heading: drop the repeated
     course title and byline, and the "N. " in front of the chapter title
     (LaTeX numbers the chapters itself).
  3. Move the exercise solutions -- code cells tagged hide-cell and the
     <details> proof blocks -- into an appendix, under their exercise
     numbers, leaving a page reference behind.
  4. Export to LaTeX via nbconvert, with
       - templates/book: a Unicode monospace font, and coloured boxes for
         the notebooks' <div class="alert alert-KIND"> blocks,
       - templates/book/alerts.lua: the pandoc filter that produces them.

Then:
  5. Assemble _book_build/book.tex with \\documentclass{book}, one chapter
     per notebook and the Solutions appendix, compile with xelatex (two
     passes for the table of contents and page references), and copy the
     result to exports/book.pdf.

Usage:  python3 nb2book.py
"""

import os, re, shutil, subprocess
from copy import deepcopy
from pathlib import Path

import nbformat
from nbconvert.exporters import LatexExporter
from nbconvert.filters import pandoc as _pandoc_mod
from nbconvert.preprocessors import ExecutePreprocessor

ROOT = Path(__file__).parent.resolve()
NB_DIR = ROOT / 'notebooks'
BUILD = ROOT / '_book_build'
OUTPUT = ROOT / 'exports' / 'book.pdf'
TEMPLATE_DIR = ROOT / 'templates'
TEMPLATE_NAME = 'book'
LUA_FILTER = TEMPLATE_DIR / TEMPLATE_NAME / 'alerts.lua'

TITLE = 'Computational Aspects of Complex Reflection Groups'
AUTHOR = r'Götz Pfeiffer \\ University of Galway'

CHAPTERS = ['orbits', 'coxeter', 'enumerate', 'linear']

CELL_TIMEOUT = 1800   # seconds; the E7 conjugacy classes take a while


# ── pandoc: alert boxes, and lists right after a paragraph ────────────────────
# nbconvert renders every markdown cell through convert_pandoc.  Add the
# alert-box filter, and enable lists_without_preceding_blankline so that a
# list starting on the line right after a paragraph renders as a list, the
# way Jupyter renders it.
_orig_convert_pandoc = _pandoc_mod.convert_pandoc

def _convert_pandoc(source, from_format, to_format, extra_args=None):
    if 'markdown' in from_format:
        from_format += '+lists_without_preceding_blankline'
        extra_args = list(extra_args or []) + [f'--lua-filter={LUA_FILTER}']
    return _orig_convert_pandoc(source, from_format, to_format, extra_args=extra_args)

_pandoc_mod.convert_pandoc = _convert_pandoc


# ── notebook transformations ──────────────────────────────────────────────────
TITLE_LINE = re.compile(r'^#\s+\d+\.\s+(.*)$', re.MULTILINE)
EXERCISE = re.compile(r'^\*\*Exercise (\d+\.\d+)')
JULIA_LOGO = re.compile(r'<img src="images/julia\.png"[^>]*>')
DETAILS = re.compile(r'<details>\s*<summary>(.*?)</summary>\s*', re.DOTALL)
IMAGE_LINE = re.compile(r'(?m)^(!\[[^\]]*\]\([^)]*\))\n(?=!\[)')


def chapter_title(nb):
    """Reduce the title cell to "# Title"; return the title."""
    first = nb.cells[0]
    m = TITLE_LINE.search(first.source)
    assert first.cell_type == 'markdown' and m, 'no "# N. Title" in the first cell'
    first.source = '# ' + m.group(1).strip()
    return m.group(1).strip()


def is_solution(cell):
    if cell.cell_type == 'code':
        return 'hide-cell' in cell.metadata.get('tags', [])
    return '<details>' in cell.source


def unhide(cell):
    """A solution as it appears in the appendix: nothing collapsed."""
    cell = deepcopy(cell)
    if cell.cell_type == 'code':
        cell.metadata.pop('jupyter', None)
        cell.metadata.pop('tags', None)
    else:
        summary = DETAILS.search(cell.source)
        label = re.sub(r'</?b>', '**', summary.group(1)) if summary else ''
        cell.source = DETAILS.sub(label + '\n', cell.source).replace('</details>', '')
    return cell


def split_solutions(nb):
    """Move solution cells out of nb; return [(exercise, [cells])]."""
    kept, solutions, exercise, noted = [], [], None, set()
    for cell in nb.cells:
        if cell.cell_type == 'markdown':
            m = EXERCISE.match(cell.source)
            if m:
                exercise = m.group(1)
        if not is_solution(cell):
            kept.append(cell)
            continue
        assert exercise, 'solution cell before the first exercise'
        if exercise not in noted:
            noted.add(exercise)
            solutions.append((exercise, []))
            kept.append(nbformat.v4.new_markdown_cell(
                f'*Solution: page \\pageref{{sol:{exercise}}}.*'))
        solutions[-1][1].append(unhide(cell))
    nb.cells = kept
    return solutions


def fix_markdown(nb):
    """Adjust the markdown for print.

    - A small julia logo in the tip boxes (pandoc would drop the <img>).
    - Images on consecutive lines form one paragraph, and so sit side by
      side, which is too wide for the page: give each its own paragraph.
    """
    for cell in nb.cells:
        if cell.cell_type == 'markdown':
            cell.source = JULIA_LOGO.sub('![](images/julia.png){width=0.6cm}', cell.source)
            cell.source = IMAGE_LINE.sub(r'\1\n\n', cell.source)


def fix_outputs(nb):
    """Make the outputs printable.

    - Drop SVG where a PNG exists: nbconvert prefers SVG for LaTeX, but
      converting it needs inkscape.  The graph plots come with both.
    - Julia elides long matrices with a vertical ellipsis, which the
      monospace font lacks: use a colon.
    """
    for cell in nb.cells:
        for out in cell.get('outputs', []):
            data = out.get('data', {})
            if 'image/png' in data:
                data.pop('image/svg+xml', None)
            if 'text/plain' in data:
                data['text/plain'] = data['text/plain'].replace('⋮', ':')
            if 'text' in out:
                out['text'] = out['text'].replace('⋮', ':')


def solutions_notebook(chapters, metadata):
    """The appendix, as one more notebook: a section per chapter."""
    cells = [nbformat.v4.new_markdown_cell('# Solutions')]
    for title, solutions in chapters:
        if not solutions:
            continue
        cells.append(nbformat.v4.new_markdown_cell(f'## {title}'))
        for exercise, sol_cells in solutions:
            cells.append(nbformat.v4.new_markdown_cell(
                f'### Exercise {exercise}\n\n\\label{{sol:{exercise}}}'))
            cells.extend(sol_cells)
    nb = nbformat.v4.new_notebook(cells=cells)
    nb.metadata = deepcopy(metadata)
    return nb


# ── LaTeX helpers (as in CS2101) ──────────────────────────────────────────────
def extract_preamble(tex):
    return tex[:tex.find(r'\begin{document}')].rstrip()

def extract_body(tex):
    m = re.search(r'\\begin\{document\}(.*?)\\end\{document\}', tex, re.DOTALL)
    return re.sub(r'^\\maketitle\s*\n?', '', m.group(1).strip(), flags=re.MULTILINE)

def adapt_preamble(preamble):
    """article -> book; drop the per-notebook title, author and date."""
    preamble = re.sub(r'\\documentclass(\[.*?\])?\{article\}', r'\\documentclass\1{book}', preamble)
    return re.sub(r'\\(title|author|date)\{[^}]*\}\n?', '', preamble)

def demote_headings(body):
    """section -> chapter, subsection -> section, subsubsection -> subsection."""
    mapping = {'subsubsection': 'subsection', 'subsection': 'section', 'section': 'chapter'}
    return re.sub(r'\\(subsubsection|subsection|section)(\*?)\{',
                  lambda m: '\\' + mapping[m.group(1)] + m.group(2) + '{', body)


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    BUILD.mkdir(exist_ok=True)
    images = BUILD / 'images'
    if images.exists():
        shutil.rmtree(images)
    shutil.copytree(NB_DIR / 'images', images)

    exporter = LatexExporter(extra_template_basedirs=[str(TEMPLATE_DIR)],
                             template_name=TEMPLATE_NAME)
    exporter.environment.filters['convert_pandoc'] = _convert_pandoc

    def to_latex(nb, name):
        tex, resources = exporter.from_notebook_node(nb)
        # nbconvert names output figures by cell index alone: prefix them with
        # the notebook name, so that the chapters don't overwrite each other's.
        for fname, data in resources.get('outputs', {}).items():
            unique = f'{name}_{fname}'
            (BUILD / unique).write_bytes(data)
            tex = tex.replace('{' + fname + '}', '{' + unique + '}')
        # pandoc labels each heading after its text, so "Setup" or "Exercises"
        # would be defined once per chapter: prefix them with the notebook name.
        # (Only the sol: labels are referred to across notebooks.)
        tex = re.sub(r'\\(?:label|ref)\{(?!sol:)|\\hyperref\[(?!sol:)',
                     lambda m: m.group(0) + name + ':', tex)
        return tex

    preamble, bodies, appendix, metadata = None, [], [], None
    for name in CHAPTERS:
        print(f'[nb2book] {name}', flush=True)
        nb = nbformat.read(NB_DIR / f'{name}.ipynb', as_version=4)
        print('  executing ...', flush=True)
        ExecutePreprocessor(timeout=CELL_TIMEOUT).preprocess(
            nb, {'metadata': {'path': str(NB_DIR)}})
        fix_outputs(nb)
        title = chapter_title(nb)
        appendix.append((title, split_solutions(nb)))
        fix_markdown(nb)
        metadata = metadata or nb.metadata
        print('  exporting to LaTeX ...', flush=True)
        tex = to_latex(nb, name)
        preamble = preamble or extract_preamble(tex)
        bodies.append(demote_headings(extract_body(tex)))

    print('[nb2book] solutions', flush=True)
    sol = solutions_notebook(appendix, metadata)
    fix_markdown(sol)
    bodies.append('\\appendix')
    bodies.append(demote_headings(extract_body(to_latex(sol, 'solutions'))))

    print('[nb2book] assembling book.tex ...', flush=True)
    (BUILD / 'book.tex').write_text('\n'.join([
        adapt_preamble(preamble),
        f'\\title{{{TITLE}}}',
        f'\\author{{{AUTHOR}}}',
        '\\date{}',
        '',
        '\\begin{document}',
        '\\maketitle',
        '\\tableofcontents',
        '',
        '\n\n'.join(bodies),
        '',
        '\\end{document}',
    ]))

    # Two passes: the table of contents and the \pageref's need the .aux file.
    for n in (1, 2):
        print(f'[nb2book] xelatex pass {n} ...', flush=True)
        r = subprocess.run(['xelatex', '-interaction=nonstopmode', 'book.tex'],
                           cwd=BUILD, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout[-3000:])
            raise SystemExit('[nb2book] xelatex failed -- see output above')

    OUTPUT.parent.mkdir(exist_ok=True)
    shutil.copy(BUILD / 'book.pdf', OUTPUT)
    print(f'[nb2book] done -> {OUTPUT.relative_to(ROOT)}', flush=True)


if __name__ == '__main__':
    main()
