import nbformat as nbf
from glob import glob

# Collect a list of all notebooks in the content folder
notebooks = glob("./Chemistry_Lessons_Jupyter_Notebooks/**/**/*.ipynb", recursive=True)
# Text to look for in adding tags
text_search_dict = {
    "# HIDE CODE": "hide_input"  # Hide the input w/ a button to show
}

# Search through each notebook and change hide_input tag to hide-input
for ipath in notebooks:
    ntbk = nbf.read(ipath, nbf.NO_CONVERT)

    for cell in ntbk.cells:
        cell_tags = cell.get('metadata', {}).get('tags', [])
        for i, tag in enumerate(cell_tags):
            if tag == "hide_input":
                cell_tags[i] = "hide-input"
        if len(cell_tags) > 0:
            cell['metadata']['tags'] = cell_tags

    nbf.write(ntbk, ipath)
