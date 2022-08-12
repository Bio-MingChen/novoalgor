# novoalgor

A series of tools for single cell assignment or cell annotation

Document to see http://39.102.38.22:8090/page/display?document_id=161

## Available tools

```shell
$novoalgor --help
Get into miniconda scanpy_env environment successfully!
Usage: novoalgor.py [OPTIONS] COMMAND [ARGS]...

  Tools for cell marker search ,celltype assignment and visualization

Options:
  --help  Show this message and exit.

Commands:
  cellvisual      visualization for celltype assignment
  fast-celltype   run fast-celltype algorithm to celltype assignment
  get-top-gene    get top n markers from marker file
  global-filter   filter specific markers from input marker file globally
  global-search   search and filter specific markers globally
  local-search    search and filter specific markers by input celltype...
  show-precision  Compute precision for prediction results and output...
```