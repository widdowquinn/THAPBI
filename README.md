# THAPBI

This repository contains code related to the LWEC/Tree Health/THAPBI project.

## Python style conventions

In this project, we're trying to keep to the Python [PEP8 style convention](https://www.python.org/dev/peps/pep-0008/), the [PEP257 docstring conventions](https://www.python.org/dev/peps/pep-0257/), and the [Zen of Python](https://www.python.org/dev/peps/pep-0020/). To help in this, a pre-commit hook script is provided in the `git_hooks` subdirectory that, if deployed in the `Git` repository, checks Python code for PEP8 correctness before permitting a `git commit` command to go to completion.

If the `pep8` module is not already present, it can be installed using `pip install pep8`

### Installing the `git hook`

To install the pre-commit hook:

1. clone the repository with `git clone https://github.com/widdowquinn/THAPBI` (you may already have done this)
2. change directory to the root of the repository with `cd THAPBI`
3. copy the pre-commit script to the `.git/hooks` directory with `cp git_hooks/pre-commit .git/hooks/`

### More information

* Git hooks (`git`): [https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks](https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks)
* Git hooks (tutorial): [http://githooks.com/](http://githooks.com/)
* PEP8: [https://www.python.org/dev/peps/pep-0008/](https://www.python.org/dev/peps/pep-0008/)
* PEP257: [https://www.python.org/dev/peps/pep-0257/](https://www.python.org/dev/peps/pep-0257/)
* Zen of Python (PEP20): [https://www.python.org/dev/peps/pep-0020/](https://www.python.org/dev/peps/pep-0020/)