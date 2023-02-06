# PRODIGY Development

## Installation

We use `poetry` to manage the dependencies and the virtual environment, so you need to install it first; check the [official documentation](https://python-poetry.org/docs/#installation) for more details.

Clone the repository and install the dependencies:

```text
git clone https://github.com/haddocking/prodigy.git && cd prodigy
poetry install
```

## Testing

To run the tests, use the following command:

```text
python -m unittest
```

## Code style

We use `trunk` as the "all-purpose" linting tool, check its [documentation](https://docs.trunk.io/docs/install).

To check for code style issues, run:

```text
trunk check
```

To automatically fix the issues, run:

```text
trunk fmt
```
