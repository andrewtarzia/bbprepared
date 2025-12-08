# List all recipes.
default:
  @just --list

# Do a dev install.
setup:
  uv sync --all-extras --dev

# Run code checks.
check:
  #!/usr/bin/env bash

  error=0
  trap error=1 ERR

  echo
  ( set -x; uv run ruff check src tests examples )

  echo
  ( set -x; uv run ruff format --check src tests examples )

  echo
  ( set -x; uv run mypy src examples tests )

  echo
  ( set -x; uv run pytest --cov=src --cov-report term-missing )

  echo
  ( set -x; uv run make -C docs doctest )

  test $error = 0

# Auto-fix code issues.
fix:
  uv run ruff format .
  uv run ruff check --fix .

# Build a release.
build:
  uv build

# Build docs.
docs:
  rm -rf docs/source/_autosummary
  uv run make -C docs html
  echo Docs are in $PWD/docs/build/html/index.html
