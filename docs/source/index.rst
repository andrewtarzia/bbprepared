.. toctree::
   :hidden:
   :caption: bbprepared
   :maxdepth: 2

   Recipes <recipes>
   Generators <generators>
   Containers <containers>
   Selectors <selectors>
   Modifiers <modifiers>
   Processing ensembles <processes>
   New building blocks <building_blocks>

.. toctree::
  :hidden:
  :maxdepth: 2
  :caption: Modules:

  Modules <modules>

============
Introduction
============

| GitHub: https://www.github.com/andrewtarzia/bbprepared

.. figure:: _static/logo.png

:mod:`bbprepared` or ``bbprep`` is a toolkit aimed at simplifying the
preparation of your building blocks for
`stk <https://stk.readthedocs.io/en/stable/>`_ construction and analysis.


Installation
============

To get :mod:`.bbprepared`, you can install it with pip:

.. code-block:: bash

  pip install bbprepared

Developer Setup
---------------

To develop with :mod:`bbprepared`, you can clone the repo and use
`just <https://github.com/casey/just>`_ and `uv <https://docs.astral.sh>`_
to setup the dev environment:

.. code-block:: bash

  just setup


Examples
========

* See `Recipes` for usage and the code tests include additional examples.
* See `gists <https://gist.github.com/andrewtarzia>`_ for usage.


How To Cite
===========

If you use ``bbprepared`` please mention the URL

  https://github.com/andrewtarzia/bbprepared

Or use the `citation file <https://github.com/andrewtarzia/bbprepared/blob/main/CITATION.cff>`_.



Acknowledgements
================

Funded by the European Union - Next Generation EU, Mission 4 Component 1
CUP E13C22002930006.

This work is a mixture of codes developed throughout my postdoc in the
`Jelfs Group <http://www.jelfs-group.org/>`_, and during my time as a developer
of `stk <https://stk.readthedocs.io/en/stable/>`_ and
`stko <https://github.com/JelfsMaterialsGroup/stko>`_ with Lukas Turcani.
