:maintainers:
  `andrewtarzia <https://github.com/andrewtarzia/>`_
:documentation: https://bbprepared.readthedocs.io/en/latest/


============
Introduction
============

``bbprepared`` is a toolkit aimed at simplifying the
preparation of your building blocks for
`stk <https://stk.readthedocs.io/en/stable/>`_ construction and analysis.



Installation
============

To get :mod:`.bbprepared`, you can install it with pip. However, there is a
distinction based on if you want the "lite" versions of ``stk`` and ``stko``.
This command will just give you ``bbprepared``, but will not provide all the
dependancies:

.. code-block:: bash

  pip install bbprepared

You can then add ``stk``/``stko``:

.. code-block:: bash

  pip install stk stko

Or ``stk-lite``/``stko-lite``:

.. code-block:: bash

  pip install stk-lite stko-lite

You can do this in one step as:

.. code-block:: bash

  pip install "bbprepared[full]" # For stk/stko
  pip install "bbprepared[lite]" # For stk-lite/stko-lite

.. warning::

  You should not mix ``stk``/``stko`` and ``stk-lite``/``stko-lite`` and
  ``bbprepared`` will raise an error if you do.

  Fix: "pip uninstall -y stk stk-lite stko stko-lite && " then
  "pip install 'bbprepared[full]'" (or 'bbprepared[lite]')


Developer Setup
---------------

To develop with :mod:`bbprepared`, you can clone the repo and use
`just <https://github.com/casey/just>`_ and `uv <https://docs.astral.sh>`_
to setup the dev environment with ``stk-lite``/``stko-lite``:

.. code-block:: bash

  just setup




Examples
========

* See ``examples/`` directory for some examples.
* See ``Recipes`` for usage and the code tests include additional examples.
* See `gists <https://gist.github.com/andrewtarzia>`_ for usage.
* In particular `this tutorial <https://gist.github.com/andrewtarzia/23e7843365bad76277322e207d95e627>`_

How To Cite
===========

If you use ``bbprepared`` please mention the URL

  https://github.com/andrewtarzia/bbprepared

Or cite using thr ``CITATION.cff`` file.


Acknowledgements
================

Funded by the European Union - Next Generation EU, Mission 4 Component 1
CUP E13C22002930006.

This work is a mixture of codes developed throughout my postdoc in the
`Jelfs Group <http://www.jelfs-group.org/>`_, and during my time as a developer
of `stk <https://stk.readthedocs.io/en/stable/>`_ and
`stko <https://github.com/JelfsMaterialsGroup/stko>`_ with Lukas Turcani.
