HPRFE2: High performance reduced finite element square
======================================================

HPRFE2 project consists of two parts: a set of utilities for bases generation and datases of a material (covered in this documentation), and the multiscale constitutive law for its use in FEA software (to be documented later).

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Quick start
===========

In this section we will follow a step-by-step process to review the several stages involved in the generation of a sample material dataset.
For this, we will create a project for the material ``COMPOSITE_01`` and sample a unit cell of this material using 9 trajectories.

To follow along, make sure the :ref:`installation <install>` is complete before proceding, no need to configure it just yet.

Let's get started.

Sampling
--------

Create root directory for our material::

  >>> mkdir COMPOSITE_01
  >>> cd COMPOSITE_01

Generate (and eventually edit) an initial configuration file::

  >>> python offline_common.py --init
  >>> ls
  configuration.json
  >>>

The configuration file contains the at least following needed parameters:

.. code-block:: json
  :emphasize-lines: 4,8

    {
      "config_data": {
        "cases_test_dataset": [ 5 ],
        "rve_data_points": [ 200, 400 ],
        "rve_data_points_range_list": [[100, 1600, 20], [1600, 2600, 100]],
        "rve_data_points_rom": true,
        "rve_data_modes": [20, 30, 40, 50, 60],
        "strain_svd_cutoff": 0.1,
      }
    }

At this point, we need a Kratos case already set up and working, which includes ``MainKratos.py``, ``model.mdpa`` (unit cell discretization), ``materials.json`` (COMPOSITE_01 constituve model and material parameters), and ``ProjectParameters.json`` (case configuration for Kratos).
For this tutorial we will use a test Kratos case bundled with the installation.

Create a ``training`` directory and populate it with our Kratos case::

  >>> pwd
  COMPOSITE_01
  >>> mkdir training
  >>> cd training
  >>> cp HPRFE2/test/template/* .

In the following step, we generate the sampling directories::

  >>> python3 offline_bases generate
  created sample_01
  created sample_02
  created sample_03
  created sample_04
  ...

This process create the file structure for each trajectory, and populates them with the Kratos case.
The only difference between them is the strain tensor value (in ``ProjectParameters.json``), which is taken from the file ``_training_strain_set.dat``.

At this point, we should have the following file structure (here showing only ``case_00``, as it is the same for the other directories)::

  PORTO_08
  ├── configuration.json
  └── training
      ├── case_00
      │   ├── MainKratos.py
      │   ├── materials.json
      │   ├── model.mdpa
      │   └── ProjectParameters.json
      ├── case_01
      ├── case_02
      ├── case_03
      ├── case_04
      ├── case_05
      ├── case_06
      ├── case_07
      ├── case_08
      ├── MainKratos.py
      ├── materials.json
      ├── model.mdpa
      ├── ProjectParameters.json
      ├── _training_strain_set.dat
      └── _training_strain_sets

We must now run every case.
In this tutorial, we just enter to each directory and run Kratos::

  >>> cd case_00
  >>> python3 MainKratos.py
  >>> cd ..

but in real-life cases we should have our own script for managing the jobs (more on this later).

Basis generation
----------------

- Create bases directory
- Generate bases 
- Generate datasets


.. _install:
Install
=======

HPRFE2 needs to be installed inside a previous Kratos installation.

.. todo::
        Make sure to add definition of required environment variables.

.. todo::
        Add reference to Kratos installation.

.. note::
        Assuming Kratos already installed

Basic steps::

   >>> cd /path/to/Kratos/applications
   >>> git checkout MultiscaleROMApplication
   >>> cd MultiscaleROMApplication


Configuration
=============

The default configuration file of each material is ``configuration.json``.
It must be present and located at the root directory of the material.

An initial configuration file with default parameters can be generated with 
>>> python offline_common.py --init

This step generates a ``configuration.json`` similar to::
  
    {
      "config_data": {
        "cases_test_dataset": [ 5 ],
        "rve_data_points": [ 200, 400 ],
        "rve_data_points_range_list": [[100, 1600, 20], [1600, 2600, 100]],
        "rve_data_points_rom": true,
        "rve_data_modes": [20, 30, 40, 50, 60],
        "strain_svd_cutoff": 0.1,
      }
    }
    
  

Usage
=====

Development
==========

.. automodule:: hprfe2


useful #1 -- auto members
=========================

This is something I want to say that is not in the docstring.

.. automodule:: hprfe2.useful_1
   :members:

useful #2 -- explicit members
=============================

This is something I want to say that is not in the docstring.

.. automodule:: hprfe2.useful_2
   :members: public_fn_with_sphinxy_docstring, _private_fn_with_docstring

.. autoclass:: MyPublicClass
   :members: get_foobar, _get_baz



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
