Low-Frequency Radio Sky Simulation Tools
========================================

Jingying WANG, Junhua GU, and Weitian LI


Introduction
------------
This repository contains the tools that were used to simulate the
low-frequency radio sky images for use in
[Wang et al. 2010](http://adsabs.harvard.edu/abs/2010ApJ...723..620W).
These tools were written by Jingying WANG and Junhua GU.
A significant part of the code has been modified by Weitian LI to
make them easier to understand and to use.


Repo Structure
--------------
* `cluster`:
  Original code written by WANG and GU to simulate the radio halos
  in galaxy clusters.
* `coscalc`:
  Cosmology calculator.
* `fio`:
  A library for easier manipulate FITS files in C++.
* `halos`:
  New code written by LI to re-implement the simulation method of
  radio halos (not finished ...).
* `ica`:
  ICA tools to separate radio halos from other foreground components.
* `pointsource`:
  Tools to simulate various types of point sources.
* `scripts`:
  Scripts to simulate the Galactic emission maps, to help analyze
  resutls, etc.
