# HOA Library

#### The high order ambisonics library.

Sound space is one of the principal dimensions of the contemporary musical thought, especially in the electroacoustic music domain but also in intermedia arts. In this context, the [CICM](http://cicm.mshparisnord.org/ "CICM") has made spatialization its principal research axis. This project's aim is to give to musician spatialization models based on high order ambisonics and sound fields synthesis. This project is developed in a part of the Paris 8 University [LABEX arts H2H](http://www.labex-arts-h2h.fr/ "LABEX arts H2H"). You can visit the official website : [HoaLibrary](http://www.mshparisnord.fr/hoalibrary/ "HoaLibrary").

![Image](http://hoalibrary.mshparisnord.fr/wp-content/themes/hoa/images/hoa-icon03.svg "Hoa-Icon")

[![Travis](https://img.shields.io/travis/CICM/HoaLibrary-Light.svg?label=travis)](https://travis-ci.org/CICM/HoaLibrary-Light)
[![Appveyor](https://img.shields.io/appveyor/ci/CICM/HoaLibrary-Light.svg?label=appveyor)](https://ci.appveyor.com/project/CICM/HoaLibrary-Light/history)
[![Coverage Status](https://coveralls.io/repos/github/CICM/HoaLibrary-Light/badge.svg?branch=master)](https://coveralls.io/github/CICM/HoaLibrary-Light?branch=master)
[![Documentation](https://img.shields.io/badge/docs-doxygen-blue.svg)](http://cicm.github.io/HoaLibrary-Light/)

## Compilation :

Requires git and [CMake](https://cmake.org/).

```shell
# clone repository
$ git clone https://github.com/CICM/HoaLibrary-Light.git
$ cd HoaLibrary-Light
# init submodules
$ git submodule update --init --recursive
# generate project
$ mkdir build && cd build
$ cmake ..
# build project
$ cmake --build . --config Release
```

## Authors :

2012-2019: Pierre Guillot, Eliott Paris and [others](https://github.com/CICM/HoaLibrary-Light/graphs/contributors).

## Implementations :

- [HoaLibrary-Max](https://github.com/CICM/HoaLibrary-Max "Max") - HoaLibrary for Max.
- [HoaLibrary-Pd](https://github.com/CICM/HoaLibrary-PD "PD") - HoaLibrary for Pure Data.
- [HoaLibrary-Unity](https://github.com/CICM/HoaLibrary-Unity "Unity") - HoaLibrary for Unity.
- [HoaLibrary-Faust](https://github.com/CICM/HoaLibrary-Faust "Faust") - HoaLibrary for Faust.
- [OfxHoa](https://github.com/CICM/ofxHoa "Open Framework") - HoaLibrary for openFrameworks.
- [Cinder-Hoa](https://github.com/saynono/Cinder-Hoa) - HoaLibrary for Cinder.

> Warning: these implementations may not be up to date with this repo.

## ThirdParty

This repo is using the following third parties:

- [Eigen](http://eigen.tuxfamily.org)

#### License :

The HOA Library in under the [GNU](http://www.gnu.org/copyleft/gpl.html "GNU Public License"). If you'd like to avoid the restrictions of the GPL and use HOA Library for a closed-source product, you contact the [CICM](http://cicm.mshparisnord.org/ "CICM").
