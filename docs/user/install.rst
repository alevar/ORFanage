.. _install:

Installation of ORF-an
========================


Using Pre-Built Binaries
-------------------------------

Several pre-compiled binaries are pre-built for every release of ORFan to make it easy to use with no additional installation requirements.

Latest stable binary release can be obtained using the following command

- On Linux: ::

	$ wget https://github.com/alevar/ORFanage/releases/latest/download/orfanage

- On MacOS: ::

	$ wget https://github.com/alevar/ORFanage/releases/latest/download/orfanage
	
Alternatively, the latest binaries may be found at the latest release page: https://github.com/alevar/ORFanage/releases/latest

Building from source
-------------------------------

If you want to build it from source, we recommend cloning the git repository as shown below.

::

    $ git clone https://github.com/alevar/ORFanage.git --recursive
    $ cd ORFanage
    $ cmake -DCMAKE_BUILD_TYPE=Release -G "Unix Makefiles" .
    $ make -j4

For a fully static build ``-DORFANAGE_STATIC_BUILD=1`` needs to be added to the list of arguments in the cmake command.

By default make install will likely require administrative privileges. To specify custom installation path ``-DCMAKE_INSTALL_PREFIX=<custom/installation/path>`` needs to be added to the list of arguments in the cmake command.

If you are using a very old version of Git (< 1.6.5) the flag ``--recursive`` does not exist. In this case you need to clone the submodule separately via ``git submodule update --init --recursive``.

**Requirements**

Operating System
GNU/Linux

Architecture
Intel/AMD platforms that support ``POPCNT``

Compiler
GCC ≥ 4.9, Clang ≥ 3.8

Build system
CMake ≥ 3.2

Language support
C++14
