This is a fork of MOODS version 1.0.1 it adds an additionaly python object
`MOODSSearch` that lets you run many search with minimal parameter 
reloading overhead so that you don't need to concatenate searches 
and then parse the results to get speed ups with `MOODS.search`

# MOODS version 1.0.1
=======================

MOODS is a suite of algorithms for matching position weight matrices (PWM) against DNA sequences. It features advanced matrix matching algorithms implemented in C++ that can be used to scan hundreds of matrices against chromosome-sized sequences in few seconds.

MOODS has been designed to be used as a library wherever PWM matching is needed. It can be used as standalone analysis tool or as a component in larger programs. It contains interfaces for BioPerl and Biopython toolkits. MOODS can thus be easily called from C++, Python and Perl programs.

MOODS is licenced under GPL version 3 license and under the Biopython license, so you may use it under terms of either of those licenses. See COPYING.GPLv3 and COPYING.BIOPYTHON files for license details.

The project web page is at http://www.cs.helsinki.fi/group/pssmfind. Please refer there for further information and documentation.

moods-py.cc is derived from Moods 1.0.1 source.  This derivation is
Copyright (C) 2014  Nick Conway Wyss Institute for Biologically Inspired Engineering
and is GPL version 3 license 

## CONTACT

MOODS has been written by Pasi Rastas, Janne Korhonen and Petri Martinmäki. The project is currently maintained by Janne Korhonen (jazkorho@cs.helsinki.fi).

for questions regarding the MOODS.moody.MOODSSearch and moods-py.cc please use issues on github

## INSTALLATION

You need to compile the C++ library before installing the Perl and Python interfaces. This can be done by running "make" in the "src" directory.

Installing extensions:
- For perl extension see perl/README.
- For python extension see python/README.md
