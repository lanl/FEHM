---
title : Installation Instructions
layout : page_getting-started
permalink: /getting-started/installation
hero_height: is-hidden
---

# Installation Instructions

## Install from Git (Recommended)

1. Obtain the FEHM repositoryfrom Github. To obtain, type the following command into a terminal:
```
git clone https://github.com/lanl/FEHM.git
```

1. Build FEHM. In a terminal, navigate to FEHM/src and type the following command:
```
make
```


## Creating the FEHM binary from source (UNIX)

On the system where FEHM is to be installed, make an installation directory, with subdirectories src and objects:

```

  mkdir fehm
  mkdir src objects
```

Copy all fehm source files (i.e., extract them from a tar file -- ``fehm_src.tar``) into the src directory. Obtain source files from [https://github.com/lanl/FEHM/releases/latest](https://github.com/lanl/FEHM/releases/latest):

```
  cd fehm/src
  tar xvf fehm_src.tar
```
A Makefile is included and should be placed in your objects directory. To compile and link FEHM, change into the objects directory and compile the code:

```
  cd fehm/objects
  make -OR- make -f Makefile
```
The makefile creates an executable called:

```
  xfehm_v2.30
```
It should be noted that FEHM uses the GZSOLVE Application (ECD-97) reuse components, solve_new, solve_rdof, and slvesu. The GZSOLVE subroutines are compiled directly into this version of FEHM.

## Installation Verification and Validation

A series of test scripts have been developed to automate the validation procedure for FEHM. They are described in more detail in the FEHM-VTP.APP, of the Validation Test Plan for the FEHM Application Version 2.30 (10086-VTP-2.30-00). See the [FEHM VTP](https://www.lanl.gov/orgs/ees/fehm/docs/FEHM_VERIFICATION_V3.3.0.pdf) for a discussion of the tests performed and their results.

In addition to an extensive test suite, FEHM also has a new simplified ``FehmPyTests`` test suite for FEHM developers to easily test new code and add new tests for existing or future functionality. FehmPyTests uses the Python unit test module and a general test method that can be called for each new test case. More information on the simplified test suite can be found [**Here**](fehmpytests).
