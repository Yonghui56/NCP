#!/bin/sh

# Script for GNU build system
# Use this to rebuild configure script

autoheader
aclocal
libtoolize
aclocal
automake
autoconf

