#!/bin/bash

echo "In octave there might be variables printed just because a semicolon is forgotten"
echo "Use this script to find these culprits"

export LC_ALL=C
#export LC_CTYPE='en_US.ISO-8859'

find $1 -maxdepth 1 -type f -name *.m -print0 | xargs -0 cat | grep -v ';' | egrep -v '(for|end|%|function|if|else|case|switch|otherwise|do|until|while|\.\.\.|plot|xlabel|ylabel|figure|hold on|hold off|title|drawnow|colormap|tic|toc|exit|return|pause|addpath)' |  sed '/^\s*$/d'
