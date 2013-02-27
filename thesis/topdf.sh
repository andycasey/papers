#!/bin/bash

# necessary to add more chapters to this list as you need them
export TEXINPUTS="`pwd`/introduction//:`pwd`/chapter1//:`pwd`/conclusions//:`pwd`/appendices//:`pwd`:"

TEX_CMD='pdflatex '
TEX_OPTS = ''

${TEX_CMD} ${TEX_OPTS} thesis.tex || exit 1
bibtex thesis || exit 1
${TEX_CMD} ${TEX_OPTS} thesis.tex || exit 1
${TEX_CMD} ${TEX_OPTS} thesis.tex || exit 1
${TEX_CMD} ${TEX_OPTS} thesis.tex || exit 1
${TEX_CMD} ${TEX_OPTS} thesis.tex || exit 1
rm *.aux
mv thesis.* output/.
mv output/thesis.tex .
cp output/thesis.pdf .

