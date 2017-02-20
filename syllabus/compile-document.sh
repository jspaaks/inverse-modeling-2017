# call this as follows:
# ./compile-document.sh <name-of-tex-file>
# ./compile-document.sh ./tex/scge2011-inverse-modeling.tex
#

latex -halt-on-error -file-line-error -interaction=nonstopmode -output-format=pdf -output-directory=out/ $1

