# call this as follows:
# ./compile-document.sh
#

BASEFILENAME=inverse-modeling

rm ./out/${BASEFILENAME}.dvi
rm ./out/${BASEFILENAME}.idx
rm ./out/${BASEFILENAME}.toc
rm ./out/${BASEFILENAME}.out
rm ./out/${BASEFILENAME}.log
rm ./out/${BASEFILENAME}.aux
rm ./out/${BASEFILENAME}.bbl
rm ./out/${BASEFILENAME}.blg

latex -halt-on-error \
      -file-line-error \
      -interaction=nonstopmode \
      -output-format=pdf \
      -output-directory=out/ \
      ./tex/${BASEFILENAME}.tex

bibtex ./out/${BASEFILENAME}.aux

latex -halt-on-error \
      -file-line-error \
      -interaction=nonstopmode \
      -output-format=pdf \
      -output-directory=out/ \
      ./tex/${BASEFILENAME}.tex

latex -halt-on-error \
      -file-line-error \
      -interaction=nonstopmode \
      -output-format=pdf \
      -output-directory=out/ \
      ./tex/${BASEFILENAME}.tex
