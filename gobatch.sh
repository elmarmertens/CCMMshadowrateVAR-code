#!/bin/bash
# shell script to launch matlab codes in background

wd=\'`pwd`\'
preamble='cd '${wd}'; display(pwd), parpool; '
# echo $preamble


os=$(uname)

if [ $os = Darwin ]
then
     alias matlab='~/Applications/MATLAB_R2022a.app/bin/matlab'
     caffeinate -iw $$ &
fi

for mfile in $@ 
do
     thismfile=`basename $mfile .m`
     # echo $thismfile
     matlab -nodisplay -batch "$preamble $thismfile"
done

