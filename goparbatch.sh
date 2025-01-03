#!/bin/zsh
# shell script to launch matlab codes in background

wd=$(pwd)
preamble="cd ${wd}; display $wd; parpool; "
# echo $preamble

os=$(uname)

if [[ $os = Darwin ]]
then
     alias matlab='~/Applications/MATLAB_R2024b.app/bin/matlab'
     caffeinate -iw $$ &
fi

for mfile in "$@" 
do
     thismfile=$(basename "$mfile" .m)
     # echo $thismfile
     matlab -nodesktop -noFigureWindows -batch "$preamble $thismfile"
done

