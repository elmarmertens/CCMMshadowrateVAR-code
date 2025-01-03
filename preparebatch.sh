#!/bin/bash
# shell script to launch matlab codes in background

wd=\'$(pwd)\'
preamble='cd '${wd}'; display(pwd), parpool; '
# echo $preamble

os=$(uname)

if [ $os = Darwin ]; then
	alias matlab='~/Applications/MATLAB_R2023b.app/bin/matlab'
	caffeinate -iw $$ &
fi

ELBbound=0.25
ELBlabel=25
p=12 # common  lag length for each VAR

MCMCdraws=1e3

foobatch='runbatch.m'
rm -f $foobatch

for mfile in $@; do
	thismfile=$(basename $mfile .m)

	for thislabel in fredsxMD20 fredsxMD20exYield fredsxMD14longyields; do

		# apply fredsxMD14longyields only to linear VAR
		if [ "$thislabel" = "fredsxMD14longyields" ] && [ "$thismfile" != "goVAR" ]; then
			continue
		fi

		# skip general model with yields
		if [ "$thismfile" = "goVARshadowrateGeneral" ] && [ "$thislabel" = "fredsxMD20" ]; then
			continue
		fi

		datalabel=${thislabel}-2022-09
		echo gobatcher:$thismfile w/$MCMCdraws $datalabel

		foofile=bot${thismfile}${thislabel}
		echo $foofile >>$foobatch

		# overwrite parameters
		sed '/SED-PARAMETERS-HERE/a\ 
        datalabel='\'$datalabel\''; \
        p='$p'; \
        ELBbound='$ELBbound'; \
		MCMCdraws='$MCMCdraws'; \
		fcstNdraws= 10 * MCMCdraws; \
        ' \
			$mfile >$foofile.m

	done # datalabel
done  # mfiles