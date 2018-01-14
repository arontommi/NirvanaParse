#!/bin/bash

#input = $1
#output = $2
#dir for nirvana = $3



dotnet	~/$3/nirvana/Nirvana.dll \
	-c ~/$3/nirvana/Data/Cache/24/GRCh37/Ensembl84 \
	--sd ~/$3/nirvana/Data/GRCh37 \
	-r ~/$3/nirvana/Data/References/5/Homo_sapiens.GRCh37.Nirvana.dat \
	--vcf \
    -i $1 \
    -o $2
