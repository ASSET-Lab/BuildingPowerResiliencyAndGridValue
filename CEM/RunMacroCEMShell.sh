#!/bin/bash

for wsGen in "70"; do
	for wsGenYr in "2030"; do
		for co2Cap in "90"; do
			for co2CapYr in "2045"; do
				for forceBldgUp in "0" "1" "5"; do
					sbatch RunMacroCEMJob.sbat $wsGen $wsGenYr $co2Cap $co2CapYr $forceBldgUp
				done
			done
		done
	done
done
