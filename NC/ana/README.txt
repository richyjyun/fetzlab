Simple functions necessary for extracting sleep times (using delta oscillations and accelerometers), performing spike sorting (using PCA and 2-window discrimination), and plotting stLFPs over time. 

RJY 06/22/2018

Hierarchy:

ParseNC (example script)
getMetaData
saveSleep
	getOscillPower
	getSleepTimes
saveSortingParams
	getSortingParams
saveSortedSpikes (requires saveSortingParams)
	sortChannel
saveSTLFP (requires saveSortedSpikes)
	callps2pdf
		ps2pdf
	