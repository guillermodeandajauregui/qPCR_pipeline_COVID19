### Exploring the eds file format

### What we know so far:

- .eds files are zip files. Within them, there is an "apldbio" dir, containing an "sds" dir
- Within the sds dir there is a lot of files and dirs

#### the file "analysis_result.txt" seems to have most of the important information

##### This file starts with this header:

Session Name	
Well	Sample Name	Detector	Task	Ct	Avg Ct	Ct SD	Delta Ct	Qty	Avg Qty	Qty SD


Then registers begin, looking like this (elipses are mine, but they extend to 40; which matches what we know of my original files):

0	Agua	ORF1ab	Target	40.0
Rn values	5.2778935	5.280504	5.2740984	5.269397	....
Delta Rn values	-0.008802516	-0.0061365096	-0.012487017	.... 
DDCT Values	0		Agua	ORF1ab	Target	40.0									false
0	Agua	RNAsaP	Target	40.0						
Rn values	0.6095734	0.61071074	0.60998964	0.61091477	.....
Delta Rn values	-0.0049007013	-0.004908542	-0.006774796	-0.0069948295	....
DDCT Values	0		Agua	RNAsaP	Target	40.0									false

### Note that the wells are zero indexed. Based on the plate image shared from the machine, well 0 is in position "A1" and well 95 is position "H12"

