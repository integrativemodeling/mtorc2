# Analysis scripts for integrative structure determination of Aβ40-sTREM2 Complex

Directory for the analysis scripts for modeling of Aβ40fibril-sTREM2 complex.

## List of files and directories:

- 00_run_replica_exchange_analysis.py 			A script to check how good was the selection of temperature range of your replica exchange simulation   
- 01_run_analysis_trajectories.py 			A trajectory analysis script to extract relevant information like score satisfaction, individual score distribution 
- 02_run_extract_models_based_on_XLs_satisfraction.py 	A script to extract models which satisfied at least 80% of the XL-MS data	
- 03_structure_based_clustering.py  			A script to calculate lowest rmsd after trying 89 possible alignments between Aβ40s from two structures 
- 04_calculate_sampling_precision.py 			A script to check sampling exhaustiveness and sampling precision
- 05_final_transformation.py				A script to align the structures based on aligned results from 03_structure_based_clustering.py and put the cluster centers and all the members in seperate clusters
- 06_calculate_cluster_precision.py  			A script to calculate cluster precision
- density.txt						A defination file to calculate localized density

## Scripts used to analyze the comformational sampling, extracting structures of the Abeta40-sTREM2 complex 

- `00_run_replica_exchange_analysis.py`, `01_run_analysis_trajectories.py`, and `02_run_extract_models_based_on_XLs_satisfraction.py`  

## Scripts used to assess sampling exhaustiveness and computing sampling and structural precision 

- `03_structure_based_clustering.py`, `04_calculate_sampling_precision.py`, `05_final_transformation.py`, and `06_calculate_cluster_precision.py` .


## Information

_Author(s)_: Dibyendu Mondal

_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Publications_:
