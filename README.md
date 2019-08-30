
 
 # refGenes - Reference Genes Analysis
   This pipeline performs an automated analysis of the three most used algorithms to verify the stability of reference genes
 
## STEP.1 -Download this repository to a preference path:<br>
	# Cloning the whole repository (Git is required)
	# Open the terminal on unixOS and navigate to your work directory
   	 $ git clone https://github.com/hanielcedraz/endoGene.git
   	 $ cd rendoGene
	 $ chmod +x endoGene.R
	
	# Downloading zip file
	   - Under the repository name, click Clone or download.
	   - Click on Download ZIP
	  
## STEP.2 - Install R<br>
	Access https://cran.r-project.org


  ## Required files:
    dataset file (default: endogenous_ct.txt)
    efficiency list (default: efficiencies_list.txt)

  ## Usage on linux or Mac OS terminal: </br>
    $ ./refGenes.R -f dataset_file.txt -e efficiency_list.txt
    
    
    
    Obs. Also you can open the refGenes.R file with your favorite R IDE and follow the script (not recomended)
    


# results
 ## Tables
   ### Bestkeeper
   |   GeneName    |   SD_Value    |
   |   ---------   |   ---------   |
   |   RPL5    |   0.666    |
   |   MRPS27  |   0.819    |
   |   HPRT1   |   0.832    |
   |   HMBS    |   0.836    |
   |   MRPS30  |   0.92    |
   |   EEF1    |   0.926    |
   |   HSF3    |   1.068    |
   |   LDAH    |   1.101    |
   |   ACTB    |   1.167    |
   |   ACTA1    |   1.176    |
   |   HSP90    |   1.356    |
   |   TFRC    |   1.458    |
   |   HSP70    |   2.443    |


   ### Genorm
   |   rank  |  Gene |  MValue   |
   |   -----   |   -----   |   -----  |
   |   1  |   MRPS27  |   0.334591317    |
   |   1  |   MRPS30  |   0.334591317    |
   |   3  |   RPL5  |   0.480490089    |
   |   4  |   EEF1  |   0.642778744    |
   |   5  |   HPRT1  |   0.803047241    |
   |   6  |   LDAH  |   0.964291421    |
   |   7  |   ACTA1  |   1.039277933    |
   |   8  |   ACTB  |   1.118329727    |
   |   9  |   HMBS  |   1.19942579    |
   |   10  |   HSF3  |   1.269338025    |
   |   11  |   TFRC  |   1.369055031    |
   |   12  |   HSP90  |   1.471023544    |
   |   13  |   HSP70  |   1.657421221    |
   
   ### NormFinder
|	GeneName	 |  GroupDif	|	GroupSD |	Stability   |
|	-----	   |	-----		|	-----  |	-----	   |
|	RPL5	|	0.43	|	0.62	|	0.29	|
|	HPRT1	|	0.34	|	0.76	|	0.31	|
|	MRPS30	|	0.34	|	0.75	|	0.31	|
|	MRPS27	|	0.4	|	0.75	|	0.32	|
|	ACTA1	|	0.06	|	1.27	|	0.32	|
|	EEF1	|	0.56	|	0.78	|	0.34	|
|	LDAH	|	0.27	|	1.11	|	0.37	|
|	HSF3	|	0.4	|	1.13	|	0.4	|
|	ACTB	|	1.47	|	0.91	|	0.49	|
|	HMBS	|	1.44	|	1.1	|	0.5	|
|	TFRC	|	0.95	|	1.65	|	0.55	|
|	HSP90	|	1.76	|	1.44	|	0.56	|
|	HSP70	|	1.11	|	2.51	|  0.78	|



   ### Final Ranking
|	Rank	|	GeneName	|
|	---------	|	---------	|
|	1	|	MRPS27	|
|	2	|	MRPS30	|
|	3	|	RPL5	|
|	4	|	HPRT1	|
|	5	|	EEF1	|
|	6	|	LDAH	|
|	7	|	ACTA1	|
|	8	|	HSF3	|
|	9	|	HMBS	|
|	10	|	ACTB	|
|	11	|	TFRC	|
|	12	|	HSP90	|
|	13	|	HSP70	|




## Plots
<a href="https://ibb.co/Fnk2QBB"><img src="https://i.ibb.co/Fnk2QBB/Rplot-gene-stability-by-Best-Keeper.png" alt="Rplot-gene-stability-by-Best-Keeper" border="0"></a> 
<a href="https://ibb.co/J3jxYV0"><img src="https://i.ibb.co/J3jxYV0/Rplot-gene-stability-by-genorm.png" alt="Rplot-gene-stability-by-genorm" border="0"></a> <br>
<a href="https://ibb.co/rtjgk3w"><img src="https://i.ibb.co/rtjgk3w/Rplot-gene-stability-by-Norm-Finder.png" alt="Rplot-gene-stability-by-Norm-Finder" border="0"></a> 
<a href="https://ibb.co/m4PsPjD"><img src="https://i.ibb.co/m4PsPjD/Rplot-gene-variation-by-genorm.png" alt="Rplot-gene-variation-by-genorm" border="0"></a>
