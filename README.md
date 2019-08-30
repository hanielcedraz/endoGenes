
 
 # refGenes - Reference Genes Analysis
   This pipeline performs an automated analysis of the three most used algorithms to verify the stability of reference genes
 
## STEP.1 -Download this repository to a preference path:<br>
	# Cloning the whole repository (Git is required)
	# Open the terminal on unixOS and navigate to your work directory
   	 $ git clone https://github.com/hanielcedraz/refGenes.git
   	 $ cd refGenes
	
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
    
    
|   GeneName    |   SD_Value    |
|   ---------   |   ---------   |
|   RPL5    |   0.666    |
|   MRPS27  |   0.819    |
|   HPRT1   |   0.832    |
|   HMBS    |   0.836    |
|   MRPS30  |   0.92     |
|   EEF1    |   0.926    |


<a href="https://ibb.co/Fnk2QBB"><img src="https://i.ibb.co/Fnk2QBB/Rplot-gene-stability-by-Best-Keeper.png" alt="Rplot-gene-stability-by-Best-Keeper" border="0"></a> <a href="https://ibb.co/J3jxYV0"><img src="https://i.ibb.co/J3jxYV0/Rplot-gene-stability-by-genorm.png" alt="Rplot-gene-stability-by-genorm" border="0"></a> <a href="https://ibb.co/rtjgk3w"><img src="https://i.ibb.co/rtjgk3w/Rplot-gene-stability-by-Norm-Finder.png" alt="Rplot-gene-stability-by-Norm-Finder" border="0"></a> <a href="https://ibb.co/m4PsPjD"><img src="https://i.ibb.co/m4PsPjD/Rplot-gene-variation-by-genorm.png" alt="Rplot-gene-variation-by-genorm" border="0"></a>
