
 
 # endoGenes- Reference Genes Analysis
   This pipeline performs an automated analysis of the three most used algorithms to verify the stability of reference genes
 
## STEP.1 -Download this repository to a preference path:<br>
	# Cloning the whole repository (Git is required)
	# Open the terminal on unixOS and navigate to your work directory
   	 $ git clone https://github.com/hanielcedraz/endoGenes.git
   	 $ cd endoGenes
	 $ chmod +x endoGenes.R
	
	# Downloading zip file
	   - Under the repository name, click in Clone or download.
	   - Click on Download ZIP
	  
## STEP.2 - Install R<br>
	Access https://cran.r-project.org


  ## Required files:
    dataset file (default: endogenous_ct.txt)
    efficiency list (default: efficiencies_list.txt)

  ## Usage on linux or Mac OS terminal: </br>
    $ ./endoGenes.R -f dataset_file.txt -e efficiency_list.txt
    
    
    
    ps. If you do not use Linux or macOS, use the file run_endoGenes_notebook.ipynb to run your analysis.
    


# Results
 ## Tables
   ### Bestkeeper
 <a href="https://imgbb.com/"><img src="https://i.ibb.co/SBcygpd/Screenshot-2019-11-21-at-15-41-06.png" alt="Screenshot-2019-11-21-at-15-41-06" border="0"></a>

   ### Genorm
<a href="https://imgbb.com/"><img src="https://i.ibb.co/GFc62DB/Screenshot-2019-11-21-at-15-33-01.png" alt="Screenshot-2019-11-21-at-15-33-01" border="0"></a>
   
   ### NormFinder
<a href="https://imgbb.com/"><img src="https://i.ibb.co/VxD8SJF/Screenshot-2019-11-21-at-15-36-55.png" alt="Screenshot-2019-11-21-at-15-36-55" border="0"></a>

   ### Final Ranking
<a href="https://imgbb.com/"><img src="https://i.ibb.co/W0HQjJy/Screenshot-2019-11-21-at-15-44-52.png" alt="Screenshot-2019-11-21-at-15-44-52" border="0"></a>




## Plots
<a href="https://ibb.co/Fnk2QBB"><img src="https://i.ibb.co/Fnk2QBB/Rplot-gene-stability-by-Best-Keeper.png" alt="Rplot-gene-stability-by-Best-Keeper" border="0"></a> 
<a href="https://ibb.co/J3jxYV0"><img src="https://i.ibb.co/J3jxYV0/Rplot-gene-stability-by-genorm.png" alt="Rplot-gene-stability-by-genorm" border="0"></a> <br>
<a href="https://ibb.co/rtjgk3w"><img src="https://i.ibb.co/rtjgk3w/Rplot-gene-stability-by-Norm-Finder.png" alt="Rplot-gene-stability-by-Norm-Finder" border="0"></a> 
<a href="https://ibb.co/m4PsPjD"><img src="https://i.ibb.co/m4PsPjD/Rplot-gene-variation-by-genorm.png" alt="Rplot-gene-variation-by-genorm" border="0"></a>
