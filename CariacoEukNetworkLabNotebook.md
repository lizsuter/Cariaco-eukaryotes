##**Cariaco Eukaryotic Paper**  
###**Lab Notebook for Network Analysis** 
*Begin Jun 2018* 

Supplements Matlab and R scripts. Workflow as follows:  
1)  "ExtractTopOTURawAbundanceTable.m"  *for extracting raw abundances from raw data*  
2) "CariacoNetwork_Rscript.R" *for setting up matrices*  
3) Python script (embedded here or in R script) for calculating network using [SparCC](https://bitbucket.org/yonatanf/sparcc/overview)  
4) Back to R script for making plot file  
5) Use exported plot file in [Gephi](https://gephi.org/) to visualize and make network figure  
  
This workflow follows that of Gerikas et al. 2018 (ISME), whose scripts are available at their [Github site](https://github.com/vaulot/Ribeiro_CARBOM_ISME_2018)  

## June 12, 2018
*I am starting this lab notebook today however many notes were already taken in the R and matlab scripts. I will do my best to combine all notes here and leave the scripts as "clean" as possible so they can be posted on Github.*  
Up to now, I was formatting the taxa tables to match the bubble plots by summing at interesting taxonomic levels (eg. Alveolata; Protalveolata; Syndiniales; Group I *or* uncultured marine labyrithulid DH417-EKD10. See bubble plots for all). However, another scan of the literature shows that other groups doing network analyses only use OTU-level taxonomic classification. Some papers use the top 100 most abundant OTUs or, for example, Gerikas Ribeiro et al. use the OTUs with >20% abundance. For us, a 20% cutoff this would drastically limit the # of OTUs used in the analysis. However, I think it's important to use a consistent cutoff (maybe 1%?) rather than just "top 100" since we are going to pull OTUs from the 3 different datasets.  
 
I am going back to Matlab script in **Documents\Manuscripts\Cariaco_Euks\Network\Cariaco>v1\ExtractTopOTURawAbundanceTable.m** and modifying the script so that it pulls the top OTUs from all 3 datasets (Euk 18S, Bac 16S, Arch 16S) for which the abundance of that OTU is >1% in at least one sample. The script will identify those top OTUs, then pull the *RAW* abundance of that OTU from each sample, without averaging across duplicate samples. This script and subsequent downstream R scripts, etc are in the folder **v2**.  
  
  
  
*Some interesting things that I noticed while doing this:*   

* One of the abundat Euk OTUs (denovo270351) is specifically a Cariaco Basin clade (an Armophorea ciliate).  
* An OTU (denovo195279) that comes out as "uncultured" even at the phylum level, but seems important in the network analysis, is closely related to an unculultured Cariaco protist from Edgcomb et al. 2011 (Accession GU820578.1). They call it an uncultired dinoflagellate (clone AA5F15RM1H02)
	* There are 3 other abundant OTUs in the Euk domain that were not assigned a phylum but blasting them shows significant matches. I left all 3 as "undefined" in the network but keep in mind:
		* denovo285757 matches rhizarian clones. Some are from Cariaco basin (see Genbank GU821081.1, which is clone BC95F13RM1A09 from Orsi et al. 2011)
		* denovo12021 matches rhizarian clones. Top hit is from Orsi et al. (Genbank GU821495.1, clone AA3F14RJ2H11)
		* denovo264290 matches "uncultured eukaryote clone"- maybe alveolate/ syndiniales? Top hit is Lie et al. 2014
		* **NOTE** why aren't the hits to our Cariaco study coming up? Is it because they are assembled so they are not in this database?
	
* The Euk OTU with no blast hit at the domain level (denovo289734) has a closest blastn relative as uncultured Syndiniales. Not sure why this didn't come out in the clustering based on SILVA but I excluded from analysis due to this.
	* Same thing happened for a Bac OTU which was in top OTUsd (rank 81). Denovo158368 showed no domain in OTU table but has a blast hit similar to uncultured alpha proteobacteria from other studies. Deleted for this analysis though.
* Making a cutoff of 0.2 for the correlation matrix for plotting the network makes a really complex network since we have a high diversity. So I increased the cutoff so that there would be less edges- only the really strong ones. However, this leads to a matrix with islands that are separate from the core of the network. I think I need to delete the nodes as well 
	* Gerikas Ribeiro et al. specifically comment in their script that when they made the cutoff 0.2, they didn't need to delete any nodes. But we may need to do that.


**Python script for building network of top Euk OTUs** after formatting matrix in R.
NOTE- it would be better to do section section in a loop. Figure that out. (*Running like this leads to some kind of error. Not every line works and doesn't create the "perm_cor) file. So I keep having to go back, see which files weren't written and rerun that line. I have no idea why this happens.
Also double check this logically against the SparrCC site and Gerikas et al.

```
python2.7 sparcc_euksonly/SparCC.py 18Sotu_matrix.txt -i 5 --cor_file=sparcc_euksonly/basis_corr/cor_sparcc.out
python2.7 sparcc_euksonly/SparCC.py 18Sotu_matrix.txt -i 5 --cor_file=sparcc_euksonly/basis_corr/cor_pearson.out -a pearson
python2.7 sparcc_euksonly/SparCC.py 18Sotu_matrix.txt -i 5 --cor_file=sparcc_euksonly/basis_corr/cor_spearman.out -a spearman

python2.7 sparcc_euksonly/MakeBootstraps.py 18Sotu_matrix.txt -n 100 -t permutation_#.txt -p sparcc_euksonly/pvals/

python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_0.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_0.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_1.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_1.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_2.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_2.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_3.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_3.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_4.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_4.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_5.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_5.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_6.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_6.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_7.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_7.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_8.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_8.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_9.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_9.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_10.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_10.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_11.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_11.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_12.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_12.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_13.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_13.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_14.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_14.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_15.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_15.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_16.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_16.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_17.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_17.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_18.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_18.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_19.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_19.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_20.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_20.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_21.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_21.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_22.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_22.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_23.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_23.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_24.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_24.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_25.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_25.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_26.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_26.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_27.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_27.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_28.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_28.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_29.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_29.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_30.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_30.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_31.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_31.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_32.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_32.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_33.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_33.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_34.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_34.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_35.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_35.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_36.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_36.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_37.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_37.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_38.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_38.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_39.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_39.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_40.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_40.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_41.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_41.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_42.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_42.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_43.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_43.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_44.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_44.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_45.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_45.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_46.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_46.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_47.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_47.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_48.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_48.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_49.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_49.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_50.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_50.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_51.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_51.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_52.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_52.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_53.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_53.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_54.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_54.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_55.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_55.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_56.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_56.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_57.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_57.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_58.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_58.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_59.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_59.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_60.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_60.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_61.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_61.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_62.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_62.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_63.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_63.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_64.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_64.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_65.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_65.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_66.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_66.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_67.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_67.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_68.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_68.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_69.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_69.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_70.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_70.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_71.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_71.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_72.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_72.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_73.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_73.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_74.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_74.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_75.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_75.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_76.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_76.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_77.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_77.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_78.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_78.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_79.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_79.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_80.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_80.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_81.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_81.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_82.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_82.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_83.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_83.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_84.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_84.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_85.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_85.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_86.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_86.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_87.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_87.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_88.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_88.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_89.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_89.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_90.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_90.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_91.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_91.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_92.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_92.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_93.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_93.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_94.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_94.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_95.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_95.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_96.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_96.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_97.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_97.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_98.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_98.txt
python2.7 sparcc_euksonly/SparCC.py sparcc_euksonly/pvals/permutation_99.txt -i 5 --cor_file=sparcc_euksonly/pvals/perm_cor_99.txt


python2.7 sparcc_euksonly/PseudoPvals.py sparcc_euksonly/basis_corr/cor_sparcc.out sparcc_euksonly/pvals/perm_cor_#.txt 100 -o sparcc_euksonly/pvals/pvals.one_sided.txt -t one_sided
python2.7 sparcc_euksonly/PseudoPvals.py sparcc_euksonly/basis_corr/cor_sparcc.out sparcc_euksonly/pvals/perm_cor_#.txt 100 -o sparcc_euksonly/pvals/pvals.two_sided.txt -t two_sided
```



**Python script for building network of top OTUs from all 3 domains** after formatting matrix in R.
NOTE- it would be better to do second section in a loop. Figure that out.   
Also double check this logically against the SparrCC site and Gerikas et al.

```
python2.7 sparcc_alldomains/SparCC.py all_domain_otu_matrix.txt -i 5 --cor_file=sparcc_alldomains/basis_corr/cor_sparcc.out
python2.7 sparcc_alldomains/SparCC.py all_domain_otu_matrix.txt -i 5 --cor_file=sparcc_alldomains/basis_corr/cor_pearson.out -a pearson
python2.7 sparcc_alldomains/SparCC.py all_domain_otu_matrix.txt -i 5 --cor_file=sparcc_alldomains/basis_corr/cor_spearman.out -a spearman

python2.7 sparcc_alldomains/MakeBootstraps.py all_domain_otu_matrix.txt -n 100 -t permutation_#.txt -p sparcc_alldomains/pvals/

python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_0.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_0.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_1.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_1.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_2.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_2.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_3.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_3.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_4.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_4.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_5.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_5.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_6.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_6.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_7.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_7.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_8.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_8.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_9.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_9.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_10.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_10.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_11.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_11.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_12.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_12.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_13.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_13.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_14.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_14.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_15.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_15.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_16.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_16.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_17.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_17.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_18.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_18.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_19.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_19.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_20.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_20.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_21.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_21.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_22.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_22.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_23.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_23.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_24.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_24.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_25.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_25.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_26.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_26.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_27.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_27.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_28.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_28.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_29.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_29.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_30.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_30.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_31.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_31.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_32.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_32.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_33.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_33.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_34.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_34.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_35.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_35.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_36.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_36.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_37.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_37.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_38.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_38.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_39.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_39.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_40.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_40.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_41.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_41.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_42.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_42.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_43.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_43.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_44.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_44.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_45.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_45.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_46.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_46.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_47.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_47.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_48.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_48.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_49.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_49.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_50.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_50.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_51.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_51.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_52.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_52.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_53.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_53.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_54.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_54.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_55.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_55.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_56.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_56.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_57.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_57.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_58.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_58.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_59.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_59.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_60.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_60.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_61.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_61.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_62.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_62.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_63.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_63.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_64.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_64.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_65.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_65.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_66.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_66.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_67.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_67.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_68.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_68.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_69.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_69.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_70.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_70.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_71.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_71.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_72.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_72.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_73.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_73.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_74.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_74.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_75.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_75.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_76.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_76.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_77.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_77.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_78.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_78.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_79.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_79.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_80.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_80.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_81.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_81.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_82.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_82.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_83.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_83.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_84.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_84.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_85.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_85.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_86.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_86.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_87.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_87.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_88.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_88.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_89.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_89.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_90.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_90.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_91.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_91.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_92.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_92.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_93.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_93.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_94.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_94.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_95.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_95.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_96.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_96.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_97.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_97.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_98.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_98.txt
python2.7 sparcc_alldomains/SparCC.py sparcc_alldomains/pvals/permutation_99.txt -i 5 --cor_file=sparcc_alldomains/pvals/perm_cor_99.txt


python2.7 sparcc_alldomains/PseudoPvals.py sparcc_alldomains/basis_corr/cor_sparcc.out sparcc_alldomains/pvals/perm_cor_#.txt 100 -o sparcc_alldomains/pvals/pvals.one_sided.txt -t one_sided
python2.7 sparcc_alldomains/PseudoPvals.py sparcc_alldomains/basis_corr/cor_sparcc.out sparcc_alldomains/pvals/perm_cor_#.txt 100 -o sparcc_alldomains/pvals/pvals.two_sided.txt -t two_sided
```



## June 14-19, 2018
*Next goals*  

* Create a metadata file to import into Gephi so that we can color code network plot by things like domain, etc.
	* Make excel file (sparcc_metadata) where I split the full OTU name; taxonomy. Made column for domain, phylum or class; and short name (where I used the most common name assoicated with that OTU and the furthest taxonomic refinement that was identified (eg. clone name, etc)
	* Remembered to take Nitrospinae out of Deltaproteibacteria, like in itags paper. Actual taxonomy is Bacteria; Nitrospinae; Nitrospinae
	* Abbreviations used in short name column
		* MG = Marine Group
		* THSG = terrestrial hot spring group
		* MCG = Miscellaneous crenarchaeaotic group
		* MBG = Marine Benthic Group
		* DSEG = deep sea Eutyarchaeotic group
		* DHVEG = deeps sea hydrothermal vent Euryarchaeota group
		* BS-GSO2 = Baltic Seaa Gammaproteobacterial Sulfur Oxidizer 2 [Thiorhodospira]
		* MCG-15 = Miscellaneous crenarchaeaotic group 15 (*changed from Group C3*)
		* GSO477 (*changed from E01-9C-26 marine group*)
		* Arctic96BD-19 (*changed from ZD0405*)
	* Include "Trophic Level" as a variable in metadata
		* I am classifying certain groups as autotrophic, heterotrophic, etc. These are **putatitive**. I can't identify all but here are which clades I am identifying
		* <u>Autotrophic</u>
			* Archaea
				* MG-I
			* Bacteria
				* BS-GSO2
				* SAR406  (MG A)
				* GSO477
				* Arctic96BD-19
				* SUP05
				* Nitrospinae
				* MG-B
				* Arctic97B-4 marine group **putative** [see Yilmaz et al. 2015]
				* SAR202  **putative** [see Yilmaz et al. 2015]
		* <u>Heterotrophic/ Fermentative</u> *Consider splitting these into aerobic heterotroph/ sulfate reducer/ fermentative*
			* Archaea
				* MG-II (according to Zhang et al. 2015)
			* Bacteria
				* SAR11
				* Sphingobacteriales;  WCHB1-69
				* Burkholderia glumae
				* Salinisphaeraceae;  ZD0417
				* Flavobacteriaceae;  Sufflavibacter (according to Kwon et al. 2007)
				* OM190 (according to Lage and Bondoso 2014)
				* Fibrobacteria;  P.palm C70
				* Bacillus sp. Asd3
				* MSBL8
				* Cytophagales;  Flammeovirgaceae;  Marinoscillum
				* Flavobacteriaceae; NS2b marine group
				* Pseudoalteromonas;  bacterium 2D702 (blasted this OTU, denovo219642, and it's closet relative is Pseudoalteromonas shioyasakiensis. See [Matsuyama et al. 2014](http://ijs.microbiologyresearch.org/content/journal/ijsem/10.1099/ijs.0.055558-0#tab2) for some characteristics.
				* Desulfobacteraceae;  SEEP-SRB1
				* SAR86 
				* Desulfarculaceae
				* Phycisphaeraceae;  JL-ETNP-F27
				* Hyd24-12 [**putative fermentative**]
				* Vibrio
				* Pla3 lineage [*this may be wrong. According to Elshahed et al. 2007, the branch previously known as the Pla3 lineage includes the anammox clades. So this could be autotrophic. However, they are mostly PA which I think suggests they are heterotrophic. Also the closest blast relative to denovo317108 is in the Chloroflexi*]
				* WS3 [**putative fermentative** see Baker et al. 2015]
				* Flavobacteriales; Cryomorphaceae;  Fluviicola [see Ch. 42, Bowman in The Prokaryotes]
				* Desulfobacteraceae;  SEEP-SRB1
				* Acidimicrobiales;  Sva0996 marine group [see Chen et al. 2016 and Connelly et al. 2014]
		* <u>Bacterial Predator</u> - **putative**
			* Sphingobacteriales;  Saprospiraceae. [See McIlroy and Nielsen 2014 in The Prokaryotes]
		* <u>Syntrophic</u> - **putative**
			* Coprothermobacters p. Dex80-4 may be fermentative with syntrophic relationship with archaeal hydrogenotrophs. See Palatsi et al. 2011. and Gagliano et al. 2015
			* Anaerolineaceae. See Liang et al. 2015 "*Bacterial clone libraries showed organisms of Anaerolineaceae (within the phylum of Chloroflexi) were predominant (45.5%), indicating syntrophically cooperation with Methanosaeta archaea was likely involved in the process of methanogenic degradation of alkanes*."
		* <u>Did not give designation to the following groups because of unclear ecological roles</u> *Come back and try to pull out methanogens*
			* Archaea
				* VC2.1Arc6
				* MCG-15
				* MCG
				* MG-III
				* CCA47
				* DSEG
				* MEG
				* MHVG
				* Thaumarchaeota terrestrial group
				* MBG-B
				* MBG-D and DHVEG-1
				* THSCG
		* <u>Look into possible interesting roles of:</u>
			* Burkholderia glumae is land plant pathogen (phytopathogen). It did not come out significantly correlated in network so didn't look into it much further.
			* Saprospiraceae family: "Isolates and in situ strains have a demonstrated ability for the hydrolysis and utilization of complex carbon sources, with the helical gliding strains also associated with predation of other bacteria and algae." (McIlroy and Nielsen 2014)
			*  SAR11 and SAR86. Heterotrophic but also genome streamlined. They exhibit carbon compound specialization.
				* See Dupont et al. 2012 





## Jan 17, 2019
*I have continued this analysis but using PR2 annotations instead of SILVA, which Ginny and Taylor suggested would be more rubust for Eukaryotes. Maria redid the annotation earlier in the year*.  

*Also in addition to the [SparCC](https://bitbucket.org/yonatanf/sparcc/overview)   website I mentioned earlier, I found this [tutorial](https://rachaellappan.github.io/16S-analysis/correlation-between-otus-with-sparcc.html#running-sparcc) helpful.*


- I have been plotting abundances and manipulating text files in R
	- See ```/Users/admin/Google Drive/Wagner/Research/Manuscripts/Cariaco_Euks/Cariaco_euks.Rproj```

	- Made text file with most abundant eukaroyotes (>0.1% in any sample) = ```Euks.001.txt```
- Run that in SparCC using python (from working directory ```/Users/admin/Google Drive/Wagner/Research/Manuscripts/Cariaco_Euks```)

```
mkdir Network_PR2_OTUs/euks.001_basis_corr

python2.7 sparcc-master/SparCC.py Network_PR2_OTUs/euks.001.txt --cor_file=Network_PR2_OTUs/euks.001_basis_corr/cor_sparcc.out
```

Output:   

```
reading data
Finished parsing table.
Table dimensions: (41,837)
**** Data has been transposed! ****
First 3 column labels are :denovo158326 ,denovo12177 ,denovo155840
First 3 row labels are :AE3b103B ,AE3b103A ,AE3b234A 

computing correlations
	Running iteration0
	Running iteration1
	Running iteration2
	Running iteration3
	Running iteration4
	Running iteration5
	Running iteration6
	Running iteration7
	Running iteration8
	Running iteration9
	Running iteration10
	Running iteration11
	Running iteration12
	Running iteration13
	Running iteration14
	Running iteration15
	Running iteration16
	Running iteration17
	Running iteration18
	Running iteration19
(837, 837)
writing results
wrote Network_PR2_OTUs/euks.001_basis_corr/cor_sparcc.out
wrote cov_mat_SparCC.out
Done!
```

- The above produces the SparCC correlations, which has certain benefits as described by Rachael Lappan in the [website](https://rachaellappan.github.io/16S-analysis/correlation-between-otus-with-sparcc.html#running-sparcc) cited above:
	- "*The issue with correlations between bacteria in microbiome data is that the data is compositional; it describes the relative abundance of bactera, not absolute like if you were to determine bacterial load with qPCR. Each OTU contributes a proportion of reads that add up to 100% in each sample. The problem with searching for correlations here is that if there is an increase in one OTU, all of the other OTUs must decrease, and this ends up looking like a negative correlation between the increased OTU and everything else when this may not be the case. A tool called SparCC determines correlations between OTUs taking this specific problem into account.*"
	- But I could also calculate a Spearman or Pearson Correlation matrix (as I did last time)
- Next step is to make the bootstrapped datasets (create 100 simulated datasets of the OTU table, shuffled randomly) in order to generate pseduo p-values
	- The [tutorial](https://rachaellappan.github.io/16S-analysis/correlation-between-otus-with-sparcc.html#running-sparcc) helps by using a loop, which I was trying to do last time I did this but couldn't figure it out:

```
mkdir Network_PR2_OTUs/euks.001_simulated_datasets/

python2.7 sparcc-master/MakeBootstraps.py Network_PR2_OTUs/euks.001.txt -n 100 -t permutation_#.txt -p Network_PR2_OTUs/euks.001_simulated_datasets/
```


- Next run SparCC on each of these simulated files

```
mkdir Network_PR2_OTUs/euks.001_boot_corr
mkdir Network_PR2_OTUs/euks.001_boot_cov

for i in `seq 0 99`; do python2.7 sparcc-master/SparCC.py Network_PR2_OTUs/euks.001_simulated_datasets/permutation_$i.txt -c Network_PR2_OTUs/euks.001_boot_corr/simulated_sparcc_$i.txt -v Network_PR2_OTUs/euks.001_boot_cov/simulated_sparcc_$i.txt >> Network_PR2_OTUs/boot_euks.001_corr_sparcc.log; done


```

- This takes a long time (It's doing the 20 SparCC iterations for each of the 100 simulated datasets). Left this running at 2:24pm and went to get lunch.
- Results of each run are in log file.

- Next calculated the pseudo pvalues
	- In the tutorial they use one-sided p values because they "Take into account the direction of the correlation." Look into this

```
mkdir Network_PR2_OTUs/euks.001_pvals

python2.7 sparcc-master/PseudoPvals.py Network_PR2_OTUs/euks.001_basis_corr/cor_sparcc.out Network_PR2_OTUs/euks.001_boot_corr/simulated_sparcc_#.txt 100 -o Network_PR2_OTUs/euks.001_pvals/pvals_one_sided.txt -t one_sided
```

- Now move back into R for plotting




## Jan 21, 2019
*Rerunning the analysis but this time for eukaryote OTUs that are at least 0.01% abundant (RA >= 0.0001)*.  
*NOTE: This will be much slower than last time since so many more OTUs. Now I retained 5466 OTUs rather than 837 like before when cutoff was 0.1%*.  



- See ```/Users/admin/Google Drive/Wagner/Research/Manuscripts/Cariaco_Euks/Cariaco_euks.Rproj```

	- And text file ```euks.0001.txt```
- Run in SparCC using python (from working directory: ```/Users/admin/Google Drive/Wagner/Research/Manuscripts/Cariaco_Euks```)

```
mkdir Network_PR2_OTUs/euks.0001_basis_corr

python2.7 sparcc-master/SparCC.py Network_PR2_OTUs/euks.0001.txt --cor_file=Network_PR2_OTUs/euks.0001_basis_corr/cor_sparcc.out
```

_ This is working through the 20 interations very slowly (~3 hours!)

Log:   

```
reading data
Finished parsing table.
Table dimensions: (41,5466)
**** Data has been transposed! ****
First 3 column labels are :denovo148708 ,denovo15209 ,denovo12896
First 3 row labels are :AE3b103B ,AE3b103A ,AE3b234A 

computing correlations
	Running iteration0
in slow
	Running iteration1
in slow
	Running iteration2
in slow
	Running iteration3
in slow
	Running iteration4
in slow
	Running iteration5
in slow
	Running iteration6
in slow
	Running iteration7
in slow
	Running iteration8
in slow
	Running iteration9
in slow
	Running iteration10
in slow
	Running iteration11
in slow
	Running iteration12
in slow
	Running iteration13
in slow
	Running iteration14
in slow
	Running iteration15
in slow
	Running iteration16
in slow
	Running iteration17
in slow
	Running iteration18
in slow
	Running iteration19
in slow
(5466, 5466)
writing results
wrote Network_PR2_OTUs/euks.0001_basis_corr/cor_sparcc.out
wrote cov_mat_SparCC.out
Done!
```

- To run the bootstrapping this is going to take forever on my computer so I will launch an AWS instance and run it from there
-  Re-started my EDAMAME instance in amazon dashboard (this one already has EBS volume attached. Been paying "rent" on it for several months now. Logged in in terminal:

```
chmod og-rwx ~/Documents/amazon-key.pem
ssh -i ~/Documents/amazon-key.pem ubuntu@ec2-54-157-4-252.compute-1.amazonaws.com

```

- Mount the EBS volume that I made last time when I was working on Cariaco N cycling manuscript. That holds all my input/output files:

```
sudo mount /dev/xvdf Cariaco
```

- Open file browser for this volume in Filezilla
	- Host is the public DNS from Amazon EC2 (ec2-3-86-89-79.compute-1.amazonaws.com), the user name is ubuntu, port is 22, and leave password empty. Load key file as password (Settings/ SFTP)
- Made new directory on EBS: ```~/Cariaco/Network_Analysis/Euk_Network```
	- Copied here the "sparcc-master" folder with the sparcc functions
	- Also copy the results of the Sparcc on my dataset (when done running) (```euks.0001_basis_corr/cor_sparcc.out``` file)
	- Also copy the initial ```euks.0001.txt``` input file
- Instance has python ```python -V``` returns ```Python 2.7.6```
	- But had to update, install numpy & pandas for the scripts to work
```
sudo apt-get update
sudo apt-get install python-numpy
sudo apt-get install python-pandas
```

- NOTE I tried runnning the initial SparCC.py on my dataset in AWS but it is taking just as long and the one on my local computer is more than halfway done. So just import results of the local one into EBS.
- Next make the bootstrapped datasets (create 100 simulated datasets of the OTU table, shuffled randomly) in order to generate pseduo p-values
	-	Run in tmux so don't lose it if connection is lost 
	
```
tmux new -s EukNetwork

mkdir euks.0001_simulated_datasets/

python2.7 sparcc-master/MakeBootstraps.py euks.0001.txt -n 100 -t permutation_#.txt -p euks.0001_simulated_datasets/
```

- To get out of tmux session:
	- ctl+b, let go, d
- To get back into tmux session:
	- ```tmux attach -t EukNetwork```

Log:  
```
Finished parsing table.
Table dimensions: (41,5466)
**** Data has been transposed! ****
First 3 column labels are :denovo148708 ,denovo15209 ,denovo12896
First 3 row labels are :AE3b103B ,AE3b103A ,AE3b234A 
```



- Next run SparCC on each of these simulated files. This will take a long time so do 5 iterations on each of the 100 bootstrapped files (instead of 20)

```
mkdir euks.0001_boot_corr
mkdir euks.0001_boot_cov

for i in `seq 0 99`; do python2.7 sparcc-master/SparCC.py euks.0001_simulated_datasets/permutation_$i.txt -i 5 -c euks.0001_boot_corr/simulated_sparcc_$i.txt -v euks.0001_boot_cov/simulated_sparcc_$i.txt >> boot_euks.0001_corr_sparcc.log; done

# Practice with just one file, one iteration- in ubuntu
python2.7 sparcc-master/SparCC.py euks.0001_simulated_datasets/permutation_0.txt -i 1 -c euks.0001_boot_corr/simulated_sparcc_0.txt -v euks.0001_boot_cov/simulated_sparcc_0.txt >> boot_euks.0001_corr_sparcc.log
# this worked (with just one iteration) but took awhile.

# Also try in local computer to see what's faster
python2.7 sparcc-master/SparCC.py  Network_PR2_OTUs/euks.0001_simulated_datasets/permutation_0.txt -i 1 -c  Network_PR2_OTUs/euks.0001_boot_corr/simulated_sparcc_0.txt -v  Network_PR2_OTUs/euks.0001_boot_cov/simulated_sparcc_0.txt >> boot_euks.0001_corr_sparcc.log
# started at 5:47- took forever
# just run in original loop in AWS- might take a few days

```
- This will take a long time. Started ~5:45pm on 1/21/19
	- Check progress by looking for the output files in ```euks.0001_boot_corr``` and ```euks.0001_boot_cov``` (there should be 100 of each generated over time)
	- As of 1/22/19 8:15 am, still running 
	- As 1/22/19 3:30pm, starting to see output files. 00, 01, and 02 only so far
	- 1/22/19 9:00pm, still only see 00, 01, and 02
	- 1/23/19 8:30am, up to file 04
	- 1/24/19 8:40am up to 07
	- 1/25 7pm up to file # 13.
	- 1/26 10:30am up to file # 16
	- 2/1 9:30am up to file # 38
	- 2/2 4:30 pm up to #42
	- 2/8 4:30 pm up to #64
	- **2/16 2:30pm up to #94....almost there**

- Next calculate the pseudo pvalues (one sided)
	
```
mkdir euks.0001_pvals

python2.7 sparcc-master/PseudoPvals.py euks.0001_basis_corr/cor_sparcc.out euks.0001_boot_corr/simulated_sparcc_#.txt 100 -o euks.0001_pvals/pvals_one_sided.txt -t one_sided
```

- Export all results to local computer. Import correlation file (```euks.001_basis_corr/cor_sparcc.out```) and p value file (```euks.0001_pvals/pvals_one_sided.txt```) into R.







## Jan 25, 2019
*Still waiting for high resolution Euk Network Correlation to run in SparCC on Amazon AWS. Meanwhile I made an OTU table with OTUs from all 3 domains. Retained OTUs with relative abundance > 0.01%. Run SparrCC on this matrix*.  


- Working directory ```/Users/admin/Google Drive/Wagner/Research/Manuscripts/Cariaco_Euks```)
- Text file with 3 domain OTU table = ```Network_PR2_OTUs/otu.table.001.txt```


```
mkdir Network_PR2_OTUs/otu.table.001_basis_corr

python2.7 sparcc-master/SparCC.py Network_PR2_OTUs/otu.table.001.txt --cor_file=Network_PR2_OTUs/otu.table.001_basis_corr/cor_sparcc.out
```

Output:   

```
reading data
Finished parsing table.
Table dimensions: (36,1327)
**** Data has been transposed! ****
First 3 column labels are :denovo180502 ,denovo80843 ,denovo217943
First 3 row labels are :A3a314A ,A2a237B ,A2b237A 

computing correlations
	Running iteration0
in slow
	Running iteration1
in slow
	Running iteration2
in slow
	Running iteration3
in slow
	Running iteration4
in slow
	Running iteration5
in slow
	Running iteration6
in slow
	Running iteration7
in slow
	Running iteration8
in slow
	Running iteration9
in slow
	Running iteration10
in slow
	Running iteration11
in slow
	Running iteration12
in slow
	Running iteration13
in slow
	Running iteration14
in slow
	Running iteration15
in slow
	Running iteration16
in slow
	Running iteration17
in slow
	Running iteration18
in slow
	Running iteration19
in slow
(1327, 1327)
writing results
wrote Network_PR2_OTUs/otu.table.001_basis_corr/cor_sparcc.out
wrote cov_mat_SparCC.out
Done!
```

- Next step is to make the bootstrapped datasets 

```
mkdir Network_PR2_OTUs/otu.table.001_simulated_datasets/

python2.7 sparcc-master/MakeBootstraps.py Network_PR2_OTUs/otu.table.001.txt -n 100 -t permutation_#.txt -p Network_PR2_OTUs/otu.table.001_simulated_datasets/
```


- Next run SparCC on each of these simulated files

```
mkdir Network_PR2_OTUs/otu.table.001_boot_corr
mkdir Network_PR2_OTUs/otu.table.001_boot_cov

for i in `seq 0 99`; do python2.7 sparcc-master/SparCC.py Network_PR2_OTUs/otu.table.001_simulated_datasets/permutation_$i.txt -c Network_PR2_OTUs/otu.table.001_boot_corr/simulated_sparcc_$i.txt -v Network_PR2_OTUs/otu.table.001_boot_cov/simulated_sparcc_$i.txt >> Network_PR2_OTUs/boot_otu.table.001_corr_sparcc.log; done


```


- Next calculated the pseudo pvalues
	- This will take a long time- probably several hours. 
	- Leave running over night ~7pm on Jan. 25th
	- 1/26 10:30am- husband closed my computer last night. Files go up to #60 but it looks like it is continuing to run now that computer is open. See if next file shows up soon.
		- Yes, it's still running!
		- Finished ~4:30pm

```
mkdir Network_PR2_OTUs/otu.table.001_pvals

python2.7 sparcc-master/PseudoPvals.py Network_PR2_OTUs/otu.table.001_basis_corr/cor_sparcc.out Network_PR2_OTUs/otu.table.001_boot_corr/simulated_sparcc_#.txt 100 -o Network_PR2_OTUs/otu.table.001_pvals/pvals_one_sided.txt -t one_sided
```

- Now move back into R for plotting



## Feb 1st-8th, 2019
**Playing with network object from R in Gephi for visualization**

- I downloaded a few packs in Gephi because I didn't like any of the default algorithms for plotting 
- Settled on the "Circle Pack Layout"
	- This moves the nodes around based on some attribute that you pick, like modularity or taxonomy
- Then I calculated the modularity so that I got 4 basic groups of nodes that had a lot of interconnections.
	- Made this the first plotting "hierarchy" in my network
	- Second hierarchy is taxonomy1 (domain), third hierarchy is taxonomy2, fourth hierarchy is taxonomy3, and 5th hierarchy is taoxnomy4 
- The nodes are color-coded by taxonomy label 4
	- The level of taxonomy divides the nodes by Syndiniales, Polycystinea, Dinophyceae, Thermoplasmatales, Flavobacteriales, ...Chromatiales, etc.
- Also filtered out some of the edges
	- Only retained edges for correlations >0.1505
		- Note- all these correlations are significant at the <0.001 level. These were filtered in R before even exporting.
	- This really highlights the interconnectdeness within the modules
	- At this point, note that the edges could be positive or negative correlations (can filter them out by type later)
- --> Export this image as ```3domain_network_circlepacklayout.pdf.svg```
	- NOTE to get this to work, had to uninstall the plugin in Gephi for exporting polygon shaped nodes. Found [here](https://github.com/gephi/gephi/issues/1759) that that was messing up the svg file export.


**Next focus in on each module to see which OTUs are there. My hypothesis is that this is dividing by depth/ oxygen.**.  

- In Gephi I can filter by Attribute>Parittion>Modularity Class
	- NOTE: From original total network figure,
		- Module 0 is center left
		- Module 1 is center right
		- Module 2 is top
		- Module 3 is bottom
- I exported each module as it's own svg image
	- Note: grey edges are negative correlations
	- Black edges are positive correlations
	    
- Start analyzing by identifying key taxa
	- Larger nodes are those with most connections
	- See which OTUs these are in Gephi and label in Inkscape
	      
- Open svg in Inkscape. Follow this [tutorial](https://wiki.cns.iu.edu/display/SCI2TUTORIAL/2.4+Saving+Visualizations+for+Publication) for making legends
	- Ungroup
	- Copy nodes of each color onto side to start making legend


**STOP THIS IS MESSY. GO BACK AND RE-RUN CORRELATIONS WITH ENV PARAMETERS**


- Working directory ```/Users/admin/Google Drive/Wagner/Research/Manuscripts/Cariaco_Euks```)
- Text file with 3 domain OTU table = ```Network_PR2_OTUs/otu.table.005.txt```


```
mkdir Network_PR2_OTUs/otu.table.005_basis_corr

python2.7 sparcc-master/SparCC.py Network_PR2_OTUs/otu.table.005.txt --cor_file=Network_PR2_OTUs/otu.table.005_basis_corr/cor_sparcc.out

mkdir Network_PR2_OTUs/otu.table.005_simulated_datasets/

python2.7 sparcc-master/MakeBootstraps.py Network_PR2_OTUs/otu.table.005.txt -n 100 -t permutation_#.txt -p Network_PR2_OTUs/otu.table.005_simulated_datasets/

mkdir Network_PR2_OTUs/otu.table.005_boot_corr
mkdir Network_PR2_OTUs/otu.table.005_boot_cov

for i in `seq 0 99`; do python2.7 sparcc-master/SparCC.py Network_PR2_OTUs/otu.table.005_simulated_datasets/permutation_$i.txt -c Network_PR2_OTUs/otu.table.005_boot_corr/simulated_sparcc_$i.txt -v Network_PR2_OTUs/otu.table.005_boot_cov/simulated_sparcc_$i.txt >> Network_PR2_OTUs/boot_otu.table.005_corr_sparcc.log; done

mkdir Network_PR2_OTUs/otu.table.005_pvals

python2.7 sparcc-master/PseudoPvals.py Network_PR2_OTUs/otu.table.005_basis_corr/cor_sparcc.out Network_PR2_OTUs/otu.table.005_boot_corr/simulated_sparcc_#.txt 100 -o Network_PR2_OTUs/otu.table.005_pvals/pvals_one_sided.txt -t one_sided


```

## Feb 9th, 2019
**Visualizing new network object from R in Gephi for visualization. These were made with OTUs >0.5% in any sample.**

- Modularity classes
	- Module 1 = center right
	- Module 7 = bottom
	- Module 0 = center left
	- Module 4 = top
	- Other modules consist of only 2, 3nodes. Won't focus on these
- Also saved screenshot of legend from Gephi
	- Color coded edges so negative correlations are grey, positice correlations are black. Thickness is weight of correlation
	- Only retained p-values <0.05

**Observations**

- Zoom in on module 7  (made svg file: ```3domain_network_circlepack_module7.svg```)
	- Env variables here are TZVS, H2S, NH4, NO2 --> more reduced forms of N and S
	- I'm going to go through manually and checking the clades in each node. In the list below, I am underlying the largest nodes, which are those with the highest degree (most connected edges)
	- Bacterial clades include 
		- Gammaproteobacteria
			- denovo130803 -B Collweliaceae
			- denovo290780 -B Oceanospirillales ZD0405
			- denovo186784 -B Salinisphaeraceae ZD0417 marine group
		- Deltaproteobacteria
			- <u>denovo293082 Bdellovibrionaceae OM27 clade</u> (Node ID 37)
			- denovo65050 Desulfobacteraceae SEEP-SRB1
		- Beta
			- denovo260142 Burkholderiaceae
		- Alpha
			- denovo216053 SAR11 clade Surface 1
		- Firmicutes
			- <u>denovo179427 Bacillales</u> (Node 161)
		- Spirochaetes
			- denovo156783 MSBL8
		- uncultured classes
			- denovo35025 Hyd24-12
			- denovo432995 Candidate WS3
		- Actinobacteria
			- denovo147540 Acidimicrobiales Sva0996 marine group
		- Chloroflexi
			- <u>denovo18923 SAR202 clade</u> (node ID 56)
			- denovo162656 SAR202 clade
		- Flavobacteria
			- denovo188771 Flavobacteriaceae Winogradskyella sp. PC-8
		-   Planctomycetes
			- denovo343001 Planctomyces   uncultured planctomycete
		-   Deferribacteres
			- denovo266202 SAR406 clade(Marine group A)
			- denovo50524 SAR406 clade(Marine group A)
	- Archaeal clades
		- Thaumarchaeota> Marine Group I
			- denovo95773
			- denovo140903
			- denovo12243
		- Thaumarchaeota > Group C3
			- denovo59484
			- denovo17773
			- denovo54443
		- Thaumarchaeota> Marine Benthic Group B
			- denovo182863
		- Euryarchaeota> Halobacteria>
			- denovo151320  Halobacteriales> Marine Hydrothermal Vent Group(MHVG)
		-   Euryarchaeota> Thermoplasmata
			- denovo189344 CCA47
			- denovo134192   Marine Group II
			- denovo70093  Marine Group II
			- denovo46022 Marine Group III
			- denovo94410 Marine Group II
			- denovo172503 Marine Group III
	- Eukaryotic Clades
		- Dinophyta> Dinophyceae *(This is the purple cluster on the left of center)*
			- denovo360862
			- denovo13939
			- denovo390752
			- <u>denovo418823</u> (Node ID 60)
			- denovo414723
			- denovo91312
			- denovo163006
			- denovo182436
			- denovo222790
			- denovo185876
			- denovo264593
		- Dinophyta> Syndiniales
			- Dino-Group-I
				- denovo264077 
				- denovo365415 
				- denovo48969 
				- denovo277382 
				- denovo396937 
			- Dino-Group-II
				- denovo239023
				- denovo96680 
				- denovo137116 
				- denovo360861 
				- denovo88790 
				- denovo392758  
				- <u>denovo280133</u> (Node ID 58)
				- denovo52459 
				- denovo227679  
				- denovo69987 
				- denovo427153 
				- denovo86282 
				- <u>denovo216884 </u> (Node ID 116)
				- denovo161618
				- denovo74320 
				- denovo60155 
				- denovo69987 
				- denovo330450 
				- denovo74320 
				- denovo84199 
				- denovo86282 
				- denovo427153 
				- denovo311118
		- Opisthokonta> Fungi> Ascomycota
			- denovo374491 Pezizomycotina> Sordariomycetes X sp.
			- denovo1692 Chytridiomycota> Chytridiomycotina> Chytridiomycetes X sp.
			- denovo39773 Chytridiomycota> Chytridiomycotina> Chytridiomycetes X sp.
				- *Quick google shows that Chytridiomycetes can be parasitic*
		- Alveolata
			- denovo210343 Perkinsea> Perkinsida
			- denovo36921 Ciliophora> Oligohymenophorea> Scuticociliatia> Cyclidiidae
		- Hacrobia
			- denovo137001 Prymnesiophyceae> Prymnesiales> Chrysochromulinaceae> Chrysochromulina> Chrysochromulina
		- Stramenopiles> Ochrophyta> Bacillariophyta> 
			- denovo447413 Raphid-pennate> Hippodonta> Hippodonta capitata
			- denovo138158 Polar-centric-Mediophyceae> Skeletonema> Skeletonema tropicum
		- Rhizaria> Radiolaria> 
			- <u>denovo342387 RAD-B> RAD-B-Group-IV</u> (Node ID 201)
			- denovo200689 Polycystinea> Spumellarida> Spumellarida-Group-I
			- denovo214671 Polycystinea> Spumellarida> Spumellarida-Group-I
			- denovo143989 Polycystinea> Spumellarida> Spumellarida-Group-I
			- denovo325994 Polycystinea> Spumellarida> Spumellarida-Group-I
		- Metazoa
			- denovo84770 Arthropoda> Crustacea> Maxillopoda> 
			- denovo442721 Ctenophora> Thalassocalyce> 
			- denovo31328 Annelida> Hirudo> Hirudo verbana
			- denovo283352 Platyhelminthes> Turbellaria> Seriata> Romankenkius
			- denovo445637 Cnidaria> Hydrozoa> Rhizophysa 
			- denovo217908 Cnidaria> Hydrozoa
			- denovo331754 Cnidaria> Hydrozoa> Vogtia pentacantha
	- I added some labels to the image, ```3domain_network_circlepack_module7_labelled.svg``` in inkscape in order to identify major clades/ clusters 
	- Next I am using the "Ego Network" filter to visually isolate the nodes with higest degrees (identified above- need to note "Node ID" in Data Lab in Gephi in order to do this filter)
		- For example, for denovo18923 SAR202 clade (node ID 56)- made image ```denovo18923_module7_egonet.svg``` that includes just this node and the nodes that connect to it, within Module 7
		- Did same for other nodes from Module 7 with high degrees, indicated above
		- *Continue with other 3 major modules but instead of noting each individual node, just identify the large ones (high degree) to save time.*
- Module 1 (svg file: ```3domain_network_circlepack_module1.svg```)
	- Environmental Variables are BNP, Fluorescence, CH4 (this is a weird mix but CH4 here is only correlated with one OTU in module)
	- Filtered by nodes with highest degrees (>=5). There are 15 nodes in this module with these high degrees:
		- denovo43029 Halobacteria> Halobacteriales (node ID 125) (Archaea)
		- denovo231149, GammaproteObacteria> Chromatiales (node ID 101)
		- denovo463083, GammaproteObacteria> Oceanosprillales (node ID 262)
		- BNP (Node ID 156)
		- denovo92826 Dinophyta> Syndiniales (node ID 129)
		- denovo138894 Dinophyta> Syndiniales (node ID 76)
		- denovo277655 Dinophyta> Syndiniales (node ID 91)
		- denovo230543 Dinophyta> Syndiniales (node ID 144)
		- denovo441155 Dinophyta> Syndiniales (node ID 115)
		- denovo376983 Dinophyta> Syndiniales (node ID 155)
		- denovo206994 Dinophyta> Syndiniales (node ID 270)
		- denovo411424 Dinophyta> Syndiniales (node ID 2)
		- denovo416105 Dinophyta> Dinophyceae (node ID 217)
		- denovo319829 Radiolaria> Polycystinea (node ID 71)
		- denovo415619 Radiolaria> Polycystinea (node ID 159)
			- Made an ego network image for each of the nodes listed above
			- Take home: maybe the polycystinea are the hosts of some of the syndiniales. There are a lot of assoications between these groups of OTUs
			- *At this point, I'm seeing drawbacks of presenting pos and neg correlations on the same graph. Consider presenting them separately (and calculation modularity separately). FInish the other 2 modules this way for now*
	



## Feb 15th 2019
Put details ofeach of the four major modules in powerpoint. But analyzing the ego network of each node in each module is getting cumbersome.

Instead look at whole network and find nodes that have high degrees:

- Filtered by nodes that have >= 5 degrees. Exported as .svg. Label with OTU ID and display in powerpoint.
	- 37	denovo293082
	- 161	denovo179427
	- 56	denovo18923
	- 60	denovo4188231
	- 58	denovo280133
	- 288 denovo239023
	- 116 denovo216884
	- 201	denovo342387
	- 120 denovo445637
	- 160 denovo442721
	- 61	denovo84770
	- 70 TZVS
	- 57	denovo17773
	- 26 denovo176118
	- 22	denovo186624
	- 326 denovo312514
	- 24 denovo339768
	- 141 Salinity
	- 80 denovo138115
	- 297 denovo2959
	- 24 denovo 339768
	- 221 denovo299729
	- 346 denovo308553
	- 125 denovo43029
	- 156 BNP
	- 101denovo231149
	- 262 denovo463083
	- 159 denovo415619
	- 71 denovo319829
	- 129 denovo92826
	- 76 denovo138894
	- 91 denovo277655
	- 144 denovo230543
	- 115 denovo441155
	- 155 denovo376983
	- 270 denovo206994
	- 2 denovo411424
	- 217 denovo416105
	- 98 denovo 348086
	- 242 denovo370279
	- 109 denovo139931
	- 135 denovo155760
	- 180 denovo260057
	- 333 denovo267690
	- 97 denovo328327
	- 171 denovo403524
	- 4 denovo422661
	- 40 denovo149934
	- 300 denovo199108
	- 183 denovo116564

- Went back to literature. 
	- Berry and Widder 2014: Deciphering microbial interactions and detecting keystone species with co-occurrence networks
		- "This very important result underscores that when co-occurrence networks are used to infer putative interactions samples should be drawn from similar environments in order to minimize the effects of habitat filtering or else the resulting network will suffer from a lack of interpretability."
		- Check Box 1 with recomendations
	- Banerjee et al. 2018: Keystone taxa as drivers of microbiome structure and functioning
		- Cite Berry and Widder and "recommend that the combined score of high mean degree, high closeness centrality and low betweenness centrality should be used as a threshold for defining keystone taxa in microbial communities" 
	- Sun et al. 2013: Core sediment bacteria drive community response to anthropogenic contamination over multiple environmental gradients
		- Assess a whole bunch of parameters for characteristics of networks- says something about connectedness of communities. Try something similarto compare across ecotype

- So I am going to re-do the analysis. Make 3 separate networks, one for each redox zone (From Suter et al 2018). and remake networks. Some of the parameters identified as important by Berry and Widder can be calculated in Gephi.




## Feb 16th 2019
- Finished making .txt files for each redox regime
- Run SparCC on each

Oxycline, input text file = ```otu.table.oxy.001.txt```

```
mkdir Network_PR2_OTUs/otu.table.oxy.001_basis_corr

python2.7 sparcc-master/SparCC.py Network_PR2_OTUs/otu.table.oxy.001.txt --cor_file=Network_PR2_OTUs/otu.table.oxy.001_basis_corr/cor_sparcc.out

mkdir Network_PR2_OTUs/otu.table.oxy.001_simulated_datasets/

python2.7 sparcc-master/MakeBootstraps.py Network_PR2_OTUs/otu.table.oxy.001.txt -n 100 -t permutation_#.txt -p Network_PR2_OTUs/otu.table.oxy.001_simulated_datasets/

mkdir Network_PR2_OTUs/otu.table.oxy.001_boot_corr
mkdir Network_PR2_OTUs/otu.table.oxy.001_boot_cov

for i in `seq 0 99`; do python2.7 sparcc-master/SparCC.py Network_PR2_OTUs/otu.table.oxy.001_simulated_datasets/permutation_$i.txt -c Network_PR2_OTUs/otu.table.oxy.001_boot_corr/simulated_sparcc_$i.txt -v Network_PR2_OTUs/otu.table.oxy.001_boot_cov/simulated_sparcc_$i.txt >> Network_PR2_OTUs/boot_otu.table.oxy.001_corr_sparcc.log; done

mkdir Network_PR2_OTUs/otu.table.oxy.001_pvals

python2.7 sparcc-master/PseudoPvals.py Network_PR2_OTUs/otu.table.oxy.001_basis_corr/cor_sparcc.out Network_PR2_OTUs/otu.table.oxy.001_boot_corr/simulated_sparcc_#.txt 100 -o Network_PR2_OTUs/otu.table.oxy.001_pvals/pvals_one_sided.txt -t one_sided

```


Shallow Anoxic, input text file = ```otu.table.shallanox.001.txt```

```
mkdir Network_PR2_OTUs/otu.table.shallanox.001_basis_corr

python2.7 sparcc-master/SparCC.py Network_PR2_OTUs/otu.table.shallanox.001.txt --cor_file=Network_PR2_OTUs/otu.table.shallanox.001_basis_corr/cor_sparcc.out

mkdir Network_PR2_OTUs/otu.table.shallanox.001_simulated_datasets/

python2.7 sparcc-master/MakeBootstraps.py Network_PR2_OTUs/otu.table.shallanox.001.txt -n 100 -t permutation_#.txt -p Network_PR2_OTUs/otu.table.shallanox.001_simulated_datasets/

mkdir Network_PR2_OTUs/otu.table.shallanox.001_boot_corr
mkdir Network_PR2_OTUs/otu.table.shallanox.001_boot_cov

for i in `seq 0 99`; do python2.7 sparcc-master/SparCC.py Network_PR2_OTUs/otu.table.shallanox.001_simulated_datasets/permutation_$i.txt -c Network_PR2_OTUs/otu.table.shallanox.001_boot_corr/simulated_sparcc_$i.txt -v Network_PR2_OTUs/otu.table.shallanox.001_boot_cov/simulated_sparcc_$i.txt >> Network_PR2_OTUs/boot_otu.table.shallanox.001_corr_sparcc.log; done

mkdir Network_PR2_OTUs/otu.table.shallanox.001_pvals

python2.7 sparcc-master/PseudoPvals.py Network_PR2_OTUs/otu.table.shallanox.001_basis_corr/cor_sparcc.out Network_PR2_OTUs/otu.table.shallanox.001_boot_corr/simulated_sparcc_#.txt 100 -o Network_PR2_OTUs/otu.table.shallanox.001_pvals/pvals_one_sided.txt -t one_sided

```



Euxinic, input text file = ```otu.table.eux.001.txt```

```
mkdir Network_PR2_OTUs/otu.table.eux.001_basis_corr

python2.7 sparcc-master/SparCC.py Network_PR2_OTUs/otu.table.eux.001.txt --cor_file=Network_PR2_OTUs/otu.table.eux.001_basis_corr/cor_sparcc.out

mkdir Network_PR2_OTUs/otu.table.eux.001_simulated_datasets/

python2.7 sparcc-master/MakeBootstraps.py Network_PR2_OTUs/otu.table.eux.001.txt -n 100 -t permutation_#.txt -p Network_PR2_OTUs/otu.table.eux.001_simulated_datasets/

mkdir Network_PR2_OTUs/otu.table.eux.001_boot_corr
mkdir Network_PR2_OTUs/otu.table.eux.001_boot_cov

for i in `seq 0 99`; do python2.7 sparcc-master/SparCC.py Network_PR2_OTUs/otu.table.eux.001_simulated_datasets/permutation_$i.txt -c Network_PR2_OTUs/otu.table.eux.001_boot_corr/simulated_sparcc_$i.txt -v Network_PR2_OTUs/otu.table.eux.001_boot_cov/simulated_sparcc_$i.txt >> Network_PR2_OTUs/boot_otu.table.eux.001_corr_sparcc.log; done

mkdir Network_PR2_OTUs/otu.table.eux.001_pvals

python2.7 sparcc-master/PseudoPvals.py Network_PR2_OTUs/otu.table.eux.001_basis_corr/cor_sparcc.out Network_PR2_OTUs/otu.table.eux.001_boot_corr/simulated_sparcc_#.txt 100 -o Network_PR2_OTUs/otu.table.eux.001_pvals/pvals_one_sided.txt -t one_sided

```



## Feb 16th-18th 2019
Finished SparCC runs, made network objects in R, then brought into Gephi. I have each network in a separate workspace in the Gephi file: ```otu.table.all3redox.001.graph.gephi```.  

- Shallow Anoxic is Workspace 9.  
- Oxycline is Workspace 3.  
- Euxinic is Workspace 4.

Start calculating some parameters of each network. Follow Berry and Widder and Sun et al. papers. I can do this easily with igraph in R (see script) and also in Gephi. 

Parameters calculated in iGraph (note that Gephi can also be used for some of these and the numbers agree):


Parameter         	| Oxycline         | Shallow Anoxic        | Euxinic               |
--------------------|------------------|-----------------------|-----------------------|
Mean degree (k) 	|2.577661                     |18.15563      |67.87257           |
Avg shortest path length (l) |8.058858         |2.533821      | 1.880004          |
Mean clustering coefficient (C)|0.002647838    |0.03653293    | 0.1338309         |
Network Diameter|19    |4    | 3         |
Mean Betweenness centrality (CB)|3782.003      |462.875       |248.9788           |
Closeness centrality (CC) |2.128513e-05        |0.00295705    |0.002641891        |


- Explanations:
	- **Mean degree** = the degree of a node counts the number of edges it has.
	- **Shortest path** = the shortest path between any two nodes is the single path with fewest edges between them. 
	- **Clustering coefficient** = (aka transitivity) the probability that the adjacent vertices of a vertex are connected
		- First 3 parameters above can be used to say something about overall network connectivity.
	- **Betweenness centrality** = the betweenness centrality of a node is equal to the number of shortest paths between any two nodes in the graph passing through that node.
		- Eg. of all the shortest paths between any two nodes, how many of those go through this node?
		- According to Berry and Widder, this should below for keystone OTUs, which seems counterintuitive.
		- It also seems like igraph and Gephi are using different algorithms for this because I get similar but different results for individiual OTUs. Each lists the same (more or less) OTUs as having the lowest betweenness values but the values are not precisely equivalent and the order is different. Use igraph results since that is what Berry and Widder use to say something about keystoneness.
	- **Closeness centrality** = the average distance of this node to any other node.
		- Eg., how many steps is required to access every other vertex from a given vertex
		- According to Berry and Widder, this should be high for keystone OTUs. 
			- This also seems counterintuitive but taken together with the point about betweenness I guess this means that while keystones are highly connected to many other nodes (high mean degree) they are not necessarily on the shortest paths between nodes and therefore have high avg distance to other nodes. Non-keystones might have low average distances to their neighbors but are not highly connected to many neighbors.



- Also looked at degree distribution and graph distance reports in Gephi. See embedded reports below:

# Oxycline Network, Degree Distribution:

![Oxy Deg Dist](file:///Users/admin/Google%20Drive/Wagner/Research/Manuscripts/Cariaco_Euks/Network_PR2_OTUs/otu.table.oxy.001.degree.distribution.pdf)

# Oxycline Network, Graph Distance Report:

![Oxy Deg Dist](file:///Users/admin/Google%20Drive/Wagner/Research/Manuscripts/Cariaco_Euks/Network_PR2_OTUs/otu.table.oxy.001.degree.distribution.pdf)

# Shallow Anoxic Network, Degree Distribution:

![ShallAnox Deg Dist](file:///Users/admin/Google%20Drive/Wagner/Research/Manuscripts/Cariaco_Euks/Network_PR2_OTUs/otu.table.shallanox.001.degree.distribution.pdf)

# Shallow Anoxic Network, Graph Distance Report:

![ShallAnox Deg Dist](file:///Users/admin/Google%20Drive/Wagner/Research/Manuscripts/Cariaco_Euks/Network_PR2_OTUs/otu.table.shallanox.001.graphdist.pdf)

# Euxinic Network, Degree Distribution:

![Eux Deg Dist](file:///Users/admin/Google%20Drive/Wagner/Research/Manuscripts/Cariaco_Euks/Network_PR2_OTUs/otu.table.eux.001.degree.distribution.pdf)

# Euxinic Network, Graph Distance Report:
![Eux Deg Dist](file:///Users/admin/Google%20Drive/Wagner/Research/Manuscripts/Cariaco_Euks/Network_PR2_OTUs/otu.table.eux.001.graphdist.pdf)







- I also extracted lists of node from each network with the highest degrees, lowest betweeness centralities, and highest closeness centralities. Then picked out OTUs which appear in all 3 lists for each network (all details in R script).
- Exported network plots for each redox regime.
	- All 3 are plotted with same algorithm (Force Atlas 2) with same parameters
		- Tolerance 1
		- Uncheck approximate repulsion
		- Approximation 1.2
		- Scaling 7
		- Stonger gravity unchecked
		- Gravity 10
		- Dissuade hubs, Linlog mode, and Prevent Overlap all unchecked
		- Edge weight influence 1.0
	- Made .svg files of networks. Tonight spend time fixing legends
-  Using suggestionsfrom Berry and Widder, I extracted candidate keystone OTUs from each network (see R script)
	- Oxycline
		- None! No single OTU satisfied all 3 parameter cutoffs
	- ShallowAnoxic
		- "denovo158333 -B" 167	 Bacteria	  Candidate division WS3 uncultured bacterium
		- "denovo64710 -E"  235	 Eukaryota Alveolata Dinophyta Syndiniales Dino-Group-I Dino-Group-I-Clade-1 Dino-Group-I-Clade-1 X Dino-Group-I-Clade-1 X sp. strain10
		- "denovo431262 -E" 253 Eukaryota	Alveolata Perkinsea Perkinsida Perkinsida X	Perkinsida XX Perkinsida XXX Perkinsida XXX sp.
		- "denovo244612 -B" 308 Bacteria Proteobacteria Gammaproteobacteria Chromatiales Ectothiorhodospiraceae Thiorhodospira uncultured bacterium
		- "denovo49812 -E" 410 	Eukaryota	Alveolata Dinophyta Dinophyceae Dinophyceae X Dinophyceae XX Dinophyceae XXX Dinophyceae XXX sp.
		- "denovo294767 -B" 520 Bacteria Proteobacteria Betaproteobacteria Burkholderiales Oxalobacteraceae Herbaspirillum uncultured bacterium
	- Euxinic
		- "denovo21694 -E" 167 	Eukaryota	Rhizaria Cercozoa Filosa-Thecofilosea Filosa-Thecofilosea X Filosa-Thecofilosea XX Filosa-Thecofilosea XXX Filosa-Thecofilosea XXX sp.
		- "denovo92288 -A" 299 denovo92288 -A Archaea Thaumarchaeota Miscellaneous Crenarchaeotic Group uncultured crenarchaeote
		- "denovo48360 -A"  76 Archaea	  Euryarchaeota Halobacteria Halobacteriales	  Deep Sea Euryarchaeotic Group(DSEG) uncultured euryarchaeote
		- "denovo161618 -E" 348 	Eukaryota Alveolata Dinophyta Syndiniales Dino-Group-II Dino-Group-II-Clade-14 Dino-Group-II-Clade-14 X Dino-Group-II-Clade-14 X sp.		- "denovo288601 -E"	380 Eukaryota	 Stramenopiles Ochrophyta Bacillariophyta Bacillariophyta X Raphid-pennate	Hippodonta Hippodonta capitata		- "denovo106090 -E" 283 Eukaryota	 Alveolata Dinophyta Syndiniales Dino-Group-II	Dino-Group-II X Dino-Group-II XX Dino-Group-II XX sp.
		- "denovo165566 -A" Archaea Thaumarchaeota Group C3 uncultured archaeon
- Next make ego networks of the keystone species in the shallow anoxic network and see who they are connected to. Use Node ID as file label (listed above- 118, 227, etc) because this is how Gephi filters ego networks
	- Also make ego network for some interesting groups
		- 128, Cariacotrichea
		- 185, 89, 576 Myxococcales (see [Munoz-Dorado et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4880591/) to read about predatory behavior)
		- 484, 420, 431, 20, 200 Bdellovibrionales

		


## Feb 21st- 22nd 2019
Open questions:

- Diversity of individual Syndiniales clades with depth/ size fraction? 
	- Bubble Plot of clades
	- PA may be active infection
	- FL could be a lot of spores
	- Orsi aruges that extent of Syn is constrained by their host, specificially Radiolarian hosts 
		- Check Euk network to pick out hosts
- Associations of ciliates with proks, indicating possible symbioses? 90% had epibionts (Orsi et al.) In addition, many had bacterial prey in vacuoles
	- Ciliate Node IDs in Shallow anoxic sample are 7, 128, 186, 241, 250, 281, 315, 370, 417, 422, 446, 491, 526, 598. Think about how to plot/ show this
	- Also specifically check for associations with methanogens (especially MG groups and Group C3 [now MCG-15])


Also- SparCC Analysis of the Euk-only network is done. Download from EC2 using Filezilla

- Open file browser for this volume in Filezilla
	- Host is the public DNS from Amazon EC2 (ec2-3-86-89-79.compute-1.amazonaws.com), the user name is ubuntu, port is 22, and leave password empty. Load key file as password (Settings/ SFTP)
- Copied files from ```~/Cariaco/Network_Analysis/Euk_Network```to my external hard drive, WD Passport, in folder called: ```CariacoEukPaper⁩``` file folders starts with ```euks.0001``` (don't confuse with files from previous analysis, ```euks.001```, which are on my computer but not on WD Passport
	- Since I have files on my hard drive, remove from the attached volume in the EC2 instance (otherwise I will be paying for these). From ```Cariaco``` directory, used the "dangerous" force remove a directory command ```rm -rf Network_Analysis```
	- Stopped the instance (didn't terminate because still need info there for nitrogen and sulfur papers)


