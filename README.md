# Hairpin
## Misc scripts/tools for working with UKB csv files
- [Running Hairpin](#Usage-for-PGS-analysis)
  * [Command Line Arguments](#command-line-arguments)
  * [Sample Calls](#sample-calls)
  * [Possible Errors](#possible-errors)
- [PGS Specifications]($PGS-specs)
- [Plotting](#plotting)
- [References](#references)

Code here written by Ophelia Dominguez, Shevaughn Holness, Andy Dahl, with help from Renée Fonseca and Manuela Costantino.

### Usage for PGS analysis
TODO: an actual description of this lol.
#### Command Line Arguments
To intialize hairpin with a new phenotype, run `workflow.sh` with the following positional arguments: 
- `new`: input `1` or `0` corresponding to a new phenotype run or not 
- `phen_name`: the name you want to call this phenotype. e.g. `height`
- `phen_id`: the phenotype id as it will show up in the .pheno file. e.g. `34.0.0`
- `phen_path`: the full path to your .pheno file
- `sum_file`: the name of your summary statistics file (TODO: make this also an absolute path thing)
- `binary`: is the phenotype binary? Input `T` or `F`
- `phen_rel_path`: the relative path of the phenotype in the extracted phenotypes folder. (TODO: get rid of this, redundant but will have to edit all the scripts)
- `stat_col`: the name of the beta statistic column in your summary statistics file.
- `ID_col`: the name of the snpID column in you summary statistics file. Put `not_avail` if not in your file.
- `chr_col`: the name of the chromosome column in your summary statistics file.
- `bp_col`: the name of the column with the variant's position in your summary statistics file. Put `not_avail` if not in your file.
- `a1`: the effect allele column in your summary statistics file.
- `a2`: the other allele column in your summary statistics file.
- `pcol`: the column containing the pvalues in your summary statistics file.
- `stat`: statistics type used for the PGS. (TODO: get rid of this? Will just be beta every time unless we want to play around?)
- `maf`: the column containing the allele frequency in your summary statistics file. Add `:0.05` to the end, e.g. `EAF:0.05`

#### Sample calls
Running Hairpin without a new phenotype:
- `sbatch workflow.sh 0`
Running Hairpin with days of physical activity as a new phenotype:
- `sbatch /gpfs/data/ukb-share/dahl/ophelia/hairpin/workflow_test/workflow.sh 1 days_phys X884.0.0 /gpfs/data/ukb-share/extracted_phenotypes/days_phys/days_phys.pheno ACTIVITY1_single_p5e-8_sumstats.txt F days_phys/days_phys.pheno BETA SNPID CHR BP EFFECT_ALLELE OTHER_ALLELE PVALUE beta EAF:0.05`

#### Possible Errors
TODO: everything lol. Writing this like it will become an actual package but if so the thing will break as soon as it's not on a HPC cluster.

### PGS Specs:

In hairpin/txt_files there are txt files for both the pcs and pvalues. These are new line delimited files that all the code reads in and is relative to. To add a new number of pcs or p-value threshold, just add a new line with the desired value. For formatting purposes please keep the numerically ascending order.
Then rerun the formatting scripts found in the txt_files directory.

Feel free to edit the PGS specifications found in `submissions/hairpin_master.sh`. The current specs are: TODO

### Plotting

TODO

### References
- Becker, J., Burik, C.A.P., Goldman, G., Wang, N., Jayashankar, H., Bennett, M., Belsky, D.W., Karlsson Linnér, R., Ahlskog, R., Kleinman, A., Hinds, D.A., 23andMe Research Group, Caspi, A., Corcoran, D.L., Moffitt, T.E., Poulton, R., Sugden, K., Williams, B.S., Harris, K.M., Steptoe, A., Ajnakina, O., Milani, L., Esko, T., Iacono, W.G., McGue, T., Magnusson, P.K.E., Mallard, T.T., Harden, K.P., Tucker-Drob, E.M., Herd, P., Freese, J., Young, A., Beauchamp, J.P., Koellinger, P.D., Oskarsson, S., Johannesson, M., Visscher, P.M., Meyer, M.N., Laibson, D., Cesarini, D., Benjamin, D.J., Turley, P., and Okbay, A. (2021). Resource Profile and User Guide of the Polygenic Index Repository. Nature Human Behaviour. Published online June 17. doi:10.1038/s41562-021-01119-3.
- Okbay, A., Wu, Y., Wang, N. et al. Polygenic prediction of educational attainment within and between families from genome-wide association analyses in 3 million individuals. Nat Genet 54, 437–449 (2022). https://doi.org/10.1038/s41588-022-01016-z
- Becker, J., Burik, C.A.P., Goldman, G. et al. Resource profile and user guide of the Polygenic Index Repository. Nat Hum Behav 5, 1744–1758 (2021). https://doi.org/10.1038/s41562-021-01119-3
- Shungin, D., Winkler, T., Croteau-Chonka, D. et al. New genetic loci link adipose and insulin biology to body fat distribution. Nature 518, 187–196 (2015). https://doi.org/10.1038/nature14132
- Global Lipids Genetics Consortium. Discovery and refinement of loci associated with lipid levels. Nat Genet 45, 1274–1283 (2013). https://doi.org/10.1038/ng.2797
- Loic Yengo, Sailaja Vedantam, Eirini Marouli, …, Yukinori Okada, Andrew R. Wood, Peter M. Visscher, Joel N. Hirschhorn. A Saturated Map of Common Genetic Variants Associated with Human Height from 5.4 Million Individuals of Diverse Ancestries. bioRxiv 2022.01.07.475305; doi: https://doi.org/10.1101/2022.01.07.475305
- Locke, A., Kahali, B., Berndt, S. et al. Genetic studies of body mass index yield new insights for obesity biology. Nature 518, 197–206 (2015). https://doi.org/10.1038/nature14177
- Choi SW, and O’Reilly PF. "PRSice-2: Polygenic Risk Score Software for Biobank-Scale Data." GigaScience 8, no. 7 (July 1, 2019). https://doi.org/10.1093/gigascience/giz082.