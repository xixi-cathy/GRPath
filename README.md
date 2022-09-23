<!-- ABOUT GRPath -->

## About GRPath

**GRPath** is a computational framework that calculates the genotype-to-phenotype regulation paths in the form of “putative causal region (pcRegion) - putative causal variant (pcVariant) - putative causal gene (pcGene) - noteworthy cell type - disease phenotype” for diseases.


**Reference:**
Xi, X., Li, H., Chen, S., Lv, T., Ma, T., Jiang, R., Zhang, P., Wong, W.H., Zhang, X., 2022. Unfolding the genotype-to-phenotype black box of cardiovascular diseases through cross-scale modeling. iScience 25, 104790. [(Link)](https://doi.org/10.1016/j.isci.2022.104790)



<img src=".\figure\fig.jpg" alt="overview" style="zoom:10%;" />



<!-- AUTHOR -->

## Author

Xi Xi xix19@mails.tsinghua.edu.cn



<!-- DEPENDENCIES -->

## Dependencies

We recommend using GRPath in Anaconda environment, where compatible Python and R packages should be installed. You may download the **requirement.txt** file to find out all the dependencies needed.



<!-- EXAMPLE -->

## Example

You may download the code and demo_data, and run:

```sh
bash predict_region_variant_egene.sh
bash predict_regulation_path.sh
```

to find the pcPaths for heart failure caused by coronary heart disease (cHF).
