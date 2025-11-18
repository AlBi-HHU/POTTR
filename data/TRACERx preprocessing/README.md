# Preprocessing the TRACERx data to run with POTTR

The phylogenetic trees from TRACERx are available at [zenodo](https://zenodo.org/record/7822002).
Download and unpack the zip file from zenodo first.
We need the files `20221109_TRACERx421_phylogenetic_trees.rds`, `20221109_TRACERx421_mutation_table.rds`, and `20221109_TRACERx421_all_patient_df.rds`.
The script to extract the trees and write them to separate files is written in R.
Execute the R script in the [data/TRACERx preprocessing](data/TRACERx%20preprocessing) directory and provide the path to the TRACERx data:

```shell
Rscript create_tree_files.r </path/to/TRACERx-figurecode> processed_data/
```

This will create the directory `processed_data` with all trees stored as xlsx files.
Next, execute [create_tracerx_trees.py](create_tracerx_trees.py):

```shell
python create_tracerx_trees.py -i processed_data/ -o data_tracerx/ -mut </path/to/TRACERx-figurecode/data/LCCE_converted_data/20221109_TRACERx421_mutation_table.rds> -s Hugo+Type -drivers
```

This will create all TRACERx trees as txt files using the MASTRO format. See the [MASTRO repository](https://github.com/VandinLab/MASTRO/tree/main) for details.

### Program arguments

The following arguments are available:

| Argument    | Argument long form      | Description                                                                                                                                                                                                                                                |
|-------------|-------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| -h          | --help                  | show help message and exit                                                                                                                                                                                                                                 |
| -i <path>   | --input <path>          | Path to input files                                                                                                                                                                                                                                        |
| -o <path>   | --output <path>         | Path to store output files                                                                                                                                                                                                                                 |
| -mut <path> | --mutation-table <path> | Path to mutation table from TRACERx (20221109_TRACERx421_mutation_table.rds)                                                                                                                                                                               |
| -s <scheme> | --scheme <scheme>       | Labelling scheme of nodes. Choose from Hugo (only Hugo symbol, first occurrence, no duplications), Hugo+Cluster (Hugo symbol and cluster number), or Hugo+AAChange (Hugo symbol and mutation level information), Hugo+Type (Hugo symbol and mutation type) |
| -drivers    | --only-drivers          | If set, mutations will be filtered to only keep driver mutations. Default, mutations will not be filtered.                                                                                                                                                 |
| -g <path>   | --gene-list <path>      | Path to tsv file containing genes for filtering                                                                                                                                                                                                            |

For our results, we used the scheme `Hugo+Type` and filtered to keep only driver mutations.

For a larger data set, instead of setting the driver gene flag, you can filter the data by known oncogenes. See for example the list of genes in [OncoKB_Cancer_Gene_List.tsv](OncoKB_Cancer_Gene_List.tsv) provided by Ahmed Shuaibi on [GitHub](https://github.com/raphael-group/dialect/tree/main).