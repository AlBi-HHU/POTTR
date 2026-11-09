## Reproducing our results

Follow these instruction to reproduce our results.
We ran the script `get_output_sizes.py`, which executes POTTR on the selected dataset for every `k=2` up to the given `max-k`.
Use `--max-k 89` for MASTRO's lung data, `--max-k 120` for MASTRO's AML data,  `--max-k 394` for the TRACERx data, and `--max-k 12` for the TLS and TLSCL data.
The script will run until the recurrent trajectory only consists of the root node.

```
python get_output_sizes.py --pottr-path <absolute_path_to_code_dir> --dags ../data/<specify_data> --max-k <k> --output-path ../data/output --pool-size 0 --cores 128
```

The script produces an output folder using the current time that has a subfolder for each run of `k`. `--pool-size 0` let's POTTR run until no further solution is found for the current `k`.
After running `get_output_sizes.py`, use [summarize_trajectories.py](<../scripts/summarize trajectories/summarize_trajectories.py>) to get a visual summary of all the unique trajectories identified by POTTR.
