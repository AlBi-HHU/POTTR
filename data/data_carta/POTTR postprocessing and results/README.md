Commands to create the cell differentiation results:

```shell
python run_POTTR.py -d ../Data/data_carta/graphs_TLS-CL_somite_collapsed -o ../output/cell_differentiation_maps/output_TLS-CL/ -k <2-11> -v -parallel -pool 5000 -c 100
```

```shell
python run_POTTR.py -d ../Data/data_carta/graphs_TLS_somite_collapsed -o ../output/cell_differentiation_maps/output_TLS/ -k <2-12> -v -parallel -pool 5000 -c 100
```

Execute for all values from k=2 up to 11 and 12 for TLSCL and TLS, respectively.

After running POTTR on the provided cell differentiation maps, use the jupyter notebook [relabel_trajectories](relabel_trajectories.ipynb) to relabel nodes with cell types. After that, use the python script [create_figures](create_figures.py) to get png figures of the trajectories.