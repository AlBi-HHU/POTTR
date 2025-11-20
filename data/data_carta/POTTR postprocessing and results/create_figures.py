import os
import re

import networkx as nx
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

# === CONFIG ===
base_dir = "path/to/Output/cell_differentiation_maps/output_TLS-CL_2/"
png_out = "path/to/Output/cell_differentiation_maps/output_TLS-CL_2/pngs_all"
output_pdf = "all_graphs.pdf"

os.makedirs(png_out, exist_ok=True)

def draw_graph(G, out_path):
    dot = nx.drawing.nx_pydot.to_pydot(G)
    dot.write_png(out_path)

# === Main ===
seen = list()
pngs = []

folders = []
for i in range(11, 1, -1):
    folders.append(base_dir + 'out_k' + str(i))

for folder in folders:
    folder_path = os.path.join(folder, 'relabeled_graphs')
    for file in os.listdir(folder_path):
        if file.endswith(".gexf") and not file.startswith('.'):
            path = os.path.join(folder_path, file)
            G = nx.read_gexf(path)
            edge_set = set(G.edges())
            if edge_set in seen:
                print('Skip ', file)
                continue
            seen.append(edge_set)
            number = re.search('out_k(\d+)', folder).group(1)
            out_file = file.replace('.gexf', '')
            out_png = os.path.join(png_out, f"{number}_{out_file}.png")
            draw_graph(G, out_png)
