import random
import re


def shuffle_da_list(a_list, zero=True):
    mutations_copy = a_list.copy()
    random.shuffle(mutations_copy)

    if zero:
        mutations_copy = ['0'] + mutations_copy
    print(mutations_copy)
    return mutations_copy


def replace_match(match):
    return replacements[match.group(0)]


graph_ids = '47 72\n'
graph = '0->-DNMT3A NPM1->-FLT3 NPM1->-NRAS DNMT3A->-NPM1'
graph_split = graph.split(' ')
mutations = set()
for edge in graph_split:
    if '-?-' in edge:
        a, b = edge.split('-?-')
        mutations.add(a)
        mutations.add(b)
    if '->-' in edge:
        a, b = edge.split('->-')
        mutations.add(a)
        mutations.add(b)

mutations = list(mutations)
if '0' in mutations:
    mutations.remove('0')

replaced_strings = []
for i in range(10):
    add_zero = ('0' in mutations)
    shuffled_list = shuffle_da_list(mutations, add_zero)
    tmp = zip(mutations, shuffled_list)
    replacements = dict(tmp)

    pattern = re.compile("|".join(re.escape(k) for k in replacements.keys()))

    result = pattern.sub(replace_match, graph)
    replaced_strings.append(result)

with open('test_graphs.txt', 'w') as file:
    file.write(graph + ' (2)\n')
    file.write(graph_ids)
    for entry in replaced_strings:
        file.write(entry + ' (2)\n')
        file.write(graph_ids)







