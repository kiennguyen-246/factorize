import json
from itertools import chain
from dimod import child_structure_dfs
from dwave.embedding import EmbeddedStructure
from dwave.system import DWaveSampler
from minorminer import minorminer
from minorminer.utils import DisconnectedChainError


def get_embedding(bqm, output):
    best_embedding = None
    source_edgelist = list(chain(bqm.quadratic, ((int(v), int(v)) for v in bqm.linear)))
    target_edgelist = child_structure_dfs(DWaveSampler()).edgelist
    min_val = 1e9
    min_len = 1e5
    trials = 0
    found = False
    while trials < 50 or not found:
        trials += 1
        embedding = minorminer.find_embedding(source_edgelist, target_edgelist, verbose=2, max_no_improvement=20,
                                              random_seed=trials, chainlength_patience=20, timeout=300)
        if len(embedding.keys()) == 0:
            print("Failed", trials)
            continue

        try:
            embed_structure = EmbeddedStructure(target_edgelist, embedding)
        except DisconnectedChainError:
            print("Failed", trials)
            continue
        found = True

        # def is_disconnected(embedding):
        #     for chain in embedding.values():
        #         for i in range(0, len(chain)):
        #             for j in range(i + 1, len(chain)):
        #                 if (chain[i], chain[j]) not in target_edgelist and (chain[j], chain[i]) not in target_edgelist:
        #                     return True
        #     return False
        #
        # if is_disconnected(embedding):
        #     print("Failed", trials)
        #     continue

        len_embedding = max(map(len, embedding.values()))
        val_embedding = sum([a for a in map(len, embedding.values())])
        print(len_embedding, val_embedding, trials)
        if (len_embedding < min_len) or (len_embedding == min_len and val_embedding < min_val):
            min_len = len_embedding
            min_val = val_embedding
            best_embedding = embedding
    print(min_len, min_val)
    embed_structure = EmbeddedStructure(target_edgelist, best_embedding)
    with open(output, 'w') as f:
        json.dump(embed_structure, f, indent=4)
    return embed_structure
