import json
import os
import time

from dwave.system import EmbeddingComposite, FixedEmbeddingComposite
from dwave.system.samplers import DWaveSampler
from dwave.embedding.chain_strength import uniform_torque_compensation
from minorminer.utils import DisconnectedChainError
from dwave.samplers import SimulatedAnnealingSampler
import dwave.inspector

from .export_embedding import get_embedding

def solve_quantum_annealing(bqm,
                            method="?_",
                            num_reads=1000):
    chain_strength_prefactor = 0.25
    annealing_time = 200
    anneal_schedule_id = -1
    chain_strength = uniform_torque_compensation(
        bqm=bqm, prefactor=chain_strength_prefactor)
    schedules = {12: [(0.0, 0.0), (40.0, 0.4), (180.0, 0.4), (200.0, 1.0)],
                 11: [(0.0, 0.0), (40.0, 0.5), (120.0, 0.5), (200.0, 1.0)],
                 13: [(0.0, 0.0), (40.0, 0.5), (130.0, 0.5), (200.0, 1.0)],
                 14: [(0.0, 0.0), (30.0, 0.5), (160.0, 0.5), (200.0, 1.0)]}
    embed_config = method + str(num_reads) + "-" + str(chain_strength)
    solver_config = method + str(num_reads) + "-" + str(chain_strength) + "s" + str(anneal_schedule_id) + "_A" + str(
        annealing_time)
    os.environ["SOLVER_CONFIG"] = solver_config
    data_name = os.getenv("INPUT_FILE")
    # data_name = ""
    output_dir = "output/" + data_name + "/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    get_embedding(bqm, output_dir + embed_config + ".json")
    with open(output_dir + embed_config + ".json", "r") as f:
        embedding = json.load(f)
    embedding = {int(k): v for k, v in embedding.items()}
    try:
        sampler = FixedEmbeddingComposite(DWaveSampler(), embedding=embedding)
    except DisconnectedChainError as e:
        raise RuntimeError(e)
    # sampler = EmbeddingComposite(DWaveSampler())

    start = time.time()
    if anneal_schedule_id == -1:
        response = sampler.sample(bqm=bqm,
                                  chain_strength=chain_strength,
                                  num_reads=num_reads,
                                  label=solver_config,
                                  annealing_time=annealing_time)
    else:
        schedule = schedules[anneal_schedule_id]
        response = sampler.sample(bqm=bqm,
                                  chain_strength=chain_strength,
                                  num_reads=num_reads,
                                  label=solver_config,
                                  # annealing_time=annealing_time,
                                  anneal_schedule=schedule)
    end = time.time()
    dwave.inspector.show(response)
    chains = response.info["embedding_context"]["embedding"].values()
    # for chain in chains:
    #     if len(chain) > 10:
    #         print(chain)

    config_dict = {
        "config": solver_config,
        "num_vars": len(bqm.variables),
        "num_qubit": sum([len(chain) for chain in chains]),
        "time_elapsed": end - start,
        "best_state": {
            "sample": response.record.sample[0].tolist(),
            "energy": response.record.energy[0],
            "chain_break_fraction": response.record.chain_break_fraction[0],
        },
        "chain_strength_prefactor": chain_strength_prefactor,
        "chain_strength": chain_strength,
        "max_chain_length": max([len(chain) for chain in chains]),
        "timing_info": response.info["timing"],
        "embedding_info": response.info["embedding_context"]
    }
    with open(output_dir + solver_config + "_1.json", "w") as f:
        json.dump(config_dict, f, indent=4)
    # file_response_output = open(output_dir + solver_config + ".txt", "w")
    # file_response_output.write("Config: " + solver_config + "\n")
    # file_response_output.write("Number of source variables: " + str(len(bqm.variables)) + "\n")
    # file_response_output.write("Number of target variables: " + str(sum([len(chain) for chain in chains])) + "\n")
    # file_response_output.write("Time Elapsed: " + str(end - start) + "\n")
    # file_response_output.write(
    #     "Best State: " + str(response.record.sample[0]) + "\n" + str(response.record.energy[0]) + "\t" + str(
    #         response.record.chain_break_fraction[0]) + "\n")
    # file_response_output.write("ChainStr/ChainLen: " + str(chain_strength) + "/" + str(
    #     max([len(chain) for chain in chains])) + "\n")
    # file_response_output.write("Info: " + str(response.info["timing"]) + "\n")
    # file_response_output.write("Embedding Info: " + str(response.info["embedding_context"]) + "\n")

    return response


def solve_simulated_annealing(bqm, method="?_", num_reads=1000):
    sampler = SimulatedAnnealingSampler()
    # beta_range = [0.1, 4]
    num_sweeps = 100
    # config = method + str(num_reads) + "-SA" + "".join(str(beta_range).split(" ")) + "s" + str(num_sweeps)
    solver_config = method + str(num_reads) + "-SA" + "s" + str(num_sweeps)
    os.environ["SOLVER_CONFIG"] = solver_config
    data_name = os.getenv("INPUT_FILE")
    output_dir = "output/" + data_name + "/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # print(config)
    start = time.time()
    response = sampler.sample(bqm,
                              num_reads=num_reads,
                              label=solver_config,
                              # beta_range=beta_range,
                              num_sweeps=num_sweeps)
    end = time.time()
    config_dict = {
        "config": solver_config,
        "num_vars": len(bqm.variables),
        "time_elapsed": end - start,
    }
    with open(output_dir + solver_config + "_1.json", "w") as f:
        json.dump(config_dict, f, indent=4)
    return response
