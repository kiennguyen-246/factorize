import json
import os
from math import log


class Formulation:
    def __init__(self) -> None:
        self.var_map = {}
        self.var_cnt = 0
        self.fixed_variables = {}
        self.fixed_var_map = {}

    def solve(self, n) -> (int, int):
        """
        The base solver for any formulations
        :param n: A safe semiprime
        :return: A tuple of two integers factorized from n
        """
        raise NotImplementedError

    # Helper functions
    def reduce_deg(self, coef=1, term=()):
        """
        Recursively reduce a HUBO term until it becomes QUBO terms.
        The formulation for optimization is: x1x2x3 = x3s + 2x1x2 - 4x1s - 4x2s + 6s
        :param term: The coefficient and a tuple of variables
        :return: A list of optimized terms
        """
        # print(coef, term)
        if len(term) <= 2:
            return {term: coef}, 0
        self.var_map[("s", self.var_cnt)] = self.var_cnt
        self.var_cnt += 1
        new_terms = {}
        new_offset = 0
        # prefix = ()
        # x = term[:-2]
        prefix = term[:-3]
        x = term[-3:-2]
        y = term[-2:-1]
        z = term[-1:]
        w = (self.var_map[("s", self.var_cnt - 1)],)
        if coef < 0:
            new_terms[x + w] = coef
            new_terms[y + w] = coef
            new_terms[z + w] = coef
            new_terms[w] = -2 * coef
        else:
            new_terms[x + w] = coef
            new_terms[y + w] = coef
            new_terms[z + w] = coef
            new_terms[w] = -1 * coef
            new_terms[x + y] = coef
            new_terms[y + z] = coef
            new_terms[z + x] = coef
            new_terms[x] = -coef
            new_terms[y] = -coef
            new_terms[z] = -coef
            new_offset = coef
        if len(term) > 3:
            new_terms = {prefix + new_term: new_terms[new_term] for new_term in new_terms}
            new_terms[prefix] = new_offset
            new_offset = 0
        # ans_terms = {}
        # ans_offset = new_offset
        # for new_term in new_terms:
        #     opt_new_terms, opt_new_offset = self.reduce_deg(new_terms[new_term], new_term)
        #     for opt_new_term in opt_new_terms:
        #         if opt_new_term not in ans_terms:
        #             ans_terms[opt_new_term] = opt_new_terms[opt_new_term]
        #         else:
        #             ans_terms[opt_new_term] += opt_new_terms[opt_new_term]
        #     ans_offset += opt_new_offset
        return new_terms, new_offset

    def analyze_response(self, response, offset=0, input_dict=None, qubo_dict=None):
        """
        Analyze the record from the solver
        :param response: The record from the solver
        :param offset: The offset of the energy
        :param input_dict: The input for the solver
        :param bqm: The binary quadratic model for the solver
        :return: A tuple of two integers factorized from n
        """
        num_read = len(response.record)
        non_zero = len(qubo_dict)
        opt_energy = self.get_opt_energy(input_dict)
        sol_count = 0
        opt_count = 0
        for i in range(0, num_read):
            energy = response.record.energy[i]
            sample = response.record.sample[i]
            num_occurences = response.record.num_occurrences[i]
            if self.is_answer(sample, input_dict):
                sol_count += num_occurences
            if energy + offset == opt_energy:
                opt_count += num_occurences
        data_name = os.getenv("INPUT_FILE")
        output_dir = "output/" + data_name + "/"
        stat_dict = {
            "num_reads": num_read,
            "non_zero": non_zero,
            "sol_pct": int(sol_count) / 1000,
            "opt_pct": int(opt_count) / 1000
        }
        solver_config = os.getenv("SOLVER_CONFIG")
        with open(output_dir + solver_config + "_2.json", "w") as f:
            json.dump(stat_dict, f, indent=4)

    def is_answer(self, sample, input_dict=None) -> bool:
        """
        Check if the sample is the answer
        :param sample: The sample from the solver
        :param input_dict: The input for the solver
        :return: True if the sample is the answer
        """
        raise NotImplementedError

    def get_answer(self, sample, input_dict=None) -> (int, int):
        """
        Get the answer from the sample
        :param sample: The sample from the solver
        :param input_dict: The input for the solver
        :return: A tuple of two integers factorized from n
        """
        raise NotImplementedError

    def get_opt_energy(self, input_dict=None) -> int:
        """
        Get the optimal energy
        :param input_dict: The input for the solver
        :return: The optimal energy
        """
        raise NotImplementedError





