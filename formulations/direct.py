from collections import defaultdict
from math import log

from dimod import BinaryQuadraticModel

from .base import Formulation
from .solvers import solve_quantum_annealing


class DirectFormulation(Formulation):
    def __init__(self):
        super().__init__()

    def solve(self, n):
        l = int(log(n, 2) + 1) // 2
        # print(l)

        # Init base variables
        for i in range(0, l * 2):
            if i == 0 or i == 4:
                continue
            self.var_map[("x", i)] = self.var_cnt
            self.var_cnt += 1
        # print(self.var_map)

        # Init pq QUBO
        pq_terms = {}
        pq_offset = 0
        for i in range(0, l):
            for j in range(0, l):
                var_tuple = (0, 0)
                if i == 0:
                    if j == 0:
                        pq_offset += 1
                    else:
                        var_tuple = (self.var_map[("x", j + l)],)
                elif j == 0:
                    var_tuple = (self.var_map[("x", i)],)
                else:
                    var_tuple = (self.var_map[("x", i)], self.var_map[("x", j + l)])
                if var_tuple != (0, 0):
                    pq_terms[var_tuple] = pow(2, i + j)
        # print(pq_terms)
        # print(pq_offset)

        # Init (n - pq)^2 HUBO
        hubo_terms = {}
        for term in pq_terms.keys():
            if term not in hubo_terms:
                hubo_terms[term] = pq_terms[term]
            else:
                hubo_terms[term] += pq_terms[term]
        for term in pq_terms.keys():
            if term not in hubo_terms:
                hubo_terms[term] = -2 * (n - pq_offset) * pq_terms[term]
            else:
                hubo_terms[term] -= 2 * (n - pq_offset) * pq_terms[term]
        for term1 in pq_terms.keys():
            for term2 in pq_terms.keys():
                var_coef = 2 * pq_terms[term1] * pq_terms[term2]
                var_set = set(list(term1) + list(term2))
                var_tuple = tuple(var_set)
                hubo_terms[var_tuple] = var_coef
        hubo_offset = pow(n - pq_offset, 2)
        # print(hubo_terms)
        print(hubo_offset)

        # Optimize HUBO
        # print(self.hubo_optimize(1, (1, 2, 4, 5)))
        opt_hubo_terms = {}
        for term in hubo_terms:
            new_opt_terms = self.hubo_optimize(hubo_terms[term], term)
            for new_opt_term in new_opt_terms:
                if new_opt_term not in opt_hubo_terms:
                    opt_hubo_terms[new_opt_term] = new_opt_terms[new_opt_term]
                else:
                    opt_hubo_terms[new_opt_term] += new_opt_terms[new_opt_term]
        # # print(opt_hubo_terms)

        # Initialize QUBO
        q = defaultdict(int)
        for term in opt_hubo_terms:
            if len(term) == 2:
                q[term] = opt_hubo_terms[term]
            else:
                q[term[0], term[0]] = opt_hubo_terms[term]
        print(q)
        bqm = BinaryQuadraticModel.from_qubo(q)

        # Solve QUBO
        response = solve_quantum_annealing(bqm=bqm, method="DI", num_reads=1000)

        # Analyze result
        energy = response.first.energy
        sample = response.first.sample
        print(energy + hubo_offset)
        ans1 = 1 + 2 * sample[0] + 4 * sample[1] + 8 * sample[2]
        ans2 = 1 + 2 * sample[3] + 4 * sample[4] + 8 * sample[5]
        print(ans1, ans2)
        return (ans1, ans2)

