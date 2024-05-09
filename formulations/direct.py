from collections import defaultdict
from math import log

from dimod import BinaryQuadraticModel

from .base import Formulation
from .solvers import solve_quantum_annealing, solve_simulated_annealing


class DirectFormulation(Formulation):
    def __init__(self):
        super().__init__()

    def solve(self, n):
        l = int(log(n, 2)) // 2 + 1
        # print(l)

        # Init base variables
        for i in range(0, l * 2):
            if i == 0 or i == l:
                continue
            self.var_map[("x", i)] = self.var_cnt
            self.var_cnt += 1
        # print(self.var_map)

        # Init pq HUBO
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
        print(pq_terms)
        print(pq_offset)

        # Init (n - pq)^2 HUBO
        hubo_terms = {}
        pq_terms[()] = pq_offset - n
        for term1 in pq_terms.keys():
            for term2 in pq_terms.keys():
                var_coef = pq_terms[term1] * pq_terms[term2]
                var_set = set(list(term1) + list(term2))
                var_tuple = tuple(var_set)
                if var_tuple not in hubo_terms:
                    hubo_terms[var_tuple] = var_coef
                else:
                    hubo_terms[var_tuple] += var_coef
        hubo_offset = hubo_terms[()]
        hubo_terms.pop(())
        print(hubo_terms)
        print(hubo_offset)

        # Optimize HUBO
        # print(self.hubo_optimize(1, (1, 2, 4, 5)))
        opt_hubo_terms = {}
        for term in hubo_terms:
            new_opt_terms, new_opt_offset = self.hubo_to_qubo(hubo_terms[term], term)
            for new_opt_term in new_opt_terms:
                if new_opt_term not in opt_hubo_terms:
                    opt_hubo_terms[new_opt_term] = new_opt_terms[new_opt_term]
                else:
                    opt_hubo_terms[new_opt_term] += new_opt_terms[new_opt_term]
            hubo_offset += new_opt_offset
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

        # # 143
        # bqm.fix_variable(self.var_map[("x", 1)], 0)
        # bqm.fix_variable(self.var_map[("x", 2)], 1)
        # bqm.fix_variable(self.var_map[("x", 3)], 1)
        # bqm.fix_variable(self.var_map[("x", 5)], 1)
        # bqm.fix_variable(self.var_map[("x", 6)], 0)
        # bqm.fix_variable(self.var_map[("x", 7)], 1)

        ## 1961
        # bqm.fix_variable(self.var_map[("x", 1)], 0)
        # bqm.fix_variable(self.var_map[("x", 2)], 1)
        # bqm.fix_variable(self.var_map[("x", 3)], 0)
        # bqm.fix_variable(self.var_map[("x", 4)], 0)
        # bqm.fix_variable(self.var_map[("x", 5)], 1)
        # bqm.fix_variable(self.var_map[("x", 7)], 0)
        # bqm.fix_variable(self.var_map[("x", 8)], 1)
        # bqm.fix_variable(self.var_map[("x", 9)], 0)
        # bqm.fix_variable(self.var_map[("x", 10)], 1)
        # bqm.fix_variable(self.var_map[("x", 11)], 1)

        # Solve QUBO
        # response = solve_simulated_annealing(bqm=bqm, method="DI", num_reads=1000)
        response = solve_quantum_annealing(bqm=bqm, method="DI", num_reads=1000)

        # Analyze result
        energy = response.first.energy
        sample = response.first.sample
        print(energy + hubo_offset)
        ans1 = 1
        ans2 = 1
        for i in range(1, l):
            ans1 *= 2
            if sample[self.var_map[("x", i)]] == 1:
                ans1 += 1
        for i in range(l + 1, l * 2):
            ans2 *= 2
            if sample[self.var_map[("x", i)]] == 1:
                ans2 += 1
        print(ans1, ans2)
        return ans1, ans2

