class Formulation:
    def __init__(self) -> None:
        self.var_map = {}
        self.var_cnt = 0

    def solve(self, n) -> (int, int):
        """
        The base solver for any formulations
        :param n: A safe semiprime
        :return: A tuple of two integers factorized from n
        """
        raise NotImplementedError

    # Helper functions
    def hubo_optimize(self, coef=1, term=()):
        """
        Recursively reduce a HUBO term until it becomes QUBO terms.
        The formulation for optimization is: x1x2x3 = x3s + 2x1x2 - 4x1s - 4x2s + 6s
        :param term: The coefficient and a tuple of variables
        :return: A list of optimized terms
        """
        # print(coef, term)
        if len(term) <= 2:
            return {term: coef}
        self.var_map[("s", self.var_cnt)] = self.var_cnt
        self.var_cnt += 1
        new_terms = {}
        tuple1 = term[:-2]
        tuple2 = term[-2:-1]
        tuple3 = term[-1:]
        slack = self.var_map[("s", self.var_cnt - 1)]
        new_terms[tuple3 + (slack,)] = coef
        new_terms[tuple1 + tuple2] = 2 * coef
        new_terms[tuple1 + (slack,)] = -4 * coef
        new_terms[tuple2 + (slack,)] = -4 * coef
        new_terms[(slack,)] = 6 * coef
        ans_terms = {}
        for new_term in new_terms:
            opt_new_terms = self.hubo_optimize(new_terms[new_term], new_term)
            for opt_new_term in opt_new_terms:
                if opt_new_term not in ans_terms:
                    ans_terms[opt_new_term] = opt_new_terms[opt_new_term]
                else:
                    ans_terms[opt_new_term] += opt_new_terms[opt_new_term]
        return ans_terms


