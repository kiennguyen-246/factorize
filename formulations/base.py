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
    def hubo_to_qubo(self, coef=1, term=()):
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
        prefix = ()
        x = term[:-2]
        # prefix = term[:-3]
        # x = term[-3:-2]
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
        # if len(term) > 3:
        #     new_terms = {prefix + new_term: new_terms[new_term] for new_term in new_terms}
        #     new_terms[prefix] = new_offset
        #     new_offset = 0
        ans_terms = {}
        ans_offset = new_offset
        for new_term in new_terms:
            opt_new_terms, opt_new_offset = self.hubo_to_qubo(new_terms[new_term], new_term)
            for opt_new_term in opt_new_terms:
                if opt_new_term not in ans_terms:
                    ans_terms[opt_new_term] = opt_new_terms[opt_new_term]
                else:
                    ans_terms[opt_new_term] += opt_new_terms[opt_new_term]
            ans_offset += opt_new_offset
        return ans_terms, ans_offset


