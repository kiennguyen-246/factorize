import os

from formulations import DirectFormulation

formulation = DirectFormulation()
with open("dataset/.inp", "r") as f:
    n = int(f.readline().strip())
    os.environ["INPUT_FILE"] = ".inp"
    ans = formulation.solve(n)
    print(ans)