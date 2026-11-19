Purpose: order cells along an inferred progression and summarize it as a graph.
Principle: PAGA abstracts the neighbour graph into cluster-level connectivity;
diffusion pseudotime (DPT) ranks cells by diffusion distance from a root cell,
so pseudotime is a relative rank, not a clock. Interpretation: choose the root
from prior biology (stem/progenitor markers), read the PAGA graph as hypotheses
about transitions, and treat inf pseudotime cells as disconnected components,
not as "late" cells. Pitfalls: parameters can manufacture the trajectory you
expected — validate against known biology and marker gradients; pseudotime is
not comparable across separate paths or datasets.
