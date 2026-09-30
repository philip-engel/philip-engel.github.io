# Collision catalogue validation — 29 September 2026

* The full configuration set equals Miranda's 379 numerical candidates minus
  the 100 exclusions in Table 2.1: **279 configurations**, not merely a count match.
* **289 marked OS/profile choices** have the expected rank, torsion, and height
  matrix, with total elliptic monodromy equal to the identity and Euler sum 12.
* The geometric J-map covering certificates were replayed, including genus-zero
  and connectedness checks, boundary extraction, and quadratic-twist signs.
* All **87 pre-existing marked models** retain exactly their matrices, section
  cocycles, height matrices, and component maps.
* **289 complete topology computations** passed with representative sections,
  original fillings, and a smooth linearization slot where needed.
* **202 further complete computations** passed on all new profiles with
  maximal-denominator twists at their additive fibers.
* The existing smoke suite and **all 99 rational homology sphere presets** passed
  with their previous integral cohomology and fundamental-group descriptions.
* The isolated deployment archive and persistent worker passed, including 6II,
  existing S6 examples, repeated requests, and invalid-input handling.
* Runtime database writes were forbidden in both catalogue-wide test runs.
* The embedded research-notebook definitions loaded the companion catalogue and
  computed the 6II profile successfully.

These tests check representative inputs and the specified models; they do not
exhaust all P, Q, and twist vectors or identify connected marked moduli spaces.
The deployed lookup adds 77,536 bytes. The 865 local-filling records are unchanged.
Full certificates and per-case outputs are retained in the research workspace.
