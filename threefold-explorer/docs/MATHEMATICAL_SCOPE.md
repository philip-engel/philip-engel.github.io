# Mathematical scope

The complete entry point is `explorer_api.compute`, which calls
`parameterized_models.explore_narrow_q` and then the exact database
Mayer–Vietoris and van Kampen implementations.

## Implemented hypotheses

- Q is globally narrow; P may be non-narrow.
- At every singular fiber, the section pair satisfies the local component
  condition inherited from the monodromy construction.
- The linearization divisor is supported on semistable I_k fibers, including
  added smooth I0 fibers. Its nonzero weights have one sign and sum to
  the height pairing of P and Q.
- Log vectors are exact rational vectors fixed by the corresponding 4 by 4
  monodromy.
- Their torsion order divides the minimal semistable-reduction degree.
- Mumford compactifications use the prescribed A2 tiling or rank-one wheel.
- Mumford models have at most 12 components and linearization order at most 12.
- Quotient compactifications use the selected minimal resolution.

## `None` and explicit zero

| Fiber | `None` | Explicit zero vector |
|---|---|---|
| I_k, including I0 | Original or prescribed Mumford filling | The same filling |
| Potentially good additive fiber | Original O(P-O) filling | Good-reduction quotient of O(P'-O') |
| Positive-index I_n* | Original O(P-O) filling | Quadratic semistable-reduction quotient of O(P'-O') |

For a quotient, zero means zero *added* twist. The canonical lifted divisor
character is still part of the construction. Integer parts of log parameters
are retained in the global clutching.

## Reported S6 result

The code computes topology from a supplied collection of marked geometric local
models. When it reports `S6_for_supplied_smooth_model = True`, the integral
cohomology is that of S6 and the computed fundamental group is trivial for that
smooth model.

## Outside the current entry point

- non-narrow Q;
- torsion points fixed modulo the lattice but not by the chosen rational lift;
- log order not dividing the semistable-reduction degree;
- linearization zeros or poles on additive fibers;
- Mumford subdivisions other than the prescribed A2 model, or models beyond
  the displayed lookup bounds;
- a general recognition theorem assigning a familiar name to every finitely
  presented fundamental group.

Detailed derivations and full validation certificates are maintained separately
and are not part of the deployment package.
