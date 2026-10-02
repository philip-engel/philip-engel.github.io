# Runtime database

The application uses 865 compact, content-addressed local models:

| Family | Records |
| --- | ---: |
| Original fillings | 85 |
| Finite good-reduction quotients | 275 |
| Quadratic I_n* quotients with narrow Q | 96 |
| Quadratic I_n* quotients with narrow P and non-narrow Q | 72 |
| Bounded Mumford fillings, including smooth wheels | 336 |
| Smooth product | 1 |

Each entry retains only the marked local cochain pair, boundary comparison,
peripheral relations, stalk summaries, marking and lookup parameters. Object
checksums are verified on reading. Derivations, expanded cellular models,
inverse certificates, historical records and private notes are excluded.

Local cochains are reduced by checked integral contractions, with their
restriction and normalization maps transported at the same time. The 865
models share 450 tables of exact finite-jet attachment coefficients in
`formulas/`; model records refer to these by checksum. The canonical marked
boundary complexes are unchanged. This keeps torsion and the full integral
clutching data while avoiding expanded translated-cell calculations for large
markings. Small comparisons use the original routine when it is faster.

The computation path is read-only. It contains no geometric model builders.
Models are generated and audited offline; request-time work consists of lookup,
full integral boundary transport, and global MV and van Kampen calculations.
Integral clutching is restored after reduction into the torsion fundamental domain.
Formula tables are decoded once per worker and cached. No formula compilation
or geometric construction takes place during a request.

All allowed component classes are present for Mumford fillings over I0 through
I9 with at most 12 components and absolute linearization order at most 12.
For I_n* the bound applies to the upstairs I_(2n) wheel, so 1 <= n <= 6.
The requested log vectors remain exactly invariant, with m dividing d except
at weight-zero I0 slots. Smooth logarithmic fillings of arbitrary order reuse
the single smooth product record: a 5x5 integral boundary marking and its
exterior powers encode the isogenous central torus and full rational clutch.
This adds no objects or formula tables to the database.

Run `sage -python tests/smoke_test.py` from the API or package directory after
unpacking the archive. Tests check table coverage and checksums, the two earlier
sphere models, all four semistable sphere presets and their untwisted controls,
mixed local narrowness, signed Mumford weights, starred twists and input errors.
