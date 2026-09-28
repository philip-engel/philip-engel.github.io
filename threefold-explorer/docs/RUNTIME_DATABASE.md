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

The computation path is read-only. It contains no geometric model builders.
Models are generated and audited offline; request-time work consists of lookup,
full integral boundary transport, and global MV and van Kampen calculations.
Integral clutching is restored after reduction into the torsion fundamental domain.

All allowed component classes are present for Mumford fillings over I0 through
I9 with at most 12 components and absolute linearization order at most 12.
For I_n* the bound applies to the upstairs I_(2n) wheel, so 1 <= n <= 6.
The requested log vectors remain exactly invariant, with m dividing d.

Run `sage -python tests/smoke_test.py` from the API or package directory after
unpacking the archive. Tests check table coverage and checksums, the two earlier
sphere models, all four semistable sphere presets and their untwisted controls,
mixed local narrowness, signed Mumford weights, starred twists and input errors.
