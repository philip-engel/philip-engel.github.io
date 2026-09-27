# Runtime database

The local-model database is the fast lookup layer for the narrow-Q topology
pipeline. Each content-addressed record contains only the information used by
the application:

- the model family and normalized parameters;
- the marked local filling and boundary cochain complexes;
- the local-to-standard boundary comparison used by Mayer–Vietoris;
- van Kampen peripheral relations and multiplicity;
- compact integral stalk summaries for verbose inspection;
- the monodromy, marking, and short geometric label needed for selection.

The index maps a normalized family/parameter key to an object checksum. Every
object is checksum-verified when read. No executable object or pickle is stored.

Expanded cell labels, raw resolutions, Smith-reduction certificates, comparison
homotopies, inverse certificates, and provenance hashes are deliberately absent.
They are required to derive or independently audit an entry, but are not inputs
to the global topology calculation.

The deployment also contains the checked parameterized constructors for the
prescribed A2 Mumford filling and quadratic I_n* quotient. If one of these
models is absent, the service constructs it, removes its derivation witnesses,
and caches the compact record for the life of the container. Mumford inputs are
bounded by 64 components and linearization order 12. Other local families stay
read-only and are selected from the finite table.
