/opt/homebrew/Library/Homebrew/cmd/shellenv.sh: line 18: /bin/ps: Operation not permitted
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

The deployment is entirely read-only. The prescribed A2 Mumford fillings are
tabulated for every allowed input with at most 12 components and linearization
order at most 12. The quadratic I_n* quotient table is complete in the stated
exactly-invariant regime. Request-time work consists only of lookup, integral
boundary re-marking, and the global Mayer--Vietoris and van Kampen calculations.
