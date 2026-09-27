# Runtime database

The local-model database is a read-only lookup table for the narrow-Q topology
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

Application requests never modify this directory. If an input requires a model
that is not indexed, the API returns an unsupported-input error. New models are
constructed and checked outside the deployed project, then exported as a new
database version.
