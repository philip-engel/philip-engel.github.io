# Fiber collision profiles

The selector includes **279 distinct Kodaira configurations**, represented by
**289 OS/profile choices**: 74 defaults and 215 alternatives. There are multiple
OS entries for some fiber configurations because their MW lattices differ.
Existing profile identifiers and example markings are preserved.

The new profiles come from genus-zero J-map covers and quadratic twists, using
[Miranda, §3](https://www.math.colostate.edu/~miranda/preprints/Perssonslist.pdf).
The classification check compares the full set with the 379 candidate
configurations minus the 100 exclusions in Miranda's Table 2.1.
[Fukae, Table 3](https://arxiv.org/pdf/math/0205062) supplies the 74 defaults.

Choose `profile="6II"` at OS1 for six cuspidal fibers. OS1 also allows one through
four II fibers and the remaining I1 fibers; five II plus two I1 is impossible.
At OS13, four III fibers are allowed, but three III plus I2 plus I1 is excluded;
that latter profile belongs to OS14. Four III forces two-torsion and cannot
belong to torsion-free OS14.

For new profiles, an explicit integral lattice isometry identifies the free
height matrix with the default OS matrix. This is a chosen marking, not a
parallel-transport computation through a specified family. Recheck the ordered
fibers and displayed local narrowness conditions when changing profile.

All matrices and marked section cocycles are precomputed in
`engine/collision_profiles.json` (about 76 KiB). Neither cover search nor lattice
isometry search runs in the applet. The local-filling database and its bounds
are unchanged. Construction certificates are retained with the research code.

Shioda's [1992 correction](https://doi.org/10.3792/pjaa.68.251), Remark (ii),
corrects entries 32 and 70; it does not add further OS entries. The engine
already uses the corrected height lattice for 32 and torsion Z/4 for 70.
