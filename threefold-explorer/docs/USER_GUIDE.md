# User guide

The interface is staged. Each stage has enough information
to validate the next one, so impossible data can be rejected near the field that
caused the problem.

## 1. Choose the rational elliptic surface

Enter an OS table number from 1 through 74. The explorer displays the ordered
Kodaira fibers, the MW group, its chosen coordinate length, its height matrix,
and the available collision profiles.

The collision selector lists complete fiber configurations. The 289 OS/profile
choices cover all 279 Persson–Miranda configurations. Selecting one revises the
ordered fibers and local section conditions. See [Collision profiles](COLLISION_PROFILES.md)
for marking conventions and the source checks.

Backend call: `describe_os_entry(os_entry, profile)`.

## 2. Enter P and Q

P and Q are integer tuples in the displayed MW basis; torsion entries come last
and are automatically reduced modulo their orders. Enter each tuple as comma-separated integers. At each fiber, all displayed congruences must hold for P or all must hold for Q. The choice may vary from fiber to fiber.

After P and Q are entered, the explorer shows their height pairing and verifies
that at least one section is narrow at each fiber, and reports which one.

Backend call: `section_and_linearization_schema(...)`.

## 3. Choose the linearization divisor

Display one horizontal slot for every ordered Kodaira fiber. Slots at additive
fibers are disabled because the current scope allows nonzero weights only at
I_k. An “add smooth fiber” button appends an I0 slot. The running sum should be
shown next to the required total degree `<P,Q>`.

The interface checks that every nonzero weight has the same sign and that their
sum equals the required degree.

## 4. Enter log modifications

Once the weights are fixed, the explorer computes each invariant lattice and
displays one log slot under the corresponding fiber. Each slot explains:

- the number of invariant coordinates;
- the integral ambient basis columns;
- the allowed denominators: arbitrary positive orders at weight-zero I0 slots,
  and divisors of the reduction order elsewhere;
- what `None` means;
- what an explicit zero means.

Exact input such as `1/2` should be accepted; decimal input should be rejected.
Additional smooth slots remain available. For advanced use, the interface may
toggle between invariant and four-coordinate ambient input.

Backend call: `log_transform_schema(...)`.

## 5. Compute and display

The result panel should lead with:

- the fundamental group description and abelianization;
- H^0 through H^6 over the integers;
- Euler characteristic;
- the integral-homology-sphere flag;
- the S6 flag for the supplied smooth geometric model.

The expandable local-model panel describes the selected fillings. See Mathematical Scope for the construction conventions.

Backend call: `compute(payload)`.

## Error presentation

Errors should be associated with the current stage: unknown collision profile,
wrong section-tuple length, both sections non-narrow at a fiber, wrong divisor sum, non-invariant log
vector, unsupported torsion order, or an uncached model exceeding the configured
component guard. Preserve the engine's mathematical explanation in a details
panel rather than replacing it with a generic “invalid input” message.
