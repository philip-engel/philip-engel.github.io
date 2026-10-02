# Mathematical scope

`explorer_api.compute` uses `fiberwise_narrow.explore_fiberwise_narrow` and the
integral database Mayer–Vietoris and van Kampen implementations.

## Section condition

The line bundle is O(P−O); Q is the translating section. At each singular
fiber, P or Q must be narrow. The choice may vary from fiber to fiber; neither
section needs to be globally narrow. For a component group with several
congruences, every congruence must hold for P or every congruence must hold for Q.
The API reports the specific fibers where neither section is narrow.

## Other hypotheses and bounds

- Linearization zeros and poles lie on I_k fibers, including added smooth I0.
  All nonzero weights have one sign and sum to the height pairing of P and Q.
- Log vectors are exact rational vectors fixed by the 4×4 monodromy, with
  torsion order dividing the minimal semistable-reduction degree, except at
  smooth I0 fibers of linearization weight zero, where every order is allowed.
- Mumford compactifications use the prescribed A2 tiling or rank-one wheel.
  Every component class satisfying the section condition is tabulated for
  I0 through I9, with max(1,k) times the absolute weight at most 12.
- Quadratic I_n* quotients are tabulated for 1 <= n <= 6: at most 12 components
  in the upstairs I_(2n) model. This includes every positive-index starred
  fiber occurring on a RES. Original fillings are also present.
- Quotient compactifications use the selected minimal resolution.

## None, zero, and integral clutching

`None` selects the original O(P−O) filling, or the prescribed Mumford filling.
At an additive fiber an explicit zero selects the reduction quotient of
O(P′−O′) with zero added twist. If P is locally narrow, these agree under the
implemented marked comparison. For non-narrow P they can differ.
At I_k, including I0, None and zero select the same filling.

Integer parts of log vectors are retained. A primitive integer vector at a
smooth fiber has torsion order one and fiber multiplicity one, while its
clutching can change the topology. The four semistable sphere presets use
this modification in the quotient direction.

## Fractional logs at smooth fibers

At an I0 slot with linearization weight zero, enter any exact rational vector
theta. Its least common denominator m is the fiber multiplicity; the reduction
is the isogenous complex torus with lattice Lambda' = Lambda + Z theta.
The filling retracts onto this torus. Its cochain ranks are (1,4,6,4,1), with
zero differentials, independently of m. No denominator-specific model is stored.

If L is a basis matrix for Lambda', the boundary homomorphism on fundamental
groups is (u,k) -> L^-1 (u + k theta). Its primitive kernel is (-m theta,m).
An integral normal-circle completion gives the boundary comparison; exterior
powers give all cochain maps. In the base-first convention the filling map
is [wedge^q(L^-t); i_theta wedge^q(L^-t)]. Van Kampen uses the relation
meridian^m = fiber^(m theta). The full integer lift of theta affects these maps.

This extension does not allow fractional twists at an I0 slot with nonzero
linearization weight: that slot has a Mumford degeneration. Add a separate
weight-zero I0 slot instead. The API schema reports `arbitrary_denominators`
and sets `allowable_denominators` to null at the new smooth log sites.

## Results

Every integral cohomology group H0 through H6 is computed from the MV cone.
Endpoint, Poincaré duality, torsion duality and van Kampen/UCT consistency are
checked afterward. Fundamental-group recognition is bounded; an unrecognized
presentation is not replaced by its abelianization.

“Topology of S6” means trivial computed fundamental group and integral
cohomology of S6 for the supplied smooth geometric model. The analytic
Mumford compactification retains the assumptions of the underlying construction.

## Outside this implementation

- Both P and Q non-narrow at the same singular fiber.
- Torsion twists fixed only modulo the lattice.
- Log order not dividing the semistable-reduction degree at singular or
  Mumford fibers (the weight-zero smooth exception is described above).
- Linearization zeros or poles on additive fibers.
- Other Mumford subdivisions or models beyond the stated bounds.
- General recognition of every finitely presented fundamental group.

Detailed derivations and validation certificates remain in the research workspace.
