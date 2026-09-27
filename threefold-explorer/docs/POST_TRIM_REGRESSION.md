# Post-trim regression

Validated on 2026-09-27 with SageMath 10.9 after exporting the compact,
read-only database.

## Historical global examples

All 45 end-to-end examples saved by the private pre-trim validation suite were
replayed through the deployment package. For every example, the newly computed
integral cohomology groups and fundamental-group description agreed exactly with
the saved result. There were no exceptions or mismatches.

The inputs remain in the private validation corpus and are not copied into this
distribution.

## Database-wide cochain check

All 487 retained records passed the following checks:

- content-address checksum verification;
- the local filling-to-boundary restriction is a cochain map;
- where present, the saved local-to-standard boundary map is an integral
  quasi-isomorphism.

The checked inventory was 85 original divisor-bundle fillings, 275 finite
good-reduction quotients, 96 positive-index starred quotients, 30 Mumford
models, and one smooth product.
