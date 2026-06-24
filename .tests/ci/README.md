# Introscape CI fixtures

The fixture files in `data/` are copied from the upstream `xin-huang/selscape`
test fixtures at commit `98e19ca636d32bbe4c40e16c54be1626d0626dc6`.

These files are used only as small integration-test inputs for the top-level
`introscape` workflow. The vendored `selscape/` directory in this repository is
not the source of truth for these fixtures.

The CI demography model in `config/demography/` is copied from the upstream
`xin-huang/sstar` example data. It is used only for the `sstar2` simulation
training step; real-data sample groups still come from the dataset metadata.
