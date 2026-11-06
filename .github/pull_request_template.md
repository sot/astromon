## Description

<!--If this fixes an issue then fill this, otherwise DELETE the line below -->
Fixes #

## Interface impacts
<!-- API changes, file format updates, coordination of changes with the community. -->

## Result versions
<!-- Tasks (@task) and cached functions (@stored_result) carry a version; results made by
another version are recomputed. DELETE the line that does not apply. -->
- [ ] Changes what a task or cached function produces, and its `version` is bumped
  (`python -m astromon.scripts.maintenance.report_stale_results --root <work> --root <archive> --bump <name>`
  shows what that will rerun)
- [ ] Does not change what any task or cached function produces

## Testing
<!-- If relevant describe any special setup for testing. -->

### Unit tests
<!-- At least one of these must be checked if unit tests exist. DELETE the unchecked/untested options. -->
- [ ] No unit tests
- [ ] Mac
- [ ] Linux
- [ ] Windows

Independent check of unit tests by [REVIEWER NAME]
- [ ] [PLATFORM]:

### Functional tests
<!-- Describe and document results of any functional tests, otherwise leave the text below -->
No functional testing.
