# End-to-end guide for contributing `VoronoiS1` to PLUMED

This guide assumes that you already copied `VoronoiS1.cpp` from the [`VoronoiCVs` prototype](https://github.com/Zhang-pchao/VoronoiCVs/blob/devel/VoronoiCV_0202/VoronoiS1.cpp)
into `src/colvar/VoronoiS1.cpp` in your working tree.  It walks through the remaining
steps required to get the feature ready for an upstream pull request.

## 1. Prepare the code

1. Run `astyle --options=astyle/options src/colvar/VoronoiS1.cpp` to align with the project style.
2. Double-check the namespaces: `VoronoiS1` must live in `namespace PLMD { namespace colvar { ... } }`.
3. Ensure the constructor is registered with `PLUMED_REGISTER_ACTION(VoronoiS1, "VORONOIS1")` so the
   dispatcher exposes the CV.
4. Confirm that `registerKeywords` forwards to `Colvar::registerKeywords` and documents every input parameter.

## 2. Embed user documentation in the source

1. Add a `//+PLUMEDOC COLVAR VORONOIS1` block at the top of the file with:
   * A concise description of the collective variable.
   * A bullet list of keywords (e.g., `ATOMS`, `NN_CUTOFF`, etc.).
   * A `\verbatim ... \endverbatim` example showing a minimal PLUMED input fragment.
2. Close the block with `//+ENDPLUMEDOC`.
3. Rebuild the manual locally with `make -C user-doc doc`.  Fix any Doxygen warnings, especially missing keyword descriptions.

## 3. Extend the user manual

1. Create or edit the corresponding RST include under `user-doc/colvar`.  Follow the naming convention
   `colvar_voronois1.dox`.  You can copy `user-doc/colvar/colvar_path.dox` as a template.
2. Inside the file, reference the generated documentation with:
   ```
   \page COLVAR_VORONOIS1 VORONOIS1
   \include src/colvar/VoronoiS1.cpp
   ```
3. Append an entry in `user-doc/colvar/index.dox` so the new page appears in the menu.
4. Re-run `make -C user-doc doc` and verify the HTML preview under `user-doc/html/colvar_voronois1.html`.

## 4. Add regression tests

1. Create a new folder `regtest/rt-voronois1-basic`.
2. Populate it with:
   * `Makefile` (copy from another CV test such as `regtest/rt-path-basic/Makefile`).
   * `plumed.dat` (exercise the CV with a small system).
   * Reference outputs (`COLVAR`, `LOG`) generated with `make` inside the test directory.
3. Register the test by appending its relative path to `regtest/rtlist`.
4. Run the suite from the repository root with `make regtest-mpi` (or `make regtest`) and ensure the new test passes.

## 5. Final validation before submission

1. Run `ctest` to confirm the C++ unit tests are unaffected.
2. Execute `make docs` to rebuild the entire documentation set.
3. Commit with a message like `Add VoronoiS1 collective variable [makedoc]` so the documentation pipeline triggers.
4. Push your branch to `github.com/Zhang-pchao/plumed2` and open a pull request targeting `plumed/plumed2`'s `master` (or
   the branch requested by maintainers).
5. In the PR description, summarize:
   * The new CV (`VoronoiS1`) and its main capability.
   * The regression test(s) that accompany it.
   * The documentation updates under `user-doc`.
6. Respond promptly to reviewer feedback and keep your branch rebased on the latest upstream changes.

Following these steps ensures that the code, documentation, and verification artifacts for `VoronoiS1` meet the
PLUMED contribution standards.
