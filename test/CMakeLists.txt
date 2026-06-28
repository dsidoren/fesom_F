# FESOM3 self-tests. add_fesom_test(<name> <nranks>) is defined in the top-level
# CMakeLists; each builds test/<name>.F90 against libfesom3 and registers a ctest.

add_fesom_test(test_params 1)
add_fesom_test(test_types  1)
add_fesom_test(test_partit 1)
add_fesom_test(test_partit 2)
add_fesom_test(test_partit 8)
add_fesom_test(test_halo   1)
add_fesom_test(test_halo   2)
add_fesom_test(test_halo   8)
add_fesom_test(test_dump   1)
add_fesom_test(test_mesh   1)
add_fesom_test(test_io_decomp 1)
add_fesom_test(test_io_decomp 2)
add_fesom_test(test_io_decomp 8)
add_fesom_test(test_vector_rotate 1)

# dump_diff.py self-check (pure python; no external data needed).
add_test(NAME dump_diff_selftest
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tools/dump_diff.py --selftest)

# End-to-end analytic driver (built from src/drivers/). Runs clean on 1 + 2 ranks.
add_test(NAME fesom_analytic_np1
         COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 1 $<TARGET_FILE:fesom_analytic>)
add_test(NAME fesom_analytic_np2
         COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 2 $<TARGET_FILE:fesom_analytic>)
