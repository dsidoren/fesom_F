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
add_fesom_test(test_bottom 1)
add_fesom_test(test_io_decomp 1)
add_fesom_test(test_io_decomp 2)
add_fesom_test(test_io_decomp 8)
add_fesom_test(test_io_decomp_gather 1)
add_fesom_test(test_io_decomp_gather 2)
add_fesom_test(test_io_posix 1)
add_fesom_test(test_vector_rotate 1)
add_fesom_test(test_io_means 1)

# dump_diff.py self-check (pure python; no external data needed).
add_test(NAME dump_diff_selftest
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tools/dump_diff.py --selftest)

# End-to-end analytic driver (built from src/drivers/). Runs clean on 1 + 2 ranks.
add_test(NAME fesom_analytic_np1
         COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 1 $<TARGET_FILE:fesom_analytic>)
add_test(NAME fesom_analytic_np2
         COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np 2 $<TARGET_FILE:fesom_analytic>)

# Conservation + no-leakage gate for the bottom-at-vertices change (docs/plans/
# 20260910-fesom3-bottom-at-vertices.md). zstar conserves to round-off, so this is a hard
# gate; the driver itself error-stops when FESOM3_CONSERVE_TOL is exceeded. Needs the pi
# mesh and the pool IC file, so it is skipped (not failed) when those are absent.
foreach(_np 1 2)
    add_test(NAME fesom_conserve_zstar_np${_np}
             COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np ${_np}
                     $<TARGET_FILE:fesom_conserve>)
    set_tests_properties(fesom_conserve_zstar_np${_np} PROPERTIES
        ENVIRONMENT "FESOM3_MESH_DIR=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi;FESOM3_IC_FILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc;FESOM3_NSTEPS=20;FESOM3_WHICH_ALE=zstar;FESOM3_CONSERVE_TOL=1e-12"
        SKIP_REGULAR_EXPRESSION "cannot open;No such file")
endforeach()
