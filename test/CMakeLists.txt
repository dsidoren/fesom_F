# FESOM3 self-tests. add_fesom_test(<name> <nranks>) is defined in the top-level
# CMakeLists; each builds test/<name>.F90 against libfesom3 and registers a ctest.

add_fesom_test(test_params 1)
add_fesom_test(test_types  1)
