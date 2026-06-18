# ======================================================================
# apply_fesom_compile_flags(<target>)
#
# Per-compiler / per-precision Fortran flags for FESOM3.
#
# The DOUBLE-PRECISION + Intel path is the **anchor build** for bit-identity.
# Its flags are transcribed *verbatim* from the byte-gate oracle
#   FESOM2 v2.7.3  ->  src/CMakeLists.txt:335-404
# NOT from the tracer_dwarf (the dwarf is a structural reference only and its
# flags differ: it adds -no-prec-sqrt/-ip to ifx and -march=core-avx2 on
# Levante, neither of which FESOM2 v2.7.3 uses on Levante — see docs/LESSONS.md).
#
# Faithful anchor (Intel classic, Levante):
#   -O3 -r8 -i4 -fp-model precise -no-prec-div -fimf-use-svml
#   -init=zero -no-wrap-margin -fpe0 -fpp           (base, both ifort/ifx)
#   -no-prec-sqrt -ip                               (ifort classic only)
#   <no -march>                                     (Levante: FESOM2's -march line
#                                                    is COMMENTED OUT -> SSE2 baseline,
#                                                    no FMA contraction. Matching it is
#                                                    required for byte-identity.)
#
# Precision is selected two ways at once (mirrors FESOM2/dwarf):
#   - default-real promotion flag (-r8 / -fdefault-real-8 ...) AND
#   - the USE_SINGLE_PRECISION / USE_HALF_PRECISION cpp define that sets WP.
# Both move together so unsuffixed literals and `real(kind=WP)` stay consistent.
# ======================================================================

# Resolve the platform strategy (set by env.sh; drives platform flag branches).
if(DEFINED ENV{FESOM_PLATFORM_STRATEGY})
    set(FESOM_PLATFORM_STRATEGY $ENV{FESOM_PLATFORM_STRATEGY})
else()
    set(FESOM_PLATFORM_STRATEGY "unknown")
endif()
message(STATUS "FESOM_PLATFORM_STRATEGY = ${FESOM_PLATFORM_STRATEGY}")

function(apply_fesom_compile_flags target_name)

    # --- precision cpp define (sets WP in mod_precision) ---
    if(USE_HALF_PRECISION)
        target_compile_definitions(${target_name} PRIVATE USE_HALF_PRECISION)
    elseif(USE_SINGLE_PRECISION)
        target_compile_definitions(${target_name} PRIVATE USE_SINGLE_PRECISION)
    endif()

    # --- default-real promotion flag, per compiler family ---
    if(USE_HALF_PRECISION)
        set(_intel_real_flag "-r2")
        set(_gnu_real_flags   "")          # no native real(2); HP is NVHPC-only
    elseif(USE_SINGLE_PRECISION)
        set(_intel_real_flag "-r4")
        set(_gnu_real_flags   "")          # WP=4 via cpp; no default-real promotion
    else()
        set(_intel_real_flag "-r8")
        set(_gnu_real_flags   "-fdefault-real-8;-fdefault-double-8")
    endif()

    # ==================================================================
    # Intel classic (ifort) and Intel LLVM (ifx)
    # ==================================================================
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR
       CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")

        if(USE_HALF_PRECISION)
            message(FATAL_ERROR "Intel Fortran has no REAL(kind=2); use NVHPC for HP.")
        endif()

        # FESOM2 v2.7.3 base flags (src/CMakeLists.txt:335). -O3 for the anchor;
        # Debug swaps to -O0 and adds run-time checks (diagnostics only; with
        # -fp-model precise these do not change FP results vs the anchor).
        target_compile_options(${target_name} PRIVATE
            $<$<NOT:$<CONFIG:Debug>>:-O3>
            $<$<CONFIG:Debug>:-O0;-g;-traceback;-check;all,noarg_temp_created,bounds,uninit>
            ${_intel_real_flag} -i4 -fp-model precise -no-prec-div -fimf-use-svml
            -init=zero -no-wrap-margin -fpe0 -fpp)

        # ifort-classic-only flags (FESOM2 v2.7.3 src/CMakeLists.txt:338-340)
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
            target_compile_options(${target_name} PRIVATE -no-prec-sqrt -ip)
        endif()

        # Platform branch — must mirror FESOM2 v2.7.3 (src/CMakeLists.txt:348-383).
        # On Levante the -march override is COMMENTED OUT in FESOM2 -> no flag here.
        if(FESOM_PLATFORM_STRATEGY STREQUAL "levante.dkrz.de")
            # (intentionally empty — matches FESOM2 v2.7.3 commented-out -march)
        else()
            target_compile_options(${target_name} PRIVATE -xHost)
        endif()

    # ==================================================================
    # GNU (gfortran)
    # ==================================================================
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

        if(USE_HALF_PRECISION)
            message(FATAL_ERROR "GNU Fortran has no native REAL(2)/FP16.")
        endif()

        # FESOM2 v2.7.3 GNU base flags (src/CMakeLists.txt:396).
        target_compile_options(${target_name} PRIVATE
            $<$<NOT:$<CONFIG:Debug>>:-O3>
            $<$<CONFIG:Debug>:-O0;-g;-fbacktrace;-fcheck=all>
            -ffloat-store -finit-local-zero -finline-functions -fimplicit-none
            ${_gnu_real_flags} -ffree-line-length-none -cpp)

        # gfortran >= 10 is strict about API rank mismatches (MPI calls).
        if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
            target_compile_options(${target_name} PRIVATE -fallow-argument-mismatch)
        else()
            target_compile_options(${target_name} PRIVATE -Wno-argument-mismatch)
        endif()

        # Platform branch (FESOM2 v2.7.3 src/CMakeLists.txt:416-434).
        # NOTE: FESOM2's GNU Levante line adds -flto, but FESOM2 links a single
        # executable whereas FESOM3 links lib+exe; -flto across our static-lib
        # boundary needs gcc-ar/gcc-ranlib and breaks plainly on gfortran 8.5.
        # GNU is the portability build (Intel is the bit-identity anchor) and
        # -flto can itself reduce FP reproducibility, so we omit it. (L3)
        if(FESOM_PLATFORM_STRATEGY STREQUAL "levante.dkrz.de")
            target_compile_options(${target_name} PRIVATE
                -march=znver3 -mtune=znver3 -ftree-vectorize)
        endif()

    else()
        message(WARNING "Unhandled Fortran compiler '${CMAKE_Fortran_COMPILER_ID}'"
                        " — no FESOM flags applied (anchor build is Intel/GNU only).")
    endif()
endfunction()
