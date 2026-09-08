###########################################################
# WANNIER90
###########################################################
add_library(qe_wannier90 INTERFACE)
qe_install_targets(qe_wannier90)
if(QE_WANNIER90_INTERNAL)
    message(STATUS "Installing Wannier90 via submodule")
    qe_git_submodule_update(external/wannier90)

    set(WANNIER90_SHARED_LIBS ${BUILD_SHARED_LIBS})
    # Wannier90's own install rules must run so that Wannier90_lib/Wannier90_post
    # end up in an export set; otherwise qe_wannier90 (which INTERFACE-links them)
    # cannot be exported as part of qeTargets.
    set(WANNIER90_INSTALL ON)
    set(WANNIER90_TEST OFF)
    add_subdirectory(wannier90)

    target_link_libraries(qe_wannier90 INTERFACE Wannier90::wannier90)

    ###########################################################
    # w90chk2chk.x
    ###########################################################
    add_executable(qe_w90chk2chk_exe wannier90/src/w90chk2chk.F90)
    set_target_properties(qe_w90chk2chk_exe PROPERTIES OUTPUT_NAME w90chk2chk.x)
    target_link_libraries(qe_w90chk2chk_exe PRIVATE Wannier90::wannier90)

    ###########################################################

    add_custom_target(w90
        DEPENDS
            qe_wannier90 Wannier90_exe Wannier90_post qe_w90chk2chk_exe
        COMMENT
            "Maximally localised Wannier Functions")

    qe_install_targets(
        # Executables
        Wannier90_exe Wannier90_post qe_w90chk2chk_exe)
else()
    find_package(Wannier90 REQUIRED)
    target_link_libraries(qe_wannier90 INTERFACE Wannier90::Wannier90)
endif()
