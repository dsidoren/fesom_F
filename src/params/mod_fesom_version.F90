module mod_fesom_version
    ! Seed library module + build/version provenance.
    implicit none
    private
    public :: fesom3_version_string

    character(len=*), parameter :: FESOM3_VERSION = "0.0.1-m0"

contains

    pure function fesom3_version_string() result(s)
        character(len=:), allocatable :: s
        s = "FESOM3 "//FESOM3_VERSION
    end function fesom3_version_string

end module mod_fesom_version
