
target("ou")
    set_languages("c99", "c++11")
    set_kind("static")

    on_load(function (target)
        import("lib.detect.find_package")
        target:add(find_package("pthread"))
    end)

    add_cxxflags("-fno-exceptions", "-fno-rtti")
    add_ldflags("-fno-exceptions", "-fno-rtti")

    -- add files
    add_files(
        "src/ou/*.cpp"
    )
    add_includedirs(
        "."
        ,"$(projectdir)/ode/src"
    )
    add_sysincludedirs(
        "$(projectdir)/include",
        "$(projectdir)/ou/include"
    )
    add_defines(
        "dOU_ENABLED"
        ,"_OU_NAMESPACE=odeou"
        ,"_OU_FEATURE_SET=_OU_FEATURE_SET_TLS"
        ,"_OU_TARGET_OS=_OU_TARGET_OS_GENUNIX"
    )

target("outest")
    set_languages("c99", "c++11")
    set_kind("binary")

    add_files(
        "test/*.cpp"
    )
    add_defines(
        "dOU_ENABLED"
        ,"_OU_NAMESPACE=odeou"
        ,"_OU_FEATURE_SET=_OU_FEATURE_SET_TLS"
        ,"_OU_TARGET_OS=_OU_TARGET_OS_GENUNIX"
    )
    add_deps("ou")
    add_sysincludedirs(
        "$(projectdir)/include",
        "$(projectdir)/ou/include"
    )

    -- add_cxxflags("-fno-exceptions", "-fno-rtti")
    -- add_ldflags("-fno-exceptions", "-fno-rtti")
