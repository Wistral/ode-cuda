
target("tests")
    set_languages("c99", "c++11")
    set_kind("binary")

    on_load(function (target)
        import("lib.detect.find_package")
        target:add(find_package("pthread"))
        -- target:add(find_package("dl"))
        -- target:add(find_package("rt"))
    end)
    add_defines("dTRIMESH_GIMPACT")
    add_sysincludedirs(
        "$(projectdir)/include",
        "$(projectdir)/ou/include"
    )
    add_includedirs(
        ".",
        "UnitTest++/src",
        "$(projectdir)/ode/src"
    )
    add_deps("ou", "ode")
    add_files(
        "*.cpp",
        "UnitTest++/src/*.cpp",
        "UnitTest++/src/Posix/*.cpp",
        "joints/*.cpp"
    )
    -- add_links("dl", "rt", "cudart_static")