
target("drawstuff")
    set_languages("c99", "c++11")
    set_kind("static")

    on_load(function (target)
        import("lib.detect.find_package")

        -- target:add(find_package("rcssnet3D", {
        --     linkdirs={"$(env SPARK_DIR)/lib/simspark", "/usr/local/lib/simspark"},
        --     includedirs={"$(env SPARK_DIR)/include/simspark", "/usr/local/include/simspark"},
        --     links={"rcssnet3D"},
        --     includes={"rcssnet/addr.hpp"}
        -- }))
        -- target:add(find_package("boost_system"))
        target:add(find_package("X11"))
        target:add(find_package("GL"))
        target:add(find_package("GLU"))
    end)

    add_defines(
        "DEFAULT_PATH_TO_TEXTURES=\"$(top_srcdir)/drawstuff/textures/\"",
        "GL_SILENCE_DEPRECATION"
    )
    -- add files
    add_files(
        "src/drawstuff.cpp",
        "src/x11.cpp"
    )
    add_includedirs(
        ".",
        "$(projectdir)/ode/src"
    )
    add_sysincludedirs(
        "$(projectdir)/include",
        "$(projectdir)/ou/include"
    )


target("dstest")
    set_languages("c99", "c++11")
    set_kind("binary")

    on_load(function (target)
        import("lib.detect.find_package")
    end)
    add_sysincludedirs(
        "$(projectdir)/include"
    )
    add_includedirs(
        ".",
        "$(projectdir)/ode/src"
    )
    add_deps("drawstuff")
    add_files(
        "dstest/dstest.cpp"
    )