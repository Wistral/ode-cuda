
target("ode")
set_languages("c99", "c++17")
set_kind("static")

on_load(function (target)
    import("lib.detect.find_package")    
    target:add(find_package("cuda"
        ,{
            linkdirs={"$(env CUDA_HOME)/lib64", "/usr/local/cuda/lib64"}
            ,includedirs={"$(env CUDA_HOME)/include", "/usr/local/cuda/include"}
            ,links={"cuda", "cudart_static"}
    }))
end)

add_defines("__ODE__","dNODEBUG")
-- add files
add_files(
    -- "nextafterf.c",
    "array.cpp",
    "box.cpp",
    "capsule.cpp",
    "collision_cylinder_box.cpp",
    "collision_cylinder_plane.cpp",
    "collision_cylinder_sphere.cpp",
    "collision_kernel.cpp",
    "collision_quadtreespace.cpp",
    "collision_sapspace.cpp",
    "collision_space.cpp",
    "collision_transform.cpp",
    "collision_trimesh_disabled.cpp",


    "collision_util.cpp",
    "convex.cpp",
    "cylinder.cpp",
    -- "default_threading.cpp",
    "error.cpp",
    "export-dif.cpp",
    -- "fastdot.cpp",
    -- "fastldltfactor.cpp",
    -- "fastldltsolve.cpp",
    -- "fastlsolve.cpp",
    -- "fastltsolve.cpp",
    -- "fastvecscale.cpp",

    "fastldlt.c",
    "fastltsolve.c", 
    "fastdot.c", 
    "fastlsolve.c",

    "heightfield.cpp",
    "lcp.cpp",
    "mass.cpp",
    "mat.cpp",
    "matrix.cpp",
    "memory.cpp",
    "misc.cpp",
    -- "objects.cpp",
    "obstack.cpp",
    "ode.cpp",
    "odeinit.cpp",
    "odemath.cpp",
    "plane.cpp",
    "quickstep.cpp",
    "ray.cpp",
    -- "resource_control.cpp",
    "rotation.cpp",
    -- "simple_cooperative.cpp",
    "sphere.cpp",
    "step.cpp",
    "stepfast.cpp",
    "timer.cpp",
    -- "threading_base.cpp",
    -- "threading_impl.cpp",
    -- "threading_pool_posix.cpp",
    -- "threading_pool_win.cpp",
    "util.cpp",
    "joints/*.cpp"

    -- ,"$(projectdir)/ou/src/ou/*.cpp"
    -- ,"$(projectdir)/ou/src/ou/customization.cpp"
    ,"odeou.cpp"
    ,"odetls.cpp"

    ,"*.cu"

)
add_cuflags("-g")
add_cxxflags("-g", "-Werror", "-Wconversion", "-Waddress")
add_cxflags("-g")
add_includedirs(
    "."
    ,"$(projectdir)/OPCODE"
)
add_sysincludedirs(
    "$(projectdir)/include",
    "$(projectdir)/ou/include"
)
add_deps("ou")
