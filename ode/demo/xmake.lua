

demos = {
    "demo_boxstack"
    ,"demo_buggy"
    ,"demo_cards"
    ,"demo_chain2"
    ,"demo_collision"
    ,"demo_convex_cd"
    ,"demo_crash"
    ,"demo_cylvssphere"
    -- ,"demo_dball"
    -- ,"demo_dhinge"
    -- ,"demo_transmission"
    ,"demo_feedback"
    ,"demo_friction"
    ,"demo_gyroscopic"
    ,"demo_gyro2"
    ,"demo_heightfield"
    ,"demo_hinge"
    ,"demo_I"
    ,"demo_jointPR"
    ,"demo_joints"
    ,"demo_jointPU"
    ,"demo_kinematic"
    ,"demo_motion"
    ,"demo_motor"
    ,"demo_ode"
    ,"demo_piston"
    ,"demo_plane2d"
    -- ,"demo_rfriction"
    ,"demo_slider"
    ,"demo_space"
    ,"demo_space_stress"
    ,"demo_step"
    -- ,"demo_tracks"
}

function _add_target(t)
    if not os.exists(t .. ".cpp") then
        return
    end
    target(t)
    set_languages("c99", "c++11")
    set_kind("binary")
    add_cxxflags("-g")
    add_cuflags("-g")
    add_ldflags(
        "-ldl", "-lrt"
    )

    add_defines(
        "DRAWSTUFF_TEXTURE_PATH=\"$(projectdir)/drawstuff/textures\""
    )
    add_sysincludedirs(
        "$(projectdir)/include"
    )
    add_deps("drawstuff", "ode")
    add_files(
        t .. ".cpp"
    )    
end

function _add_target_c(t)
    if not os.exists(t .. ".c") then
        return
    end
    target(t)
    set_languages("c99", "c++11")
    set_kind("binary")

    add_defines(
        "DRAWSTUFF_TEXTURE_PATH=\"$(projectdir)/drawstuff/textures\""
    )
    add_sysincludedirs(
        "$(projectdir)/include"
    )
    add_deps("drawstuff", "ode")
    add_files(
        t .. ".c"
    )
end

target("demos")
for i,v in ipairs(demos) do
    -- print(v)
    _add_target(v)
end

_add_target_c("demo_chain1")

_add_target("cuda_demo_ode")
-- add_files("$(projectdir)/ode/src/cuda_demo_helper.c")
add_ldflags(
    "-ldl", "-lrt"
)
