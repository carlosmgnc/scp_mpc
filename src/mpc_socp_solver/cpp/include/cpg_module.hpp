
/*
Auto-generated by CVXPYgen on May 26, 2025 at 22:53:23.
Content: Declarations for Python binding with pybind11.
*/

// User-defined parameters
struct mpc_socp__CPG_Params_cpp_t {
    std::array<double, 476> param3972;
    std::array<double, 14> param3965;
    std::array<double, 490> param3966;
    std::array<double, 6664> param3969;
    std::array<double, 1428> param3970;
    std::array<double, 102> param3967;
    std::array<double, 34> param3968;
};

// Flags for updated user-defined parameters
struct mpc_socp__CPG_Updated_cpp_t {
    bool param3972;
    bool param3965;
    bool param3966;
    bool param3969;
    bool param3970;
    bool param3967;
    bool param3968;
};

// Primal solution
struct mpc_socp__CPG_Prim_cpp_t {
    std::array<double, 490> var3963;
    std::array<double, 102> var3964;
};

// Dual solution
struct mpc_socp__CPG_Dual_cpp_t {
    std::array<double, 14> d0;
    std::array<double, 14> d1;
    double d2;
    double d3;
    double d4;
    std::array<double, 14> d5;
    double d6;
    double d7;
    double d8;
    std::array<double, 14> d9;
    double d10;
    double d11;
    double d12;
    std::array<double, 14> d13;
    double d14;
    double d15;
    double d16;
    std::array<double, 14> d17;
    double d18;
    double d19;
    double d20;
    std::array<double, 14> d21;
    double d22;
    double d23;
    double d24;
    std::array<double, 14> d25;
    double d26;
    double d27;
    double d28;
    std::array<double, 14> d29;
    double d30;
    double d31;
    double d32;
    std::array<double, 14> d33;
    double d34;
    double d35;
    double d36;
    std::array<double, 14> d37;
    double d38;
    double d39;
    double d40;
    std::array<double, 14> d41;
    double d42;
    double d43;
    double d44;
    std::array<double, 14> d45;
    double d46;
    double d47;
    double d48;
    std::array<double, 14> d49;
    double d50;
    double d51;
    double d52;
    std::array<double, 14> d53;
    double d54;
    double d55;
    double d56;
    std::array<double, 14> d57;
    double d58;
    double d59;
    double d60;
    std::array<double, 14> d61;
    double d62;
    double d63;
    double d64;
    std::array<double, 14> d65;
    double d66;
    double d67;
    double d68;
    std::array<double, 14> d69;
    double d70;
    double d71;
    double d72;
    std::array<double, 14> d73;
    double d74;
    double d75;
    double d76;
    std::array<double, 14> d77;
    double d78;
    double d79;
    double d80;
    std::array<double, 14> d81;
    double d82;
    double d83;
    double d84;
    std::array<double, 14> d85;
    double d86;
    double d87;
    double d88;
    std::array<double, 14> d89;
    double d90;
    double d91;
    double d92;
    std::array<double, 14> d93;
    double d94;
    double d95;
    double d96;
    std::array<double, 14> d97;
    double d98;
    double d99;
    double d100;
    std::array<double, 14> d101;
    double d102;
    double d103;
    double d104;
    std::array<double, 14> d105;
    double d106;
    double d107;
    double d108;
    std::array<double, 14> d109;
    double d110;
    double d111;
    double d112;
    std::array<double, 14> d113;
    double d114;
    double d115;
    double d116;
    std::array<double, 14> d117;
    double d118;
    double d119;
    double d120;
    std::array<double, 14> d121;
    double d122;
    double d123;
    double d124;
    std::array<double, 14> d125;
    double d126;
    double d127;
    double d128;
    std::array<double, 14> d129;
    double d130;
    double d131;
    double d132;
    std::array<double, 14> d133;
    double d134;
    double d135;
    double d136;
};

// Solver information
struct mpc_socp__CPG_Info_cpp_t {
    double obj_val;
    int iter;
    int status;
    double pri_res;
    double dua_res;
    double time;
};

// Solution and solver information
struct mpc_socp__CPG_Result_cpp_t {
    mpc_socp__CPG_Prim_cpp_t prim;
    mpc_socp__CPG_Dual_cpp_t dual;
    mpc_socp__CPG_Info_cpp_t info;
};

// Main solve function
mpc_socp__CPG_Result_cpp_t mpc_socp__solve_cpp(struct mpc_socp__CPG_Updated_cpp_t& CPG_Updated_cpp, struct mpc_socp__CPG_Params_cpp_t& CPG_Params_cpp);
