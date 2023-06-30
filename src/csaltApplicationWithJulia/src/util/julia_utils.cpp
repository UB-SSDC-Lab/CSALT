#include "julia_utils.hpp"

// =================================================
//  void setup_julia()
// =================================================
//
// Initializes Julia, activates the CSALT/Julia 
// environment (define in jlsrc/Project.toml), and 
// includes all of the Julia source files.
// 
// =================================================

void setup_julia()
{
    // Setup Julia context
    jl_init();

    // Instantiate string to store Julia code
    std::string c = "";

    // Activate julia csalt environment
    c += "using Pkg; Pkg.activate(\"" + JL_SRC + "\")";
    jl_eval_string(c.c_str());

    // Include source code in jlsrc
    c = "include(\"" + JL_SRC + "/include.jl\")";
    jl_eval_string(c.c_str());
}


// =================================================
//  void stop_julia()
// =================================================
//
// Pops all stored arrays to be garbage collected,
// handles Julia exceptions if any, and ends the 
// Julia process with exit ret val n.
// 
// =================================================

void stop_julia(int n)
{
    JL_GC_POP();
    handle_julia_exception();
    jl_atexit_hook(n);
}

void handle_julia_exception(void) {
    if (jl_value_t *ex = jl_exception_occurred()) {
        jl_printf(jl_stderr_stream(), "Exception: ");
        jl_call2(
            jl_get_function(jl_base_module, "showerror"),
            jl_stderr_obj(),
            ex
        );
        stop_julia(1);
        exit(1);
    }
}

jl_value_t* handle_eval_string(const char* code) {
    jl_value_t* result = jl_eval_string(code);
    handle_julia_exception();
    assert(result && "Missing return value but no exception occurred!");
    return result;
}

// ======================================
// JLvector
// ======================================

JLvector::JLvector(size_t n) : n(n)
{
    // Get array type
    jl_value_t* type = jl_apply_array_type((jl_value_t*) jl_float64_type, 1);

    // Allocate vector
    vecptr = jl_alloc_array_1d(type, n);
    vec    = (double *) jl_array_data(vecptr);

    // Add vector to GC
    JL_GC_PUSH1(&vecptr);
}

void JLvector::MakeZeroVector() 
{
    for (size_t i = 0; i < n; i++)
        vec[i] = 0.0;
}