#include "julia_utils.hpp"

using namespace jluna;

// Function for initializing jluna and the CSALT Julia module
void initialize_jluna() {
    // Initialize jluna with --threads auto flag
    //initialize(JULIA_NUM_THREADS_AUTO);

    // Initialize jluna with single thread
    initialize();

    // Load CSALT/julia source code
    Main.include(JLSRC_INCLUDE);
};