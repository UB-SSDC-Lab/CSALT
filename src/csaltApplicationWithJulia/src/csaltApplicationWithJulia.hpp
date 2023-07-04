
#ifndef CSALT_APPLICATION_WITH_JULIA_HPP
#define CSALT_APPLICATION_WITH_JULIA_HPP

#include "csalt.hpp"
#include "jluna.hpp"
#include "julia_utils.hpp"
#include "drivers/CsaltDriver.hpp"
#include "drivers/DebrisDeorbitDriver.hpp"
#include "drivers/AveragedOrbitalElementsDriver.hpp"

// Makes Julia code run fast
//JULIA_DEFINE_FAST_TLS

int main(int argc, char *argv[]);
int unsafe_main(int argc, char *argv[]);

#endif