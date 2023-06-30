
#include "AveragedOrbitalElementsDriver.hpp"

AveragedOrbitalElementsDriver::AveragedOrbitalElementsDriver() : 
    CsaltDriver("AveragedOrbitalElements")
{
    SetParameters();
}

AveragedOrbitalElementsDriver::AveragedOrbitalElementsDriver(const std::string &testName, const Integer snoptConsoleOutputLevel) :
    CsaltDriver(testName, snoptConsoleOutputLevel)
{
    SetParameters();
}

AveragedOrbitalElementsDriver::~AveragedOrbitalElementsDriver()
{

}

void AveragedOrbitalElementsDriver::SetParameters()
{
    // Define problem parameters
    LU      = 6378.0;
    mu      = 3.986e5;
    TU      = GmatMathUtil::Sqrt(LU*LU*LU / mu);
    m0      = 1200.0; 
    P       = 5.0e3;
    Isp     = 1800.0;
    g0      = 9.80664;
    eta     = 0.55;
    tMax    = 2.0*eta*P / (g0 * Isp);

    // Set number of points use in Gauss quadrature
    n       = 10;

    // Set initial state vector constraint
    x0_con.SetSize(6);
    x0_con[0] = 11359.07 / LU;
    x0_con[1] = 0.7306;
    x0_con[2] = 0.0;
    x0_con[3] = 0.2539676;
    x0_con[4] = 0.0;
    x0_con[5] = m0;

    // Set final state vector constraint
    xf_con.SetSize(5);
    xf_con[0] = 42165.0 / LU;
    xf_con[1] = 0.0;
    xf_con[2] = 0.0;
    xf_con[3] = 0.0;
    xf_con[4] = 0.0;

    // Scale variables
    mu      *= TU*TU / (LU*LU*LU);
    Isp     *= 1.0 / TU; 
    g0      *= TU*TU / (LU*1000.0);
    tMax    *= TU*TU / (LU*1000.0);
}

void AveragedOrbitalElementsDriver::SetPointPathAndProperties()
{

}

void AveragedOrbitalElementsDriver::SetupPhases()
{

}
