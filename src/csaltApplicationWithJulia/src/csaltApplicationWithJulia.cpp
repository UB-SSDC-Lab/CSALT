
#include "csaltApplicationWithJulia.hpp"

// int main(int argc, char *argv[])
// {
//     initialize_jluna();

//     return 0;
// }

int main(int argc, char *argv[])
{
    // Initialize jluna
    initialize_jluna();

    // Instantiate some requirements
    std::string msg;
    std::string prob;
    std::string inFile = "";

    // Instantiate CsaltDriver null pointer
    CsaltDriver *driver = NULL;

    // If we have command line arguments, we need to process them
    // to determine what we need to do.
    if (argc > 1)
    {
        for (int i = 1; i < argc; i++)
        {
            std::string argument = argv[i];

            if (argument == "--run" || argument == "-r")
            {
                // Check that i != argc - 1
                if (i == (argc - 1))
                {
                    std::cout << "Error: No problem specified for --run argument." << std::endl;
                    return 1;
                }

                prob = argv[i + 1];
                i++;
            }
            else if (argument == "--input" || argument == "-i")
            {
                // Check that i != argc - 1
                if (i == (argc - 1))
                {
                    std::cout << "Error: No file specified for --input argument." << std::endl;
                    return 1;
                }

                inFile = argv[i + 1];
                i++;
            }
            else if (argument == "--help" || argument == "-h")
            {
  				msg = "";
				msg += "\n********************************************\n";
				msg += "CSALTApplicationWithJulia\n";
				msg += "********************************************\n\n";
				msg += "This program runs CSALT problems. If used\n";
				msg += "without any command line arguments, it will offer\n";
				msg += "a list of problems to choose from.\n";
				msg += "\n";
				msg += "The following command line arguments are \n";
				msg += "available:\n";
				msg += "\n";

				msg += "[-r | --run]: Run the specified problem.\n";
				msg += "Currently, the following problems are \n";
				msg += "available:\n";
				msg += "   1.  DebrisDeorbit\n";
                msg += "   2.  AveragedOrbitalElements";
				msg += "\n";

                msg += "[-i | --input]: Set input file for problem.\n";
                msg += "If none specified, will used default input file.\n";
                msg += "\n";

				msg += "[-h | --help]: Displays help menu.\n";

				std::cout << msg << std::endl;              
            }
        }
    }
    else 
    {
        msg = "";
        msg += "\n********************************************\n";
        msg += "CSALTApplicationWithJulia\n";
        msg += "********************************************\n\n";
        msg += "Select a problem:\n";
        msg += "  1.  DebrisDeorbit\n";
        msg += "  2.  AveragedOrbitalElements\n";

        std::cout << msg << std::endl;
        std::cout << "Input desired problem (name or number): ";
        //std::cin >> prob;
        prob = "2";
    }

    // Create driver if valid problem was specified
    if (prob == "1" || prob == "DebrisDeorbit") {
        driver = new DebrisDeorbitDriver("DebrisDeorbit", 1);

        // Set input file
        if (inFile == "")
            driver->SetInputFile("./../data/csalt/ddo_input.txt");
        else
            driver->SetInputFile(inFile);
    }
    else if (prob == "2" || prob == "AveragedOrbitalElements") {
        driver = new AveragedOrbitalElementsDriver();
    }
    else {
        std::cout << "Error: Invalid problem specified.\n";
        return 1;
    }

    // Run driver
    driver->Run();

    return 0;
}