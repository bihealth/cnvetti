#include "cnvetti/version.h"

#include <iostream>

int main()
{
    std::cerr << "GIT_VERSION: " << GIT_VERSION << "\n"
              << "GIT_VERSION_SHORT: " << GIT_VERSION_SHORT << "\n";
    return 0;
}
