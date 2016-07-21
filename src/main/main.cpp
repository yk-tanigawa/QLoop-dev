//
// Created by yosuke on 7/21/16.
//

#include <iostream>

#include "config.h"

using namespace std;

int version() {
    std::cout << "version " << VERSION_MAJOR << "." << VERSION_MINOR << std::endl;
}


int main() {
    version();
}
