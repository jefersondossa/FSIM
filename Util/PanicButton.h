/**
 * Defines a class to be called every time somethingh is not implemented yet.
 */
#ifndef PANICBUTTON
#define PANICBUTTON
#include <iostream>
#include <stdexcept>
#include <execinfo.h>
#include <stdio.h>

/// @brief Class to check consistency in the code. Returns an error code when called.
class PanicButton
{
public:
    //Error 
    PanicButton(){
        std::cout << "Put a break point here!" << '\n';

        void* callstack[128];
        int i, frames = backtrace(callstack, 128);
        char** strs = backtrace_symbols(callstack, frames);
        for (i = 0; i < frames; ++i) {
            printf("%s\n", strs[i]);
        }
        free(strs);

        throw std::bad_exception();
    };
};


#endif
