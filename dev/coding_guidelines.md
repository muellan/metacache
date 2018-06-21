MetaCache Coding Guidelines
===========================

Project
------------------
 - DON'T introduce breaking changes to the command line interface
 - DON'T introduce breaking changes to the (default) output format
 - try to avoid breaking database format changes
 - define all entities within namespace ```mc``` 


C++ Language Usage
------------------
 - C++14

 - NO C++17 features (yet)

 - NO **owning** raw pointers
   => Every resource must be cleaned up if it's owner is destroyed.

 - NO explicit ```new``` and ```delete``` 
   (except they are well *encapsulated* inside allocators / data structures)
 
 - NO ```using std::``` statements in *headers*!
 - DON'T use macros to ```#define``` constants or "functions" - EVER!
 - DON'T use ```#pragma once```, use include guards
 - DON'T use ```typedef A B;```, use ```using B = A;```

 - DON'T return meaningless pairs or tuples from functions
 - DON'T use **unscoped** enums, use **scoped** enums: ```enum class { ... }```
 - avoid functions with more than 5 parameters
 - avoid out-parameters (non-const reference function parameters)

 - std::vector should be your default container choice
 - prefer std::unordered_map/set over std::map/set
 - prefer range-based for loops
 - prefer std:: algorithms 
 - prefer free standing functions (over member functions)
 - prefer templates over std::function


Coding Style
------------------
 - indentation: 4 SPACES, NO Tabs
 - try to keep line lengths under 80-100 characters
 - file extensions: headers: ".h", TUs: ".cpp"
 - naming:
    - localVariables
    - memberVariables_
    - function_names        (should be verbs)
    - class_names           (should be nouns)
    - TemplateParameters    
    - DON'T use "_" at the begginning of names;
      these are reserved for std:: library entities.

```cpp
#include <vector>

namespace mc {

/**
 * @brief doxygen-style comment
 */
template<class MyTemplateParam>
class my_class_name
{
public:
    struct simple_data {
        int x;
        int y;
    };

    void my_function_name(int myFunctionParam) {
        if(myFunctionParam > 0) {  
            for(const auto& l: myMember_)  {
                ...
            }
        }
    }

private:
    //member names **end** in "_", 
    std::vector<simple_data> myMember_;
};


template<class TParam>
void foo() {
...
}


} //namespace mc
```

