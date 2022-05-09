#!/bin/bash

# approx 52.000 lines of code

find . \( -name \*.cxx -or -name \*.h -or -name \*.hpp -or -name \*.cpp \) ! -path "./ThirdParty/*" ! -path "./build/*" ! -path "./release/*" ! -path "./build_Debug/*" ! -path "./sps/*" ! -path "./build_debug/*" | xargs -i -t less {} | wc -l
