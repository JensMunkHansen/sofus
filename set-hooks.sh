#!/bin/bash
cd .git/hooks
ln -sf ../../bin/githook-astyle.sh ./pre-commit
cd -

