# Laboratory reaction-event bridge reproduction

See ../../REACTION_EVENT.md for energy conventions and the scope boundary.
The C++ test exercises400 events across5 channels x2 conventions.

Build reaction_parent_reference_driver.cpp against ../../../include and
/path/to/build/libpb11.a, naming the executable reaction_parent_reference_driver
in this directory. Run reference_reaction_parent.py with mpmath1.3.0.
It generates50 deterministic input cases and compares the parent to80-digit
direct invariant formulas, writing reaction-parent-reference.json. Canonical
mass/Q values are inputs to this bridge test; their separate nuclear-data
validation is not replaced by it.

The four *-test.log files are final full-suite passes. consumer/ is a separate
CMake project using PB11_PREFIX and ifx -r8; build and run both executables.
No complete reaction probability distribution is selected in these event tests.
