# Mailing list
Communication about this project occurs on the lib2geom
[mailing list](https://lists.sourceforge.net/lists/listinfo/lib2geom-devel).


# Help Wanted
We greatly appreciate contributions.  You don't need to be a `math-whiz` or
`über-hacker` (though these are definitely appreciated :) ) to help.  The tasks
of code cleanup, consistency, testing, documentation, and toys mostly just
require perseverance, and benefit the project greatly.

As far as very specialized skill, we are always in need of mathy people, even if
it is just for their insight on problems and techniques (as opposed to coding).

# Coding Style
Please refer to the
[Coding Style Guidelines](http://www.inkscape.org/doc/coding_style.php)
if you have specific questions on the style to use for code. If reading style
guidelines doesn't interest you, just follow the general style of the
surrounding code, so that it is at least internally consistent.

# Compiling
For Windows instructions, see [README.win32.md](README.win32.md)

For Debian-like platforms, the following packages are required:
 - cairo v1.1.7 or later (Debian package libcairo2-dev)
 - cmake
 - make
 - libboost-dev
 - libgsl0-dev (though eventually it will only be required in tests)
 - refblas3* on dapper

To compile, use
```bash
cmake .
Make
```

If you have problems, just ask on the mailing list.

# Running tests
For Debian-like platforms, after compiling lib2geom, issue this command line
```bash
make test
```

# Adding a unit test
Make sure you write it using GTest syntax - look at e.g.
[tests/affine-test.cpp](tests/affine-test.cpp).

To add the test to the build, add the test to
[tests/CMakeLists.txt](tests/CMakeLists.txt) in under
`SET(2GEOM_GTESTS_SRC)`.
