# wideSplitBB
A branch-and-bound code based on wide split disjunctions.

---
This code implements the algorithm introduced in Section 4.3 of my [PhD thesis](http://amsdottorato.unibo.it/7612/1/wiese_sven_tesi.pdf). Note that there is no cut separation in this repository's version. This has been tested to compile on Mac OS X but should work on any UNIX system with a working CPLEX installation. Adjust the CPLEX path in the Makefile in order to compile. After compilation, type
~~~
./branch-and-hole
~~~
in order to display the set of possible command line options.

In instances, there is an example instance, composed of an .mps file and a .txt file in custom format that contains information on valid simple wide split disjuctions.
