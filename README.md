# SpTRSV
A project that implements naive and optimized algorithms to solve sparse triangular systems with single and multicore approaches.

This project uses OpenMp and the Boost Library (preferably boost_1_81_0) so they have to be installed on your system.

To run this project, first run: 
```
make
```
Note that the make file assumes that the boost library is located int the following path so feel free to change it if yours is installed in another path.
```
/usr/local/include/boost/
```

Next run:

```
./main pathToMatrix pathToRightVector

Note that the matrix and the vector have to be sparse .mtx files.
