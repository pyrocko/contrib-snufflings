HYPOSAT-snuffling
-----------------

This snuffling is a wrapper to the *HYPOSAT* localization program by Johannes 
Schweitzer. Before using this snuffling, download the source codes from the 
NORSAR ftp server:
[ftp://ftp.norsar.no/pub/outgoing/johannes/hyposat/](ftp://ftp.norsar.no/pub/outgoing/johannes/hyposat/)

Follow the installation instructions as described in the manual 
`hyposat_man.pdf` in Chapter 2. If you are using a linux
derivative such as Debian or Ubuntu, you most likely want to download one of
the compressed `hyposat_l.XXX` files (where "l" stands for linux).

#### Note

In the Makefile the first non-commented line by default sets G77 to use
as Fortran compiler. You might want to change this to gfortran or f77 on more
recent Debian/Ubuntu versions.

#### One Final Hint

*HYPOSAT* is capable of using more than just *P* or *S* phases. To exploit 
that you can change/extend the `phase_key_mapping` in your snuffler's 
configuration file located at `~/.pyrocko/snuffler.pf`.


