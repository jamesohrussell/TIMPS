# Tracked IMERG Precipitation Systems (TIPS) development codes
C++, Python, and shell scripts for the development of a set of tracked precipitation systems derived from gridded precipitation fields.

## C++ codes. Developed by Gregor Skok.
Forward in Time (FiT) tracking algorithm codes. Initially developed for NASA TRMM TMPA data. Has been used to track features in NASA IMERG data. Can be used to track features in any generic gridded variable.

## Python codes. Developed by James Russell.
These scripts take NASA IMERG data and with the FiT algorithm, can be used to identify and calculate various variables for features in IMERG. The python scripts take IMERG data and convert it into a format that can be read by the FiT algorithm. They then take the original IMERG data and the output of the FiT algorithm to produce files that represent precipitation features.

## Instructions
Instructions on how to use these codes and an example script to run all scripts can be found in run_TIPS.sh
