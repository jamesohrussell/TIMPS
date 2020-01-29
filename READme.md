# Tracked IMERG Precipitation Systems (TIPS) development codes
C++ and Python codes for the development of a set of tracked features derived from gridded precipitation fields

## C++ codes
Forward in Time (FiT) tracking algorithm codes developed by Gregor Skok. Initially developed for NASA TRMM TMPA data. Has been used to track features in NASA IMERG data. Can be used to track features in any generic gridded variable.

## Python codes
These scripts take NASA IMERG data and with the FiT algorithm, can be used to identify and calculate various variables for features in IMERG. Developed by James Russell. The python scripts take IMERG data and convert it into a format that can be read by the FiT algorithm. They then take the original IMERG data and the output of the FiT algorithm to produce files that represent precipitation features.
