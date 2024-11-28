Developers comments for gbnsr6 

All my modification, or most of them are labeled
by "B.Aguilar". Here are some comments of my modifications:

- I added a new file called "newtop" to the block /files/. The
the parameters chunk, necks, and facts, required by Molecular dynamics
simulartions with AR6, are printed in the
file referenced by "newtop". (this functionality is not used now 
will be important for igb9 that has been tested.)

- In order to read the flag '-md' and the file name for "newtop" 
we added a couple of lines in the source file "gbnsr6.F90".
The default name for "newtop" is ' '. This will be used to check
if the user wants to generate a new topology file for a md simulation
with AR6.

- In "egb.F90" we check if "newtop" is equal to ' '. If that is
the case then we do nothing, otherwise, if the  user provided a good 
file name after the flag "-md", we compute chunks factors and print it
in the the file name in "newtop".


   

