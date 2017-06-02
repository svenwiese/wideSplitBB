// Authors:	Sven Wiese
//-----------------------------------------------------
#ifndef holes_H
#define holes_H

int
readholes (const char *filename,
           int        *nvarswithholes_p,
           double     **lbvarswithholes_p,
           double     **ubvarswithholes_p,
           int        **nholesvarswithholes_p,
           char       ***namevarswithholes_p,
           double     ***lbholes_p,
           double     ***ubholes_p);

#endif

