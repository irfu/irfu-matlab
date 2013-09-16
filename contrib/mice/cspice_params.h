/*
-Procedure

   cspice_params.h - parameter values from SPICELIB

-Abstract

   Taken from pool.f version:

   SPICELIB Version 8.3.0, 22-DEC-2004 (NJB)


      MAXVAR      is the maximum number of variables that the
                  kernel pool may contain at any one time.
                  MAXVAR should be a prime number.

      MAXLEN      is the maximum length of the variable names
                  that can be stored in the kernel pool.

      MAXVAL      is the maximum number of distinct values that
                  may belong to the variables in the kernel pool.
                  Each variable must have at least one value, and
                  may have any number, so long as the total number
                  does not exceed MAXVAL. MAXVAL must be at least
                  as large as MAXVAR.

      MXNOTE      is the maximum number of distinct variable-agents
                  pairs that can be maintained by the kernel pool.
                  (A variable is "paired" with an agent, if that agent
                  is to be notified whenever the variable is updated.)

      MAXAGT      is the maximum number of agents that can be kept
                  on the distribution list for notification of updates
                  to kernel variables.

      MAXCHR      is the maximum number of characters that can be
                  stored in a component of a string valued kernel
                  variable.

      MAXLIN      is the maximum number of character strings that
                  can be stored as data for kernel pool variables.


   Taken from zzrvar.f version:

   SPICELIB Version 1.6.0, 06-AUG-2002 (BVS)

      LINLEN      is the maximum length of a line in the kernel file.

-Disclaimer

   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
   CALIFORNIA  INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE 
   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED 
   "AS-IS" TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING
   ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR
   A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE 
   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.

   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, 
   OR NASA BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, 
   BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF 
   ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY 
   AND LOST PROFITS, REGARDLESS OF WHETHER CALTECH, JPL, OR 
   NASA BE ADVISED, HAVE REASON TO KNOW, OR, IN FACT, SHALL 
   KNOW OF THE POSSIBILITY.

   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE 
   OF THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO 
   INDEMNIFY CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING 
   FROM THE ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.

-Required_Reading

   None.

-Keywords

   paramters

-Brief_I/O

    None.

-Detailed_Input

    None.

-Detailed_Output

    None.

-Parameters

    None.

-Exceptions

    None.

-Files

    None.

-Particulars

    KEEP THIS FILE SYNCHED WITH POOL.F and ZZRVAR.F.

-Examples

    None.

-Restrictions

    None.

-Literature_References

    None.

-Author_and_Institution

   E. D. Wright    (JPL)

-Version

   -CSPICE Version 1.2.0, 24-MAY-2010 (EDW) (NJB)

      Increased MAXVAL to 200000.
        
   -CSPICE Version 1.1.0, 23-FEB-2009 (EDW)

      Added LINLEN parameter.

   -CSPICE Version 1.0.0, 27-APR-2006 (EDW)

-Index_Entries

-&
*/

#define         MAXVAR         5003 
 
#define         MAXVAL         200000 
 
#define         MAXLIN         4000 
 
#define         MAXCHR         80 
 
#define         MXNOTE         2000 
 
#define         MAXLEN         32 
 
#define         MAXAGT         1000 

#define         LINLEN         132  


