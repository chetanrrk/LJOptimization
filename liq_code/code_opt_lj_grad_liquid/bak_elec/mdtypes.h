/*
 * Copyright (C) 2003-2004 by David J. Hardy.  All rights reserved.
 *
 * mdtypes.h - data types and constants for MDAPI and related codes
 */

#ifndef MDTYPES_H
#define MDTYPES_H

#include <limits.h>
#include <math.h>

  /*
   * Version number string
   */
#define MDAPI_VERSION  "0.9.0"

  /*
   * Name string for MDAPI
   */
#define MDAPI_STRING  "MDAPI v" MDAPI_VERSION

  /*
   * The predefined MDAPI types provide a compact way of describing
   * CHARMM-style force field parameters, topology, and trajectory
   * information for a macromolecular system.
   *
   * These data types are intended as elements in a collection of
   * interrelated arrays that describe the molecular system.  The force
   * field parameter data types define the fundamental constants used
   * in the CHARMM force field.  The topology data structures indicate
   * the bonded connections between particular atoms.  The "atom" array of
   * type MD_Atom gives the implicit ordering for all N atoms of the system.
   * The system trajectory is given in N-length arrays "pos" (position)
   * and "vel" (velocity) with the same ordering as the "atom" array.
   *
   * The "prm" field in the topology data types is an index mapping onto
   * the corresponding parameter array:
   *
   *     MD_Atom "atom"     corresponds with       MD_AtomPrm "atomprm"
   *     MD_Bond "bond"                            MD_BondPrm "bondprm"
   *     MD_Angle "angle"                          MD_AnglePrm "angleprm"
   *     MD_Tors "dihed"                           MD_TorsPrm "dihedprm"
   *     MD_Tors "impr"                            MD_TorsPrm "imprprm"
   *
   * CHARMM force field distinguishes between two kinds of torsion angles:
   * dihedral and improper.  These are maintained in separate lists.
   * 
   * MD_Excl "excl" does not have a corresponding force field parameter data 
   * type.  Also, MD_NbfixPrm "nbfixprm" is not indexed by any topology data 
   * type but is instead used for modifying particular pairs of nonbonded
   * parameter entries in MD_AtomPrm.  
   *
   * Each parameter data type occurs on a per atom-type (or per bond-type)
   * basis.  
   *
   *   MD_AtomPrm "atomprm"   ---  Params given on a per atom-type basis.  
   *   MD_BondPrm "bondprm"   ---  Params given on a per bond-type basis.  
   *   MD_AnglePrm "angleprm" ---  Params given on a per angle-type basis.  
   *   MD_TorsPrm "dihedprm"  ---  Params given on a per dihedral-type basis.
   *   MD_TorsPrm "imprprm"   ---  Params given on a per improper-type basis.
   *   MD_NbfixPrm "nbrixprm" ---  These params override the van der Waals 
   *                               constants for a pair of atom-types.
   *                               Contains pair of indices into MD_AtomPrm.
   *
   * The topology data types occur on a per atom (or per bond) basis.  
   * The MD_Atom array provides the implicit ordering of atoms in the system.  
   * The "atom[]" member of MD_Bond, MD_Angle, MD_Tors, and MD_Excl types
   * identifies the connection between specific atoms by indexing the MD_Atom
   * array.  The "prm" field gives the index into the corresponding force
   * field parameter array.  
   *
   *   MD_Atom "atom"   ---  Constants defined on a per atom basis.
   *   MD_Bond "bond"   ---  Linear bonds involve two atoms.  
   *   MD_Angle "angle" ---  Angular bonds involve three atoms.  
   *   MD_Tors "dihed"  ---  Dihedral bonds involve four atoms.  
   *   MD_Tors "impr"   ---  Improper bonds involve four atoms.  
   *   MD_Excl "excl"   ---  Pairs of atoms to exclude from nonbonded 
   *                         interactions.  
   *
   * Dihedral and improper bond multiplicity is represented by storing
   * the sets of parameters consecutively in MD_TorsPrm array.  The "mult"
   * member gives the count of parameters, as listed in subarray order.
   * For instance, suppose the "dihedprm" array has in its kth entry a set
   * of dihedral parameters with multiplicity 1.  Then
   *
   *   dihedprm[k] == { ... , mult == 1, ... }
   *
   * Now suppose the next array element begins a set of parameters for a
   * dihedral interaction with multiplicity 4.  This means that the
   * parameters would be stored in the next four entries with
   *
   *   dihedprm[k+1] == { ... , mult == 1, ... }
   *   dihedprm[k+2] == { ... , mult == 2, ... }
   *   dihedprm[k+3] == { ... , mult == 3, ... }
   *   dihedprm[k+4] == { ... , mult == 4, ... }
   *
   * In case the dihedral from the mth entry of the MD_Tors "dihed" array
   * uses this particular multiplicity 4 parameter set, the entry would
   * index the last entry of the set,
   *
   *   dihed[m] == { ... , prm == k+4, ... }
   *
   * The trajectory information, such as position and velocity, is stored 
   * in MD_Dvec or MD_Fvec arrays.  The ordering of these array elements
   * should agree with the ordering of the MD_Atom array.  
   *
   * The "name" and "type" members are intended for diagnostic purposes
   * only and should not really be necessary for correct functioning of
   * the engine.  However, the front end does have the responsibility to
   * make sure that the "prm" and "atom" members correctly index between
   * arrays, so that every topology array element correctly shows which
   * atoms are bonded together and the corresponding set of force field
   * parameters.
   */


  /*
   * Common units:
   *   time        - picoseconds (ps, 1 ps = 10^{-12} seconds)
   *   timestep    - femtoseconds (fs, 1 fs = 10^{-3} ps = 10^{-15} seconds)
   *   length      - Angstroms (A)
   *   energy      - kilocalories per mole (kcal/mol)
   *   mass        - atomic mass units (AMU)
   *   charge      - electron charge (e)
   *   temperature - degrees Kelvin (K)
   *
   * MDAPI units:
   *   timestep    - fs
   *   numsteps    - (unitless integer)
   *   position    - A
   *   velocity    - A/fs
   *   force       - kcal/mol/A
   *   energy      - kcal/mol
   *   temperature - K
   *
   * (force field parameter constants have their units listed below)
   */

  /* convert time from fs to ps */
#define MD_PICOSEC  0.001

  /* convert time from ps to fs */
#define MD_FEMTOSEC  1000.0

  /*
   * convert velocity from A/ps to A/fs
   * (use to convert typical file storage units to MDAPI veloc units)
   */
#define MD_ANGSTROM_FS  MD_PICOSEC

  /*
   * convert velocity from A/fs to A/ps
   * (use to convert MDAPI veloc units back to file storage units)
   */
#define MD_ANGSTROM_PS  MD_FEMTOSEC

  /* convert energy from kcal/mol into derived units AMU*A^2/fs^2 */
#define MD_ENERGY_CONST  0.0004184

  /* convert force from kcal/mol/A into derived units AMU*A/fs^2 */
#define MD_FORCE_CONST  MD_ENERGY_CONST

  /* convert energy from derived units AMU*A^2/fs^2 into kcal/mol */
#define MD_KCAL_MOL  (1.0 / MD_ENERGY_CONST)

  /* convert force from derived units AMU*A/fs^2 into kcal/mol/A */
#define MD_KCAL_MOL_A  MD_KCAL_MOL

  /* Coulomb's constant for electrostatics, units kcal*A/mol/e^2 */
//#define MD_COULOMB  332.0636

  /* Boltzmann constant for temperature, units kcal/mol/K */
#define MD_BOLTZMAN  0.001987191
/* #define MD_BOLTZMAN  0.0019872067 */

  /* pi */
#ifdef M_PI
#define MD_PI  M_PI
#else
#define MD_PI  3.1415926535897932384626433832795028841971693993751
#endif

  /* convert from angle degrees to radians */
#define MD_RADIANS  (MD_PI / 180.0)

  /* convert from radians to angle degrees */
#define MD_DEGREES  (180.0 / MD_PI)


  /*
   * Primary data types:
   *   char (1 byte)
   *   int32 (4 bytes = 32 bits)
   *   float (4 bytes) 
   *   double (8 bytes)
   *   MD_Fvec (3 floats = 12 bytes)
   *   MD_Dvec (3 doubles = 24 bytes)
   *   MD_Name (8-char array = 8 bytes)
   *   MD_String (64-char array = 64 bytes)      -- deprecate?
   *   MD_Message (512-char array = 512 bytes)   -- deprecate?
   */

  /* number of primary types */
  enum {
    MD_NUMBER_PRIMARY_TYPES = 9
  };

  /* predefined char array lengths */
  enum {
    MD_NAME_SIZE = 8,
    MD_STRING_SIZE = 64,
    MD_MESSAGE_SIZE = 512
  };

  /*
   * ensure that int32 is 32-bit (4 byte) integer (uint32 is unsigned)
   */
#if   INT_MAX == 2147483647
  typedef int int32;
  typedef unsigned int uint32;
  enum {
    MD_INT32_MAX = INT_MAX,
    MD_INT32_MIN = INT_MIN
  };
#elif INT_MAX == 32767 && LONG_MAX == 2147483647
  typedef long int32;
  typedef unsigned long uint32;
  enum {
    MD_INT32_MAX = LONG_MAX,
    MD_INT32_MIN = LONG_MIN
  };
#elif INT_MAX == 9223372036854775807L && SHRT_MAX == 2147483647
  typedef short int32;
  typedef unsigned short uint32;
  enum {
    MD_INT32_MAX = SHRT_MAX,
    MD_INT32_MIN = SHRT_MIN
  };
#endif

  typedef struct MD_Fvec_tag { float  x, y, z; } MD_Fvec;
  typedef struct MD_Dvec_tag { double x, y, z; } MD_Dvec;
  typedef char MD_Name    [ MD_NAME_SIZE    ];
  typedef char MD_String  [ MD_STRING_SIZE  ];
  typedef char MD_Message [ MD_MESSAGE_SIZE ];


  /*
   * Derived data types:
   *   force field parameter data types
   *   topology data types
   */

  /* force field parameter data types */

  typedef struct MD_AtomPrm_tag {
    double emin;     /* Lennard-Jones energy min (kcal/mol) */
    double rmin;     /* Lennard-Jones distance for emin (A) */
    double emin14;   /* modified 1-4 energy min (kcal/mol) */
    double rmin14;   /* modified 1-4 distance for emin14 (A) */
    MD_Name type;    /* string to identify atom type */
  } MD_AtomPrm; 
  
  typedef struct MD_BondPrm_tag {
    double k;        /* spring coefficient (kcal/mol/A^2) */
    double r0;       /* equilibrium length (A) */
    MD_Name type[2]; /* strings to identify atom types */
  } MD_BondPrm;

  typedef struct MD_AnglePrm_tag {
    double k_theta;  /* coefficient for theta (kcal/mol/rad^2) */
    double theta0;   /* equilibrium angle (radians) */
    double k_ub;     /* coef for Urey-Bradley term (kcal/mol/A^2) */
    double r_ub;     /* equil length for Urey-Bradley term (A) */
    MD_Name type[3]; /* strings to identify atom types */
  } MD_AnglePrm;

  typedef struct MD_TorsPrm_tag {
    double k_tor;    /* torsion coef (kcal/mol for n>0 OR kcal/mol/rad^2) */
    double phi;      /* phase shift (radians) */
    int32 n;         /* periodicity */
    int32 mult;      /* multiplicity of torsion */
    MD_Name type[4]; /* strings to identify atom types */
  } MD_TorsPrm;

  typedef struct MD_NbfixPrm_tag {
    double emin;     /* Lennard-Jones energy min (kcal/mol) */
    double rmin;     /* Lennard-Jones distance for emin (A) */
    double emin14;   /* modified 1-4 energy min (kcal/mol) */
    double rmin14;   /* modified 1-4 distance for emin14 (A) */
    int32 prm[2];    /* index MD_AtomPrm array */
    MD_Name type[2]; /* strings to identify atom types */
  } MD_NbfixPrm;


  /* topology data types */

  typedef struct MD_Atom_tag {
    double m;        /* mass (AMU) */
    double q;        /* charge (e) */
    int32 prm;       /* index MD_AtomPrm array */
    char is_hydrogen;  /* set TRUE for light atoms (having m < 1.5) */
    char is_water;     /* set TRUE for the oxygen in a water molecule */
    char is_fixed;     /* set TRUE for fixed atom */
    char is_other;     /* set TRUE for Drude particle */
    MD_Name name;    /* string to identify atom name */
    MD_Name type;    /* string to identify atom type name */
  } MD_Atom;

  typedef struct MD_Bond_tag {
    int32 atom[2];   /* index MD_Atom array */
    int32 prm;       /* index MD_BondPrm array */
  } MD_Bond;

  typedef struct MD_Angle_tag {
    int32 atom[3];   /* index MD_Atom array */
    int32 prm;       /* index MD_AnglePrm array */
  } MD_Angle;

  typedef struct MD_Tors_tag {
    int32 atom[4];   /* index MD_Atom array */
    int32 prm;       /* index MD_TorsPrm array */
  } MD_Tors;

  typedef struct MD_Excl_tag {
    int32 atom[2];   /* index MD_Atom array */
  } MD_Excl;


  /*
   * type ID numbers used by MDAPI
   */

  /* bit shifting constants for type encoding */
  enum { 
    MD_SHIFT_FLIP = 16,
    MD_SHIFT_TYPE = 19
  };

  /* type size mask (based on shifting values above) */
  enum {
    MD_SIZEOF_MASK = 0xFFFF
  };

  /* encoded flags to identify byte reordering information */
  enum {
    MD_FLIP_NONE    = 0,
    MD_FLIP_4BYTE   = 1 << MD_SHIFT_FLIP,
    MD_FLIP_8BYTE   = 2 << MD_SHIFT_FLIP,
    MD_FLIP_12BYTE  = 3 << MD_SHIFT_FLIP,
    MD_FLIP_24BYTE  = 4 << MD_SHIFT_FLIP,
    MD_FLIP_DERIVED = 5 << MD_SHIFT_FLIP,
    MD_FLIP_MASK    = 7 << MD_SHIFT_FLIP
  };

  /*
   * constants to uniquely identify each data type
   * (includes size information in bytes)
   */
  enum {
    MD_CHAR
      =  (0 << MD_SHIFT_TYPE) | MD_FLIP_NONE    | (int32)sizeof(char),
    MD_INT32
      =  (1 << MD_SHIFT_TYPE) | MD_FLIP_4BYTE   | (int32)sizeof(int32),
    MD_FLOAT
      =  (2 << MD_SHIFT_TYPE) | MD_FLIP_4BYTE   | (int32)sizeof(float),
    MD_DOUBLE
      =  (3 << MD_SHIFT_TYPE) | MD_FLIP_8BYTE   | (int32)sizeof(double),
    MD_FVEC
      =  (4 << MD_SHIFT_TYPE) | MD_FLIP_12BYTE  | (int32)sizeof(MD_Fvec),
    MD_DVEC
      =  (5 << MD_SHIFT_TYPE) | MD_FLIP_24BYTE  | (int32)sizeof(MD_Dvec),
    MD_NAME
      =  (6 << MD_SHIFT_TYPE) | MD_FLIP_NONE    | (int32)sizeof(MD_Name),
    MD_STRING
      =  (7 << MD_SHIFT_TYPE) | MD_FLIP_NONE    | (int32)sizeof(MD_String),
    MD_MESSAGE
      =  (8 << MD_SHIFT_TYPE) | MD_FLIP_NONE    | (int32)sizeof(MD_Message),
    MD_ATOMPRM
      =  (9 << MD_SHIFT_TYPE) | MD_FLIP_DERIVED | (int32)sizeof(MD_AtomPrm),
    MD_BONDPRM
      = (10 << MD_SHIFT_TYPE) | MD_FLIP_DERIVED | (int32)sizeof(MD_BondPrm),
    MD_ANGLEPRM
      = (11 << MD_SHIFT_TYPE) | MD_FLIP_DERIVED | (int32)sizeof(MD_AnglePrm),
    MD_TORSPRM
      = (12 << MD_SHIFT_TYPE) | MD_FLIP_DERIVED | (int32)sizeof(MD_TorsPrm),
    MD_NBFIXPRM
      = (13 << MD_SHIFT_TYPE) | MD_FLIP_DERIVED | (int32)sizeof(MD_NbfixPrm), 
    MD_ATOM
      = (14 << MD_SHIFT_TYPE) | MD_FLIP_DERIVED | (int32)sizeof(MD_Atom),
    MD_BOND
      = (15 << MD_SHIFT_TYPE) | MD_FLIP_DERIVED | (int32)sizeof(MD_Bond),
    MD_ANGLE
      = (16 << MD_SHIFT_TYPE) | MD_FLIP_DERIVED | (int32)sizeof(MD_Angle),
    MD_TORS
      = (17 << MD_SHIFT_TYPE) | MD_FLIP_DERIVED | (int32)sizeof(MD_Tors),
    MD_EXCL
      = (18 << MD_SHIFT_TYPE) | MD_FLIP_DERIVED | (int32)sizeof(MD_Excl)
  };

  /* macro to return byte size on type constants */
#define MD_SIZEOF(type)  ((type) & MD_SIZEOF_MASK)


#endif  /* MDTYPES_H */
