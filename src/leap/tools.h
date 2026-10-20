/*
 *      File:   tools.h
 *
 ************************************************************************
 *                                                                      *
 *     Designed by:    Christian Schafmeister                           *
 *     Author:         Christian Schafmeister                           *
 *                                                                      *
 *     VERSION: 1.0                                                     *
 *     Programmers:                                                     *
 *             Christian Schafmeister                                   *
 *             David Rivkin                                             *
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *      Description:
 *              This file contains several routines that
 *              perform different functions.
 *              These functions are currently called from
 *              the command line interpretor.
 */



                /* Different functions that can be performed during */
                /* a distance search */

#define DISTANCE_SEARCH_CREATE_BONDS    1
#define DISTANCE_SEARCH_PRINT_WARNINGS  2

#define DEFAULT_DISTANCE_SEARCH         2.0


                /* NOSHELL is used to specify a ridiculously thick */
                /* solvent shell, so that NO solvent is removed */
                /* because of being too far from the solute */

#define NOSHELL         999999.0

extern bool     bToolGeometricCenter( OBJEKT oObjekt, VECTOR *vPCenter );
extern int      iToolDistanceSearch( CONTAINER cCont, double dCloseness, 
                        bool bAbsoluteDistance, int iOperation );
extern LIST     lToolListOfResidues( UNIT uUnit, LIST lResidues );
extern void     ToolCenterUnitByRadii( UNIT uUnit, bool bOrient );
extern void     ToolInitSolventPotential( UNIT uUnit, VARARRAY vaSolvent,
                        int *iPMinPotRes, int *iPMaxPotRes );
extern void     ToolOctBoxCheck( UNIT uSolute, double *dPBuf, bool bMsg, bool bIsotropic );
extern void     ToolOrientPrincipleAxisAlongCoordinateAxis( UNIT uUnit );
extern void     ToolReplaceSolvent( UNIT uUnit, VARARRAY vaSolvent, int iSolv,
                        UNIT uIon, double dCharge, 
                        int *iPMinPotRes, int *iPMaxPotRes );
extern void     ToolSanityCheckBox( UNIT uUnit );
extern UNIT     zToolSetupSolvent( UNIT uSolvent );
extern void     zToolSolvateAndShell( UNIT uSolute, UNIT uSolvent,
                        double dXW, double dYW, double dZW, double dCloseness,
                        double dFarness, bool bShell, bool bClip, bool bOct,
                        bool bIsotropic );
extern void     ToolSolvateCell( UNIT uSolute, UNIT uSolvent, double dCloseness);
extern void     BuildFractionalTransforms( UNIT uUnit, MATRIX M, MATRIX Mi );

extern void     zToolSolvateInSphere( UNIT uSolute, UNIT uSolvent, 
                        VECTOR *vPCenter, double dRadius, double dCloseness );
extern void     ToolSetUnitBoxByCenters( UNIT uUnit );

