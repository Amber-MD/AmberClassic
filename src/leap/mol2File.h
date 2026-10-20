/*
 *      File:   mol2File.h
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
 *        mol2File.h      based on pdbFile.h
 * 
 *      Christine Cezard (2007)
 *
 *     Universite de Picardie - Jules Verne, Amiens
 *     http://q4md-forcefieldtools.org
 *
*/

extern void     Mol2Write(FILE *fOut, UNIT uUnit, int choice);
extern UNIT     uTriposReadUnit(FILE *fIn);


