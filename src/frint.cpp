///////////////////////////////////////////////////////////////////////////////
// =========================================================================
// =========================================================================
// ===============FRACTIONAL INTEGRATION TOOLBOX============================
// =========================================================================
// =========================================================================
//  FRACTIONAL INTEGRATION TOOLBOX
//  RELEASE 1.0
//  JUNE 2013
//
//  A toolbox from the Santamaria Laboratory at the
//  University of Texas at San Antonio.
//  ( http://utsa.edu/Santamarialab/index.htm )
//  In collaboration
//  with the Computational Systems Biology Core/Computational Biology
//  Initiative( www.cbi.utsa.edu )
//
//  Grant Acknowledgments:
//  "This work received computational support from
//  Computational System Biology Core, funded by the National Institute
//  on Minority Health and Health Disparities (G12MD007591) from the
//  National Institutes of Health."
//
//  This work was supported by grants EF-1137897 and HDR-0932339 from the
//  National Science Foundation
//
//  Software Developer(s):
//  Toma.Marinov@utsa.edu
//  Nelson.Ramirez@utsa.edu
//
//  Project PI:
//  Fidel.Santamaria@utsa.edu
//
//  Notes:
//  Evaluate fractional integrals at half intervals.
//  Multi-core enabled via OpenMP, processing each signal independently.
//  First order corrections applied, 2nd order can be added easily.
//  Double and Single Precision enabled.
//  Interpolation only mode for GPU interfacing enabled.
//
//  For Matlab(mcc) to recognize C++ files, it must have the .cpp
//  extension.
//
//  This C++ code should compile from within the Matlab command line
//  as follows:
//   compileLinux.m
//   compileWindows.m
//
//  Key Interpolation Algorithms:
//  CUBIC SPLINE:
//  HERMITE SPLINE:
//
//  Additional references:
//  http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
//  http://en.wikipedia.org/wiki/Spline_interpolation
//  http://en.wikipedia.org/wiki/Monotone_cubic_interpolation
//  Both NC methods are described in "Algorithms for the Fractional Calculus:
//  A Selection of Numerical Methods", K. Diethlem, N.J. Ford, A.D. Dreed,
//  Yu. Luchko, Comp.Methods Appl. Mech. Engrg. 194(2005)
//
///////////////////////////////////////////////////////////////////////////////
//Fractional Integration Toolbox 1.0
//UTSA RESEARCH LICENSE (SOURCE CODE)
//The University of Texas at San Antonio has developed certain software and
//documentation that it desires to make available without charge to anyone for
//academic, research, experimental or personal use. This license is designed
//to guarantee freedom to use the software for these purposes. If you wish to
//distribute or make other use of the software, you may purchase a license to
//do so from the University of Texas.
//The accompanying source code is made available to you under the terms of
//this UT Research License (this "UTRL"). By clicking the "ACCEPT" button,
//or by installing or using the code, you are consenting to be bound by
//this UTRL. If you do not agree to the terms and conditions of this license,
//do not click the "ACCEPT" button, and do not install or use any part
//of the code.
//The terms and conditions in this UTRL not only apply to the source code
//made available by UT, but also to any improvements to, or derivative works
//of, that source code made by you and to any object code compiled from such
//source code, improvements or derivative works.
//1. DEFINITIONS.
//1.1 "Commercial Use" shall mean use of Software or Documentation by Licensee
//for direct or indirect financial, commercial or strategic gain or advantage,
//including without limitation: (a) bundling or integrating the Software with
//any hardware product or another software product for transfer, sale or license
//to a third party (even if distributing the Software on separate media and not
//charging for the Software); (b) providing customers with a link to the
//Software or a copy of the Software for use with hardware or another software
//product purchased by that customer; or (c) use in connection with the
//performance of services for which Licensee is compensated.
//1.2 "Derivative Products" means any improvements to, or other derivative
//works of, the Software made by Licensee.
//1.3 "Documentation" shall mean all manuals, user documentation, and other
//related materials pertaining to the Software that are made available to
//Licensee in connection with the Software.
//1.4 "Licensor" shall mean The University of Texas.
//1.5 "Licensee" shall mean the person or entity that has agreed to the
//terms hereof and is exercising rights granted hereunder.
//1.6 "Software" shall mean the computer program(s) referred to as
//"Fractional Integration Toolbox 1.0" made available under this UTRL in
//source code form, including any error corrections, bug fixes, patches,
//updates or other modifications that Licensor may in its sole discretion
//make available to Licensee from time to time, and any object code compiled
//from such source code.
//2. GRANT OF RIGHTS.
//Subject to the terms and conditions hereunder, Licensor hereby grants to
//Licensee a worldwide, non- transferable, non-exclusive license to (a) install,
//use and reproduce the Software for academic, research,
//experimental and personal use (but specifically excluding Commercial Use);
//(b) use and modify the Software to create Derivative Products, subject
//to Section 3.2; and (c) use the Documentation, if any, solely in connection
//with Licensee's authorized use of the Software.
//3. RESTRICTIONS; COVENANTS.
//3.1 Licensee may not: (a) distribute, sub-license or otherwise transfer copies
//or rights to the Software (or any portion thereof) or the Documentation;
//(b) use the Software (or any portion thereof) or Documentation for
//Commercial Use, or for any other use except as described in Section 2;
//(c) copy the Software or Documentation other than for archival and
//backup purposes; or (d) remove any product identification, copyright,
//proprietary notices or labels from the Software and Documentation.
//This UTRL confers no rights upon Licensee except those expressly granted herein.
//3.2 Licensee hereby agrees that it will provide a copy of all Derivative
//Products to Licensor and that its use of the Derivative Products will be
//subject to all of the same terms, conditions, restrictions and limitations
//on use imposed on the Software under this UTRL. Licensee hereby grants
//Licensor a worldwide, non- exclusive, royalty-free license to reproduce,
//prepare derivative works of, publicly display, publicly
//perform, sublicense and distribute Derivative Products. Licensee also hereby
//grants Licensor a worldwide, non-exclusive, royalty-free patent license to
//make, have made, use, offer to sell, sell, import and otherwise transfer the
//Derivative Products under those patent claims licensable by Licensee that
//are necessarily infringed by the Derivative Products.
//4. PROTECTION OF SOFTWARE.
//4.1 Confidentiality. The Software and Documentation are the confidential and
//proprietary information of Licensor. Licensee agrees to take adequate steps to
//protect the Software and Documentation from unauthorized disclosure or use.
//Licensee agrees that it will not disclose the Software or Documentation to
//any third party.
//4.2 Proprietary Notices. Licensee shall maintain and place on any copy of
//Software or Documentation that it reproduces for internal use all notices
//as are authorized and/or required hereunder. Licensee shall include a copy
//of this UTRL and the following notice, on each copy of the Software and
//Documentation. Such license and notice shall be embedded in each copy of
//the Software, in the video screen display, on the physical medium embodying
//the Software copy and on any Documentation:
//Copyright © 2013, The University of Texas at San Antonio. All rights reserved.
//UNIVERSITY EXPRESSLY DISCLAIMS ANY AND ALL WARRANTIES CONCERNING THIS SOFTWARE
//AND DOCUMENTATION, INCLUDING ANY WARRANTIES OF MERCHANTABILITY, FITNESS FOR
//ANY PARTICULAR PURPOSE, NON- INFRINGEMENT AND WARRANTIES OF PERFORMANCE, AND
//ANY WARRANTY THAT MIGHT OTHERWISE ARISE FROM COURSE OF DEALING OR USAGE OF
//TRADE. NO WARRANTY IS EITHER EXPRESS OR IMPLIED WITH RESPECT TO THE USE OF
//THE SOFTWARE OR DOCUMENTATION. Under no circumstances shall University be
//liable for incidental, special, indirect, direct or consequential damages
//or loss of profits, interruption of business, or related expenses which
//may arise from use of Software or Documentation, including but not limited
//to those resulting from defects in Software and/or Documentation, or loss
//or inaccuracy of data of any kind.
//5. WARRANTIES.
//5.1 Disclaimer of Warranties. TO THE EXTENT PERMITTED BY APPLICABLE LAW,
//THE SOFTWARE AND DOCUMENTATION ARE BEING PROVIDED ON AN "AS IS" BASIS WITHOUT
//ANY WARRANTIES OF ANY KIND RESPECTING THE SOFTWARE OR DOCUMENTATION, EITHER
//EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO ANY WARRANTY OF DESIGN,
//MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR NON-INFRINGEMENT.
//5.2 Limitation of Liability. UNDER NO CIRCUMSTANCES UNLESS REQUIRED BY
//APPLICABLE LAW SHALL LICENSOR BE LIABLE FOR INCIDENTAL, SPECIAL, INDIRECT,
//DIRECT OR CONSEQUENTIAL DAMAGES OR LOSS OF PROFITS, INTERRUPTION OF BUSINESS,
//OR RELATED EXPENSES WHICH MAY ARISE AS A RESULT OF THIS LICENSE OR OUT OF
//THE USE OR ATTEMPT OF USE OF SOFTWARE OR DOCUMENTATION INCLUDING BUT NOT
//LIMITED TO THOSE RESULTING FROM DEFECTS IN SOFTWARE AND/OR DOCUMENTATION,
//OR LOSS OR INACCURACY OF DATA OF ANY KIND. THE FOREGOING EXCLUSIONS AND
//LIMITATIONS WILL APPLY TO ALL CLAIMS AND ACTIONS OF ANY KIND, WHETHER
//BASED ON CONTRACT, TORT (INCLUDING, WITHOUT LIMITATION, NEGLIGENCE), OR
//ANY OTHER GROUNDS.
//6. INDEMNIFICATION.
//Licensee shall indemnify, defend and hold harmless Licensor, the
//University of Texas System, their Regents, and their officers,
//agents and employees from and against any claims, demands, or
//causes of action whatsoever caused by, or arising out of, or resulting
//from, the exercise or practice of the license granted hereunder by
//Licensee, its officers, employees, agents or representatives.
//7. TERMINATION.
//If Licensee breaches this UTRL, Licensee's right to use the Software
//and Documentation will terminate immediately without notice, but all
//provisions of this UTRL except Section 2 will survive termination and
//continue in effect. Upon termination, Licensee must destroy all copies
//of the Software and Documentation.
//8. GOVERNING LAW; JURISDICTION AND VENUE.
//The validity, interpretation, construction and performance of this UTRL
//shall be governed by the laws of the State of Texas. The Texas state
//courts of Travis County, Texas (or, if there is exclusive federal
//jurisdiction,the United States District Court for the Central
//District of Texas) shall have exclusive jurisdiction and venue
//over any dispute arising out of this UTRL, and Licensee consents
//to the jurisdiction of such courts. Application of the
//United Nations Convention on Contracts for the International
//Sale of Goods is expressly excluded.
//9. EXPORT CONTROLS.
//This license is subject to all applicable export restrictions. Licensee
//must comply with all export and import laws and restrictions and
//regulations of any United States or foreign agency or authority relating
//to the Software and its use.
//10. U.S. GOVERNMENT END-USERS.
//The Software is a "commercial item," as that term is defined in
//48 C.F.R. 2.101, consisting of "commercial computer software" and
//"commercial computer software documentation," as such terms are
//used in 48 C.F.R. 12.212 (Sept. 1995) and 48 C.F.R. 227.7202 (June 1995).
//Consistent with 48 C.F.R. 12.212, 48 C.F.R. 27.405(b)(2) (June 1998)
//and 48 C.F.R. 227.7202, all U.S. Government End Users acquire the
//Software with only those rights as set forth herein.
//11. MISCELLANEOUS
//If any provision hereof shall be held illegal, invalid or unenforceable,
//in whole or in part, such provision shall be modified to the minimum extent
//necessary to make it legal, valid and enforceable, and the legality, validity
//and enforceability of all other provisions of this UTRL shall not be affected
//thereby. Licensee may not assign this UTRL in whole or in part, without
//Licensor's prior written consent. Any attempt to assign this UTRL without
//such consent will be null and void. This UTRL is the complete and exclusive
//statement between Licensee and Licensor relating to the subject matter
//hereof and supersedes all prior oral and written and all contemporaneous
//oral negotiations, commitments and understandings of the parties, if any.
//Any waiver by either party of any default or breach hereunder shall
//not constitute a waiver of any provision of this UTRL or of any
//subsequent default or breach of the same or a different kind.
//END OF LICENSE
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <cassert>
/*=============================================================================
 *=============================================================================
 *===============COMPILE CONTROL===============================================
 *===============THESE CAN BE CONTROLLED AT THE COMPILE COMMAND================
 *=============================================================================
 *===========================================================================*/
//#define MAININTERFACE
//#define MATLABINTERFACE
//#define FRDEBUG
//#define OPENMP
//#define SINGLEPRECISION
//#define DOUBLEPRECISION
//#define INTERPOLATEONLY
/*=============================================================================
 *=============================================================================
 *===============DATATYPE======================================================
 *=============================================================================
 *===========================================================================*/
#ifdef SINGLEPRECISION
#define DATATYPE float
#endif

#ifdef DOUBLEPRECISION
#define DATATYPE double
#endif

/*=============================================================================
 *=============================================================================
 *===============OPENMP========================================================
 *=============================================================================
 *===========================================================================*/
#ifdef OPENMP
#include <omp.h>
#endif

using namespace std;

/*=============================================================================
 *=============================================================================
 *===============MATLAB INTERFACING============================================
 *=============================================================================
 *===========================================================================*/
#ifdef MATLABINTERFACE
#include "mex.h"
#endif

/* Input Arguments */
#define DATAMATRIX      prhs[0]
#define DELTAT        prhs[1]
#define ALPHA       prhs[2]
#define C1          prhs[3]
#define C2          prhs[4]
#define NUMCORES    prhs[5]

/* Output Arguments */
#define OUTPUT  plhs[0]


/*=============================================================================
 *=============================================================================
 *=============================HERMITE=========================================
 *=============================================================================
 *===========================================================================*/
#ifdef HERMITE
///////////////////////////////////////////////////////////////////////////
/// Novel Code:
/// HERMITE INTERPOLATION:
///
///
/// Algorithmic Reference:
/// 
/// [1] "Monotone Piecewise Cubic Interpolation" F.N. Fritsch, R.E. Carlson, SIAM J. Number. Anal.Vol 17, No.2, April 1980.
///
///http://epubs.siam.org/doi/abs/10.1137/0717021
///http://epubs.siam.org/doi/pdf/10.1137/0905021
///http://epubs.siam.org/doi/pdf/10.1137/0717021
///http://www.reading.ac.uk/web/FILES/maths/Tomos_Roberts.pdf
///http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19910011517_1991011517.pdf
///
/// 08/13/12, Toma.Marinov, Nelson Ramirez
///
/// Copyright 2012-2013.  UTSA. All Rights Reserved. 
///////////////////////////////////////////////////////////////////////////
int hermite(DATATYPE * dataMatrixPtr, int numRows, int numCols, DATATYPE* deltatPtr, DATATYPE* alphaPtr, DATATYPE* c1Ptr, DATATYPE* c2Ptr, DATATYPE* numCoresPtr,DATATYPE* outputPtr  ) {

   
    DATATYPE alpha = alphaPtr[0];
    DATATYPE c1 = c1Ptr[0];
    DATATYPE c2 = c2Ptr[0];
    DATATYPE numCores = numCoresPtr[0];
    int rows = numRows;  // maps to first input of hermite
    int cols = numCols;
    int N = numRows; // this is the time dimention
    vector<DATATYPE> weightsZeroCorrection(N,0); // N
    vector<DATATYPE> weightsFirstCorrection(2*N-1,0); //2N-1
    weightsZeroCorrection.at(0) = ((1+alpha)*pow(N,alpha))-pow(N,1+alpha)+pow(N-1,1+alpha);
    weightsZeroCorrection.at(N-1) = 1;
   
    for ( int i = 2; i < weightsZeroCorrection.size(); i++ ) {
        DATATYPE arg1 = pow(N-i+1,1+alpha);
        DATATYPE arg2 = 2*pow((N-i), 1+alpha) ;
        DATATYPE arg3 = pow(N-i-1,1+alpha);
        weightsZeroCorrection.at(i-1) =arg1-arg2+arg3;
    }
   
    int nFirst = weightsFirstCorrection.size();  /* the length of the first order correction vector */
    weightsFirstCorrection.at(0) = (1+alpha)*pow(nFirst,alpha) - pow(nFirst,1+alpha) + pow(nFirst-1,1+alpha);
    weightsFirstCorrection.at(nFirst-1) = 1;
   
    for ( int i = 2; i < weightsFirstCorrection.size(); i++ ) {
        DATATYPE arg1 = pow(nFirst-i+1,1+alpha);
        DATATYPE arg2 = 2*pow((nFirst-i), 1+alpha) ;
        DATATYPE arg3 = pow(nFirst-i-1,1+alpha);
        weightsFirstCorrection.at(i-1) =arg1-arg2+arg3;
    }
   
    ///////////////////////////////////////////////////////////////////////
    // Constants per signal
    ///////////////////////////////////////////////////////////////////////
    DATATYPE startT = 0; // should always be zero, comment out???
    DATATYPE deltaT = deltatPtr[0];   // ??? Notes: pass in only T0,TN=N*deltaT, N or T0, deltaT, N
    DATATYPE oneOverDt = 1/deltaT;
    DATATYPE oneOverTwoDt = 1/(2*deltaT);
    DATATYPE oneEighthDt = 0.125*deltaT;	
    int numIntervals = rows -1;
    int twoTimesNumIntervals = 2*numIntervals;
    int rowsMinusTwo = rows -2;
    ///////////////////////////////////////////////////////////////////////
    // column by column processing of the signals
    // a signal is defined at each column
    ///////////////////////////////////////////////////////////////////////
#ifdef OPENMP
int numCoresInt = (int)numCores;  // convert from float/double to integer
omp_set_num_threads(numCoresInt);  // requires integer
#endif

///////////////////////////////////////////////////////////////////////////
// Loop through all input signals.
// Signals should be in columns.
///////////////////////////////////////////////////////////////////////////
#pragma omp parallel for
    for ( int i = 0; i < cols; i++ ) {
#ifdef FRDEBUG
#ifdef OPENMP
uint nThreads, tid;
tid = omp_get_thread_num();
nThreads = omp_get_num_threads();
cout<<"tid="<<tid<<"nThreads="<<nThreads<<endl;
#endif
#endif

        DATATYPE * y = (DATATYPE*)&(dataMatrixPtr[i*rows]); /* C++ to C integration */
      
        ///////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////
        // Generate Spline
        //
        // If we assume that delta K's that are equal to 0 are a small
        // fraction of the cases, having an efficient way to avoid
        // doing this within the main for loop is important.
        ///////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////
        vector<double> dely(numIntervals,0);
        vector<double> d(rows,0);
 
        ///////////////////////////////////////////////////////////////////
        //
        // Calculate dely and d
        //
        ///////////////////////////////////////////////////////////////////
        // starting values
        dely[0]= (y[1]-y[0])*oneOverDt; // calculate the first secant
        d[0]= dely[0]; // calculate the first derivative
        // both dely and d have the same upper limit
     
        // while parallelizable, avoid...
        for ( int j = 1; j <= rowsMinusTwo; j++ ){
            dely[j]=(y[j+1]-y[j])*oneOverDt;  // secant
            d[j]= (y[j+1]-y[j-1])*oneOverTwoDt;  // derivative, updated on v0.0.2
        }
        // ending values
        d[numIntervals] = dely[rowsMinusTwo];
       
        ///////////////////////////////////////////////////////////////////
        //
        // Go from k = 0 to k = rows - 2
        //
        // ??? Perhaps in next release, investigate OpenMP 
	    // ??? parallelism of the loop.  There is a potential race
	    // ??? condition that needs to be addressed: k and k+1
        //
        ///////////////////////////////////////////////////////////////////
        for ( int k = 0 ; k <= rowsMinusTwo; k++ ) {
            if ( dely[k] < 1.0e-30 ) {
                d[k] = 0;
                d[k+1] = 0;  // invalidates parallelism at k index level
            }
            else {
                double alpha = d[k]/dely[k];
                double beta = d[k+1]/dely[k];
                bool boolA = (alpha < 0 );
                bool boolB = (beta < 0 );
                if ( boolA || boolB ) {
                    d[k] = 0;
                    d[k+1] = 0;
                }
                else {
                   
                    double alphabeta = alpha*alpha + beta*beta;
                   
                    ///////////////////////////////////////////////////////
                    // Monotonicity check: This varation produces
                    // best visual results.
                    // L2:
                    ///////////////////////////////////////////////////////
                    if ( alphabeta > 9  )  {
                        double tau = 3/pow(alphabeta,0.5);
                        d[k] = alpha * tau * dely[k];
                        d[k+1] = beta* tau* dely[k];
                       
                    }
                    ///////////////////////////////////////////////////////
                    // Monotonicity check: This variation does
                    // not produce best visual results.
                    ///////////////////////////////////////////////////////
                    //if ( alpha > 3 ) {
                    //  d[k] = 3* dely[k];
                    //} // end alpha > 3 if check
                    //if ( beta > 3 ) {
                    //  d[k+1] = 3*dely[k+1];
                    //} // end beta > 3 if check
                } // end if block
            } // end else block
        }  // end for loop

        ///////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////
        // End Generate Spline
        ///////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////

     
       
        // Handle special cases( last element )
	// We need zeroTotal and oneTotal to be local thread variables.
	// Do not move out of the loop.
        DATATYPE zeroTotal = y[numIntervals]*weightsZeroCorrection[numIntervals];
        DATATYPE oneTotal = y[numIntervals]*weightsFirstCorrection[twoTimesNumIntervals];

	// do not parallelize, however, perhaps later investigate
	// effect of loop unrolling on performance.
        for ( int k = 0; k < numIntervals; k++ ) {
            zeroTotal = zeroTotal+(y[k]*weightsZeroCorrection[k]);
            DATATYPE yInt = (y[k]+y[k+1])*0.5  + (d[k]-d[k+1])*oneEighthDt;
            oneTotal=oneTotal+((y[k]*weightsFirstCorrection[2*k])) + ((yInt)*weightsFirstCorrection[2*k+1]);   // do not parallelize
        }
        ///////////////////////////////////////////////////////////////////
        // CALCULATE RESULT FOR THIS SIGNAL
        ///////////////////////////////////////////////////////////////////
        DATATYPE result = (c1 * zeroTotal) + (c2 * oneTotal);
        outputPtr[i] = result;

 }  // end loop through each signal( there is 1 signal per matlab column )
return 0;
} // end function

///////////////////////////////////////////////////////////////////////////
/// Novel Code:
/// HERMITE INTERPOLATION:
///
///
/// Algorithmic Reference:
/// 
/// [1] "Monotone Piecewise Cubic Interpolation" F.N. Fritsch, R.E. Carlson, SIAM J. Number. Anal.Vol 17, No.2, April 1980.
///
///http://epubs.siam.org/doi/abs/10.1137/0717021
///http://epubs.siam.org/doi/pdf/10.1137/0905021
///http://epubs.siam.org/doi/pdf/10.1137/0717021
///http://www.reading.ac.uk/web/FILES/maths/Tomos_Roberts.pdf
///http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19910011517_1991011517.pdf
///
/// 08/13/12-08/22/12, Toma.Marinov, Nelson Ramirez
///
/// Copyright 2012.  UTSA. All Rights Reserved. 
///////////////////////////////////////////////////////////////////////////
int hermiteinterpolateonly(DATATYPE * dataMatrixPtr, int numRows, int numCols, DATATYPE* deltatPtr, DATATYPE* alphaPtr, DATATYPE* c1Ptr, DATATYPE* c2Ptr, DATATYPE* numCoresPtr,DATATYPE* outputPtr  ) {

   
    DATATYPE alpha = alphaPtr[0];
    DATATYPE c1 = c1Ptr[0];
    DATATYPE c2 = c2Ptr[0];
    DATATYPE numCores = numCoresPtr[0];
    int rows = numRows;  // maps to first input of hermite
    int cols = numCols;
    int N = numRows; // this is the time dimention
  
    int totalInterpolatedRowCount = (2*rows)-1;
    DATATYPE startT = 0; // should always be 0, comment out???
    DATATYPE deltaT = deltatPtr[0];
    DATATYPE oneOverDt = 1/deltaT;
    DATATYPE oneOverTwoDt = 1/(2*deltaT);
    DATATYPE oneEighthDt = 0.125*deltaT;		
    int numIntervals = rows -1;
    int twoTimesNumIntervals = 2*numIntervals;
    int rowsMinusTwo = rows -2;
    int twoTimesRowsMinusOne = 2*rows -1;
   
 
    ///////////////////////////////////////////////////////////////////////
    // column by column processing of the signals
    // a signal is defined at each column
    ///////////////////////////////////////////////////////////////////////
#ifdef OPENMP
int numCoresInt = (int)numCores;  // convert from float/double to integer
omp_set_num_threads(numCoresInt);  // requires integer
#endif

///////////////////////////////////////////////////////////////////////////
// Loop through all input signals.
// Signals should be in columns.
///////////////////////////////////////////////////////////////////////////
#pragma omp parallel for
    for ( int i = 0; i < cols; i++ ) {
#ifdef FRDEBUG
#ifdef OPENMP
uint nThreads, tid;
tid = omp_get_thread_num();
nThreads = omp_get_num_threads();
cout<<"tid="<<tid<<"nThreads="<<nThreads<<endl;
#endif
#endif

        DATATYPE * y = (DATATYPE*)&(dataMatrixPtr[i*rows]); /* C++ to C integration */
      
        ///////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////
        // Generate Spline
        //
        // If we assume that delta K's that are equal to 0 are a small
        // fraction of the cases, having an efficient way to avoid
        // doing this within the main for loop is important.
        ///////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////
        vector<double> dely(numIntervals,0);
        vector<double> d(rows,0);
 
        ///////////////////////////////////////////////////////////////////
        //
        // Calculate dely and d
        //
        ///////////////////////////////////////////////////////////////////
        // starting values
        dely[0]= (y[1]-y[0])*oneOverDt; // calculate the first secant
        d[0]= dely[0]; // calculate the first derivative
        // both dely and d have the same upper limit
     
        // while parallelizable, avoid...
        for ( int j = 1; j <= rows-2; j++ ){
            dely[j]=(y[j+1]-y[j])*oneOverDt;  // secant
            d[j]= (y[j+1]-y[j-1])*oneOverTwoDt;  // derivative
        }
        // ending values
        d[rows-1] = dely[rows-2];
       
        ///////////////////////////////////////////////////////////////////
        //
        // Go from k = 0 to k = rows - 2
        //
        //
        //
        ///////////////////////////////////////////////////////////////////
        for ( int k = 0 ; k <= rowsMinusTwo; k++ ) {
            if ( dely[k] < 1.0e-30 ) {
                d[k] = 0;
                d[k+1] = 0;  // invalidates parallelism at k index level
            }
            else {
                double alpha = d[k]/dely[k];
                double beta = d[k+1]/dely[k];
                bool boolA = (alpha < 0 );
                bool boolB = (beta < 0 );
                if ( boolA || boolB ) {
                    d[k] = 0;
                    d[k+1] = 0;
                }
                else {
                   
                    double alphabeta = alpha*alpha + beta*beta;
                   
                    ///////////////////////////////////////////////////////
                    // Monotonicity check: This varation produces
                    // best visual results.
                    // L2:
                    ///////////////////////////////////////////////////////
                    if ( alphabeta > 9  )  {
                        double tau = 3/pow(alphabeta,0.5);
                        d[k] = alpha * tau * dely[k];
                        d[k+1] = beta* tau* dely[k];
                       
                    }
                    ///////////////////////////////////////////////////////
                    // Monotonicity check: This variation does
                    // not produce best visual results.
                    ///////////////////////////////////////////////////////
                    //if ( alpha > 3 ) {
                    //  d[k] = 3* dely[k];
                    //} // end alpha > 3 if check
                    //if ( beta > 3 ) {
                    //  d[k+1] = 3*dely[k+1];
                    //} // end beta > 3 if check
                   
                } // end if block
            } // end else block
        }  // end for loop

        ///////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////
        // End Generate Spline
        ///////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////
       

        // get pointer to ith signal original data
        DATATYPE * dataPtr = (DATATYPE*)&(dataMatrixPtr[i*rows]); /* C++ to C integration */
        // get pointer to ith output signal
        DATATYPE * dataOutputPtr = (DATATYPE*)&(outputPtr[i*(twoTimesRowsMinusOne)]); /* C++ to C integration */
       
        int counter =0;
        for (int k = 1; k < twoTimesNumIntervals; k=k+2 ) {
            dataOutputPtr[k-1] = dataPtr[counter]; 
            DATATYPE yInt = (y[counter]+y[counter+1])*0.5  + (d[counter]-d[counter+1])*oneEighthDt;
            dataOutputPtr[k] = yInt;
            counter++;
        }
        
        dataOutputPtr[twoTimesNumIntervals] = dataPtr[numIntervals];
       
     
 }  // end loop through each signal( there is 1 signal per matlab column )
return 0;
} // end function



#endif

/*=============================================================================
 *=============================================================================
 *=============================CUBIC===========================================
 *=============================================================================
 *===========================================================================*/
#ifdef CUBIC

///////////////////////////////////////////////////////////////////////////
/// Novel Code:
/// Cubic Interpolation  with Tridiagonal Solver:
///
/// All input vectors must be pre-allocated
/// All input vectors must be the same length
/// lowerDiagonal and upperDiagonal must be the same length as
/// centralDiagonal must be the same length as d
/// centralDiagonal must be the same length as result
/// Algorithmic Reference:
/// 
/// Reference:
/// [1] Carl de Boor, A Practical Guide to Splines, Springer Verlag, 1978.
/// [2] Conte, S.D., and deBoor, C. (1972). Elementary Numerical Analysis. McGraw-Hill, New York.
/// [3] http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_%28Thomas_algorithm%29
///
///
/// 08/17/12-8/22/12, Toma.Marinov, Nelson Ramirez
///
/// Copyright 2012.  UTSA. All Rights Reserved. 
///////////////////////////////////////////////////////////////////////////
int cubic(DATATYPE * dataMatrixPtr, int numRows, int numCols, DATATYPE* deltatPtr, DATATYPE* alphaPtr, DATATYPE* c1Ptr, DATATYPE* c2Ptr, DATATYPE* numCoresPtr,DATATYPE* outputPtr ) {

    DATATYPE alpha = alphaPtr[0];
    DATATYPE c1 = c1Ptr[0];
    DATATYPE c2 = c2Ptr[0];
    DATATYPE numCores = numCoresPtr[0];
    int rows = numRows;
    int cols = numCols;
    int N = numRows; // this is the time dimention
    int numIntervals = rows -1;
    int rowsMinusTwo = rows -2;
    int twoTimesNumIntervals = 2*numIntervals;

    // The zero order correction is always calculated
    vector<DATATYPE> weightsZeroCorrection(N,0); // N
    vector<DATATYPE> weightsFirstCorrection(2*N-1,0); //2N-1
  
    // Leading weights
    weightsZeroCorrection.at(0) = ((1+alpha)*pow(N,alpha))-pow(N,1+alpha)+pow(N-1,1+alpha);
    weightsZeroCorrection.at(N-1) = 1;
    for ( int i = 2; i < weightsZeroCorrection.size(); i++ ) {
        DATATYPE arg1 = pow(N-i+1,1+alpha);
        DATATYPE arg2 = 2*pow((N-i), 1+alpha) ;
        DATATYPE arg3 = pow(N-i-1,1+alpha);
        weightsZeroCorrection.at(i-1) =arg1-arg2+arg3;
    }
   
    // First order weights
    int nFirst = weightsFirstCorrection.size();  // the length of the first order correction vector
    weightsFirstCorrection.at(0) = (1+alpha)*pow(nFirst,alpha) - pow(nFirst,1+alpha) + pow(nFirst-1,1+alpha);
    weightsFirstCorrection.at(nFirst-1) = 1;
   
    for ( int i = 2; i < weightsFirstCorrection.size(); i++ ) {
        DATATYPE arg1 = pow(nFirst-i+1,1+alpha);
        DATATYPE arg2 = 2*pow((nFirst-i), 1+alpha) ;
        DATATYPE arg3 = pow(nFirst-i-1,1+alpha);
        weightsFirstCorrection.at(i-1) =arg1-arg2+arg3;
    }
   
    // Time and deltaT information
    DATATYPE startT = 0;  // this should always be zero, do we need this???
    DATATYPE deltaT = deltatPtr[0];
    DATATYPE twoDeltaT = 2*deltaT;
    DATATYPE oneOverDeltaT = 1/deltaT;
    DATATYPE threeTimesOneOverDeltaTSquared = 3/(deltaT*deltaT);
    DATATYPE deltaTSquaredOverSix = (deltaT*deltaT)/6;
    DATATYPE deltaTSquaredOverEight = (deltaT*deltaT)/8;
    DATATYPE deltaTSquaredOverFortyEight = (deltaT*deltaT)/48;	
    ///////////////////////////////////////////////////////////////////////
    // column by column processing of the signals
    // a signal is defined at each column
    ///////////////////////////////////////////////////////////////////////
#ifdef OPENMP
int numCoresInt = (int)numCores;  // convert from float/DATATYPE to integer
omp_set_num_threads(numCoresInt);  // requires integer
#endif

///////////////////////////////////////////////////////////////////////////
// Loop through all input signals.
// Signals should be in columns.
///////////////////////////////////////////////////////////////////////////
#pragma omp parallel for
    for ( int i = 0; i < cols; i++ ) {
#ifdef FRDEBUG
#ifdef OPENMP
uint nThreads, tid;
tid = omp_get_thread_num();
nThreads = omp_get_num_threads();
cout<<"tid="<<tid<<"nThreads="<<nThreads<<endl;
#endif
#endif

        DATATYPE * y = (DATATYPE*)&(dataMatrixPtr[i*rows]); // C++ to C integration

        ///////////////////////////////////////////////////////////////////
        // Create the cubic spline coefficients
        ///////////////////////////////////////////////////////////////////
        vector<double> secondDerivatives(rows,0);  // maps to x in the linear system
       
        ///////////////////////////////////////////////////////////////////
        // Bulirsh, Stoer. Introduction to numerical analysis(2nd ed, Springer, 1993), pg. 98-99
        // Non-periodic case.
        // Case C:
        // lambda = upper diagonal
        //     h1/ hn+h1 = deltaT / ( 2*deltaT ) = 0.5
        // mu = lower diagonal
        //     1-lambda = hn / hn+h1 = 0.5
        // 2 = central diagonal
        ///////////////////////////////////////////////////////////////////
        vector<double> lowerDiagonal(rows,0.5);   // maps to mu in the linear system ( of everthing but last is 0.5 )
        vector<double> centralDiagonal(rows,2); // maps to b in the linear system
        vector<double> upperDiagonal(rows,0.5);   // maps to lambda in the linear system ( of anything but 0 is 0.5)
        vector<double> d(rows,0); // maps to d in the linear system
       
        ///////////////////////////////////////////////////////////////////
        // Setting up the boundary conditions
        ///////////////////////////////////////////////////////////////////
        d[0] = 0; // d0= 6/h1( ( (y1-y0) /h1) -y0' )
        upperDiagonal[0] = 0.0; // lambda at 0 is 0
       
        // This loop goes from [1 to rows-2]
        for ( int jj = 1; jj < numIntervals; jj++ ){
            ///////////////////////////////////////////////////////////////
            // Evaluate the d = 6/(hj+hj+1) { (yj+1-yj)/(hj+1) - (yj - yj-1)/(hj) }
            //  = 6/( 2* (deltaT*deltaT) )
            //  = (3/(deltaT*deltaT)) * ( y[jj+1] - 2*y[jj] + y[jj-1] )
            ///////////////////////////////////////////////////////////////
            d[jj] = threeTimesOneOverDeltaTSquared * ( y[jj+1] - 2*y[jj] + y[jj-1] ) ;
        }
        d[numIntervals] = 0.0;  // dn = 6/hn*(yn'- (yn-yn-1)/hn)
        lowerDiagonal[numIntervals] = 0.0; // mu at rows -1 = 0
      
       
        ///////////////////////////////////////////////////////////////////
        // Calculate the spline coefficients
        // and solve for the second derivative using
        // [1] Bulirsh, Stoer. Introduction to numerical analysis(2nd ed, Springer, 1993), pg. 102
        // [2] Wikipedia: TDMA algorithm
        // [3] ( deBoor, 1972 )
        //
        ///////////////////////////////////////////////////////////////////
        for ( int ij = 1; ij < rows; ij++ ) {
            double temp = lowerDiagonal[ij] / centralDiagonal[ij-1];
            centralDiagonal[ij] = centralDiagonal[ij] - temp* upperDiagonal[ij-1];
            d[ij] = d[ij] - temp * d[ij-1];
        }
        // The last element
        secondDerivatives[numIntervals] = d[numIntervals] / centralDiagonal[rows-1];
        for ( int ij = rowsMinusTwo; ij >= 0; ij-- ) {
            secondDerivatives[ij] = (d[ij]- upperDiagonal[ij]*secondDerivatives[ij+1])/centralDiagonal[ij];
        }
       
        ///////////////////////////////////////////////////////////////////
        // Evaluate fractional integral
        ///////////////////////////////////////////////////////////////////
       
        // Handle special cases( last element )
        DATATYPE zeroTotal = y[numIntervals]*weightsZeroCorrection[numIntervals];
        DATATYPE oneTotal = y[numIntervals]*weightsFirstCorrection[twoTimesNumIntervals];
        // this loop goes from j = [0 to row-2]
        for ( int j = 0; j < numIntervals; j++ ) {
             zeroTotal = zeroTotal+(y[j]*weightsZeroCorrection[j]);
             //////////////////////////////////////////////////////////////
             ////////////////////EVALUATE SPLINE///////////////////////////
             //////////////////////////////////////////////////////////////
            
             //////////////////////////////////////////////////////////////
             // We still need to cleanup this expression
             //
             // In order to have the diagonal elements always equal to 2
             // The interpolation equation becomes the following:
             //
             // Bulirsh, Stoer. Introduction to numerical analysis(2nd ed, Springer, 1993), pg. 110
             // yInt = ....
             // alpha[j] = y[j]
             // gamma[j] = Mj /2;
             // beta[j] = (y[j+1] - y[j]) / hj+1  - (2*Mj +Mj+1) /6 * hj+1
             // delta = (Mj+1-Mj) / ( 6*hj+1)
             //
             //
             // yInt = alpha(j) + Beta(j)/
             // 5 terms
             // 1)
             //
             //////////////////////////////////////////////////////////////
             DATATYPE yInt = 0;
             yInt =  y[j]
                     + 0.5*( y[j+1] - y[j] )
                     - ( secondDerivatives[j]  + 0.5*secondDerivatives[j+1]) *deltaTSquaredOverSix
                     + (deltaTSquaredOverEight * secondDerivatives[j])
                     + deltaTSquaredOverFortyEight*(secondDerivatives[j+1] - secondDerivatives[j] ) ;

            
             //cout<<"yInt = "<<yInt<<endl;
             //////////////////////////////////////////////////////////////
             //////////////////////////////////////////////////////////////
             //////////////////////////////////////////////////////////////
            oneTotal=oneTotal+((y[j]*weightsFirstCorrection[2*j])) + ((yInt)*weightsFirstCorrection[2*j+1]);
 
        } // end loop through all the intervals
        ///////////////////////////////////////////////////////////////////
        // Key update is to avoid searching for an interval, since
        // our processing is already at an interval level, there
        // is no need to
        ///////////////////////////////////////////////////////////////////
       
        ///////////////////////////////////////////////////////////////////
        // CALCULATE RESULT FOR THIS SIGNAL
        ///////////////////////////////////////////////////////////////////
        DATATYPE result = (c1 * zeroTotal) + (c2 * oneTotal);
 
        outputPtr[i] = result;

 }  // end loop through each signal( there is 1 signal per matlab column )

    return 0;
}  // end  function

///////////////////////////////////////////////////////////////////////////
/// Novel Code:
/// Cubic interpolation only.
/// Cubic Interpolation with Tridiagonal Solver:
///
/// All input vectors must be pre-allocated
/// All input vectors must be the same length
/// lowerDiagonal and upperDiagonal must be the same length as
/// centralDiagonal must be the same length as d
/// centralDiagonal must be the same length as result
/// Algorithmic Reference:
/// 
/// Reference:
/// [1] Carl de Boor, A Practical Guide to Splines, Springer Verlag, 1978.
/// [2] Conte, S.D., and deBoor, C. (1972). Elementary Numerical Analysis. McGraw-Hill, New York.
/// [3] http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_%28Thomas_algorithm%29
///
///
/// 08/17/12-8/22/12, Toma.Marinov, Nelson Ramirez
///
/// Copyright 2012.  UTSA. All Rights Reserved. 
///////////////////////////////////////////////////////////////////////////
int cubicinterpolateonly(DATATYPE * dataMatrixPtr, int numRows, int numCols, DATATYPE* deltatPtr, DATATYPE* alphaPtr, DATATYPE* c1Ptr, DATATYPE* c2Ptr, DATATYPE* numCoresPtr,DATATYPE* outputPtr ) {


    DATATYPE alpha = alphaPtr[0];
    DATATYPE c1 = c1Ptr[0];
    DATATYPE c2 = c2Ptr[0];
    DATATYPE numCores = numCoresPtr[0];
    int rows = numRows;
    int cols = numCols;
    int N = numRows; // this is the time dimention
    int numIntervals = rows -1;
    int twoTimesRowsMinusOne = 2*rows-1;
    int twoTimesNumIntervals = 2*numIntervals;
    int rowsMinusTwo = rows -2;
    DATATYPE deltaT = deltatPtr[0];
    DATATYPE dtOver2 = deltaT/2;
    DATATYPE twoDeltaT = 2*deltaT;
    DATATYPE threeTimesOneOverDeltaTSquared = 3/(deltaT*deltaT);
    DATATYPE deltaTSquaredOverSix = (deltaT*deltaT)/6;
    DATATYPE deltaTSquaredOverEight = (deltaT*deltaT)/8;
    DATATYPE deltaTSquaredOverFortyEight = (deltaT*deltaT)/48;	

    ///////////////////////////////////////////////////////////////////////
    // column by column processing of the signals
    // a signal is defined at each column
    ///////////////////////////////////////////////////////////////////////
#ifdef OPENMP
int numCoresInt = (int)numCores;  // convert from float/DATATYPE to integer
omp_set_num_threads(numCoresInt);  // requires integer
#endif

///////////////////////////////////////////////////////////////////////////
// Loop through all input signals.
// Signals should be in columns.
///////////////////////////////////////////////////////////////////////////
#pragma omp parallel for
    for ( int i = 0; i < cols; i++ ) {
#ifdef FRDEBUG
#ifdef OPENMP
uint nThreads, tid;
tid = omp_get_thread_num();
nThreads = omp_get_num_threads();
cout<<"tid="<<tid<<"nThreads="<<nThreads<<endl;
#endif
#endif
        DATATYPE * y = (DATATYPE*)&(dataMatrixPtr[i*rows]); /* C++ to C integration */

        vector<double> secondDerivatives(rows,0);  // maps to x in the linear system
       
        ///////////////////////////////////////////////////////////////////
        // Bulirsh, Stoer. Introduction to numerical analysis(2nd ed, Springer, 1993), pg. 98-99
        // Non-periodic case.
        // Case C:
        // lambda = upper diagonal
        //     h1/ hn+h1 = deltaT / ( 2*deltaT ) = 0.5
        // mu = lower diagonal
        //     1-lambda = hn / hn+h1 = 0.5
        // 2 = central diagonal
        ///////////////////////////////////////////////////////////////////
        vector<double> lowerDiagonal(rows,0.5);   // maps to mu in the linear system ( of everthing but last is 0.5 )
        vector<double> centralDiagonal(rows,2); // maps to b in the linear system
        vector<double> upperDiagonal(rows,0.5);   // maps to lambda in the linear system ( of anything but 0 is 0.5)
        vector<double> d(rows,0); // maps to d in the linear system
       
        ///////////////////////////////////////////////////////////////////
        // Setting up the boundary conditions
        ///////////////////////////////////////////////////////////////////
        d[0] = 0; // d0= 6/h1( ( (y1-y0) /h1) -y0' )
        upperDiagonal[0] = 0.0; // lambda at 0 is 0
       
        // This loop goes from [1 to rows-2]
        for ( int jj = 1; jj < numIntervals; jj++ ){
            ///////////////////////////////////////////////////////////////
            // Evaluate the d = 6/(hj+hj+1) { (yj+1-yj)/(hj+1) - (yj - yj-1)/(hj) }
            //  = 6/( 2* (deltaT*deltaT) )
            //  = (3/(deltaT*deltaT)) * ( y[jj+1] - 2*y[jj] + y[jj-1] )
            ///////////////////////////////////////////////////////////////
            d[jj] = threeTimesOneOverDeltaTSquared * ( y[jj+1] - 2*y[jj] + y[jj-1] ) ;
        }
        d[numIntervals] = 0.0;  // dn = 6/hn*(yn'- (yn-yn-1)/hn)
        lowerDiagonal[numIntervals] = 0.0; // mu at rows -1 = 0
      
       
        ///////////////////////////////////////////////////////////////////
        // Calculate the spline coefficients
        // and solve for the second derivative using
        // [1] Bulirsh, Stoer. Introduction to numerical analysis(2nd ed, Springer, 1993), pg. 102
        // [2] Wikipedia: TDMA algorithm
        // [3] ( deBoor, 1972 )
        //
        ///////////////////////////////////////////////////////////////////
        for ( int ij = 1; ij < rows; ij++ ) {
            double temp = lowerDiagonal[ij] / centralDiagonal[ij-1];
            centralDiagonal[ij] = centralDiagonal[ij] - temp* upperDiagonal[ij-1];
            d[ij] = d[ij] - temp * d[ij-1];
        }
        // The last element
        secondDerivatives[numIntervals] = d[numIntervals] / centralDiagonal[numIntervals];
        for ( int ij = rowsMinusTwo; ij >= 0; ij-- ) {
            secondDerivatives[ij] = (d[ij]- upperDiagonal[ij]*secondDerivatives[ij+1])/centralDiagonal[ij];
        }
       
        ///////////////////////////////////////////////////////////////////
        // Evaluate fractional integral
        ///////////////////////////////////////////////////////////////////
      
        // Handle special cases( last element )
        
        // get pointer to ith signal original data
        DATATYPE * dataPtr = (DATATYPE*)&(dataMatrixPtr[i*rows]); /* C++ to C integration */
        // get pointer to ith output signal
        DATATYPE * dataOutputPtr = (DATATYPE*)&(outputPtr[i*(twoTimesRowsMinusOne)]); /* C++ to C integration */
  
        int counter =0;
        for (int k = 1; k <twoTimesNumIntervals; k=k+2 ) {
            dataOutputPtr[k-1] = dataPtr[counter]; 
            
            DATATYPE yInt = y[counter]
                     + 0.5*( y[counter+1] - y[counter] )
                     - ( secondDerivatives[counter]  + 0.5*secondDerivatives[counter+1]) *deltaTSquaredOverSix
                     + (deltaTSquaredOverEight * secondDerivatives[counter])
                     + deltaTSquaredOverFortyEight*(secondDerivatives[counter+1] - secondDerivatives[counter] ) ;



            dataOutputPtr[k] = yInt;
            counter++;
        }
        
        dataOutputPtr[twoTimesNumIntervals] = dataPtr[numIntervals];
        
        

        // make sure to release all dynamically allocated memory
        
       
 }  // end loop through each signal( there is 1 signal per matlab column )


    return 0;
}  // end function


#endif



#ifdef MATLABINTERFACE
/*=============================================================================
 *=============================================================================
 *===============MATLAB INTERFACING============================================
 *=============================================================================
 *===========================================================================*/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray*prhs[] )
    
{
    ///////////////////////////////////////////////////////////////////////
    // It is important to note, that Matlab stores the
    // data in column-major format.
    ///////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////
    // These variables hold the size of each input matrix
    ///////////////////////////////////////////////////////////////////////
    mwSize m,n;
    mwSize mT, nT;
    mwSize mAlpha, nAlpha;
    mwSize mc1,nc1;
    mwSize mc2,nc2;
    mwSize mNumCores, nNumCores;
   
    mwSize mOutput, nOutput;
   
    // Fractional Integration Input Parameters
    DATATYPE * dataMatrixPtr;
    DATATYPE * deltatPtr;
    DATATYPE * alphaPtr;
    DATATYPE * c1Ptr;
    DATATYPE * c2Ptr;
    DATATYPE * numCoresPtr;
 
    // Fractional Integration Output Parameters
    DATATYPE * outputPtr;
   

    // Check for proper number of arguments
    if (nrhs != 6) {
       
           
    ///////////////////////////////////////////////////////////////////////
    // SPLINE CASES
    ///////////////////////////////////////////////////////////////////////     
    mexErrMsgTxt("Six input arguments required.");      
    return;
   
    } else if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }

    // FRACTIONAL INTEGRATION
    // Parameter # 1, DataMatrix
    // The DataMatrix must have the "signals" along the columns...
    // For performance reasons.
    m = mxGetM(DATAMATRIX); // this should map to the time axis
    n = mxGetN(DATAMATRIX); // this is the spatial dimension
   
#ifdef SINGLEPRECISION
    if (!mxIsSingle(DATAMATRIX) || mxIsComplex(DATAMATRIX ) ) {
        mexErrMsgTxt("FRACTINT requires that DATAMATRIX be a non-complex mxn matrix of single datatype.");
    }
#endif
   

#ifdef DOUBLEPRECISION
    if (!mxIsDouble(DATAMATRIX) || mxIsComplex(DATAMATRIX ) ) {
        mexErrMsgTxt("FRACTINT requires that DATAMATRIX be a non-complex mxn matrix of double datatype.");
    }
#endif
   
  
#ifdef FRDEBUG
    cout<<"m="<<m<<endl;
    cout<<"n="<<n<<endl;
#endif

    dataMatrixPtr = (DATATYPE*)mxGetPr(DATAMATRIX);

#ifdef FRDEBUG
    printf("dataMatrixPtr= %x\n",dataMatrixPtr);
    printf("dataMatrixPtr= %f\n",dataMatrixPtr[0]);

    int offsetCount = 0;
    // i maps the rows
    for ( int i = 0; i < m; i++ ) {
        for ( int j = 0; j < n; j++ ) {
            printf("dataMatrixPtr= %x\n",&(dataMatrixPtr[offsetCount]));
            printf("dataMatrixPtr= %f\n",dataMatrixPtr[offsetCount]);
            offsetCount=offsetCount+1;
        }
    }
#endif
    // delta t 
    mT = mxGetM(DELTAT); // this gets number of rows
    nT = mxGetN(DELTAT); // this gets the number of columns
    deltatPtr = (DATATYPE*)mxGetPr(DELTAT);
#ifdef FRDEBUG
    printf("deltaT = %u\n",&(deltatPtr[0]));
    printf("deltatPtr = %f\n",(deltatPtr[0]));
#endif
    
 

    // alpha
    mAlpha = mxGetM(ALPHA);
    nAlpha = mxGetN(ALPHA);
    alphaPtr = (DATATYPE*)mxGetPr(ALPHA);

#ifdef FRDEBUG
    printf("alpha = %u\n",&(alphaPtr[0]));
    printf("alpha = %f\n",(alphaPtr[0]));
#endif
    ///////////////////////////////////////////////////////////////////////
    // we need to pass in coeffOne (C1)
    // (-1/3) * ( (dt^alpha) / ( gamma(2+alpha))  )
    //
    // we need to pass in coeffTwo (C2)
    // (-4/(2^alpha))*C1
    ///////////////////////////////////////////////////////////////////////
    mc1 = mxGetM(C1);
    nc1 = mxGetN(C1);
    c1Ptr = (DATATYPE*)mxGetPr(C1);

#ifdef FRDEBUG
    printf("c1 = %u\n",&(c1Ptr[0]));
    printf("c1 = %f\n",(c1Ptr[0]));
#endif

    mc2 = mxGetM(C2);
    nc2 = mxGetN(C2);
    c2Ptr = (DATATYPE*)mxGetPr(C2);

#ifdef FRDEBUG
    printf("c2 = %u\n",&(c2Ptr[0]));
    printf("c2 = %f\n",(c2Ptr[0]));
#endif
   
    ///////////////////////////////////////////////////////////////////////
    // Process number of cores
    ///////////////////////////////////////////////////////////////////////
    mNumCores = mxGetM(NUMCORES);
    nNumCores = mxGetN(NUMCORES);
    numCoresPtr = (DATATYPE*)mxGetPr(NUMCORES);

#ifdef FRDEBUG
    printf("numCores = %u\n",&(numCoresPtr[0]));
    printf("numCores = %f\n",(numCoresPtr[0]));
#endif
   
    ///////////////////////////////////////////////////////////////////////
    // Create a matrix for the return argument
    // The output matrix should have one item for
    // each spatial vector.
    ///////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////
    // Handle regular fractional integration case
    ///////////////////////////////////////////////////////////////////////
#ifndef INTERPOLATEONLY
#ifdef SINGLEPRECISION   
    // The last argument equal 0 indicates it is non-complex.
    //cout<<"Regular integration mode, single precision"<<endl;
    OUTPUT = mxCreateNumericMatrix(1, n, mxSINGLE_CLASS, 0);
#endif

#ifdef DOUBLEPRECISION
   // cout<<"Regular integration mode, double precision"<<endl;
    OUTPUT = mxCreateDoubleMatrix(1, n, mxREAL);
#endif
#endif

    ///////////////////////////////////////////////////////////////////////
    // Handle interpolation only case
    ///////////////////////////////////////////////////////////////////////
#ifdef INTERPOLATEONLY
#ifdef SINGLEPRECISION   
    // The last argument equal 0 indicates it is non-complex.
    // cout<<"Interpolate only, single precision"<<endl;
    OUTPUT = mxCreateNumericMatrix(2*m-1, n, mxSINGLE_CLASS, 0);
#endif

#ifdef DOUBLEPRECISION
    // cout<<"Interpolate only,double precision"<<endl;
    OUTPUT = mxCreateDoubleMatrix(2*m-1, n, mxREAL);
#endif
#endif
    ///////////////////////////////////////////////////////////////////////
    // Assign pointers to the various parameters
    ///////////////////////////////////////////////////////////////////////
    outputPtr = (DATATYPE*)mxGetPr(OUTPUT);

#ifdef FRDEBUG
    cout<<"address of outputPtr="<<&(outputPtr[0])<<endl;
    cout<<"outputPtr = "<<outputPtr[0]<<endl;
#endif
   
 
///////////////////////////////////////////////////////////////////////
// Call the HERMITE Custom Fully Proprietary code
// Works for both single and double precision.          
///////////////////////////////////////////////////////////////////////
int rc = 0; // keep track of the return code

#ifdef HERMITE
rc = hermite(dataMatrixPtr,m,n,deltatPtr,alphaPtr, c1Ptr, c2Ptr,numCoresPtr,outputPtr);
#endif

#ifdef INTERPOLATEONLY
#ifdef HERMITE
rc = hermiteinterpolateonly(dataMatrixPtr,m,n,deltatPtr,alphaPtr, c1Ptr, c2Ptr,numCoresPtr,outputPtr);
#endif
#endif

 
///////////////////////////////////////////////////////////////////////
// Call the SPLINE Custom Fully Proprietary code ( Cubic Splines )
// Works for both single and double precision.          
///////////////////////////////////////////////////////////////////////
#ifdef CUBIC
rc = cubic(dataMatrixPtr,m,n,deltatPtr,alphaPtr, c1Ptr, c2Ptr,numCoresPtr,outputPtr);
///////////////////////////////////////////////////////////////////////////
//  In certain situations, specially for the cubic spline case,
//  we may be unable to integrate a signal.
//  Add a way to check whether the call succeeded or failed.
///////////////////////////////////////////////////////////////////////////
if ( rc < 0 ) {
    cout<<"Cubic Splines results might be highly inaccurate due to coefficient matrix that is not diagonally dominant"<<endl;

}
#endif

#ifdef INTERPOLATEONLY
#ifdef CUBIC
rc = cubicinterpolateonly(dataMatrixPtr,m,n,deltatPtr,alphaPtr, c1Ptr, c2Ptr,numCoresPtr,outputPtr);
///////////////////////////////////////////////////////////////////////////
//  In certain situations, specially for the cubic spline case,
//  we may be unable to integrate a signal.
//  Add a way to check whether the call succeeded or failed.
///////////////////////////////////////////////////////////////////////////
if ( rc < 0 ) {
    cout<<"Cubic Splines results might be highly inaccurate due to coefficient matrix that is not diagonally dominant"<<endl;

}
#endif
#endif



return;
}
#endif

#ifdef MAININTERFACE
int main( int argc, char** argv) {
    ///////////////////////////////////////////////////////////////////////
    // This code is meant to be run within Matlab.
    // The main interface is primarily for testing purposes.
    //
    // Notes, on X86, single precision is emulated.
    // Using single precision will allow handling larger matrices
    // since only half the storage is used, however with
    // a slight performance degradation on architectures
    // that do not directly support single precision.
    //
    // todo items:
    // 2) handle signal lengths < 3 ( handled: switch case... )
    // 3) windows compiling.. ( tricky we can leave for a bit later )
    // 4) hermite  ( done )
    // 5) GSL integration ( later )
    // 6) GPU ( done in a preliminary manner. )
    // 7) MDCS ( just a script wrapping the C/C++, GPU code )
    //      - To prototype the fully distributed algorithm.
    //      - Specific to each application domain ( diffusion,... )
    // 8) Zero -order only cases. ( done)
    //    - Not worth using C version for zero-order cases.
    // 10) Add single precision for all cases
    // 11) Investigate Fully Distributed Algorithm
    //    - Part of MDCS ...[ major undertaking ]
    // 12) Investigate Real-Time Mode... according to meeting
    //     on 08/07/12. ( later )
    // 13) Implement custom code for all third-party library code.
    //     being used. ( e.g. Tri-diagonal solver, [do..] )
    // 14) Code cleanup, ensure unnecessary code is removed for
    //     enabling performance analysis.
    // 15) Memory scaling.  Using cells of signals. ( later )
    // 16) Log scale algorithm. ( later )
    // 17) Add second order corrections ( later )
    // 18) Add third order corrections ( later )
    // 19) Single precision stability issue (investigate further )
    //     - Perhaps consider whether to create the weights
    //       always in double precision?
    // 20) Windows compiling issue ( OpenMP support requires Visual
    //     Studio 2010 Professional, using trial version confirms this )
    //     (Compiles with OpenMP on Windows 7)
    ///////////////////////////////////////////////////////////////////////
    cout<<"FRINT TOOLBOX v1.0.0"<<endl;
    return 0;
}
#endif

