%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FRACTIONAL INTEGRATION TOOLBOX
% RELEASE 1.0
% JUNE 2013
% 
% A toolbox from the Santamaria Laboratory at the 
% University of Texas at San Antonio. 
% ( http://utsa.edu/Santamarialab/index.htm )
% In collaboration
% with the Computational Systems Biology Core/Computational Biology
% Initiative( www.cbi.utsa.edu )
% 
% Grant Acknowledgments:
% "This work received computational support from 
% Computational System Biology Core, funded by the National Institute 
% on Minority Health and Health Disparities (G12MD007591) from the 
% National Institutes of Health."
%
% This work was supported by grants EF-1137897 and HDR-0932339 from the
% National Science Foundation
%
% Software Developer(s): 
% Toma.Marinov@utsa.edu
% Nelson.Ramirez@utsa.edu
%
% Project PI:
% Fidel.Santamaria@utsa.edu
%
%  Fractional Integration Toolbox 1.0
%  UTSA RESEARCH LICENSE (SOURCE CODE)
%  The University of Texas at San Antonio has developed certain software and
%  documentation that it desires to make available without charge to anyone for
%  academic, research, experimental or personal use. This license is designed
%  to guarantee freedom to use the software for these purposes. If you wish to
%  distribute or make other use of the software, you may purchase a license to
%  do so from the University of Texas.
%  The accompanying source code is made available to you under the terms of
%  this UT Research License (this "UTRL"). By clicking the "ACCEPT" button,
%  or by installing or using the code, you are consenting to be bound by
%  this UTRL. If you do not agree to the terms and conditions of this license,
%  do not click the "ACCEPT" button, and do not install or use any part
%  of the code.
%  The terms and conditions in this UTRL not only apply to the source code
%  made available by UT, but also to any improvements to, or derivative works
%  of, that source code made by you and to any object code compiled from such
%  source code, improvements or derivative works.
%  1. DEFINITIONS.
%  1.1 "Commercial Use" shall mean use of Software or Documentation by Licensee
%  for direct or indirect financial, commercial or strategic gain or advantage,
%  including without limitation: (a) bundling or integrating the Software with
%  any hardware product or another software product for transfer, sale or license
%  to a third party (even if distributing the Software on separate media and not
%  charging for the Software); (b) providing customers with a link to the
%  Software or a copy of the Software for use with hardware or another software
%  product purchased by that customer; or (c) use in connection with the
%  performance of services for which Licensee is compensated.
%  1.2 "Derivative Products" means any improvements to, or other derivative
%  works of, the Software made by Licensee.
%  1.3 "Documentation" shall mean all manuals, user documentation, and other
%  related materials pertaining to the Software that are made available to
%  Licensee in connection with the Software.
%  1.4 "Licensor" shall mean The University of Texas.
%  1.5 "Licensee" shall mean the person or entity that has agreed to the
%  terms hereof and is exercising rights granted hereunder.
%  1.6 "Software" shall mean the computer program(s) referred to as
%  "Fractional Integration Toolbox 1.0" made available under this UTRL in
%  source code form, including any error corrections, bug fixes, patches,
%  updates or other modifications that Licensor may in its sole discretion
%  make available to Licensee from time to time, and any object code compiled
%  from such source code.
%  2. GRANT OF RIGHTS.
%  Subject to the terms and conditions hereunder, Licensor hereby grants to
%  Licensee a worldwide, non- transferable, non-exclusive license to (a) install,
%  use and reproduce the Software for academic, research,
%  experimental and personal use (but specifically excluding Commercial Use);
%  (b) use and modify the Software to create Derivative Products, subject
%  to Section 3.2; and (c) use the Documentation, if any, solely in connection
%  with Licensee's authorized use of the Software.
%  3. RESTRICTIONS; COVENANTS.
%  3.1 Licensee may not: (a) distribute, sub-license or otherwise transfer copies
%  or rights to the Software (or any portion thereof) or the Documentation;
%  (b) use the Software (or any portion thereof) or Documentation for
%  Commercial Use, or for any other use except as described in Section 2;
%  (c) copy the Software or Documentation other than for archival and
%  backup purposes; or (d) remove any product identification, copyright,
%  proprietary notices or labels from the Software and Documentation.
%  This UTRL confers no rights upon Licensee except those expressly granted herein.
%  3.2 Licensee hereby agrees that it will provide a copy of all Derivative
%  Products to Licensor and that its use of the Derivative Products will be
%  subject to all of the same terms, conditions, restrictions and limitations
%  on use imposed on the Software under this UTRL. Licensee hereby grants
%  Licensor a worldwide, non- exclusive, royalty-free license to reproduce,
%  prepare derivative works of, publicly display, publicly
%  perform, sublicense and distribute Derivative Products. Licensee also hereby
%  grants Licensor a worldwide, non-exclusive, royalty-free patent license to
%  make, have made, use, offer to sell, sell, import and otherwise transfer the
%  Derivative Products under those patent claims licensable by Licensee that
%  are necessarily infringed by the Derivative Products.
%  4. PROTECTION OF SOFTWARE.
%  4.1 Confidentiality. The Software and Documentation are the confidential and
%  proprietary information of Licensor. Licensee agrees to take adequate steps to
%  protect the Software and Documentation from unauthorized disclosure or use.
%  Licensee agrees that it will not disclose the Software or Documentation to
%  any third party.
%  4.2 Proprietary Notices. Licensee shall maintain and place on any copy of
%  Software or Documentation that it reproduces for internal use all notices
%  as are authorized and/or required hereunder. Licensee shall include a copy
%  of this UTRL and the following notice, on each copy of the Software and
%  Documentation. Such license and notice shall be embedded in each copy of
%  the Software, in the video screen display, on the physical medium embodying
%  the Software copy and on any Documentation:
%  Copyright © 2013, The University of Texas at San Antonio. All rights reserved.
%  UNIVERSITY EXPRESSLY DISCLAIMS ANY AND ALL WARRANTIES CONCERNING THIS SOFTWARE
%  AND DOCUMENTATION, INCLUDING ANY WARRANTIES OF MERCHANTABILITY, FITNESS FOR
%  ANY PARTICULAR PURPOSE, NON- INFRINGEMENT AND WARRANTIES OF PERFORMANCE, AND
%  ANY WARRANTY THAT MIGHT OTHERWISE ARISE FROM COURSE OF DEALING OR USAGE OF
%  TRADE. NO WARRANTY IS EITHER EXPRESS OR IMPLIED WITH RESPECT TO THE USE OF
%  THE SOFTWARE OR DOCUMENTATION. Under no circumstances shall University be
%  liable for incidental, special, indirect, direct or consequential damages
%  or loss of profits, interruption of business, or related expenses which
%  may arise from use of Software or Documentation, including but not limited
%  to those resulting from defects in Software and/or Documentation, or loss
%  or inaccuracy of data of any kind.
%  5. WARRANTIES.
%  5.1 Disclaimer of Warranties. TO THE EXTENT PERMITTED BY APPLICABLE LAW,
%  THE SOFTWARE AND DOCUMENTATION ARE BEING PROVIDED ON AN "AS IS" BASIS WITHOUT
%  ANY WARRANTIES OF ANY KIND RESPECTING THE SOFTWARE OR DOCUMENTATION, EITHER
%  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO ANY WARRANTY OF DESIGN,
%  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR NON-INFRINGEMENT.
%  5.2 Limitation of Liability. UNDER NO CIRCUMSTANCES UNLESS REQUIRED BY
%  APPLICABLE LAW SHALL LICENSOR BE LIABLE FOR INCIDENTAL, SPECIAL, INDIRECT,
%  DIRECT OR CONSEQUENTIAL DAMAGES OR LOSS OF PROFITS, INTERRUPTION OF BUSINESS,
%  OR RELATED EXPENSES WHICH MAY ARISE AS A RESULT OF THIS LICENSE OR OUT OF
%  THE USE OR ATTEMPT OF USE OF SOFTWARE OR DOCUMENTATION INCLUDING BUT NOT
%  LIMITED TO THOSE RESULTING FROM DEFECTS IN SOFTWARE AND/OR DOCUMENTATION,
%  OR LOSS OR INACCURACY OF DATA OF ANY KIND. THE FOREGOING EXCLUSIONS AND
%  LIMITATIONS WILL APPLY TO ALL CLAIMS AND ACTIONS OF ANY KIND, WHETHER
%  BASED ON CONTRACT, TORT (INCLUDING, WITHOUT LIMITATION, NEGLIGENCE), OR
%  ANY OTHER GROUNDS.
%  6. INDEMNIFICATION.
%  Licensee shall indemnify, defend and hold harmless Licensor, the
%  University of Texas System, their Regents, and their officers,
%  agents and employees from and against any claims, demands, or
%  causes of action whatsoever caused by, or arising out of, or resulting
%  from, the exercise or practice of the license granted hereunder by
%  Licensee, its officers, employees, agents or representatives.
%  7. TERMINATION.
%  If Licensee breaches this UTRL, Licensee's right to use the Software
%  and Documentation will terminate immediately without notice, but all
%  provisions of this UTRL except Section 2 will survive termination and
%  continue in effect. Upon termination, Licensee must destroy all copies
%  of the Software and Documentation.
%  8. GOVERNING LAW; JURISDICTION AND VENUE.
%  The validity, interpretation, construction and performance of this UTRL
%  shall be governed by the laws of the State of Texas. The Texas state
%  courts of Travis County, Texas (or, if there is exclusive federal
%  jurisdiction,the United States District Court for the Central
%  District of Texas) shall have exclusive jurisdiction and venue
%  over any dispute arising out of this UTRL, and Licensee consents
%  to the jurisdiction of such courts. Application of the
%  United Nations Convention on Contracts for the International
%  Sale of Goods is expressly excluded.
%  9. EXPORT CONTROLS.
%  This license is subject to all applicable export restrictions. Licensee
%  must comply with all export and import laws and restrictions and
%  regulations of any United States or foreign agency or authority relating
%  to the Software and its use.
%  10. U.S. GOVERNMENT END-USERS.
%  The Software is a "commercial item," as that term is defined in
%  48 C.F.R. 2.101, consisting of "commercial computer software" and
%  "commercial computer software documentation," as such terms are
%  used in 48 C.F.R. 12.212 (Sept. 1995) and 48 C.F.R. 227.7202 (June 1995).
%  Consistent with 48 C.F.R. 12.212, 48 C.F.R. 27.405(b)(2) (June 1998)
%  and 48 C.F.R. 227.7202, all U.S. Government End Users acquire the
%  Software with only those rights as set forth herein.
%  11. MISCELLANEOUS
%  If any provision hereof shall be held illegal, invalid or unenforceable,
%  in whole or in part, such provision shall be modified to the minimum extent
%  necessary to make it legal, valid and enforceable, and the legality, validity
%  and enforceability of all other provisions of this UTRL shall not be affected
%  thereby. Licensee may not assign this UTRL in whole or in part, without
%  Licensor's prior written consent. Any attempt to assign this UTRL without
%  such consent will be null and void. This UTRL is the complete and exclusive
%  statement between Licensee and Licensor relating to the subject matter
%  hereof and supersedes all prior oral and written and all contemporaneous
%  oral negotiations, commitments and understandings of the parties, if any.
%  Any waiver by either party of any default or breach hereunder shall
%  not constitute a waiver of any provision of this UTRL or of any
%  subsequent default or breach of the same or a different kind.
%  END OF LICENSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MinGW64 Compiling Process on Windows 64 bit systems( e.g. Windows 7):
% Compiling Setup Notes:
% Step 1: Install minGW64
% (http://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/)
% Step 2: Set the Windows Path environment variable to the
%  minGW64 bin directory: C:\MinGW64\bin
% Step 3: Make sure to restart Matlab.
% Step 4: edit(fullfile(prefdir,'mexopts.bat'))
% Refer to:  http://stackoverflow.com/questions/8552580/using-gccmingw-as-matlabs-mex-compiler   
% A sample mexopts.bat file:
% @echo off
% set MINGWPATH="C:\MinGW64"
% set PATH=%MINGWPATH%\bin;%PATH%
% set COMPILER=x86_64-w64-mingw32-g++
% set COMPFLAGS=-c -m64 -I"%MATLAB%\extern\include" -DMATLAB_MEX_FILE -Wall -fopenmp 
% set OPTIMFLAGS=-O3 -DNDEBUG
% set DEBUGFLAGS=-g
% set NAME_OBJECT=-o
% set LINKER=x86_64-w64-mingw32-g++
% set LINKFLAGS=-shared -L"%MATLAB%\bin\win64" -L"%MATLAB%\extern\lib\win64" -lmex -lmx -leng -lmat -lmwlapack -lmwblas -fopenmp
% set NAME_OUTPUT=-o "%OUTDIR%%MEX_NAME%%MEX_EXT%"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HERMITE ( With First Correction( First Order ) ): Double Precision
mex -v frint.cpp -output frint_hermite_double COMPILER="x86_64-w64-mingw32-g++" COMPFLAGS="$COMPFLAGS -O3 -fopenmp -DMATLABINTERFACE -DHERMITE -DOPENMP -DDOUBLEPRECISION -DFIRSTORDER -DNDEBUG" LINKER="x86_64-w64-mingw32-g++" LINKFLAGS="$LINKFLAGS -fopenmp" DEBUGFLAGS="" OPTIMFLAGS="-O3 -DNDEBUG"

% HERMITE Interpolation Only: Double Precision
mex -v frint.cpp -output interpolate_hermite_double COMPILER="x86_64-w64-mingw32-g++" COMPFLAGS="$COMPFLAGS -O3 -fopenmp -DMATLABINTERFACE -DINTERPOLATEONLY -DHERMITE -DOPENMP -DDOUBLEPRECISION -DNDEBUG" LINKER="x86_64-w64-mingw32-g++" LINKFLAGS="$LINKFLAGS -fopenmp" DEBUGFLAGS="" OPTIMFLAGS="-O3 -DNDEBUG"

% CUBIC ( With First Correction( First Order ) ): Double Precision
mex -v frint.cpp -output frint_cubic_double COMPILER="x86_64-w64-mingw32-g++" COMPFLAGS="$COMPFLAGS -O3 -fopenmp -DMATLABINTERFACE -DCUBIC -DOPENMP -DDOUBLEPRECISION -DFIRSTORDER -DNDEBUG" LINKER="x86_64-w64-mingw32-g++" LINKFLAGS="$LINKFLAGS -fopenmp" DEBUGFLAGS="" OPTIMFLAGS="-O3 -DNDEBUG"

% CUBIC Interpolation Only: Double Precision
mex -v frint.cpp -output interpolate_cubic_double COMPILER="x86_64-w64-mingw32-g++" COMPFLAGS="$COMPFLAGS -O3 -fopenmp -DMATLABINTERFACE -DINTERPOLATEONLY -DCUBIC -DOPENMP -DDOUBLEPRECISION -DNDEBUG" LINKER="x86_64-w64-mingw32-g++" LINKFLAGS="$LINKFLAGS -fopenmp" DEBUGFLAGS="" OPTIMFLAGS="-O3 -DNDEBUG"
