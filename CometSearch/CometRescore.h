/*
   Copyright 2012 University of Washington

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include "CometDataInternal.h"

struct RescoreThreadData
{
    int iQueryIndex;

    RescoreThreadData()
    {
       iQueryIndex = -1;
    }

    RescoreThreadData(int iQueryIndex_in)
    {
       iQueryIndex = iQueryIndex_in;
    }
};

class CometRescore
{
public:
    CometRescore();
    ~CometRescore();
    static bool Rescore(int minNumThreads, int maxNumThreads);
    static void RescoreThreadProc(RescoreThreadData *pThreadData);
    static void RescorePeptides(int iQueryIndex);
private:
    static vector<double> DeltamassCompare(double dDeltaXcorrMass = 0.0);
    
    //static double        _pdAAforwardClosest[MAX_PEPTIDE_LEN];
    //static double        _pdAAreverseClosest[MAX_PEPTIDE_LEN];
    //static unsigned int  _uiBinnedIonMassesClosest[MAX_FRAGMENT_CHARGE+1][9][MAX_PEPTIDE_LEN]; // Comet-PTM
    
    //struct iQuery
    //{
       //int        iQueryIndex;
    //};
};

//unsigned int  CometRescore::_uiBinnedIonMassesClosest[MAX_FRAGMENT_CHARGE+1][9][MAX_PEPTIDE_LEN]; // Comet-PTM
