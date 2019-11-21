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

#include "Common.h"
#include "CometDataInternal.h"
#include "ThreadPool.h"
#include "CometRescore.h"
#include "CometMassSpecUtils.h"
#include "CometStatus.h"

CometRescore::CometRescore()
{
}

CometRescore::~CometRescore()
{
}

bool CometRescore::Rescore(int minNumThreads,
                           int maxNumThreads)
{
   bool bSucceeded = true;

   // Create the thread pool containing g_staticParams.options.iNumThreads,
   // each hanging around and sleeping until asked to do a rescore.
   ThreadPool<RescoreThreadData *> *pRescoreThreadPool = new ThreadPool<RescoreThreadData*>(RescoreThreadProc,
         minNumThreads, maxNumThreads);

   for (int i=0; i<(int)g_pvQuery.size(); i++)
   {
      RescoreThreadData *pThreadData = new RescoreThreadData(i);
      pRescoreThreadPool->Launch(pThreadData);

      bSucceeded = !g_cometStatus.IsError() && !g_cometStatus.IsCancel();
      if (!bSucceeded)
      {
         break;
      }
   }

   // Wait for active post analysis threads to complete processing.
   //pRescoreThreadPool->WaitForThreads();
   
   delete pRescoreThreadPool;
   pRescoreThreadPool = NULL;

   // Check for errors one more time since there might have been an error
   // while we were waiting for the threads.
   if (bSucceeded)
   {
      bSucceeded = !g_cometStatus.IsError() && !g_cometStatus.IsCancel();
   }

   return bSucceeded;
}

void CometRescore::RescoreThreadProc(RescoreThreadData *pThreadData)
{
   int iQueryIndex = pThreadData->iQueryIndex;

   //iQuery iRescore;
   //iRescore.iQueryIndex = iQueryIndex;
   RescorePeptides(iQueryIndex);
   //iRescore.RescorePeptides(iRescore);

   delete pThreadData;
   pThreadData = NULL;
}

vector<double> CometRescore::DeltamassCompare(double dDeltaXcorrMass)
{
	// Gets all theoretical deltamasses within range

	std::vector<double> vMasses = g_staticParams.vectorDeltaPeaks;
	std::vector<double> vClosest;

	for (int k = 0; k < vMasses.size(); k++)
	{
		double dMass = vMasses[k];
		double dDistance = abs(dMass - dDeltaXcorrMass);
		if (dDistance <= g_staticParams.tolerances.dDeltaTolerance) // within tolerance
			vClosest.push_back(dMass);
	}

	if (vClosest.empty())
		vClosest.push_back(-9999);

	return vClosest;
}

void CometRescore::RescorePeptides(int iQueryIndex)
{
    // REPEAT FOR DECOYS! (DECOYS HAS A DIFFERENT SIZE!)

    Query* pQuery = g_pvQuery.at(iQueryIndex);
    int iSize = pQuery->iMatchPeptideCount;
    if (iSize > g_staticParams.options.iNumStored)
      iSize = g_staticParams.options.iNumStored;

    for (int i=0; i<iSize; i++)
    {
        Threading::LockMutex(pQuery->accessMutex); // Do we need a mutex here?
        
       if (pQuery->_pResults[i].bDoRecom == 1) // Rescore only marked candidates
       {
            int bestPosClosest = -1;
            double dXcorrClosest;
            double dBestXcorr;
            
            // Vectors to store the results for each deltamass
            std::vector<double> vdXcorrClosest;
            std::vector<double> vdBestXcorr;
            //std::vector<std::string> vszXcorrProfileClosest;
            //std::vector<std::string> vszXcorrProfileBClosest;
            //std::vector<std::string> vszXcorrProfileYClosest;
            std::vector<std::string> vszXcorrType;
            std::vector<int> viBestPosClosest;

            //std::string szXcorrProfileClosest = "";
            //std::string szXcorrProfileBClosest = "";
            //std::string szXcorrProfileYClosest = "";
            std::string szXcorrType = "";

            // Get information from this Result
            int iStartPos = pQuery->_pResults[i].iStartPos;
            int iEndPos = pQuery->_pResults[i].iEndPos;
            int iLenMinus1 = iEndPos - iStartPos; // Equals iLenPeptide minus 1.
            int iProteinSeqLengthMinus1 = pQuery->_pResults[i].vcProteinSeq.size() - 1;
            int iModPos = pQuery->_pResults[i].iModPos;
            float fXcorr = pQuery->_pResults[i].fXcorr;
            double dDeltaXcorrMass = pQuery->_pResults[i].dDeltaXcorrMass;
            std::vector<char> vzProteinSeq = pQuery->_pResults[i].vcProteinSeq;
            char* szProteinSeq = reinterpret_cast<char*>(vzProteinSeq.data());
            
            //std::vector<std::vector<std::vector<std::vector<double> > > > vctXcorrs = pQuery->_pResults[i].vctXcorrs;

            // Get deltamasses within range
            std::vector<double> vClosest;
            vClosest = DeltamassCompare(dDeltaXcorrMass);
            
           if (vClosest[0] > -9999) // No deltamasses found (TODO: do this in a more elegant way)
           {
            // Calculate new score for each deltamass
            for (int s = 0; s < vClosest.size(); s++) // For each theoretical deltamass within range
            {
                std::vector<std::vector<std::vector<std::vector<double> > > > vctXcorrs = pQuery->_pResults[i].vctXcorrs;
                // Reset variables
                double dDeltaXcorrMassClosest = vClosest[s];
                unsigned int  _uiTempBinnedIonMassesClosest[MAX_FRAGMENT_CHARGE+1][9][MAX_PEPTIDE_LEN];
                unsigned int (*p_uiTempBinnedIonMassesClosest)[MAX_FRAGMENT_CHARGE+1][9][MAX_PEPTIDE_LEN];
                double _pdAAforwardClosest[MAX_PEPTIDE_LEN];
                double _pdAAreverseClosest[MAX_PEPTIDE_LEN];
                int i;
                double dBion = g_staticParams.precalcMasses.dNtermProton;
                double dYion = g_staticParams.precalcMasses.dCtermOH2Proton;

                bool* pbTempDuplFragmentClosest = pQuery->pbDuplFragment;

                if (iStartPos == 0)
                    dBion += g_staticParams.staticModifications.dAddNterminusProtein;
                if (iEndPos == iProteinSeqLengthMinus1)
                    dYion += g_staticParams.staticModifications.dAddCterminusProtein;

                int iPos;
                for (i = iStartPos; i < iEndPos; i++)
                {
                    iPos = i - iStartPos;

                    dBion += g_staticParams.massUtility.pdAAMassFragment[(int)szProteinSeq[i]];
                    double dBionDeltaClosest = dBion + dDeltaXcorrMassClosest;
                    if (dBionDeltaClosest >= 0.0)
                        _pdAAforwardClosest[iPos] = dBionDeltaClosest;
                    else
                        _pdAAforwardClosest[iPos] = 0.0;

                    dYion += g_staticParams.massUtility.pdAAMassFragment[(int)szProteinSeq[iEndPos - iPos]];
                    double dYionDeltaClosest = dYion + dDeltaXcorrMassClosest;
                    if (dYionDeltaClosest >= 0.0)
                        _pdAAreverseClosest[iPos] = dYionDeltaClosest;
                    else
                        _pdAAreverseClosest[iPos] = 0.0;
                }

                // Now get the set of binned fragment ions once to compare this peptide against all matching spectra.
                for (int ctCharge = 1; ctCharge <= g_massRange.iMaxFragmentCharge; ctCharge++)
                {
                   for (int ctIonSeries = 0; ctIonSeries < g_staticParams.ionInformation.iNumIonSeriesUsed; ctIonSeries++)
                   {
                      int iWhichIonSeries = g_staticParams.ionInformation.piSelectedIonSeries[ctIonSeries];
                      for (int ctLen = 0; ctLen < iLenMinus1; ctLen++)
                      {
                          pbTempDuplFragmentClosest[BIN(CometMassSpecUtils::GetFragmentIonMass(iWhichIonSeries, ctLen, ctCharge, _pdAAforwardClosest, _pdAAreverseClosest))] = false; // TODO this crashes sometimes ?
                      }
                   }
                }
                
                for (int ctCharge = 1; ctCharge <= g_massRange.iMaxFragmentCharge; ctCharge++)
                {
                   for (int ctIonSeries = 0; ctIonSeries < g_staticParams.ionInformation.iNumIonSeriesUsed; ctIonSeries++)
                   {
                      int iWhichIonSeries = g_staticParams.ionInformation.piSelectedIonSeries[ctIonSeries];

                      // As both _pdAAforwardClosest and _pdAAreverseClosest are increasing, loop through
                      // iLenPeptide-1 to complete set of internal fragment ions.
                      for (int ctLen = 0; ctLen < iLenMinus1; ctLen++)
                      {
                         int iValDeltaClosest = BIN(CometMassSpecUtils::GetFragmentIonMass(iWhichIonSeries, ctLen, ctCharge, _pdAAforwardClosest, _pdAAreverseClosest));
                         unsigned int  _uiBinnedIonMassesClosest[MAX_FRAGMENT_CHARGE+1][9][MAX_PEPTIDE_LEN];
                         _uiTempBinnedIonMassesClosest[ctCharge][ctIonSeries][ctLen] = _uiBinnedIonMassesClosest[ctCharge][ctIonSeries][ctLen]; // reset

                         if (pbTempDuplFragmentClosest[iValDeltaClosest] == false)
                         {
                            _uiTempBinnedIonMassesClosest[ctCharge][ctIonSeries][ctLen] = iValDeltaClosest;
                            pbTempDuplFragmentClosest[iValDeltaClosest] = true;
                         }
                         else
                            _uiTempBinnedIonMassesClosest[ctCharge][ctIonSeries][ctLen] = 0;
                      }
                   }
                }

                if (!g_staticParams.variableModParameters.bRequireVarMod)
                {
                   // spiros: XcorrScore now takes into account deltaXcorr

                   int ctLen,
                       ctIonSeries,
                       ctCharge;
                   double dXcorrScale = 0.005; // Scale intensities to 50 and divide score by 1E4.
                   int iLenPeptideMinus1 = iEndPos - iStartPos;
                
                   // Pointer to either regular or decoy uiBinnedIonMasses[][][].
                   p_uiTempBinnedIonMassesClosest = &_uiTempBinnedIonMassesClosest; 
                   
                   int iWhichIonSeries;
                   bool bUseNLPeaks = false;

                   float **ppSparseFastXcorrData;              // use this if bSparseMatrix
                   float *pFastXcorrData;                      // use this if not using SparseMatrix

                   dXcorrClosest = 0.0;

                   //int iMax = pQuery->_spectrumInfoInternal.iArraySize / SPARSE_MATRIX_SIZE + 1;

                   // spiros: store all the fragment Xcorrs and DeltaXcorrs in this vector to calculate where the modification most likely is
//                    std::vector<std::vector<std::vector<std::vector<double> > > > vctXcorrs(pQuery->_spectrumInfoInternal.iMaxFragCharge,   // per charge
//                       std::vector<std::vector<std::vector<double> > >(g_staticParams.ionInformation.iNumIonSeriesUsed,                 // per series
//                          std::vector<std::vector<double> >(3, std::vector<double>(iLenPeptideMinus1, 0))));   // one for Xcorr, one for ClosestDeltaXcorr
                   //std::vector<std::vector<std::vector<std::vector<double> > > > vctXcorrs = pQuery->_pResults[i].vctXcorrs;

                   const short XCORR_ID = 0;
                   const short DELTA_XCORR_ID = 1;
                   const short DELTA_XCORR_CLOSEST_ID = 2;

                   for (ctCharge = 1; ctCharge <= pQuery->_spectrumInfoInternal.iMaxFragCharge; ctCharge++)
                   {
                      for (ctIonSeries = 0; ctIonSeries<g_staticParams.ionInformation.iNumIonSeriesUsed; ctIonSeries++)
                      {
                         iWhichIonSeries = g_staticParams.ionInformation.piSelectedIonSeries[ctIonSeries];

                         if (g_staticParams.ionInformation.bUseNeutralLoss && (iWhichIonSeries == ION_SERIES_A || iWhichIonSeries == ION_SERIES_B || iWhichIonSeries == ION_SERIES_Y))
                            bUseNLPeaks = true;
                         else
                            bUseNLPeaks = false;

                         if (ctCharge == 1 && bUseNLPeaks)
                         {
                            ppSparseFastXcorrData = pQuery->ppfSparseFastXcorrDataNL;
                            pFastXcorrData = pQuery->pfFastXcorrDataNL;
                         }
                         else
                         {
                            ppSparseFastXcorrData = pQuery->ppfSparseFastXcorrData;
                            pFastXcorrData = pQuery->pfFastXcorrData;
                         }
                         int binDeltaClosest, xDeltaClosest, yDeltaClosest;
                         //bool bDeltaValid = false;
                         //bool bDeltaClosestValid = false;

                         for (ctLen = 0; ctLen<iLenPeptideMinus1; ctLen++)
                         {
                            //MH: newer sparse matrix converts bin to sparse matrix bin
                            //bDeltaClosestValid = false;
                            if (g_staticParams.options.bUseDeltaXcorr)
                            {
                               binDeltaClosest = *(*(*(*p_uiTempBinnedIonMassesClosest + ctCharge) + ctIonSeries) + ctLen);
                               xDeltaClosest = binDeltaClosest / SPARSE_MATRIX_SIZE;

                               if (xDeltaClosest < pQuery->iFastXcorrData && ppSparseFastXcorrData[xDeltaClosest] != NULL)
                               {
                                  yDeltaClosest = binDeltaClosest - (xDeltaClosest*SPARSE_MATRIX_SIZE);
                                  //bDeltaClosestValid = true;
                                  dXcorrClosest += ppSparseFastXcorrData[xDeltaClosest][yDeltaClosest]; //andrea deltamass
                                  //vctXcorrs[ctCharge - 1][ctIonSeries][DELTA_XCORR_ID][ctLen] = ppSparseFastXcorrData[xDeltaClosest][yDeltaClosest];
                                  vctXcorrs[ctCharge - 1][ctIonSeries][DELTA_XCORR_CLOSEST_ID][ctLen] = ppSparseFastXcorrData[xDeltaClosest][yDeltaClosest];
                               }
                            }
                         }
                      }
                   }

                   // comet-ptm spiros start: make this a separate function
                   std::vector<double> vctXcorrAtPosClosest(iLenPeptideMinus1 + 1, 0);
                   std::vector<double> vctXcorrAtPosBClosest(iLenPeptideMinus1 + 1, 0);
                   std::vector<double> vctXcorrAtPosYClosest(iLenPeptideMinus1 + 1, 0);
                   std::vector<double>::iterator bestClosest;

                   // spiros
                   // encodes b and y series modifications, b series at [0], y series at [1]
                   // std::vector<std::bitset<MAX_PEPTIDE_LEN> > vctBitMods(2, std::bitset<64>());

                   if (g_staticParams.options.bUseDeltaXcorr && iLenPeptideMinus1 > 0)
                   {
                      // calculate the best Xcorr: for each residue, calculate the Xcorr with the modification at that residue
                      // marco fixed bug (previously it was just "<", so the last amino acid was never checked
                      for (int j = 0; j <= iLenPeptideMinus1; ++j)
                      {
                         for (int i = 0; i < iLenPeptideMinus1; ++i)
                         {
                            for (ctIonSeries = 0; ctIonSeries < g_staticParams.ionInformation.iNumIonSeriesUsed; ctIonSeries++)
                            {
                               iWhichIonSeries = g_staticParams.ionInformation.piSelectedIonSeries[ctIonSeries];

                               if (iWhichIonSeries == ION_SERIES_B)   // b series
                               {
                                  if (j - i > 0)
                                  {
                                     for (ctCharge = 1; ctCharge <= pQuery->_spectrumInfoInternal.iMaxFragCharge; ctCharge++)
                                     {
                                        vctXcorrAtPosClosest[j] += vctXcorrs[ctCharge - 1][ctIonSeries][XCORR_ID][i];    // deltaXcorrClosest
                                        vctXcorrAtPosBClosest[j] += vctXcorrs[ctCharge - 1][ctIonSeries][XCORR_ID][i];    // deltaXcorrClosest only for B series
                                     }
                                  }
                                  else
                                  {
                                     for (ctCharge = 1; ctCharge <= pQuery->_spectrumInfoInternal.iMaxFragCharge; ctCharge++)
                                     {
                                        vctXcorrAtPosClosest[j] += vctXcorrs[ctCharge - 1][ctIonSeries][DELTA_XCORR_CLOSEST_ID][i];    // deltaXcorrClosest
                                        vctXcorrAtPosBClosest[j] += vctXcorrs[ctCharge - 1][ctIonSeries][DELTA_XCORR_CLOSEST_ID][i];    // deltaXcorrClosest only for B series
                                     }
                                  }
                               }
                               else if (iWhichIonSeries == ION_SERIES_Y) // y series
                               {
                                  if (iLenPeptideMinus1 - j - i > 0)
                                  {
                                     for (ctCharge = 1; ctCharge <= pQuery->_spectrumInfoInternal.iMaxFragCharge; ctCharge++)
                                     {
                                        vctXcorrAtPosClosest[j] += vctXcorrs[ctCharge - 1][ctIonSeries][XCORR_ID][i];    // deltaXcorrClosest
                                        vctXcorrAtPosYClosest[j] += vctXcorrs[ctCharge - 1][ctIonSeries][XCORR_ID][i];    // deltaXcorrClosest only for Y series
                                     }
                                  }
                                  else
                                  {
                                     for (ctCharge = 1; ctCharge <= pQuery->_spectrumInfoInternal.iMaxFragCharge; ctCharge++)
                                     {
                                        vctXcorrAtPosClosest[j] += vctXcorrs[ctCharge - 1][ctIonSeries][DELTA_XCORR_CLOSEST_ID][i];    // deltaXcorrClosest
                                        vctXcorrAtPosYClosest[j] += vctXcorrs[ctCharge - 1][ctIonSeries][DELTA_XCORR_CLOSEST_ID][i];    // deltaXcorrClosest only for Y series
                                     }
                                  }
                               }
                            }
                         }
                      }
                      
                      if (g_staticParams.options.bPrintXcorrProfiles)
                      {
                         // filling the strings with the Xcorr profile
                         for (int j = 0; j <= iLenPeptideMinus1; ++j)
                         {
                            // andrea deltamass start
                            //szXcorrProfileClosest += std::to_string(vctXcorrAtPosClosest[j] * dXcorrScale);
                            //szXcorrProfileBClosest += std::to_string(vctXcorrAtPosBClosest[j] * dXcorrScale);
                            //szXcorrProfileYClosest += std::to_string(vctXcorrAtPosYClosest[j] * dXcorrScale);
                            if (j < iLenPeptideMinus1)
                            {
                               //szXcorrProfileClosest += ", ";
                               //szXcorrProfileBClosest += ", ";
                               //szXcorrProfileYClosest += ", ";
                            }
                         }
                      }

                      bestClosest = std::max_element(vctXcorrAtPosClosest.begin(), vctXcorrAtPosClosest.end());
                      bestPosClosest = std::distance(vctXcorrAtPosClosest.begin(), bestClosest);
                   }

                   if (bestPosClosest != -1)
                   dXcorrClosest = *bestClosest;

                   if (dXcorrClosest < XCORR_CUTOFF)
                      dXcorrClosest = XCORR_CUTOFF;
                   else
                      dXcorrClosest *= dXcorrScale;

                   if (!g_staticParams.options.bUseDeltaClosest || dDeltaXcorrMassClosest == -9999)
                   {
                      dXcorrClosest = 0;
                      //szXcorrProfileClosest = "0";
                      //szXcorrProfileBClosest = "0";
                      //szXcorrProfileYClosest = "0";
                   }
                   
                   if (fXcorr > dXcorrClosest)
                   {
                      dBestXcorr = fXcorr;
                      szXcorrType = "experimental";
                      viBestPosClosest.push_back(iModPos); // Keep previous position
                   }
                   else if (fXcorr == dXcorrClosest)
                   {
                      dBestXcorr = fXcorr;
                      szXcorrType = "same";
                      viBestPosClosest.push_back(iModPos); // Keep previous position
                   }
                   else if (fXcorr < dXcorrClosest)
                   {
                      dBestXcorr = dXcorrClosest;
                      szXcorrType = "theoretical";
                      viBestPosClosest.push_back(bestPosClosest); // Update position
                   }
                }

                // Add results from XcorrScore() to the output vectors
                vdXcorrClosest.push_back(dXcorrClosest);
                vdBestXcorr.push_back(dBestXcorr);
                //vszXcorrProfileClosest.push_back(szXcorrProfileClosest);
                //vszXcorrProfileBClosest.push_back(szXcorrProfileBClosest);
                //vszXcorrProfileYClosest.push_back(szXcorrProfileYClosest);
                vszXcorrType.push_back(szXcorrType);
            }

            if (!vdXcorrClosest.empty())
            {
               // Get index of highest scoring theoretical deltamass
               std::vector<double>::iterator best_score = max_element(vdXcorrClosest.begin(), vdXcorrClosest.end());
               int best_index = std::distance(vdXcorrClosest.begin(), best_score);

                // Get all results for the candidate with that deltamass
               double dDeltaXcorrMassClosest = vClosest[best_index];
               double dXcorrClosest = vdXcorrClosest[best_index];
               double dBestXcorr = vdBestXcorr[best_index];
               //std::string szXcorrProfileClosest = vszXcorrProfileClosest[best_index];
               //std::string szXcorrProfileBClosest = vszXcorrProfileBClosest[best_index];
               //std::string szXcorrProfileYClosest = vszXcorrProfileYClosest[best_index];
               std::string szXcorrType = vszXcorrType[best_index];
               int iBestPos = viBestPosClosest[best_index];
               
               // Calculate corrected Xcorr
               double dXcorrCorrClosest;
               if (g_staticParams.options.bUseXcorrCorr)
               {
                  if (pQuery->_spectrumInfoInternal.iChargeState < 3) // R = 1
                  {
                     dXcorrCorrClosest = log10(fXcorr) / log10(2*(iEndPos - iStartPos + 1)/110);
                  }
                  else // R = 1.22
                  {
                     float fXcorr_R = fXcorr/1.22;
                     dXcorrCorrClosest = log10(fXcorr_R) / log10(2*(iEndPos - iStartPos + 1)/110);
                  }
               }
               
               if (g_staticParams.options.iDecoySearch != 2)
               {
                  // Store results
                  pQuery->_pResults[i].dDeltaXcorrMassClosest = dDeltaXcorrMassClosest;
                  pQuery->_pResults[i].fXcorrClosest = (float)dXcorrClosest;
                  if (g_staticParams.options.bUseXcorrCorr)
                     pQuery->_pResults[i].fXcorrCorrClosest = (float)dXcorrCorrClosest;
                  pQuery->_pResults[i].fBestXcorr = (float)dBestXcorr;
                  pQuery->_pResults[i].iModPos = iBestPos;
                  //pQuery->_pResults[i].szXcorrProfileClosest = szXcorrProfileClosest;
                  //pQuery->_pResults[i].szXcorrProfileBClosest = szXcorrProfileBClosest;
                  //pQuery->_pResults[i].szXcorrProfileYClosest = szXcorrProfileYClosest;
                  strcpy(pQuery->_pResults[i].szXcorrType, const_cast<char*>(szXcorrType.c_str()));
                }
            }
            
           }
           else
           {
              pQuery->_pResults[i].dDeltaXcorrMassClosest = -9999;
              pQuery->_pResults[i].fXcorrClosest = 0;
              pQuery->_pResults[i].fBestXcorr = pQuery->_pResults[i].fXcorr;
              strcpy(pQuery->_pResults[i].szXcorrType, "experimental");
           }
           
           if (pQuery->_pResults[i].fNonModXcorr > pQuery->_pResults[i].fBestXcorr)
           {
              pQuery->_pResults[i].fBestXcorr = pQuery->_pResults[i].fNonModXcorr;
              strcpy(pQuery->_pResults[i].szXcorrType, "nonmod");   
           }
       }
       else
       {
          pQuery->_pResults[i].dDeltaXcorrMassClosest = -9999;
          pQuery->_pResults[i].fXcorrClosest = 0;
          pQuery->_pResults[i].fBestXcorr = pQuery->_pResults[i].fXcorr;
          strcpy(pQuery->_pResults[i].szXcorrType, "not rescored");
       }
        
       Threading::UnlockMutex(pQuery->accessMutex);
    }
}
