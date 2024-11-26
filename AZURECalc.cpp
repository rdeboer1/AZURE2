#include "AZURECalc.h"
#include "Config.h"
#include "CNuc.h"
#include "EData.h"
#include <iostream>
#include <iomanip>

double AZURECalc::operator()(const vector_r& p) const {

  int thisIteration=data()->Iterations();
  data()->Iterate();
  bool isFit=data()->IsFit();

  CNuc * localCompound = NULL;
  EData *localData = NULL;
  if(isFit) {
    localCompound = compound()->Clone();
    localData = data()->Clone();
  } else {
    localCompound = compound();
    localData = data();
  }

  //Fill Compound Nucleus From Minuit Parameters
  localCompound->FillCompoundFromParams(p);
  localData->FillNormsFromParams(p);
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());
  
  //loop over segments and points
  double chiSquared=0.0;
  double segmentChiSquared=0.0;
  ESegmentIterator firstSumIterator = localData->GetSegments().end();
  ESegmentIterator lastSumIterator = localData->GetSegments().end();
  for(EDataIterator data=localData->begin();data!=localData->end();data++) {
    if(data.segment()->GetPoints().begin()==data.point()) {
      segmentChiSquared=0.0;
      if(data.segment()->IsTotalCapture()) {
	firstSumIterator=data.segment();
	lastSumIterator=data.segment()+data.segment()->IsTotalCapture()-1;
      } 
    }
    if(!data.point()->IsMapped()) data.point()->Calculate(localCompound,configure());
    if(firstSumIterator!=localData->GetSegments().end()&&
       data.segment()!=lastSumIterator) continue;
    double fitCrossSection=data.point()->GetFitCrossSection();
    ESegmentIterator thisSegment = data.segment();
    if(data.segment()==lastSumIterator) {
      int pointIndex=data.point()-data.segment()->GetPoints().begin()+1;
      for(ESegmentIterator it=firstSumIterator;it<data.segment();it++) 
	fitCrossSection+=it->GetPoint(pointIndex)->GetFitCrossSection();
      thisSegment = firstSumIterator;
    }
    double dataNorm=thisSegment->GetNorm();
    double CrossSection=data.point()->GetCMCrossSection()*dataNorm;
    double CrossSectionError=data.point()->GetCMCrossSectionError()*dataNorm;
    double chi=(fitCrossSection-CrossSection)/CrossSectionError;
    double pointChiSquared=pow(chi,2.0);
//  enable for alternate goodness of fit function
/*
    if(pointChiSquared > 1.0e-15){  
      pointChiSquared = -1.0*log((1.0-exp(-1.0*pointChiSquared/2.0))/pointChiSquared);
    }
*/
    segmentChiSquared+=pointChiSquared;
    if(data.segment()->GetPoints().end()-1==data.point()) {
      if(!isFit) thisSegment->SetSegmentChiSquared(segmentChiSquared);
      if(data.segment()==lastSumIterator) {
	firstSumIterator=localData->GetSegments().end();
	lastSumIterator=localData->GetSegments().end();
      }
      double dataNormNominal=thisSegment->GetNominalNorm();
      double dataNormError=dataNormNominal/100.*thisSegment->GetNormError();
      if(dataNormError!=0.) {
        double normChiSquared = pow((dataNorm-dataNormNominal)/dataNormError,2.0);
        if(!isFit) thisSegment->SetNormChiSquared(normChiSquared);
	segmentChiSquared += normChiSquared;
/*
        if(segmentChiSquared > 1.0e-15){
           double normChiSquared = pow((dataNorm-dataNormNominal)/dataNormError,2.0);
           segmentChiSquared -= normChiSquared;	
	   segmentChiSquared += -1.0*log((1.0-exp(-1.0*normChiSquared/2.0))/normChiSquared);
	}
*/
      }
      chiSquared+=segmentChiSquared;
    }
  }

/*  for(int j=1;j<=localCompound->NumJGroups();j++) {
    for(int la=1;la<=localCompound->GetJGroup(j)->NumLevels();la++) {
      for(int ch=1;ch<=localCompound->GetJGroup(j)->NumChannels();ch++) {
        if(gammaW != 0.0 && abs(p.at(i)) > gammaW){
          wignerChiSquared = pow((abs(p.at(i)) - gammaW)/p.at(i),2.0);

//      std::cout << "Param value: " << p.at(i) << ", gW: " << gammaW << ", c2: " << wignerChiSquared << std::endl;     
        }
      }
    }
  }
*/

  //use chiflag.dat to control chi2 options
  std::ifstream chioption;
  chioption.open("chiflag.dat");
  int chiflag=0; // 0 normal chi2, 1 Wigner limit
  if(chioption){
    chioption >> chiflag;
    chioption.close();
  }
//  std::cout<<chiflag<<std::endl;
  if(chiflag==1){
//  std::cout << "Param size: " << p.size() << ", Wigner size: " << localCompound->GetWignerLimit().size() << std::endl;

    for (size_t i=0;i<localCompound->GetWignerLimit().size();i++){
  
      double gammaW = localCompound->GetWignerLimit().at(i);
      double wignerChiSquared = 0.0;

      if(gammaW != 0.0 && abs(p.at(i)) > gammaW){
        wignerChiSquared = pow((abs(p.at(i)) - gammaW)/gammaW,2.0);
//      if(paramChiSquared > 1.0e-15){
//        paramChiSquared = -1.0*log((1.0-exp(-1.0*paramChiSquared/2.0))/paramChiSquared); 
//      std::cout << "Param value: " << p.at(i) << ", gW: " << gammaW << ", c2: " << wignerChiSquared << std::endl;     
        }
        chiSquared += wignerChiSquared;
    }  
  }

  if(!localData->IsErrorAnalysis()&&thisIteration!=0) {
    if(thisIteration%100==0) configure().outStream
			       << "\r\tIteration: " << std::setw(6) << thisIteration
			       << " Chi-Squared: " << chiSquared;  configure().outStream.flush();

    if(thisIteration%1000==0) {
      localData->WriteOutputFiles(configure(),isFit);
      localCompound->TransformOut(configure());
      localCompound->PrintTransformParams(configure());
    }
  }
  if(isFit) {
    delete localCompound;
    delete localData;
  }
  if(configure().stopFlag&&isFit) return 0.;
  else return chiSquared;
}
