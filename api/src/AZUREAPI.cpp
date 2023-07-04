#include "AZUREAPI.h"
#include "AZUREParams.h"

#include "GSLException.h"

#include "Config.h"
#include "CNuc.h"
#include "EData.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>

bool AZUREAPI::Initialize( ){

  configure().paramMask |= Config::USE_EXTERNAL_CAPTURE;

  std::string file;
  if( configure().paramMask & Config::CALCULATE_WITH_DATA ) file = configure().outputdir + "intEC.dat";
  else file = configure().outputdir + "intEC.extrap";

  std::ifstream in(file.c_str());
  if( !in ) configure().paramMask &= ~Config::USE_PREVIOUS_INTEGRALS;
  else configure().paramMask |= Config::USE_PREVIOUS_INTEGRALS;

  configure().integralsfile=file;

  azureMain_ = new AZUREMain( configure( ) );

  data_ = azureMain_->data();
  compound_ = azureMain_->compound( );

  //configure().outStream << "Filling Compound Nucleus..." << std::endl;
  if(compound()->Fill(configure())==-1) {
    //configure().outStream << "Could not fill compound nucleus from file." << std::endl;
    return -1;
  } else if(compound()->NumPairs()==0 || compound()->NumJGroups()==0) {
    //configure().outStream << "No nuclear data exists. Calculation not possible." << std::endl; 
    return -1;
  } 
  if((configure().screenCheckMask|configure().fileCheckMask) & 
     Config::CHECK_COMPOUND_NUCLEUS) compound()->PrintNuc(configure());

  if(!(configure().paramMask & Config::CALCULATE_REACTION_RATE)) {
    //Fill the data object from the segments and data file
    //  Compound object is passed to the function for pair key verification and
    //  center of mass conversions, s-factor conversions, etc.
    //configure().outStream << "Filling Data Structures..." << std::endl;
    if(configure().paramMask & Config::CALCULATE_WITH_DATA) {
      if(data()->Fill(configure(),compound())==-1) {
	//configure().outStream << "Could not fill data object from file." << std::endl;
	return -1;
      } else if(data()->NumSegments()==0) {
	//configure().outStream << "There is no data provided." << std::endl;
	return -1;
      }
    } else {
      if(data()->MakePoints(configure(),compound())==-1) {
	//configure().outStream << "Could not fill data object from file." << std::endl;
	return -1;
      } else if(data()->NumSegments()==0) {
	//configure().outStream << "Extrapolation segments produce no data." << std::endl;
	return -1;
      }
    } 
    if((configure().fileCheckMask|configure().screenCheckMask) & Config::CHECK_DATA)
      data()->PrintData(configure());
  } else {
    if(!compound()->IsPairKey(configure().rateParams.entrancePair)||!compound()->IsPairKey(configure().rateParams.exitPair)) {
      //configure().outStream << "Reaction rate pairs do not exist in compound nucleus." << std::endl;
      return -1;
    } else {
      compound()->GetPair(compound()->GetPairNumFromKey(configure().rateParams.entrancePair))->SetEntrance();
    }
  }

  //Initialize compound nucleus object
  try {
    compound()->Initialize(configure());
  } catch (GSLException e) {
    configure().outStream << e.what() << std::endl;
    configure().outStream << std::endl
			  << "Calculation was aborted." << std::endl;
    return -1;
  }

  if(data()->Initialize(compound(),configure())==-1) return -1;

  return 0;
  
}

bool AZUREAPI::UpdateParameters( ) {

  names_.clear( );
  fixed_.clear( );
  values_.clear( );
  transform_.clear( );

  AZUREParams params;
  compound()->FillMnParams(params.GetMinuitParams());
  data()->FillMnParams(params.GetMinuitParams());

  compound()->FillCompoundFromParams(params.GetMinuitParams( ).Params( ));

  compound()->CalcShiftFunctions( configure() );
  compound()->TransformOut( configure() );

  for(int i = 0; i < params.GetMinuitParams().Params().size(); i++){
    names_.push_back( params.GetMinuitParams().Parameter(i).GetName() );
    fixed_.push_back( params.GetMinuitParams().Parameter(i).IsFixed() );
    values_.push_back( params.GetMinuitParams().Parameter(i).Value() );
  }

  transform_ = compound()->GetTransformParams( configure() );

  return true;

}

bool AZUREAPI::UpdateSegments(vector_r& p) {

  calculatedSegments_.clear( );

  CNuc* localCompound = NULL;
  EData* localData = NULL;
  localCompound = compound()->Clone();
  localData = data()->Clone();

  AZUREParams params;
  localCompound->FillCompoundFromParamsPhysical(p);
  bool isValid = localCompound->TransformIn( configure( ) );

  if( !isValid ){ 
    calculatedSegments_.push_back( std::numeric_limits<double>::infinity() );
    return false;
  }

  localCompound->FillMnParams(params.GetMinuitParams());
  localData->FillMnParams(params.GetMinuitParams());

  localCompound->FillCompoundFromParams(params.GetMinuitParams( ).Params( ));

  //Fill Compound Nucleus From Minuit Parameters
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());

  //loop over segments and points
  ESegmentIterator firstSumIterator = localData->GetSegments().end();
  ESegmentIterator lastSumIterator = localData->GetSegments().end();
  for(EDataIterator data=localData->begin();data!=localData->end();data++) {
    if(data.segment()->GetPoints().begin()==data.point()) {

      if(data.segment()->IsTotalCapture()) {
	      firstSumIterator=data.segment();
	      lastSumIterator=data.segment()+data.segment()->IsTotalCapture()-1;
      }

    }
    
    if(!data.point()->IsMapped()) data.point()->Calculate(localCompound,configure());
    

    if(firstSumIterator!=localData->GetSegments().end()&& data.segment()!=lastSumIterator) continue;
    
    double fitCrossSection=data.point()->GetFitCrossSection();
    
    ESegmentIterator thisSegment = data.segment();
    if(data.segment()==lastSumIterator) {
      int pointIndex=data.point()-data.segment()->GetPoints().begin()+1;
      for(ESegmentIterator it=firstSumIterator;it<data.segment();it++) 
	      fitCrossSection+=it->GetPoint(pointIndex)->GetFitCrossSection();
        thisSegment = firstSumIterator;
    }

    calculatedSegments_.push_back( fitCrossSection );

  }

  return true;

}

bool AZUREAPI::CalculateExternalCapture( ){

  configure().paramMask |= Config::USE_EXTERNAL_CAPTURE;

  configure().paramMask &= ~Config::USE_PREVIOUS_INTEGRALS;
  data()->CalculateECAmplitudes( compound( ), configure( ) );
  configure().paramMask |= Config::USE_PREVIOUS_INTEGRALS;

  return true;

}

bool AZUREAPI::UpdateData( ) {

  dataEnergies_.clear( );
  dataSegments_.clear( );
  dataSegmentsErrors_.clear( );

  CNuc* localCompound = NULL;
  EData* localData = NULL;
  localCompound = compound()->Clone();
  localData = data()->Clone();

  //loop over segments and points
  ESegmentIterator firstSumIterator = localData->GetSegments().end();
  ESegmentIterator lastSumIterator = localData->GetSegments().end();
  for(EDataIterator data=localData->begin();data!=localData->end();data++) {
    if(data.segment()->GetPoints().begin()==data.point()) {

      if(data.segment()->IsTotalCapture()) {
	      firstSumIterator=data.segment();
	      lastSumIterator=data.segment()+data.segment()->IsTotalCapture()-1;
      }

    }    

    if(firstSumIterator!=localData->GetSegments().end()&& data.segment()!=lastSumIterator) continue;
    
    double energy=data.point()->GetCMEnergy( );
    double crossSection=data.point()->GetCMCrossSection();
    double crossSectionError=data.point()->GetCMCrossSectionError();

    dataEnergies_.push_back( energy );
    dataSegments_.push_back( crossSection );
    dataSegmentsErrors_.push_back( crossSectionError );

  }

  return true;

}