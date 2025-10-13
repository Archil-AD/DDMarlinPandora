/**
 *  @file   DDMarlinPandora/src/DDCaloHitCreator.cc
 * 
 *  @brief  Implementation of the calo hit creator class.
 * 
 *  $Log: $
 */

#include "DDCaloHitCreatorALLEGRO.h"

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "UTIL/CellIDDecoder.h"


#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/DetType.h>
#include <DD4hep/DetectorSelector.h>
#include <DD4hep/Detector.h>

#include <algorithm>
#include <cmath>
#include <limits>

//forward declarations. See in DDPandoraPFANewProcessor.cc

dd4hep::rec::LayeredCalorimeterData * getExtension(unsigned int includeFlag, unsigned int excludeFlag=0);


DDCaloHitCreatorALLEGRO::DDCaloHitCreatorALLEGRO(const Settings &settings, const pandora::Pandora *const pPandora) :
    DDCaloHitCreator(settings, pPandora),
    m_hcalBarrelSegmentation(nullptr),
    m_hcalEndcapSegmentation(nullptr)
{
  dd4hep::Detector & theDetector = dd4hep::Detector::getInstance();

  for (unsigned int iSys = 0; iSys < m_settings.m_readoutNames.size(); iSys++) {
    // Check if readout exists
    streamlog_out(DEBUG) << "Readout: " << m_settings.m_readoutNames[iSys] << std::endl;
    if (theDetector.readouts().find(m_settings.m_readoutNames[iSys]) ==
        theDetector.readouts().end()) {
      streamlog_out(ERROR) << "Readout <<" << m_settings.m_readoutNames[iSys] << ">> does not exist." << std::endl;
      throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    // get segmentation
    dd4hep::DDSegmentation::Segmentation* aSegmentation =
        theDetector.readout(m_settings.m_readoutNames[iSys]).segmentation().segmentation();
    if (aSegmentation == nullptr) {
      streamlog_out(ERROR) << "Segmentation does not exist." << std::endl;
      throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    std::string segmentationType = aSegmentation->type();
    streamlog_out(DEBUG) << "Segmentation type : " << segmentationType << std::endl;

    if(m_settings.m_hcalBarrelSystemId == m_settings.m_readoutSystemId[iSys])
      m_hcalBarrelSegmentation = std::shared_ptr<dd4hep::DDSegmentation::Segmentation>(aSegmentation,[](dd4hep::DDSegmentation::Segmentation*) {
          // no-op deleter: don't delete the pointer
          });
    else
      m_hcalEndcapSegmentation = std::shared_ptr<dd4hep::DDSegmentation::Segmentation>(aSegmentation,[](dd4hep::DDSegmentation::Segmentation*) {
          // no-op deleter: don't delete the pointer
          });
  }

  // make sure that we got readout segmentation objects
  if(!m_hcalBarrelSegmentation || !m_hcalEndcapSegmentation)
  {
    streamlog_out(ERROR) << "Unable to get segmentation objects." << std::endl;
    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDCaloHitCreatorALLEGRO::~DDCaloHitCreatorALLEGRO()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDCaloHitCreatorALLEGRO::CreateCaloHits(const EVENT::LCEvent *const pLCEvent)
{
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateECalCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateMuonCaloHits(pLCEvent));

    return pandora::STATUS_CODE_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------------------------------------

void DDCaloHitCreatorALLEGRO::GetCommonCaloHitProperties(const EVENT::CalorimeterHit *const pCaloHit, PandoraApi::CaloHit::Parameters &caloHitParameters) const
{
    const float *pCaloHitPosition(pCaloHit->getPosition());
    const pandora::CartesianVector positionVector(pCaloHitPosition[0], pCaloHitPosition[1], pCaloHitPosition[2]);

    // AD: HCAL cells are rectangular. If the hit is in the HCAL get the dimensions from the segmentation class
    if(caloHitParameters.m_hitType.Get() == pandora::HCAL)
      caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
    else caloHitParameters.m_cellGeometry = pandora::POINTING_THETAPHI;

    caloHitParameters.m_positionVector = positionVector;
    caloHitParameters.m_expectedDirection = positionVector.GetUnitVector();
    caloHitParameters.m_pParentAddress = (void*)pCaloHit;
    caloHitParameters.m_inputEnergy = pCaloHit->getEnergy();
    caloHitParameters.m_time = pCaloHit->getTime();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDCaloHitCreatorALLEGRO::GetEndCapCaloHitProperties(const EVENT::CalorimeterHit *const pCaloHit, const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer> &layers,
    PandoraApi::CaloHit::Parameters &caloHitParameters, float &absorberCorrection) const
{
    caloHitParameters.m_hitRegion = pandora::ENDCAP;

    //FIXME! WHAT DO WE DO HERE?
    const int physicalLayer(std::min(static_cast<int>(caloHitParameters.m_layer.Get()), static_cast<int>(layers.size()-1)));

    // AD: HCAL cells are rectangular. if the hit is in the HCAL get the dimensions from the segmentation class
    if(caloHitParameters.m_hitType.Get() == pandora::HCAL)
    {
      std::vector<double>  cellSize = m_hcalEndcapSegmentation->cellDimensions(pCaloHit->getCellID0());
      caloHitParameters.m_cellSize0 = cellSize[0]/dd4hep::mm;
      caloHitParameters.m_cellSize1 = cellSize[1]/dd4hep::mm;
    }
    else
    {
      caloHitParameters.m_cellSize0 = layers[physicalLayer].cellSize0;
      caloHitParameters.m_cellSize1 = layers[physicalLayer].cellSize1;
    }
    double thickness = (layers[physicalLayer].inner_thickness+layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;
    double nRadLengths = layers[physicalLayer].inner_nRadiationLengths;
    double nIntLengths = layers[physicalLayer].inner_nInteractionLengths;
    double layerAbsorberThickness = (layers[physicalLayer].inner_thickness-layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;

    if(physicalLayer>0){
        thickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;
        nRadLengths += layers[physicalLayer-1].outer_nRadiationLengths;
        nIntLengths += layers[physicalLayer-1].outer_nInteractionLengths;
        layerAbsorberThickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;

    }

    caloHitParameters.m_cellThickness = thickness;
    caloHitParameters.m_nCellRadiationLengths = nRadLengths;
    caloHitParameters.m_nCellInteractionLengths = nIntLengths;

    if (caloHitParameters.m_nCellRadiationLengths.Get() < std::numeric_limits<float>::epsilon() || caloHitParameters.m_nCellInteractionLengths.Get() < std::numeric_limits<float>::epsilon())
    {
        streamlog_out(WARNING) << "CaloHitCreator::GetEndCapCaloHitProperties Calo hit has 0 radiation length or interaction length: \
            not creating a Pandora calo hit." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }


    /*
    //FIXME! do we need this?
    absorberCorrection = 1.;
    for (unsigned int i = 0, iMax = layers.size(); i < iMax; ++i)
    {
        float absorberThickness((layers[i].inner_thickness - layers[i].sensitive_thickness/2.0 )/dd4hep::mm);

        if (i>0)
            absorberThickness += (layers[i-1].outer_thickness - layers[i-1].sensitive_thickness/2.0)/dd4hep::mm;

        if (absorberThickness < std::numeric_limits<float>::epsilon())
            continue;

        if (layerAbsorberThickness > std::numeric_limits<float>::epsilon())
            absorberCorrection = absorberThickness / layerAbsorberThickness;

        break;
    }
    */

    caloHitParameters.m_cellNormalVector = (pCaloHit->getPosition()[2] > 0) ? pandora::CartesianVector(0, 0, 1) :
        pandora::CartesianVector(0, 0, -1);

//     streamlog_out(DEBUG) <<" GetEndCapCaloHitProperties: physLayer: "<<physicalLayer <<" layer: "<<caloHitParameters.m_layer.Get()<<" nX0: "<<    caloHitParameters.m_nCellRadiationLengths.Get() <<" nLambdaI: "<<    caloHitParameters.m_nCellInteractionLengths.Get()<<" thickness: "<<caloHitParameters.m_cellThickness.Get()<<std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDCaloHitCreatorALLEGRO::GetBarrelCaloHitProperties( const EVENT::CalorimeterHit *const pCaloHit,
                                                   const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer> &layers,
                                                   unsigned int barrelSymmetryOrder,
                                                   PandoraApi::CaloHit::Parameters &caloHitParameters,
                                                   FloatVector const& normalVector,
                                                   float &absorberCorrection ) const
{
    caloHitParameters.m_hitRegion = pandora::BARREL;

    // GM for debug
    // streamlog_out(DEBUG) << "In GetBarrelCaloHitProperties" << std::endl;
    // streamlog_out(DEBUG) << "Nlayers: " << layers.size() << std::endl;
    // streamlog_out(DEBUG) << "Hit layer: " << caloHitParameters.m_layer.Get() << std::endl;
    // streamlog_out(DEBUG) << "Hit type: " << caloHitParameters.m_hitType.Get() << std::endl;

    // GM check if we need to handle the case of projective cells rather than fixed size..
    const int physicalLayer(std::min(static_cast<int>(caloHitParameters.m_layer.Get()), static_cast<int>(layers.size()-1)));

    // AD: HCAL cells are rectangular. If the hit is in the HCAL get the dimensions from the segmentation class
    if(caloHitParameters.m_hitType.Get() == pandora::HCAL)
    {
      std::vector<double>  cellSize = m_hcalBarrelSegmentation->cellDimensions(pCaloHit->getCellID0());
      caloHitParameters.m_cellSize0 = cellSize[0]/dd4hep::mm;
      caloHitParameters.m_cellSize1 = cellSize[1]/dd4hep::mm;
    }
    else
    {
      caloHitParameters.m_cellSize0 = layers[physicalLayer].cellSize0;
      caloHitParameters.m_cellSize1 = layers[physicalLayer].cellSize1;
    }

    /*
    streamlog_out(MESSAGE) << "In GetBarrelCaloHitProperties" << std::endl;
    streamlog_out(MESSAGE) << "Hit type: " << caloHitParameters.m_hitType.Get() << std::endl;
    streamlog_out(MESSAGE) << "Hit cellGeometry: " << caloHitParameters.m_cellGeometry.Get() << std::endl;
    streamlog_out(MESSAGE) << "Hit cellSize0: " << caloHitParameters.m_cellSize0.Get() << std::endl;
    streamlog_out(MESSAGE) << "Hit cellSize1: " << caloHitParameters.m_cellSize1.Get() << std::endl;
    */

    double thickness = (layers[physicalLayer].inner_thickness+layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;
    double nRadLengths = layers[physicalLayer].inner_nRadiationLengths;
    double nIntLengths = layers[physicalLayer].inner_nInteractionLengths;

    double layerAbsorberThickness = (layers[physicalLayer].inner_thickness-layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;
    if(physicalLayer>0){
        thickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;
        nRadLengths += layers[physicalLayer-1].outer_nRadiationLengths;
        nIntLengths += layers[physicalLayer-1].outer_nInteractionLengths;
        layerAbsorberThickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;
    }

    caloHitParameters.m_cellThickness = thickness;
    caloHitParameters.m_nCellRadiationLengths = nRadLengths;
    caloHitParameters.m_nCellInteractionLengths = nIntLengths;

    if (caloHitParameters.m_nCellRadiationLengths.Get() < std::numeric_limits<float>::epsilon() || caloHitParameters.m_nCellInteractionLengths.Get() < std::numeric_limits<float>::epsilon())
    {
        streamlog_out(WARNING) << "CaloHitCreator::GetBarrelCaloHitProperties Calo hit has 0 radiation length or interaction length: \
            not creating a Pandora calo hit." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    //FIXME! do we need this?
    absorberCorrection = 1.;
    for (unsigned int i = 0, iMax = layers.size(); i < iMax; ++i)
    {
        float absorberThickness((layers[i].inner_thickness - layers[i].sensitive_thickness/2.0 )/dd4hep::mm);

        if (i>0)
            absorberThickness += (layers[i-1].outer_thickness - layers[i-1].sensitive_thickness/2.0)/dd4hep::mm;

        if (absorberThickness < std::numeric_limits<float>::epsilon())
            continue;

        if (layerAbsorberThickness > std::numeric_limits<float>::epsilon())
            absorberCorrection = absorberThickness / layerAbsorberThickness;

        break;
    }

    if (barrelSymmetryOrder > 2)
    {

      if ( pCaloHit->getCellID0() != 0 ) {

        auto staveDetElement = m_volumeManager.lookupDetElement( pCaloHit->getCellID0() );
        dd4hep::Position local1(0.0, 0.0, 0.0);
        dd4hep::Position local2(normalVector[0],normalVector[1],normalVector[2]);
        dd4hep::Position global1(0.0, 0.0, 0.0);
        dd4hep::Position global2(0.0, 0.0, 0.0);
        staveDetElement.nominal().localToWorld( local1, global1 );
        staveDetElement.nominal().localToWorld( local2, global2 );
        dd4hep::Position normal( global2-global1 );

        streamlog_out(DEBUG6) << "   detelement: " << staveDetElement.name()
                              << "   parent: " << staveDetElement.parent().name()
                              << "   grandparent: " << staveDetElement.parent().parent().name()
                              << "   cellID: " << pCaloHit->getCellID0()
                              << "   PhiLoc:"  << atan2( global1.y(), global1.x() )*180/M_PI
                              << "   PhiNor:"  << atan2( normal.y(), normal.x() )*180/M_PI
                              << " normal vector "
                              << std::setw(15) << normal.x()
                              << std::setw(15) << normal.y()
                              << std::setw(15) << normal.z()
                              << std::endl;

        caloHitParameters.m_cellNormalVector = pandora::CartesianVector( normal.x(), normal.y(), normal.z() );
      } else {
        const double phi = atan2( pCaloHit->getPosition()[1], pCaloHit->getPosition()[0] );

        streamlog_out( WARNING )  << "This hit does not have any cellIDs set, will use phi-direction for normal vector "
                                  << " phi:" << std::setw(15) << phi*180/M_PI
                                  << std::endl;

        caloHitParameters.m_cellNormalVector = pandora::CartesianVector( std::cos(phi), std::sin(phi) , 0.0 );
      }



    }
    else
    {
        const float *pCaloHitPosition( pCaloHit->getPosition() );
        const float phi = std::atan2( pCaloHitPosition[1], pCaloHitPosition[0] );
        caloHitParameters.m_cellNormalVector = pandora::CartesianVector(std::cos(phi), std::sin(phi), 0);
    }

//     streamlog_out(DEBUG)<<" GetBarrelCaloHitProperties: physLayer: "<<physicalLayer <<" layer: "<<caloHitParameters.m_layer.Get()<<" nX0: "<<    caloHitParameters.m_nCellRadiationLengths.Get() <<" nLambdaI: "<<    caloHitParameters.m_nCellInteractionLengths.Get()<<" thickness: "<<caloHitParameters.m_cellThickness.Get()<<std::endl;
}

