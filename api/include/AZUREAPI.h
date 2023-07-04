#ifndef AZUREAPI_H
#define AZUREAPI_H

#include "AZUREMain.h"

#include "Constants.h"
#include <vector>

class Config;
class EData;
class CNuc;

///A function class to perform the calculation of the chi-squared value

/*!
 * The AZUREAPI function class calculates the cross section based on a 
 * parameter set for all available data, and returns a chi-squared value.
 * This function class is what Minuit calls repeatedly during the fitting
 * process to perform the minimization.
 */

class AZUREAPI {
 public:
  /*!
   * The AZUREAPI object is created with reference to an EData and CNuc object.
   *. The runtime configurations are also passed through a Config structure.
   */
  AZUREAPI(Config& configure) : configure_(configure) { };
  
  ~AZUREAPI() {};

  bool Initialize( );

  // Update data objects
  bool UpdateData( );
  // Calculate the external capture for data
  bool CalculateExternalCapture( );
  // Reads the parameters values
  bool UpdateParameters( );
  // Calculated segments values
  bool UpdateSegments(vector_r& p);
  
  /*!
   * Returns a reference to the Config structure.
   */
  Config &configure() const {return configure_;};
  /*!
   * Returns a pointer to the EData object.
   */
  EData *data() const {return data_;};
  /*!
   * Returns a pointer to the CNuc object.
   */
  CNuc *compound() const {return compound_;};
  /*!
   * Returns a pointer to the parameter values object.
   */
  vector_r params_values() const {return transform_;};
  /*!
   * Returns a pointer to the parameter names object.
   */
  std::vector<std::string> params_names() const {return names_;};
  /*!
   * Returns a pointer to the parameter names object.
   */
  std::vector<bool> params_fixed() const {return fixed_;};
  /*!
   * Returns a pointer to the calculated segments object.
   */
  vector_r data_energies() const {return dataEnergies_;};
  /*!
   * Returns a pointer to the calculated segments object.
   */
  vector_r data_segments() const {return dataSegments_;};
  /*!
   * Returns a pointer to the calculated segments object.
   */
  vector_r data_segments_errors() const {return dataSegmentsErrors_;};
  /*!
   * Returns a pointer to the calculated segments object.
   */
  vector_r calculated_segments() const {return calculatedSegments_;};
 
 private:

  // Configuration
  Config &configure_;
  AZUREMain* azureMain_;
  EData *data_;
  CNuc *compound_;

  // Parameters
  std::vector<std::string> names_;
  std::vector<bool> fixed_;
  vector_r values_, transform_;

  // Data
  vector_r dataEnergies_;
  vector_r dataSegments_;
  vector_r dataSegmentsErrors_;
  vector_r calculatedSegments_;

};

#endif