#ifndef GENMATRIXFUNC_H
#define GENMATRIXFUNC_H

#include "Constants.h"
#include "Config.h"

class EPoint;
class CNuc;

/// A temporaray T-Matrix structure

/*!
 * The TempTMatrix structure is used to coherently add T-matrix elements from 
 * pathways  with like \f$ J,l,l' \f$ values for the calculation of angle integrated cross section.
 * This is primarly used to facilitate the interference between internal and external pathways.
 */

struct TempTMatrix {
  ///Total spin value of temporary matrix element
  double jValue;
  /// Entrance orbital angular momentum for temporary matrix element
  int lValue;
  /// Exit orbital angular momentum for temporary matrix element
  int lpValue;
  ///Value of temporary matrix element
  complex TMatrix;
};


///A generalized function class to calculate cross sections

/*!
 * The GenMatrixFunc function class is the general form of the function used 
 * to calculate cross section from R-Matrix parameters.  It is the parent class of
 * AMatrixFunc and RMatrixFunc.
 */

class GenMatrixFunc {
 public:
  GenMatrixFunc() {};
  virtual ~GenMatrixFunc(){};
  /*!
   *This virtual function in implemented in the child class.
   */
  virtual void ClearMatrices()=0;
  /*!
   *This virtual function in implemented in the child class.
   */
  virtual void FillMatrices(EPoint*)=0;
  /*!
   *This virtual function in implemented in the child class.
   */
  virtual void InvertMatrices()=0;
  /*!
   *This virtual function in implemented in the child class.
   */
  virtual void CalculateTMatrix(EPoint*)=0;
  void CalculateCrossSection(EPoint*);
  void NewTempTMatrix(TempTMatrix);
  void NewTempTMatrixE1(TempTMatrix);
  void NewTempTMatrixE2(TempTMatrix);
  void AddToTempTMatrix(int,complex);
  void AddToTempTMatrixE1(int,complex);
  void AddToTempTMatrixE2(int,complex);
  void ClearTempTMatrices();
  void ClearTempTMatricesE1();
  void ClearTempTMatricesE2();
  void AddTMatrixElement(int,int,complex,int decayNum=1);
  void AddECTMatrixElement(int,int,complex);
  int IsTempTMatrix(double,int,int);
  int IsTempTMatrixE1(double,int,int);
  int IsTempTMatrixE2(double,int,int);
  int NumTempTMatrices() const;
  int NumTempTMatricesE1() const;
  int NumTempTMatricesE2() const;
  TempTMatrix *GetTempTMatrix(int);
  TempTMatrix *GetTempTMatrixE1(int);
  TempTMatrix *GetTempTMatrixE2(int);
  complex GetTMatrixElement(int,int,int decayNum=1) const;
  complex GetECTMatrixElement(int,int) const;

  /*!
   *This virtual function in implemented in the child class.
   */
  virtual CNuc *compound() const = 0;
  /*!
   *This virtual function in implemented in the child class.
   */
  virtual const Config& configure() const = 0;
 protected:
  ///Vector of internal T-matrix elements accessable to child class.
  std::vector<matrix_c> tmatrix_;
  ///Vector of external T-matrix elements accessable to child class.
  matrix_c ec_tmatrix_;
 private:
  std::vector<TempTMatrix> temp_t_matrices_;
  std::vector<TempTMatrix> temp_t_matrices_E1_;
  std::vector<TempTMatrix> temp_t_matrices_E2_;
};

#endif
