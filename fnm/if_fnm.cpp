#include <sps/memory> // Figure out why it is needed to import this before <memory>
#include <if_fnm.h>
#include <fnm.hpp>

int ApertureCreate(Aperture** obj) {
  *obj = (Aperture*) new fnm::Aperture<float>();
  return 0;
}

int ApertureDestroy(Aperture* obj) {
  auto p = (fnm::Aperture<float>*) obj;
  delete p;
  return 0;
}

int ApertureLinearCreate(Aperture** obj, size_t nElements, float width, float kerf, float height) {
  *obj = (Aperture*) new fnm::Aperture<float>(nElements, width, kerf, height);
  return 0;
}

#ifdef FNM_CLOSURE_FUNCTIONS

int ApertureRwFloatParamSet(Aperture* obj, int fsel, const float* f, size_t nDim, ...) {
  auto p = (fnm::Aperture<float>*) obj;

  va_list args;
  va_start(args, nDim);
  int retval = p->RwFloatParamSet(fsel, f, nDim, args);
  va_end(args);
  return retval;
}

int ApertureRwFloatParamGet(Aperture* obj, int fsel, float** f, size_t nDim, ...) {
  auto p = (fnm::Aperture<float>*) obj;

  va_list args;
  va_start(args, nDim);
  int retval = p->RwFloatParamGet(fsel, f, nDim, args);
  va_end(args);
  return retval;
}

int FreeCArray(void* pData) {
  if (pData) {
    free(pData);
    return 0;
  }
  return -1;
}

#endif

int ApertureNDivWSet(Aperture* obj, size_t nDiv) {
  auto p = (fnm::Aperture<float>*) obj;
  p->NDivWSet(nDiv);
  return 0;
}

size_t ApertureNDivWGet(Aperture* obj) {
  auto p = (fnm::Aperture<float>*) obj;
  return p->NDivWGet();
}

int ApertureNDivHSet(Aperture* obj, size_t nDiv) {
  auto p = (fnm::Aperture<float>*) obj;
  p->NDivHSet(nDiv);
  return 0;
}

size_t ApertureNDivHGet(Aperture* obj) {
  auto p = (fnm::Aperture<float>*) obj;
  return p->NDivHGet();
}

int ApertureNThreadsSet(Aperture* obj, size_t nthreads) {
  auto p = (fnm::Aperture<float>*) obj;
  p->NThreadsSet(nthreads);
  return 0;
}

size_t ApertureNThreadsGet(Aperture* obj) {
  auto p = (fnm::Aperture<float>*) obj;
  return p->NThreadsGet();
}

int ApertureFocusingTypeSet(Aperture* obj, int ftype) {
  auto p = (fnm::Aperture<float>*) obj;
  p->FocusingTypeSet(ftype);
  return 0;
}

int ApertureCalcCwFieldRef(Aperture* obj, float* pos, size_t nPositions, size_t nDim,
                           SPS_FCOMPLEX** odata, size_t* nOutPositions) {
  auto p = (fnm::Aperture<float>*) obj;
  p->CalcCwFieldRef(pos, nPositions, nDim,
                    (std::complex<float>**)odata, nOutPositions);
  return 0;
}

int ApertureCalcCwFieldFast(Aperture* obj, float* pos, size_t nPositions, size_t nDim,
                            SPS_FCOMPLEX** odata, size_t* nOutPositions) {
  auto p = (fnm::Aperture<float>*) obj;
  p->CalcCwFast(pos, nPositions, nDim,
                (std::complex<float>**)odata, nOutPositions);
  return 0;
}
