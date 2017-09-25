/**
 * @file   circular.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Mon Jul 24 15:13:32 2017
 *
 * @brief
 *
 *
 */

#pragma once

/*
 *  This file is part of SOFUS.
 *
 *  SOFUS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SOFUS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with SOFUS.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <sps/config.h>
#include <fnm/config.h>
#include <fnm/fnm_export.h>
#include <fnm/fnm_types.hpp> // Sysparm

#include <cstddef> // size_t
#include <complex>

namespace fnm {
  template<class T>
  struct CircularApertureData;
}

namespace fnm {
  template <class T>
  class FNM_EXPORT CircularAperture {
  public:
    CircularAperture();
    CircularAperture(T radius);
    ~CircularAperture();

    /**
     * Get max number of angular abcissas
     *
     *
     * @return # of abcissa in the width dimension
     */
    const size_t& NMaxDivAGet() const;

    const T& DensityGet() const;
    void DensitySet(const T& rho);

    void RadiusSet(const T& radius);
    const T& RadiusGet() const;

    void DelaySet(const T& delay);

    const T& DelayGet() const;
    /**
     * Set max number of angular abcissas
     *
     * @param nDivW
     */
    int NMaxDivASet(const size_t& nDivW);
    /**
     * Get center frequency
     *
     * @return
     */
    const T& F0Get() const;

    /**
     * Set center frequency
     *
     * @param f0
     */
    void F0Set(const T& f0);

    void FsSet(const T& fs);

    const T& FsGet() const;

    void CSet(const T& fs);

    const T& CGet() const;

    /**
     * Get width of pulse [s]
     *
     * @return
     */
    const T& WGet() const;

    /**
     * Set width of pulse [s]. TODO: Figure this out
     *
     * @param w
     */
    void WSet(const T& w);

    const T& GridSectorScaleGet() const;

    /**
     * Set center frequency
     *
     * @param gss
     */
    void GridSectorScaleSet(const T& gss);

    /**
     *
     * \f$
     * p(r,z;\omega) = v(w) rho_0 * a * c / pi * int
     * \f$
     *
     * \f$
     * H(r,z;k) = a / (ik\pi) int
     * \f$
     *
     * \f$
     * p(r,z,k) = -i v(\omega) * \rho_0 \omega H
     * \f$
     *
     *
     * @param pos
     * @param nPositions
     * @param nDim
     * @param odata
     * @param nOutPositions
     */
    int CalcCwFieldRef(const T* pos, const size_t nPositions, const size_t nDim,
                       std::complex<T>** odata, size_t* nOutPositions);

    /**
     * Compute transient field at the specified positions. The pulse is a simple
     * tone-burst specified by setting the length w and center frequency f0.
     *
     * auto a = CircularAperture<float>(radius);
     * a.F0Set(1e6);
     * a.WSet(nCycles / 1e6)
     *
     * a.CalcTransientFieldRef(pPositions, nPositions, 3, ppOutData, pNSignals, pNSamples)
     *
     * @param pos
     * @param nPositions
     * @param nDim
     * @param odata
     * @param nSignals
     * @param nSamples
     *
     * @return
     */
    T CalcTransientFieldRef(const T* pos, const size_t nPositions, const size_t nDim,
                            T** odata, size_t* nSignals, size_t* nSamples);


  private:
    ///< Data
    CircularApertureData<T>* m_data;

    ///< Sysparm
    sysparm_t<T>* m_sysparm;

  };
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
