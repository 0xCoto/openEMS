/*
*	Copyright (C) 2025 Gadi Lahav (gadi@rfwithcare.com)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*/

#ifndef ENGINE_EXT_MODEABSORB_H
#define ENGINE_EXT_MODEABSORB_H

#include <vector>

#include "engine_extension.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"
#include "engine_extension_dispatcher.h"

class Operator_Ext_ModeAbsorb;

// File-based modal absorber engine. Mirrors Engine_Ext_WaveguideAbsorber's
// matched-modal-source kernel: store H^n at the H-plane in the pre-current
// hook, then in the post-current hook compute the modal projections,
// directionally decompose into a_bwd, and subtractively write back to E
// and H.

class Engine_Ext_ModeAbsorb : public Engine_Extension
{
public:
	Engine_Ext_ModeAbsorb(Operator_Ext_ModeAbsorb* op_ext);
	virtual ~Engine_Ext_ModeAbsorb();

	virtual void DoPreCurrentUpdates() {Engine_Ext_ModeAbsorb::DoPreCurrentUpdates(0);}
	virtual void DoPreCurrentUpdates(int threadID);

	virtual void Apply2Voltages() {Engine_Ext_ModeAbsorb::Apply2Voltages(0);}
	virtual void Apply2Voltages(int threadID);

	virtual void Apply2Current() {Engine_Ext_ModeAbsorb::Apply2Current(0);}
	virtual void Apply2Current(int threadID);

protected:
	Operator_Ext_ModeAbsorb* m_Op_MA;

	template <typename EngType>
	void DoPreCurrentUpdatesImpl(EngType* eng, int threadID);

	template <typename EngType>
	void Apply2VoltagesImpl(EngType* eng, int threadID);

	template <typename EngType>
	void Apply2CurrentImpl(EngType* eng, int threadID);

	int m_ny, m_nyP, m_nyPP;
	int m_dir;

	unsigned int m_posStart[3];

	unsigned int m_numLines_E[2];
	double** m_E_OverlapW[2];
	double** m_E_SubtractW[2];

	unsigned int m_numLines_H[2];
	double** m_H_OverlapW[2];
	double** m_H_SubtractW[2];

	unsigned int m_E_abs_line;
	unsigned int m_H_abs_line;

	// H^n storage for temporal averaging (single layer): [comp][posP][posPP]
	double** m_H_stored[2];
	// 2026-05-11 Phase J: Yee stagger correction.  H lives at z = (k+½)·dz
	// while E lives at z = k·dz.  Sampling H at a single plane gives a dz/2
	// spatial offset → tan(β·dz/2) phase error in V/I formula → ~−13 dB S11
	// floor at 2.5 GHz.  We average H at index m_H_abs_line and index
	// (m_H_abs_line − 1) to interpolate H to the E plane.  m_H_stored_alt
	// holds the previous-timestep H at the adjacent plane.
	double** m_H_stored_alt[2];

	double m_Zw;

	// --- subtract-known-incident additions (CST-style) ---
	// During the first N_CAL_STEPS, no source-emission subtraction is applied
	// and (Emm, Signal[n]) pairs are accumulated to derive K_eff via
	// least-squares. After calibration completes, the absorber subtracts the
	// known incident from Emm before the V/I decomposition:
	//      a_inc_known = 0.5 * Signal[n] * K_eff * port_amplitude
	//      Emm <- Emm - a_inc_known
	// removing the source's immediate forward+backward emission. The V/I
	// decomposition then operates on the residual (mostly true reflections /
	// non-modal noise) and gives a better-behaved a_bwd.
	double       m_PortAmplitude;
	unsigned int m_TS;            // running timestep count (incremented per Apply2Current)
	bool         m_Calibrated;    // true after K_eff fit completes
	double       m_K_eff;         // calibrated incident-to-Emm conversion
	double       m_calib_num;     // sum(Emm * Signal)
	double       m_calib_den;     // sum(Signal * Signal)

	// --- state-space recursion (CST CBBPortOperatorStateSpace path) ---
	// Persistent ψ (size m_Op_MA->m_SS_n) holding the state vector across
	// time steps.  Update per CN with u_{n+1}≈u_n approximation:
	//      ψ_new = P · ψ + Q · u_n
	//      y_n   = C · ψ_new + D · u_n
	//      correction = (u_n − y_n) · e_mode   (subtracted from E at Π)
	std::vector<double> m_psi;
	std::vector<double> m_psi_new;

	// --- Luo-Chen 1-D Yee per-mode absorber (Phase C, 2026-05-10) ---
	// Klein-Gordon update for modal dispersion:
	//      V[k]^{n+1} = 2·V[k]^n − V[k]^{n-1}
	//                 + A·(V[k+1]−2V[k]+V[k-1])^n − B·V[k]^n
	// where A = (c·dt/dz_line)², B = (ω_c·dt)²
	// gives ω² = c²·β² + ω_c² (correct waveguide TE dispersion).
	// Drive: V[0]^n = Emm at start of step.  V[-1] = V_minus1 from 3-D.
	// Mur ABC at k=N-1.  Replace V_3D modal V with V[0]^{n+1}.
	bool   m_LC_Active;
	int    m_LC_N;
	double m_LC_dz_line;       // 1-D line cell pitch (independent of 3-D)
	double m_LC_A;             // (c·dt/dz_line)² — KG spatial coefficient
	double m_LC_B;             // (ω_c·dt)²       — KG cutoff coefficient
	double m_LC_mur_coef;      // (c·dt − dz_line)/(c·dt + dz_line)
	std::vector<double> m_LC_V;        // V[k]^n           (length N)
	std::vector<double> m_LC_V_prev;   // V[k]^{n-1}       (length N)
	double m_LC_V_prev_far;            // V[N-1]^{n-1} (for Mur)
	double m_LC_V_prev_near;           // V[N-2]^{n-1} (for Mur)

	// --- CFS-CPML modal aux fields (Phase O, 2026-05-11) ---
	// Modal CFS-CPML at cells inside WG just before port plane.
	// For each cell k=0..N_cpml-1 (k=0 at port plane, k=N-1 deep inside):
	//   ψ_V_modal[k] tracks the past ∂I_modal/∂z via:
	//     ψ_V_modal[k]^(n+1) = b_z[k]·ψ_V_modal[k]^n + a_z[k]·∂I_modal/∂z
	//   V_modal correction: V_modal[k] -= (dt/ε)·ψ_V_modal[k]^(n+1)
	// Similarly for ψ_I_modal[k] with ∂V_modal/∂z.
	// b_z, a_z computed from σ_z(k), κ_z(k), α_z(k) profiles.
	// Enable via MODE_CPML_ENABLE=1; default OFF.
	bool                  m_CPML_Active;
	int                   m_CPML_N;       // number of CPML cells
	std::vector<double>   m_CPML_b_z;     // size N_cpml
	std::vector<double>   m_CPML_a_z;     // size N_cpml
	std::vector<double>   m_CPML_kappa;   // size N_cpml
	std::vector<double>   m_psi_V_modal;  // [k]: aux for V_modal (E side)
	std::vector<double>   m_psi_I_modal;  // [k]: aux for I_modal (H side)

	// --- Wang 2022 wave-equation CFS-PML for 1-D Klein-Gordon line ---
	// Path 1 implementation (2026-05-11).  Replaces the Mur ABC at the
	// far end of the Luo-Chen 1-D KG line with a CFS-PML region of N_pml
	// cells.  Paper: IEEE MWCL 2022 "Convolutional Implementation and
	// Analysis of the CFS-PML ABC for the FDTD Method Based on Wave
	// Equation".  Optimal params: m=4, κ_max=4, α_max=0.1,
	// σ_opt=(m+1)/(150π·Δs).
	//
	// PML occupies cells k=[N − N_pml, N − 1] of the 1-D line.  Each PML
	// cell holds:
	//   ξ[k+1/2]: first-level aux for stretched ∂V/∂z at staggered z+½dz
	//   φ[k]    : second-level aux for stretched ∂²V/∂z² at cell centre
	//   σ_z[k], α_z[k], κ_z[k]: profiles
	//   b_z[k], c_z[k]: derived recursion coefficients
	// PEC (V=0) at k=N−1 closes the line.  Enable via USE_WANG_CFS_PML=1.
	bool                  m_WP_Active;
	int                   m_WP_N_pml;     // CFS-PML cells at far end
	double                m_WP_kappa_max;
	double                m_WP_alpha_max;
	double                m_WP_sigma_max;
	int                   m_WP_p_order;
	std::vector<double>   m_WP_kappa_half;   // κ at k+½ (size N_pml+1)
	std::vector<double>   m_WP_kappa_cell;   // κ at k   (size N_pml)
	std::vector<double>   m_WP_b_half;       // b at k+½
	std::vector<double>   m_WP_c_half;       // c at k+½
	std::vector<double>   m_WP_b_cell;       // b at k
	std::vector<double>   m_WP_c_cell;       // c at k
	std::vector<double>   m_WP_xi;           // ξ[k+½] (size N_pml+1)
	std::vector<double>   m_WP_phi;          // φ[k]   (size N_pml)
};

#endif // ENGINE_EXT_MODEABSORB_H
