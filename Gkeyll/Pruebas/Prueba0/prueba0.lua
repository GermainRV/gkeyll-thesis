-- Gkyl ------------------------------------------------------------------------
-- 1x1v simulation of collisionless damping of an electron Langmuir wave and
-- ion acustic wave, using kinetic ions and electrons. We use normalize units.
--------------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()
-- Pseudo-random number generator from SciLua package for initial conditions.
-- Specific rng is the mrg32k3a (multiple recursive generator) of Ecuyer99.
local prng = require "sci.prng"
local rng = prng.mrg32k3a(123456)
--------------------------------------------------------------------------------
-- Parameters for plasma simulation
--------------------------------------------------------------------------------
-- UNIVERSAL CONSTANTS
epsilon_0 = 1.0			-- Permitivity of free space
mu_0 	  = 1			-- Permeability of free space
k_B 	  = 1			-- Boltzmann constant
e 		  = 1.0			-- Elementary charge
c 		  = 1			-- Speed of light
-- PARAMETERS RATIO
denRatio  = 1			-- Ion-electron density ratio
massRatio = 100			-- Ion-electron mass ratio
tempRatio = 0.1			-- Ion-electron temperature ratio
velRatio  = 0.05			-- Electron thermal velocity-light speed ratio
-- PLASMA PARAMETERS
q_Elc 	 = -e							-- Electron charge
q_Ion 	 = e							-- Ion charge
m_Elc, m_Ion = 1, massRatio				-- Electron and Ion mass
n_Elc, n_Ion = 1, denRatio				-- Electron and Ion number density
T_Elc 	 = 1							-- Electron temperature
T_Ion 	 = tempRatio*T_Elc				-- Ion temperature
vth_Elc  = math.sqrt(k_B*T_Elc/m_Elc)	-- Electron thermal speed
vth_Ion	 = math.sqrt(k_B*T_Ion/m_Ion)	-- Ion thermal speed
w_pe	 = math.sqrt((q_Elc^2)*n_Elc/(epsilon_0*m_Elc))	-- Electron plasma frequency
w_pi	 = math.sqrt((q_Ion^2)*n_Ion/(epsilon_0*m_Ion))	-- Ion plasma frequency
lambda_De = vth_Elc/w_pe					-- Electron debye length
ux_Elc, uy_Elc, uz_Elc = 0.0, 0.0, 0.0	-- Bulk flow velocity of the electrons
ux_Ion, uy_Ion, uz_Ion = 0.0, 0.0, 0.0	-- Bulk flow velocity of the ions
-- SIMULATION PARAMETERS
Lx = 500*lambda_De		-- Plasma length
Nx = 64					-- Number of points in x
Nv = 64					-- Number of points in v 
dx = 2*Lx/(Nx-1)    	-- Step size in x
--dk = 1/(Nx*dx) 			-- Step size in k
dk = 2*math.pi/(2*Lx)

-- Maxwellian in (r,v)-space, given the density (den), bulk flow
-- velocity (u_r), and thermal speed (vth).
DIMV = 1
local function maxwellian1D1V(x, vx, ux, den, vth)
	local v2   = (vx - ux)^2
	return den/(math.sqrt(2*math.pi*vth^2)^DIMV)*math.exp(-v2/(2*vth^2))
end

A = 1e-4 -- Amplitude of sinusoidal perturbation. E = 0.1 mV/m
-- Building random array to excite all wavemodes
-- Perturbing the first N wave modes with random amplitudes and phases.
denk   = 1
Nmodes = 54*denk--32*denK
P={}
for i = 1,Nmodes,1 do
   P[i] = rng:sample()
end

plasmaApp = Plasma.App {
	logToFile = true,
	tEnd         = 5000/w_pe,	-- End time.
	nFrame       = 2500,		-- Number of output frames.
	lower        = {-Lx},		-- Lower boundary of configuration space.
	upper        = {Lx},		-- Upper boundary of configuration space.
	cells        = {Nx},		-- Configuration space cells.
	polyOrder    = 2,           -- Polynomial order.
	periodicDirs = {1},         -- Periodic directions.
	decompCuts = {1},         	-- cuts in each configuration direction
   
   elc = Plasma.Species {
		charge = q_Elc, mass = mElc,
		lower = {-6.0*vth_Elc},     -- Velocity space lower boundary.
		upper = { 6.0*vth_Elc},     -- Velocity space upper boundary.
		cells = {Nv},              	-- Number of cells in velocity space.
		init = function (t, xn)    	-- Initial conditions.
			local x, vx = xn[1], xn[2]
			local rho_0 = 0
			for i=1,Nmodes,1 do
				local k = i*dk/denk
				rho_0 = rho_0 + A*(epsilon_0*k/n_Elc)*math.cos(k*x + 2*math.pi*P[i])
			end
			return maxwellian1D1V(x,vx,ux_Elc,n_Elc,vth_Elc)*(1 + rho_0)
		end,
		evolve = true, 				-- Evolve species
		diagnostics = { "M0", "M1i", "M2" },
   },

   ion = Plasma.Species {
		charge = q_Ion, mass = m_Ion,
		lower = {-6.0*vth_Ion},     -- Velocity space lower boundary.
		upper = { 6.0*vth_Ion},     -- Velocity space upper boundary.
		cells = {Nv},              	-- Number of cells in velocity space.
		init  = function (t, xn)   	-- Initial conditions.
			local x, vx = xn[1], xn[2]
			return maxwellian1D1V(x,vx,ux_Ion,n_Ion,vth_Ion)
		end,
		evolve = true, 				-- Evolve species
   },

   field = Plasma.Field {
		epsilon0 = epsilon_0, mu0 = mu_0,
		init = function (t, xn)		-- Initial conditions.
			local x, vx = xn[1], xn[2]
			local Ex = 0.0
			for i=1,Nmodes,1 do
				-- local k = 2*math.pi*i*dk/denK
				local k = i*dk/denk
				--Ex = Ex + A*q_Elc*n_Elc/(k*epsilon_0)*math.sin(k*x + 2*math.pi*P[i])/Nmodes
				Ex = Ex + A*math.sin(k*x + 2*math.pi*P[i])
			end
			return Ex, 0.0, 0.0, 0.0, 0.0, 0.0
		end,
		evolve = true, 				-- Evolve field
   },
}
-- Run application.
plasmaApp:run()
