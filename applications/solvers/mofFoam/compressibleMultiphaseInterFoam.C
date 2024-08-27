/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    compressibleMultiphaseInterFoam

Description
    Solver for n compressible, non-isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e.  laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "compressibleMultiphaseMixture.H"
#include "compressibleMomentumTransportModels.H"
#include "pimpleControl.H"
#include "fvModels.H"
#include "fvConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    volScalarField& p = mixture.p();
    volScalarField& T = mixture.T();

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;


        // --- Map of the massoffluidFoam solver
        /**
         * 
         * Overview of solution algorithm within a time step 𝑘. 
            - for each time-step 𝑘, do: 
                - 0 Update time increment $𝛥𝑡^𝑘$  based on Eq. (45) 
                - 1 Phase and energy advection 
                    - 1a Solve Eq. (46) for $𝜌^𝑘_𝑖,𝑐𝑜𝑛𝑣$ using $𝒖^{𝑘−1}$  and $𝛼_𝑖^{𝑘−1}$ 
                    - 1b Solve Eq. (48) for $𝐻_{𝑖,𝑐𝑜𝑛𝑣}^𝑘$ using $𝒖^{𝑘−1}$, $𝑄^{𝑘−1}_{𝑎𝑏𝑠,𝑖}$, and $𝑆^{𝑘−1}_{𝑝,𝑖}$  
                    - 1c Calculate $𝑇^𝑘_{𝑐𝑜𝑛𝑣}$ from $𝐻^𝑘_{𝑐𝑜𝑛𝑣}$ as laid out in Section 3.2 
                    - 1d Phase-norming by minimizing Eq. (58) for $𝛼^𝑘_{𝑖,𝑐𝑜𝑛𝑣}$
                - 2 Beam propagation 
                    - 2a Solve steady-state RTEs (39) and (40) iteratively using $𝛼^𝑘_{𝑖,𝑐𝑜𝑛𝑣}$ , where in each iteration $j$, (39) is solved for $𝐼^𝑔$ using $𝐼^{𝑗−1}_{𝑐𝑜𝑛𝑑}$  and subsequently, (40) is solved for $𝐼_{𝑐𝑜𝑛𝑑}$ using $𝐼^𝑗_𝑔$ , until $|(𝐼_𝑔+𝐼_{𝑐𝑜𝑛𝑑})^𝑗−(𝐼_𝑔+𝐼_{𝑐𝑜𝑛𝑑})_{𝑗−1} |/ (𝐼_𝑔+𝐼_{𝑐𝑜𝑛𝑑})_{𝑗−1} < 𝛿$ 
                    - 2b Calculate $𝑄^𝑘_{𝑎𝑏𝑠,𝐼}$ from $𝐼^𝑘_{𝑐𝑜𝑛𝑑}$ 
                    - 2c Calculate $𝑄^𝑘_{𝑎𝑏𝑠,𝑀𝑅}$ via ray tracing algorithm using $𝑄^𝑘_{𝑎𝑏𝑠,𝐼}$ and $𝛼^𝑘_{𝑖,𝑐𝑜𝑛𝑣}$ 
                    - 2d Calculate absorbed laser energy $𝑄^𝑘_{𝑎𝑏𝑠}=𝑄^𝑘_{𝑎𝑏𝑠,𝐼}+𝑄^𝑘_{𝑎𝑏𝑠,𝑀𝑅}$ 
                - 3 Heat conduction 
                    - 3a Solve Eq. (20) for $𝑇^𝑘_{𝑐𝑜𝑛𝑑𝑢𝑐𝑡}$ from $𝑇^𝑘_{𝑐𝑜𝑛𝑣}$ 
                    - 3b Update $𝐻^𝑘_{𝑖,𝑐𝑜𝑛𝑑𝑢𝑐𝑡}$ from $𝑇^𝑘_{𝑐𝑜𝑛𝑑𝑢𝑐𝑡}$ via Eq. (50) 
                - 4 Phase change 
                    - 4a Evaluate Eqs. (65)–(68) using $𝑇^𝑘_{𝑐𝑜𝑛𝑑𝑢𝑐𝑡}$ and $𝐻^𝑘_{𝑖,𝑐𝑜𝑛𝑑𝑢𝑐𝑡}$ 
                    - 4b Obtain $𝜌^𝑘_𝑖$ via Eq. (63) 
                    - 4c Obtain $𝐻^𝑘_𝑖$ and thus $𝐻^𝑘$ via Eq. (70) and obtain $𝑇^𝑘$ from $𝐻^𝑘$ as laid out in Section 3.2 
                    - 4d Phase-norming by minimizing Eq. (58) for $𝛼^𝑘_𝑖$ 
                    - 4e Calculate $𝛹^𝑘$ from Eq. (62) using $𝛼^𝑘_𝑖$ 
                - 5 PISO-Loop 
                    - 5a Evaluate momentum predictor (75) using current values 
                    - 5b Solve pressure equation (77) 
                    - 5c Calculate $𝒖^𝑘$ from Eq. (76) , calculate $𝑆^𝑘_{𝑃,𝑖}$ for use in time-step $𝑘 + 1$ 
                - end time-step
         * 
         */ 

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            fvModels.correct();

            mixture.solve();

            #include "contErr.H"

            #include "UEqn.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
